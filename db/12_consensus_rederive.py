#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 12_consensus_rederive.py — re-derive per-copy variants against the population-majority consensus.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
12_consensus_rederive.py

Re-derive the entire variant table against the POPULATION-MAJORITY CONSENSUS
(not CHM13), indel-aware, with masking of SNPs adjacent to indels.

Two passes over each haplotype's full-unit MAFFT alignment + mafft --add T2T:

  PASS 1 — build global consensus:
    For every T2T position (1..2168), tally the aligned base across ALL copies
    of ALL haplotypes; consensus[p] = majority base. Cache the T2T-augmented
    alignment to disk so PASS 2 does not re-run mafft --add.

  PASS 2 — re-derive each copy vs consensus:
    * SNP:  copy base != consensus base (both non-gap)
    * del:  copy gap where T2T/consensus has a base (anchor = preceding T2T pos)
    * ins:  copy base where T2T has a gap   (anchor = preceding T2T pos)
    * mask: any SNP whose T2T position is within MASK_DIST of an indel anchor
    Stored as alignment_source='consensus_t2t' with columns:
      consensus_pos, ref(=consensus base), alt(=copy base),
      chm13_base, var_type('snp'|'ins'|'del'), masked(0|1), region

Existing gene_unit_t2t rows are left intact as a fallback.

Paths are read from environment variables:
  FIVES_DB    path to 5S_rDNA.db
  FIVES_DATA  input derived-data directory (sequences/, cache, consensus output)
  FIVES_REFS  reference directory (consensus FASTA)

Usage:
  python3 12_consensus_rederive.py [--dry-run] [--sample HG00097] [--jobs 20]
  python3 12_consensus_rederive.py --pass1-only   # just build + save consensus
"""

import sys, os, re, json, subprocess, tempfile, sqlite3
from pathlib import Path
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

# ── Paths (from environment variables) ────────────────────────────────────────

ROOT = Path(os.environ.get("FIVES_DATA", "data"))

HPRC     = ROOT
DB_PATH  = Path(os.environ.get("FIVES_DB", "5S_rDNA.db"))
T2T_FA   = Path(os.environ.get("FIVES_REFS", "refs")) / "5S_t2t_consensus.fa"
SEQ_DIR  = HPRC / "sequences"
CONS_OUT = HPRC / "consensus_reference.json"   # saved consensus + per-pos table

MAFFT = "mafft"

# Also handle the alternate CHM13/HG002 alignment file locations
ALT_ALN = {
    ("CHM13", "hap1"):      ROOT / "sequences/CHM13/5S_array_aligned.fa",
    ("HG002_GIAB", "hap1"): ROOT / "sequences/HG002/hg002v1_MATERNAL_gene_aln_full.fa",
    ("HG002_GIAB", "hap2"): ROOT / "sequences/HG002/hg002v1_PATERNAL_gene_aln_full.fa",
}

COMP = str.maketrans("ACGTacgt", "TGCAtgca")
MASK_DIST            = 5     # bp; SNP within this of an indel anchor → masked
ORIENT_SNP_THRESHOLD = 100
MAJORITY_MIN_COUNT   = 10    # min aligned copies to call a consensus base (else use CHM13)


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args():
    a = sys.argv[1:]
    return {
        "dry_run":    "--dry-run" in a,
        "pass1_only": "--pass1-only" in a,
        "sample":     _opt(a, "--sample"),
        "jobs":       int(_opt(a, "--jobs") or 20),
    }

def _opt(args, flag):
    if flag in args:
        i = args.index(flag)
        if i + 1 < len(args): return args[i+1]
    return None


# ── FASTA (no BioPython) ──────────────────────────────────────────────────────

def parse_fasta(path) -> dict:
    seqs, h, buf = {}, None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if h: seqs[h] = "".join(buf)
                h, buf = line[1:].split()[0], []
            elif h:
                buf.append(line.upper())
    if h: seqs[h] = "".join(buf)
    return seqs

def parse_fasta_str(text) -> dict:
    seqs, h, buf = {}, None, []
    for line in text.splitlines():
        line = line.rstrip()
        if line.startswith(">"):
            if h: seqs[h] = "".join(buf)
            h, buf = line[1:].split()[0], []
        elif h:
            buf.append(line.upper())
    if h: seqs[h] = "".join(buf)
    return seqs


# ── mafft --add T2T, cached ───────────────────────────────────────────────────

def aln_path_for(sample_id, hap_label):
    alt = ALT_ALN.get((sample_id, hap_label))
    if alt and alt.exists():
        return alt
    p = SEQ_DIR / sample_id / f"{sample_id}_{hap_label}_gene_aln.fa"
    return p if p.exists() else None


def augmented_alignment(sample_id, hap_label, t2t_seq, cache_dir):
    """
    Run mafft --add T2T on the haplotype alignment, auto-detecting orientation.
    Returns (t2t_row, {copy_num: aligned_seq}, orientation).
    Caches the augmented FASTA to cache_dir for reuse across passes.
    """
    cache = cache_dir / f"{sample_id}_{hap_label}_t2taug.fa"
    orient_cache = cache_dir / f"{sample_id}_{hap_label}_orient.txt"

    if cache.exists() and orient_cache.exists():
        aug = parse_fasta(cache)
        orient = orient_cache.read_text().strip()
        t2t_row = aug.pop("T2T_ref", None) or aug.pop("T2T_rc", None)
        copies = {int(re.match(r"copy(\d+)", k).group(1)): v
                  for k, v in aug.items() if k.startswith("copy")}
        return t2t_row, copies, orient

    aln = aln_path_for(sample_id, hap_label)
    if aln is None:
        raise FileNotFoundError(f"no alignment for {sample_id} {hap_label}")

    def run_add(ref_seq, ref_id):
        with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as tmp:
            tmp.write(f">{ref_id}\n{ref_seq}\n"); fn = tmp.name
        try:
            r = subprocess.run([MAFFT, "--add", fn, "--quiet", str(aln)],
                               capture_output=True, text=True, check=True)
        finally:
            os.unlink(fn)
        return parse_fasta_str(r.stdout)

    aug = run_add(t2t_seq, "T2T_ref")
    t2t_row = aug.get("T2T_ref")
    copies = {int(re.match(r"copy(\d+)", k).group(1)): v
              for k, v in aug.items() if k.startswith("copy")}
    sample = list(copies.values())[:10]
    avg = sum(sum(1 for a, b in zip(t2t_row, s) if a != "-" and b != "-" and a != b)
              for s in sample) / max(len(sample), 1)
    orient = "fwd"
    if avg > ORIENT_SNP_THRESHOLD:
        t2t_rc = t2t_seq.translate(COMP)[::-1]
        aug = run_add(t2t_rc, "T2T_rc")
        t2t_row = aug.get("T2T_rc")
        copies = {int(re.match(r"copy(\d+)", k).group(1)): v
                  for k, v in aug.items() if k.startswith("copy")}
        orient = "rc"

    # cache
    with open(cache, "w") as f:
        rid = "T2T_rc" if orient == "rc" else "T2T_ref"
        f.write(f">{rid}\n{t2t_row}\n")
        for cn, s in copies.items():
            f.write(f">copy{cn:03d}\n{s}\n")
    orient_cache.write_text(orient)
    return t2t_row, copies, orient


# ── Column → T2T position map ─────────────────────────────────────────────────

def column_to_t2t(t2t_row, t2t_len, orient):
    """
    Return dict {alignment_col: t2t_pos_1based} for non-gap T2T columns.
    For rc orientation, positions are mapped back to forward coordinates.
    """
    col2pos = {}
    c = 0
    for col, base in enumerate(t2t_row):
        if base != "-":
            c += 1
            if orient == "rc":
                col2pos[col] = t2t_len - c + 1
            else:
                col2pos[col] = c
    return col2pos


def fwd_base(base, orient):
    """Convert an aligned base to forward-strand orientation."""
    if base == "-":
        return "-"
    return base.translate(COMP) if orient == "rc" else base


# ── PASS 1: accumulate consensus tallies ──────────────────────────────────────

def pass1_tally(task):
    sample_id, hap_label, t2t_seq, cache_dir = (
        task["sample_id"], task["hap_label"], task["t2t_seq"], task["cache_dir"])
    T2T_LEN = len(t2t_seq)
    try:
        t2t_row, copies, orient = augmented_alignment(
            sample_id, hap_label, t2t_seq, cache_dir)
    except Exception as e:
        return {"err": f"{sample_id} {hap_label}: {e}"}

    col2pos = column_to_t2t(t2t_row, T2T_LEN, orient)
    # tally[pos][base] += 1 ; gaps counted as '-' so we can detect majority-deletions
    tally = defaultdict(lambda: defaultdict(int))
    for cn, seq in copies.items():
        for col, pos in col2pos.items():
            b = fwd_base(seq[col], orient)
            if b in "ACGT":
                tally[pos][b] += 1
            elif b == "-":
                tally[pos]["-"] += 1
    # convert to plain dict for pickling
    return {"tally": {p: dict(d) for p, d in tally.items()}}


# ── PASS 2: re-derive vs consensus, indel-aware ───────────────────────────────

def pass2_derive(task):
    sample_id  = task["sample_id"]
    hap_label  = task["hap_label"]
    hap_id     = task["haplotype_id"]
    cid_map    = task["copy_id_map"]
    t2t_seq    = task["t2t_seq"]
    consensus  = task["consensus"]      # list, index 0 = pos1
    gap_pos    = task.get("gap_positions", set())  # positions where ref = deletion
    cache_dir  = task["cache_dir"]
    T2T_LEN = len(t2t_seq)

    out = {"sample_id": sample_id, "hap_label": hap_label, "haplotype_id": hap_id,
           "status": "ok", "rows": [], "orient": "?", "n_snp": 0,
           "n_indel": 0, "n_masked": 0, "error": None}
    try:
        t2t_row, copies, orient = augmented_alignment(
            sample_id, hap_label, t2t_seq, cache_dir)
        out["orient"] = orient
        col2pos = column_to_t2t(t2t_row, T2T_LEN, orient)
        ncols = len(t2t_row)

        def region_of(p):
            if p < 630: return "nts_pre"
            if p <= 748: return "gene"
            return "nts_post"

        for cn, seq in copies.items():
            cid = cid_map.get(cn)
            if cid is None:
                continue
            # Walk the alignment in forward T2T order.
            # Build ordered list of (col, pos) for T2T-anchored columns.
            # For insertion columns (T2T gap), anchor to the nearest preceding pos.
            variants = []           # (pos, ref, alt, vtype, region)
            indel_positions = []    # T2T positions adjacent to an indel

            # iterate columns in alignment order; track running t2t position
            # For rc we still iterate columns left→right; positions decrease, so
            # we collect then handle masking by absolute pos distance afterward.
            last_pos = 0
            ins_run_anchor = None
            ins_bases = []

            # We process per-column. Deletions: T2T base present, copy gap.
            # Insertions: T2T gap (col not in col2pos), copy base present.
            for col in range(ncols):
                t2t_b = fwd_base(t2t_row[col], orient)
                cp_b  = fwd_base(seq[col], orient)
                if col in col2pos:
                    pos = col2pos[col]
                    cons_b = consensus[pos - 1] if 1 <= pos <= T2T_LEN else t2t_b
                    chm_b  = t2t_seq[pos - 1]
                    if pos in gap_pos:
                        # Reference state here is DELETION (majority of copies lack base).
                        # Invert polarity: a copy WITH a base is the insertion variant;
                        # a copy with a gap matches the reference (no variant). No SNPs.
                        if cp_b in "ACGT":
                            variants.append((pos, "-", cp_b, "ins", region_of(pos)))
                            indel_positions.append(pos)
                        last_pos = pos
                    elif cp_b == "-":
                        # deletion at this T2T position
                        variants.append((pos, cons_b, "-", "del", region_of(pos)))
                        indel_positions.append(pos)
                        last_pos = pos
                    elif cp_b in "ACGT" and cp_b != cons_b:
                        variants.append((pos, cons_b, cp_b, "snp", region_of(pos)))
                        last_pos = pos
                    else:
                        last_pos = pos
                else:
                    # insertion column (T2T gap)
                    if cp_b in "ACGT":
                        anchor = last_pos if last_pos > 0 else 1
                        variants.append((anchor, "-", cp_b, "ins", region_of(anchor)))
                        indel_positions.append(anchor)

            # Mask SNPs within MASK_DIST of any indel position
            indel_set = set(indel_positions)
            for (pos, ref, alt, vtype, region) in variants:
                masked = 0
                if vtype == "snp":
                    if any(abs(pos - ip) <= MASK_DIST for ip in indel_set):
                        masked = 1
                chm_b = t2t_seq[pos - 1] if 1 <= pos <= T2T_LEN else "N"
                out["rows"].append(
                    (cid, "consensus_t2t", int(pos), ref, alt, region,
                     chm_b, vtype, masked))
                if vtype == "snp":
                    out["n_snp"] += 1
                    out["n_masked"] += masked
                else:
                    out["n_indel"] += 1

    except Exception as e:
        out["status"] = "error"
        out["error"] = str(e)
    return out


# ── DB helpers ────────────────────────────────────────────────────────────────

def db_connect():
    con = sqlite3.connect(DB_PATH)
    con.execute("PRAGMA foreign_keys = ON")
    con.execute("PRAGMA journal_mode = WAL")
    return con

def ensure_columns(con):
    cols = {r[1] for r in con.execute("PRAGMA table_info(variant)").fetchall()}
    if "chm13_base" not in cols:
        con.execute("ALTER TABLE variant ADD COLUMN chm13_base TEXT")
    if "var_type" not in cols:
        con.execute("ALTER TABLE variant ADD COLUMN var_type TEXT DEFAULT 'snp'")
    if "masked" not in cols:
        con.execute("ALTER TABLE variant ADD COLUMN masked INTEGER DEFAULT 0")
    con.commit()

def get_haplotypes(con, sample_filter):
    where = "a.cohort IN ('HPRC_Year1','HPRC_Release2','CHM13','HG002_GIAB')"
    params = []
    if sample_filter:
        where += " AND a.sample_id = ?"; params.append(sample_filter)
    rows = con.execute(f"""
        SELECT a.sample_id, h.hap_label, h.haplotype_id
        FROM assembly a JOIN haplotype h USING(assembly_id)
        WHERE {where} ORDER BY a.sample_id, h.hap_label
    """, params).fetchall()
    haps = []
    for sid, hl, hid in rows:
        cmap = dict(con.execute(
            "SELECT copy_number, copy_id FROM copy WHERE haplotype_id=?", (hid,)
        ).fetchall())
        haps.append({"sample_id": sid, "hap_label": hl,
                     "haplotype_id": hid, "copy_id_map": cmap})
    return haps


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()
    mode = "[DRY-RUN] " if args["dry_run"] else ""
    print(f"{mode}Consensus re-derivation (indel-aware, MASK_DIST={MASK_DIST})")
    print(f"  ROOT={ROOT}")

    t2t_seq = list(parse_fasta(T2T_FA).values())[0]
    T2T_LEN = len(t2t_seq)
    print(f"  T2T: {T2T_LEN} bp")

    con = db_connect()
    ensure_columns(con)
    haps = get_haplotypes(con, args["sample"])
    con.close()
    print(f"  Haplotypes: {len(haps)}")

    cache_dir = HPRC / "t2t_aug_cache"
    cache_dir.mkdir(exist_ok=True)

    # ── PASS 1 ────────────────────────────────────────────────────────────────
    print("\n=== PASS 1: build consensus ===")
    global_tally = defaultdict(lambda: defaultdict(int))
    tasks1 = [{"sample_id": h["sample_id"], "hap_label": h["hap_label"],
               "t2t_seq": t2t_seq, "cache_dir": cache_dir} for h in haps]

    with ProcessPoolExecutor(max_workers=args["jobs"]) as pool:
        for i, res in enumerate(as_completed([pool.submit(pass1_tally, t) for t in tasks1])):
            r = res.result()
            if "err" in r:
                print(f"  PASS1 error: {r['err']}")
                continue
            for pos, d in r["tally"].items():
                for b, n in d.items():
                    global_tally[pos][b] += n
            if (i + 1) % 50 == 0:
                print(f"  ...{i+1}/{len(tasks1)} haplotypes tallied", flush=True)

    # Build consensus sequence (T2T frame). The majority state at a position may
    # be a gap ('-') → the consensus is "deleted" there (most copies lack the base).
    # We keep the 2168 T2T coordinate frame but flag such positions so PASS 2
    # inverts the indel polarity (copies WITH a base become insertions).
    consensus = list(t2t_seq)
    n_flipped = 0
    flip_table = []
    consensus_gap_positions = []   # positions where deletion is the majority state
    for pos in range(1, T2T_LEN + 1):
        d = global_tally.get(pos, {})
        total = sum(d.values())
        if total < MAJORITY_MIN_COUNT:
            continue  # keep CHM13 base
        major = max(d, key=d.get)
        if major == "-":
            # Majority of copies have a deletion here → reference state is gap.
            consensus_gap_positions.append({
                "pos": pos, "chm13": t2t_seq[pos-1],
                "del_frac": d["-"] / total, "total": total,
            })
            # keep t2t base in `consensus` string as context placeholder
            continue
        if major != t2t_seq[pos - 1]:
            consensus[pos - 1] = major
            n_flipped += 1
            flip_table.append({
                "pos": pos, "chm13": t2t_seq[pos-1], "consensus": major,
                "major_frac": d[major] / total, "total": total,
            })
    consensus = "".join(consensus)
    print(f"  Consensus built: {n_flipped} SNP positions differ from CHM13")
    for ft in sorted(flip_table, key=lambda x: -x["major_frac"])[:12]:
        print(f"    pos {ft['pos']:4d}: CHM13 {ft['chm13']} → consensus {ft['consensus']} "
              f"({ft['major_frac']*100:.0f}% of {ft['total']} copies)")
    print(f"  Majority-deletion positions (reference = gap): {len(consensus_gap_positions)}")
    for gp in sorted(consensus_gap_positions, key=lambda x: -x["del_frac"]):
        print(f"    pos {gp['pos']:4d}: CHM13 {gp['chm13']} → consensus GAP "
              f"({gp['del_frac']*100:.0f}% of {gp['total']} copies deleted)")

    # Save consensus
    CONS_OUT.write_text(json.dumps({
        "consensus": consensus, "chm13": t2t_seq,
        "n_flipped": n_flipped, "flip_table": flip_table,
        "consensus_gap_positions": [g["pos"] for g in consensus_gap_positions],
        "consensus_gap_detail": consensus_gap_positions,
        "mask_dist": MASK_DIST,
    }, indent=2))
    print(f"  Saved consensus → {CONS_OUT}")
    gap_set = set(g["pos"] for g in consensus_gap_positions)

    if args["pass1_only"]:
        print("\n[pass1-only] stopping.")
        return

    # ── PASS 2 ────────────────────────────────────────────────────────────────
    print("\n=== PASS 2: re-derive vs consensus ===")
    tasks2 = [{"sample_id": h["sample_id"], "hap_label": h["hap_label"],
               "haplotype_id": h["haplotype_id"], "copy_id_map": h["copy_id_map"],
               "t2t_seq": t2t_seq, "consensus": consensus,
               "gap_positions": gap_set, "cache_dir": cache_dir}
              for h in haps]

    con = db_connect()
    # Clear any prior consensus_t2t rows
    if not args["dry_run"]:
        con.execute("DELETE FROM variant WHERE alignment_source='consensus_t2t'")
        con.commit()

    tot_snp = tot_indel = tot_masked = ok = err = 0
    with ProcessPoolExecutor(max_workers=args["jobs"]) as pool:
        for res in as_completed([pool.submit(pass2_derive, t) for t in tasks2]):
            r = res.result()
            if r["status"] != "ok":
                print(f"  PASS2 error {r['sample_id']} {r['hap_label']}: {r['error']}")
                err += 1
                continue
            if not args["dry_run"] and r["rows"]:
                con.executemany("""
                    INSERT INTO variant
                      (copy_id, alignment_source, consensus_pos, ref, alt, region,
                       chm13_base, var_type, masked)
                    VALUES (?,?,?,?,?,?,?,?,?)
                """, r["rows"])
                con.commit()
            tot_snp += r["n_snp"]; tot_indel += r["n_indel"]; tot_masked += r["n_masked"]
            ok += 1

    con.close()
    print(f"\n{mode}Done. OK={ok} err={err}")
    print(f"  SNPs: {tot_snp:,}  (masked near indel: {tot_masked:,})")
    print(f"  Indels: {tot_indel:,}")
    print(f"  Clean SNPs (unmasked): {tot_snp - tot_masked:,}")


if __name__ == "__main__":
    main()
