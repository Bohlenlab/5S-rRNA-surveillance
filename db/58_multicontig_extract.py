#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 58_multicontig_extract.py — extract the full 5S array for haplotypes split across multiple contigs.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
58 — Multi-contig 5S array extraction with flank-anchored ordering.

Recovers the full 5S array for haplotypes whose array is split across several
assembly contigs. Processes all array-fragment contigs, orders/orients them
against the CHM13 chr1q42 5S locus, and emits a complete, ordered per-haplotype
TSV.

Definitions
  array-fragment contig : >=3 5S copies in a tandem run (median spacing in
                          [1800,2600] bp). Others (1-2 hits) = dispersed/orphan
                          5S, recorded with array_member=0, category='dispersed_5S'.
  edge contig           : array abuts ONE end (break), large unique flank on the
                          other -> anchors a 5' or 3' array terminus.
  interior contig       : array abuts BOTH ends (no flank) -> goes between edges.

Ordering / orientation
  minimap2 each edge contig's flank to CHM13v2_5S_region.fa; the flank mapping to
  the CHM13 5' flank (<50 kb) marks the 5' terminus, the 3' flank (>327 kb) the 3'.
  order_resolved:
    1 contig                 -> 'single'
    2 edges                  -> 'resolved'
    2 edges + 1 interior     -> 'resolved' (one interior slot; orientation flagged)
    2 edges + 2 interiors    -> 'partial'  (interiors placed in one orientation, flagged)
    anything unanchorable    -> 'unresolved' (kept, flagged; copies still counted)

Usage: 58_multicontig_extract.py SAMPLE HAP ASSEMBLY[.gz]
Imports helper functions and constants from the companion assembly-analysis
module (02_assembly_analysis_external.py).

Paths and the tool-binary directory are read from environment variables:
  FIVES_DATA  input/derived-data directory (blast/, sequences/, databases_mc/, tmp/)
  FIVES_REFS  reference directory (consensus FASTA, CHM13 5S region FASTA)
  FIVES_BIN   directory prepended to PATH for external tools
"""
import os, sys, re, subprocess, tempfile, importlib.util
from pathlib import Path
import numpy as np
import pandas as pd

SCRIPTS = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("ext", SCRIPTS / "02_assembly_analysis_external.py")
ext = importlib.util.module_from_spec(spec); spec.loader.exec_module(ext)

T2T   = Path(os.environ.get("FIVES_DATA", "data"))
HPRC  = T2T
CONS  = Path(os.environ.get("FIVES_REFS", "refs")) / "5S_t2t_consensus.fa"
CHM13_REGION = Path(os.environ.get("FIVES_REFS", "refs")) / "CHM13v2_5S_region.fa"
CHM13_ARRAY_LO, CHM13_ARRAY_HI = 50001, 327609   # 5S array within CHM13v2_5S_region.fa
TMPDIR = HPRC / "tmp"; TMPDIR.mkdir(parents=True, exist_ok=True)
ENVBIN = os.environ.get("FIVES_BIN", "")
os.environ["PATH"] = ENVBIN + ":" + os.environ.get("PATH", "")

QUERY_LEN, NTS_PRE_BP, NTS_POST_BP = ext.QUERY_LEN, ext.NTS_PRE_BP, ext.NTS_POST_BP
BLAST_COLS, BLAST_FMT = ext.BLAST_COLS, ext.BLAST_FMT
MIN_LEN, MIN_PID = ext.BLAST_MIN_LEN, ext.BLAST_MIN_PIDENT
SPACING_LO, SPACING_HI = 1800, 2600
FLANK_PROBE = 60000          # bp of flank to align for orientation
MIN_ARRAY_COPIES = 3


def run(cmd):
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")


def blast_assembly(tmp_fa, blast_out):
    if blast_out.exists() and blast_out.stat().st_size > 50:
        return
    with tempfile.TemporaryDirectory(dir=TMPDIR) as td:
        db = Path(td) / "db"
        run(f"makeblastdb -in '{tmp_fa}' -dbtype nucl -out '{db}' -logfile /dev/null")
        run(f"blastn -query '{CONS}' -db '{db}' -outfmt '{BLAST_FMT}' "
            f"-word_size 7 -evalue 1e-10 -perc_identity 90 -num_threads 6 -out '{blast_out}'")


def contig_geometry(g, slen):
    """g = hits on one contig (already lo/hi). Returns geometry dict."""
    lo, hi = int(g.gene_lo.min()), int(g.gene_hi.max())
    left, right = lo, slen - hi
    sp = np.diff(np.sort(g.gene_lo.values))
    med_sp = int(np.median(sp)) if len(sp) else 0
    return dict(n=len(g), slen=slen, array_lo=lo, array_hi=hi,
                left_flank=left, right_flank=right,
                big_flank=max(left, right), break_gap=min(left, right),
                med_spacing=med_sp, strand=g.sstrand.value_counts().idxmax())


def anchor_flank(tmp_fa, contig, geom):
    """minimap2 the contig's larger flank to CHM13 region; return '5prime'/'3prime'/None."""
    left, right = geom["left_flank"], geom["right_flank"]
    if max(left, right) < 20000:
        return None
    with tempfile.TemporaryDirectory(dir=TMPDIR) as td:
        flank_fa = Path(td) / "flank.fa"
        # probe the unique sequence ADJACENT to the array (not the distal contig end)
        if left >= right:           # flank upstream of the array
            lo = max(1, geom['array_lo'] - FLANK_PROBE); hi = max(1, geom['array_lo'] - 1)
        else:                        # flank downstream of the array
            lo = min(geom['slen'], geom['array_hi'] + 1); hi = min(geom['slen'], geom['array_hi'] + FLANK_PROBE)
        reg = f"{contig}:{lo}-{hi}"
        run(f"samtools faidx '{tmp_fa}' '{reg}' > '{flank_fa}'")
        paf = subprocess.run(
            f"minimap2 -x asm20 -t 4 '{CHM13_REGION}' '{flank_fa}' 2>/dev/null",
            shell=True, capture_output=True, text=True, executable="/bin/bash").stdout
    best = None; best_mlen = 0
    for ln in paf.splitlines():
        f = ln.split("\t")
        if len(f) < 9:
            continue
        tstart, tend, mlen = int(f[7]), int(f[8]), int(f[9])
        if mlen > best_mlen:
            best_mlen = mlen
            mid = (tstart + tend) / 2
            best = "5prime" if mid < CHM13_ARRAY_LO else ("3prime" if mid > CHM13_ARRAY_HI else None)
    # which array side this flank sits on, on the contig, matters for join orientation
    side = "upstream" if left >= right else "downstream"
    return dict(anchor=best, flank_side=side, mlen=best_mlen)


def order_contigs(geoms, anchors):
    """Return ordered list of contig names (5'->3') + order_resolved tag + per-contig flip."""
    names = list(geoms)
    if len(names) == 1:
        return names, "single", {names[0]: False}
    edges = [c for c in names if geoms[c]["big_flank"] >= 20000 and geoms[c]["break_gap"] < geoms[c]["med_spacing"] * 3]
    interiors = [c for c in names if c not in edges]
    five = [c for c in edges if anchors.get(c, {}).get("anchor") == "5prime"]
    three = [c for c in edges if anchors.get(c, {}).get("anchor") == "3prime"]
    resolved = "resolved"
    # need exactly one 5' and one 3' edge for a clean order
    if len(five) == 1 and len(three) == 1:
        c5, c3 = five[0], three[0]
        mids = [c for c in names if c not in (c5, c3)]
        if len(mids) >= 2:
            resolved = "partial"      # >=2 interiors: relative order/orientation arbitrary
        ordered = [c5] + mids + [c3]
    else:
        # ambiguous anchoring (e.g. 4-contig with 2 same-end anchors): best-effort
        # order — 5'-anchored first, 3'-anchored last, interiors/ambiguous between
        # by copy count. 'partial' if any anchor exists (flag for order-sensitive
        # analyses), else 'unresolved'.
        any_anchor = any(anchors.get(c, {}).get("anchor") for c in names)
        def keyf(c):
            a = anchors.get(c, {}).get("anchor")
            return (0 if a == "5prime" else 2 if a == "3prime" else 1, -geoms[c]["n"])
        ordered = sorted(names, key=keyf)
        resolved = "partial" if any_anchor else "unresolved"
    # orientation per contig: make each contig's copies run 5'->3' along the array.
    # 5' edge: array should be at its 3' (downstream) end -> flip if flank is downstream.
    flip = {}
    for rank, c in enumerate(ordered):
        a = anchors.get(c, {})
        if c == ordered[0] and len(ordered) > 1:
            flip[c] = (a.get("flank_side") == "downstream")   # 5' edge: flank upstream
        elif c == ordered[-1] and len(ordered) > 1:
            flip[c] = (a.get("flank_side") == "upstream")     # 3' edge: flank downstream
        else:
            flip[c] = (geoms[c]["strand"] == "minus")         # interior: best-effort by strand
    return ordered, resolved, flip


def extract_oriented(tmp_fa, items, out_fa):
    """items: list of (region, cid, revcomp_bool). Extract each region in the given
    orientation so every copy ends up in CONSENSUS orientation regardless of the
    copy's genomic strand."""
    if not items:
        Path(out_fa).write_text(""); return
    out = {}
    for rev in (True, False):
        regs = [it[0] for it in items if it[2] == rev]
        ids  = [it[1] for it in items if it[2] == rev]
        if not regs: continue
        cmd = ["samtools", "faidx"] + (["-i"] if rev else []) + [str(tmp_fa)] + regs
        r = subprocess.run(cmd, capture_output=True, text=True)
        for (_, seq), cid in zip(ext.read_fasta_ordered(r.stdout), ids):
            out[cid] = seq
    with open(out_fa, "w") as fh:
        for _, cid, _ in items:
            if cid in out: fh.write(f">{cid}\n{out[cid]}\n")


def analyse(sample_id, hap, assembly):
    out_dir = HPRC / "sequences" / sample_id; out_dir.mkdir(parents=True, exist_ok=True)
    blast_out = HPRC / "blast" / f"{sample_id}_{hap}_blast.txt"
    db_out    = HPRC / "databases_mc" / f"{sample_id}_{hap}.tsv"
    db_out.parent.mkdir(parents=True, exist_ok=True)
    gene_fa = out_dir / f"{sample_id}_{hap}_mc_gene.fa"
    pre_fa  = out_dir / f"{sample_id}_{hap}_mc_pre.fa"
    post_fa = out_dir / f"{sample_id}_{hap}_mc_post.fa"
    gene_aln= out_dir / f"{sample_id}_{hap}_mc_gene_aln.fa"
    pre_aln = out_dir / f"{sample_id}_{hap}_mc_pre_aln.fa"
    post_aln= out_dir / f"{sample_id}_{hap}_mc_post_aln.fa"

    with tempfile.TemporaryDirectory(dir=TMPDIR) as td:
        tmp_gz = Path(td) / "asm.fa.gz"; tmp_fa = Path(td) / "asm.fa"
        run(f"cp '{assembly}' '{tmp_gz}'")
        run(f"gunzip -c '{tmp_gz}' > '{tmp_fa}'")
        run(f"samtools faidx '{tmp_fa}'")
        blast_assembly(tmp_fa, blast_out)

        df = pd.read_csv(blast_out, sep="\t", names=BLAST_COLS)
        hits = df[(df.length >= MIN_LEN) & (df.pident >= MIN_PID)].copy()
        hits["gene_lo"] = hits[["sstart", "send"]].min(axis=1).astype(int)
        hits["gene_hi"] = hits[["sstart", "send"]].max(axis=1).astype(int)

        # classify contigs
        geoms, orphan_rows = {}, []
        for ctg, g in hits.groupby("sseqid"):
            slen = int(g.slen.iloc[0]); g = g.sort_values("gene_lo")
            if len(g) >= MIN_ARRAY_COPIES:
                gm = contig_geometry(g, slen)
                if SPACING_LO <= gm["med_spacing"] <= SPACING_HI:
                    geoms[ctg] = (gm, g); continue
            for _, r in g.iterrows():    # orphan / dispersed
                orphan_rows.append((ctg, int(r.gene_lo), int(r.gene_hi), float(r.pident), str(r.sstrand)))
        if not geoms:
            print(f"  {sample_id} {hap}: no array-fragment contig"); return

        # anchor + order
        anchors = {c: anchor_flank(tmp_fa, c, gm) or {} for c, (gm, _) in geoms.items()}
        gm_only = {c: gm for c, (gm, _) in geoms.items()}
        ordered, resolved, flip = order_contigs(gm_only, anchors)
        if os.environ.get("MC_DEBUG"):
            for c in geoms:
                g = gm_only[c]; a = anchors.get(c, {})
                print(f"  [dbg] {c}: n={g['n']} slen={g['slen']} L={g['left_flank']} R={g['right_flank']} "
                      f"bigflank={g['big_flank']} break_gap={g['break_gap']} medsp={g['med_spacing']} "
                      f"-> anchor={a.get('anchor')} side={a.get('flank_side')} mlen={a.get('mlen')}")
            print(f"  [dbg] ordered={ordered} resolved={resolved}")

        # build global copy list 5'->3'
        copies = []   # (contig, contig_rank, gene_lo, gene_hi, pident, mismatch, gapopen, strand)
        for rank, ctg in enumerate(ordered, 1):
            gm, g = geoms[ctg]
            g2 = g.sort_values("gene_lo", ascending=not flip[ctg])
            for _, r in g2.iterrows():
                copies.append(dict(contig=ctg, rank=rank, gene_lo=int(r.gene_lo), gene_hi=int(r.gene_hi),
                                   pident=float(r.pident), mismatch=int(r.mismatch),
                                   gapopen=int(r.gapopen), strand=str(r.sstrand)))
        N = len(copies)

        # extract per-copy regions in CONSENSUS orientation (strand-aware). The 5S
        # array consensus matches the minus strand of most contigs; for minus copies
        # the consensus 5' (NTS-pre) lies DOWNSTREAM of the gene and the seq must be
        # reverse-complemented. For plus copies (minority contigs) it is the mirror:
        # NTS-pre lies UPSTREAM and no revcomp. This is applied per-copy according
        # to each copy's genomic strand.
        fai = {l.split("\t")[0]: int(l.split("\t")[1]) for l in (Path(str(tmp_fa)+".fai")).read_text().splitlines()}
        gene_items, pre_items, post_items = [], [], []
        for i, cp in enumerate(copies, 1):
            cid = f"copy{i:03d}"; clen = fai[cp["contig"]]; ctg = cp["contig"]
            lo, hi = cp["gene_lo"], cp["gene_hi"]; minus = (cp["strand"] == "minus")
            gene_items.append((f"{ctg}:{lo}-{hi}", cid, minus))
            if minus:
                if hi + 1 <= clen: pre_items.append((f"{ctg}:{hi+1}-{min(clen, hi+NTS_PRE_BP)}", cid, True))
                if lo > 1:         post_items.append((f"{ctg}:{max(1, lo-NTS_POST_BP)}-{lo-1}", cid, True))
            else:  # plus-strand copy: consensus 5'/3' map to the opposite genomic sides, no revcomp
                if lo > 1:         pre_items.append((f"{ctg}:{max(1, lo-NTS_PRE_BP)}-{lo-1}", cid, False))
                if hi + 1 <= clen: post_items.append((f"{ctg}:{hi+1}-{min(clen, hi+NTS_POST_BP)}", cid, False))
        extract_oriented(tmp_fa, gene_items, gene_fa)
        extract_oriented(tmp_fa, pre_items,  pre_fa)
        extract_oriented(tmp_fa, post_items, post_fa)

    # MAFFT (all copies of the haplotype together)
    for fa, aln in [(gene_fa, gene_aln), (pre_fa, pre_aln), (post_fa, post_aln)]:
        if fa.stat().st_size < 10:
            Path(aln).write_text(""); continue
        with open(aln, "w") as o:
            subprocess.run(["mafft", "--auto", "--quiet", "--thread", "4", str(fa)],
                           stdout=o, stderr=subprocess.DEVNULL)

    def parse_safe(p):
        return ext.parse_alignment(p) if Path(p).exists() and Path(p).stat().st_size > 10 \
               else ([], np.zeros((0,0)), "", np.zeros((0,0)), [])
    g_cids,g_mat,g_cons,_,g_ng = parse_safe(gene_aln)
    p_cids,p_mat,p_cons,_,p_ng = parse_safe(pre_aln)
    o_cids,o_mat,o_cons,_,o_ng = parse_safe(post_aln)
    gi={c:i for i,c in enumerate(g_cids)}; pi={c:i for i,c in enumerate(p_cids)}; oi={c:i for i,c in enumerate(o_cids)}

    rows = []
    for i, cp in enumerate(copies, 1):
        cid = f"copy{i:03d}"
        gv  = ext.call_variants(g_mat[gi[cid]], g_cons, g_ng) if cid in gi else []
        prv = ext.call_variants(p_mat[pi[cid]], p_cons, p_ng) if cid in pi else []
        pov = ext.call_variants(o_mat[oi[cid]], o_cons, o_ng) if cid in oi else []
        gv_gene = [v for v in gv if NTS_PRE_BP <= int(v.split(':')[0]) < NTS_PRE_BP+QUERY_LEN]
        nxt = copies[i] if i < N else None
        same_contig = nxt and nxt["contig"] == cp["contig"]
        sp = (nxt["gene_lo"] - cp["gene_lo"]) if same_contig else None
        if i == 1: border = "5-prime_array_border"
        elif i == N: border = "3-prime_array_border"
        elif not same_contig: border = "junction_border"
        elif copies[i-2]["contig"] != cp["contig"]: border = "junction_border"
        else: border = "interior"
        rows.append(dict(copy_id=i, sample_id=sample_id, haplotype=hap,
            array_chrom=cp["contig"], source_contig=cp["contig"], contig_rank=cp["rank"],
            order_resolved=resolved, array_member=1, category="array",
            gene_lo_local=cp["gene_lo"], gene_hi_local=cp["gene_hi"],
            unit_start_local=cp["gene_lo"]-1420,
            unit_end_local=(sp+cp["gene_lo"]-1421) if sp else cp["gene_hi"]+NTS_PRE_BP,
            strand=cp["strand"], gene_pct_identity=cp["pident"],
            gene_mismatches=cp["mismatch"], gene_gaps=cp["gapopen"],
            category2="identical" if not gv_gene else "highly_similar",
            unit_length_bp=(sp if sp else None), spacing_to_next_bp=sp,
            n_snv_gene=len(gv), n_snv_5s_gene=len(gv_gene),
            n_snv_nts_pre=len(prv), n_snv_nts_post=len(pov),
            gene_variants="; ".join(gv) or "none",
            nts_pre_variants="; ".join(prv) or "none",
            nts_post_variants="; ".join(pov) or "none",
            border_note=border))
    # orphan rows (flagged, not in array)
    for ctg, lo, hi, pid, strand in orphan_rows:
        rows.append(dict(copy_id=None, sample_id=sample_id, haplotype=hap,
            array_chrom=ctg, source_contig=ctg, contig_rank=None, order_resolved=resolved,
            array_member=0, category="dispersed_5S", gene_lo_local=lo, gene_hi_local=hi,
            unit_start_local=None, unit_end_local=None, strand=strand, gene_pct_identity=pid,
            gene_mismatches=None, gene_gaps=None, category2="dispersed",
            unit_length_bp=None, spacing_to_next_bp=None, n_snv_gene=None, n_snv_5s_gene=None,
            n_snv_nts_pre=None, n_snv_nts_post=None, gene_variants="none",
            nts_pre_variants="none", nts_post_variants="none", border_note="dispersed"))
    out = pd.DataFrame(rows)
    out.to_csv(db_out, sep="\t", index=False)
    nc = (out.array_member == 1).sum(); no = (out.array_member == 0).sum()
    print(f"  {sample_id} {hap}: {len(ordered)} array contig(s) [{resolved}], "
          f"{nc} array copies + {no} dispersed -> {db_out.name}")


if __name__ == "__main__":
    analyse(sys.argv[1], sys.argv[2], sys.argv[3])
