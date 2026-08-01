#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 01_annotate_45S_assemblies.py — Annotates 45S rDNA array locations and unique
# flanks in HPRC diploid assemblies.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
01_annotate_45S_assemblies.py

Annotate where the 45S rDNA arrays sit in the HPRC year-1 diploid assemblies —
contigs, coordinates, and (since the arrays themselves are not fully resolved) the UNIQUE
FLANK assembled adjacent to each rDNA stretch, plus which CHM13 NOR (chr13/14/15/21/22)
each flank belongs to.

Per (sample, hap):
  1. stream the genbank assembly from S3 to a temp file (no-sign-request)
  2. minimap2 -cx asm20 the rDNA unit (U13369.1) onto the assembly  -> rDNA-bearing contigs
     + the target span covered by rDNA hits on each contig
  3. for each rDNA contig, take the flank sequence on each side of the rDNA span and
     minimap2 it to CHM13v2.0 -> assign NOR + DJ/PJ side, and measure how much unique
     (non-rDNA) flank is actually assembled.

Outputs (append-as-you-go):
  <FIVES_DATA>/annotation/rDNA_contigs_<sample>.tsv     per rDNA-bearing contig
  <FIVES_DATA>/annotation/rDNA_assembly_inventory.tsv   merged, all samples

Paths are read from environment variables (FIVES_REFS, FIVES_DATA;
MINIMAP2, SAMTOOLS, AWS).

Usage:
  python3 01_annotate_45S_assemblies.py --samples ALL --jobs 6
  python3 01_annotate_45S_assemblies.py --samples HG00438,HG00735        # pilot
"""
import sys, os, argparse, tempfile, subprocess, shutil, traceback
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd

H45  = Path(os.environ.get("FIVES_DATA", "data"))
ANN  = H45 / "annotation"
REF  = H45 / "reference"
MINIMAP2, SAMTOOLS, AWS = os.environ.get("MINIMAP2", "minimap2"), os.environ.get("SAMTOOLS", "samtools"), os.environ.get("AWS", "aws")
RDNA_UNIT = REF / "human_rDNA_U13369.fa"
CHM13     = Path(os.environ.get("FIVES_REFS", "refs")) / "chm13.asm20.mmi"   # prebuilt asm20 index (fast reuse)
INV       = H45 / "data_inventory.tsv"

# CHM13 NOR (rDNA array) spans from chm13v2.0_censat_v2.1.bed
NOR = {"chr13": (5_770_548, 9_348_041), "chr14": (2_099_537, 2_817_811),
       "chr15": (2_506_442, 4_707_485), "chr21": (3_108_298, 5_612_715),
       "chr22": (4_793_794, 5_720_650)}
MIN_FLANK = 5000          # report flank only if >= this much unique sequence assembled
THREADS   = 4

def sh(cmd, **kw):
    return subprocess.run(cmd, shell=True, executable="/bin/bash", **kw)

def paf(cmd):
    return sh(cmd, capture_output=True, text=True).stdout

def classify(chrom, pos):
    """Classify a CHM13 alignment midpoint: (nor, zone) where zone is
    'array' (inside an acrocentric rDNA span), 'DJ' (distal/telomere-side unique flank,
    lower coord than the array), 'PJ' (proximal/centromere-side, higher coord), or
    ('','unique') for any other (non-acrocentric or far) unique sequence."""
    if chrom not in NOR:
        return ("", "unique")
    lo, hi = NOR[chrom]
    if pos < lo:  return (chrom, "DJ")
    if pos > hi:  return (chrom, "PJ")
    return (chrom, "array")

def process(sample, hap, s3):
    rows = []
    wdir = Path(tempfile.mkdtemp(prefix=f"{sample}_{hap}_", dir="/tmp"))
    asm = wdir / "asm.fa"
    try:
        if not s3 or s3 == "nan":
            return rows, f"SKIP {sample} {hap}: no S3 path"
        sh(f"'{AWS}' s3 cp --no-sign-request '{s3}' - 2>/dev/null | gzip -dc > '{asm}'")
        if not asm.exists() or asm.stat().st_size < 1_000_000:
            return rows, f"FAIL {sample} {hap}: download/gunzip failed"
        sh(f"'{SAMTOOLS}' faidx '{asm}'")
        clen = {l.split('\t')[0]: int(l.split('\t')[1])
                for l in (asm.with_suffix('.fa.fai')).read_text().splitlines()}
        # 1. identify rDNA-bearing contigs: align the rDNA unit (keep secondaries so
        #    a collapsed multi-unit array is still flagged) and total the matched bp.
        p = paf(f"'{MINIMAP2}' -cx asm20 -t {THREADS} -p 0.1 -N 50 '{asm}' '{RDNA_UNIT}' 2>/dev/null")
        unit_bp = {}
        for L in p.splitlines():
            f = L.split("\t")
            if len(f) < 11: continue
            tgt, ml = f[5], int(f[10])
            if ml < 1000: continue
            unit_bp[tgt] = unit_bp.get(tgt, 0) + ml
        if not unit_bp:
            return rows, f"OK {sample} {hap}: 0 rDNA contigs"
        # 2. per rDNA contig: map the WHOLE contig to CHM13 and partition query bp into
        #    array vs unique DJ/PJ flank; record the array<->flank boundary (anchor point).
        for tgt in sorted(unit_bp, key=unit_bp.get, reverse=True):
            L = clen[tgt]
            ctg = wdir / "ctg.fa"
            sh(f"'{SAMTOOLS}' faidx '{asm}' '{tgt}' > '{ctg}'")
            cp = paf(f"'{MINIMAP2}' -cx asm20 -t {THREADS} --secondary=no '{CHM13}' '{ctg}' 2>/dev/null")
            zb = {"array": 0, "DJ": 0, "PJ": 0, "unique": 0}
            nor_bp = {}                       # nor -> array bp (dominant NOR)
            arr_qs, arr_qe = [], []           # contig coords of array-aligned blocks
            flank_blocks = {"DJ": [], "PJ": []}
            for Lh in cp.splitlines():
                g = Lh.split("\t")
                if len(g) < 11: continue
                qs, qe, ch, ts, te, ml = int(g[2]), int(g[3]), g[5], int(g[7]), int(g[8]), int(g[10])
                if ml < 2000: continue
                nor, zone = classify(ch, (ts + te) // 2)
                qbp = qe - qs
                if zone == "array":
                    zb["array"] += qbp; arr_qs.append(qs); arr_qe.append(qe)
                    if nor: nor_bp[nor] = nor_bp.get(nor, 0) + qbp
                elif zone in ("DJ", "PJ"):
                    zb[zone] += qbp; flank_blocks[zone].append((qs, qe, nor))
                else:
                    zb["unique"] += qbp
            dom_nor = max(nor_bp, key=nor_bp.get) if nor_bp else ""
            arr_lo = min(arr_qs) if arr_qs else -1
            arr_hi = max(arr_qe) if arr_qe else -1
            # unique flank adjacent to the array, by side, restricted to the dominant NOR
            def side_flank(side):
                bl = [(qs, qe) for qs, qe, n in flank_blocks[side] if not dom_nor or n == dom_nor]
                return sum(qe - qs for qs, qe in bl), bl
            dj_bp, dj_bl = side_flank("DJ")
            pj_bp, pj_bl = side_flank("PJ")
            rows.append((sample, hap, tgt, L, unit_bp[tgt], zb["array"], dom_nor,
                         arr_lo, arr_hi, dj_bp, pj_bp, zb["unique"]))
        return rows, (f"OK {sample} {hap}: {len(unit_bp)} rDNA contig(s), "
                      f"top NOR flank DJ/PJ recorded")
    except Exception as e:
        return rows, f"FAIL {sample} {hap}: {e}\n{traceback.format_exc()[:300]}"
    finally:
        shutil.rmtree(wdir, ignore_errors=True)

HDR = ("sample\thap\tcontig\tcontig_len\tunit_match_bp\tarray_bp\tnor\tarray_qstart\t"
       "array_qend\tDJ_flank_bp\tPJ_flank_bp\tother_unique_bp\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", default="ALL")
    ap.add_argument("--jobs", type=int, default=6)
    a = ap.parse_args()
    ANN.mkdir(parents=True, exist_ok=True)
    inv = pd.read_csv(INV, sep="\t", dtype=str)
    inv = inv[inv.has_ont_meth.astype(str).isin({"True", "true", "1"})]
    if a.samples != "ALL":
        want = set(a.samples.split(","))
        inv = inv[inv.sample_id.isin(want)]
    jobs = []
    for _, r in inv.iterrows():
        jobs.append((r.sample_id, "hap1", r.assembly_hap1_s3))
        jobs.append((r.sample_id, "hap2", r.assembly_hap2_s3))
    merged = ANN / "rDNA_assembly_inventory.tsv"
    if not merged.exists():
        merged.write_text(HDR)
    print(f"{len(jobs)} (sample,hap) jobs, {a.jobs} concurrent", flush=True)
    done = 0
    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        futs = {ex.submit(process, s, h, s3): (s, h) for s, h, s3 in jobs}
        for fut in as_completed(futs):
            rows, msg = fut.result()
            if rows:
                with open(merged, "a") as f:
                    for row in rows:
                        f.write("\t".join(map(str, row)) + "\n")
            done += 1
            print(f"  [{done}/{len(jobs)}] {msg}", flush=True)
    print(f"DONE -> {merged}", flush=True)

if __name__ == "__main__":
    main()
