#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 02_ont_45S_flank_stream.py — Genome-competitive ONT read placement to call
# 45S rDNA flank and array-edge methylation per NOR.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
02_ont_45S_flank_stream.py  <sample> [--cells N]

45S rDNA flank methylation by GENOME-COMPETITIVE ONT read placement across the five NORs
(chr13/14/15/21/22), using ONT Dorado 5mCG_5hmCG reads.

Why genome competition does the work: the five rDNA arrays are near-identical, so a read that
sits purely inside an array maps ambiguously across all five NORs -> low MAPQ -> dropped. Only
reads ANCHORED in unique DJ/PJ flank map confidently to a single NOR (MAPQ>=20); their array
portion then yields edge-relative methylation. This isolates exactly the array EDGE copies.

Two passes:
  pass 1  all reads -> NOR-flank ref (arrays masked to N) -> flank-anchored candidate reads
  pass 2  candidates -> whole CHM13 genome -> primary in a NOR region & MAPQ>=20 -> modkit
Output: <FIVES_DATA>/methylation/calls_{sample}.tsv  (read_id, chrom, ref_position, meth)

Paths are read from environment variables (FIVES_REFS, FIVES_DATA;
AWS, SAMTOOLS, MINIMAP2, MODKIT).
"""
import sys, os, argparse, subprocess, shutil
from pathlib import Path
import pandas as pd

os.environ.setdefault("TMPDIR", "/tmp")
Path(os.environ["TMPDIR"]).mkdir(parents=True, exist_ok=True)
H45  = Path(os.environ.get("FIVES_DATA", "data"))
WORK = H45 / "work"
AWS, SAM, MM2, MODKIT = os.environ.get("AWS", "aws"), os.environ.get("SAMTOOLS", "samtools"), os.environ.get("MINIMAP2", "minimap2"), os.environ.get("MODKIT", "modkit")
PASS1_REF = H45 / "reference" / "human_45S_transcribed_NR046235.fa"   # 45S transcribed region only
#   (13.4 kb: 5'ETS-18S-ITS1-5.8S-ITS2-28S-3'ETS; RefSeq NR_046235). rRNA-gene sequence is
#   rDNA-specific and carries no dispersed repeats, so this is a clean, cheap pre-filter.
#   Every genuine rDNA / array-edge / junction read contains a gene -> anchors here;
#   pure-IGS-only reads (rare) are skipped. Unique single-NOR placement of hybrids happens
#   in pass 2 vs the whole CHM13 (MAPQ>=20).
GENOME    = Path(os.environ.get("FIVES_REFS", "refs")) / "chm13.map-ont.mmi"
GENOME_FA = Path(os.environ.get("FIVES_REFS", "refs")) / "chm13v2.0.fa"
CFG  = Path(os.environ.get("FIVES_DATA", "data")) / "aws_config"   # multipart config
S3ROOT = "s3://human-pangenomics/working/HPRC"
NOR = {"chr13": (5_770_548, 9_348_041), "chr14": (2_099_537, 2_817_811),
       "chr15": (2_506_442, 4_707_485), "chr21": (3_108_298, 5_612_715),
       "chr22": (4_793_794, 5_720_650)}
W = 250_000                                              # flank half-window (matches pass1 ref)
MOD_LO, MOD_HI, MAPQ = 0.2, 0.8, 20
AWSENV = f"AWS_CONFIG_FILE='{CFG}'"
MMT, FQT, MODT = 12, 8, 6

def run(cmd): subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
def sh(cmd):  return subprocess.run(cmd, shell=True, capture_output=True, text=True, executable="/bin/bash").stdout

def cleanup(wdir, keep_dl):
    """Remove scratch. With keep_dl, preserve the downloaded cell*_dl.bam(s) for re-runs."""
    if not keep_dl:
        shutil.rmtree(wdir, ignore_errors=True); return
    for p in wdir.glob("*"):
        if "_dl.bam" in p.name: continue
        try: p.unlink()
        except (IsADirectoryError, PermissionError): shutil.rmtree(p, ignore_errors=True)

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("sample"); ap.add_argument("--cells", type=int, default=0)
    ap.add_argument("--keep-dl", action="store_true",
                    help="reuse an already-downloaded cell{i}_dl.bam if present, and do not delete it")
    ap.add_argument("--stream", action="store_true",
                    help="pipe each BAM in-flight from S3 (no local BAM write/read) — use for the cohort batch")
    ap.add_argument("--resume", action="store_true",
                    help="reuse an existing non-empty cand.fq and skip pass 1 (re-run pass 2 only)")
    ap.add_argument("--prefix", default="",
                    help="S3 prefix of the donor's ONT meth BAMs (from inventory ont_meth_prefix); "
                         "REQUIRED for Release-2 donors whose path differs from the batch-1 default")
    a = ap.parse_args(); s = a.sample
    wdir = WORK / s; wdir.mkdir(parents=True, exist_ok=True)
    pref = a.prefix if a.prefix and a.prefix != "nan" else f"{S3ROOT}/{s}/raw_data/nanopore/sup5.0.0_5mCG_5hmCG/"
    if not pref.endswith("/"): pref += "/"
    bams = [l.split()[-1] for l in sh(f"'{AWS}' s3 ls --no-sign-request '{pref}'").splitlines()
            if l.strip().endswith(".bam")]
    if not bams: sys.exit(f"no ONT meth bam for {s}")
    if a.cells > 0: bams = bams[:a.cells]
    print(f"[{s}] {len(bams)} ONT meth bam(s): {bams}", flush=True)
    cand = wdir / "cand.fq"
    if a.resume and cand.exists() and cand.stat().st_size > 0:
        print(f"  [{s}] --resume: reusing cand.fq ({cand.stat().st_size/1e9:.2f} GB), skipping pass 1", flush=True)
        bams = []                                        # skip the pass-1 loop entirely
    else:
        open(cand, "w").close()
    # pass-1 tail: fastq stream -> rDNA-gene pre-filter -> append primary-mapped candidates
    P1TAIL = (f" | '{MM2}' -ax map-ont -y -t {MMT} '{PASS1_REF}' - 2>/dev/null"
              f" | '{SAM}' view -b -F 0x904 2>/dev/null"
              f" | '{SAM}' fastq -@4 -T MM,ML - 2>/dev/null >> '{cand}'")
    for i, bf in enumerate(bams):
        if a.stream:
            # IN-FLIGHT: pipe the BAM straight from S3, never written to local disk;
            # a network-bound single stream that avoids the local BAM write+read.
            run(f"{AWSENV} '{AWS}' s3 cp --no-sign-request '{pref}{bf}' - 2>/dev/null"
                f" | '{SAM}' fastq -@{FQT} -F 0x900 -T MM,ML - 2>/dev/null" + P1TAIL)
        else:
            tmp = wdir / f"cell{i}_dl.bam"
            try:
                if a.keep_dl and tmp.exists() and tmp.stat().st_size > 1_000_000_000:
                    print(f"  [{s}] reusing existing {tmp.name} ({tmp.stat().st_size/1e9:.1f} GB)", flush=True)
                else:
                    run(f"{AWSENV} '{AWS}' s3 cp --no-sign-request '{pref}{bf}' '{tmp}' --quiet")
                run(f"'{SAM}' fastq -@{FQT} -F 0x900 -T MM,ML '{tmp}' 2>/dev/null" + P1TAIL)
            finally:
                if tmp.exists() and not a.keep_dl: os.remove(tmp)
        print(f"  [{s}] cell {i} pass1 done", flush=True)
    ncand = sh(f"grep -c '^@' '{cand}'").strip()
    print(f"  [{s}] {ncand} candidate reads -> genome", flush=True)
    out = H45 / "methylation" / f"calls_{s}.tsv"
    fb = wdir / "pass2.bam"
    # --secondary=no: rDNA UL reads hit all 5 near-identical arrays; computing+emitting those
    # secondaries is the pass-2 bottleneck and we discard them anyway (saved BAM uses -F 0x900,
    # win uses -F 0x904). MAPQ is still derived from the suppressed second-best chain, so the
    # genome-competition / MAPQ>=20 filter is unaffected.
    run(f"'{MM2}' -ax map-ont -y --secondary=no -t {MMT} '{GENOME}' '{cand}' 2>/dev/null"
        f" | '{SAM}' sort -@{FQT} -o '{fb}' 2>/dev/null"); run(f"'{SAM}' index '{fb}'")
    regions = " ".join(f"{ch}:{lo-W}-{hi+W}" for ch, (lo, hi) in NOR.items())
    # SAVE all NOR-region reads (array-interior + flank), ALL MAPQ, primary only, MM/ML
    # preserved -> raw material for later within-unit array methylation work. Array reads
    # map ambiguously across the 5 near-identical NORs so most carry MAPQ 0; we deliberately
    # do NOT apply the MAPQ>=20 floor here (that floor is only for the edge analysis below).
    saved = H45 / "methylation" / f"array_reads_{s}.bam"
    run(f"'{SAM}' view -b -F 0x900 '{fb}' {regions} > '{saved}'"); run(f"'{SAM}' index '{saved}'")
    nsav = sh(f"'{SAM}' view -c '{saved}'").strip()
    print(f"  [{s}] saved {nsav} NOR-region reads (all MAPQ) -> {saved.name}", flush=True)
    # EDGE analysis: flank-anchored reads only -> MAPQ>=20 (confident single-NOR placement)
    win = wdir / "win.bam"
    run(f"'{SAM}' view -b -q {MAPQ} -F 0x904 '{fb}' {regions} > '{win}'"); run(f"'{SAM}' index '{win}'")
    nwin = sh(f"'{SAM}' view -c '{win}'").strip()
    print(f"  [{s}] {nwin} reads primary-in-NOR MAPQ>={MAPQ} (flank-anchored, for edge calls)", flush=True)
    mod = wdir / "modkit.tsv"
    run(f"'{MODKIT}' extract full --cpg --mapped-only --reference '{GENOME_FA}' -t {MODT} '{win}' '{mod}' 2>/dev/null")
    if not mod.exists() or mod.stat().st_size < 100:
        pd.DataFrame(columns=["read_id", "chrom", "ref_position", "meth"]).to_csv(out, sep="\t", index=False)
        print(f"  [{s}] no modkit output"); cleanup(wdir, a.keep_dl); return
    m = pd.read_csv(mod, sep="\t", usecols=["read_id", "ref_position", "mod_qual", "mod_code", "chrom"])
    keep = False
    for ch, (lo, hi) in NOR.items():
        keep = keep | ((m.chrom == ch) & (m.ref_position >= lo - W) & (m.ref_position < hi + W))
    m = m[keep]
    m = m[(m.mod_code == "m") & ((m.mod_qual <= MOD_LO) | (m.mod_qual >= MOD_HI))].copy()
    m["meth"] = (m.mod_qual >= MOD_HI).astype(int)
    m[["read_id", "chrom", "ref_position", "meth"]].to_csv(out, sep="\t", index=False)
    nor_counts = m.groupby("chrom").read_id.nunique().to_dict()
    print(f"[{s}] DONE: {m.read_id.nunique()} flank reads, {len(m)} confident CpG calls; per-NOR reads {nor_counts}", flush=True)
    cleanup(wdir, a.keep_dl)

if __name__ == "__main__":
    main()
