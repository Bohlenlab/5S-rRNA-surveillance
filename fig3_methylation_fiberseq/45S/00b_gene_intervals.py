#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 00b_gene_intervals.py — Locate the transcribed 45S region in each CHM13 rDNA
# copy and write its genomic intervals; input to the 45S edge-methylation plot.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
Align the 45S transcribed unit (NR_046235) to each CHM13 NOR region and write the
matched genomic intervals to methylation/gene_intervals.tsv (chrom, start, end).
"""
import subprocess, os, tempfile
from pathlib import Path

BASE = Path(os.environ.get("FIVES_DATA", "data"))
REFFA = str(Path(os.environ.get("FIVES_REFS", "refs")) / "chm13v2.0.fa")
GENE = BASE / "reference/human_45S_transcribed_NR046235.fa"
SAMTOOLS = os.environ.get("SAMTOOLS", "samtools")
MINIMAP2 = os.environ.get("MINIMAP2", "minimap2")
NOR = {"chr13": (5_770_548, 9_348_041), "chr14": (2_099_537, 2_817_811), "chr15": (2_506_442, 4_707_485),
       "chr21": (3_108_298, 5_612_715), "chr22": (4_793_794, 5_720_650)}
TMP = Path(tempfile.gettempdir())
regfa = TMP / "reg.fa"
out = open(BASE / "methylation/gene_intervals.tsv", "w"); out.write("chrom\tstart\tend\n"); ncopies = 0
for ch, (lo, hi) in NOR.items():
    subprocess.run(f"{SAMTOOLS} faidx {REFFA} {ch}:{lo-260000}-{hi+260000} > {regfa}", shell=True)
    paf = subprocess.run(f"{MINIMAP2} -cx asm20 -N 400 -p 0.02 {regfa} {GENE} 2>/dev/null",
                         shell=True, capture_output=True, text=True).stdout
    off = lo - 260000; iv = []
    for L in paf.splitlines():
        f = L.split("\t")
        if len(f) >= 11 and int(f[10]) >= 8000: iv.append((off + int(f[7]), off + int(f[8])))
    iv.sort(); m = []
    for a, b in iv:
        if m and a <= m[-1][1] + 200: m[-1] = (m[-1][0], max(m[-1][1], b))
        else: m.append((a, b))
    for a, b in m: out.write(f"{ch}\t{a}\t{b}\n")
    ncopies += len(m)
out.close(); print(f"wrote methylation/gene_intervals.tsv ({ncopies} transcribed copies)")
