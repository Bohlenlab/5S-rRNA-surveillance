#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 00_build_pass1_ref.py — Builds the candidate-capture reference (rDNA unit plus
# satellite-masked NOR flanks) for the ONT 45S streamer.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
00_build_pass1_ref.py

Build the pass-1 candidate-capture reference for the ONT 45S streamer:
  - the rDNA unit (U13369)  -> captures array-interior + junction-spanning reads
  - unique flank on each side of each NOR array, with satellite blocks masked to N
    (censat annotation only; SINE/LINE left intact) to provide flank anchors for
    edge/flank coverage without satellite over-inclusion.

Output: <FIVES_DATA>/reference/pass1_ref.fa
Prints per-window unmasked (usable anchor) bp per NOR.

Paths are read from environment variables (FIVES_REFS, FIVES_DATA, SAMTOOLS).
"""
import os
import subprocess
from pathlib import Path

REF = Path(os.environ.get("FIVES_REFS", "refs"))
H45 = Path(os.environ.get("FIVES_DATA", "data"))
SAM = os.environ.get("SAMTOOLS", "samtools")
FA = REF / "chm13v2.0.fa"
CENSAT = REF / "chm13v2.0_censat_v2.1.bed"
UNIT = H45 / "reference" / "human_rDNA_U13369.fa"
NOR = {"chr13": (5_770_548, 9_348_041), "chr14": (2_099_537, 2_817_811),
       "chr15": (2_506_442, 4_707_485), "chr21": (3_108_298, 5_612_715),
       "chr22": (4_793_794, 5_720_650)}
W = 50_000

# load censat intervals per chrom (BED 0-based). Mask only true tandem satellites
# (hsat/bsat/gsat/alpha-HOR/mon/ACRO_composite/CER/SATR/...); keep 'ct_' (centromeric-
# transition / junction sequence = the DJ/PJ anchor) unmasked.
sat = {}
for L in CENSAT.read_text().splitlines():
    f = L.split("\t")
    if len(f) < 4: continue
    if f[3].startswith("ct_"): continue          # keep transition/junction sequence
    sat.setdefault(f[0], []).append((int(f[1]), int(f[2])))

def faidx(region):
    return "".join(subprocess.run(f"{SAM} faidx {FA} {region}", shell=True,
                                  capture_output=True, text=True).stdout.splitlines()[1:])

out = open(H45 / "reference" / "pass1_ref.fa", "w")
for ch, (lo, hi) in NOR.items():
    for side, a, b in [("DJ", max(1, lo - W), lo), ("PJ", hi, hi + W)]:
        seq = list(faidx(f"{ch}:{a+1}-{b}"))          # window, 1-based inclusive
        # mask any base inside a censat satellite interval
        masked = 0
        for s0, e0 in sat.get(ch, []):
            if e0 <= a or s0 >= b: continue
            for i in range(max(s0, a) - a, min(e0, b) - a):
                if 0 <= i < len(seq) and seq[i] != "N":
                    seq[i] = "N"; masked += 1
        s = "".join(seq); usable = len(s) - s.count("N")
        out.write(f">{ch}_{side}_flank {ch}:{a+1}-{b} usable_bp:{usable}\n")
        for i in range(0, len(s), 60): out.write(s[i:i+60] + "\n")
        print(f"{ch} {side}: window {len(s)} bp, satellite-masked {masked}, USABLE anchor {usable} bp")
# append the rDNA unit verbatim
out.write(UNIT.read_text() if UNIT.read_text().startswith(">") else "")
out.close()
subprocess.run(f"{SAM} faidx {H45}/reference/pass1_ref.fa", shell=True)
print("wrote reference/pass1_ref.fa")
