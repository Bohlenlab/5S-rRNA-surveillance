#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# folding_energy.py — Per-variant ViennaRNA RNAfold folding-energy changes (ddG)
# for all single-nucleotide variants of the 5S rRNA sense (transcript) strand.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
5S rRNA folding energies on the SENSE (transcript) strand.

For every single-nucleotide variant at every position of the 119-nt 5S rRNA
transcript, computes the whole-molecule minimum free energy (MFE) via ViennaRNA
RNAfold (default parameters, 37 C) and the folding-energy change
    ddG = MFE(mutant) - MFE(WT)   [kcal/mol]
where positive ddG is destabilizing and negative ddG is stabilizing. Positions
are numbered 5'->3' along the transcript. Sense-strand outputs carry a
_CORRECTED suffix.

Outputs (under FIVES_OUT):
  - 5S-rRNA_sense_CORRECTED.fa                    WT sense-strand FASTA
  - RNAfold_sense_per_variant_CORRECTED.csv       per-variant MFE/ddG table
  - Folding_energies_sense_matrix_CORRECTED.csv   position x base ddG matrix
"""
import RNA, pandas as pd, os

OUT = os.environ.get("FIVES_OUT", "output")
# 5S rRNA sense (transcript) strand, 119 nt
SENSE = "GTCTACGGCCATACCACCCTGAACGCGCCCGATCTCGTCTGATCTCGGAAGCTAAGCAGGGTCGGGCCTGGTTAGTACTTGGATGGGAGACCGCCTGGGAATACCGGGTGCTGTAGGCT"

wt_struct, wt_mfe = RNA.fold(SENSE.replace("T", "U"))
print(f"SENSE WT: len={len(SENSE)}  MFE={wt_mfe:.2f} kcal/mol")
print(f"WT structure: {wt_struct}")

rows = []
for i, ref in enumerate(SENSE):
    pos = i + 1
    for alt in "ACGT":
        if alt == ref:
            continue
        mut = SENSE[:i] + alt + SENSE[i+1:]
        _, mfe = RNA.fold(mut.replace("T", "U"))
        rows.append(dict(Transcript_ID=f"5S_sense-{pos}-{ref}-{alt}", Pos=pos, Ref=ref,
                         Alt=alt, MFE_wt=round(wt_mfe, 2), MFE_mut=round(mfe, 2),
                         ddG=round(mfe - wt_mfe, 2)))
d = pd.DataFrame(rows)

# WT sense-strand FASTA
with open(f"{OUT}/5S-rRNA_sense_CORRECTED.fa", "w") as f:
    f.write(f">5S_rRNA sense strand (119 nt)\n{SENSE}\n")

# per-variant MFE/ddG table (sense strand)
d.to_csv(f"{OUT}/RNAfold_sense_per_variant_CORRECTED.csv", index=False)

# position x base ddG matrix
mat = d.pivot(index="Pos", columns="Alt", values="ddG")
mat.to_csv(f"{OUT}/Folding_energies_sense_matrix_CORRECTED.csv")

print(f"\nwrote {len(d)} variants across {d.Pos.nunique()} positions")
print(f"ddG range: {d.ddG.min():.1f} to {d.ddG.max():.1f}  (median {d.ddG.median():.2f})")
print(d.head(6).to_string(index=False))
