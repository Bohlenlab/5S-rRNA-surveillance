#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# folding_validate.py — Cross-algorithm validation of sense-strand ddG by
# re-folding every SNV with RNAstructure and correlating against RNAfold ddG.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""Independent cross-algorithm validation of the sense-strand ΔΔG:
re-fold WT + every SNV with RNAstructure (Mathews-lab Turner-2004 implementation),
compute ddG = MFE(mut) - MFE(WT), and correlate with the ViennaRNA RNAfold ddG.
High correlation => the folding energies are algorithm-independent (validated)."""
import subprocess, os, re, tempfile, pandas as pd
from scipy.stats import spearmanr, pearsonr
OUT = os.environ.get("FIVES_OUT", "output")
os.environ["DATAPATH"] = os.environ.get("RNASTRUCTURE_DATAPATH", os.environ.get("DATAPATH", ""))
FOLD = os.environ.get("RNASTRUCTURE_FOLD", "Fold")
SENSE = "GTCTACGGCCATACCACCCTGAACGCGCCCGATCTCGTCTGATCTCGGAAGCTAAGCAGGGTCGGGCCTGGTTAGTACTTGGATGGGAGACCGCCTGGGAATACCGGGTGCTGTAGGCT"
td = tempfile.mkdtemp()
def rs_mfe(seq):
    fa, ct = f"{td}/s.fa", f"{td}/s.ct"
    open(fa, "w").write(">s\n" + seq + "\n")
    subprocess.run(["arch", "-x86_64", FOLD, "-mfe", fa, ct], capture_output=True)
    m = re.search(r"ENERGY = (-?\d+\.?\d*)", open(ct).readline())
    return float(m.group(1)) if m else float("nan")

wt = rs_mfe(SENSE); print(f"RNAstructure WT MFE = {wt:.2f}")
d = pd.read_csv(f"{OUT}/RNAfold_sense_per_variant.csv")
rs = []
for k, (p, a) in enumerate(zip(d.Pos, d.Alt)):
    mut = SENSE[:p-1] + a + SENSE[p:]
    rs.append(round(rs_mfe(mut) - wt, 2))
    if k % 50 == 0: print(f"  {k}/{len(d)}", flush=True)
d["ddG_RNAstructure"] = rs
d.to_csv(f"{OUT}/RNAfold_sense_per_variant.csv", index=False)
sr = spearmanr(d.ddG, d.ddG_RNAstructure); pr = pearsonr(d.ddG, d.ddG_RNAstructure)
print(f"\n[VALIDATION] RNAfold ΔΔG vs RNAstructure ΔΔG (independent algorithm):")
print(f"   Spearman ρ={sr[0]:.3f} (p={sr[1]:.1e})   Pearson r={pr[0]:.3f}   n={len(d)}")
