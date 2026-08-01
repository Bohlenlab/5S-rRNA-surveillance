# -----------------------------------------------------------------------------
# de_annotate.py — build per-donor 5S-variant group labels (carrier / expresser /
# tertile / continuous dose) used by the de_v2 DE contrasts.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""Build per-donor 5S-variant group labels for every de_v2 DE contrast, from the RNA-only pileups.

Source (C.RNAVAF): results/eqtl/extreme/donor_variant_rnavaf.tsv
  donor x variant table; wgs_vaf = DNA VAF, rna_vaf = editing/FP-corrected RNA VAF.
Variant universe (C.GEN, from de_common.py): the gene-region functional 5S variants, where the
  editing/false-positive rule has already dropped any allele called in >0.5% of non-carriers (NONC=0.5%).
A donor is a GEN carrier iff it appears with ANY GEN variant; NC (non-carrier) = no GEN variant.
The donor universe is every donor in the counts matrix (C.H5AD obs.donor), so NC donors are included.

Output: de_v2/groups_donor.tsv, one row per donor with several label columns, each selected by a
different contrast downstream (de_pertissue.py picks the column for the requested contrast):
  g_3group   NC / silent (carrier, RNA<cut) / expr (carrier, RNA>=cut)
  g_dna      carrier vs NC             (pure genotype/carriage effect, ignores RNA)
  g_rna      expr vs nonexpr           (anyone reaching the RNA expresser threshold)
  cont_rna, cont_rna_log   continuous per-donor max RNA-VAF (NC=0; log10 for the design)
  g_tertile  NC / none / low / high    (RNA-expresser tertiles among carriers)

No CLI args; run directly. Paths/env come from de_common.py. Emits group counts to stdout.
"""
import numpy as np, pandas as pd
import de_common as C

EXPR_CUT = 0.007  # RNA-VAF at/above which a carrier counts as a high-expresser (recalibrated, GSEA-robust)
DETECT = 0.003    # RNA-VAF detection floor; carriers below this are the tertile 'none' bin

d = pd.read_csv(C.RNAVAF, sep="\t")
g = d[d.variant.isin(C.GEN)].copy()                       # keep only GEN (gene-region functional) variant rows
# one row per carrier donor: strongest DNA and RNA VAF across their carried GEN variants, and how many
per = g.groupby("donor").agg(max_wgs=("wgs_vaf", "max"), max_rna=("rna_vaf", "max"),
                             n_gen=("variant", "nunique")).reset_index()

# all donors present in the counts matrix (defines the NC universe: donors absent from `per` are non-carriers)
A_obs_donors = pd.read_hdf if False else None                          # (unused; no-op kept for parity)
import anndata as ad
donors = pd.Index(sorted(set(ad.read_h5ad(C.H5AD, backed="r").obs["donor"].astype(str))))

lab = pd.DataFrame({"donor": donors})
lab = lab.merge(per, on="donor", how="left")              # left join: non-carriers get NaN aggregates
lab["is_carrier"] = lab["n_gen"].notna().astype(int)      # carrier iff the donor had >=1 GEN variant
lab["max_rna"] = lab["max_rna"].fillna(0.0)               # non-carriers -> 0 RNA-VAF
lab["max_wgs"] = lab["max_wgs"].fillna(0.0)               # non-carriers -> 0 DNA-VAF

# --- contrast columns (each is a candidate `group` label for the per-tissue DESeq2 model) ---
# 3-group: NC / silent (carrier that does NOT reach the expresser cutoff) / expr (carrier at/above it)
def three(r):
    if not r.is_carrier: return "NC"
    return "expr" if r.max_rna >= EXPR_CUT else "silent"
lab["g_3group"] = lab.apply(three, axis=1)

# DNA-only: carrier vs NC -- pure carriage/genotype effect, independent of whether the variant is expressed
lab["g_dna"] = np.where(lab.is_carrier == 1, "carrier", "NC")

# RNA-only: expr vs nonexpr -- anyone reaching the RNA expresser threshold, regardless of carrier status
lab["g_rna"] = np.where(lab.max_rna >= EXPR_CUT, "expr", "nonexpr")

# continuous RNA dose: carriers carry their max RNA-VAF, NC=0. log10(x + 1e-4) pseudocount for the design.
lab["cont_rna"] = lab["max_rna"]
lab["cont_rna_log"] = np.log10(lab["max_rna"] + 1e-4)

# RNA-expresser tertiles among CARRIERS with detectable RNA: split the detected (>=DETECT) group at its
# median into low/high; carriers below DETECT are 'none'; non-carriers are 'NC'. (median of an empty set -> inf)
det = lab[(lab.is_carrier == 1) & (lab.max_rna >= DETECT)]
med = det.max_rna.median() if len(det) else np.inf
def tert(r):
    if not r.is_carrier: return "NC"
    if r.max_rna < DETECT: return "none"
    return "high" if r.max_rna >= med else "low"
lab["g_tertile"] = lab.apply(tert, axis=1)

lab.to_csv(f"{C.ROOT}/de_v2/groups_donor.tsv", sep="\t", index=False)

print(f"donors total={len(lab)}  GEN carriers={lab.is_carrier.sum()}  NC={(lab.is_carrier==0).sum()}")
print("\n3-group :", lab.g_3group.value_counts().to_dict())
print("DNA     :", lab.g_dna.value_counts().to_dict())
print(f"RNA(>= {EXPR_CUT}):", lab.g_rna.value_counts().to_dict())
print("tertile :", lab.g_tertile.value_counts().to_dict(), f"(detected-median max_rna={med:.4f})")
print(f"\nexpresser donors (rna>= {EXPR_CUT}): {(lab.max_rna>=EXPR_CUT).sum()}  "
      f"| detected (>= {DETECT}): {(lab.max_rna>=DETECT).sum()}")
print("wrote de_v2/groups_donor.tsv")
