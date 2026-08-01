#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 95_icd10_positional_enrichment.py — positional Manhattan plot of FDR-significant
# ICD10 associations along the 5S rDNA region with a gene-vs-NTS enrichment test.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
95_icd10_positional_enrichment.py

Positional plot of FDR-significant ICD10 associations along the covered 5S rDNA
region, highlighting gene-body variants, plus a Fisher test of gene-variant
enrichment for significant phenotype associations versus the surrounding NTS.

Input : 81_results/per_variant_results.csv  (variant x ICD10-block tests)
Output: HPRC/figures/02_variant_calling_qc/95_icd10_positional/
          icd10_positional_enrichment.pdf / .png   (manhattan, gene highlighted)
          icd10_positional_enrichment_data.tsv      (per-variant min-P / n_sig)
          icd10_sig_hits.tsv                         (the FDR<0.05 associations)
          icd10_enrichment_test.txt                  (gene-vs-NTS Fisher result)

Significance = Benjamini-Hochberg FDR < 0.05 (column `fdr`, ICD10 block level).
Gene body: T2T pos 630-748 (UKBB POS 669-787).  Coordinates plotted in T2T.
"""
import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

T2T   = Path(os.environ.get("FIVES_DATA", "data"))
SRC   = T2T / "81_results/per_variant_results.csv"
OUT   = Path(os.environ.get("FIVES_OUT", "output")) / "02_variant_calling_qc/95_icd10_positional"
OUT.mkdir(parents=True, exist_ok=True)

GENE_LO, GENE_HI = 630, 748          # T2T gene body
WIN_LO,  WIN_HI  = 467, 967          # covered window (T2T)
FDR_CUT          = 0.05
GENE_COLOR       = "#ffe14d"
C_GENE, C_NTS    = "#e67e22", "#9aa0a6"

# ── load: ICD10 block-level, converged tests ────────────────────────────────
df = pd.read_csv(SRC)
d  = df[(df.phenotype_level == "block") & (df.converged)].copy()
d["grp"] = np.where(d.region == "gene", "gene", "nts")
d["sig"] = d.fdr < FDR_CUT
d["nlp"] = -np.log10(d.pval.clip(lower=1e-300))

# p-value corresponding to the FDR cutoff (largest p still passing BH)
p_line = d.loc[d.sig, "pval"].max() if d.sig.any() else np.nan

# ── per-variant summary (for table) ─────────────────────────────────────────
gv = (d.groupby(["variant_id", "t2t_pos", "region"])
        .agg(n_tested=("pval", "size"),
             n_sig=("sig", "sum"),
             min_p=("pval", "min"),
             min_fdr=("fdr", "min")).reset_index())
gv["grp"] = np.where(gv.region == "gene", "gene", "nts")
gv["hit"] = gv.n_sig >= 1
gv.sort_values("t2t_pos").to_csv(OUT / "icd10_positional_enrichment_data.tsv",
                                 sep="\t", index=False)

# ── the FDR-significant associations themselves ─────────────────────────────
hits = (d[d.sig][["variant_id", "t2t_pos", "ukbb_pos", "region", "n_carriers",
                  "phenotype", "phenotype_desc", "n_carrier_cases",
                  "or_", "ci_lo", "ci_hi", "pval", "fdr"]]
        .sort_values(["region", "pval"]))
hits.to_csv(OUT / "icd10_sig_hits.tsv", sep="\t", index=False)

# ── enrichment test: gene vs NTS ────────────────────────────────────────────
# variant-level: fraction of variants with >=1 FDR-sig association
ctv = pd.crosstab(gv.grp, gv.hit).reindex(index=["gene", "nts"],
                                          columns=[True, False], fill_value=0)
OR_v, p_v = stats.fisher_exact(ctv.values)
# test-level: fraction of all tests that are FDR-sig
ctt = pd.crosstab(d.grp, d.sig).reindex(index=["gene", "nts"],
                                        columns=[True, False], fill_value=0)
OR_t, p_t = stats.fisher_exact(ctt.values)

gene_rate = gv[gv.grp == "gene"].hit.mean()
nts_rate  = gv[gv.grp == "nts"].hit.mean()

# one-sided permutation test (primary): shuffle gene/NTS region labels across
# variants, null distribution of the number of gene variants with >=1 FDR-sig
# association. Exact, properly calibrated for the small hit count.
NPERM = 20000
rng    = np.random.default_rng(42)
hitvar = gv.hit.values
n_gene = int((gv.grp == "gene").sum())
obs    = int(gv[gv.grp == "gene"].hit.sum())
idx    = np.arange(len(gv))
null   = np.array([hitvar[rng.choice(idx, n_gene, replace=False)].sum()
                   for _ in range(NPERM)])
p_perm = (np.sum(null >= obs) + 1) / (NPERM + 1)

with open(OUT / "icd10_enrichment_test.txt", "w") as fh:
    fh.write("Gene-variant enrichment for FDR<0.05 ICD10 associations (vs surrounding NTS)\n")
    fh.write(f"Significance: BH-FDR < {FDR_CUT}; ICD10 block level; converged logistic tests\n\n")
    fh.write("Variant-level (variant carries >=1 FDR-sig association):\n")
    fh.write(f"  gene: {int(ctv.loc['gene',True])}/{int(ctv.loc['gene'].sum())} "
             f"({gene_rate*100:.1f}%)   NTS: {int(ctv.loc['nts',True])}/{int(ctv.loc['nts'].sum())} "
             f"({nts_rate*100:.1f}%)\n")
    fh.write(f"  Permutation (1-sided, {NPERM} perms): obs={obs} gene hits, "
             f"null mean={null.mean():.2f}  P={p_perm:.4g}\n")
    fh.write(f"  Fisher exact (2-sided): OR={OR_v:.2f}  P={p_v:.3g}\n\n")
    fh.write("Test-level (per variant x phenotype association):\n")
    fh.write(f"  gene sig {int(ctt.loc['gene',True])}/{int(ctt.loc['gene'].sum())}  "
             f"NTS sig {int(ctt.loc['nts',True])}/{int(ctt.loc['nts'].sum())}\n")
    fh.write(f"  Fisher exact: OR={OR_t:.2f}  P={p_t:.3g}\n")

# ── figure: manhattan along covered region ──────────────────────────────────
CM = 1 / 2.54
fig, ax = plt.subplots(figsize=(9 * CM, 4.2 * CM))

ax.axvspan(GENE_LO, GENE_HI, color=GENE_COLOR, alpha=0.35, zorder=0)

# all (non-sig) tests, faint, colored by region
ns = d[~d.sig]
for grp, col in [("nts", C_NTS), ("gene", C_GENE)]:
    s = ns[ns.grp == grp]
    ax.scatter(s.t2t_pos, s.nlp, s=2.5, c=col, alpha=0.18,
               linewidths=0, zorder=2, rasterized=True)
# FDR-significant tests, emphasised
sg = d[d.sig]
for grp, col in [("nts", "#5f6368"), ("gene", "#d35400")]:
    s = sg[sg.grp == grp]
    ax.scatter(s.t2t_pos, s.nlp, s=7, c=col, edgecolors="none",
               linewidths=0, zorder=5, rasterized=True,
               label=f"{'Gene' if grp=='gene' else 'NTS'} FDR<0.05 ({len(s)} assoc.)")

if np.isfinite(p_line):
    ax.axhline(-np.log10(p_line), color="#c0392b", lw=0.6, ls="--", zorder=3)
    ax.text(WIN_HI, -np.log10(p_line), " FDR 0.05", fontsize=4.5,
            color="#c0392b", va="bottom", ha="right")

ax.set_xlim(WIN_LO - 5, WIN_HI + 5)
ax.set_xlabel("Position on 5S rDNA repeat unit (T2T)", fontsize=6)
ax.set_ylabel("ICD10 association\n$-$log$_{10}$ P", fontsize=6)
ax.tick_params(labelsize=5.5)
ax.spines[["top", "right"]].set_visible(False)
ax.legend(fontsize=4.6, loc="upper left", handletextpad=0.2,
          borderpad=0.3, labelspacing=0.2, framealpha=0.85)

# enrichment annotation
ax.text(0.985, 0.97,
        f"Gene vs NTS enrichment\n{obs} gene hits vs {null.mean():.1f} exp., "
        f"OR={OR_v:.1f}\nP={p_perm:.3f} (perm), {p_v:.3f} (Fisher)",
        transform=ax.transAxes, fontsize=4.6, va="top", ha="right",
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#cccccc", lw=0.4))

fig.tight_layout(pad=0.3)
fig.savefig(OUT / "icd10_positional_enrichment.pdf", dpi=400, bbox_inches="tight")
fig.savefig(OUT / "icd10_positional_enrichment.png", dpi=300,
            bbox_inches="tight", facecolor="white")

print(f"FDR<{FDR_CUT}: {int(d.sig.sum())} associations / {gv.hit.sum()} variants")
print(f"gene vs NTS: {obs} gene hits vs {null.mean():.2f} expected | "
      f"permutation P={p_perm:.4g} | Fisher OR={OR_v:.2f} P={p_v:.3g}")
print(f"Saved figure + tables -> {OUT}")
