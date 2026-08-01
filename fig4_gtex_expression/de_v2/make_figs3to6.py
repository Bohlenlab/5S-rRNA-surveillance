# -----------------------------------------------------------------------------
# make_figs3to6.py — summary panels of the all-ancestry bulk analyses: GSEA bar,
# pathway-NES heatmap, and per-tissue RP-regulon and DNA-RNA concordance bars.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""Figures 3-6: summary panels of the all-ancestry bulk analyses.
Inputs: GSEA key table and per-tissue summary table.
"""
import os, numpy as np, pandas as pd
from pathlib import Path
import de_figstyle as S
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
S.setup()
FIG = str(Path(os.environ.get("FIVES_OUT", "output")))
GREEN, RED, BLUE, ORANGE, GREY = "#2ca02c", "#d62728", "#1f77b4", "#ff7f0e", "#bdbdbd"

g = pd.read_csv(f"{Path(os.environ.get('FIVES_DATA','data'))}/gsea_key.tsv", sep="\t")
for c in ["NES", "FDR q-val", "NOM p-val"]:
    g[c] = pd.to_numeric(g[c], errors="coerce")
pt = pd.read_csv(f"{Path(os.environ.get('FIVES_DATA','data'))}/pertissue.tsv", sep="\t")
CULT = {"Cells - Cultured fibroblasts", "Cells - EBV-transformed lymphocytes"}

# ===== Figure 3: expresser-axis GSEA =====
sub = g[g.contrast == "expr_vs_silent"].copy()
# top pathways by |NES| among FDR<0.25, plus force-in the custom RP sets
top = sub[(sub["FDR q-val"] < 0.25)].reindex(sub["NES"].abs().sort_values(ascending=False).index)
keep = pd.concat([top.head(14), sub[sub.Term.isin(["Cyto_Ribosomal_Proteins", "Mito_Ribosomal_Proteins"])]]).drop_duplicates("Term")
keep = keep.sort_values("NES")
fig, ax = plt.subplots(figsize=(S.PANEL * 1.6, S.PANEL * 1.4))
cols = [RED if n < 0 else BLUE for n in keep.NES]
ax.barh(range(len(keep)), keep.NES, color=cols, edgecolor="k", lw=0.5)
ax.set_yticks(range(len(keep))); ax.set_yticklabels([t[:38] for t in keep.Term], fontsize=7)
ax.axvline(0, color="k", lw=1.0)
ax.set_xlabel("NES (expresser vs silent)")
ax.set_title("Expresser-vs-silent GSEA\n(blue = up in expressers, red = down)")
S.save(fig, f"{FIG}/Figure3_expresser_GSEA.pdf")

# ===== Figure 4: NES heatmap across contrasts =====
PATHS = ["Cyto_Ribosomal_Proteins", "Mito_Ribosomal_Proteins",
         "TNF-alpha Signaling via NF-kB", "Interferon Gamma Response", "Inflammatory Response",
         "Oxidative Phosphorylation", "MYC Targets V1", "mTORC1 Signaling", "E2F Targets",
         "p53 Pathway", "Epithelial Mesenchymal Transition", "Hypoxia"]
CONS = ["expr_vs_silent", "rna_expr_vs_not", "RNA_dose_AGG", "dose_high_vs_low", "DNA_carrier"]
M = pd.DataFrame(index=PATHS, columns=CONS, dtype=float)
for p in PATHS:
    for c in CONS:
        r = g[(g.Term == p) & (g.contrast == c)]
        if len(r): M.loc[p, c] = r.NES.iloc[0]
fig, ax = plt.subplots(figsize=(S.PANEL * 1.5, S.PANEL * 1.6))
vmax = np.nanmax(np.abs(M.values.astype(float)))
im = ax.imshow(M.values.astype(float), cmap="RdBu_r", norm=TwoSlopeNorm(0, -vmax, vmax), aspect="auto")
ax.set_xticks(range(len(CONS))); ax.set_xticklabels(CONS, rotation=40, ha="right", fontsize=7)
ax.set_yticks(range(len(PATHS))); ax.set_yticklabels(PATHS, fontsize=7)
for i in range(len(PATHS)):
    for j in range(len(CONS)):
        v = M.values[i, j]
        if not np.isnan(v): ax.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=6,
                                    color="white" if abs(v) > vmax * 0.6 else "k")
fig.colorbar(im, ax=ax, shrink=0.5, label="NES")
ax.set_title("Pathway NES across contrasts")
S.save(fig, f"{FIG}/Figure4_GSEA_heatmap.pdf")

# ===== Figure 5: per-tissue expresser RP-regulon effect =====
d = pt.dropna(subset=["RP_mean_lfc"]).sort_values("RP_mean_lfc")
fig, ax = plt.subplots(figsize=(S.PANEL * 1.4, S.PANEL * 2))
cols = [ORANGE if t in CULT else (RED if v < 0 else BLUE) for t, v in zip(d.tissue, d.RP_mean_lfc)]
ax.barh(range(len(d)), d.RP_mean_lfc * 100, color=cols, edgecolor="k", lw=0.3)
ax.set_yticks(range(len(d))); ax.set_yticklabels([t[:34] for t in d.tissue], fontsize=6)
ax.axvline(0, color="k", lw=1.0)
ax.set_xlabel("mean cyto-RP log2FC ×100  (expresser vs silent)")
ax.set_title("Per-tissue ribosomal-protein regulon shift\n(orange = cultured cells)")
S.save(fig, f"{FIG}/Figure5_pertissue_RP.pdf")

# ===== Figure 6: per-tissue DNA->RNA concordance =====
d = pt.dropna(subset=["concordance_rho"]).sort_values("concordance_rho")
fig, ax = plt.subplots(figsize=(S.PANEL * 1.4, S.PANEL * 2))
cols = [ORANGE if t in CULT else (BLUE if v > 0 else GREY) for t, v in zip(d.tissue, d.concordance_rho)]
ax.barh(range(len(d)), d.concordance_rho, color=cols, edgecolor="k", lw=0.3)
ax.set_yticks(range(len(d))); ax.set_yticklabels([f"{t[:30]} (n={n})" for t, n in zip(d.tissue, d.n_carrier)], fontsize=6)
ax.axvline(0, color="k", lw=1.0)
ax.set_xlabel("per-tissue concordance ρ (β_DNA carrier × β_RNA dose, top-300)")
ax.set_title("Per-tissue DNA-RNA concordance\n(orange = cultured cells)")
S.save(fig, f"{FIG}/Figure6_pertissue_concordance.pdf")
print("done figs 3-6")
