#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 62_border_methylation_figure.py — Methylation of 5S array border copies versus interior copies (HiFi).
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
62_border_methylation_figure.py

Methylation of 5S array border copies (first/last copy of each haplotype array)
versus interior copies, from the HiFi data.

Inputs:
  border_copy_methylation_hifi.tsv   per-border-copy
  border_copy_pos_hifi.tsv           gene-anchored 15 bp bins
Interior baseline: copy_methylation_hifi / copy_meth_pos_hifi (5S_rDNA.db)

Output: <FIVES_OUT>/10_border_methylation/border_methylation_hifi.pdf

Paths are read from environment variables (see repository README):
    FIVES_DB    path to 5S_rDNA.db
    FIVES_DATA  directory with input TSVs
    FIVES_OUT   output directory
"""
import os
import sqlite3, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon
from pathlib import Path

DB    = Path(os.environ.get("FIVES_DB", "5S_rDNA.db"))
HPRC  = Path(os.environ.get("FIVES_DATA", "data"))
OUT   = Path(os.environ.get("FIVES_OUT", "output")) / "10_border_methylation"
OUT.mkdir(parents=True, exist_ok=True)

GENE_S = 630          # gene start in interior unit coords (copy_meth_pos wpos)
HYPO   = 65.0
COL_INT, COL_5, COL_3 = "#4c72b0", "#D6604D", "#E8B92E"
pctf = plt.FuncFormatter(lambda v, _: f"{v:.0f}%")

# ── border data ───────────────────────────────────────────────────────────────
pc = pd.read_csv(HPRC / "border_copy_methylation_hifi.tsv", sep="\t")
pc["pct"] = pc["mean_meth"] * 100
b5 = pc[pc.border_note.str.startswith("5") & (pc.n_conf >= 10)]
b3 = pc[pc.border_note.str.startswith("3") & (pc.n_conf >= 10)]

# ── interior baseline from DB ─────────────────────────────────────────────────
con = sqlite3.connect(DB)
interior = pd.read_sql("""
    SELECT a.sample_id AS sample, h.hap_label AS hap, cm.mean_meth*100 AS pct
    FROM copy_methylation_hifi cm JOIN copy c ON cm.copy_id=c.copy_id
    JOIN haplotype h ON c.haplotype_id=h.haplotype_id
    JOIN assembly a ON h.assembly_id=a.assembly_id
    WHERE c.border_note='interior' AND cm.n_conf_calls>=10
""", con)
int_hap = interior.groupby(["sample", "hap"])["pct"].mean().reset_index(name="interior_pct")
# interior gene-anchored positional profile: rel = wpos_bin - GENE_S
ipos = pd.read_sql("""
    SELECT p.wpos_bin, SUM(p.n_meth) nm, SUM(p.n_conf) nc
    FROM copy_meth_pos_hifi p JOIN copy_methylation_hifi cm ON p.copy_id=cm.copy_id
    JOIN copy c ON cm.copy_id=c.copy_id
    WHERE c.border_note='interior' GROUP BY p.wpos_bin
""", con)
con.close()
ipos["rel"] = ipos["wpos_bin"] - GENE_S
ipos["pct"] = ipos["nm"] / ipos["nc"] * 100

# ── figure ────────────────────────────────────────────────────────────────────
fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(15, 4.6))
fig.suptitle("5S array border-copy CpG methylation (HiFi)", fontsize=11, y=1.02)

# A: overall methylation distribution by copy class
groups = [("Interior", interior["pct"].values, COL_INT),
          (f"5′ border\n(n={len(b5)})", b5["pct"].values, COL_5),
          (f"3′ border\n(n={len(b3)})", b3["pct"].values, COL_3)]
parts = axA.violinplot([g[1] for g in groups], showextrema=False, points=200)
for pcb, (_, _, c) in zip(parts["bodies"], groups):
    pcb.set_facecolor(c); pcb.set_alpha(0.55); pcb.set_edgecolor(c)
for i, (lab, vals, c) in enumerate(groups, 1):
    q25, q50, q75 = np.percentile(vals, [25, 50, 75])
    axA.plot([i, i], [q25, q75], color=c, lw=2)
    axA.plot(i, q50, "o", color="white", mec=c, mew=1.4, ms=6, zorder=5)
axA.axhline(HYPO, color="grey", ls="--", lw=0.8)
axA.set_xticks([1, 2, 3]); axA.set_xticklabels([g[0] for g in groups])
axA.set_ylim(0, 100); axA.yaxis.set_major_formatter(pctf)
axA.set_ylabel("CpG methylation (%)")
axA.set_title("(A) Overall methylation by copy class")
axA.grid(axis="y", lw=0.3, alpha=0.4)

# B: gene-anchored positional profile, 5' border vs interior
bp = pd.read_csv(HPRC / "border_copy_pos_hifi.tsv", sep="\t")
p5 = (bp[bp.border_note.str.startswith("5")].groupby("rel_bin")
      .agg(nm=("n_meth", "sum"), nc=("n_conf", "sum")).reset_index())
p5 = p5[p5.nc >= 20]; p5["pct"] = p5.nm / p5.nc * 100
ip = ipos[(ipos.rel >= -1450) & (ipos.rel <= 1500) & (ipos.nc >= 50)]
axB.axvspan(0, 119, color="#aec6cf", alpha=0.5, zorder=0)
axB.text(60, 102, "5S gene", ha="center", fontsize=7, color="#444")
axB.plot(ip["rel"], ip["pct"], color=COL_INT, lw=2, label="Interior copies")
axB.plot(p5["rel_bin"], p5["pct"], color=COL_5, lw=2, label="5′ border copies")
axB.set_xlim(-1450, 1500); axB.set_ylim(0, 105)
axB.yaxis.set_major_formatter(pctf)
axB.set_xlabel("Position relative to 5S gene start (bp)")
axB.set_ylabel("CpG methylation (%)")
axB.set_title("(B) Gene-anchored profile: 5′ border vs interior")
axB.legend(fontsize=8, loc="center right"); axB.grid(lw=0.3, alpha=0.4)

# C: paired 5' border vs own-haplotype interior mean
m = b5.merge(int_hap, on=["sample", "hap"], how="inner")
axC.scatter(m["interior_pct"], m["pct"], s=18, alpha=0.6, color=COL_5, edgecolor="none")
axC.plot([0, 100], [0, 100], "k--", lw=0.8, alpha=0.6)
w, p = wilcoxon(m["pct"], m["interior_pct"])
axC.set_xlim(0, 100); axC.set_ylim(0, 100)
axC.xaxis.set_major_formatter(pctf); axC.yaxis.set_major_formatter(pctf)
axC.set_xlabel("Interior mean methylation (same haplotype)")
axC.set_ylabel("5′ border methylation")
axC.set_title(f"(C) Paired per-haplotype (n={len(m)})\n"
              f"5′<interior in {(m['pct']<m['interior_pct']).sum()}/{len(m)}  ·  "
              f"Wilcoxon p={p:.1e}")
axC.grid(lw=0.3, alpha=0.4)

plt.tight_layout()
out = OUT / "border_methylation_hifi.pdf"
fig.savefig(out, dpi=200, bbox_inches="tight")
plt.close(fig)
# data tables
(OUT / "data").mkdir(exist_ok=True)
pc.to_csv(OUT / "data/border_copy_methylation_hifi.tsv", sep="\t", index=False)
m.round(2).to_csv(OUT / "data/border5_vs_interior_paired.tsv", sep="\t", index=False)
print(f"Saved: {out}")
print(f"5′ border mean={b5['pct'].mean():.1f}%  3′ border mean={b3['pct'].mean():.1f}%  "
      f"interior mean={interior['pct'].mean():.1f}%")
print(f"paired n={len(m)}, 5′<interior {(m['pct']<m['interior_pct']).sum()}, p={p:.2e}")
