#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 41_methylation_overview_hifi.py — Four-panel overview of per-copy CpG methylation across the 5S array (HiFi).
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
41_methylation_overview_hifi.py

Four-panel overview of per-copy CpG methylation across the 5S array
(215 HPRC probands, HiFi).

Panel A  — per-copy methylation distribution (High / Intermediate / Low)
Panel B  — along-array gradient: median + IQR + 10–90th pctile (hap1+hap2 pooled)
Panel C  — n low-meth copies per donor vs. total copy number
Panel D  — regional methylation by copy class
           If alu_n/alu_meth columns are populated:
             three bars: 5S gene | ALU SINE | other NTS
           Else (3-region mode):
             NTS-pre | gene | NTS-post

Input : 5S_rDNA.db (tables copy_methylation_hifi, copy, haplotype, assembly).
Output: <FIVES_OUT>/03_methylation_full215/00_methylation_overview_hifi.pdf

Paths are read from environment variables (see repository README):
    FIVES_DB   path to 5S_rDNA.db
    FIVES_OUT  output directory
"""

import os
import sqlite3, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats as sc
from pathlib import Path

DB     = Path(os.environ.get("FIVES_DB", "5S_rDNA.db"))
OUTDIR = Path(os.environ.get("FIVES_OUT", "output")) / "03_methylation_full215"

# ── load data ─────────────────────────────────────────────────────────────────
con = sqlite3.connect(DB)

# detect ALU columns
has_alu = con.execute(
    "SELECT COUNT(*) FROM pragma_table_info('copy_methylation_hifi') "
    "WHERE name='alu_n'"
).fetchone()[0] > 0
alu_populated = has_alu and con.execute(
    "SELECT SUM(alu_n) FROM copy_methylation_hifi").fetchone()[0] > 0

df = pd.read_sql("""
    SELECT cm.copy_id, cm.n_conf_calls, cm.mean_meth,
           cm.nts_pre_n, cm.nts_pre_meth,
           cm.gene_n, cm.gene_meth,
           cm.nts_post_n, cm.nts_post_meth,
           COALESCE(cm.alu_n,   0) AS alu_n,
           COALESCE(cm.alu_meth,0) AS alu_meth,
           c.copy_number, h.hap_label, h.haplotype_id, h.array_order_resolved,
           a.sample_id, a.superpopulation
    FROM copy_methylation_hifi cm
    JOIN copy c ON cm.copy_id=c.copy_id
    JOIN haplotype h ON c.haplotype_id=h.haplotype_id
    JOIN assembly a ON h.assembly_id=a.assembly_id
    WHERE c.border_note='interior' AND cm.n_conf_calls>=10
      AND a.cohort IN ('HPRC_Year1','HPRC_Release2')
""", con)
con.close()

def safe_frac(n, d):
    return np.where(d > 0, n / d * 100, np.nan)

df["meth_pct"]     = df["mean_meth"] * 100
df["nts_pre_pct"]  = safe_frac(df["nts_pre_meth"],  df["nts_pre_n"])
df["gene_pct"]     = safe_frac(df["gene_meth"],      df["gene_n"])
df["nts_post_pct"] = safe_frac(df["nts_post_meth"],  df["nts_post_n"])
df["alu_pct"]      = safe_frac(df["alu_meth"],       df["alu_n"])
# other NTS = nts_pre + (nts_post minus ALU)
df["other_nts_n"]    = df["nts_pre_n"] + (df["nts_post_n"] - df["alu_n"]).clip(lower=0)
df["other_nts_meth"] = df["nts_pre_meth"] + (df["nts_post_meth"] - df["alu_meth"]).clip(lower=0)
df["other_nts_pct"]  = safe_frac(df["other_nts_meth"], df["other_nts_n"])

df["copy_class"] = pd.cut(df["meth_pct"], bins=[-1, 35, 65, 101],
                           labels=["Low", "Intermediate", "High"])
# along-array position: rank by copy_number (cross-contig order).
# Exclude partial-order haplotypes.
df.loc[df["array_order_resolved"] == "partial", "copy_number"] = np.nan
df["pct_pos"] = (df.groupby("haplotype_id")["copy_number"]
                   .rank(method="first").sub(1)
                 / df.groupby("haplotype_id")["copy_number"]
                      .transform("count").sub(1) * 100)

n_samples = df["sample_id"].nunique()
n_copies  = len(df)
n_calls   = df["n_conf_calls"].sum()

CLS_COLORS = {"High": "#2166AC", "Intermediate": "#F7CA18", "Low": "#D6604D"}

# ── figure layout ─────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(16, 12))
gs  = gridspec.GridSpec(2, 2, figure=fig, hspace=0.42, wspace=0.35)
axA = fig.add_subplot(gs[0, 0])
axB = fig.add_subplot(gs[0, 1])
axC = fig.add_subplot(gs[1, 0])
axD = fig.add_subplot(gs[1, 1])

pct_fmt = plt.FuncFormatter(lambda x, _: f"{x:.0f}%")
fig.suptitle(
    f"5S rDNA CpG methylation — overview\n"
    f"HPRC Year-1 + Release-2  ·  {n_samples} probands, 428 haplotypes, "
    f"{n_copies:,} interior copies, {n_calls/1e6:.0f}M confident calls",
    fontsize=10, y=1.01
)

# ── A: bimodal distribution ────────────────────────────────────────────────────
bins_hist = np.linspace(0, 100, 51)
for cls, color in CLS_COLORS.items():
    sub = df[df["copy_class"] == cls]["meth_pct"]
    axA.hist(sub, bins=bins_hist, color=color, alpha=0.75, label=cls,
             edgecolor="none", density=False)
axA.axvline(35, color="k", ls="--", lw=0.9, alpha=0.6)
axA.axvline(65, color="k", ls="--", lw=0.9, alpha=0.6)
pct_high = (df["copy_class"] == "High").mean() * 100
pct_low  = (df["copy_class"] == "Low").mean() * 100
pct_int  = (df["copy_class"] == "Intermediate").mean() * 100
axA.set_xlabel("Per-copy mean CpG methylation (%)", fontsize=9)
axA.set_ylabel("Number of copies", fontsize=9)
axA.set_title(
    f"(A) Per-copy methylation distribution\n"
    f"High≥65%: {pct_high:.0f}%  |  Low<35%: {pct_low:.0f}%  |  Inter: {pct_int:.0f}%",
    fontsize=9)
axA.xaxis.set_major_formatter(pct_fmt)
axA.legend(fontsize=8, framealpha=0.85)
axA.tick_params(labelsize=8)
axA.grid(lw=0.3, alpha=0.4, axis="y")

# ── B: along-array gradient ────────────────────────────────────────────────────
N_BINS = 20
bin_edges = np.linspace(0, 100, N_BINS + 1)
PCTILES   = [10, 25, 50, 75, 90]
mids=[]; q10s=[]; q25s=[]; q50s=[]; q75s=[]; q90s=[]
for i in range(len(bin_edges) - 1):
    msk = (df["pct_pos"] >= bin_edges[i]) & (df["pct_pos"] < bin_edges[i+1])
    sub = df.loc[msk, "meth_pct"].dropna().values
    if len(sub) < 5:
        continue
    mids.append((bin_edges[i] + bin_edges[i+1]) / 2)
    q10, q25, q50, q75, q90 = np.percentile(sub, PCTILES)
    q10s.append(q10); q25s.append(q25); q50s.append(q50)
    q75s.append(q75); q90s.append(q90)
mids = np.array(mids)
q10s=np.array(q10s); q25s=np.array(q25s); q50s=np.array(q50s)
q75s=np.array(q75s); q90s=np.array(q90s)
axB.fill_between(mids, q10s, q90s, color="#4c72b0", alpha=0.14, label="10th–90th pctile")
axB.fill_between(mids, q25s, q75s, color="#4c72b0", alpha=0.35, label="IQR (25th–75th)")
axB.plot(mids, q50s, color="#4c72b0", lw=2.4, label="Median")
axB.set_xlim(0, 100); axB.set_ylim(0, 100)
axB.set_xlabel("Position along array (copy rank 0–100%)", fontsize=9)
axB.set_ylabel("CpG methylation per copy (%)", fontsize=9)
axB.set_title("(B) Along-array gradient\nhap1+hap2 pooled  |  median + IQR + 10–90th pctile",
              fontsize=9)
axB.xaxis.set_major_formatter(pct_fmt); axB.yaxis.set_major_formatter(pct_fmt)
axB.tick_params(labelsize=8); axB.grid(lw=0.3, alpha=0.4)
axB.legend(fontsize=8, framealpha=0.85, loc="upper right")

# ── C: set-point ─────────────────────────────────────────────────────────────
sp = (df.groupby("sample_id")
       .agg(n_low=("copy_class", lambda x: (x == "Low").sum()),
            n_total=("copy_id", "count"))
       .reset_index())
axC.scatter(sp["n_total"], sp["n_low"], s=18, alpha=0.55, color="#555555", zorder=3)
slope, intercept, r, p, _ = sc.linregress(sp["n_total"], sp["n_low"])
xs = np.array([sp["n_total"].min(), sp["n_total"].max()])
axC.plot(xs, xs * slope + intercept, "k--", lw=1.2, label=f"fit (slope={slope:.2f})")
axC.axhline(sp["n_low"].mean(), color="#D6604D", lw=1.5,
            label=f"mean = {sp['n_low'].mean():.0f}")
axC.set_xlabel("Total interior copies (≥10 calls)", fontsize=9)
axC.set_ylabel("Low-methylation copies (<35%)", fontsize=9)
axC.set_title(
    f"(C) Low-methylation copies vs total copy number\n"
    f"mean={sp['n_low'].mean():.0f}, CV={sp['n_low'].std()/sp['n_low'].mean():.2f}  "
    f"|  r={r:.2f}, slope={slope:.2f}",
    fontsize=9)
axC.tick_params(labelsize=8); axC.grid(lw=0.3, alpha=0.4)
axC.legend(fontsize=8, framealpha=0.85)

# ── D: regional methylation by copy class ─────────────────────────────────────
if alu_populated:
    regions = ["5S gene\n(119 bp)", "ALU SINE\n(280 bp)", "Other NTS\n(1769 bp)"]
    col_map = {"5S gene\n(119 bp)": "gene_pct",
               "ALU SINE\n(280 bp)": "alu_pct",
               "Other NTS\n(1769 bp)": "other_nts_pct"}
    title_suffix = "5S gene | ALU SINE | other NTS"
else:
    regions = ["NTS-pre\n(629 bp)", "5S gene\n(119 bp)", "NTS-post\n(1419 bp)"]
    col_map = {"NTS-pre\n(629 bp)": "nts_pre_pct",
               "5S gene\n(119 bp)": "gene_pct",
               "NTS-post\n(1419 bp)": "nts_post_pct"}
    title_suffix = "NTS-pre | gene | NTS-post"

reg_data = {}
for cls in ["High", "Intermediate", "Low"]:
    sub = df[df["copy_class"] == cls]
    reg_data[cls] = {r: sub[col_map[r]].dropna().mean() for r in regions}

x = np.arange(len(regions)); width = 0.25
for i, cls in enumerate(["High", "Intermediate", "Low"]):
    vals = [reg_data[cls][r] for r in regions]
    axD.bar(x + (i - 1) * width, vals, width,
            color=CLS_COLORS[cls], alpha=0.85, label=cls,
            edgecolor="white", lw=0.5)
    for j, v in enumerate(vals):
        if not np.isnan(v):
            axD.text(x[j] + (i - 1) * width, v + 1.5, f"{v:.0f}%",
                     ha="center", va="bottom", fontsize=6.5, color="#333333")

axD.set_xticks(x)
axD.set_xticklabels(regions, fontsize=8)
axD.set_ylabel("Mean CpG methylation (%)", fontsize=9)
axD.set_ylim(0, 100)
axD.yaxis.set_major_formatter(pct_fmt)
axD.set_title(f"(D) Regional methylation by copy class\n{title_suffix}", fontsize=9)
axD.legend(fontsize=8, framealpha=0.85)
axD.tick_params(labelsize=8)
axD.grid(lw=0.3, alpha=0.4, axis="y")

plt.tight_layout()
out = OUTDIR / "00_methylation_overview_hifi.pdf"
fig.savefig(out, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"Saved: {out}")
print(f"Panel D mode: {'ALU split' if alu_populated else 'legacy 3-region (re-run script 37 for ALU)'}")
