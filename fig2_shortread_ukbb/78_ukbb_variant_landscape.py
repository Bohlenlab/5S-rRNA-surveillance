#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 78_ukbb_variant_landscape.py — descriptive figure of the UK Biobank 5S rDNA
# variant landscape at a 0.30% VAF threshold (position, spectrum, rarity).
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
78_ukbb_variant_landscape.py

Descriptive landscape of UKBB 5S rDNA variant calls at primary threshold 0.30% VAF,
with comparisons to neighboring thresholds (0.15%, 0.20%, 0.50%).

Output: figures/02_variant_calling_qc/78_ukbb_variant_landscape.pdf
"""

import os
import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D

# ── paths ─────────────────────────────────────────────────────────────────────

T2T    = Path(os.environ.get("FIVES_DATA", "data"))
DB     = Path(os.environ.get("FIVES_DB", "5S_rDNA.db"))
OUTDIR = Path(os.environ.get("FIVES_OUT", "output")) / "02_variant_calling_qc"
OUTDIR.mkdir(parents=True, exist_ok=True)

WIN_LO, WIN_HI = 467, 967
GENE_LO, GENE_HI = 630, 748
N_TOTAL = 490_075

THRESH_MAIN = 0.003
THRESHOLDS  = [0.001, 0.0015, 0.002, 0.003, 0.005, 0.010]
THRESH_LABS = ["0.10%", "0.15%", "0.20%", "0.30%", "0.50%", "1.00%"]

GENE_COLOR = "#fff3b0"
REG_COLORS = {"gene": "#e67e22", "nts_pre": "#2980b9", "nts_post": "#27ae60"}
REG_ORDER  = ["gene", "nts_pre", "nts_post"]
REG_XLABS  = ["Gene\n(630–748)", "NTS pre\n(467–629)", "NTS post\n(749–967)"]

# ── load data ─────────────────────────────────────────────────────────────────

print("Loading UKBB data …", flush=True)
con = sqlite3.connect(DB)

raw = con.execute("""
    SELECT t2t_pos, t2t_ref, t2t_alt, region,
           n_carriers_ad1, n_carriers_ad5, mean_vaf, median_vaf, vaf_array
    FROM ukbb_population_variants
    ORDER BY t2t_pos
""").fetchall()

depth_df = pd.DataFrame(con.execute("""
    SELECT t2t_pos, median_dp, p5_dp, p25_dp, p75_dp, p95_dp
    FROM ukbb_depth_profile ORDER BY t2t_pos
""").fetchall(), columns=["t2t_pos","median_dp","p5_dp","p25_dp","p75_dp","p95_dp"])

# True positives: union of assembly-confirmed and HiFi-confirmed variants
tp_rows = con.execute("""
    SELECT DISTINCT consensus_pos, ref, alt FROM (
        SELECT v.consensus_pos, v.ref, v.alt
        FROM variant v
        JOIN copy c ON v.copy_id = c.copy_id
        JOIN haplotype h ON c.haplotype_id = h.haplotype_id
        JOIN assembly a ON h.assembly_id = a.assembly_id
        WHERE v.alignment_source = 'gene_unit_t2t'
          AND a.cohort = 'HPRC_Year1'
          AND v.consensus_pos BETWEEN ? AND ?
        UNION
        SELECT rv.consensus_pos, rv.ref, rv.alt
        FROM read_variant rv
        JOIN assembly a ON rv.assembly_id = a.assembly_id
        WHERE rv.modality = 'hifi'
          AND rv.vaf IS NOT NULL
          AND a.cohort = 'HPRC_Year1'
          AND rv.consensus_pos BETWEEN ? AND ?
    )
""", (WIN_LO, WIN_HI, WIN_LO, WIN_HI)).fetchall()
tp_set = {(int(r[0]), r[1], r[2]) for r in tp_rows}
con.close()
print(f"  True positives (assembly ∪ HiFi): {len(tp_set):,}")

records = []
all_vafs_main = []   # pooled VAF values for panel F

for pos, ref, alt, reg, n1, n5, mv, mdv, blob in raw:
    vafs = np.frombuffer(blob, dtype=np.float32)
    rec  = dict(t2t_pos=int(pos), ref=ref, alt=alt, region=reg or "unknown",
                n_ad1=int(n1), n_ad5=int(n5),
                mean_vaf=float(mv or 0), median_vaf=float(mdv or 0))
    for lab, T in zip(THRESH_LABS, THRESHOLDS):
        nc = int(len(vafs) - np.searchsorted(vafs, T))
        rec[f"n_{lab}"] = nc
    records.append(rec)
    # collect VAFs above primary threshold for distribution plot
    above = vafs[vafs >= THRESH_MAIN]
    if len(above):
        all_vafs_main.append((reg or "unknown", above))

df = pd.DataFrame(records)
df["is_tp"] = [
    (int(r["t2t_pos"]), r["ref"], r["alt"]) in tp_set for r in records
]
df_called = df[df["n_0.30%"] >= 1].copy()
all_vafs_main_arr = np.concatenate([a for _, a in all_vafs_main]) if all_vafs_main else np.array([])

print(f"  {len(df):,} total variants, {len(df_called):,} called at 0.30%")
print(f"  TP variants present in UKBB: {df['is_tp'].sum():,}  "
      f"({df_called['is_tp'].sum():,} also called at 0.30%)")
for lab in THRESH_LABS:
    n = (df[f"n_{lab}"] >= 1).sum()
    print(f"    {lab}: {n:,} variants")

# ── substitution types ────────────────────────────────────────────────────────

_COMP = {"A":"T","T":"A","C":"G","G":"C"}

def canonical_snv(ref, alt):
    if len(ref) != 1 or len(alt) != 1 or ref == alt:
        return "other"
    if ref in ("C","T"):
        return f"{ref}>{alt}"
    return f"{_COMP.get(ref,'?')}>{_COMP.get(alt,'?')}"

df_snv = df_called.copy()
df_snv["snv_type"] = [canonical_snv(r, a) for r, a in zip(df_snv["ref"], df_snv["alt"])]
df_snv = df_snv[df_snv["snv_type"] != "other"]

SNV_ORDER  = ["C>A","C>G","C>T","T>A","T>C","T>G"]
SNV_COLORS = ["#1f77b4","#000000","#d62728","#cccccc","#2ca02c","#ff7f0e"]

# ── rarity buckets at 0.30% ───────────────────────────────────────────────────

BUCK_EDGES  = [1, 2, 11, 101, 1001, N_TOTAL + 1]
BUCK_LABELS = ["1 (singleton)", "2–10", "11–100", "101–1000", ">1000"]
BUCK_COLORS = ["#e41a1c","#ff7f00","#4daf4a","#377eb8","#984ea3"]

def bucket(n):
    for i in range(len(BUCK_EDGES) - 1):
        if BUCK_EDGES[i] <= n < BUCK_EDGES[i+1]:
            return i
    return len(BUCK_EDGES) - 2

df_called["bucket"] = df_called["n_0.30%"].apply(bucket)

# ── figure layout ─────────────────────────────────────────────────────────────

print("Generating figure …", flush=True)

fig = plt.figure(figsize=(21, 17))
gs  = gridspec.GridSpec(3, 3, figure=fig, hspace=0.52, wspace=0.38,
                        top=0.93, bottom=0.06, left=0.06, right=0.97)

fig.suptitle(
    "UKBB 5S rDNA variant landscape  ·  n = 490,075 samples  ·  Primary threshold: 0.30% VAF",
    fontsize=12, fontweight="bold")

# ── A: position scatter + depth (full width) ─────────────────────────────────

ax_A = fig.add_subplot(gs[0, :])

for reg in REG_ORDER:
    sub     = df_called[(df_called["region"] == reg) & ~df_called["is_tp"]]
    sub_tp  = df_called[(df_called["region"] == reg) &  df_called["is_tp"]]
    for subset, edge, lw, zo in [(sub, "none", 0, 3), (sub_tp, "red", 0.8, 4)]:
        sizes = np.clip(np.log10(subset["n_0.30%"].values + 1) * 12, 4, 80)
        ax_A.scatter(subset["t2t_pos"], subset["n_0.30%"],
                     s=sizes, c=REG_COLORS[reg], alpha=0.70, zorder=zo,
                     edgecolors=edge, linewidths=lw,
                     label=f"{reg} (n={len(sub)+len(sub_tp):,})" if edge == "none" else None)

# TP outline legend entry
ax_A.scatter([], [], s=20, c="grey", edgecolors="red", linewidths=0.8,
             label=f"assembly ∪ HiFi confirmed (n={df_called['is_tp'].sum():,})")

ax_A.axvspan(GENE_LO, GENE_HI, color=GENE_COLOR, alpha=0.5, zorder=1)
ax_A.set_yscale("log")
ax_A.set_xlim(WIN_LO - 5, WIN_HI + 5)
ax_A.set_xlabel("T2T consensus position")
ax_A.set_ylabel("Carriers at VAF ≥ 0.30%  (log)")
ax_A.set_title(
    f"A  Variant map: {len(df_called):,} variants × {N_TOTAL:,} samples  "
    f"(marker size ∝ log₁₀ carrier count;  black outline = assembly ∪ HiFi confirmed)",
    fontsize=10, loc="left")
ax_A.legend(fontsize=8, title="Region", title_fontsize=8, loc="upper right")
ax_A.grid(True, lw=0.3, alpha=0.4)

ax_A2 = ax_A.twinx()
ax_A2.fill_between(depth_df["t2t_pos"], depth_df["p25_dp"], depth_df["p75_dp"],
                   alpha=0.12, color="grey")
ax_A2.plot(depth_df["t2t_pos"], depth_df["median_dp"],
           color="grey", lw=1, alpha=0.5, label="Median depth")
ax_A2.set_ylabel("Read depth (right axis)", color="grey", fontsize=8)
ax_A2.tick_params(axis="y", labelcolor="grey", labelsize=7)
ax_A2.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v/1e3:.1f}k"))

# ── B: allele frequency spectrum ──────────────────────────────────────────────

ax_B = fig.add_subplot(gs[1, 0])

bins_afs = np.logspace(0, np.log10(N_TOTAL), 55)

style_map = {
    "0.15%": dict(histtype="step",  lw=1.5, color="#1f77b4"),
    "0.30%": dict(histtype="stepfilled", alpha=0.45, lw=2.0, color="#d62728"),
    "0.50%": dict(histtype="step",  lw=1.5, color="#2ca02c"),
}
for lab, kw in style_map.items():
    vals = df[df[f"n_{lab}"] >= 1][f"n_{lab}"].values
    ax_B.hist(vals, bins=bins_afs, label=f"{lab} (n={len(vals):,})", **kw)

ax_B.set_xscale("log")
ax_B.set_xlabel("Carriers (log)")
ax_B.set_ylabel("Variants")
ax_B.set_title("B  Allele frequency spectrum", fontsize=10, loc="left")
ax_B.legend(fontsize=8); ax_B.grid(True, lw=0.3, alpha=0.4)
ax_B.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{int(v):,}"))

# ── C: per-region breakdown ───────────────────────────────────────────────────

ax_C = fig.add_subplot(gs[1, 1])

var_counts   = [len(df_called[df_called["region"] == r])            for r in REG_ORDER]
med_carriers = [df_called[df_called["region"] == r]["n_0.30%"].median() for r in REG_ORDER]
# Variants per nt window length
win_lengths  = [(GENE_HI - GENE_LO + 1), (GENE_LO - WIN_LO), (WIN_HI - GENE_HI)]

x = np.arange(3)
bars_C = ax_C.bar(x, var_counts, 0.5,
                  color=[REG_COLORS[r] for r in REG_ORDER], alpha=0.85)
for bar, cnt, wl in zip(bars_C, var_counts, win_lengths):
    ax_C.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 3,
              f"{cnt}\n({cnt/wl:.1f}/nt)", ha="center", va="bottom", fontsize=7.5)

ax_C.set_xticks(x); ax_C.set_xticklabels(REG_XLABS, fontsize=8)
ax_C.set_ylabel("Variants called at 0.30%")
ax_C.set_title("C  Variant count per region  (rate per nt in labels)", fontsize=10, loc="left")

ax_C2 = ax_C.twinx()
ax_C2.plot(x, med_carriers, "D--", color="black", lw=1.5, ms=7, zorder=5)
ax_C2.set_ylabel("Median carrier count  (◆)", color="black", fontsize=8)
ax_C2.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
ax_C.grid(True, lw=0.3, alpha=0.4, axis="y")

# ── D: substitution spectrum ──────────────────────────────────────────────────

ax_D = fig.add_subplot(gs[1, 2])

# count by snv type × region
snv_by_reg = {reg: {} for reg in REG_ORDER}
for reg in REG_ORDER:
    sub = df_snv[df_snv["region"] == reg]
    for t in SNV_ORDER:
        snv_by_reg[reg][t] = int((sub["snv_type"] == t).sum())

x_D = np.arange(len(SNV_ORDER))
bar_w = 0.25
for i, reg in enumerate(REG_ORDER):
    counts = [snv_by_reg[reg].get(t, 0) for t in SNV_ORDER]
    ax_D.bar(x_D + (i - 1) * bar_w, counts, bar_w,
             color=REG_COLORS[reg], alpha=0.85, label=reg)

ax_D.set_xticks(x_D); ax_D.set_xticklabels(SNV_ORDER, fontsize=9)
ax_D.set_xlabel("Canonical SNV type (pyrimidine reference)")
ax_D.set_ylabel("Variants called at 0.30%")
ax_D.set_title("D  Substitution spectrum by region", fontsize=10, loc="left")
ax_D.legend(fontsize=7); ax_D.grid(True, lw=0.3, alpha=0.4, axis="y")

# ── E: rarity distribution at 0.30% ──────────────────────────────────────────

ax_E = fig.add_subplot(gs[2, 0])

bucket_counts = [int((df_called["bucket"] == i).sum()) for i in range(len(BUCK_LABELS))]
total_called  = len(df_called)

bars_E = ax_E.bar(range(len(BUCK_LABELS)), bucket_counts, color=BUCK_COLORS, alpha=0.85)
for bar, cnt in zip(bars_E, bucket_counts):
    pct = 100 * cnt / total_called
    ax_E.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
              f"{cnt}\n({pct:.1f}%)", ha="center", va="bottom", fontsize=7.5)

ax_E.set_xticks(range(len(BUCK_LABELS)))
ax_E.set_xticklabels(BUCK_LABELS, fontsize=8)
ax_E.set_xlabel("Carrier count at 0.30%")
ax_E.set_ylabel("Variants")
ax_E.set_title("E  Rarity distribution at 0.30% VAF", fontsize=10, loc="left")
ax_E.grid(True, lw=0.3, alpha=0.4, axis="y")

# ── F: VAF distribution of carriers (pooled) ─────────────────────────────────

ax_F = fig.add_subplot(gs[2, 1])

bins_vaf = np.linspace(THRESH_MAIN, min(0.50, float(all_vafs_main_arr.max()) + 0.01), 80)
for reg, color in REG_COLORS.items():
    these = np.concatenate([a for r, a in all_vafs_main if r == reg]) \
            if any(r == reg for r, _ in all_vafs_main) else np.array([])
    if len(these):
        ax_F.hist(these, bins=bins_vaf, color=color, alpha=0.55,
                  label=f"{reg} (n={len(these):,})", density=True)

ax_F.axvline(THRESH_MAIN, color="red",  lw=1.5, ls="--", label="0.30% threshold")
ax_F.axvline(0.0015,      color="blue", lw=1.0, ls="--", alpha=0.7, label="0.15% UKB AD5")
ax_F.set_xlabel("VAF per carrier  (VAF ≥ 0.30%)")
ax_F.set_ylabel("Density")
ax_F.set_title("F  Per-carrier VAF distribution (pooled across all variants)", fontsize=10, loc="left")
ax_F.legend(fontsize=7); ax_F.grid(True, lw=0.3, alpha=0.4)
ax_F.set_xlim(0, min(0.55, bins_vaf[-1] + 0.02))
ax_F.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v*100:.1f}%"))

# ── G: top 25 most-frequent variants ─────────────────────────────────────────

ax_G = fig.add_subplot(gs[2, 2])

top25 = df_called.nlargest(25, "n_0.30%").sort_values("n_0.30%")
labels_G = [f"pos{r.t2t_pos} {r.ref}>{r.alt}  [{r.region}]" for r in top25.itertuples()]
colors_G = [REG_COLORS.get(r.region, "grey") for r in top25.itertuples()]
n_main_col = top25["n_0.30%"].values
ax_G.barh(range(len(top25)), n_main_col, color=colors_G, alpha=0.85)
ax_G.set_yticks(range(len(top25)))
ax_G.set_yticklabels(labels_G, fontsize=6.5)
ax_G.set_xlabel("Carriers at VAF ≥ 0.30%")
ax_G.set_title("G  Top 25 most frequent variants at 0.30%", fontsize=10, loc="left")
ax_G.grid(True, lw=0.3, alpha=0.4, axis="x")
max_n = int(n_main_col.max())
for i, (n, pct) in enumerate(zip(n_main_col, 100 * n_main_col / N_TOTAL)):
    ax_G.text(n + max_n * 0.01, i, f" {pct:.2f}%", va="center", fontsize=6)

plt.savefig(OUTDIR / "78_ukbb_variant_landscape.pdf", dpi=180, bbox_inches="tight")
plt.close()
print(f"\nFigure → {OUTDIR}/78_ukbb_variant_landscape.pdf")

# ── print top 25 to stdout ─────────────────────────────────────────────────────
print("\n── Top 25 variants at 0.30% VAF ────────────────────────────────────────────")
top25_out = df_called.nlargest(25, "n_0.30%")[
    ["t2t_pos","ref","alt","region","n_0.30%","n_0.15%","n_0.50%","mean_vaf"]].copy()
top25_out["freq_%"] = (top25_out["n_0.30%"] / N_TOTAL * 100).round(3)
print(top25_out.to_string(index=False))
