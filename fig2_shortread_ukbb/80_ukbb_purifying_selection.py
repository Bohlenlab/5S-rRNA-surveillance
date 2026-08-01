#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 80_ukbb_purifying_selection.py — carrier-frequency association tests between
# incorporation efficiency and UK Biobank 5S variant carrier counts at three VAF
# thresholds, with a statistical summary table.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
80_ukbb_purifying_selection.py

Tests whether purifying selection against incorporation-defective 5S rRNA
variants is detectable in the UKBB population at three VAF thresholds.

Four complementary tests:
  1. Carrier count (MW, binary defective/normal split)
  2. Spearman correlation (continuous incorp score vs carrier count)
  3. Conditional carrier count (MW, detected variants only)
  4. Mean VAF per carrier (MW, detected variants only)

Layout:
  Row 0 (A/B/C): boxplot+strip of carrier count at 0.30% / 1.00% / 3.00%
  Row 1 (D/E)  : binned-median trend of incorp vs carrier count at 0.30% and 3.00%
  Row 1 (F)    : statistical summary table

Output: figures/02_variant_calling_qc/80_ukbb_purifying_selection.pdf
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
from matplotlib.patches import Rectangle
from scipy.stats import mannwhitneyu, spearmanr

# ── paths ─────────────────────────────────────────────────────────────────────

T2T    = Path(os.environ.get("FIVES_DATA", "data"))
DB     = Path(os.environ.get("FIVES_DB", "5S_rDNA.db"))
OUTDIR = Path(os.environ.get("FIVES_OUT", "output")) / "02_variant_calling_qc"
OUTDIR.mkdir(parents=True, exist_ok=True)

N_TOTAL    = 490_075
INCORP_CUT = 0.50
THRESHOLDS = [("0.30%", 0.003), ("1.00%", 0.010), ("3.00%", 0.030)]
COL_DEF    = "#d62728"
COL_NOR    = "#1f77b4"

# ── load ──────────────────────────────────────────────────────────────────────

print("Loading …", flush=True)
con = sqlite3.connect(DB)
rows = con.execute("""
    SELECT fa.consensus_pos, fa.ref_base, fa.alt_base,
           fa.incorp_60s_mean, fa.rna_expr_mean, uv.vaf_array
    FROM functional_annotation fa
    JOIN ukbb_population_variants uv
      ON uv.t2t_pos = fa.consensus_pos
     AND uv.t2t_ref = fa.ref_base
     AND uv.t2t_alt = fa.alt_base
    WHERE fa.incorp_60s_mean IS NOT NULL
""").fetchall()
con.close()

data = []
for pos, ref, alt, incorp, expr, blob in rows:
    vafs   = np.frombuffer(blob, dtype=np.float32).copy()
    is_def = incorp < INCORP_CUT
    rec    = dict(pos=int(pos), ref=ref, alt=alt,
                  incorp=float(incorp), is_defective=is_def)
    for lab, T in THRESHOLDS:
        above = vafs[vafs >= T]
        rec[f"n_{lab}"]        = len(above)
        rec[f"mean_vaf_{lab}"] = float(above.mean()) if len(above) else np.nan
    data.append(rec)

df     = pd.DataFrame(data)
df_def = df[df["is_defective"]]
df_nor = df[~df["is_defective"]]
print(f"  {len(df_def)} defective  {len(df_nor)} normal")

# ── compute all stats (for summary table) ─────────────────────────────────────

results = {}   # (lab, test) → dict
for lab, T in THRESHOLDS:
    col   = f"n_{lab}"
    d_all = df_def[col].values.astype(float)
    n_all = df_nor[col].values.astype(float)

    # Test 1: MW on full carrier count
    _, p1    = mannwhitneyu(d_all, n_all, alternative="less")
    med_d, med_n = np.median(d_all), np.median(n_all)
    depl     = 100*(med_n - med_d)/med_n if med_n > 0 else np.nan

    # Test 2: Spearman (continuous)
    rho, p2  = spearmanr(df["incorp"].values, df[col].values)

    # Test 3: conditional MW (detected only, n > 0)
    det      = df[df[col] > 0]
    d_det    = det[det["is_defective"]][col].values.astype(float)
    n_det    = det[~det["is_defective"]][col].values.astype(float)
    if len(d_det) > 1 and len(n_det) > 1:
        _, p3 = mannwhitneyu(d_det, n_det, alternative="less")
    else:
        p3 = np.nan

    # Test 4: mean VAF per carrier
    d_vaf = df_def[f"mean_vaf_{lab}"].dropna().values.astype(float)
    n_vaf = df_nor[f"mean_vaf_{lab}"].dropna().values.astype(float)
    if len(d_vaf) > 1 and len(n_vaf) > 1:
        _, p4 = mannwhitneyu(d_vaf, n_vaf, alternative="less")
    else:
        p4 = np.nan

    # Test 5: Fisher's exact on zero rate (are defective variants more likely undetected?)
    from scipy.stats import fisher_exact
    n_zeros_d = int((d_all == 0).sum())
    n_zeros_n = int((n_all == 0).sum())
    _, p5 = fisher_exact(
        [[n_zeros_d, len(d_all) - n_zeros_d],
         [n_zeros_n, len(n_all) - n_zeros_n]],
        alternative="greater")   # H1: defective have MORE zeros

    results[lab] = dict(
        d_all=d_all, n_all=n_all,
        med_d=med_d, med_n=med_n, depl=depl, p1=p1,
        rho=rho, p2=p2,
        d_det=d_det, n_det=n_det,
        med_det_d=np.median(d_det) if len(d_det) else np.nan,
        med_det_n=np.median(n_det) if len(n_det) else np.nan,
        p3=p3,
        d_vaf=d_vaf, n_vaf=n_vaf,
        med_vaf_d=np.nanmedian(d_vaf) if len(d_vaf) else np.nan,
        med_vaf_n=np.nanmedian(n_vaf) if len(n_vaf) else np.nan,
        p4=p4,
        n_zeros_d=n_zeros_d, n_zeros_n=n_zeros_n, p5=p5,
    )

# ── helpers ───────────────────────────────────────────────────────────────────

rng = np.random.default_rng(42)

def sig_stars(p):
    if np.isnan(p):  return "n.s."
    if p < 0.001:    return "***"
    if p < 0.01:     return "**"
    if p < 0.05:     return "*"
    return "n.s."

def p_fmt(p):
    if np.isnan(p):   return "n/a"
    if p < 0.001:     return f"{p:.2e}"
    return f"{p:.3f}"

# Y-axis ticks for log1p plots
LOG1P_TICKS    = [0, 1, 5, 10, 50, 100, 500, 1000, 5000, 10000, 50000, 100000, 300000]
LOG1P_TICKVALS = [np.log1p(t) for t in LOG1P_TICKS]

def set_log1p_yticks(ax, y_max_orig):
    ticks = [(v, lv) for v, lv in zip(LOG1P_TICKS, LOG1P_TICKVALS)
             if v <= y_max_orig * 1.05]
    if ticks:
        ax.set_yticks([lv for _, lv in ticks])
        ax.set_yticklabels([f"{v:,}" for v, _ in ticks], fontsize=7)

# ── figure layout ─────────────────────────────────────────────────────────────

fig = plt.figure(figsize=(19, 18))
gs  = gridspec.GridSpec(3, 3, figure=fig, hspace=0.52, wspace=0.35,
                        top=0.94, bottom=0.04, left=0.07, right=0.97)

fig.suptitle(
    "Purifying selection against incorporation-defective 5S rRNA variants  ·  UKBB n=490,075\n"
    f"284 annotated gene-body variants: {len(df_def)} defective (incorp<0.5) "
    f"vs {len(df_nor)} normal (incorp≥0.5)",
    fontsize=11)

# ── Row 0: boxplot + strip of carrier count ───────────────────────────────────

for col_i, (lab, T) in enumerate(THRESHOLDS):
    ax  = fig.add_subplot(gs[0, col_i])
    res = results[lab]

    y_max_orig = max(res["d_all"].max(), res["n_all"].max())
    y_max_t    = np.log1p(y_max_orig) * 1.12
    ax.set_ylim(-0.25, y_max_t)

    for arr, col_c, xc in [(res["d_all"], COL_DEF, 0), (res["n_all"], COL_NOR, 1)]:
        yt = np.log1p(arr)
        # jittered strip
        jit = rng.uniform(-0.22, 0.22, len(yt))
        ax.scatter(xc + jit, yt, s=6, alpha=0.30, c=col_c, zorder=2)
        # box
        q25, med, q75 = np.percentile(yt, [25, 50, 75])
        iqr = q75 - q25
        ax.add_patch(Rectangle((xc - 0.22, q25), 0.44, q75 - q25,
                                fc=col_c, alpha=0.22, zorder=3, lw=0))
        # whiskers (1.5 IQR)
        lo_w = max(float(yt.min()), q25 - 1.5*iqr)
        hi_w = min(float(yt.max()), q75 + 1.5*iqr)
        ax.plot([xc, xc], [lo_w, q25], c=col_c, lw=1.2, zorder=4)
        ax.plot([xc, xc], [q75, hi_w], c=col_c, lw=1.2, zorder=4)
        # median bar
        ax.plot([xc - 0.28, xc + 0.28], [med, med],
                color="black", lw=2.5, zorder=5)
        # median label just above bar (in original units)
        orig_med = int(np.round(np.expm1(med)))
        ax.text(xc, med + 0.06, f"{orig_med:,}",
                ha="center", va="bottom", fontsize=8, fontweight="bold")

    set_log1p_yticks(ax, y_max_orig)
    ax.set_xlim(-0.5, 1.5)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Defective\n(incorp<0.5)", "Normal\n(incorp≥0.5)"], fontsize=8.5)
    ax.set_ylabel("UKBB carriers  (log₁₊ₓ scale)" if col_i == 0 else "")
    ax.grid(True, lw=0.3, alpha=0.4, axis="y")

    # stat annotation — placed at fixed y position near top, no overlap
    stars = sig_stars(res["p1"])
    depl_s = f"{res['depl']:.0f}% depletion" if not np.isnan(res["depl"]) else "med = 0 both"
    fc = "#fff3b0" if stars != "n.s." else "#f0f0f0"
    ec = "goldenrod" if stars != "n.s." else "grey"
    ax.text(0.5, 0.99,
            f"{stars}  {depl_s}   MW p={p_fmt(res['p1'])}",
            transform=ax.transAxes, ha="center", va="top", fontsize=8,
            bbox=dict(boxstyle="round,pad=0.25", fc=fc, ec=ec, alpha=0.9))

    # zero-count annotation below x-axis
    n_zeros_d = int((res["d_all"] == 0).sum())
    n_zeros_n = int((res["n_all"] == 0).sum())
    if n_zeros_d + n_zeros_n > 0:
        ax.text(0,  -0.22, f"{n_zeros_d} zeros", ha="center", fontsize=7, color="grey")
        ax.text(1,  -0.22, f"{n_zeros_n} zeros", ha="center", fontsize=7, color="grey")

    panel = ["A","B","C"][col_i]
    ax.set_title(f"{panel}  Carrier count  —  VAF ≥ {lab}", fontsize=10, loc="left")

# ── Row 1 D&E: binned-median trend (incorp vs carrier count) ─────────────────

for col_i, (lab, T) in enumerate(THRESHOLDS):
    ax  = fig.add_subplot(gs[1, col_i])
    col = f"n_{lab}"
    res = results[lab]

    # Clip incorp to 0–2 for readability (1 outlier at 4.51)
    mask    = df["incorp"] <= 2.0
    df_clip = df[mask]
    n_excl  = int((~mask).sum())

    y_max_orig = float(df_clip[col].max())
    y_max_t    = np.log1p(y_max_orig) * 1.12
    ax.set_ylim(-0.15, y_max_t)

    # scatter: coloured by group
    for sub, col_c, s, alpha in [
            (df_clip[df_clip["is_defective"]],  COL_DEF, 15, 0.55),
            (df_clip[~df_clip["is_defective"]], COL_NOR,  8, 0.30)]:
        jit = rng.uniform(-0.01, 0.01, len(sub))
        ax.scatter(sub["incorp"] + jit,
                   np.log1p(sub[col].values.astype(float)),
                   s=s, alpha=alpha, c=col_c, zorder=2)

    # Binned median line (6 equal-freq bins across full incorp range, not clipped)
    incorp_full = df["incorp"].values
    q_edges     = np.quantile(incorp_full, np.linspace(0, 1, 7))
    bin_x, bin_med, bin_q25, bin_q75 = [], [], [], []
    for lo, hi in zip(q_edges[:-1], q_edges[1:]):
        mask_b = (incorp_full >= lo) & (incorp_full <= hi)
        vals   = np.log1p(df.loc[mask_b, col].values.astype(float))
        if vals.size >= 3:
            bin_x.append(float((lo + hi) / 2))
            bin_med.append(float(np.median(vals)))
            bin_q25.append(float(np.percentile(vals, 25)))
            bin_q75.append(float(np.percentile(vals, 75)))

    bx = np.array(bin_x); bm = np.array(bin_med)
    bq25 = np.array(bin_q25); bq75 = np.array(bin_q75)
    # Only plot within visible x range
    vis = bx <= 2.0
    ax.fill_between(bx[vis], bq25[vis], bq75[vis],
                    alpha=0.18, color="black", zorder=3)
    ax.plot(bx[vis], bm[vis], "o-", color="black", lw=2, ms=5,
            zorder=4, label="Binned median (IQR)")

    ax.axvline(INCORP_CUT, color="grey", lw=1, ls="--", alpha=0.7,
               label="Defective cutoff (0.5)")
    ax.set_xlim(0, 2.02)
    ax.set_xlabel("Incorporation efficiency (normalised to WT)")
    ax.set_ylabel("UKBB carriers  (log₁₊ₓ)" if col_i == 0 else "")
    ax.grid(True, lw=0.3, alpha=0.4)

    set_log1p_yticks(ax, y_max_orig)

    # legend for groups
    from matplotlib.lines import Line2D
    handles = [
        Line2D([0],[0], marker="o", color="w", markerfacecolor=COL_DEF, ms=7,
               label=f"Defective (n={int(df_clip['is_defective'].sum())})"),
        Line2D([0],[0], marker="o", color="w", markerfacecolor=COL_NOR, ms=5,
               label=f"Normal (n={int((~df_clip['is_defective']).sum())})"),
        Line2D([0],[0], color="black", lw=2, label="Binned median ± IQR"),
    ]
    ax.legend(handles=handles, fontsize=7.5, loc="upper left")

    if n_excl:
        ax.text(0.99, 0.02, f"({n_excl} outlier(s) with incorp>2 excluded from plot)",
                transform=ax.transAxes, ha="right", va="bottom",
                fontsize=6.5, color="grey", style="italic")

    # stat annotation — show zero-rate decomposition + conditional test
    rho, p_sp  = res["rho"], res["p2"]
    p_cond     = res["p3"]
    p_zero     = res["p5"]
    n_det_d    = len(res["d_det"]); n_det_n = len(res["n_det"])
    pct_d      = 100 * n_det_d / len(res["d_all"])
    pct_n      = 100 * n_det_n / len(res["n_all"])
    med_det_d  = res["med_det_d"]; med_det_n = res["med_det_n"]

    # choose annotation background based on strongest signal
    best_p = min(p for p in [p_sp, p_cond, p_zero] if not np.isnan(p))
    stars_best = sig_stars(best_p)
    fc = "#fff3b0" if stars_best != "n.s." else "#f0f0f0"
    ec = "goldenrod" if stars_best != "n.s." else "grey"

    lines = [
        f"Spearman ρ={rho:.3f} (p={p_fmt(p_sp)}) {sig_stars(p_sp)}",
        f"Zero rate: def {pct_d:.0f}% vs nor {pct_n:.0f}%",
        f"  Fisher p={p_fmt(p_zero)} {sig_stars(p_zero)}",
        f"Detected-only (n={n_det_d}/{n_det_n}):",
        f"  med={med_det_d:.0f} vs {med_det_n:.0f}  MW p={p_fmt(p_cond)} {sig_stars(p_cond)}",
    ]
    ax.text(0.99, 0.99, "\n".join(lines),
            transform=ax.transAxes, ha="right", va="top", fontsize=7.5,
            bbox=dict(boxstyle="round,pad=0.3", fc=fc, ec=ec, alpha=0.9))

    panel = ["D","E","F"][col_i]
    ax.set_title(f"{panel}  Incorp. efficiency vs carriers  —  VAF ≥ {lab}",
                 fontsize=10, loc="left")

# ── Row 1 F: statistical summary table ───────────────────────────────────────

ax_F = fig.add_subplot(gs[2, :])
ax_F.axis("off")
ax_F.set_title("G  Statistical summary", fontsize=10, loc="left", pad=8)

def row_color(p):
    if np.isnan(p): return "#f7f7f7"
    if p < 0.001:   return "#c8e6c9"   # green
    if p < 0.01:    return "#dcedc8"
    if p < 0.05:    return "#fff9c4"   # yellow
    return "#f7f7f7"                   # grey

# Build table rows
LABS = ["0.30%", "1.00%", "3.00%"]
col_labels = ["Test", "Metric", "0.30%", "1.00%", "3.00%"]

def _zr(lab):
    r = results[lab]
    return (f"def {r['n_zeros_d']}/{len(r['d_all'])} "
            f"nor {r['n_zeros_n']}/{len(r['n_all'])}")

table_data = [
    ["1. MW carrier count",
     "Median depletion",
     f"{results['0.30%']['depl']:.0f}%",
     f"{results['1.00%']['depl']:.0f}%",
     "med=0 both"],
    ["",
     "p-value",
     f"{p_fmt(results['0.30%']['p1'])} {sig_stars(results['0.30%']['p1'])}",
     f"{p_fmt(results['1.00%']['p1'])} {sig_stars(results['1.00%']['p1'])}",
     f"{p_fmt(results['3.00%']['p1'])} {sig_stars(results['3.00%']['p1'])}"],
    ["2. Spearman ρ",
     "ρ",
     f"{results['0.30%']['rho']:.3f}",
     f"{results['1.00%']['rho']:.3f}",
     f"{results['3.00%']['rho']:.3f}"],
    ["",
     "p-value",
     f"{p_fmt(results['0.30%']['p2'])} {sig_stars(results['0.30%']['p2'])}",
     f"{p_fmt(results['1.00%']['p2'])} {sig_stars(results['1.00%']['p2'])}",
     f"{p_fmt(results['3.00%']['p2'])} {sig_stars(results['3.00%']['p2'])}"],
    ["3. Conditional MW",
     "med def / nor",
     (f"{results['0.30%']['med_det_d']:.0f} / {results['0.30%']['med_det_n']:.0f}  "
      f"(n={len(results['0.30%']['d_det'])}/{len(results['0.30%']['n_det'])})"),
     (f"{results['1.00%']['med_det_d']:.0f} / {results['1.00%']['med_det_n']:.0f}  "
      f"(n={len(results['1.00%']['d_det'])}/{len(results['1.00%']['n_det'])})"),
     (f"{results['3.00%']['med_det_d']:.0f} / {results['3.00%']['med_det_n']:.0f}  "
      f"(n={len(results['3.00%']['d_det'])}/{len(results['3.00%']['n_det'])})")],
    ["",
     "p-value",
     f"{p_fmt(results['0.30%']['p3'])} {sig_stars(results['0.30%']['p3'])}",
     f"{p_fmt(results['1.00%']['p3'])} {sig_stars(results['1.00%']['p3'])}",
     f"{p_fmt(results['3.00%']['p3'])} {sig_stars(results['3.00%']['p3'])}"],
    ["4. Fisher zero rate",
     "zeros def / nor",
     _zr("0.30%"), _zr("1.00%"), _zr("3.00%")],
    ["  (more zeros in def?)",
     "p-value",
     f"{p_fmt(results['0.30%']['p5'])} {sig_stars(results['0.30%']['p5'])}",
     f"{p_fmt(results['1.00%']['p5'])} {sig_stars(results['1.00%']['p5'])}",
     f"{p_fmt(results['3.00%']['p5'])} {sig_stars(results['3.00%']['p5'])}"],
    ["5. Mean VAF/carrier",
     "med def / nor",
     f"{results['0.30%']['med_vaf_d']*100:.2f}% / {results['0.30%']['med_vaf_n']*100:.2f}%",
     f"{results['1.00%']['med_vaf_d']*100:.2f}% / {results['1.00%']['med_vaf_n']*100:.2f}%",
     f"{results['3.00%']['med_vaf_d']*100:.2f}% / {results['3.00%']['med_vaf_n']*100:.2f}%"],
    ["",
     "p-value",
     f"{p_fmt(results['0.30%']['p4'])} {sig_stars(results['0.30%']['p4'])}",
     f"{p_fmt(results['1.00%']['p4'])} {sig_stars(results['1.00%']['p4'])}",
     f"{p_fmt(results['3.00%']['p4'])} {sig_stars(results['3.00%']['p4'])}"],
]

tbl = ax_F.table(
    cellText=table_data,
    colLabels=col_labels,
    loc="upper center", cellLoc="center",
    bbox=[0, 0.35, 0.62, 0.60])
tbl.auto_set_font_size(False)
tbl.set_fontsize(8.5)
tbl.scale(1, 1.6)

# Color p-value rows
p_row_map = {1: ("p1","p1","p1"), 3: ("p2","p2","p2"),
             5: ("p3","p3","p3"), 7: ("p5","p5","p5"), 9: ("p4","p4","p4")}
for row_i, keys in p_row_map.items():
    for col_j, (lab, key) in enumerate(zip(LABS, keys)):
        p_val = results[lab][key]
        c = row_color(p_val)
        for ci in range(5):
            tbl[row_i + 1, ci].set_facecolor(c)

# Header row
for ci in range(5):
    tbl[0, ci].set_facecolor("#e3f2fd")
    tbl[0, ci].set_text_props(fontweight="bold")

# Legend for the significance shading below the table
note = (
    "Cell shading reflects the test p-value for each VAF threshold:\n"
    "  green   p < 0.01\n"
    "  yellow  p < 0.05\n"
    "  grey    n.s.\n"
    "Tests compare incorporation-defective (incorp<0.5) vs normal (incorp>=0.5)\n"
    "gene-body variants on carrier count and mean VAF per carrier."
)
ax_F.text(0.65, 0.95, note, transform=ax_F.transAxes,
          va="top", fontsize=8, family="monospace",
          bbox=dict(boxstyle="round,pad=0.4", fc="#f5f5f5", ec="grey", alpha=0.9))

# ── save ─────────────────────────────────────────────────────────────────────

out = OUTDIR / "80_ukbb_purifying_selection.pdf"
plt.savefig(out, dpi=180, bbox_inches="tight")
plt.close()
print(f"Figure → {out}")
