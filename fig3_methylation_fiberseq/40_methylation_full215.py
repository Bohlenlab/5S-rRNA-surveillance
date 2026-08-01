#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 40_methylation_full215.py — Per-copy CpG methylation figures across the 5S array for the full 215-proband cohort (ONT), from the SQLite database.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
40_methylation_full215.py

Per-copy CpG methylation figures across the 5S array for the full 215-proband
cohort (HPRC Year 1 + Release 2), from the SQLite database (ONT).

Data source: 5S_rDNA.db  ·  copy_methylation table
  — per-copy regional summaries: n_conf_calls, n_meth, mean_meth,
    nts_pre_n/meth, gene_n/meth, nts_post_n/meth
  — joined to copy (variants), haplotype, assembly (superpopulation)

Within-unit profiles are summarised at 3-region resolution
(NTS-pre / gene / NTS-post); between-copy and between-sample comparisons
use the per-copy regional means.

Array position (pct_pos) is computed from copy rank in ascending
unit_start_local order per haplotype; no orientation flip is applied.

Input : 5S_rDNA.db (tables copy_methylation, copy, variant, haplotype, assembly).
Output: <FIVES_OUT>/03_methylation_full215/{figname}.pdf
Data:   <FIVES_OUT>/03_methylation_full215/data/{figname}.tsv  (wide format)

Paths are read from environment variables (see repository README):
    FIVES_DB   path to 5S_rDNA.db
    FIVES_OUT  output directory
"""

import os
import sqlite3
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from scipy import stats
from scipy.ndimage import gaussian_filter1d
import statsmodels.formula.api as smf
import warnings
warnings.filterwarnings("ignore")
from pathlib import Path

# ── paths ──────────────────────────────────────────────────────────────────────
DB_PATH = Path(os.environ.get("FIVES_DB", "5S_rDNA.db"))
OUTDIR  = Path(os.environ.get("FIVES_OUT", "output")) / "03_methylation_full215"
DATADIR = OUTDIR / "data"
OUTDIR.mkdir(exist_ok=True)
DATADIR.mkdir(exist_ok=True)

# ── constants ──────────────────────────────────────────────────────────────────
MIN_CALLS   = 10
METH_CUTOFF = 0.50   # methylated vs unmethylated
HYPO_THR    = 65.0   # hypomethylated threshold (%)
HI_CUT      = 0.65   # copy class High
LO_CUT      = 0.35   # copy class Low

GENE_START, GENE_END = 630, 748   # within 2168 bp unit

POP_COLORS = {
    "AFR": "#e6194b", "AMR": "#3cb44b", "EAS": "#4363d8",
    "SAS": "#f58231", "EUR": "#911eb4", "unknown": "#808080",
}
HAP_COLORS = {"hap1": "#1f77b4", "hap2": "#ff7f0e"}
REGION_COLORS = {"NTS-pre": "#f4a582", "gene": "#aec6cf", "NTS (excl. ALU)": "#92c5de"}
pct_fmt = plt.FuncFormatter(lambda v, _: f"{v:.0f}%")

COHORTS = ("HPRC_Year1", "HPRC_Release2")


# ═══════════════════════════════════════════════════════════════════════════════
# 1. LOAD DATA
# ═══════════════════════════════════════════════════════════════════════════════

print("Loading data from SQLite …", flush=True)

con = sqlite3.connect(DB_PATH)

df = pd.read_sql_query(f"""
    SELECT
        cm.copy_id,
        cm.n_conf_calls, cm.n_meth, cm.mean_meth,
        cm.nts_pre_n, cm.nts_pre_meth,
        cm.gene_n,    cm.gene_meth,
        cm.nts_post_n,cm.nts_post_meth,
        COALESCE(cm.alu_n,   0) AS alu_n,
        COALESCE(cm.alu_meth,0) AS alu_meth,
        c.unit_start_local,
        c.n_snv_gene, c.n_snv_5s_gene, c.n_snv_nts_pre, c.n_snv_nts_post,
        h.hap_label, h.haplotype_id,
        a.sample_id, a.superpopulation, a.cohort
    FROM copy_methylation cm
    JOIN copy c ON cm.copy_id = c.copy_id
    JOIN haplotype h ON c.haplotype_id = h.haplotype_id
    JOIN assembly a ON h.assembly_id = a.assembly_id
    WHERE c.border_note = 'interior'
      AND cm.n_conf_calls >= {MIN_CALLS}
      AND a.cohort IN ({','.join(repr(c) for c in COHORTS)})
""", con)

# also load all interior copy counts per haplotype (for copy-number axis)
hap_ncopy = pd.read_sql_query(f"""
    SELECT h.haplotype_id, COUNT(*) as n_copies_interior
    FROM copy c
    JOIN haplotype h ON c.haplotype_id = h.haplotype_id
    JOIN assembly a ON h.assembly_id = a.assembly_id
    WHERE c.border_note = 'interior'
      AND a.cohort IN ({','.join(repr(c) for c in COHORTS)})
    GROUP BY h.haplotype_id
""", con)

alu_snv = pd.read_sql_query("""
    SELECT copy_id, COUNT(*) AS n_snv_alu
    FROM variant
    WHERE consensus_pos >= 787 AND consensus_pos < 1066 AND masked = 0
    GROUP BY copy_id
""", con)
con.close()

df = df.merge(alu_snv, on="copy_id", how="left")
df["n_snv_alu"] = df["n_snv_alu"].fillna(0).astype(int)

# ── derived per-copy columns ──────────────────────────────────────────────────
df["meth_pct"] = df["mean_meth"] * 100.0

def safe_frac(meth, n):
    return np.where(n > 0, meth / n * 100.0, np.nan)

df["nts_pre_pct"]  = safe_frac(df["nts_pre_meth"],  df["nts_pre_n"])
df["gene_pct"]     = safe_frac(df["gene_meth"],      df["gene_n"])
df["nts_post_pct"] = safe_frac(df["nts_post_meth"],  df["nts_post_n"])
df["alu_pct"]      = safe_frac(df["alu_meth"],       df["alu_n"])
# "other NTS" = nts_pre + (nts_post minus ALU slice) — excludes gene and ALU
df["other_nts_n"]     = df["nts_pre_n"]    + (df["nts_post_n"]    - df["alu_n"]).clip(lower=0)
df["other_nts_meth"]  = df["nts_pre_meth"] + (df["nts_post_meth"] - df["alu_meth"]).clip(lower=0)
df["other_nts_pct"]   = safe_frac(df["other_nts_meth"], df["other_nts_n"])
df["n_snv_other_nts"] = (df["n_snv_nts_pre"] + (df["n_snv_nts_post"] - df["n_snv_alu"])
                         .clip(lower=0)).astype(int)

df["is_methylated"] = (df["meth_pct"] > METH_CUTOFF * 100).astype(int)
df["is_hypo"]       = (df["meth_pct"] < HYPO_THR).astype(int)

def classify(m):
    if m >= HI_CUT * 100: return "High"
    if m >= LO_CUT * 100: return "Intermediate"
    return "Low"

df["copy_class"] = df["meth_pct"].apply(classify)

# per-haplotype mean methylation (for delta_meth)
hap_mean = df.groupby("haplotype_id")["meth_pct"].transform("mean")
df["delta_meth"] = df["meth_pct"] - hap_mean

# array position (% of array length, based on copy rank in sorted order)
def pct_positions(group):
    ranked = group["unit_start_local"].rank(method="first") - 1
    n = len(ranked)
    return ranked / max(n - 1, 1) * 100.0

df["pct_pos"] = df.groupby("haplotype_id", group_keys=False).apply(pct_positions)

# join haplotype copy counts
df = df.merge(hap_ncopy, on="haplotype_id", how="left")

# superpopulation fallback
df["superpopulation"] = df["superpopulation"].fillna("unknown")

pops_present = sorted(df["superpopulation"].unique())

print(f"  {df['sample_id'].nunique()} individuals, "
      f"{df.groupby(['sample_id','hap_label']).ngroups} haplotypes, "
      f"{len(df):,} copies with ≥{MIN_CALLS} confident calls", flush=True)


# ═══════════════════════════════════════════════════════════════════════════════
# HELPERS
# ═══════════════════════════════════════════════════════════════════════════════

def save_fig(fig, name):
    path = OUTDIR / f"{name}.pdf"
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {path.name}", flush=True)


def save_data(df_out, name, **kwargs):
    path = DATADIR / f"{name}.tsv"
    df_out.to_csv(path, sep="\t", index=False, **kwargs)


def cohort_title():
    n_y1 = df[df["cohort"] == "HPRC_Year1"]["sample_id"].nunique()
    n_r2 = df[df["cohort"] == "HPRC_Release2"]["sample_id"].nunique()
    return (f"Full cohort: {df['sample_id'].nunique()} probands  "
            f"({n_y1} Year-1 + {n_r2} Release-2)")


def pop_legend_handles():
    return [mpatches.Patch(color=POP_COLORS.get(p, "#808080"), label=p)
            for p in pops_present if p in POP_COLORS]


def add_region_shading(ax, alpha=0.30):
    """Light gene-body shading on a within-unit position axis."""
    ax.axvspan(GENE_START, GENE_END, color="#aec6cf", alpha=alpha, zorder=0)


def t_ci(sample_means, alpha=0.95):
    """Mean ± t-distribution 95% CI for an array of sample means."""
    finite = np.asarray(sample_means, dtype=float)
    finite = finite[np.isfinite(finite)]
    n = len(finite)
    if n < 3:
        m = float(np.nanmean(sample_means)) if len(sample_means) else np.nan
        return m, np.nan, np.nan
    m  = float(np.mean(finite))
    se = float(np.std(finite, ddof=1) / np.sqrt(n))
    t  = float(stats.t.ppf((1 + alpha) / 2, df=n - 1))
    return m, m - t * se, m + t * se


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 16 — METHYLATION PROFILES (regional resolution)
# ═══════════════════════════════════════════════════════════════════════════════

print("Fig 16 …", flush=True)

# ── per-sample, per-haplotype, per-region means ────────────────────────────────
reg_long_rows = []
for region, col in [("NTS-pre", "nts_pre_pct"), ("gene", "gene_pct"), ("NTS (excl. ALU)", "other_nts_pct")]:
    tmp = (df.dropna(subset=[col])
             .groupby(["sample_id", "hap_label", "superpopulation"])[col]
             .mean().reset_index()
             .rename(columns={col: "meth_pct"}))
    tmp["region"] = region
    reg_long_rows.append(tmp)
reg_long = pd.concat(reg_long_rows, ignore_index=True)

# per-sample (both haps averaged) per-region mean
reg_sample = (reg_long.groupby(["sample_id", "superpopulation", "region"])
              ["meth_pct"].mean().reset_index())

# per-sample per-hap overall mean
sample_mean = (df.groupby(["sample_id", "hap_label", "superpopulation"])
               ["meth_pct"].mean().reset_index())

# wide for hap1 vs hap2
hap_wide = (sample_mean.groupby(["sample_id", "superpopulation", "hap_label"])
            ["meth_pct"].mean()
            .unstack("hap_label").reset_index()
            .dropna(subset=["hap1", "hap2"]))

# per-copy means
copy_mean_df = df[["sample_id", "hap_label", "superpopulation", "copy_id", "meth_pct"]].copy()

# ── figure ─────────────────────────────────────────────────────────────────────
REGION_ORDER = ["NTS-pre", "gene", "NTS (excl. ALU)"]
REGION_X     = {"NTS-pre": 0, "gene": 1, "NTS (excl. ALU)": 2}

fig, axes = plt.subplots(2, 2, figsize=(14, 11))
fig.suptitle(f"5S rRNA CpG methylation — {cohort_title()}\n"
             "Regional summaries (NTS-pre / gene / NTS excl. ALU) from SQLite copy_methylation table",
             fontsize=9, y=1.01)

# Panel A: per-sample 3-region profiles
ax = axes[0, 0]
pop_traces = {p: {r: [] for r in REGION_ORDER} for p in pops_present}
for (sid, hap, spop), grp in reg_long.groupby(["sample_id", "hap_label", "superpopulation"]):
    col = POP_COLORS.get(spop, "#808080")
    xs = [REGION_X[r] for r in REGION_ORDER if r in grp["region"].values]
    ys = [grp.loc[grp["region"] == r, "meth_pct"].values[0]
          for r in REGION_ORDER if r in grp["region"].values]
    if len(xs) == 3:
        ax.plot(xs, ys, color=col, lw=0.6, alpha=0.18, marker="o", ms=2, zorder=2)
        for r, y in zip(REGION_ORDER, ys):
            pop_traces[spop][r].append(y)
for pop in pops_present:
    vals = [np.nanmean(pop_traces[pop][r]) for r in REGION_ORDER
            if len(pop_traces[pop][r]) > 0]
    if len(vals) == 3:
        ax.plot(range(3), vals, color=POP_COLORS.get(pop, "#808080"),
                lw=2.2, marker="o", ms=5, alpha=0.95, zorder=4, label=pop)
ax.set_xticks([0, 1, 2])
ax.set_xticklabels(REGION_ORDER, fontsize=9)
ax.set_ylabel("Mean CpG methylation (%)", fontsize=8)
ax.set_title("(A) Per-sample region profiles\n(thin = individual haps, thick = pop. mean)",
             fontsize=9)
ax.yaxis.set_major_formatter(pct_fmt)
ax.set_xlim(-0.3, 2.3)
ax.tick_params(labelsize=7)
ax.grid(axis="y", lw=0.3, alpha=0.4)
ax.legend(handles=pop_legend_handles(), fontsize=7, loc="best", framealpha=0.85)

# Panel B: per-region boxplot
ax = axes[0, 1]
for ri, region in enumerate(REGION_ORDER):
    sub = reg_sample[reg_sample["region"] == region]
    bp = ax.boxplot(sub["meth_pct"].values, positions=[ri * 1.5], widths=0.55,
                    patch_artist=True, notch=False,
                    medianprops=dict(color="black", lw=1.5),
                    boxprops=dict(facecolor=REGION_COLORS[region], alpha=0.5),
                    whiskerprops=dict(lw=0.8), flierprops=dict(marker=""),
                    showfliers=False)
    np.random.seed(42 + ri)
    jitter = np.random.uniform(-0.18, 0.18, len(sub))
    for i, row in enumerate(sub.itertuples()):
        ax.scatter(ri * 1.5 + jitter[i], row.meth_pct,
                   s=18, color=POP_COLORS.get(row.superpopulation, "#808080"),
                   alpha=0.7, zorder=4, edgecolors="white", linewidths=0.3)
ax.set_xticks([0, 1.5, 3.0])
ax.set_xticklabels(REGION_ORDER, fontsize=9)
ax.set_ylabel("Mean methylation (%)", fontsize=8)
ax.set_title("(B) Per-region methylation\n(each dot = one proband, haps averaged)",
             fontsize=9)
ax.yaxis.set_major_formatter(pct_fmt)
ax.tick_params(labelsize=7)
ax.set_xlim(-0.6, 3.6)
ax.grid(axis="y", lw=0.3, alpha=0.4)
ax.legend(handles=pop_legend_handles(), fontsize=7, framealpha=0.85)

# Panel C: hap1 vs hap2 scatter
ax = axes[1, 0]
ax.plot([0, 100], [0, 100], "k--", lw=0.8, alpha=0.4, zorder=1)
for _, row in hap_wide.iterrows():
    ax.scatter(row["hap1"], row["hap2"],
               s=30, color=POP_COLORS.get(row["superpopulation"], "#808080"),
               alpha=0.8, zorder=3, edgecolors="white", linewidths=0.4)
r_hap = np.corrcoef(hap_wide["hap1"], hap_wide["hap2"])[0, 1]
ax.text(0.97, 0.04, f"r = {r_hap:.2f}", ha="right", va="bottom",
        fontsize=7, transform=ax.transAxes, color="dimgrey")
ax.set_xlabel("Haplotype 1 mean methylation (%)", fontsize=8)
ax.set_ylabel("Haplotype 2 mean methylation (%)", fontsize=8)
ax.set_title("(C) Hap1 vs hap2 mean methylation per proband", fontsize=9)
ax.yaxis.set_major_formatter(pct_fmt)
ax.xaxis.set_major_formatter(pct_fmt)
ax.set_xlim(0, 105); ax.set_ylim(0, 105)
ax.tick_params(labelsize=7)
ax.grid(lw=0.3, alpha=0.35)
ax.legend(handles=pop_legend_handles(), fontsize=7, framealpha=0.85)

# Panel D: per-copy distribution per sample
ax = axes[1, 1]
sample_order = (hap_wide.sort_values(["superpopulation", "sample_id"])
                ["sample_id"].tolist())
y_pos = {sid: i for i, sid in enumerate(sample_order)}
ax.set_yticks(range(len(sample_order)))
ax.set_yticklabels(sample_order, fontsize=4.0)
prev_pop = None
for i, sid in enumerate(sample_order):
    spop = df.loc[df["sample_id"] == sid, "superpopulation"].iloc[0]
    if spop != prev_pop:
        if prev_pop is not None:
            ax.axhline(i - 0.5, color="grey", lw=0.6, ls="--", alpha=0.5)
        prev_pop = spop
np.random.seed(0)
for sid in sample_order:
    sub = copy_mean_df[copy_mean_df["sample_id"] == sid]
    spop = sub["superpopulation"].iloc[0]
    col  = POP_COLORS.get(spop, "#808080")
    y    = y_pos[sid]
    jitter = np.random.uniform(-0.35, 0.35, len(sub))
    ax.scatter(sub["meth_pct"].values, y + jitter,
               s=3, color=col, alpha=0.25, linewidths=0, zorder=2)
    q25, q75 = np.percentile(sub["meth_pct"], [25, 75])
    ax.plot([q25, q75], [y, y], color=col, lw=2, alpha=0.85, zorder=3)
ax.set_xlabel("Per-copy mean methylation (%)", fontsize=8)
ax.set_title("(D) Per-copy methylation distribution per proband\n"
             "(dots = copies; bar = IQR)", fontsize=9)
ax.xaxis.set_major_formatter(pct_fmt)
ax.set_xlim(0, 110)
ax.tick_params(labelsize=5.5)
ax.grid(axis="x", lw=0.3, alpha=0.35)

plt.tight_layout()
save_fig(fig, "16_methylation_profiles")

# ── data table 16 ─────────────────────────────────────────────────────────────
t16 = (reg_sample.pivot(index=["sample_id", "superpopulation"],
                        columns="region", values="meth_pct")
       .rename(columns={"NTS-pre": "NTS_pre_pct", "gene": "gene_pct",
                        "NTS-post": "NTS_post_pct"})
       .reset_index())
h1 = (sample_mean[sample_mean["hap_label"] == "hap1"]
      .rename(columns={"meth_pct": "hap1_overall_pct"})[["sample_id", "hap1_overall_pct"]])
h2 = (sample_mean[sample_mean["hap_label"] == "hap2"]
      .rename(columns={"meth_pct": "hap2_overall_pct"})[["sample_id", "hap2_overall_pct"]])
t16 = t16.merge(h1, on="sample_id", how="left").merge(h2, on="sample_id", how="left")
save_data(t16.round(2), "16_methylation_profiles")


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 17 — METHYLATION ALONG THE ARRAY
# ═══════════════════════════════════════════════════════════════════════════════

print("Fig 17 …", flush=True)

N_BINS       = 20
SMOOTH_SIGMA = 10

def smooth_trace(pct_arr, meth_arr, n_pts=200, sigma=SMOOTH_SIGMA):
    if len(pct_arr) < 3:
        return None, None
    idx = np.argsort(pct_arr)
    x   = np.linspace(0, 100, n_pts)
    y   = np.interp(x, pct_arr[idx], meth_arr[idx])
    return x, gaussian_filter1d(y, sigma=sigma)

def bin_stats(hap_df, bin_edges):
    mids, means, lo, hi = [], [], [], []
    for i in range(len(bin_edges) - 1):
        msk = (hap_df["pct_pos"] >= bin_edges[i]) & (hap_df["pct_pos"] < bin_edges[i+1])
        sub = hap_df[msk]
        mid = (bin_edges[i] + bin_edges[i+1]) / 2
        sm  = sub.groupby("sample_id")["meth_pct"].mean().values
        m, l, h = t_ci(sm)
        mids.append(mid); means.append(m); lo.append(l); hi.append(h)
    return np.array(mids), np.array(means), np.array(lo), np.array(hi)

bin_edges   = np.linspace(0, 100, N_BINS + 1)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

fig, axes = plt.subplots(1, 3, figsize=(18, 6))
fig.suptitle(
    f"5S array — CpG methylation along the array\n"
    f"{cohort_title()}  |  copy rank → 0–100% (no orientation correction)",
    fontsize=9, y=1.02
)

for panel_idx, hap in enumerate(("hap1", "hap2")):
    ax       = axes[panel_idx]
    hap_data = df[df["hap_label"] == hap]
    pop_traces = {p: [] for p in pops_present}
    for sid, grp in hap_data.groupby("sample_id"):
        spop  = grp["superpopulation"].iloc[0]
        x_sm, y_sm = smooth_trace(grp["pct_pos"].values, grp["meth_pct"].values)
        if x_sm is None:
            continue
        ax.plot(x_sm, y_sm, color=POP_COLORS.get(spop, "#808080"),
                lw=0.6, alpha=0.15, zorder=2)
        pop_traces[spop].append(y_sm)
    x_grid = np.linspace(0, 100, 200)
    for pop in pops_present:
        traces = np.array(pop_traces[pop])
        if len(traces) == 0:
            continue
        ax.plot(x_grid, np.nanmean(traces, axis=0),
                color=POP_COLORS.get(pop, "#808080"),
                lw=2.2, alpha=0.95, zorder=4, label=pop)
    ax.set_xlim(0, 100); ax.set_ylim(0, 100)
    ax.set_xlabel("Position along array (0% = first copy in DB order)", fontsize=8)
    ax.set_ylabel("Mean CpG methylation per copy (%)", fontsize=8)
    pl = "A" if hap == "hap1" else "B"
    ax.set_title(f"({pl}) {hap}  |  {hap_data['sample_id'].nunique()} probands",
                 fontsize=9)
    ax.yaxis.set_major_formatter(pct_fmt)
    ax.xaxis.set_major_formatter(pct_fmt)
    ax.tick_params(labelsize=7)
    ax.grid(lw=0.3, alpha=0.4)
    ax.legend(handles=pop_legend_handles(), fontsize=7, loc="lower right",
              framealpha=0.85)

ax = axes[2]
# panel C: pool hap1+hap2, show per-copy quantile distribution
PCTILES = [10, 25, 50, 75, 90]
rows17c = []
mids_c = []; q10s = []; q25s = []; q50s = []; q75s = []; q90s = []
means_c = []; ns_c = []
for i in range(len(bin_edges) - 1):
    msk = (df["pct_pos"] >= bin_edges[i]) & (df["pct_pos"] < bin_edges[i+1])
    sub = df.loc[msk, "meth_pct"].dropna().values
    if len(sub) < 5:
        continue
    mid = (bin_edges[i] + bin_edges[i+1]) / 2
    q10, q25, q50, q75, q90 = np.percentile(sub, PCTILES)
    mids_c.append(mid); q10s.append(q10); q25s.append(q25); q50s.append(q50)
    q75s.append(q75); q90s.append(q90)
    means_c.append(float(sub.mean())); ns_c.append(len(sub))
    rows17c.append({"bin_center_pct": round(mid, 1), "n_copies": len(sub),
                    "mean": round(float(sub.mean()), 2),
                    "q10": round(q10, 2), "q25": round(q25, 2),
                    "q50_median": round(q50, 2),
                    "q75": round(q75, 2), "q90": round(q90, 2)})
mids_c = np.array(mids_c); q10s = np.array(q10s); q25s = np.array(q25s)
q50s = np.array(q50s); q75s = np.array(q75s); q90s = np.array(q90s)
ax.fill_between(mids_c, q10s, q90s, color="#4c72b0", alpha=0.14,
                label="10th–90th pctile")
ax.fill_between(mids_c, q25s, q75s, color="#4c72b0", alpha=0.35,
                label="IQR (25th–75th)")
ax.plot(mids_c, q50s, color="#4c72b0", lw=2.4, label="Median")
ax.set_xlim(0, 100); ax.set_ylim(0, 100)
ax.set_xlabel("Position along array (%)", fontsize=8)
ax.set_ylabel("CpG methylation per copy (%)", fontsize=8)
n_all_c = df["sample_id"].nunique()
ax.set_title(
    f"(C) hap1+hap2 pooled — per-copy distribution per {100//N_BINS}% bin\n"
    f"n = {n_all_c} probands, {sum(ns_c):,} copy-observations  |  median + IQR + 10–90th",
    fontsize=9)
ax.yaxis.set_major_formatter(pct_fmt)
ax.xaxis.set_major_formatter(pct_fmt)
ax.tick_params(labelsize=7)
ax.grid(lw=0.3, alpha=0.4)
ax.legend(fontsize=8, framealpha=0.85)

plt.tight_layout()
save_fig(fig, "17_methylation_along_array")

# ── data table 17: per-bin quantile summary (wide) ────────────────────────────
save_data(pd.DataFrame(rows17c), "17_methylation_along_array")


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 19 — COPY CLASSIFICATION + COPY NUMBER VS METHYLATION
# ═══════════════════════════════════════════════════════════════════════════════

print("Fig 19 …", flush=True)

# per-haplotype stats
hap_stats = (df.groupby(["sample_id", "hap_label", "superpopulation", "haplotype_id"])
             .agg(n_callable=("copy_id", "count"),
                  n_meth=("is_methylated", "sum"),
                  mean_meth=("meth_pct", "mean"),
                  n_copies_interior=("n_copies_interior", "first"))
             .reset_index())
hap_stats["n_unmeth"]  = hap_stats["n_callable"] - hap_stats["n_meth"]
hap_stats["pct_meth"]  = hap_stats["n_meth"] / hap_stats["n_callable"] * 100

wide19 = hap_stats.pivot(index=["sample_id", "superpopulation"],
                         columns="hap_label",
                         values=["pct_meth", "n_meth", "n_unmeth",
                                 "n_callable", "n_copies_interior",
                                 "mean_meth"])
wide19.columns = ["_".join(c) for c in wide19.columns]
wide19 = wide19.reset_index().dropna(subset=["pct_meth_hap1", "pct_meth_hap2"])
wide19["imbalance"] = (wide19["pct_meth_hap1"] - wide19["pct_meth_hap2"]).abs()
wide19["mean_pct_meth"] = (wide19["pct_meth_hap1"] + wide19["pct_meth_hap2"]) / 2
r19, p19 = stats.pearsonr(wide19["pct_meth_hap1"], wide19["pct_meth_hap2"])

fig = plt.figure(figsize=(16, 13))
gs  = GridSpec(2, 2, hspace=0.40, wspace=0.32, figure=fig)
ax_hist   = fig.add_subplot(gs[0, 0])
ax_bar    = fig.add_subplot(gs[0, 1])
ax_scat   = fig.add_subplot(gs[1, 0])
ax_imbal  = fig.add_subplot(gs[1, 1])

fig.suptitle(
    f"5S rRNA per-copy methylation classification  (cutoff > {METH_CUTOFF*100:.0f}%)\n"
    f"{cohort_title()}  ·  {len(wide19)} probands with both haplotypes  ·  "
    f"{len(df):,} callable interior copies",
    fontsize=10, y=1.01
)

# A: distribution + cutoff
ax = ax_hist
ax.hist(df["meth_pct"], bins=50, color="#4d94c9", alpha=0.8, edgecolor="none")
ax.axvline(METH_CUTOFF * 100, color="#e6194b", lw=1.8, ls="--",
           label=f"cutoff {METH_CUTOFF*100:.0f}%")
pct_below = 100 * (df["meth_pct"] <= METH_CUTOFF * 100).mean()
ax.text(15, ax.get_ylim()[1] * 0.85,
        f"≤50%: {pct_below:.1f}%\n(unmethylated)",
        fontsize=8, ha="center")
ax.text(80, ax.get_ylim()[1] * 0.85,
        f">50%: {100-pct_below:.1f}%\n(methylated)",
        fontsize=8, ha="center")
ax.set_xlabel("Mean CpG methylation per copy (%)", fontsize=9)
ax.set_ylabel("Number of copies", fontsize=9)
ax.set_title(f"(A) Per-copy methylation  (n={len(df):,} copies, "
             f"≥{MIN_CALLS} calls)", fontsize=9)
ax.xaxis.set_major_formatter(pct_fmt)
ax.legend(fontsize=8)
ax.grid(axis="y", lw=0.3, alpha=0.4)

# B: stacked bar per sample
ax = ax_bar
sample_order_19 = wide19.sort_values("mean_pct_meth")["sample_id"].values
x_pos = np.arange(len(sample_order_19)) * 2.2
for i, sid in enumerate(sample_order_19):
    row = wide19[wide19["sample_id"] == sid].iloc[0]
    spop  = row["superpopulation"]
    color = POP_COLORS.get(spop, "#808080")
    for j, hap in enumerate(["hap1", "hap2"]):
        xc = x_pos[i] + j * 0.9
        nm  = row.get(f"n_meth_{hap}", 0) or 0
        nu  = row.get(f"n_unmeth_{hap}", 0) or 0
        if nm + nu == 0:
            continue
        ax.bar(xc, nm, width=0.8, color=color, alpha=0.85, bottom=0)
        ax.bar(xc, nu, width=0.8, color=color, alpha=0.28,
               bottom=nm, hatch="///", linewidth=0)
ax.set_xticks(x_pos + 0.45)
ax.set_xticklabels(sample_order_19, rotation=90, fontsize=3.8, ha="right")
ax.set_ylabel("Callable interior copies", fontsize=9)
ax.set_title("(B) Methylated (solid) / unmethylated (hatched) per haplotype\n"
             "hap1=left, hap2=right; sorted by ascending mean methylation",
             fontsize=8.5)
ax.tick_params(labelsize=7)
ax.grid(axis="y", lw=0.3, alpha=0.4)
ax.legend(handles=pop_legend_handles(), fontsize=7, loc="upper left",
          framealpha=0.85)

# C: hap1 vs hap2 % methylated copies
ax = ax_scat
for _, row in wide19.iterrows():
    ax.scatter(row["pct_meth_hap1"], row["pct_meth_hap2"],
               color=POP_COLORS.get(row["superpopulation"], "#808080"),
               s=30, alpha=0.8, linewidths=0, zorder=3)
ax.plot([0, 100], [0, 100], color="grey", lw=0.8, ls="--", alpha=0.5)
ax.set_xlim(0, 105); ax.set_ylim(0, 105)
ax.set_xlabel("% methylated copies — hap1", fontsize=9)
ax.set_ylabel("% methylated copies — hap2", fontsize=9)
ax.set_title(f"(C) Hap1 vs. hap2 methylated copy fraction\n"
             f"Pearson r = {r19:.3f}, p = {p19:.3g}  (n={len(wide19)} probands)",
             fontsize=9)
ax.xaxis.set_major_formatter(pct_fmt)
ax.yaxis.set_major_formatter(pct_fmt)
ax.tick_params(labelsize=8)
ax.grid(lw=0.3, alpha=0.4)
ax.legend(handles=pop_legend_handles(), fontsize=7, framealpha=0.85)

# D: imbalance distribution
ax = ax_imbal
ax.hist(wide19["imbalance"], bins=25, color="#7d5ba6", alpha=0.8, edgecolor="none")
ax.axvline(wide19["imbalance"].median(), color="#f58231", lw=1.8, ls="--",
           label=f"median {wide19['imbalance'].median():.1f}%")
ax.axvline(20, color="#e6194b", lw=1.2, ls=":", alpha=0.8, label="20% reference")
ax.set_xlabel("|hap1 % methylated − hap2 % methylated|", fontsize=9)
ax.set_ylabel("Number of probands", fontsize=9)
ax.set_title(f"(D) Inter-haplotype imbalance\n"
             f"(>20%: {(wide19['imbalance'] > 20).sum()} / {len(wide19)} probands)",
             fontsize=9)
ax.xaxis.set_major_formatter(pct_fmt)
ax.legend(fontsize=8)
ax.grid(axis="y", lw=0.3, alpha=0.4)

plt.tight_layout()
save_fig(fig, "19_methylation_copy_classification")

# companion: copy number vs methylation
print("Fig 19b …", flush=True)
fig2, ax2 = plt.subplots(1, 1, figsize=(7, 5))
for _, row in hap_stats.iterrows():
    col = POP_COLORS.get(row["superpopulation"], "#808080")
    hap_style = "o" if row["hap_label"] == "hap1" else "^"
    ax2.scatter(row["n_copies_interior"], row["mean_meth"],
                color=col, s=18, alpha=0.55, marker=hap_style,
                linewidths=0, zorder=3)
r_cn, p_cn = stats.pearsonr(hap_stats["n_copies_interior"].dropna(),
                             hap_stats.loc[hap_stats["n_copies_interior"].notna(),
                                           "mean_meth"])
m_cn, b_cn, *_ = stats.linregress(hap_stats["n_copies_interior"].dropna(),
                                    hap_stats.loc[hap_stats["n_copies_interior"].notna(),
                                                  "mean_meth"])
xs2 = np.array([hap_stats["n_copies_interior"].min(),
                hap_stats["n_copies_interior"].max()])
ax2.plot(xs2, m_cn * xs2 + b_cn, color="black", lw=1.5, ls="--", alpha=0.7,
         label=f"r={r_cn:.2f}, p={p_cn:.2g}")
ax2.set_xlabel("Interior copy number per haplotype", fontsize=9)
ax2.set_ylabel("Mean CpG methylation (%)", fontsize=9)
ax2.set_title(
    f"Copy number vs. mean methylation per haplotype\n"
    f"{cohort_title()}  |  ○=hap1, △=hap2",
    fontsize=9)
ax2.yaxis.set_major_formatter(pct_fmt)
ax2.tick_params(labelsize=8)
ax2.grid(lw=0.3, alpha=0.4)
pop_handles2 = pop_legend_handles()
pop_handles2 += [Line2D([0], [0], marker="o", color="k", ls="", ms=6, label="hap1"),
                 Line2D([0], [0], marker="^", color="k", ls="", ms=6, label="hap2")]
ax2.legend(handles=pop_handles2, fontsize=7, framealpha=0.85, ncol=2)
plt.tight_layout()
save_fig(fig2, "19_copynumber_vs_methylation")

# ── data tables 19 ────────────────────────────────────────────────────────────
save_data(wide19.round(2), "19_methylation_copy_classification")
save_data(hap_stats[["sample_id", "hap_label", "superpopulation",
                      "n_copies_interior", "mean_meth",
                      "n_callable", "n_meth", "n_unmeth", "pct_meth"]].round(2),
          "19_copynumber_vs_methylation")


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 20 — METHYLATION ALONG ARRAY BY REGION
# ═══════════════════════════════════════════════════════════════════════════════

print("Fig 20 …", flush=True)

REGIONS20 = [
    ("NTS-pre",         "nts_pre_pct"),
    ("gene",            "gene_pct"),
    ("NTS (excl. ALU)", "other_nts_pct"),
]

fig, axes = plt.subplots(3, 3, figsize=(18, 14))
fig.suptitle(
    f"Methylation along array — by region\n"
    f"{cohort_title()}  |  no orientation correction",
    fontsize=10, y=1.01
)

rows20 = []
for ri, (region, col) in enumerate(REGIONS20):
    region_data = df.dropna(subset=[col]).copy()
    region_data["region_meth"] = region_data[col]

    for panel_col, hap in enumerate(["hap1", "hap2"]):
        ax = axes[ri, panel_col]
        hap_data = region_data[region_data["hap_label"] == hap]
        pop_traces = {p: [] for p in pops_present}
        for sid, grp in hap_data.groupby("sample_id"):
            spop = grp["superpopulation"].iloc[0]
            x_sm, y_sm = smooth_trace(grp["pct_pos"].values,
                                       grp["region_meth"].values)
            if x_sm is None:
                continue
            ax.plot(x_sm, y_sm, color=POP_COLORS.get(spop, "#808080"),
                    lw=0.6, alpha=0.15)
            pop_traces[spop].append(y_sm)
        x_grid = np.linspace(0, 100, 200)
        for pop in pops_present:
            traces = np.array(pop_traces[pop])
            if len(traces) == 0:
                continue
            ax.plot(x_grid, np.nanmean(traces, axis=0),
                    color=POP_COLORS.get(pop, "#808080"), lw=2.2, label=pop)
        ax.set_xlim(0, 100); ax.set_ylim(0, 100)
        ax.xaxis.set_major_formatter(pct_fmt)
        ax.yaxis.set_major_formatter(pct_fmt)
        pl = chr(65 + ri * 3 + panel_col)
        ax.set_title(f"({pl}) {region} / {hap}", fontsize=9)
        ax.set_ylabel("Meth (%)", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(lw=0.3, alpha=0.4)
        if ri == 0 and panel_col == 0:
            ax.legend(handles=pop_legend_handles(), fontsize=7, framealpha=0.85)

        # bin stats for table
        mids, means, lo, hi = bin_stats(
            hap_data.rename(columns={"region_meth": "meth_pct"}), bin_edges)
        for mid, m, l, h in zip(mids, means, lo, hi):
            rows20.append({"region": region, "hap": hap,
                            "bin_center_pct": mid,
                            "mean_meth_pct": m, "ci_lo": l, "ci_hi": h})

    # summary panel (col 2)
    ax = axes[ri, 2]
    for hap in ["hap1", "hap2"]:
        hap_data = region_data[region_data["hap_label"] == hap].rename(
            columns={"region_meth": "meth_pct"})
        mids, means, lo, hi = bin_stats(hap_data, bin_edges)
        color = HAP_COLORS[hap]
        valid = ~np.isnan(means)
        ax.plot(mids[valid], means[valid], color=color, lw=2.2, label=f"{hap}")
        ax.fill_between(mids[valid], lo[valid], hi[valid],
                        color=color, alpha=0.20)
    ax.set_xlim(0, 100); ax.set_ylim(0, 100)
    ax.xaxis.set_major_formatter(pct_fmt)
    ax.yaxis.set_major_formatter(pct_fmt)
    pl = chr(65 + ri * 3 + 2)
    ax.set_title(f"({pl}) {region} — summary ± 95% CI", fontsize=9)
    ax.tick_params(labelsize=7)
    ax.grid(lw=0.3, alpha=0.4)
    ax.legend(fontsize=8, framealpha=0.85)
    if ri == 2:
        ax.set_xlabel("Position along array (%)", fontsize=8)

plt.tight_layout()
save_fig(fig, "20_methylation_along_array_by_region")
save_data(pd.DataFrame(rows20).round(2), "20_methylation_along_array_by_region")


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 21 — VARIANTS VS METHYLATION
# ═══════════════════════════════════════════════════════════════════════════════

print("Fig 21 …", flush=True)

GENE_BINS = [(-1, 0, "0"), (0, 1, "1"), (1, 2, "2"), (2, 200, "≥3")]

def gene_bin_label(n):
    if n == 0: return "0"
    if n == 1: return "1"
    if n == 2: return "2"
    return "≥3"

df["gene_snv_bin"] = df["n_snv_5s_gene"].apply(gene_bin_label)
GENE_BIN_ORDER = ["0", "1", "2", "≥3"]

fig, axes = plt.subplots(2, 3, figsize=(16, 10))
fig.suptitle(
    f"Sequence variants vs. CpG methylation\n{cohort_title()}",
    fontsize=10, y=1.01
)

def decile_boxplot(ax, groups, data_col, order, color, ylabel):
    """Box-and-whisker with 10th/90th percentile whiskers, median bar."""
    positions = range(len(order))
    for i, label in enumerate(order):
        sub = groups.get(label, pd.Series(dtype=float))
        if len(sub) == 0:
            continue
        q10, q25, q50, q75, q90 = np.percentile(sub, [10, 25, 50, 75, 90])
        ax.bar(i, q75 - q25, bottom=q25, width=0.5, color=color, alpha=0.5)
        ax.plot([i - 0.25, i + 0.25], [q50, q50], color="black", lw=2)
        ax.plot([i, i], [q10, q25], color="grey", lw=1)
        ax.plot([i, i], [q75, q90], color="grey", lw=1)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, fontsize=9)
    ax.yaxis.set_major_formatter(pct_fmt)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.grid(axis="y", lw=0.3, alpha=0.4)

# A: gene snv groups vs raw meth
ax = axes[0, 0]
grps_raw = {lb: df.loc[df["gene_snv_bin"] == lb, "meth_pct"] for lb in GENE_BIN_ORDER}
decile_boxplot(ax, grps_raw, "meth_pct", GENE_BIN_ORDER, "#4d94c9", "Raw methylation (%)")
ax.set_title("(A) Gene SNV groups vs raw methylation\n"
             "(box=IQR, line=median, whiskers=10/90th)", fontsize=8.5)
ax.set_xlabel("n_snv_5s_gene", fontsize=8)

# B: gene snv groups vs delta_meth
ax = axes[0, 1]
grps_delta = {lb: df.loc[df["gene_snv_bin"] == lb, "delta_meth"] for lb in GENE_BIN_ORDER}
decile_boxplot(ax, grps_delta, "delta_meth", GENE_BIN_ORDER, "#e07b3a",
               "Δ methylation (copy − hap mean, %)")
ax.set_title("(B) Gene SNV groups vs haplotype-centred Δmeth\n"
             "(removes inter-haplotype confounder)", fontsize=8.5)
ax.set_xlabel("n_snv_5s_gene", fontsize=8)
ax.axhline(0, color="grey", lw=0.8, ls="--", alpha=0.5)

# C: total variants vs raw meth scatter
ax = axes[0, 2]
total_snv = df["n_snv_gene"].clip(0, 20)
ax.scatter(total_snv + np.random.default_rng(0).uniform(-0.3, 0.3, len(df)),
           df["meth_pct"], s=3, alpha=0.12, color="#4d94c9", linewidths=0)
slope_raw, b_raw, r_raw, p_raw, _ = stats.linregress(df["n_snv_gene"], df["meth_pct"])
xs_r = np.linspace(0, 20, 50)
ax.plot(xs_r, slope_raw * xs_r + b_raw, color="#c0392b", lw=2, zorder=5,
        label=f"r={r_raw:.3f}\np={p_raw:.2g}")
ax.set_xlabel("Total repeat variants (n_snv_gene, clipped at 20)", fontsize=8)
ax.set_ylabel("Raw methylation (%)", fontsize=8)
ax.set_title("(C) Total variants vs raw methylation", fontsize=9)
ax.yaxis.set_major_formatter(pct_fmt)
ax.legend(fontsize=8); ax.grid(lw=0.3, alpha=0.35)

# D: total variants vs delta_meth
ax = axes[1, 0]
ax.scatter(total_snv + np.random.default_rng(1).uniform(-0.3, 0.3, len(df)),
           df["delta_meth"], s=3, alpha=0.12, color="#e07b3a", linewidths=0)
slope_d, b_d, r_d, p_d, _ = stats.linregress(df["n_snv_gene"], df["delta_meth"])
ax.plot(xs_r, slope_d * xs_r + b_d, color="#c0392b", lw=2, zorder=5,
        label=f"r={r_d:.3f}\np={p_d:.2g}")
ax.axhline(0, color="grey", lw=0.8, ls="--", alpha=0.5)
ax.set_xlabel("Total repeat variants (clipped at 20)", fontsize=8)
ax.set_ylabel("Δ methylation (copy − hap mean, %)", fontsize=8)
ax.set_title("(D) Total variants vs Δmeth", fontsize=9)
ax.yaxis.set_major_formatter(pct_fmt)
ax.legend(fontsize=8); ax.grid(lw=0.3, alpha=0.35)

# E: proportion methylated (>50%) per gene SNV group
ax = axes[1, 1]
bars_e = [df.loc[df["gene_snv_bin"] == lb, "is_methylated"].mean() * 100
          for lb in GENE_BIN_ORDER]
ns_e   = [df["gene_snv_bin"].eq(lb).sum() for lb in GENE_BIN_ORDER]
ax.bar(range(4), bars_e, color="#4d94c9", alpha=0.8, edgecolor="black", lw=0.5)
for i, (b, n) in enumerate(zip(bars_e, ns_e)):
    ax.text(i, b + 0.8, f"n={n:,}", ha="center", fontsize=7)
ax.set_xticks(range(4)); ax.set_xticklabels(GENE_BIN_ORDER, fontsize=9)
ax.set_xlabel("n_snv_5s_gene", fontsize=8)
ax.set_ylabel("% methylated copies (>50%)", fontsize=8)
ax.set_title("(E) Proportion methylated per gene-SNV group", fontsize=9)
ax.yaxis.set_major_formatter(pct_fmt)
ax.grid(axis="y", lw=0.3, alpha=0.4)

# F: proportion methylated per total-variant bin
ax = axes[1, 2]
bin_edges_f = [-1, 1, 4, 7, 200]
bin_labels_f = ["0–1", "2–4", "5–7", "≥8"]
df["total_snv_bin"] = pd.cut(df["n_snv_gene"], bins=bin_edges_f,
                              labels=bin_labels_f)
bars_f = [df.loc[df["total_snv_bin"] == lb, "is_methylated"].mean() * 100
          for lb in bin_labels_f]
ns_f   = [df["total_snv_bin"].eq(lb).sum() for lb in bin_labels_f]
ax.bar(range(4), bars_f, color="#e07b3a", alpha=0.8, edgecolor="black", lw=0.5)
for i, (b, n) in enumerate(zip(bars_f, ns_f)):
    ax.text(i, b + 0.8, f"n={n:,}", ha="center", fontsize=7)
ax.set_xticks(range(4)); ax.set_xticklabels(bin_labels_f, fontsize=9)
ax.set_xlabel("n_snv_gene (full repeat bins)", fontsize=8)
ax.set_ylabel("% methylated copies (>50%)", fontsize=8)
ax.set_title("(F) Proportion methylated per total-variant bin", fontsize=9)
ax.yaxis.set_major_formatter(pct_fmt)
ax.grid(axis="y", lw=0.3, alpha=0.4)

plt.tight_layout()
save_fig(fig, "21_variants_vs_methylation")

# data table 21
t21 = df[["copy_id", "sample_id", "hap_label", "superpopulation",
          "meth_pct", "delta_meth", "is_methylated",
          "n_snv_gene", "n_snv_5s_gene", "n_snv_nts_pre", "n_snv_nts_post",
          "gene_snv_bin"]].copy()
save_data(t21.round(3), "21_variants_vs_methylation")


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 22 — REGION VARIANTS VS METHYLATION
# ═══════════════════════════════════════════════════════════════════════════════

print("Fig 22 …", flush=True)

REGIONS22 = [
    ("5S rRNA gene (119 bp)",    "n_snv_5s_gene",   "meth_pct",
     [-1, 0, 1, 200], ["0", "1", "≥2"]),
    ("NTS-pre (629 bp)",         "n_snv_nts_pre",   "meth_pct",
     [-1, 0, 1, 3, 200], ["0", "1", "2–3", "≥4"]),
    ("NTS excl. ALU (1769 bp)",  "n_snv_other_nts", "meth_pct",
     [-1, 0, 2, 5, 200], ["0", "1–2", "3–5", "≥6"]),
    ("Full repeat (2168 bp)",    "n_snv_gene",      "meth_pct",
     [-1, 1, 4, 7, 200], ["0–1", "2–4", "5–7", "≥8"]),
]

fig, axes = plt.subplots(4, 2, figsize=(12, 16))
fig.suptitle(
    f"Methylation vs. regional variant count\n{cohort_title()}",
    fontsize=10, y=1.01
)

rows22 = []
for ri, (region_name, var_col, meth_col, bin_edges22, bin_labels22) in enumerate(REGIONS22):
    sub22 = df.dropna(subset=[meth_col]).copy()
    sub22["var_bin"] = pd.cut(sub22[var_col], bins=bin_edges22,
                               labels=bin_labels22)
    hap_mean22 = sub22.groupby("haplotype_id")[meth_col].transform("mean")
    sub22["delta22"] = sub22[meth_col] - hap_mean22

    for ci, (y_col, y_label, color) in enumerate(
            [(meth_col, "CpG methylation (%)", "#4d94c9"),
             ("delta22", "Δ (copy − hap mean, %)", "#e07b3a")]):
        ax = axes[ri, ci]
        grps22 = {lb: sub22.loc[sub22["var_bin"] == lb, y_col]
                  for lb in bin_labels22}
        decile_boxplot(ax, grps22, y_col, bin_labels22, color, y_label)
        if ci == 1:
            ax.axhline(0, color="grey", lw=0.8, ls="--", alpha=0.5)
        ax.set_title(f"{'ABCDEFGH'[ri*2+ci]}  {region_name}\n"
                     f"{'raw meth' if ci==0 else 'Δmeth (hap-centred)'}",
                     fontsize=8)
        ax.tick_params(labelsize=8)

        for lb in bin_labels22:
            vals = grps22[lb].dropna()
            rows22.append({"region": region_name, "var_col": var_col,
                            "bin": lb, "metric": "raw" if ci == 0 else "delta",
                            "n": len(vals),
                            "median": np.median(vals) if len(vals) else np.nan,
                            "q25": np.percentile(vals, 25) if len(vals) else np.nan,
                            "q75": np.percentile(vals, 75) if len(vals) else np.nan})

plt.tight_layout()
save_fig(fig, "22_region_variants_vs_methylation")
save_data(pd.DataFrame(rows22).round(3), "22_region_variants_vs_methylation")


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 23 — VARIANTS HYPOMETHYLATION STATS
# ═══════════════════════════════════════════════════════════════════════════════

print("Fig 23 …", flush=True)

REGIONS23 = [
    ("5S gene\n(119 bp)",       "n_snv_5s_gene",  [-1,0,1,200],   ["0","1","≥2"]),
    ("NTS-pre\n(629 bp)",       "n_snv_nts_pre",  [-1,0,1,3,200], ["0","1","2–3","≥4"]),
    ("NTS excl. ALU\n(1769 bp)", "n_snv_other_nts", [-1,0,2,5,200], ["0","1–2","3–5","≥6"]),
    ("Full repeat\n(2168 bp)",  "n_snv_gene",     [-1,1,4,7,200], ["0–1","2–4","5–7","≥8"]),
]

def wilson_ci(k, n, z=1.96):
    if n == 0: return np.nan, np.nan, np.nan
    p = k / n
    denom = 1 + z**2/n
    centre = (p + z**2/(2*n)) / denom
    half   = z * np.sqrt(p*(1-p)/n + z**2/(4*n**2)) / denom
    return p * 100, (centre - half) * 100, (centre + half) * 100

fig23, axes23 = plt.subplots(2, 4, figsize=(18, 9))
fig23.suptitle(
    f"Hypomethylation (<{HYPO_THR:.0f}%) vs. variant count, position confound check\n"
    f"{cohort_title()}",
    fontsize=10, y=1.01
)

rows23 = []
for ri, (region_name, var_col, bin_edges23, bin_labels23) in enumerate(REGIONS23):
    ax = axes23[0, ri]
    df["var_bin23"] = pd.cut(df[var_col], bins=bin_edges23, labels=bin_labels23)
    proportions, lo_ci, hi_ci = [], [], []
    ns23 = []
    for lb in bin_labels23:
        grp23 = df[df["var_bin23"] == lb]
        n = len(grp23)
        k = grp23["is_hypo"].sum()
        p, lo23, hi23 = wilson_ci(k, n)
        proportions.append(p)
        lo_ci.append(lo23)
        hi_ci.append(hi23)
        ns23.append(n)
        rows23.append({"region": region_name, "var_col": var_col, "bin": lb,
                       "n": n, "n_hypo": k, "pct_hypo": p,
                       "ci_lo": lo23, "ci_hi": hi23})
    xs23 = range(len(bin_labels23))
    ax.bar(xs23, proportions, color="#4d94c9", alpha=0.75, edgecolor="black", lw=0.5)
    for i, (p, lo23, hi23) in enumerate(zip(proportions, lo_ci, hi_ci)):
        if not (np.isnan(lo23) or np.isnan(hi23)):
            ax.plot([i, i], [lo23, hi23], color="#333333", lw=1.5, zorder=5)
        ax.text(i, (hi23 or p or 0) + 0.8, f"{ns23[i]:,}", ha="center",
                fontsize=6.5, rotation=30)

    # chi-squared
    conts = [[df.loc[df["var_bin23"] == lb, "is_hypo"].sum(),
              (df["var_bin23"] == lb).sum() - df.loc[df["var_bin23"] == lb, "is_hypo"].sum()]
             for lb in bin_labels23]
    try:
        chi2, pchi2, *_ = stats.chi2_contingency(np.array(conts))
        ax.set_title(f"{'ABCD'[ri]}1  {region_name}\n"
                     f"χ² p={pchi2:.2g}", fontsize=8.5)
    except Exception:
        ax.set_title(f"{'ABCD'[ri]}1  {region_name}", fontsize=8.5)
    ax.set_xticks(list(xs23))
    ax.set_xticklabels(bin_labels23, fontsize=9)
    ax.yaxis.set_major_formatter(pct_fmt)
    ax.set_ylabel(f"% hypomethylated (<{HYPO_THR:.0f}%)", fontsize=8)
    ax.set_xlabel(var_col, fontsize=7)
    ax.grid(axis="y", lw=0.3, alpha=0.4)

# Row 1: B = n_snv_5s_gene vs array position quintile
ax_b = axes23[1, 0]
df["pos_quintile"] = pd.qcut(df["pct_pos"], q=5,
                              labels=["Q1\n(0–20%)", "Q2", "Q3", "Q4", "Q5\n(80–100%)"])
for qi, q_label in enumerate(["Q1\n(0–20%)", "Q2", "Q3", "Q4", "Q5\n(80–100%)"]):
    grp_q = df.loc[df["pos_quintile"] == q_label, "n_snv_5s_gene"]
    if len(grp_q) == 0:
        continue
    q25q, q50q, q75q = np.percentile(grp_q, [25, 50, 75])
    ax_b.bar(qi, q75q - q25q, bottom=q25q, width=0.5, color="#7d5ba6", alpha=0.65)
    ax_b.plot([qi - 0.25, qi + 0.25], [q50q, q50q], color="black", lw=2)
ax_b.set_xticks(range(5))
ax_b.set_xticklabels(["Q1\n(0–20%)", "Q2", "Q3", "Q4", "Q5\n(80–100%)"], fontsize=8)
ax_b.set_ylabel("n_snv_5s_gene", fontsize=8)
ax_b.set_title("B  Gene-variant count vs.\narray-position quintile", fontsize=9)
ax_b.grid(axis="y", lw=0.3, alpha=0.4)

# C: logistic OR — n_snv_5s_gene before / after adding position covariate
ax_c = axes23[1, 1]
df["dist_from_centre"] = (df["pct_pos"] - 50).abs()
or_vals, ci_lo_or, ci_hi_or, labels_or = [], [], [], []
for var_col_c in ["n_snv_5s_gene", "n_snv_nts_pre", "n_snv_nts_post", "n_snv_gene"]:
    for label_c, formula_c in [
        (f"{var_col_c}\n(no pos)", f"is_hypo ~ {var_col_c}"),
        (f"{var_col_c}\n+position", f"is_hypo ~ {var_col_c} + dist_from_centre"),
    ]:
        try:
            mod = smf.logit(formula_c, data=df).fit(disp=False)
            coef = mod.params[var_col_c]
            lo_c = mod.conf_int().loc[var_col_c, 0]
            hi_c = mod.conf_int().loc[var_col_c, 1]
            or_vals.append(np.exp(coef))
            ci_lo_or.append(np.exp(lo_c))
            ci_hi_or.append(np.exp(hi_c))
        except Exception:
            or_vals.append(np.nan); ci_lo_or.append(np.nan); ci_hi_or.append(np.nan)
        labels_or.append(label_c)

xs_c = range(len(or_vals))
colors_c = ["#4d94c9", "#e07b3a"] * (len(or_vals) // 2)
ax_c.bar(xs_c, or_vals, color=colors_c, alpha=0.75, edgecolor="black", lw=0.5)
for i, (lo_c, hi_c) in enumerate(zip(ci_lo_or, ci_hi_or)):
    if not np.isnan(lo_c):
        ax_c.plot([i, i], [lo_c, hi_c], color="#333", lw=1.5, zorder=5)
ax_c.axhline(1.0, color="grey", lw=1, ls="--", alpha=0.6)
ax_c.set_xticks(list(xs_c))
ax_c.set_xticklabels(labels_or, fontsize=5.5, rotation=30, ha="right")
ax_c.set_ylabel("OR for hypomethylation", fontsize=8)
ax_c.set_title("C  Logistic OR — before/after\nposition adjustment", fontsize=9)
ax_c.grid(axis="y", lw=0.3, alpha=0.4)

for col_c in [2, 3]:
    axes23[1, col_c].axis("off")

plt.tight_layout()
save_fig(fig23, "23_variants_hypo_stats")
save_data(pd.DataFrame(rows23).round(3), "23_variants_hypo_stats")


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 29 — METHYLATION PROFILE BY COPY CLASS (3-region)
# ═══════════════════════════════════════════════════════════════════════════════

print("Fig 29 …", flush=True)

CLASS_ORDER  = ["High", "Intermediate", "Low"]
CLASS_COLORS = {"High": "#2c7fb8", "Intermediate": "#f0a202", "Low": "#d7191c"}
CLASS_LABELS = {"High": f"High ≥{HI_CUT*100:.0f}%",
                "Intermediate": f"Intermediate {LO_CUT*100:.0f}–{HI_CUT*100:.0f}%",
                "Low": f"Low <{LO_CUT*100:.0f}%"}

# 4 non-overlapping regions: NTS-pre, gene, ALU SINE, NTS excl. ALU (post-ALU spacer)
REGION_PAIRS = [
    ("NTS-pre (630 bp)",          "nts_pre_pct",   0),
    ("5S gene (119 bp)",          "gene_pct",      1),
    ("ALU SINE (279 bp)",         "alu_pct",       2),
    ("NTS excl. ALU (1769 bp)",   "other_nts_pct", 3),
]

fig, (ax_line, ax_bar) = plt.subplots(1, 2, figsize=(12, 5))
fig.suptitle(
    f"Within-unit methylation profile by copy class\n{cohort_title()}",
    fontsize=9, y=1.02
)

rows29 = []
for cls in CLASS_ORDER:
    sub29 = df[df["copy_class"] == cls]
    xs29 = [x for _, _, x in REGION_PAIRS]
    ys29 = [sub29[col].mean() for _, col, _ in REGION_PAIRS]
    col29 = CLASS_COLORS[cls]
    ax_line.plot(xs29, ys29, color=col29, lw=2.5, marker="o", ms=7,
                 label=f"{CLASS_LABELS[cls]}  (n={len(sub29):,})", zorder=4)
    for rname, col_name, _ in REGION_PAIRS:
        rows29.append({"copy_class": cls, "region": rname,
                       "n_copies": sub29[col_name].notna().sum(),
                       "mean_meth_pct": sub29[col_name].mean()})

ax_line.axvspan(1.5, 2.5, color="#d6604d", alpha=0.10, lw=0)  # ALU highlight
ax_line.set_xlim(-0.5, 3.5)
ax_line.set_xticks([0, 1, 2, 3])
ax_line.set_xticklabels([r for r, _, _ in REGION_PAIRS], fontsize=8)
ax_line.set_ylabel("Mean CpG methylation (%)", fontsize=9)
ax_line.set_title("(A) Regional methylation profile by copy class", fontsize=9)
ax_line.yaxis.set_major_formatter(pct_fmt)
ax_line.legend(fontsize=8, loc="best", framealpha=0.85)
ax_line.grid(axis="y", lw=0.3, alpha=0.4)

# bar version for clarity
x_pos29 = np.arange(4)
width = 0.25
for ci, cls in enumerate(CLASS_ORDER):
    sub29 = df[df["copy_class"] == cls]
    ys29  = [sub29[col].mean() for _, col, _ in REGION_PAIRS]
    ax_bar.bar(x_pos29 + (ci-1) * width, ys29, width * 0.9,
               color=CLASS_COLORS[cls], alpha=0.8,
               label=CLASS_LABELS[cls])
ax_bar.set_xticks(x_pos29)
ax_bar.set_xticklabels([r for r, _, _ in REGION_PAIRS], fontsize=8)
ax_bar.set_ylabel("Mean CpG methylation (%)", fontsize=9)
ax_bar.set_title("(B) Regional methylation by copy class (grouped bar)", fontsize=9)
ax_bar.yaxis.set_major_formatter(pct_fmt)
ax_bar.legend(fontsize=8, framealpha=0.85)
ax_bar.grid(axis="y", lw=0.3, alpha=0.4)

plt.tight_layout()
save_fig(fig, "29_methylation_profile_by_copy_class")
save_data(pd.DataFrame(rows29).round(2), "29_methylation_profile_by_copy_class")


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 30 — METHYLATION HEATMAP (3-region)
# ═══════════════════════════════════════════════════════════════════════════════

print("Fig 30 …", flush=True)

# ── try to build a truly positional heatmap from copy_meth_pos ────────────────
con = sqlite3.connect(DB_PATH)
_has_pos_table = con.execute(
    "SELECT COUNT(*) FROM sqlite_master WHERE type='table' AND name='copy_meth_pos'"
).fetchone()[0] > 0

GENE_S, GENE_E = 630, 748
ALU_S,  ALU_E  = 787, 1066
BIN_SIZE = 15
REPEAT_LEN_HM = 2168
N_BINS_HM = REPEAT_LEN_HM // BIN_SIZE   # 144 bins

if _has_pos_table:
    print("  copy_meth_pos table found — building positional heatmap …")
    from scipy.ndimage import convolve as _convolve
    from matplotlib.gridspec import GridSpec as _GS

    pos_data = pd.read_sql_query("""
        SELECT p.copy_id, p.wpos_bin,
               CAST(p.n_meth AS REAL) / p.n_conf AS bin_meth
        FROM copy_meth_pos p
        JOIN copy_methylation cm USING(copy_id)
        JOIN copy c ON p.copy_id=c.copy_id
        WHERE c.border_note='interior' AND cm.n_conf_calls>=10
          AND p.n_conf >= 2
    """, con)
    valid_ids = df["copy_id"].values
    pos_data  = pos_data[pos_data["copy_id"].isin(valid_ids)]

    bin_cols = np.arange(0, REPEAT_LEN_HM, BIN_SIZE)
    pivot = (pos_data.pivot_table(index="copy_id", columns="wpos_bin",
                                  values="bin_meth", aggfunc="mean")
                     .reindex(columns=bin_cols))
    sort_order = (df.set_index("copy_id")["meth_pct"]
                    .reindex(pivot.index).sort_values().index)
    pivot      = pivot.loc[sort_order]
    mat_pos    = pivot.values * 100        # (n_copies, 144), NaN where no CpG data

    n_cells   = mat_pos.size
    n_missing = int(np.isnan(mat_pos).sum())
    col_x     = bin_cols + BIN_SIZE / 2   # bin centres for x-axis
    col_mean_obs = np.nanmean(mat_pos, axis=0)

    # ── neighbour imputation ──
    def _impute(mat):
        valid   = (~np.isnan(mat)).astype(float)
        filled0 = np.where(np.isnan(mat), 0.0, mat)
        kernel  = np.ones((5, 7))           # 5 copy-rows × 7 position-cols
        ssum = _convolve(filled0, kernel, mode="nearest")
        cnt  = _convolve(valid,   kernel, mode="nearest")
        neigh = np.divide(ssum, cnt, out=np.full_like(ssum, np.nan), where=cnt > 0)
        out = np.where(np.isnan(mat), neigh, mat)
        if np.isnan(out).any():
            tmp = pd.DataFrame(out).interpolate(axis=1, limit_direction="both")
            col_fill = np.nanmean(out, axis=0)
            tmp = tmp.apply(lambda r: r.fillna(
                pd.Series(col_fill, index=tmp.columns)), axis=1)
            out = tmp.values
        return out

    mat_imp = _impute(mat_pos)
    n_residual = int(np.isnan(mat_imp).sum())

    # ── shared draw function ──────────────────────────────────────────────────
    cmap_hm = plt.get_cmap("RdYlGn").copy()
    cmap_hm.set_bad("#cccccc", alpha=0.6)   # NaN → light grey

    def _draw(mat, interp, imputed):
        n_copies = mat.shape[0]
        fig2 = plt.figure(figsize=(12, 13))
        gs2  = _GS(2, 1, figure=fig2, height_ratios=[1, 9], hspace=0.04)
        ax_top = fig2.add_subplot(gs2[0])
        ax_hm  = fig2.add_subplot(gs2[1])

        # top profile (always shows observed means, not imputed)
        ax_top.axvspan(GENE_S, GENE_E, color="#aec6cf", alpha=0.35, zorder=0)
        ax_top.axvspan(ALU_S,  ALU_E,  color="#c8a0e8", alpha=0.25, zorder=0)
        ax_top.plot(col_x, col_mean_obs, color="#222", lw=1.2)
        ax_top.set_xlim(0, REPEAT_LEN_HM); ax_top.set_ylim(0, 100)
        ax_top.yaxis.set_major_formatter(pct_fmt)
        ax_top.set_ylabel("Column\nmean", fontsize=8)
        ax_top.set_xticklabels([])
        ax_top.grid(axis="y", lw=0.3, alpha=0.4)
        imp_note = (f"neighbour-imputed  ({100*n_missing/n_cells:.0f}% cells filled"
                    + (f", {n_residual} residual" if n_residual else "")
                    + ")" if imputed else f"raw  ({100*n_missing/n_cells:.0f}% cells NaN → grey)")
        smooth_note = "bilinear display" if interp == "bilinear" else "no display smoothing"
        ax_top.set_title(
            f"Per-copy positional CpG methylation  ·  {n_copies:,} copies  ·  {BIN_SIZE} bp bins\n"
            f"{cohort_title()}  ·  {imp_note}  ·  {smooth_note}",
            fontsize=9)

        # heatmap
        im = ax_hm.imshow(mat, aspect="auto", cmap=cmap_hm, vmin=0, vmax=100,
                          extent=[0, REPEAT_LEN_HM, 0, n_copies],
                          origin="lower", interpolation=interp)
        for start, end, label in [(GENE_S, GENE_E, "5S gene"),
                                   (ALU_S,  ALU_E,  "ALU SINE")]:
            ax_hm.axvline(start, color="white", lw=0.8, ls="--", alpha=0.55)
            ax_hm.axvline(end,   color="white", lw=0.8, ls="--", alpha=0.55)
            ax_hm.text((start+end)/2, n_copies*0.995, label,
                       ha="center", va="top", fontsize=7,
                       color="white", fontweight="bold")
        plt.colorbar(im, ax=ax_hm, label="CpG methylation (%)",
                     fraction=0.02, pad=0.01)
        ax_hm.set_xlabel("Position within repeat unit (bp)", fontsize=9)
        ax_hm.set_ylabel("Individual copies (sorted by mean methylation ↑)", fontsize=9)
        ax_hm.set_xlim(0, REPEAT_LEN_HM)
        plt.tight_layout()
        return fig2

    # ── produce four outputs ──────────────────────────────────────────────────
    # 30: raw, nearest (no smoothing)
    save_fig(_draw(mat_pos, "nearest",  imputed=False), "30_methylation_heatmap")
    # 30 PNG: raw, bilinear (display smoothing)
    fig_png = _draw(mat_pos, "bilinear", imputed=False)
    fig_png.savefig(OUTDIR / "30_methylation_heatmap.png", dpi=200, bbox_inches="tight")
    plt.close(fig_png)
    # 30b: imputed, nearest
    save_fig(_draw(mat_imp, "nearest",  imputed=True), "30b_methylation_heatmap_imputed")
    # 30b PNG: imputed, bilinear
    fig_imp_png = _draw(mat_imp, "bilinear", imputed=True)
    fig_imp_png.savefig(OUTDIR / "30b_methylation_heatmap_imputed.png",
                        dpi=200, bbox_inches="tight")
    plt.close(fig_imp_png)

    # data tables
    hm_wide = pivot.copy() * 100
    hm_wide.columns = [f"bp{int(c)}" for c in hm_wide.columns]
    hm_wide.insert(0, "mean_meth_pct",
                   df.set_index("copy_id")["meth_pct"].reindex(hm_wide.index))
    hm_wide = hm_wide.reset_index()
    save_data(hm_wide.round(2), "30_methylation_heatmap")

    hm_wide_imp = pd.DataFrame(mat_imp, index=pivot.index,
                                columns=[f"bp{int(c)}" for c in bin_cols])
    hm_wide_imp.insert(0, "mean_meth_pct",
                       df.set_index("copy_id")["meth_pct"].reindex(hm_wide_imp.index))
    hm_wide_imp = hm_wide_imp.reset_index()
    save_data(hm_wide_imp.round(2), "30b_methylation_heatmap_imputed")

    print(f"  → raw: {len(mat_pos):,} copies × {N_BINS_HM} bins  "
          f"({100*n_missing/n_cells:.0f}% NaN)")
    print(f"  → imputed: {n_residual} residual NaN after neighbourhood fill")

    # ── Fig 30c: within-copy positional profiles by methylation class ─────────
    print("Fig 30c …", flush=True)

    METH_CUT = 65.0   # >65% = methylated, <65% = hypomethylated
    PCTILES_30C = [10, 25, 50, 75, 90]

    # annotate each copy_id with its class
    copy_cls = (df.set_index("copy_id")["meth_pct"]
                  .reindex(pos_data["copy_id"])
                  .values)
    pos_data_c = pos_data.copy()
    pos_data_c["cls"] = np.where(copy_cls >= METH_CUT, "high", "low")
    pos_data_c["meth_pct_bin"] = pos_data_c["bin_meth"] * 100

    groups = {
        f"Methylated (≥{METH_CUT:.0f}%)": ("high", "#2166AC"),
        f"Hypomethylated (<{METH_CUT:.0f}%)": ("low",  "#D6604D"),
    }
    n_hi = (pos_data_c["cls"] == "high").sum()  # copy×bin rows, not copies
    n_lo = (pos_data_c["cls"] == "low").sum()
    n_copies_hi = pos_data_c[pos_data_c["cls"]=="high"]["copy_id"].nunique()
    n_copies_lo = pos_data_c[pos_data_c["cls"]=="low"]["copy_id"].nunique()

    rows30c = []
    _CM = 1 / 2.54   # inches per cm
    _PX = 5 * _CM    # 5 cm panel edge
    fig30c, axes30c = plt.subplots(1, 2, figsize=(_PX*2 + 1.4, _PX + 1.0),
                                   sharey=True, constrained_layout=True)
    fig30c.suptitle(
        f"Within-copy positional CpG methylation profile by copy class\n"
        f"{cohort_title()}  ·  {BIN_SIZE} bp bins  ·  median + IQR + 10–90th pctile",
        fontsize=8
    )

    for ax, (label, (cls_key, color)) in zip(axes30c, groups.items()):
        sub = pos_data_c[pos_data_c["cls"] == cls_key]
        n_copies_cls = sub["copy_id"].nunique()

        stats_rows = []
        for bstart in bin_cols:
            bdata = sub.loc[sub["wpos_bin"] == bstart, "meth_pct_bin"].dropna()
            if len(bdata) < 3:
                continue
            q10, q25, q50, q75, q90 = np.percentile(bdata, PCTILES_30C)
            stats_rows.append({
                "wpos_bin": int(bstart),
                "wpos_centre": int(bstart) + BIN_SIZE / 2,
                "n_copies": len(bdata),
                "mean": float(bdata.mean()),
                "q10": q10, "q25": q25, "q50_median": q50,
                "q75": q75, "q90": q90,
                "class": label,
            })
            rows30c.append(stats_rows[-1])

        if not stats_rows:
            continue
        st = pd.DataFrame(stats_rows)
        xs = st["wpos_centre"].values
        ax.fill_between(xs, st["q10"], st["q90"],
                        color=color, alpha=0.15, label="10th–90th pctile")
        ax.fill_between(xs, st["q25"], st["q75"],
                        color=color, alpha=0.35, label="IQR (25th–75th)")
        ax.plot(xs, st["q50_median"], color=color, lw=1.6, label="Median")

        for fstart, fend, fname, fcolor in [
                (GENE_S, GENE_E, "5S gene", "#aec6cf"),
                (ALU_S,  ALU_E,  "ALU SINE", "#c8a0e8")]:
            ax.axvspan(fstart, fend, color=fcolor, alpha=0.25, zorder=0)
            ax.text((fstart + fend) / 2, 104, fname,
                    ha="center", va="bottom", fontsize=6.5, color="#444")

        ax.set_xlim(0, REPEAT_LEN_HM)
        ax.set_ylim(-2, 108)
        ax.yaxis.set_major_formatter(pct_fmt)
        ax.set_xlabel("Position within repeat unit (bp)", fontsize=8)
        if ax is axes30c[0]:
            ax.set_ylabel("CpG methylation (%)", fontsize=8)
        ax.set_title(f"{label}  ({n_copies_cls:,} copies)", fontsize=8)
        ax.legend(fontsize=7, framealpha=0.85, loc="lower right")
        ax.tick_params(labelsize=8)
        ax.grid(lw=0.3, alpha=0.4, axis="y")
        ax.set_box_aspect(1)

    save_fig(fig30c, "30c_methylation_profile_by_class_positional")
    tbl30c = pd.DataFrame(rows30c)
    # pivot to wide: one row per wpos_bin, columns for high/low quantiles
    wide30c = tbl30c.pivot_table(
        index="wpos_bin", columns="class",
        values=["n_copies", "mean", "q10", "q25", "q50_median", "q75", "q90"])
    wide30c.columns = ["_".join(str(c) for c in col).strip()
                       for col in wide30c.columns]
    wide30c = wide30c.reset_index()
    save_data(wide30c.round(2), "30c_methylation_profile_by_class_positional")
    print(f"  → 30c: {n_copies_hi:,} methylated + {n_copies_lo:,} hypo copies")

else:
    # ── fallback: 3-region heatmap (no copy_meth_pos table yet) ───────────────
    print("  copy_meth_pos not found — falling back to 3-region heatmap.")
    print("  Run scripts 37b + 38 (server then local) to get positional version.")
    hm = (df[["copy_id", "meth_pct", "nts_pre_pct", "gene_pct", "other_nts_pct"]]
          .dropna(subset=["nts_pre_pct", "gene_pct", "other_nts_pct"])
          .sort_values("meth_pct").copy())
    mat = hm[["nts_pre_pct", "gene_pct", "other_nts_pct"]].values

    def _make_heatmap_3reg(mat, interp):
        fig_h, axes_hm = plt.subplots(1, 1, figsize=(6, 10))
        im = axes_hm.imshow(mat, aspect="auto", cmap="RdYlGn", vmin=0, vmax=100,
                             extent=[-0.5, 2.5, 0, len(mat)],
                             origin="lower", interpolation=interp)
        plt.colorbar(im, ax=axes_hm, label="CpG methylation (%)", fraction=0.04, pad=0.02)
        axes_hm.set_xticks([0, 1, 2])
        axes_hm.set_xticklabels(["NTS-pre", "gene", "NTS (excl. ALU)"], fontsize=10)
        axes_hm.set_ylabel("Individual copies (sorted ↑)", fontsize=9)
        axes_hm.set_title(
            f"Per-copy methylation heatmap (3-region fallback)\n"
            f"{cohort_title()}  ·  {len(mat):,} copies", fontsize=9)
        plt.tight_layout()
        return fig_h

    save_fig(_make_heatmap_3reg(mat, "nearest"), "30_methylation_heatmap")
    fig_png = _make_heatmap_3reg(mat, "bilinear")
    fig_png.savefig(OUTDIR / "30_methylation_heatmap.png", dpi=200, bbox_inches="tight")
    plt.close(fig_png)
    hm_wide = hm[["copy_id", "meth_pct", "nts_pre_pct", "gene_pct", "other_nts_pct"]].copy()
    hm_wide.columns = ["copy_id", "mean_meth_pct", "NTS_pre_pct", "gene_pct", "NTS_exclALU_pct"]
    save_data(hm_wide.round(2), "30_methylation_heatmap")

con.close()


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 36 — LOW-METH COPIES DONOR STABILITY
# ═══════════════════════════════════════════════════════════════════════════════

print("Fig 36 …", flush=True)

# per-donor: all donors, all hap labels (hap1/hap2 and mat/pat), no minimum filter
# group by sample_id only so all 215 probands with any methylation data appear
donor36 = (df.groupby("sample_id")
             .agg(n_low=("is_hypo", "sum"),
                  n_total=("copy_id", "count"))
             .reset_index())
donor36["n_high"]   = donor36["n_total"] - donor36["n_low"]
donor36["frac_low"] = donor36["n_low"] / donor36["n_total"]
N36 = len(donor36)
print(f"  {N36} donors  (min {donor36['n_total'].min()} copies, "
      f"max {donor36['n_total'].max()} copies)", flush=True)

nlow, ntot, nhigh, frac = (donor36["n_low"].values, donor36["n_total"].values,
                            donor36["n_high"].values, donor36["frac_low"].values)

def cv36(x): return np.std(x, ddof=1) / np.mean(x)

r36, p36 = stats.pearsonr(ntot, nlow)
slope36, b36, *_ = stats.linregress(ntot, nlow)
pred_num  = np.full(N36, nlow.mean())
pred_frac = frac.mean() * ntot
rss_num   = float(np.sum((nlow - pred_num) ** 2))
rss_frac  = float(np.sum((nlow - pred_frac) ** 2))

print(f"  {N36} donors  | n_low mean={nlow.mean():.1f} CV={cv36(nlow):.2f}  "
      f"| set-point RSS={rss_num:.0f} vs prop RSS={rss_frac:.0f}",
      flush=True)


def _draw_fig36(donor_df, title_suffix, figname):
    nlo = donor_df["n_low"].values
    nto = donor_df["n_total"].values
    nhi = donor_df["n_high"].values
    fr  = donor_df["frac_low"].values
    N   = len(donor_df)
    r, _ = stats.pearsonr(nto, nlo)
    pred_n = np.full(N, nlo.mean())
    pred_f = fr.mean() * nto
    rn = float(np.sum((nlo - pred_n)**2))
    rf = float(np.sum((nlo - pred_f)**2))

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(
        f"Lowly-methylated (<{HYPO_THR:.0f}%) copies per donor\n"
        f"{cohort_title()}  ·  {title_suffix}",
        fontsize=10, y=1.01
    )

    ax = axes[0, 0]
    ax.hist(nlo, bins=np.arange(max(0, nlo.min()-0.5), nlo.max()+1.5),
            color="#2c7fb8", edgecolor="white", lw=0.4)
    ax.axvline(nlo.mean(), color="#c0392b", lw=2, label=f"mean={nlo.mean():.1f}")
    ax.axvspan(nlo.mean()-nlo.std(ddof=1), nlo.mean()+nlo.std(ddof=1),
               color="#c0392b", alpha=0.12, label=f"±1 SD ({nlo.std(ddof=1):.1f})")
    ax.set_xlabel("Low-meth copies per donor (<65%)", fontsize=9)
    ax.set_ylabel("Donors", fontsize=9)
    ax.set_title(f"(A) Per-donor count of lowly-meth copies  (CV={cv36(nlo):.2f})",
                 fontsize=9, loc="left")
    ax.legend(fontsize=8); ax.grid(axis="y", lw=0.3, alpha=0.4)

    ax = axes[0, 1]
    sc = ax.scatter(nto, nlo, s=35, c=nto, cmap="viridis_r", alpha=0.8,
                    edgecolors="white", lw=0.3, zorder=3)
    plt.colorbar(sc, ax=ax, label="Total copies", fraction=0.04, pad=0.02)
    xs = np.linspace(nto.min(), nto.max(), 50)
    ax.axhline(nlo.mean(), color="#c0392b", lw=2,
               label=f"fixed number, RSS={rn:.0f}")
    ax.plot(xs, fr.mean()*xs, color="#777", lw=2, ls="--",
            label=f"fixed fraction ({fr.mean()*100:.0f}%), RSS={rf:.0f}")
    ax.set_xlabel("Total classified copies per donor", fontsize=9)
    ax.set_ylabel("Lowly-methylated copies per donor", fontsize=9)
    ax.set_title(f"(B) Low-meth copies vs total copies  r={r:.2f}", fontsize=9, loc="left")
    ax.legend(fontsize=8, loc="upper left"); ax.grid(lw=0.3, alpha=0.4)

    ax = axes[1, 0]
    cv_vals = {"n_low\n(<65%)": cv36(nlo), "n_high\n(≥65%)": cv36(nhi),
               "n_total\ncopies": cv36(nto), "fraction\nlow": cv36(fr)}
    colors36 = ["#2c7fb8", "#41ab5d", "#888888", "#f0a202"]
    ax.bar(range(len(cv_vals)), list(cv_vals.values()), color=colors36,
           edgecolor="black", lw=0.5)
    for i, v in enumerate(cv_vals.values()):
        ax.text(i, v+0.004, f"{v:.2f}", ha="center", va="bottom", fontsize=9)
    ax.set_xticks(range(len(cv_vals)))
    ax.set_xticklabels(list(cv_vals.keys()), fontsize=9)
    ax.set_ylabel("Coefficient of variation (SD/mean)", fontsize=9)
    ax.set_title("(C) Coefficient of variation by quantity",
                 fontsize=9, loc="left")
    ax.grid(axis="y", lw=0.3, alpha=0.4)

    ax = axes[1, 1]
    ord_ = np.argsort(nto)
    xx   = np.arange(N)
    ax.bar(xx, nlo[ord_],  color="#2c7fb8", label=f"lowly-meth (<{HYPO_THR:.0f}%)")
    ax.bar(xx, nhi[ord_],  bottom=nlo[ord_], color="#dddddd",
           label=f"highly-meth (≥{HYPO_THR:.0f}%)")
    ax.plot(xx, nlo[ord_], color="#c0392b", lw=1.2, marker="o", ms=2.5,
            label="low count")
    ax.set_xlabel("Donor (sorted by total copy number →)", fontsize=9)
    ax.set_ylabel("Copies per donor", fontsize=9)
    ax.set_title("(D) Per-donor low- and high-meth copies, sorted by total",
                 fontsize=9, loc="left")
    ax.legend(fontsize=8, loc="upper left"); ax.set_xticks([])
    ax.grid(axis="y", lw=0.3, alpha=0.35)

    plt.tight_layout()
    save_fig(fig, figname)


# Fig 36: all donors (no filter)
_draw_fig36(donor36,
            f"{N36} probands — all donors, all haplotype labels",
            "36_lowmeth_copies_donor_stability")
save_data(donor36[["sample_id", "n_low", "n_high", "n_total", "frac_low"]].round(4),
          "36_lowmeth_copies_donor_stability")

# Fig 36b: quality-filtered (both hap1+hap2 with ≥20 covered copies each)
MIN_HAP_METH = 20
hap36b = (df.groupby(["sample_id", "hap_label"])
           .agg(n_low=("is_hypo", "sum"), n_total=("copy_id", "count"))
           .reset_index())
hap36b["n_high"] = hap36b["n_total"] - hap36b["n_low"]
hap36b = hap36b[hap36b["n_total"] >= MIN_HAP_METH]
g36b   = hap36b.groupby("sample_id").filter(
    lambda d: set(d["hap_label"]) >= {"hap1", "hap2"})
donor36b = g36b.groupby("sample_id").agg(
    n_low=("n_low", "sum"), n_high=("n_high", "sum"),
    n_total=("n_total", "sum")).reset_index()
donor36b["frac_low"] = donor36b["n_low"] / donor36b["n_total"]
print(f"  36b: {len(donor36b)} donors (hap1+hap2, ≥{MIN_HAP_METH} covered copies per hap)",
      flush=True)
_draw_fig36(donor36b,
            f"{len(donor36b)} probands — both hap1+hap2, ≥{MIN_HAP_METH} copies/hap",
            "36b_lowmeth_copies_donor_stability_filtered")
save_data(donor36b[["sample_id", "n_low", "n_high", "n_total", "frac_low"]].round(4),
          "36b_lowmeth_copies_donor_stability_filtered")


print("\nAll done.", flush=True)
print(f"Figures: {OUTDIR}")
print(f"Data:    {DATADIR}")
