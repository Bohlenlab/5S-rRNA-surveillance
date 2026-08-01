#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 94_heritability_vaf.py — twin-sharing and trio-inheritance heritability of 5S
# variant counts across a VAF-threshold sweep (Falconer midparent regression).
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
94_heritability_vaf.py

Computes twin-sharing and trio-inheritance analyses using a VAF-based
filter (VAF ≥ threshold) instead of AD read-count thresholds.

Uses the optimized 0.003 (0.3%) VAF cutoff as the primary threshold, and
sweeps across [0.001, 0.002, 0.003, 0.005, 0.01, 0.02, 0.05] to show the
full curve — directly analogous to the AD>0/1/2/5/10/20/50 sweep.

Inputs:
  carriers : results_500k/carriers_AD1_with_dp.tsv  (POS, REF, ALT, SAMPLE_ID, AD, DP, VAF)
  twins    : pairs.tsv                              (ID1, ID2)
  trios    : input/trios_annotate.tsv               (FID, IID, parent1_ID, parent2_ID, ...)
  scram    : input/trios_annotate_scram.tsv          (same format, scrambled IID)

Outputs:
  results_500k/plots_vaf/twin_sharing_vaf.pdf/.png
  results_500k/plots_vaf/trio_inheritance_vaf.pdf/.png
  results_500k/plots_vaf/twin_sharing_vaf_data.tsv
  results_500k/plots_vaf/trio_inheritance_vaf_data.tsv
"""

import os, random
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator
import numpy as np
import pandas as pd

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE      = Path(os.environ.get("FIVES_DATA", "data"))
CARRIERS  = BASE / "results_500k/carriers_AD1_with_dp.tsv"
PAIRS     = BASE / "pairs.tsv"
TRIOS     = BASE / "input/trios_annotate.tsv"
SCRAM     = BASE / "input/trios_annotate_scram.tsv"
OUT_DIR   = Path(os.environ.get("FIVES_OUT", "output")) / "02_variant_calling_qc/94_heritability"
OUT_DIR.mkdir(parents=True, exist_ok=True)

VAF_THRESHOLDS    = [0.001, 0.002, 0.003, 0.005, 0.01, 0.02, 0.05]
VAF_LABELS        = ["0.1%","0.2%","0.3%","0.5%","1%","2%","5%"]
OPTIMIZED_VAF     = 0.003   # highlighted cutoff
COHORT_FREQ_UPPER = 0.30    # exclude variants present in >30% of all cohort samples
MIN_VAF_FOR_FREQ  = 0.001   # VAF floor when computing cohort frequencies
SEED              = 42

# ── Load twin pairs and trio annotations ─────────────────────────────────────
print("Loading twin pairs and trio annotations …")
twins_df = pd.read_csv(PAIRS, sep="\t", dtype=str)
twins_pairs = [(r.ID1.strip(), r.ID2.strip()) for _, r in twins_df.iterrows()]

trios_df = pd.read_csv(TRIOS, sep="\t", dtype=str)[["IID","parent1_ID","parent2_ID"]].dropna()
scram_df = pd.read_csv(SCRAM, sep="\t", dtype=str)[["IID","parent1_ID","parent2_ID"]].dropna()

twin_ids = {id_ for pair in twins_pairs for id_ in pair}
trio_ids = (set(trios_df["IID"]) | set(trios_df["parent1_ID"]) | set(trios_df["parent2_ID"]) |
            set(scram_df["IID"]) | set(scram_df["parent1_ID"]) | set(scram_df["parent2_ID"]))
all_ids  = twin_ids | trio_ids

print(f"  {len(twins_pairs)} twin pairs ({len(twin_ids)} unique IDs)")
print(f"  {len(trios_df)} trios ({len(trio_ids)} unique IDs)")

# ── Load carriers — single pass for both cohort-frequency table and twin/trio data ─
print(f"Streaming {CARRIERS.name} for cohort freq + {len(all_ids)} twin/trio IDs …")
variant_sample_sets = defaultdict(set)   # (pos, ref, alt) → set of sample_ids
all_sample_ids_seen = set()
rel_chunks = []

for chunk in pd.read_csv(CARRIERS, sep="\t", chunksize=1_000_000,
                          dtype={"SAMPLE_ID": str, "POS": int,
                                 "AD": int, "DP": int, "VAF": float}):
    all_sample_ids_seen.update(chunk["SAMPLE_ID"].unique())

    # Cohort-frequency table: all rows at VAF >= MIN_VAF_FOR_FREQ
    sub_freq = chunk[chunk["VAF"] >= MIN_VAF_FOR_FREQ][["POS","REF","ALT","SAMPLE_ID"]]
    for key, grp in sub_freq.groupby(["POS","REF","ALT"]):
        variant_sample_sets[key].update(grp["SAMPLE_ID"].tolist())

    # Twin/trio carrier rows
    sub = chunk[chunk["SAMPLE_ID"].isin(all_ids)]
    if len(sub):
        rel_chunks.append(sub[["SAMPLE_ID","POS","REF","ALT","VAF"]])

total_samples = len(all_sample_ids_seen)
print(f"  Total unique samples in cohort: {total_samples:,}")

# Build exclusion set: variants present in >10% of cohort
common_variants = {var for var, sids in variant_sample_sets.items()
                   if len(sids) / total_samples > COHORT_FREQ_UPPER}
print(f"  {len(common_variants)} variants present in >{COHORT_FREQ_UPPER*100:.0f}% of samples — excluded")
del variant_sample_sets, all_sample_ids_seen  # free memory

carriers = pd.concat(rel_chunks, ignore_index=True)
print(f"  {len(carriers):,} twin/trio rows, {carriers.SAMPLE_ID.nunique()} unique IDs")

# Build per-sample dict: sample_id → list of (pos, ref, alt, vaf)
sample_variants = defaultdict(list)
for row in carriers.itertuples(index=False):
    sample_variants[row.SAMPLE_ID].append((row.POS, row.REF, row.ALT, row.VAF))
print(f"  Variant dict built for {len(sample_variants)} samples")

# ── Helper functions ──────────────────────────────────────────────────────────
def vset(sample_id, min_vaf):
    """Variant set with cohort-common variants excluded — used for twin sharing."""
    return {(p, r, a) for p, r, a, v in sample_variants.get(sample_id, [])
            if v >= min_vaf and (p, r, a) not in common_variants}

def vset_all(sample_id, min_vaf):
    """Variant set without cohort filter — used for trio h² regression."""
    return {(p, r, a) for p, r, a, v in sample_variants.get(sample_id, [])
            if v >= min_vaf}

def mean_parent_sharing_all(child_id, p1_id, p2_id, min_vaf):
    """Mean of pct_in_parent1 and pct_in_parent2 (30% cohort filter) — for beeswarm."""
    sc = vset(child_id, min_vaf)
    if not sc:
        return np.nan
    sp1 = vset(p1_id, min_vaf)
    sp2 = vset(p2_id, min_vaf)
    return (len(sc & sp1) / len(sc) + len(sc & sp2) / len(sc)) / 2

def symmetric_sharing(id1, id2, min_vaf):
    """Symmetric sharing fraction: mean(|s1∩s2|/|s1|, |s1∩s2|/|s2|), NaN if either empty."""
    s1 = vset(id1, min_vaf)
    s2 = vset(id2, min_vaf)
    fracs = []
    if s1: fracs.append(len(s1 & s2) / len(s1))
    if s2: fracs.append(len(s1 & s2) / len(s2))
    return np.mean(fracs) if fracs else np.nan

def estimate_h2(df, min_vaf, n_boot=1000, rng=None):
    """Narrow-sense h² via offspring-on-midparent regression (Falconer & Mackay).
    Phenotype = per-individual variant count at VAF >= min_vaf (all variants, no cohort filter).
    Returns (h2, ci_lo, ci_hi, n_trios, y_arr, x_arr).
    """
    if rng is None:
        rng = np.random.default_rng(42)
    y, x = [], []
    for row in df.itertuples(index=False):
        child, p1, p2 = row.IID, row.parent1_ID, row.parent2_ID
        y.append(len(vset_all(child, min_vaf)))
        x.append((len(vset_all(p1, min_vaf)) + len(vset_all(p2, min_vaf))) / 2)
    y = np.array(y, dtype=float)
    x = np.array(x, dtype=float)
    n = len(y)
    if np.var(x, ddof=1) == 0:
        return np.nan, np.nan, np.nan, n, y, x
    h2 = np.cov(y, x, ddof=1)[0, 1] / np.var(x, ddof=1)
    boot = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        xb, yb = x[idx], y[idx]
        vb = np.var(xb, ddof=1)
        if vb > 0:
            boot.append(np.cov(yb, xb, ddof=1)[0, 1] / vb)
    if len(boot) < 10:
        return h2, np.nan, np.nan, n, y, x
    ci_lo, ci_hi = np.percentile(boot, [2.5, 97.5])
    return h2, ci_lo, ci_hi, n, y, x

# ── Twin sharing analysis ─────────────────────────────────────────────────────
print("\nComputing twin sharing (real and scrambled) …")
random.seed(SEED)
id2s = [id2 for _, id2 in twins_pairs]
random.shuffle(id2s)
scram_pairs = [(id1, sid2) for (id1, _), sid2 in zip(twins_pairs, id2s)]

twin_real  = {v: [] for v in VAF_THRESHOLDS}   # list of (id1, id2, sharing)
twin_scram = {v: [] for v in VAF_THRESHOLDS}
missing = 0
for (id1, id2), (sid1, sid2) in zip(twins_pairs, scram_pairs):
    if id1 not in sample_variants or id2 not in sample_variants:
        missing += 1
        continue
    for vaf in VAF_THRESHOLDS:
        f = symmetric_sharing(id1, id2, vaf)
        if not np.isnan(f):
            twin_real[vaf].append((id1, id2, f))
        fs = symmetric_sharing(sid1, sid2, vaf)
        if not np.isnan(fs):
            twin_scram[vaf].append((sid1, sid2, fs))

if missing:
    print(f"  Warning: {missing} twin pairs missing carrier data")

print("  VAF threshold  |  Real (n, median)  |  Scrambled (n, median)")
print("  " + "-"*55)
for vaf, lbl in zip(VAF_THRESHOLDS, VAF_LABELS):
    r, s = twin_real[vaf], twin_scram[vaf]
    print(f"  ≥{lbl:<6}  |  n={len(r):<4} med={np.median([x[2] for x in r]):.3f}  "
          f"|  n={len(s):<4} med={np.median([x[2] for x in s]):.3f}")

# ── Trio heritability (offspring-on-midparent regression) ─────────────────────
print("\nComputing trio h² (offspring-on-midparent regression, 1000-rep bootstrap) …")
print("  (no cohort-frequency filter for h² — common variants cancel in regression)")
_rng = np.random.default_rng(SEED)
trio_h2_real  = []   # list of (h2, ci_lo, ci_hi, n, y, x) per threshold
trio_h2_scram = []

for min_vaf, lbl in zip(VAF_THRESHOLDS, VAF_LABELS):
    h2r, lo_r, hi_r, nr, yr, xr = estimate_h2(trios_df, min_vaf, rng=_rng)
    h2s, lo_s, hi_s, ns, ys, xs = estimate_h2(scram_df, min_vaf, rng=_rng)
    trio_h2_real.append((h2r, lo_r, hi_r, nr, yr, xr))
    trio_h2_scram.append((h2s, lo_s, hi_s, ns, ys, xs))
    print(f"  ≥{lbl:<6}  real h²={h2r:.3f} [{lo_r:.3f}, {hi_r:.3f}] n={nr}  |  "
          f"scram h²={h2s:.3f} [{lo_s:.3f}, {hi_s:.3f}]")

# ── Trio per-child sharing (beeswarm — no cohort filter, ~25 variants/child) ──
print("\nComputing trio per-child sharing (mean pct in parent1 + parent2) …")
trio_share_real  = {v: [] for v in VAF_THRESHOLDS}   # list of (child, p1, p2, sharing)
trio_share_scram = {v: [] for v in VAF_THRESHOLDS}

for df_rows, results in [(trios_df, trio_share_real), (scram_df, trio_share_scram)]:
    for row in df_rows.itertuples(index=False):
        child, p1, p2 = row.IID, row.parent1_ID, row.parent2_ID
        for vaf in VAF_THRESHOLDS:
            f = mean_parent_sharing_all(child, p1, p2, vaf)
            if not np.isnan(f):
                results[vaf].append((child, p1, p2, f))

print("  VAF threshold  |  Real (n, median)  |  Scrambled (n, median)")
print("  " + "-"*55)
for vaf, lbl in zip(VAF_THRESHOLDS, VAF_LABELS):
    r, s = trio_share_real[vaf], trio_share_scram[vaf]
    print(f"  ≥{lbl:<6}  |  n={len(r):<4} med={np.median([x[3] for x in r]):.3f}  "
          f"|  n={len(s):<4} med={np.median([x[3] for x in s]):.3f}")

# ── Save data tables ──────────────────────────────────────────────────────────
rows_twin = []
for vaf, lbl in zip(VAF_THRESHOLDS, VAF_LABELS):
    for id1, id2, v in twin_real[vaf]:
        rows_twin.append({"type":"real",  "vaf_threshold":vaf, "label":lbl, "id1":id1, "id2":id2, "sharing":v})
    for id1, id2, v in twin_scram[vaf]:
        rows_twin.append({"type":"scram", "vaf_threshold":vaf, "label":lbl, "id1":id1, "id2":id2, "sharing":v})
pd.DataFrame(rows_twin).to_csv(OUT_DIR / "twin_sharing_vaf_data.tsv", sep="\t", index=False)

rows_trio = []
for (h2r, lo_r, hi_r, nr, *_), (h2s, lo_s, hi_s, ns, *_), vaf, lbl in zip(
        trio_h2_real, trio_h2_scram, VAF_THRESHOLDS, VAF_LABELS):
    rows_trio.append({"type":"real",  "vaf_threshold":vaf, "label":lbl,
                      "h2":h2r, "ci_lo":lo_r, "ci_hi":hi_r, "n":nr})
    rows_trio.append({"type":"scram", "vaf_threshold":vaf, "label":lbl,
                      "h2":h2s, "ci_lo":lo_s, "ci_hi":hi_s, "n":ns})
pd.DataFrame(rows_trio).to_csv(OUT_DIR / "trio_heritability_vaf_data.tsv", sep="\t", index=False)

rows_trio_share = []
for vaf, lbl in zip(VAF_THRESHOLDS, VAF_LABELS):
    for child, p1, p2, v in trio_share_real[vaf]:
        rows_trio_share.append({"type":"real",  "vaf_threshold":vaf, "label":lbl,
                                "child_id":child, "parent1_id":p1, "parent2_id":p2, "mean_parent_sharing":v})
    for child, p1, p2, v in trio_share_scram[vaf]:
        rows_trio_share.append({"type":"scram", "vaf_threshold":vaf, "label":lbl,
                                "child_id":child, "parent1_id":p1, "parent2_id":p2, "mean_parent_sharing":v})
pd.DataFrame(rows_trio_share).to_csv(OUT_DIR / "trio_inheritance_vaf_data.tsv", sep="\t", index=False)
print("\nData tables saved.")

# ── Beeswarm helper ───────────────────────────────────────────────────────────
def beeswarm_x(values, center, width=0.32, n_bins=40):
    values = np.array(values)
    if len(values) == 0:
        return np.array([])
    y_min, y_max = values.min(), values.max()
    if y_max == y_min:
        return np.full(len(values), center)
    bins = np.linspace(y_min - 1e-9, y_max + 1e-9, n_bins + 1)
    xs = np.empty(len(values))
    for i in range(n_bins):
        mask = (values >= bins[i]) & (values < bins[i+1])
        n = mask.sum()
        if n == 0:
            continue
        offsets = np.linspace(-width/2, width/2, n) if n > 1 else np.array([0.0])
        rng = np.random.default_rng(i)
        rng.shuffle(offsets)
        xs[mask] = center + offsets
    return xs

def draw_panel(ax, real_dict, scram_dict, vaf_thresholds, vaf_labels,
               ylabel, title, optimized_vaf, n_pairs):
    gap, sep = 1.2, 0.55
    xtick_pos, xtick_labels = [], []

    for i, (vaf, lbl) in enumerate(zip(vaf_thresholds, vaf_labels)):
        center  = i * gap
        x_scram = center - sep / 2
        x_real  = center + sep / 2

        rv = np.array([x[-1] for x in real_dict[vaf]])
        sv = np.array([x[-1] for x in scram_dict[vaf]])

        s_real  = 18 if len(rv) < 200 else 4
        s_scram = 18 if len(sv) < 200 else 4

        xs_r = beeswarm_x(rv, x_real)
        xs_s = beeswarm_x(sv, x_scram)

        ax.scatter(xs_s, sv, s=s_scram, color="black",   alpha=0.5, linewidths=0, zorder=2)
        ax.scatter(xs_r, rv, s=s_real,  color="#2166ac",  alpha=0.6, linewidths=0, zorder=2)

        bar_w = 0.22
        ax.plot([x_scram - bar_w, x_scram + bar_w], [np.median(sv)]*2,
                color="#e31a1c", linewidth=2.5, zorder=4, solid_capstyle="round")
        ax.plot([x_real  - bar_w, x_real  + bar_w], [np.median(rv)]*2,
                color="#e31a1c", linewidth=2.5, zorder=4, solid_capstyle="round")

        # Highlight optimized VAF column
        if abs(vaf - optimized_vaf) < 1e-9:
            ax.axvspan(center - gap/2 + 0.05, center + gap/2 - 0.05,
                       color="#ffffcc", alpha=0.55, zorder=0)
            ax.text(center, 1.03, "★ 0.3%", ha="center", va="bottom",
                    fontsize=8.5, color="#996600", fontweight="bold",
                    transform=ax.get_xaxis_transform())

        xtick_pos   += [x_scram, x_real]
        xtick_labels += [lbl, lbl]

    ax.set_xticks(xtick_pos)
    ax.set_xticklabels(xtick_labels, fontsize=9)
    ax.set_xlim(-gap * 0.6, (len(vaf_thresholds) - 1) * gap + gap * 0.6)
    ax.set_ylim(-0.02, 1.08)
    ax.yaxis.set_major_locator(MultipleLocator(0.25))
    ax.set_xlabel("VAF threshold (left dot = scrambled, right = real)", fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(f"{title}  (n={n_pairs})", fontsize=11)
    ax.grid(axis="y", alpha=0.25, linestyle="--")
    ax.spines[["top","right"]].set_visible(False)

    handles = [
        ax.scatter([], [], s=25, color="#2166ac", label="Real"),
        ax.scatter([], [], s=25, color="black",   label="Scrambled"),
        Line2D([0],[0], color="#e31a1c", linewidth=2.5, label="Median"),
    ]
    ax.legend(handles=handles, fontsize=9, framealpha=0.9, loc="upper left")

# ── Figure 1: Twin sharing ────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(14, 7))
draw_panel(ax, twin_real, twin_scram, VAF_THRESHOLDS, VAF_LABELS,
           ylabel="Fraction of variants shared\nbetween twin pair",
           title=f"Identical twin variant sharing — real vs scrambled",
           optimized_vaf=OPTIMIZED_VAF,
           n_pairs=len(twins_pairs))
fig.suptitle("Twin variant sharing by VAF threshold\n"
             f"(UKBB 500k; cohort-common variants excluded [>{COHORT_FREQ_UPPER*100:.0f}% of samples]; optimized cutoff ★ = 0.3%)",
             fontsize=12, fontweight="bold", y=1.02)
fig.tight_layout()
fig.savefig(OUT_DIR / "twin_sharing_vaf.pdf", bbox_inches="tight")
fig.savefig(OUT_DIR / "twin_sharing_vaf.png", dpi=200, bbox_inches="tight", facecolor="white")
plt.close()
print(f"Saved: {OUT_DIR}/twin_sharing_vaf.pdf")

# ── Figure 2: Trio heritability line plot ─────────────────────────────────────
x_pos = np.arange(len(VAF_THRESHOLDS))

h2r_vals  = np.array([t[0] for t in trio_h2_real])
lo_r_vals = np.array([t[1] for t in trio_h2_real])
hi_r_vals = np.array([t[2] for t in trio_h2_real])
h2s_vals  = np.array([t[0] for t in trio_h2_scram])
lo_s_vals = np.array([t[1] for t in trio_h2_scram])
hi_s_vals = np.array([t[2] for t in trio_h2_scram])

fig, ax = plt.subplots(figsize=(10, 6))

ax.plot(x_pos, h2r_vals, "o-", color="#2166ac", lw=2, ms=7, zorder=3, label="Real trios")
mask_r = ~np.isnan(lo_r_vals)
ax.fill_between(x_pos[mask_r], lo_r_vals[mask_r], hi_r_vals[mask_r],
                color="#2166ac", alpha=0.2, zorder=2)

ax.plot(x_pos, h2s_vals, "s--", color="black", lw=1.5, ms=5, zorder=3, label="Scrambled trios")
mask_s = ~np.isnan(lo_s_vals)
ax.fill_between(x_pos[mask_s], lo_s_vals[mask_s], hi_s_vals[mask_s],
                color="gray", alpha=0.15, zorder=2)

ax.axhline(0, color="#888", lw=0.8, ls=":")
ax.axhline(1, color="#888", lw=0.8, ls=":", label="h²=1 (full heritability)")

opt_idx = VAF_THRESHOLDS.index(OPTIMIZED_VAF)
ax.axvspan(opt_idx - 0.4, opt_idx + 0.4, color="#ffffcc", alpha=0.55, zorder=0)
ymax_annot = np.nanmax(np.concatenate([h2r_vals, [1.0]])) * 1.05
ax.text(opt_idx, ymax_annot, "★ 0.3%", ha="center", va="bottom",
        fontsize=9, color="#996600", fontweight="bold")

ax.set_xticks(x_pos)
ax.set_xticklabels(VAF_LABELS, fontsize=10)
ax.set_xlabel("VAF threshold", fontsize=11)
ax.set_ylabel("Narrow-sense heritability (h²)\noffspring ~ midparent regression", fontsize=11)
ax.set_title(f"5S rDNA variant count heritability — trios  (n={len(trios_df)})", fontsize=11)
ax.legend(fontsize=10, framealpha=0.9)
ax.spines[["top", "right"]].set_visible(False)
ax.grid(axis="y", alpha=0.25, ls="--")

fig.suptitle("Trio heritability by VAF threshold (Falconer offspring-on-midparent regression, 95% CI)\n"
             f"(UKBB 500k; cohort-common variants excluded [>{COHORT_FREQ_UPPER*100:.0f}% of samples]; optimized cutoff ★ = 0.3%)",
             fontsize=11, fontweight="bold", y=1.02)
fig.tight_layout()
fig.savefig(OUT_DIR / "trio_heritability_vaf.pdf", bbox_inches="tight")
fig.savefig(OUT_DIR / "trio_heritability_vaf.png", dpi=200, bbox_inches="tight", facecolor="white")
plt.close()
print(f"Saved: {OUT_DIR}/trio_heritability_vaf.pdf")

# ── Figure 3: Per-child scatter at optimized VAF ──────────────────────────────
opt_idx = VAF_THRESHOLDS.index(OPTIMIZED_VAF)
_, _, _, _, yr_opt, xr_opt = trio_h2_real[opt_idx]
_, _, _, _, ys_opt, xs_opt = trio_h2_scram[opt_idx]
h2_opt = h2r_vals[opt_idx]

# Regression line
x_line = np.linspace(0, max(xr_opt.max(), xs_opt.max()), 100)
y_line = np.mean(yr_opt) + h2_opt * (x_line - np.mean(xr_opt))

fig, axes = plt.subplots(1, 2, figsize=(13, 6), sharey=True)
for ax, xv, yv, color, label, is_real in [
        (axes[0], xr_opt, yr_opt, "#2166ac", "Real trios", True),
        (axes[1], xs_opt, ys_opt, "black",   "Scrambled trios", False)]:
    ax.scatter(xv, yv, s=8, alpha=0.35, color=color, linewidths=0)
    if is_real:
        ax.plot(x_line, y_line, color="#e31a1c", lw=1.8, ls="-",
                label=f"h²={h2_opt:.3f} [{lo_r_vals[opt_idx]:.3f}, {hi_r_vals[opt_idx]:.3f}]")
    else:
        h2s_opt = h2s_vals[opt_idx]
        y_scram = np.mean(ys_opt) + h2s_opt * (x_line - np.mean(xs_opt))
        ax.plot(x_line, y_scram, color="#e31a1c", lw=1.8, ls="--",
                label=f"h²={h2s_opt:.3f} [{lo_s_vals[opt_idx]:.3f}, {hi_s_vals[opt_idx]:.3f}]")
    ax.set_xlabel("Mid-parent variant count", fontsize=11)
    ax.set_ylabel("Offspring variant count", fontsize=11)
    ax.set_title(f"{label}  (n={len(xv)})", fontsize=11)
    ax.legend(fontsize=9.5, framealpha=0.9)
    ax.spines[["top","right"]].set_visible(False)
    ax.grid(alpha=0.2, ls="--")

fig.suptitle(f"Trio heritability scatter — VAF ≥ {OPTIMIZED_VAF*100:.1f}% (Falconer offspring-on-midparent regression)\n"
             f"(UKBB 500k; all variants included; optimized cutoff ★ = 0.3%)",
             fontsize=11, fontweight="bold", y=1.02)
fig.tight_layout()
fig.savefig(OUT_DIR / "trio_heritability_scatter.pdf", bbox_inches="tight")
fig.savefig(OUT_DIR / "trio_heritability_scatter.png", dpi=200, bbox_inches="tight", facecolor="white")
plt.close()
print(f"Saved: {OUT_DIR}/trio_heritability_scatter.pdf")

# ── Figure 4: Trio beeswarm — per-child mean parent sharing ──────────────────
print("\nGenerating trio inheritance beeswarm …")
fig, ax = plt.subplots(figsize=(14, 7))
draw_panel(ax, trio_share_real, trio_share_scram, VAF_THRESHOLDS, VAF_LABELS,
           ylabel="Mean fraction of child variants\nfound in parent1 or parent2",
           title=f"Trio inheritance — real vs scrambled  (mean sharing per parent)",
           optimized_vaf=OPTIMIZED_VAF,
           n_pairs=len(trios_df))
fig.suptitle("Trio inheritance by VAF threshold\n"
             f"(UKBB 500k; all variants included, no cohort filter; optimized cutoff ★ = 0.3%)",
             fontsize=12, fontweight="bold", y=1.02)
fig.tight_layout()
fig.savefig(OUT_DIR / "trio_inheritance_vaf.pdf", bbox_inches="tight")
fig.savefig(OUT_DIR / "trio_inheritance_vaf.png", dpi=200, bbox_inches="tight", facecolor="white")
plt.close()
print(f"Saved: {OUT_DIR}/trio_inheritance_vaf.pdf")
print("\nDone.")
