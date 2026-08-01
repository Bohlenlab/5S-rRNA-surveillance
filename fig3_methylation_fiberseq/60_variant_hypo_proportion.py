#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 60_variant_hypo_proportion.py — Methylation distributions and proportion-low bars binned by variant count for the 5S gene, ALU, and NTS-post regions (ONT).
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
60_variant_hypo_proportion.py

Truncated violins and proportion-low bars for 5S gene, ALU SINE, and NTS
(excl. ALU) regions. Variant bins for gene and ALU: 0 / 1 / 2 / ≥3;
NTS bins: 0 / 1 / 2-3 / 4-5 / ≥6.

Each panel 5×5 cm · 300 DPI · font 8 pt · lines 1 pt

Output: <FIVES_OUT>/03_methylation_full215/24b_variant_hypo_proportion.pdf

Paths are read from environment variables (see repository README):
    FIVES_DB   path to 5S_rDNA.db
    FIVES_OUT  output directory
"""

import os
import sqlite3, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats as sc
from scipy.stats import kruskal, mannwhitneyu
from pathlib import Path

# ── global style ──────────────────────────────────────────────────────────────
FS   = 8      # font size pt
LW   = 1.0    # line width pt
_CM  = 1/2.54
_PX  = 5 * _CM   # 5 cm panel edge

matplotlib.rcParams.update({
    "font.size":        FS,
    "axes.titlesize":   FS,
    "axes.labelsize":   FS,
    "xtick.labelsize":  FS,
    "ytick.labelsize":  FS,
    "legend.fontsize":  FS,
    "lines.linewidth":  LW,
    "patch.linewidth":  LW,
    "axes.linewidth":   LW,
    "xtick.major.width": LW,
    "ytick.major.width": LW,
    "pdf.fonttype":     42,
    "font.family":      "sans-serif",
})

# ── paths ─────────────────────────────────────────────────────────────────────
DB     = Path(os.environ.get("FIVES_DB", "5S_rDNA.db"))
OUTDIR = Path(os.environ.get("FIVES_OUT", "output")) / "03_methylation_full215"
DDIR   = OUTDIR / "data"

HYPO_THR  = 65.0
MIN_CALLS = 10
COHORTS   = ("HPRC_Year1", "HPRC_Release2")

# ── load data ─────────────────────────────────────────────────────────────────
con = sqlite3.connect(DB)
df = pd.read_sql(f"""
    SELECT cm.copy_id, cm.mean_meth,
           COALESCE(cm.alu_n,0)    AS alu_n,
           COALESCE(cm.alu_meth,0) AS alu_meth,
           cm.gene_n, cm.gene_meth,
           cm.nts_pre_n, cm.nts_pre_meth,
           cm.nts_post_n, cm.nts_post_meth,
           c.n_snv_5s_gene, c.n_snv_nts_pre, c.n_snv_nts_post
    FROM copy_methylation cm
    JOIN copy c  ON cm.copy_id       = c.copy_id
    JOIN haplotype h ON c.haplotype_id = h.haplotype_id
    JOIN assembly  a ON h.assembly_id  = a.assembly_id
    WHERE c.border_note='interior' AND cm.n_conf_calls>={MIN_CALLS}
      AND a.cohort IN ({','.join(repr(c) for c in COHORTS)})
""", con)
alu_snv = pd.read_sql("""
    SELECT copy_id, COUNT(*) AS n_snv_alu
    FROM variant
    WHERE consensus_pos>=787 AND consensus_pos<1066
      AND alignment_source='gene_unit_t2t' AND var_type='snp'
    GROUP BY copy_id
""", con)
con.close()

df = df.merge(alu_snv, on="copy_id", how="left")
df["n_snv_alu"] = df["n_snv_alu"].fillna(0).astype(int)

def safe_pct(meth, n): return np.where(n > 0, meth / n * 100.0, np.nan)

df["meth_pct"]     = df["mean_meth"] * 100.0
df["alu_pct"]      = safe_pct(df["alu_meth"],      df["alu_n"])
df["nts_post_pct"] = safe_pct(df["nts_post_meth"], df["nts_post_n"])
# "other NTS" = NTS-pre + (NTS-post minus ALU) — all non-gene, non-ALU sequence
df["other_nts_n"]     = df["nts_pre_n"]    + (df["nts_post_n"]    - df["alu_n"]).clip(lower=0)
df["other_nts_meth"]  = df["nts_pre_meth"] + (df["nts_post_meth"] - df["alu_meth"]).clip(lower=0)
df["other_nts_pct"]   = safe_pct(df["other_nts_meth"], df["other_nts_n"])
df["n_snv_other_nts"] = (df["n_snv_nts_pre"] + (df["n_snv_nts_post"] - df["n_snv_alu"])
                         .clip(lower=0)).astype(int)

n_total = len(df)

# ── region definitions ────────────────────────────────────────────────────────
# Identical bins for gene and ALU: 0 / 1 / 2 / ≥3
# NTS-post: 0 / 1 / 2-3 / 4-5 / ≥6
SHARED_BINS   = [-1, 0, 1, 2, 200]
SHARED_LABELS = ["0", "1", "2", "≥3"]

REGIONS = [
    ("5S gene (119 bp)",   "n_snv_5s_gene",  "meth_pct",
     SHARED_BINS, SHARED_LABELS, "#2166AC"),
    ("ALU SINE (279 bp)",  "n_snv_alu",       "meth_pct",
     SHARED_BINS, SHARED_LABELS, "#D6604D"),
    ("NTS excl. ALU (~1769 bp)", "n_snv_other_nts", "meth_pct",
     [-1, 0, 1, 3, 5, 200], ["0", "1", "2–3", "4–5", "≥6"], "#4DAC26"),
]

# ── Wilson CI ────────────────────────────────────────────────────────────────
def wilson_ci(k, n, z=1.96):
    if n == 0: return np.nan, np.nan, np.nan
    p = k / n
    denom  = 1 + z**2 / n
    centre = (p + z**2 / (2*n)) / denom
    half   = z * np.sqrt(p*(1-p)/n + z**2/(4*n**2)) / denom
    return p*100, max(0,(centre-half)*100), min(100,(centre+half)*100)

# ── violin helper ─────────────────────────────────────────────────────────────
def draw_violin(ax, groups, order, color):
    data_per_bin = [groups.get(lb, pd.Series(dtype=float)).dropna().values
                    for lb in order]
    non_empty = [(i, d) for i, d in enumerate(data_per_bin) if len(d) > 1]
    positions = [i+1 for i, _ in non_empty]
    data_vals = [d for _, d in non_empty]

    if data_vals:
        parts = ax.violinplot(data_vals, positions=positions,
                              showmedians=False, showextrema=False,
                              points=300, widths=0.7)
        for pc in parts["bodies"]:
            pc.set_facecolor(color); pc.set_alpha(0.45)
            pc.set_edgecolor(color); pc.set_linewidth(LW)

    for pos, d in zip(positions, data_vals):
        q5, q25, q50, q75, q95 = np.percentile(d, [5, 25, 50, 75, 95])
        ax.plot([pos, pos], [q5, q95], color=color, lw=LW,
                solid_capstyle="round")
        ax.plot([pos-.15, pos+.15], [q25, q25], color=color, lw=LW)
        ax.plot([pos-.15, pos+.15], [q75, q75], color=color, lw=LW)
        ax.plot(pos, q50, "o", color="white", ms=3.5, zorder=5,
                markeredgecolor=color, markeredgewidth=LW)
        ax.text(pos, -8, f"n={len(d):,}", ha="center", va="top",
                fontsize=FS-1, color="#666")

    ax.axhline(HYPO_THR, color="#c0392b", lw=LW, ls="--", alpha=0.75)
    ax.set_xticks(range(1, len(order)+1))
    ax.set_xticklabels(order)
    ax.set_ylim(-13, 105)
    ax.yaxis.set_major_formatter(
        plt.FuncFormatter(lambda v, _: f"{v:.0f}%" if v >= 0 else ""))
    ax.grid(axis="y", lw=0.5, alpha=0.35)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_box_aspect(1)

    # Kruskal-Wallis across all bins + Mann-Whitney each bin vs bin-0
    all_groups = [groups.get(lb, pd.Series(dtype=float)).dropna().values for lb in order]
    non_empty  = [d for d in all_groups if len(d) > 1]
    ref_group  = all_groups[0]  # 0-variant bin
    if len(non_empty) >= 2:
        kw_stat, kw_p = kruskal(*non_empty)
        kw_str = f"KW p={kw_p:.2g}" if kw_p >= 1e-4 else "KW p<1e-4"
        ax.set_title(kw_str, loc="right", style="italic", color="#444", fontsize=FS-1)
    if len(ref_group) > 1:
        for pos, lb in enumerate(order[1:], start=2):
            grp = groups.get(lb, pd.Series(dtype=float)).dropna().values
            if len(grp) > 1:
                _, mw_p = mannwhitneyu(ref_group, grp, alternative="two-sided")
                stars = ("***" if mw_p < 0.001 else
                         "**"  if mw_p < 0.01  else
                         "*"   if mw_p < 0.05  else "ns")
                ax.text(pos, 97, stars, ha="center", va="top",
                        fontsize=FS-1, color="#333")

# ── proportion bar helper ─────────────────────────────────────────────────────
def draw_proportion(ax, groups, order, color):
    ps, lo_errs, hi_errs, ns = [], [], [], []
    for lb in order:
        sub = groups.get(lb, pd.Series(dtype=float)).dropna()
        k   = (sub < HYPO_THR).sum()
        p, lo, hi = wilson_ci(k, len(sub))
        ps.append(p); ns.append(len(sub))
        lo_errs.append(0 if np.isnan(p) else p - lo)
        hi_errs.append(0 if np.isnan(p) else hi - p)

    xs = np.arange(len(order))
    valid = [not np.isnan(p) for p in ps]
    ax.bar([x for x,v in zip(xs,valid) if v],
           [p for p,v in zip(ps,valid)  if v],
           color=color, alpha=0.75, edgecolor="white", lw=LW, zorder=3, width=0.6)
    ax.errorbar([x for x,v in zip(xs,valid) if v],
                [p for p,v in zip(ps,valid)  if v],
                yerr=[[e for e,v in zip(lo_errs,valid) if v],
                      [e for e,v in zip(hi_errs,valid) if v]],
                fmt="none", color="#333", lw=LW, capsize=2.5, zorder=4)
    hi_err_max = max((e for e,v in zip(hi_errs,valid) if v), default=3)
    for x, p, n, v in zip(xs, ps, ns, valid):
        if v:
            ax.text(x, p + hi_err_max*0.2 + 1.5, f"{p:.1f}%",
                    ha="center", va="bottom", fontsize=FS-1,
                    fontweight="bold", color="#333")
        else:
            ax.text(x, 2, "n=0", ha="center", va="bottom",
                    fontsize=FS-1, color="#aaa", style="italic")

    ax.set_xticks(xs); ax.set_xticklabels(order)
    p_vals = [p for p,v in zip(ps,valid) if v]
    ax.set_ylim(0, min(100, max(p_vals + [5]) * 1.6) if p_vals else 10)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.0f}%"))
    ax.grid(axis="y", lw=0.5, alpha=0.35)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_box_aspect(1)

    # chi-square (skip empty bins)
    obs = [(groups.get(lb, pd.Series(dtype=float)).dropna() < HYPO_THR).sum()
           for lb in order]
    tot = [len(groups.get(lb, pd.Series(dtype=float)).dropna()) for lb in order]
    obs_f = [o for o,t in zip(obs,tot) if t > 0]
    tot_f = [t for t in tot if t > 0]
    exp_frac = sum(obs_f) / max(sum(tot_f), 1)
    try:
        chi2, p_chi = sc.chisquare(obs_f, [exp_frac*t for t in tot_f])
        p_str = f"p={p_chi:.2g}" if p_chi >= 1e-4 else "p<1e-4"
        ax.set_title(f"χ²={chi2:.1f}  {p_str}", loc="right",
                     style="italic", color="#444")
    except Exception:
        pass

# ── build figure ──────────────────────────────────────────────────────────────
n_rows = len(REGIONS)
# each panel ~5×5 cm; add room for labels and inter-panel gaps
fig_w = 2 * _PX + 1.5 * _CM    # two panels + inter-column gap + margins
fig_h = n_rows * _PX + 2.5 * _CM  # rows + suptitle + bottom margin

fig, axes = plt.subplots(n_rows, 2,
                         figsize=(fig_w, fig_h),
                         constrained_layout=True)
fig.suptitle(
    f"Proportion <{HYPO_THR:.0f}% methylated by variant count\n"
    f"HPRC · {n_total:,} interior copies",
    fontsize=FS)

letters = "ABCDEFGHI"
for ri, (region_name, var_col, meth_col, bin_edges, bin_labels, color) in enumerate(REGIONS):
    sub = df.dropna(subset=[meth_col]).copy()
    sub["var_bin"] = pd.cut(sub[var_col], bins=bin_edges, labels=bin_labels)
    groups = {lb: sub.loc[sub["var_bin"] == lb, meth_col] for lb in bin_labels}

    ax_v, ax_p = axes[ri, 0], axes[ri, 1]

    draw_violin(ax_v, groups, bin_labels, color)
    ax_v.set_ylabel("CpG methylation (%)")
    ax_v.set_xlabel("Number of sequence variants")
    ax_v.set_title(f"({letters[ri*2]}) {region_name}", loc="left")
    if ri == 0:
        ax_v.text(0.98, HYPO_THR + 2, f"<{HYPO_THR:.0f}%",
                  transform=ax_v.get_yaxis_transform(), ha="right",
                  fontsize=FS-1, color="#c0392b")

    draw_proportion(ax_p, groups, bin_labels, color)
    ax_p.set_ylabel(f"% copies <{HYPO_THR:.0f}%")
    ax_p.set_xlabel("Number of sequence variants")
    ax_p.set_title(f"({letters[ri*2+1]}) {region_name}", loc="left")

out = OUTDIR / "24b_variant_hypo_proportion.pdf"
fig.savefig(out, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"Saved: {out}")

# ── update data tables ────────────────────────────────────────────────────────
def make_wide(series_dict, col_prefix):
    max_n = max((len(v.dropna()) for v in series_dict.values()), default=0)
    out   = {}
    for label, vals in series_dict.items():
        arr = vals.dropna().values
        out[f"{col_prefix}={label}"] = np.concatenate(
            [arr, np.full(max_n - len(arr), np.nan)])
    return pd.DataFrame(out)

# gene: 0/1/2/≥3 → meth_pct
gene_dict = {lb: df.loc[pd.cut(df["n_snv_5s_gene"],
             bins=SHARED_BINS, labels=SHARED_LABELS) == lb, "meth_pct"]
             for lb in SHARED_LABELS}
make_wide(gene_dict, "n_snv_5s_gene").to_csv(
    DDIR / "24b_gene_variants_meth_wide.tsv", sep="\t", index=False)

# ALU: 0/1/2/≥3 → meth_pct
alu_dict = {lb: df.loc[pd.cut(df["n_snv_alu"],
            bins=SHARED_BINS, labels=SHARED_LABELS) == lb, "meth_pct"]
            for lb in SHARED_LABELS}
make_wide(alu_dict, "n_snv_alu").to_csv(
    DDIR / "24b_alu_variants_meth_wide.tsv", sep="\t", index=False)

# NTS excl. ALU: 0/1/2-3/4-5/≥6 → meth_pct
nts_labels = ["0","1","2–3","4–5","≥6"]
nts_dict   = {lb: df.loc[pd.cut(df["n_snv_other_nts"],
              bins=[-1,0,1,3,5,200], labels=nts_labels) == lb, "meth_pct"]
              for lb in nts_labels}
make_wide(nts_dict, "n_snv_other_nts").to_csv(
    DDIR / "24b_nts_variants_meth_wide.tsv", sep="\t", index=False)

print("Data tables updated.")

# ── summary ───────────────────────────────────────────────────────────────────
print(f"\n{'Region':<22} {'Bin':<6} {'n':>7} {'%<65':>7}  95% CI")
print("-" * 55)
for region_name, var_col, meth_col, bin_edges, bin_labels, _ in REGIONS:
    sub = df.dropna(subset=[meth_col]).copy()
    sub["var_bin"] = pd.cut(sub[var_col], bins=bin_edges, labels=bin_labels)
    for lb in bin_labels:
        s = sub.loc[sub["var_bin"] == lb, meth_col].dropna()
        k = (s < HYPO_THR).sum()
        p, lo, hi = wilson_ci(k, len(s))
        ci = f"[{lo:.1f}–{hi:.1f}%]" if not np.isnan(p) else "—"
        pstr = f"{p:.1f}%" if not np.isnan(p) else "n/a"
        print(f"  {region_name[:20]:<20} {lb:<6} {len(s):>7,} {pstr:>7}  {ci}")
    print()
