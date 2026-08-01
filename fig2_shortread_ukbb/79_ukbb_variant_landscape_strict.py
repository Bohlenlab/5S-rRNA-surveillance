#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 79_ukbb_variant_landscape_strict.py — UK Biobank 5S rDNA variant-landscape
# figure regenerated at two stricter VAF threshold ranges.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
79_ukbb_variant_landscape_strict.py

Seven-panel variant-landscape layout regenerated for two stricter threshold ranges:
  78b  —  primary 0.5%,  AFS comparison 0.5% / 1% / 2%
  78c  —  primary 1.0%,  AFS comparison 1% / 3% / 5%

Loads UKBB blobs once, calls make_figure() twice.
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

GENE_COLOR = "#fff3b0"
REG_COLORS = {"gene": "#e67e22", "nts_pre": "#2980b9", "nts_post": "#27ae60"}
REG_ORDER  = ["gene", "nts_pre", "nts_post"]
REG_XLABS  = ["Gene\n(630–748)", "NTS pre\n(467–629)", "NTS post\n(749–967)"]

# All thresholds we need across both figures
ALL_THRESH = {
    "0.50%": 0.005,
    "1.00%": 0.010,
    "2.00%": 0.020,
    "3.00%": 0.030,
    "5.00%": 0.050,
}

# ── load data (once) ──────────────────────────────────────────────────────────

print("Loading UKBB data …", flush=True)
con = sqlite3.connect(DB)

raw = con.execute("""
    SELECT t2t_pos, t2t_ref, t2t_alt, region,
           n_carriers_ad1, n_carriers_ad5, mean_vaf, median_vaf, vaf_array
    FROM ukbb_population_variants ORDER BY t2t_pos
""").fetchall()

depth_df = pd.DataFrame(con.execute("""
    SELECT t2t_pos, median_dp, p5_dp, p25_dp, p75_dp, p95_dp
    FROM ukbb_depth_profile ORDER BY t2t_pos
""").fetchall(), columns=["t2t_pos","median_dp","p5_dp","p25_dp","p75_dp","p95_dp"])

tp_rows = con.execute("""
    SELECT DISTINCT consensus_pos, ref, alt FROM (
        SELECT v.consensus_pos, v.ref, v.alt
        FROM variant v
        JOIN copy c ON v.copy_id = c.copy_id
        JOIN haplotype h ON c.haplotype_id = h.haplotype_id
        JOIN assembly a ON h.assembly_id = a.assembly_id
        WHERE v.alignment_source = 'gene_unit_t2t' AND a.cohort = 'HPRC_Year1'
          AND v.consensus_pos BETWEEN ? AND ?
        UNION
        SELECT rv.consensus_pos, rv.ref, rv.alt
        FROM read_variant rv
        JOIN assembly a ON rv.assembly_id = a.assembly_id
        WHERE rv.modality = 'hifi' AND rv.vaf IS NOT NULL
          AND a.cohort = 'HPRC_Year1'
          AND rv.consensus_pos BETWEEN ? AND ?
    )
""", (WIN_LO, WIN_HI, WIN_LO, WIN_HI)).fetchall()
tp_set = {(int(r[0]), r[1], r[2]) for r in tp_rows}
con.close()

# Build dataframe + store raw vaf arrays for Panel F
records   = []
vaf_store = {}   # (pos, ref, alt) → sorted float32 array

for pos, ref, alt, reg, n1, n5, mv, mdv, blob in raw:
    vafs = np.frombuffer(blob, dtype=np.float32)
    key  = (int(pos), ref, alt)
    vaf_store[key] = vafs
    rec  = dict(t2t_pos=int(pos), ref=ref, alt=alt, region=reg or "unknown",
                n_ad1=int(n1), n_ad5=int(n5),
                mean_vaf=float(mv or 0), median_vaf=float(mdv or 0),
                is_tp=(key in tp_set))
    for lab, T in ALL_THRESH.items():
        rec[f"n_{lab}"] = int(len(vafs) - np.searchsorted(vafs, T))
    records.append(rec)

df = pd.DataFrame(records)
print(f"  {len(df):,} variants loaded;  {len(tp_set):,} TPs")

# ── helpers ───────────────────────────────────────────────────────────────────

_COMP = {"A":"T","T":"A","C":"G","G":"C"}

def canonical_snv(ref, alt):
    if len(ref) != 1 or len(alt) != 1 or ref == alt:
        return "other"
    if ref in ("C","T"):
        return f"{ref}>{alt}"
    return f"{_COMP.get(ref,'?')}>{_COMP.get(alt,'?')}"

SNV_ORDER  = ["C>A","C>G","C>T","T>A","T>C","T>G"]
SNV_COLORS = ["#1f77b4","#000000","#d62728","#cccccc","#2ca02c","#ff7f0e"]

BUCK_EDGES  = [1, 2, 11, 101, 1001, N_TOTAL + 1]
BUCK_LABELS = ["1", "2–10", "11–100", "101–1000", ">1000"]
BUCK_COLORS = ["#e41a1c","#ff7f00","#4daf4a","#377eb8","#984ea3"]

def bucket(n):
    for i in range(len(BUCK_EDGES) - 1):
        if BUCK_EDGES[i] <= n < BUCK_EDGES[i + 1]:
            return i
    return len(BUCK_EDGES) - 2

# ── figure factory ─────────────────────────────────────────────────────────────

def make_figure(thresh_main_lab, thresh_b_labs, outname):
    """
    thresh_main_lab  : e.g. "0.50%"
    thresh_b_labs    : list of 3 labels for Panel B comparison, e.g. ["0.50%","1.00%","2.00%"]
    outname          : output PDF stem
    """
    T_main = ALL_THRESH[thresh_main_lab]
    col_main = f"n_{thresh_main_lab}"

    df_called = df[df[col_main] >= 1].copy()
    df_called["bucket"] = df_called[col_main].apply(bucket)
    n_tp_called = int(df_called["is_tp"].sum())

    # pooled VAFs above main threshold by region (for Panel F)
    vafs_by_reg = {reg: [] for reg in REG_ORDER}
    for rec in records:
        if rec[col_main] < 1:
            continue
        key  = (rec["t2t_pos"], rec["ref"], rec["alt"])
        vafs = vaf_store[key]
        above = vafs[vafs >= T_main]
        if len(above):
            vafs_by_reg[rec["region"]].append(above)
    vafs_by_reg = {reg: np.concatenate(v) if v else np.array([], dtype=np.float32)
                   for reg, v in vafs_by_reg.items()}

    print(f"\n{outname}  —  primary {thresh_main_lab}:  "
          f"{len(df_called):,} variants called,  {n_tp_called:,} TPs", flush=True)

    # ── layout ────────────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(21, 17))
    gs  = gridspec.GridSpec(3, 3, figure=fig, hspace=0.52, wspace=0.38,
                            top=0.93, bottom=0.06, left=0.06, right=0.97)

    fig.suptitle(
        f"UKBB 5S rDNA variant landscape  ·  n=490,075 samples  ·  "
        f"Primary threshold: {thresh_main_lab} VAF  "
        f"(comparison: {' / '.join(thresh_b_labs)})",
        fontsize=12, fontweight="bold")

    # ── A: position scatter ───────────────────────────────────────────────────
    ax_A = fig.add_subplot(gs[0, :])

    for reg in REG_ORDER:
        sub    = df_called[(df_called["region"] == reg) & ~df_called["is_tp"]]
        sub_tp = df_called[(df_called["region"] == reg) &  df_called["is_tp"]]
        for subset, edge, lw, zo in [(sub, "none", 0, 3), (sub_tp, "red", 0.8, 4)]:
            sizes = np.clip(np.log10(subset[col_main].values + 1) * 12, 4, 80)
            ax_A.scatter(subset["t2t_pos"], subset[col_main],
                         s=sizes, c=REG_COLORS[reg], alpha=0.70, zorder=zo,
                         edgecolors=edge, linewidths=lw,
                         label=f"{reg} (n={len(sub)+len(sub_tp):,})"
                               if edge == "none" else None)

    ax_A.scatter([], [], s=20, c="grey", edgecolors="red", linewidths=0.8,
                 label=f"assembly ∪ HiFi confirmed ({n_tp_called:,})")
    ax_A.axvspan(GENE_LO, GENE_HI, color=GENE_COLOR, alpha=0.5, zorder=1)
    ax_A.set_yscale("log")
    ax_A.set_xlim(WIN_LO - 5, WIN_HI + 5)
    ax_A.set_xlabel("T2T consensus position")
    ax_A.set_ylabel(f"Carriers at VAF ≥ {thresh_main_lab}  (log)")
    ax_A.set_title(
        f"A  Variant map: {len(df_called):,} variants  "
        f"(marker size ∝ log₁₀ carrier count;  red outline = assembly ∪ HiFi confirmed)",
        fontsize=10, loc="left")
    ax_A.legend(fontsize=8, title="Region", title_fontsize=8, loc="upper right")
    ax_A.grid(True, lw=0.3, alpha=0.4)

    ax_A2 = ax_A.twinx()
    ax_A2.fill_between(depth_df["t2t_pos"], depth_df["p25_dp"], depth_df["p75_dp"],
                       alpha=0.12, color="grey")
    ax_A2.plot(depth_df["t2t_pos"], depth_df["median_dp"],
               color="grey", lw=1, alpha=0.5)
    ax_A2.set_ylabel("Read depth (right axis)", color="grey", fontsize=8)
    ax_A2.tick_params(axis="y", labelcolor="grey", labelsize=7)
    ax_A2.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v/1e3:.1f}k"))

    # ── B: allele frequency spectrum ──────────────────────────────────────────
    ax_B = fig.add_subplot(gs[1, 0])

    # Dynamically set upper bin bound from actual max carrier count
    max_n = max((df[f"n_{lab}"].max() for lab in thresh_b_labs), default=1)
    bins_afs = np.logspace(0, np.log10(max(max_n, 2)), 50)

    styles = [
        dict(histtype="step",       lw=1.5, color="#1f77b4"),
        dict(histtype="stepfilled", lw=2.0, color="#d62728", alpha=0.45),
        dict(histtype="step",       lw=1.5, color="#2ca02c"),
    ]
    for lab, kw in zip(thresh_b_labs, styles):
        vals = df[df[f"n_{lab}"] >= 1][f"n_{lab}"].values
        ax_B.hist(vals, bins=bins_afs, label=f"{lab} (n={len(vals):,})", **kw)

    ax_B.set_xscale("log")
    ax_B.set_xlabel("Carriers (log)")
    ax_B.set_ylabel("Variants")
    ax_B.set_title("B  Allele frequency spectrum", fontsize=10, loc="left")
    ax_B.legend(fontsize=8); ax_B.grid(True, lw=0.3, alpha=0.4)
    ax_B.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{int(v):,}"))

    # ── C: per-region breakdown ───────────────────────────────────────────────
    ax_C = fig.add_subplot(gs[1, 1])

    var_counts   = [len(df_called[df_called["region"] == r])               for r in REG_ORDER]
    med_carriers = [df_called[df_called["region"] == r][col_main].median() for r in REG_ORDER]
    win_lengths  = [(GENE_HI - GENE_LO + 1), (GENE_LO - WIN_LO), (WIN_HI - GENE_HI)]

    x = np.arange(3)
    bars_C = ax_C.bar(x, var_counts, 0.5,
                      color=[REG_COLORS[r] for r in REG_ORDER], alpha=0.85)
    for bar, cnt, wl in zip(bars_C, var_counts, win_lengths):
        ax_C.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                  f"{cnt}\n({cnt/wl:.1f}/nt)", ha="center", va="bottom", fontsize=7.5)

    ax_C.set_xticks(x); ax_C.set_xticklabels(REG_XLABS, fontsize=8)
    ax_C.set_ylabel(f"Variants called at {thresh_main_lab}")
    ax_C.set_title("C  Variant count per region  (rate per nt in labels)",
                   fontsize=10, loc="left")

    ax_C2 = ax_C.twinx()
    ax_C2.plot(x, med_carriers, "D--", color="black", lw=1.5, ms=7, zorder=5)
    ax_C2.set_ylabel("Median carrier count  (◆)", color="black", fontsize=8)
    ax_C2.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax_C.grid(True, lw=0.3, alpha=0.4, axis="y")

    # ── D: substitution spectrum ──────────────────────────────────────────────
    ax_D = fig.add_subplot(gs[1, 2])

    df_snv = df_called.copy()
    df_snv["snv_type"] = [canonical_snv(r, a)
                          for r, a in zip(df_snv["ref"], df_snv["alt"])]
    df_snv = df_snv[df_snv["snv_type"] != "other"]

    snv_by_reg = {}
    for reg in REG_ORDER:
        sub = df_snv[df_snv["region"] == reg]
        snv_by_reg[reg] = {t: int((sub["snv_type"] == t).sum()) for t in SNV_ORDER}

    x_D = np.arange(len(SNV_ORDER))
    bar_w = 0.25
    for i, reg in enumerate(REG_ORDER):
        counts = [snv_by_reg[reg].get(t, 0) for t in SNV_ORDER]
        ax_D.bar(x_D + (i - 1) * bar_w, counts, bar_w,
                 color=REG_COLORS[reg], alpha=0.85, label=reg)
    ax_D.set_xticks(x_D); ax_D.set_xticklabels(SNV_ORDER, fontsize=9)
    ax_D.set_xlabel("Canonical SNV type (pyrimidine reference)")
    ax_D.set_ylabel(f"Variants at {thresh_main_lab}")
    ax_D.set_title("D  Substitution spectrum by region", fontsize=10, loc="left")
    ax_D.legend(fontsize=7); ax_D.grid(True, lw=0.3, alpha=0.4, axis="y")

    # ── E: rarity distribution ────────────────────────────────────────────────
    ax_E = fig.add_subplot(gs[2, 0])

    bucket_counts = [int((df_called["bucket"] == i).sum()) for i in range(len(BUCK_LABELS))]
    total_called  = len(df_called)

    bars_E = ax_E.bar(range(len(BUCK_LABELS)), bucket_counts,
                      color=BUCK_COLORS, alpha=0.85)
    for bar, cnt in zip(bars_E, bucket_counts):
        if cnt == 0:
            continue
        pct = 100 * cnt / total_called
        ax_E.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                  f"{cnt}\n({pct:.1f}%)", ha="center", va="bottom", fontsize=7.5)
    ax_E.set_xticks(range(len(BUCK_LABELS)))
    ax_E.set_xticklabels(BUCK_LABELS, fontsize=8)
    ax_E.set_xlabel(f"Carrier count at {thresh_main_lab}")
    ax_E.set_ylabel("Variants")
    ax_E.set_title(f"E  Rarity distribution at {thresh_main_lab}", fontsize=10, loc="left")
    ax_E.grid(True, lw=0.3, alpha=0.4, axis="y")

    # ── F: per-carrier VAF distribution ──────────────────────────────────────
    ax_F = fig.add_subplot(gs[2, 1])

    all_above = np.concatenate(list(vafs_by_reg.values())) \
                if any(len(v) for v in vafs_by_reg.values()) else np.array([])
    vmax = float(all_above.max()) if len(all_above) else 1.0
    bins_vaf = np.linspace(T_main, min(1.0, vmax + 0.01), 70)

    for reg, color in REG_COLORS.items():
        v = vafs_by_reg.get(reg, np.array([]))
        if len(v):
            ax_F.hist(v, bins=bins_vaf, color=color, alpha=0.55,
                      label=f"{reg} (n={len(v):,})", density=True)

    ax_F.axvline(T_main, color="red", lw=1.5, ls="--",
                 label=f"{thresh_main_lab} threshold")
    ax_F.set_xlabel(f"VAF per carrier  (VAF ≥ {thresh_main_lab})")
    ax_F.set_ylabel("Density")
    ax_F.set_title("F  Per-carrier VAF distribution (pooled)", fontsize=10, loc="left")
    ax_F.legend(fontsize=7); ax_F.grid(True, lw=0.3, alpha=0.4)
    ax_F.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v*100:.1f}%"))

    # ── G: top 25 most frequent variants ─────────────────────────────────────
    ax_G = fig.add_subplot(gs[2, 2])

    top25    = df_called.nlargest(25, col_main).sort_values(col_main)
    n_vals   = top25[col_main].values
    labels_G = [f"pos{r.t2t_pos} {r.ref}>{r.alt}  [{r.region}]"
                for r in top25.itertuples()]
    colors_G = [REG_COLORS.get(r.region, "grey") for r in top25.itertuples()]
    edges_G  = ["red" if r.is_tp else "none" for r in top25.itertuples()]
    lws_G    = [0.8 if r.is_tp else 0.0 for r in top25.itertuples()]

    bars_G = ax_G.barh(range(len(top25)), n_vals, color=colors_G, alpha=0.85)
    for bar, ec, lw in zip(bars_G, edges_G, lws_G):
        bar.set_edgecolor(ec); bar.set_linewidth(lw)
    ax_G.set_yticks(range(len(top25)))
    ax_G.set_yticklabels(labels_G, fontsize=6.5)
    ax_G.set_xlabel(f"Carriers at VAF ≥ {thresh_main_lab}")
    ax_G.set_title(f"G  Top 25 most frequent variants at {thresh_main_lab}",
                   fontsize=10, loc="left")
    ax_G.grid(True, lw=0.3, alpha=0.4, axis="x")
    max_n = int(n_vals.max()) if len(n_vals) else 1
    for i, (n, pct) in enumerate(zip(n_vals, 100 * n_vals / N_TOTAL)):
        ax_G.text(n + max_n * 0.01, i, f" {pct:.2f}%", va="center", fontsize=6)

    out = OUTDIR / f"{outname}.pdf"
    plt.savefig(out, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"  → {out}")

# ── generate both figures ─────────────────────────────────────────────────────

make_figure("0.50%", ["0.50%", "1.00%", "2.00%"], "78b_ukbb_variant_landscape_mid")
make_figure("1.00%", ["1.00%", "3.00%", "5.00%"], "78c_ukbb_variant_landscape_strict")
