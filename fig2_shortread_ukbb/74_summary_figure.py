#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 74_summary_figure.py — six-panel summary of short-read 5S rDNA variant-calling
# performance (VAF distributions, cross-platform concordance, sensitivity, precision).
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
74_summary_figure.py

Six-panel summary figure for the HPRC 5S rDNA SR variant calling method.

A  SR VAF by copy count — why low thresholds are needed (violin)
B  Call landscape — stacked log histogram by group
C  Cross-platform validation — SR VAF vs HiFi VAF scatter
D  Sensitivity by copy count at VAF >= 0.30%
E  Precision by SR VAF bin (any-LR TP definition)
F  Precision-Recall curve at operating point
"""

import os
import sqlite3
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
plt.rcParams.update({
    "font.size": 8,
    "axes.titlesize": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "figure.dpi": 300,
})
from pathlib import Path

DB     = Path(os.environ.get("FIVES_DB", "5S_rDNA.db"))
BAM    = Path(os.environ.get("FIVES_DATA", "data")) / "bam"
OUTDIR = Path(os.environ.get("FIVES_OUT", "output")) / "02_variant_calling_qc"

COHORT          = "HPRC_Year1"
BOUNDARY_EXCL    = set(range(1, 90)) | set(range(2034, 2169))
EXCL_ALL  = {"HG02486"}   # no usable HiFi or SR data
EXCL_HIFI = {"HG02818"}   # HiFi alignment failed (48× SR/HiFi ratio); SR kept
SR_MIN_AD       = 3
LR_AD_MIN       = 5
LR_VAF_MIN_HIFI = 0.003   # 0.3%
VAF_CALL        = 0.003   # reporting threshold
FOLD_MAX        = 3       # SR/HiFi VAFs must be within this many-fold to count as confirmed

NC_STRATA = [
    ("nc=1",     1,  1),
    ("nc=2",     2,  2),
    ("nc=3",     3,  3),
    ("nc=4",     4,  4),
    ("nc=5",     5,  5),
    ("nc=6–10",  6, 10),
    ("nc=11–20", 11, 20),
    ("nc>20",    21, 9999),
]

VAF_THRESHOLDS = np.concatenate([
    np.linspace(0.0005, 0.003, 26),
    np.linspace(0.004,  0.010, 13),
    np.linspace(0.012,  0.030,  5),
])

# Coarse VAF bins for precision panel
PREC_EDGES  = [0, 0.001, 0.002, 0.003, 0.005, 0.01, 0.02, 0.05, 0.20, 1.01]
PREC_LABELS = ["<0.1%", "0.1–0.2%", "0.2–0.3%", "0.3–0.5%",
               "0.5–1%", "1–2%", "2–5%", "5–20%", ">20%"]

COLORS = {
    "asm_multi":  "#2ca25f",
    "asm_sing":   "#74c476",
    "lr_rescued": "#fd8d3c",
    "fp_noise":   "#bdbdbd",
}

# Panel B fine-grained nc breakdown (lightest → darkest = fewest → most copies)
COLORS_B = {
    "asm_nc1":  "#c7e9c0",
    "asm_nc2":  "#74c476",
    "asm_nc3":  "#31a354",
    "asm_nc4p": "#006d2c",
    "lr_rescued": "#fd8d3c",
    "fp_noise":   "#bdbdbd",
}

# ── helpers ───────────────────────────────────────────────────────────────────

def parse_tsv(path, min_ad=1):
    calls = {}
    with open(path) as fh:
        for line in fh:
            p = line.rstrip().split("\t")
            if len(p) < 5:
                continue
            pos = int(p[0]); ref = p[1]
            alts = [a for a in p[2].split(",") if a and a != "<*>"]
            try:
                dp = int(p[3])
            except ValueError:
                continue
            if dp == 0:
                continue
            ads = p[4].split(",")
            for i, alt in enumerate(alts):
                try:
                    ad = int(ads[i + 1])
                except (IndexError, ValueError):
                    continue
                if ad < min_ad:
                    continue
                calls[(pos, ref, alt)] = (ad, ad / dp)
    return calls

def lr_conf_set(calls, ad_min, vaf_min):    return {k for k, (ad, vaf) in calls.items() if ad >= ad_min and vaf >= vaf_min}

# ── load assembly GT ──────────────────────────────────────────────────────────

con = sqlite3.connect(DB)
gt_rows = con.execute("""
    SELECT a.sample_id, v.consensus_pos, v.ref, v.alt, COUNT(*) nc
    FROM variant v JOIN copy c ON v.copy_id=c.copy_id
    JOIN haplotype h ON c.haplotype_id=h.haplotype_id
    JOIN assembly a ON h.assembly_id=a.assembly_id
    WHERE a.cohort=? AND v.alignment_source='gene_unit_t2t'
    GROUP BY a.sample_id, v.consensus_pos, v.ref, v.alt
""", (COHORT,)).fetchall()
con.close()

gt_map = {}
for sid, pos, ref, alt, nc in gt_rows:
    pos = int(pos)
    if pos not in BOUNDARY_EXCL:
        gt_map.setdefault(sid, {})[(pos, ref, alt)] = int(nc)

# ── data collection ───────────────────────────────────────────────────────────

# Panel A / D
nc_sr_vafs   = {s[0]: [] for s in NC_STRATA}
strata_det   = {s[0]: np.zeros(len(VAF_THRESHOLDS), dtype=int) for s in NC_STRATA}
strata_total = {s[0]: 0 for s in NC_STRATA}
# per-donor lists for D error bars: sens value per donor (only donors with ≥1 GT in stratum)
strata_sens_per_donor      = {s[0]: [] for s in NC_STRATA}
# D stacked: HiFi-confirmed fraction of detected variants, per donor
strata_hifi_sens_per_donor = {s[0]: [] for s in NC_STRATA}

# Panels B, C, E, F
hist_vafs   = {"asm_multi": [], "asm_sing": [], "lr_rescued": [], "fp_noise": []}
# Panel B fine-grained nc split
hist_vafs_B = {"asm_nc1": [], "asm_nc2": [], "asm_nc3": [], "asm_nc4p": [],
               "lr_rescued": [], "fp_noise": []}
scatter_rec = []   # (sr_vaf, hifi_vaf, group)
prec_tp      = np.zeros(len(PREC_EDGES) - 1, dtype=int)
prec_tp_asm  = np.zeros(len(PREC_EDGES) - 1, dtype=int)  # assembly GT only (no rescue)
prec_tot     = np.zeros(len(PREC_EDGES) - 1, dtype=int)
# per-donor lists for E error bars: (tp_count, tot_count) per donor per bin
prec_tp_per_donor     = [[] for _ in range(len(PREC_EDGES) - 1)]
prec_tp_asm_per_donor = [[] for _ in range(len(PREC_EDGES) - 1)]
prec_tot_per_donor    = [[] for _ in range(len(PREC_EDGES) - 1)]
pr_rec   = []
gt_n_asm = 0

n_samples = 0
for sid, gt_nc in sorted(gt_map.items()):
    if sid in EXCL_ALL: continue
    sr_tsv   = BAM / sid / f"{sid}_illumina.tsv"
    hifi_tsv = BAM / sid / f"{sid}_hifi_variants.tsv"
    if not sr_tsv.exists() or sr_tsv.stat().st_size == 0:
        continue

    sr_calls   = parse_tsv(sr_tsv,   SR_MIN_AD)
    hifi_calls = parse_tsv(hifi_tsv, 1) if (hifi_tsv.exists() and sid not in EXCL_HIFI) else {}
    hifi_conf  = lr_conf_set(hifi_calls, LR_AD_MIN, LR_VAF_MIN_HIFI)

    gt_all   = set(gt_nc.keys())
    # HiFi confirmation requires SR/HiFi VAF within FOLD_MAX (3x) — drop VAF-discordant co-calls
    def _fold_ok(k):
        if k not in sr_calls or k not in hifi_calls: return False
        hv = hifi_calls[k][1]
        return hv > 0 and (1.0/FOLD_MAX) <= sr_calls[k][1]/hv <= FOLD_MAX
    lr_any   = {k for k in (hifi_conf - gt_all) if _fold_ok(k)}
    tp_hifi = gt_all | lr_any
    gt_n_asm += len(gt_all)

    # per-donor precision accumulators
    donor_prec_tp      = np.zeros(len(PREC_EDGES) - 1, dtype=int)
    donor_prec_tp_asm  = np.zeros(len(PREC_EDGES) - 1, dtype=int)
    donor_prec_tot     = np.zeros(len(PREC_EDGES) - 1, dtype=int)

    for key, (ad, vaf) in sr_calls.items():
        if key in gt_all:
            grp = "asm_multi" if gt_nc[key] >= 2 else "asm_sing"
        elif key in lr_any:
            grp = "lr_rescued"
        else:
            grp = "fp_noise"

        hist_vafs[grp].append(vaf)
        if key in gt_all:
            _nc = gt_nc[key]
            grp_B = "asm_nc1" if _nc == 1 else "asm_nc2" if _nc == 2 else "asm_nc3" if _nc == 3 else "asm_nc4p"
        else:
            grp_B = grp   # "lr_rescued" or "fp_noise"
        hist_vafs_B[grp_B].append(vaf)
        hifi_vaf = hifi_calls[key][1] if key in hifi_calls else 0.0
        scatter_rec.append((vaf, hifi_vaf, grp))

        bi = min(np.searchsorted(PREC_EDGES[1:], vaf), len(prec_tot) - 1)
        prec_tp[bi]         += 1 if key in tp_hifi else 0
        prec_tp_asm[bi]     += 1 if key in gt_all  else 0
        prec_tot[bi]        += 1
        donor_prec_tp[bi]     += 1 if key in tp_hifi else 0
        donor_prec_tp_asm[bi] += 1 if key in gt_all  else 0
        donor_prec_tot[bi]    += 1

        pr_rec.append({"vaf": vaf, "is_asm": key in gt_all, "is_anylr": key in tp_hifi})

    for bi in range(len(PREC_EDGES) - 1):
        prec_tp_per_donor[bi].append(donor_prec_tp[bi])
        prec_tp_asm_per_donor[bi].append(donor_prec_tp_asm[bi])
        prec_tot_per_donor[bi].append(donor_prec_tot[bi])

    # per-donor sensitivity per stratum
    donor_det      = {s[0]: 0 for s in NC_STRATA}
    donor_det_hifi = {s[0]: 0 for s in NC_STRATA}
    donor_total    = {s[0]: 0 for s in NC_STRATA}
    for key, nc in gt_nc.items():
        strat = next((s[0] for s in NC_STRATA if s[1] <= nc <= s[2]), None)
        if strat is None:
            continue
        strata_total[strat] += 1
        donor_total[strat]  += 1
        if key in sr_calls:
            _, vaf = sr_calls[key]
            nc_sr_vafs[strat].append(vaf)
            strata_det[strat] += (VAF_THRESHOLDS <= vaf).astype(int)
            if vaf >= VAF_CALL:
                donor_det[strat] += 1
                if key in hifi_conf and _fold_ok(key):   # fold-concordant HiFi confirmation
                    donor_det_hifi[strat] += 1

    for s in NC_STRATA:
        if donor_total[s[0]] > 0:
            strata_sens_per_donor[s[0]].append(
                donor_det[s[0]] / donor_total[s[0]] * 100)
            strata_hifi_sens_per_donor[s[0]].append(
                donor_det_hifi[s[0]] / donor_total[s[0]] * 100)

    n_samples += 1

print(f"Samples: {n_samples}")
for g, v in hist_vafs.items():
    print(f"  {g}: {len(v):,}")

for g in hist_vafs:
    hist_vafs[g] = np.array(hist_vafs[g])
sc_df    = pd.DataFrame(scatter_rec, columns=["sr_vaf", "hifi_vaf", "grp"])

# ── DATA TABLES into the figure folder ────────────────────────────────────────
_DD = OUTDIR / "data"; _DD.mkdir(parents=True, exist_ok=True)
sc_df.assign(sr_vaf_pct=100*sc_df.sr_vaf, hifi_vaf_pct=100*sc_df.hifi_vaf)[
    ["sr_vaf_pct","hifi_vaf_pct","grp"]].to_csv(_DD/"Fig74_calls_by_group.tsv",sep="\t",index=False)
pd.DataFrame({"vaf_bin_low_pct":[100*PREC_EDGES[i] for i in range(len(PREC_EDGES)-1)],
   "vaf_bin_high_pct":[100*PREC_EDGES[i+1] for i in range(len(PREC_EDGES)-1)],
   "n_calls":prec_tot,"n_TP_anyLR":prec_tp,"n_TP_assembly":prec_tp_asm,
   "precision_anyLR_pct":np.where(prec_tot>0,100*prec_tp/prec_tot,np.nan),
   "precision_assembly_pct":np.where(prec_tot>0,100*prec_tp_asm/prec_tot,np.nan)}).to_csv(
   _DD/"Fig74_precision_by_vaf.tsv",sep="\t",index=False)
pd.DataFrame([{g:int((sc_df.grp==g).sum()) for g in ["asm_multi","asm_sing","lr_rescued","fp_noise"]}]).to_csv(
   _DD/"Fig74_group_counts.tsv",sep="\t",index=False)
# stacked histogram of short-read calls by SR-VAF, true-pos split by assembly CN (1/2/>=3)
_cmap={"asm_nc1":"tp_cn1","asm_nc2":"tp_cn2","asm_nc3":"tp_cn3plus","asm_nc4p":"tp_cn3plus",
       "lr_rescued":"rescued","fp_noise":"false_pos"}
_crows=[]
for _g,_lbl in _cmap.items():
    for _v in hist_vafs_B[_g]:
        _crows.append({"sr_vaf_pct":100*_v,"group":_lbl})
pd.DataFrame(_crows).to_csv(_DD/"Fig74_calls_by_group_CN.tsv",sep="\t",index=False)
# sensitivity by CN stratum (panel D): per-donor mean +/- SD, HiFi-confirmed split
_sens=[]
for s in NC_STRATA:
    v=np.array(strata_sens_per_donor[s[0]]); h=np.array(strata_hifi_sens_per_donor[s[0]])
    if len(v)==0: continue
    _sens.append(dict(nc_stratum=s[0], n_donors=len(v),
        sensitivity_mean_pct=round(float(v.mean()),1), sensitivity_sd=round(float(v.std()),1),
        sensitivity_sem=round(float(v.std()/np.sqrt(len(v))),1),
        hifi_confirmed_mean_pct=round(float(h.mean()),1),
        sr_only_mean_pct=round(float(v.mean()-h.mean()),1)))
pd.DataFrame(_sens).to_csv(_DD/"Fig74_sensitivity_by_CN.tsv",sep="\t",index=False)
print("74 data tables ->", _DD)
pr_df    = pd.DataFrame(pr_rec)
prec_bin = np.where(prec_tot > 0, prec_tp / prec_tot * 100, np.nan)

# PR curve
pr_rows = []
for vaf_min in VAF_THRESHOLDS:
    called = pr_df[pr_df.vaf >= vaf_min]
    if len(called) == 0:
        continue
    tp_asm = int(called.is_asm.sum())
    tp     = int(called.is_anylr.sum())
    sens   = tp_asm / gt_n_asm
    prec   = tp / len(called)
    f1     = 2 * sens * prec / (sens + prec) if (sens + prec) > 0 else 0.0
    pr_rows.append({"vaf": vaf_min, "sens": sens, "prec": prec, "f1": f1})
pr_curve = pd.DataFrame(pr_rows)
idx_call = int((pr_curve.vaf - VAF_CALL).abs().idxmin())
pt_call  = pr_curve.loc[idx_call]

# Log bins for Panel B
LOG_EDGES = np.logspace(np.log10(0.0003), np.log10(0.25), 38)  # 37 bins, 0.03%–25%
N_BINS = len(LOG_EDGES) - 1
idx_thresh_B = int(np.searchsorted(LOG_EDGES[:-1], VAF_CALL))
tick_vafs_pct = [0.05, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 5.0, 15.0]
tick_pos_B = [int(np.searchsorted(LOG_EDGES[:-1], v / 100)) for v in tick_vafs_pct]

# ═══════════════════════════════════════════════════════════════════════════════
# Figure
# ═══════════════════════════════════════════════════════════════════════════════

fig, axes = plt.subplots(2, 3, figsize=(7.2, 6.0))
fig.subplots_adjust(left=0.09, right=0.97, top=0.88, bottom=0.22,
                    hspace=0.55, wspace=0.38)
fig.suptitle(
    "5S rDNA SR variant calling — HPRC Year 1 (38 samples)\n"
    f"SR threshold: AD≥{SR_MIN_AD} & VAF≥{VAF_CALL*100:.2g}%  ·  "
    f"HiFi confirmation: AD≥{LR_AD_MIN} & VAF≥{LR_VAF_MIN_HIFI*100:.1f}%",
    fontsize=9, y=1.01
)
(ax_A, ax_B, ax_C), (ax_D, ax_E, ax_F) = axes

# ── Panel A: copy count → SR VAF (violin) ────────────────────────────────────
strat_ok_A = [s[0] for s in NC_STRATA
              if strata_total[s[0]] > 0 and len(nc_sr_vafs[s[0]]) >= 3]
vaf_data_A = [np.array(nc_sr_vafs[lbl]) * 100 for lbl in strat_ok_A]

vp = ax_A.violinplot(vaf_data_A, positions=range(len(strat_ok_A)),
                     showmedians=True, showextrema=False, widths=0.75)
for body in vp["bodies"]:
    body.set_facecolor("#74c476")
    body.set_edgecolor("#2ca25f")
    body.set_linewidth(0.6)
    body.set_alpha(0.65)
vp["cmedians"].set_color("#2ca25f")
vp["cmedians"].set_linewidth(2.0)

ax_A.axhline(VAF_CALL * 100, color="black", lw=1.5, ls="--", alpha=0.8,
             label=f"0.30% reporting threshold")
ax_A.set_yscale("log")
ax_A.set_ylim(0.06, 35)
ax_A.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:.2g}%"))
ax_A.set_xticks(range(len(strat_ok_A)))
ax_A.set_xticklabels(strat_ok_A, rotation=38, ha="right")
ax_A.set_xlabel("Assembly copy-count stratum")
ax_A.set_ylabel("Observed SR VAF (log scale)")
ax_A.set_title("A  SR detection VAF scales with copy count\n"
               "(GT variants detected by SR; n per violin = detected count)",
               loc="left")
ax_A.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2)
ax_A.grid(True, axis="y", lw=0.3, alpha=0.4, which="both")

# Annotate detection rate above each violin
idx_call_A = int(np.argmin(np.abs(VAF_THRESHOLDS - VAF_CALL)))
for i, lbl in enumerate(strat_ok_A):
    n   = strata_total[lbl]
    det = strata_det[lbl][idx_call_A]
    pct = det / n * 100
    col = "#2ca25f" if pct >= 90 else "#e6550d"
    ax_A.text(i, 22, f"{pct:.0f}%", ha="center", va="bottom",
              fontsize=6.5, color=col, fontweight="bold")
ax_A.text(-0.6, 22, "det.\n@0.3%:", ha="right", va="bottom", fontsize=6.5, color="grey")

# ── Panel B: stacked call landscape (log-binned histogram) ───────────────────
grp_order_B  = ["asm_nc1", "asm_nc2", "asm_nc3", "asm_nc4p", "lr_rescued", "fp_noise"]
grp_labels_B = ["Asm GT  nc=1", "Asm GT  nc=2", "Asm GT  nc=3",
                "Asm GT  nc≥4", "HiFi-rescued", "FP — no LR support"]
bottom_B   = np.zeros(N_BINS, dtype=float)

for grp, lbl in zip(grp_order_B, grp_labels_B):
    counts, _ = np.histogram(hist_vafs_B[grp], bins=LOG_EDGES)
    ax_B.bar(range(N_BINS), counts, bottom=bottom_B, color=COLORS_B[grp],
             label=f"{lbl}  (n={len(hist_vafs_B[grp]):,})",
             edgecolor="white", linewidth=0.25, width=1.0)
    bottom_B += counts

ax_B.axvline(idx_thresh_B - 0.5, color="black", lw=1.8, ls="--", alpha=0.85,
             label="0.30% threshold")
ax_B.set_xticks(tick_pos_B)
ax_B.set_xticklabels([f"{v:.2g}%" for v in tick_vafs_pct],
                     rotation=40, ha="right")
ax_B.set_xlabel("SR observed VAF")
ax_B.set_ylabel("Number of SR calls")
ax_B.set_title("B  Call landscape across all VAF levels\n"
               "(all SR calls from 38 samples; stacked by classification)",
               loc="left")
ax_B.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2)
ax_B.grid(True, axis="y", lw=0.3, alpha=0.4)

# Shade below-threshold region
ax_B.axvspan(0, idx_thresh_B - 0.5, color="grey", alpha=0.07, zorder=0)
ax_B.text(idx_thresh_B / 2, ax_B.get_ylim()[1] * 0.92, "below\nthreshold",
          ha="center", va="top", fontsize=6.5, color="grey", style="italic")

# ── Panel C: SR VAF vs HiFi VAF scatter (cross-platform validation) ───────────
EPS = 0.00008   # floor for log scale — "no HiFi support"
rng = np.random.default_rng(42)

for grp, col, zorder, alpha, s in [
    ("fp_noise",   COLORS["fp_noise"],   1, 0.12,  2.5),
    ("lr_rescued", COLORS["lr_rescued"], 3, 0.55,  5.5),
    ("asm_sing",   COLORS["asm_sing"],   4, 0.65,  5.5),
    ("asm_multi",  COLORS["asm_multi"],  4, 0.65,  5.5),
]:
    sub = sc_df[sc_df.grp == grp].copy()
    if grp == "fp_noise" and len(sub) > 2000:
        sub = sub.sample(2000, random_state=42)
    x = sub.sr_vaf.values
    y = sub.hifi_vaf.values
    # Zero-HiFi calls → jitter near EPS for log-scale visibility
    y_plot = np.where(y < EPS, EPS * rng.uniform(0.4, 1.0, len(y)), y)
    n_full = len(sc_df[sc_df.grp == grp])
    lbl = (f"Asm GT nc=1 (n={n_full:,})"    if grp == "asm_sing"   else
           f"Asm GT nc≥2 (n={n_full:,})"    if grp == "asm_multi"  else
           f"LR-rescued (n={n_full:,})"      if grp == "lr_rescued" else
           f"FP noise (n={n_full:,}, subsampled)")
    ax_C.scatter(x * 100, y_plot * 100, s=s, color=col, alpha=alpha,
                 zorder=zorder, label=lbl, linewidths=0)

ax_C.set_xscale("log"); ax_C.set_yscale("log")
ax_C.set_xlim(0.04, 150); ax_C.set_ylim(0.006, 150)
ax_C.axvline(VAF_CALL * 100, color="black", lw=1.2, ls="--", alpha=0.7,
             label=f"SR 0.30%")
ax_C.axhline(LR_VAF_MIN_HIFI * 100, color="#2ca25f", lw=1.2, ls=":",
             alpha=0.85, label=f"HiFi {LR_VAF_MIN_HIFI*100:.1f}% conf. threshold")
# y=x reference line
ref = np.logspace(np.log10(0.04), np.log10(150), 50)
ax_C.plot(ref, ref, color="grey", lw=0.8, ls="-", alpha=0.4, zorder=0)
ax_C.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:.2g}%"))
ax_C.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:.2g}%"))
ax_C.set_xlabel("SR VAF (%)")
ax_C.set_ylabel("HiFi VAF (%)")
ax_C.set_title("C  Cross-platform concordance\n"
               "(SR vs HiFi VAF per call; diagonal = perfect agreement)",
               loc="left")
ax_C.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), markerscale=2.5)
ax_C.grid(True, lw=0.3, alpha=0.3)
# Annotate the "no HiFi support" floor
ax_C.text(0.3, 0.01, "← no HiFi support", ha="left", va="center", fontsize=6.5,
          color="grey", style="italic", transform=ax_C.get_xaxis_transform())

# ── Panel D: sensitivity by copy count (stacked: HiFi-confirmed + SR-only) ────
strat_ok_D = [s[0] for s in NC_STRATA if strata_total[s[0]] > 0]
strat_ns_D = [strata_total[lbl] for lbl in strat_ok_D]

from scipy import stats as _stats
sens_D_vals      = [np.array(strata_sens_per_donor[lbl])      for lbl in strat_ok_D]
sens_D_hifi_vals = [np.array(strata_hifi_sens_per_donor[lbl]) for lbl in strat_ok_D]
sens_D_mean      = np.array([v.mean() for v in sens_D_vals])
sens_D_hifi_mean = np.array([v.mean() if len(v) else 0.0 for v in sens_D_hifi_vals])
sens_D_sronly_mean = sens_D_mean - sens_D_hifi_mean
sens_D_ci        = np.array([
    _stats.t.ppf(0.975, len(v) - 1) * v.std(ddof=1) / np.sqrt(len(v))
    if len(v) > 1 else 0.0 for v in sens_D_vals])

xs = np.arange(len(strat_ok_D))
# Bottom: HiFi-confirmed fraction
ax_D.bar(xs, sens_D_hifi_mean,   color="#2ca25f", edgecolor="white", linewidth=0.5,
         label="SR + HiFi confirmed")
# Top: SR-detected but not HiFi-confirmed
ax_D.bar(xs, sens_D_sronly_mean, bottom=sens_D_hifi_mean,
         color="#74c476", edgecolor="white", linewidth=0.5,
         label="SR only (no HiFi)")
ax_D.errorbar(xs, sens_D_mean, yerr=sens_D_ci,
              fmt="none", color="black", capsize=4, capthick=1.2, lw=1.2, zorder=5)

for j, (n, s, ci) in enumerate(zip(strat_ns_D, sens_D_mean, sens_D_ci)):
    ax_D.text(j, 2, f"n={n:,}", ha="center", va="bottom", fontsize=6.5,
              color="white" if s > 30 else "black", fontweight="bold")
    ax_D.text(j, s + ci + 1.5, f"{s:.0f}%", ha="center", va="bottom")

ax_D.axhline(100, color="grey", lw=0.8, ls=":", alpha=0.5)
ax_D.set_ylim(0, 118)
ax_D.set_xticks(xs)
ax_D.set_xticklabels(strat_ok_D, rotation=38, ha="right")
ax_D.set_xlabel("Assembly copy-count stratum")
ax_D.set_ylabel("Sensitivity at VAF ≥ 0.30% (%)")
ax_D.set_title("D  Detection sensitivity by copy count\n"
               "(mean ± 95% CI; stacked: HiFi-confirmed vs SR-only detections)",
               loc="left")
ax_D.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2)
ax_D.grid(True, axis="y", lw=0.3, alpha=0.4)

# ── Panel E: precision by VAF bin (stacked: assembly GT + HiFi-rescued) ────────
show_idx = list(range(1, len(PREC_LABELS)))
tot_E    = prec_tot[show_idx]
lbl_E    = [PREC_LABELS[i] for i in show_idx]

# per-donor precision: total (asm GT + HiFi-rescued) and assembly-only components
prec_E_mean     = np.full(len(show_idx), np.nan)
prec_E_asm_mean = np.full(len(show_idx), np.nan)
prec_E_ci       = np.full(len(show_idx), np.nan)
for j, i in enumerate(show_idx):
    tp_arr      = np.array(prec_tp_per_donor[i])
    tp_asm_arr  = np.array(prec_tp_asm_per_donor[i])
    tot_arr     = np.array(prec_tot_per_donor[i])
    mask = tot_arr > 0
    if mask.sum() == 0:
        continue
    vals     = tp_arr[mask]     / tot_arr[mask] * 100
    vals_asm = tp_asm_arr[mask] / tot_arr[mask] * 100
    prec_E_mean[j]     = vals.mean()
    prec_E_asm_mean[j] = vals_asm.mean()
    prec_E_ci[j] = (_stats.t.ppf(0.975, len(vals) - 1) * vals.std(ddof=1) / np.sqrt(len(vals))
                    if len(vals) > 1 else 0.0)

prec_E_rescue_mean = np.where(
    np.isnan(prec_E_mean) | np.isnan(prec_E_asm_mean),
    0.0, prec_E_mean - prec_E_asm_mean)
prec_E_asm_plot = np.where(np.isnan(prec_E_asm_mean), 0.0, prec_E_asm_mean)

xs_E = np.arange(len(show_idx))
# Bottom: assembly GT precision (without HiFi rescue)
ax_E.bar(xs_E, prec_E_asm_plot,    color="#2ca25f", edgecolor="white", linewidth=0.5,
         label="Assembly GT TP")
# Top: additional from HiFi rescue
ax_E.bar(xs_E, prec_E_rescue_mean, bottom=prec_E_asm_plot,
         color="#fd8d3c", edgecolor="white", linewidth=0.5,
         label="HiFi-rescued TP")
ax_E.errorbar(xs_E, prec_E_mean, yerr=prec_E_ci,
              fmt="none", color="black", capsize=4, capthick=1.2, lw=1.2, zorder=5)

for j, (p, se, n) in enumerate(zip(prec_E_mean, prec_E_ci, tot_E)):
    if np.isnan(p) or n == 0:
        continue
    ax_E.text(j, p + se + 1.5, f"n={n:,}", ha="center", va="bottom",
              fontsize=6.5, rotation=45)

ax_E.axhline(80, color="grey", lw=0.9, ls=":", alpha=0.55, label="80%")
ax_E.axhline(50, color="grey", lw=0.6, ls=":", alpha=0.4,  label="50%")
cut_E = next(j for j, i in enumerate(show_idx) if PREC_EDGES[i] >= VAF_CALL) - 0.5
ax_E.axvline(cut_E, color="black", lw=1.5, ls="--", alpha=0.8, label="0.30% cut")
ax_E.set_xticks(xs_E)
ax_E.set_xticklabels(lbl_E, rotation=45, ha="right")
ax_E.set_xlabel("SR observed VAF bin")
ax_E.set_ylabel("Precision — % of calls that are TP")
ax_E.set_title("E  Precision by SR VAF bin\n"
               "(mean ± 95% CI; stacked: assembly GT TP vs HiFi-rescued TP)",
               loc="left")
ax_E.set_ylim(0, 118)
ax_E.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2)
ax_E.grid(True, axis="y", lw=0.3, alpha=0.4)
ax_E.axvspan(-0.5, cut_E, color="grey", alpha=0.07, zorder=0)

# ── Panel F: PR curve ─────────────────────────────────────────────────────────
ax_F.plot(pr_curve.prec * 100, pr_curve.sens * 100,
          color="#d62728", lw=2.2, zorder=2, label="HiFi GT")

# Mark key operating points
for vaf_pct, ls, mk in [(0.20, ".", "o"), (0.30, "--", "D"),
                         (0.50, ".", "o"), (1.00, ".", "o")]:
    r = pr_curve.iloc[(pr_curve.vaf - vaf_pct/100).abs().argsort().values[0]]
    s_val = 50 if vaf_pct == 0.30 else 20
    ax_F.scatter([r.prec * 100], [r.sens * 100], color="black",
                 s=s_val, zorder=5, marker=mk)
    ax_F.annotate(f"  {vaf_pct:.2g}%",
                  (r.prec * 100, r.sens * 100),
                  textcoords="offset points", xytext=(4, 2),
                  fontsize=6.5, color="black" if vaf_pct == 0.30 else "grey")

# Summary box
summary = (f"At 0.30% VAF:\n"
           f"  Sensitivity = {pt_call.sens*100:.1f}%\n"
           f"  Precision   = {pt_call.prec*100:.1f}%\n"
           f"  F1          = {pt_call.f1:.3f}")
ax_F.text(0.04, 0.07, summary, transform=ax_F.transAxes, fontsize=6.5,
          verticalalignment="bottom",
          bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                    edgecolor="#aaaaaa", alpha=0.9))

ax_F.set_xlabel("Precision (%)")
ax_F.set_ylabel("Sensitivity (%)")
ax_F.set_title("F  Precision-Recall curve\n"
               "(sensitivity over assembly GT denominator)",
               loc="left")
ax_F.set_xlim(28, 102)
ax_F.set_ylim(18, 95)
ax_F.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2)
ax_F.grid(True, lw=0.3, alpha=0.4)

# ── save ──────────────────────────────────────────────────────────────────────
out = OUTDIR / "74_summary_figure.pdf"
plt.savefig(out, dpi=300, bbox_inches="tight")
plt.close()
print(f"Figure → {out}")
