#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 75_platform_comparison.py — cross-platform (short-read vs assembly vs HiFi)
# VAF concordance figure for 5S rDNA variant calls.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
75_platform_comparison.py

1 × 3 cross-platform VAF concordance figure (HiFi only; ONT dropped).

A  SR vs Assembly expected VAF  (coloured by copy count)
B  SR vs HiFi — assembly GT variants  (coloured by copy count)
C  SR vs HiFi — HiFi-rescued variants  (coloured by HiFi confirmation)
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
BOUNDARY_EXCL   = set(range(1, 90)) | set(range(2034, 2169))
EXCL_ALL  = {"HG02486"}   # no usable HiFi or SR data
EXCL_HIFI = {"HG02818"}   # HiFi alignment failed (48× SR/HiFi ratio); SR kept
SR_MIN_AD       = 3
LR_AD_MIN       = 5
LR_VAF_MIN_HIFI = 0.003   # 0.3%
VAF_CALL        = 0.003

EPS = 0.00008   # plotting floor for "no reads" on log scale

# (label, lo, hi, green_col, gray_col)
NC_CATS = [
    ("nc=1",  1,    1, "#c7e9c0", "#e0e0e0"),
    ("nc=2",  2,    2, "#74c476", "#bdbdbd"),
    ("nc=3",  3,    3, "#31a354", "#737373"),
    ("nc≥4",  4, 9999, "#006d2c", "#252525"),
]

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

def lr_conf_set(calls, ad_min, vaf_min):
    return {k for k, (ad, vaf) in calls.items() if ad >= ad_min and vaf >= vaf_min}

# ── load GT and total copies ──────────────────────────────────────────────────

con = sqlite3.connect(DB)
gt_rows = con.execute("""
    SELECT a.sample_id, v.consensus_pos, v.ref, v.alt, COUNT(*) nc
    FROM variant v JOIN copy c ON v.copy_id=c.copy_id
    JOIN haplotype h ON c.haplotype_id=h.haplotype_id
    JOIN assembly a ON h.assembly_id=a.assembly_id
    WHERE a.cohort=? AND v.alignment_source='gene_unit_t2t'
    GROUP BY a.sample_id, v.consensus_pos, v.ref, v.alt
""", (COHORT,)).fetchall()

n_total_rows = con.execute("""
    SELECT a.sample_id, COUNT(*) as n_copies
    FROM copy c
    JOIN haplotype h ON c.haplotype_id=h.haplotype_id
    JOIN assembly a ON h.assembly_id=a.assembly_id
    WHERE a.cohort=?
    GROUP BY a.sample_id
""", (COHORT,)).fetchall()
con.close()

gt_map = {}
for sid, pos, ref, alt, nc in gt_rows:
    pos = int(pos)
    if pos not in BOUNDARY_EXCL:
        gt_map.setdefault(sid, {})[(pos, ref, alt)] = int(nc)

n_total_map = {sid: int(n) for sid, n in n_total_rows}

# ── data collection ───────────────────────────────────────────────────────────

# Assembly GT × SR concordant (not HiFi confirmed): green nc-graded in panel A
# Assembly GT × SR × HiFi confirmed: orange in panel A
# cols: sr_vaf, asm_vaf, hifi_vaf, nc_cat, nc_col, hifi_ok
asm_rec = []

# FN: assembly GT variants missed by SR — col: asm_vaf (fraction)
fn_rec  = []

# FP: SR calls not in assembly GT and not HiFi-confirmed — col: sr_vaf (fraction)
fp_rec  = []

# HiFi-rescued SR calls (not in assembly GT)
# cols: sr_vaf, hifi_vaf
lr_rec = []

# SR-FP ∩ HiFi-conf-FP overlap enrichment vs chance
FOLD_MAX = 3   # SR and HiFi VAFs must be within this many-fold to count as confirmed
NONBOUND_SITES = (2168 - len(BOUNDARY_EXCL)) * 3   # possible non-boundary SNV (pos×alt) space
enr_obs = 0.0; enr_exp = 0.0

n_samples = 0
for sid, gt_nc in sorted(gt_map.items()):
    if sid in EXCL_ALL: continue
    sr_tsv   = BAM / sid / f"{sid}_illumina.tsv"
    hifi_tsv = BAM / sid / f"{sid}_hifi_variants.tsv"
    if not sr_tsv.exists() or sr_tsv.stat().st_size == 0:
        continue

    sr_calls   = {k: v for k, v in parse_tsv(sr_tsv, SR_MIN_AD).items()
                  if k[0] not in BOUNDARY_EXCL}
    hifi_calls = {k: v for k, v in (parse_tsv(hifi_tsv, 1) if (hifi_tsv.exists() and sid not in EXCL_HIFI) else {}).items()
                  if k[0] not in BOUNDARY_EXCL}
    hifi_conf  = lr_conf_set(hifi_calls, LR_AD_MIN, LR_VAF_MIN_HIFI)

    gt_all = set(gt_nc.keys())
    lr_any = hifi_conf - gt_all
    n_tot  = n_total_map.get(sid)

    # A valid HiFi confirmation requires HiFi-conf AND SR/HiFi VAF within FOLD_MAX (3x).
    def _conf(key):
        if key not in hifi_conf or key not in sr_calls or key not in hifi_calls:
            return False
        hv = hifi_calls[key][1]
        if hv <= 0:
            return False
        r = sr_calls[key][1] / hv
        return (1.0 / FOLD_MAX) <= r <= FOLD_MAX

    # enrichment: observed = fold-CONCORDANT SR-FP ∩ HiFi-FP co-calls; exp = random co-occurrence
    # over the possible non-boundary SNV site space (platform-independent null).
    sr_fp_set   = set(sr_calls) - gt_all
    hifi_fp_set = hifi_conf - gt_all
    enr_obs += sum(1 for k in (sr_fp_set & hifi_fp_set) if _conf(k))
    enr_exp += len(sr_fp_set) * len(hifi_fp_set) / NONBOUND_SITES

    # Assembly GT variants — concordant or FN; hifi_ok requires fold-concordance
    if n_tot:
        for key, nc in gt_nc.items():
            cat = next((c for c in NC_CATS if c[1] <= nc <= c[2]), NC_CATS[-1])
            if key in sr_calls:
                asm_rec.append((
                    sr_calls[key][1], nc / n_tot,
                    hifi_calls[key][1] if key in hifi_calls else 0.0,
                    cat[0], cat[3],
                    _conf(key)                                         # hifi_ok (fold-concordant)
                ))
            else:
                fn_rec.append((nc / n_tot, cat[0]))

    # FP: SR calls not in assembly and NOT a fold-concordant HiFi confirmation
    # (VAF-discordant HiFi co-calls now count as FP, not rescued)
    for key, (ad, vaf) in sr_calls.items():
        if key not in gt_all and not _conf(key) and vaf >= VAF_CALL:
            fp_rec.append(vaf)

    # HiFi-rescued: ALL SR-FP co-called by HiFi (panel B colours in-fold vs out-fold);
    # the COUNTED rescue (precision/enrichment) is the fold-concordant subset only.
    for key in lr_any:
        if key not in sr_calls:
            continue
        lr_rec.append((sr_calls[key][1],
                       hifi_calls[key][1] if key in hifi_calls else 0.0))

    n_samples += 1

n_hifi_ok = sum(r[5] for r in asm_rec)
print(f"Samples: {n_samples}")
print(f"Asm GT concordant (panels A+B): {len(asm_rec):,}  "
      f"(HiFi confirmed: {n_hifi_ok:,})")
print(f"FN (assembly GT, SR missed): {len(fn_rec):,}")
print(f"FP (SR only, no assembly/HiFi): {len(fp_rec):,}")
print(f"HiFi-rescued (panel C): {len(lr_rec):,}")

asm_df = pd.DataFrame(asm_rec,
                      columns=["sr_vaf", "asm_vaf", "hifi_vaf", "nc_cat", "nc_col", "hifi_ok"])
lr_df  = pd.DataFrame(lr_rec, columns=["sr_vaf", "hifi_vaf"])

# ── DATA TABLES into the figure folder (VAF as %, matching the axes) ───────────
DATADIR = OUTDIR / "data"; DATADIR.mkdir(parents=True, exist_ok=True)
# Panel A: SR observed vs Assembly expected VAF, classified TP / FN / FP
pa = [dict(status="TP", sr_observed_vaf_pct=100*r[0], assembly_expected_vaf_pct=100*r[1],
           hifi_vaf_pct=100*r[2], nc_category=r[3], hifi_confirmed=int(r[5])) for r in asm_rec]
pa += [dict(status="FN", sr_observed_vaf_pct=None, assembly_expected_vaf_pct=100*av,
            hifi_vaf_pct=None, nc_category=cat, hifi_confirmed=None) for av, cat in fn_rec]
pa += [dict(status="FP", sr_observed_vaf_pct=100*v, assembly_expected_vaf_pct=None,
            hifi_vaf_pct=None, nc_category=None, hifi_confirmed=None) for v in fp_rec]
pd.DataFrame(pa).to_csv(DATADIR/"Fig2A_SR_vs_assembly_VAF.tsv", sep="\t", index=False)
# Panel B: SR observed vs HiFi observed VAF for HiFi-rescued calls (assembly errors)
(lr_df.assign(sr_observed_vaf_pct=100*lr_df.sr_vaf, hifi_observed_vaf_pct=100*lr_df.hifi_vaf)
   [["sr_observed_vaf_pct","hifi_observed_vaf_pct"]]
   .to_csv(DATADIR/"Fig2B_SR_vs_HiFi_rescued.tsv", sep="\t", index=False))
# SR∩HiFi FP-overlap enrichment vs chance + Poisson p-value (for panel B annotation)
import math as _math
from scipy import stats as _stats
ENR = enr_obs / enr_exp if enr_exp else float("nan")
_log10p = _stats.poisson.logsf(enr_obs - 1, enr_exp) / _math.log(10) if enr_exp else float("nan")
ENR_P_TXT = ("p < 10$^{-300}$" if _log10p < -300 else f"p = 10$^{{{_log10p:.0f}}}$")
print(f"SR∩HiFi FP overlap (fold-concordant): obs={enr_obs:.0f} exp={enr_exp:.1f}  "
      f"enrichment={ENR:.1f}x  log10p={_log10p:.0f}")
print(f"Rescued: co-called={len(lr_rec):,}  fold-concordant(counted)={int(enr_obs):,}")
# headline counts
pd.DataFrame([dict(n_TP=len(asm_rec), n_TP_hifi_confirmed_fold=int(n_hifi_ok), n_FN=len(fn_rec),
   n_FP=len(fp_rec),
   n_rescued_cocalled=len(lr_rec), n_rescued_fold_concordant=int(enr_obs),
   sensitivity_pct=round(100*len(asm_rec)/(len(asm_rec)+len(fn_rec)),1),
   precision_assembly_pct=round(100*len(asm_rec)/(len(asm_rec)+len(fp_rec)),1),
   rescue_expected_by_chance=round(enr_exp,1), rescue_enrichment_x=round(ENR,2),
   rescue_log10_p=round(_log10p,1), fold_max=FOLD_MAX,
   SR_min_AD=SR_MIN_AD, HiFi_min_AD=LR_AD_MIN, VAF_call_pct=VAF_CALL*100,
   n_samples=n_samples)]).to_csv(DATADIR/"Fig2AB_summary_counts.tsv", sep="\t", index=False)
print("data tables ->", DATADIR)

# ── shared helpers ────────────────────────────────────────────────────────────

rng = np.random.default_rng(42)

def jitter_log(arr, sigma=0.025):
    return arr * np.exp(rng.normal(0, sigma, len(arr)))

def floor_log(arr):
    return np.where(arr < EPS,
                    EPS * rng.uniform(0.3, 0.95, len(arr)),
                    arr)

def setup_ax(ax, xlabel, ylabel, title,
             x_thr=None, y_thr=None, y_thr_label=None,
             xlim=(0.04, 150), ylim=(0.006, 150)):
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(*xlim); ax.set_ylim(*ylim)
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:.3g}%"))
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:.3g}%"))
    ax.grid(True, lw=0.3, alpha=0.3)
    ref = np.array([min(xlim[0], ylim[0]), max(xlim[1], ylim[1])])
    ax.plot(ref, ref, color="grey", lw=0.9, ls="-", alpha=0.35, zorder=0,
            label="y = x")
    if x_thr is not None:
        ax.axvline(x_thr * 100, color="black", lw=1.4, ls="--", alpha=0.75,
                   label=f"SR {x_thr*100:.2g}%")
    if y_thr is not None:
        lbl = y_thr_label or f"HiFi {y_thr*100:.3g}%"
        ax.axhline(y_thr * 100, color="#d62728", lw=1.4, ls=":",
                   alpha=0.85, label=lbl)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc="left")

# ── fold-concordance threshold for panel B ────────────────────────────────────
FOLD_MAX = 3   # SR and HiFi VAFs must be within this many-fold of each other

# orange VAF bins for rescued variants (panel B)
ORANGE_BINS = [
    (VAF_CALL, 0.010, "#fdd0a2", "0.3–1%"),
    (0.010,    0.030, "#fdae6b", "1–3%"),
    (0.030,    0.100, "#fd8d3c", "3–10%"),
    (0.100,    1.010, "#d94801", ">10%"),
]

# ═══════════════════════════════════════════════════════════════════════════════
# Figure  (2 × 1 stacked)
# ═══════════════════════════════════════════════════════════════════════════════

CM = 1 / 2.54
_pw, _ph = 4 * CM, 4 * CM      # 4×4 cm panel area
_lm, _rm  = 0.90, 2.80          # left (ylabel+ticks), right (legend)
_bm, _tm  = 0.65, 0.80          # bottom (xlabel), top (suptitle + panel title)
_gap      = 0.55                 # vertical gap between panels
_fw = _lm + _pw + _rm
_fh = _bm + _ph + _gap + _ph + _tm

fig = plt.figure(figsize=(_fw, _fh))
ax_A = fig.add_axes([_lm/_fw, (_bm + _ph + _gap)/_fh, _pw/_fw, _ph/_fh])
ax_B = fig.add_axes([_lm/_fw,  _bm/_fh,               _pw/_fw, _ph/_fh])

fig.suptitle(
    "5S rDNA SR variant calling — HPRC Year 1 (41 samples)\n"
    f"SR: AD≥{SR_MIN_AD} & VAF≥{VAF_CALL*100:.2g}%  ·  "
    f"HiFi conf.: AD≥{LR_AD_MIN} & VAF≥{LR_VAF_MIN_HIFI*100:.1f}%  ·  "
    f"Fold tolerance: {FOLD_MAX}×",
    fontsize=9, y=0.98
)

# ── A — SR vs Assembly  (concordant · FN left strip · FP bottom strip) ────────
_SEP_Y  = 0.07   # y-separator between main zone and FP strip (%)
_FN_X   = 0.048  # centre x of FN strip (%) — left of SR threshold
_FP_Y   = 0.048  # centre y of FP strip (%) — below assembly VAF range

ax_A.set_xscale("log"); ax_A.set_yscale("log")
ax_A.set_xlim(0.035, 150); ax_A.set_ylim(0.035, 150)
ax_A.xaxis.set_major_formatter(
    mticker.FuncFormatter(lambda v, _: f"{v:.3g}%" if v >= 0.1 else ""))
ax_A.yaxis.set_major_formatter(
    mticker.FuncFormatter(lambda v, _: f"{v:.3g}%" if v >= 0.1 else ""))
ax_A.grid(True, lw=0.3, alpha=0.3)

_ref = np.logspace(np.log10(0.1), np.log10(150), 50)
ax_A.plot(_ref, _ref, color="grey", lw=0.9, ls="-", alpha=0.35, zorder=0,
          label="y = x")
ax_A.axvline(VAF_CALL * 100, color="black",   lw=1.4, ls="--", alpha=0.75,
             label=f"SR {VAF_CALL*100:.2g}%")
ax_A.axhline(VAF_CALL * 100, color="#d62728", lw=1.4, ls=":",  alpha=0.85,
             label="asm {:.2g}%".format(VAF_CALL * 100))
ax_A.axhline(_SEP_Y, color="#aaaaaa", lw=0.7, ls=":", alpha=0.6)

# green nc-graded: assembly GT × SR concordant
for cat_lbl, lo, hi, gcol, _gray in reversed(NC_CATS):
    mask = (asm_df.nc_cat == cat_lbl).values
    if mask.sum() == 0:
        continue
    ax_A.scatter(
        jitter_log(asm_df.sr_vaf.values[mask]  * 100, sigma=0.02),
        floor_log( asm_df.asm_vaf.values[mask]) * 100,
        s=4.5, color=gcol, alpha=0.55, linewidths=0, zorder=3,
        label=f"{cat_lbl}  (n={mask.sum():,})"
    )

# FN strip — left edge, gray shades mirroring nc gradation
if fn_rec:
    fn_vafs  = np.array([v for v, _ in fn_rec]) * 100
    fn_cats  = [c for _, c in fn_rec]
    fn_x_all = _FN_X * np.exp(rng.normal(0, 0.12, len(fn_vafs)))
    first_fn = True
    for cat_lbl, lo, hi, _gcol, gray in reversed(NC_CATS):
        idx = np.array([i for i, c in enumerate(fn_cats) if c == cat_lbl])
        if len(idx) == 0:
            continue
        ax_A.scatter(fn_x_all[idx], fn_vafs[idx],
                     s=4.5, color=gray, alpha=0.65, linewidths=0, zorder=3,
                     label=(f"FN {cat_lbl}  (n={len(idx):,})"
                            if first_fn else f"    {cat_lbl}  (n={len(idx):,})"))
        first_fn = False
    ax_A.text(_FN_X * 1.18, 110, f"FN\nn={len(fn_vafs):,}",
              ha="left", va="top", fontsize=6.5, color="#636363")

# FP strip — bottom edge, uniform gray
if fp_rec:
    fp_arr = np.array(fp_rec) * 100
    fp_y   = _FP_Y * np.exp(rng.normal(0, 0.15, len(fp_arr)))
    ax_A.scatter(fp_arr, fp_y,
                 s=4.5, color="#bdbdbd", alpha=0.55, linewidths=0, zorder=3,
                 label=f"FP  (n={len(fp_arr):,})")
    ax_A.text(0.12, _FP_Y * 1.7, f"FP  n={len(fp_arr):,}",
              ha="left", va="bottom", fontsize=6.5, color="#636363")

ax_A.set_xlabel("SR observed VAF (%)")
ax_A.set_ylabel("Assembly expected VAF (%)")
ax_A.set_title("A  SR vs Assembly  (concordant · FN · FP)", loc="left")
ax_A.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), markerscale=2.5, ncol=1)

# ── B — SR vs HiFi  (HiFi-rescued variants, fold-concordance filter) ──────────
# Only SR-called (VAF ≥ threshold) rescued variants are shown.
# Orange shades = VAF level; gray = within detection but outside fold tolerance.

lr_plot = lr_df[lr_df.sr_vaf >= VAF_CALL].copy()
with np.errstate(divide="ignore", invalid="ignore"):
    ratio = np.where(lr_plot.hifi_vaf.values > 0,
                     lr_plot.sr_vaf.values / lr_plot.hifi_vaf.values,
                     np.inf)
in_fold  = (ratio >= 1 / FOLD_MAX) & (ratio <= FOLD_MAX)
out_fold = ~in_fold

ax_B.set_xscale("log"); ax_B.set_yscale("log")
ax_B.set_xlim(0.035, 150); ax_B.set_ylim(0.035, 150)
ax_B.xaxis.set_major_formatter(
    mticker.FuncFormatter(lambda v, _: f"{v:.3g}%" if v >= 0.1 else ""))
ax_B.yaxis.set_major_formatter(
    mticker.FuncFormatter(lambda v, _: f"{v:.3g}%" if v >= 0.1 else ""))
ax_B.grid(True, lw=0.3, alpha=0.3)

# reference diagonal and fold-band lines
_ref2 = np.logspace(np.log10(0.035), np.log10(150), 60)
ax_B.plot(_ref2, _ref2, color="grey", lw=0.9, ls="-", alpha=0.35, zorder=0,
          label="y = x")
_lo = np.logspace(np.log10(0.035), np.log10(150 / FOLD_MAX), 60)
ax_B.plot(_lo,           _lo * FOLD_MAX, color="#888888", lw=0.9, ls="--",
          alpha=0.50, label=f"±{FOLD_MAX}× fold")
ax_B.plot(_lo * FOLD_MAX, _lo,           color="#888888", lw=0.9, ls="--",
          alpha=0.50)

ax_B.axvline(VAF_CALL * 100, color="black",   lw=1.4, ls="--", alpha=0.75,
             label=f"SR {VAF_CALL*100:.2g}%")
ax_B.axhline(LR_VAF_MIN_HIFI * 100, color="#d62728", lw=1.4, ls=":", alpha=0.85,
             label=f"HiFi {LR_VAF_MIN_HIFI*100:.1f}%")

# enrichment of fold-concordant SR∩HiFi rescue over chance
ax_B.text(0.035, 0.965,
    f"≤{FOLD_MAX}× rescue: {ENR:.1f}× over chance\n{ENR_P_TXT}",
    transform=ax_B.transAxes, ha="left", va="top", fontsize=5.6,
    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#999999", alpha=0.9), zorder=12)

# outside fold: light gray
if out_fold.sum() > 0:
    sub = lr_plot[out_fold]
    ax_B.scatter(
        jitter_log(sub.sr_vaf.values   * 100, sigma=0.02),
        floor_log( sub.hifi_vaf.values) * 100,
        s=4.5, color="#d9d9d9", alpha=0.45, linewidths=0, zorder=3,
        label=f">{FOLD_MAX}× discordant  (n={out_fold.sum():,})"
    )

# within fold: orange shades by SR VAF level
for lo_v, hi_v, col, lbl in reversed(ORANGE_BINS):
    mask = in_fold & (lr_plot.sr_vaf.values >= lo_v) & (lr_plot.sr_vaf.values < hi_v)
    if mask.sum() == 0:
        continue
    sub = lr_plot[mask]
    ax_B.scatter(
        jitter_log(sub.sr_vaf.values   * 100, sigma=0.02),
        floor_log( sub.hifi_vaf.values) * 100,
        s=5.5, color=col, alpha=0.70, linewidths=0, zorder=4,
        label=f"SR {lbl}  (n={mask.sum():,})"
    )

ax_B.set_xlabel("SR observed VAF (%)")
ax_B.set_ylabel("HiFi observed VAF (%)")
ax_B.set_title(
    f"B  SR vs HiFi — rescued variants  "
    f"(assembly−, SR+, HiFi+  ·  ≤{FOLD_MAX}× fold tolerance)",
    fontsize=7, loc="left"
)
ax_B.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), markerscale=2.5, ncol=1)

out = OUTDIR / "75_platform_comparison.pdf"
plt.savefig(out, dpi=300)
plt.close()
print(f"Figure → {out}")
