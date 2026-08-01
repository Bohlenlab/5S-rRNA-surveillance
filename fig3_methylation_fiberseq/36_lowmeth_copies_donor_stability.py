#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 36_lowmeth_copies_donor_stability.py — Per-donor stability of the number of lowly-methylated 5S copies relative to total copy number.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
36_lowmeth_copies_donor_stability.py

Per-donor stability of the number of lowly-methylated 5S copies (summed across
both haplotypes) relative to total copy number.

Per copy: mean methylation from confident calls (mod_qual<=0.2 or >=0.8),
minimum MIN_CALLS_PER_COPY calls; a copy is low if mean <65%. Per donor (both
haplotypes required): n_low = lowly-meth copies, n_total = classified copies,
frac_low = n_low / n_total.

Analyses:
  - dispersion (CV) of n_low across donors versus that of n_total and n_high
  - n_low versus total copy number, comparing a fixed-number model (horizontal
    set-point) with a fixed-fraction model (through-origin, proportional) by
    residual sum of squares.

Output: <FIVES_OUT>/03_methylation/36_lowmeth_copies_donor_stability.pdf

Paths are read from environment variables (see repository README):
    FIVES_DATA  directory with methylation/ and databases/ inputs
    FIVES_OUT   output directory
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats
from pathlib import Path

METH      = Path(os.environ.get("FIVES_DATA", "data")) / "methylation"
DATABASES = Path(os.environ.get("FIVES_DATA", "data")) / "databases"
OUTDIR    = Path(os.environ.get("FIVES_OUT", "output")) / "03_methylation"

MOD_HI, MOD_LO = 0.8, 0.2
MIN_CALLS_PER_COPY = 10
HI_CUT = 0.65


def per_hap_counts():
    rows = []
    for meth_dir in sorted(METH.iterdir()):
        if not meth_dir.is_dir():
            continue
        sample = meth_dir.name
        for hap in ("hap1", "hap2"):
            ann = meth_dir / f"{sample}_{hap}_modkit_annotated.tsv"
            db  = DATABASES / f"{sample}_{hap}.tsv"
            if not ann.exists() or ann.stat().st_size < 500 or not db.exists():
                continue
            d = pd.read_csv(ann, sep="\t", usecols=["mod_qual", "copy_id"])
            d = d[(d["mod_qual"] <= MOD_LO) | (d["mod_qual"] >= MOD_HI)].copy()
            if d.empty:
                continue
            d["is_meth"] = (d["mod_qual"] >= MOD_HI).astype(np.int8)
            cs = d.groupby("copy_id")["is_meth"].agg(n="count", m="mean")
            cs = cs[cs["n"] >= MIN_CALLS_PER_COPY]
            if len(cs) == 0:
                continue
            nt = len(cs); nl = int((cs["m"] < HI_CUT).sum())
            rows.append({"sample": sample, "hap": hap,
                         "n_total": nt, "n_low": nl, "n_high": nt - nl})
    return pd.DataFrame(rows)


def cv(x):
    return np.std(x, ddof=1) / np.mean(x)


def main():
    print("Counting ...", flush=True)
    h = per_hap_counts()
    g = h.groupby("sample").filter(lambda d: set(d["hap"]) >= {"hap1", "hap2"})
    donor = g.groupby("sample").agg(n_low=("n_low", "sum"),
                                    n_high=("n_high", "sum"),
                                    n_total=("n_total", "sum")).reset_index()
    donor["frac_low"] = donor["n_low"] / donor["n_total"]
    N = len(donor)
    print(f"{N} donors with both haplotypes")

    nlow, ntot, nhigh, frac = (donor["n_low"].values, donor["n_total"].values,
                               donor["n_high"].values, donor["frac_low"].values)

    # dispersion
    stats_tbl = {
        "n_low":   (nlow.mean(), nlow.std(ddof=1), cv(nlow)),
        "n_high":  (nhigh.mean(), nhigh.std(ddof=1), cv(nhigh)),
        "n_total": (ntot.mean(), ntot.std(ddof=1), cv(ntot)),
        "frac_low":(frac.mean(), frac.std(ddof=1), cv(frac)),
    }
    print("\nDispersion across donors (mean / SD / CV):")
    for k, (m, s, c) in stats_tbl.items():
        print(f"  {k:9s}: mean={m:.3f}  SD={s:.3f}  CV={c:.2f}")

    # n_low vs total copy number
    r, p = stats.pearsonr(ntot, nlow)
    slope, intercept, rr, pp, se = stats.linregress(ntot, nlow)
    print(f"\nn_low ~ n_total: Pearson r={r:.2f} (p={p:.1e}); "
          f"slope={slope:.3f} copies per total copy (SE {se:.3f})")

    # competing models on n_low:
    #   fixed-NUMBER  : predict mean(n_low)            (set-point)
    #   fixed-FRACTION: predict mean(frac_low)*n_total (proportional)
    pred_num = np.full(N, nlow.mean())
    pred_frac = frac.mean() * ntot
    rss_num = np.sum((nlow - pred_num) ** 2)
    rss_frac = np.sum((nlow - pred_frac) ** 2)
    print(f"\nResidual SS — fixed-number model: {rss_num:.0f}   "
          f"fixed-fraction model: {rss_frac:.0f}  "
          f"(lower = better fit)")
    better = "FIXED NUMBER (set-point)" if rss_num < rss_frac else "FIXED FRACTION (proportional)"
    print(f"Better-fitting model: {better}")

    # ── figure ─────────────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(15, 11))
    gs = gridspec.GridSpec(2, 2, hspace=0.30, wspace=0.27)

    # A: distribution of per-donor n_low
    ax = fig.add_subplot(gs[0, 0])
    ax.hist(nlow, bins=np.arange(nlow.min()-0.5, nlow.max()+1.5),
            color="#2c7fb8", edgecolor="white", linewidth=0.4)
    ax.axvline(nlow.mean(), color="#c0392b", lw=2,
               label=f"mean={nlow.mean():.1f}")
    ax.axvspan(nlow.mean()-nlow.std(ddof=1), nlow.mean()+nlow.std(ddof=1),
               color="#c0392b", alpha=0.12, label=f"±1 SD ({nlow.std(ddof=1):.1f})")
    ax.set_xlabel("Lowly-methylated copies per donor (hap1+hap2, <65%)")
    ax.set_ylabel("donors")
    ax.set_title(f"A   Per-donor count of lowly-meth copies  (CV={cv(nlow):.2f})",
                 fontweight="bold", loc="left", fontsize=10)
    ax.legend(fontsize=8)

    # B: n_low vs total copy number — set-point vs proportional
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.scatter(ntot, nlow, s=42, color="#2c7fb8", alpha=0.8,
                edgecolors="white", linewidths=0.4, zorder=3)
    xs = np.linspace(ntot.min(), ntot.max(), 50)
    ax2.axhline(nlow.mean(), color="#c0392b", lw=2,
                label=f"fixed number (set-point), RSS={rss_num:.0f}")
    ax2.plot(xs, frac.mean()*xs, color="#777", lw=2, ls="--",
             label=f"fixed fraction ({frac.mean()*100:.0f}%), RSS={rss_frac:.0f}")
    ax2.set_xlabel("Total classified copies per donor (hap1+hap2)")
    ax2.set_ylabel("Lowly-methylated copies per donor (<65%)")
    ax2.set_title(f"B   n_low vs total copy number  r={r:.2f}, slope={slope:.3f}",
                  fontweight="bold", loc="left", fontsize=10)
    ax2.legend(fontsize=8, loc="upper left")
    ax2.set_ylim(bottom=0)

    # C: CV comparison
    ax3 = fig.add_subplot(gs[1, 0])
    keys = ["n_low", "n_high", "n_total", "frac_low"]
    cvs = [stats_tbl[k][2] for k in keys]
    cols = ["#2c7fb8", "#41ab5d", "#888888", "#f0a202"]
    ax3.bar(range(len(keys)), cvs, color=cols, edgecolor="black", linewidth=0.5)
    for i, c in enumerate(cvs):
        ax3.text(i, c+0.005, f"{c:.2f}", ha="center", va="bottom", fontsize=9)
    ax3.set_xticks(range(len(keys)))
    ax3.set_xticklabels(["# low\n(<65%)", "# high\n(≥65%)", "# total\ncopies",
                         "fraction\nlow"], fontsize=8)
    ax3.set_ylabel("Coefficient of variation (SD/mean)")
    ax3.set_title("C   Coefficient of variation by quantity",
                  fontweight="bold", loc="left", fontsize=10)

    # D: per-donor stacked low/high, sorted by total
    ax4 = fig.add_subplot(gs[1, 1])
    order = np.argsort(ntot)
    xx = np.arange(N)
    ax4.bar(xx, nlow[order], color="#2c7fb8", label="lowly-meth (<65%)")
    ax4.bar(xx, nhigh[order], bottom=nlow[order], color="#dddddd",
            label="highly-meth (≥65%)")
    ax4.plot(xx, nlow[order], color="#c0392b", lw=1.2, marker="o", ms=2.5,
             label="low count")
    ax4.set_xlabel("donor (sorted by total copy number →)")
    ax4.set_ylabel("copies per donor")
    ax4.set_title("D   Per-donor low/high copies, sorted by total copy number",
                  fontweight="bold", loc="left", fontsize=10)
    ax4.legend(fontsize=8, loc="upper left")
    ax4.set_xticks([])

    fig.suptitle(
        f"Per-donor count of lowly-methylated 5S copies vs total copy number  "
        f"({N} donors, HPRC Year 1 ONT)",
        fontsize=12, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = OUTDIR / "36_lowmeth_copies_donor_stability.pdf"
    fig.savefig(out, bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"\nSaved {out}")


if __name__ == "__main__":
    main()
