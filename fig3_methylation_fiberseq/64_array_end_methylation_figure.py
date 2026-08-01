#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 64_array_end_methylation_figure.py — End-anchored per-copy methylation across the first and last copies of each array (HiFi).
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
64_array_end_methylation_figure.py

End-anchored (unscaled) per-copy CpG methylation, oriented telomere ->
centromere: the first 20 copies from the telomere-proximal end and the last 20
from the centromere-proximal end, one copy per bin.

Orientation: every HiFi-methylation array is plus-strand / reference-oriented
(335/335 plus); the array lies at chr1q42 (~228 Mb), ~20 Mb from the q-telomere.
The high-coordinate array end (rank5) is labelled telomere-proximal and the
low-coordinate end (rank3) centromere-proximal.

Input:  array_end_methylation_hifi.tsv
Output: <FIVES_OUT>/10_border_methylation/array_end_methylation_hifi.pdf

Paths are read from environment variables (see repository README):
    FIVES_DATA  directory with input TSVs
    FIVES_OUT   output directory
"""
import os
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

HPRC = Path(os.environ.get("FIVES_DATA", "data"))
OUT  = Path(os.environ.get("FIVES_OUT", "output")) / "10_border_methylation"
NK   = 20            # copies to show from each end
MIN_CALLS = 10
pctf = plt.FuncFormatter(lambda v, _: f"{v:.0f}%")

d = pd.read_csv(HPRC / "array_end_methylation_hifi.tsv", sep="\t")
d["pct"] = d["mean_meth"] * 100
d = d[(d["n_conf"] >= MIN_CALLS) & (d["M"] >= 40)]   # clean 5'/3' separation
interior_ref = d[(d["rank5"] > NK) & (d["rank3"] > NK)]["pct"].median() \
               if ((d["rank5"] > NK) & (d["rank3"] > NK)).any() else d["pct"].median()

MIN_N = 30   # drop ranks with too few haplotypes (unreliable terminal/border copies)

def boxes(ax, rank_col, ranks, invert=False):
    data, cents, meds, ns = [], [], [], []
    for r in ranks:
        v = d.loc[d[rank_col] == r, "pct"].dropna().values
        if len(v) < MIN_N:          # unreliable — skip (e.g. sparse short 3' border)
            continue
        x = (-r if invert else r)
        data.append(v); cents.append(x); meds.append(np.median(v)); ns.append(len(v))
    bp = ax.boxplot(data, positions=cents, widths=0.62, showfliers=False,
                    patch_artist=True,
                    boxprops=dict(facecolor="white", color="#333", lw=0.7),
                    medianprops=dict(color="#c0392b", lw=1.3),
                    whiskerprops=dict(color="#333", lw=0.6),
                    capprops=dict(color="#333", lw=0.6))
    ax.plot(cents, meds, color="#c0392b", lw=1.4, alpha=0.7, zorder=4)
    # annotate n under each box so coverage is transparent
    for x, n in zip(cents, ns):
        ax.text(x, 1.5, str(n), ha="center", va="bottom", fontsize=5.2, color="#777",
                rotation=90)
    ax.axhline(interior_ref, color="#4c72b0", ls="--", lw=1.0)
    ax.axhline(65, color="grey", ls=":", lw=0.8)
    ax.set_ylim(0, 100); ax.yaxis.set_major_formatter(pctf)
    return cents, ns

fig, (axL, axR) = plt.subplots(1, 2, figsize=(15, 5), sharey=True,
                               gridspec_kw=dict(wspace=0.06))
fig.suptitle("Per-copy CpG methylation oriented telomere → centromere (HiFi)",
             fontsize=11, y=0.99)

# Telomere-proximal end = high-coordinate end = rank5 (telomere terminal at left)
xs5, ns5 = boxes(axL, "rank5", range(1, NK + 1))
axL.set_xlim(0.3, NK + 0.7)
axL.set_xticks(range(1, NK + 1, 2))
axL.set_xlabel("Copy number from telomere (1qter) end")
axL.set_ylabel("Per-copy CpG methylation (%)")
axL.set_title("(A) Telomere-proximal end")
axL.text(1, 3, "terminal", ha="center", fontsize=6.5, color="#a33", rotation=90, va="bottom")
axL.grid(axis="y", lw=0.3, alpha=0.4)

# Centromere-proximal end = low-coordinate end = rank3 (centromere terminal at far right)
xs3, ns3 = boxes(axR, "rank3", range(1, NK + 1), invert=True)
axR.set_xlim(-(NK + 0.7), -0.3)
axR.set_xticks([-r for r in range(1, NK + 1, 2)])
axR.set_xticklabels([str(r) for r in range(1, NK + 1, 2)])
axR.set_xlabel("Copy number from centromere-proximal end")
axR.set_title("(B) Centromere-proximal end")
axR.grid(axis="y", lw=0.3, alpha=0.4)
axR.text(0.98, interior_ref + 1.5, f"interior median {interior_ref:.0f}%",
         transform=axR.get_yaxis_transform(), ha="right", fontsize=7, color="#4c72b0")

plt.tight_layout()
(OUT / "data").mkdir(parents=True, exist_ok=True)
out = OUT / "array_end_methylation_hifi.pdf"
fig.savefig(out, dpi=200, bbox_inches="tight")
plt.close(fig)

# data table: per-rank summary both ends
rows = []
for end, col in [("telomere", "rank5"), ("centromere", "rank3")]:
    for r in range(1, NK + 1):
        v = d.loc[d[col] == r, "pct"].dropna().values
        if len(v) < 5:
            continue
        q = np.percentile(v, [10, 25, 50, 75, 90])
        rows.append({"array_end": end, "copy_from_end": r, "n_haplotypes": len(v),
                     "median": round(q[2], 1), "mean": round(v.mean(), 1),
                     "q10": round(q[0], 1), "q25": round(q[1], 1),
                     "q75": round(q[3], 1), "q90": round(q[4], 1),
                     "pct_below65": round((v < 65).mean() * 100, 1)})
pd.DataFrame(rows).to_csv(OUT / "data/array_end_methylation_hifi.tsv", sep="\t", index=False)
print(f"Saved: {out}")
print(f"interior reference median: {interior_ref:.1f}%")
