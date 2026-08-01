#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 66_array_end_methylation_ont_vs_hifi.py — End-anchored per-copy methylation, ONT versus HiFi, from the database methylation tables.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
66_array_end_methylation_ont_vs_hifi.py

End-anchored (telomere -> centromere) per-copy methylation, ONT versus HiFi,
built from the database copy_methylation tables (interior copies; border copies
are not in the methylation tables). Mirrors the layout of
array_end_methylation_hifi.pdf so ONT and HiFi can be compared directly.

Rank 1 = terminal interior copy at that end (the border copy is excluded).

Input : 5S_rDNA.db
Output: <FIVES_OUT>/10_border_methylation/array_end_methylation_ONT_vs_HiFi.pdf

Paths are read from environment variables (see repository README):
    FIVES_DB   path to 5S_rDNA.db
    FIVES_OUT  output directory
"""
import os
import sqlite3, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import warnings; warnings.simplefilter("ignore")
import matplotlib.pyplot as plt
from pathlib import Path

DB  = Path(os.environ.get("FIVES_DB", "5S_rDNA.db"))
OUT = Path(os.environ.get("FIVES_OUT", "output")) / "10_border_methylation"
NK, MIN_N = 20, 30
pctf = plt.FuncFormatter(lambda v, _: f"{v:.0f}%")

con = sqlite3.connect(DB)
def load(tbl):
    d = pd.read_sql(f"""
        SELECT a.sample_id s, h.hap_label hp, c.copy_number cn,
               cm.mean_meth*100 pct, cm.n_conf_calls nc
        FROM {tbl} cm JOIN copy c ON cm.copy_id=c.copy_id
        JOIN haplotype h ON c.haplotype_id=h.haplotype_id
        JOIN assembly a ON h.assembly_id=a.assembly_id
        WHERE c.border_note='interior' AND cm.n_conf_calls>=10
          AND a.cohort IN ('HPRC_Year1','HPRC_Release2')
          AND h.array_order_resolved != 'partial'""", con)   # order-dependent: drop unreliable order
    g = d.groupby(["s", "hp"])
    d = d.assign(mincn=g.cn.transform("min"), maxcn=g.cn.transform("max"),
                 M=g.cn.transform("count"))
    d = d[d.M >= 40].copy()
    d["rank_tel"] = d.cn - d.mincn + 1
    d["rank_cen"] = d.maxcn - d.cn + 1
    return d
data = {"ONT": load("copy_methylation"), "HiFi": load("copy_methylation_hifi")}
con.close()

ref = {k: v[(v.rank_tel > NK) & (v.rank_cen > NK)].pct.median() for k, v in data.items()}

def boxes(ax, d, rank_col, invert):
    cents, dat, meds, ns = [], [], [], []
    for r in range(1, NK + 1):
        v = d.loc[d[rank_col] == r, "pct"].dropna().values
        if len(v) < MIN_N:
            continue
        x = -r if invert else r
        cents.append(x); dat.append(v); meds.append(np.median(v)); ns.append(len(v))
    ax.boxplot(dat, positions=cents, widths=0.62, showfliers=False, patch_artist=True,
               boxprops=dict(facecolor="white", color="#333", lw=0.7),
               medianprops=dict(color="#c0392b", lw=1.3),
               whiskerprops=dict(color="#333", lw=0.6), capprops=dict(color="#333", lw=0.6))
    ax.plot(cents, meds, color="#c0392b", lw=1.4, alpha=0.7, zorder=4)
    for x, n in zip(cents, ns):
        ax.text(x, 1.5, str(n), ha="center", va="bottom", fontsize=4.6, color="#888", rotation=90)
    ax.set_ylim(0, 100); ax.yaxis.set_major_formatter(pctf)

fig, axes = plt.subplots(2, 2, figsize=(15, 9), sharey=True,
                         gridspec_kw=dict(wspace=0.05, hspace=0.32))
fig.suptitle("End-anchored per-copy CpG methylation (telomere → centromere), interior copies — "
             "ONT vs HiFi", fontsize=11, y=0.98)
for row, plat in enumerate(["ONT", "HiFi"]):
    d = data[plat]; rmed = ref[plat]
    axL, axR = axes[row, 0], axes[row, 1]
    boxes(axL, d, "rank_tel", invert=False)
    axL.set_xlim(0.3, NK + 0.7); axL.set_xticks(range(1, NK + 1, 2))
    axL.set_ylabel(f"{plat}\nCpG methylation (%)")
    axL.set_title(f"({'AC'[row]}) {plat} — telomere-proximal end")
    boxes(axR, d, "rank_cen", invert=True)
    axR.set_xlim(-(NK + 0.7), -0.3); axR.set_xticks([-r for r in range(1, NK + 1, 2)])
    axR.set_xticklabels([str(r) for r in range(1, NK + 1, 2)])
    axR.set_title(f"({'BD'[row]}) {plat} — centromere-proximal end")
    for ax in (axL, axR):
        ax.axhline(rmed, color="#4c72b0", ls="--", lw=1.0)
        ax.axhline(65, color="grey", ls=":", lw=0.8); ax.grid(axis="y", lw=0.3, alpha=0.4)
    if row == 1:
        axL.set_xlabel("Interior copy # from telomere (1qter) end")
        axR.set_xlabel("Interior copy # from centromere-proximal end")
    axR.text(0.98, rmed + 1.5, f"interior median {rmed:.0f}%",
             transform=axR.get_yaxis_transform(), ha="right", fontsize=7, color="#4c72b0")

plt.tight_layout()
out = OUT / "array_end_methylation_ONT_vs_HiFi.pdf"
fig.savefig(out, dpi=200, bbox_inches="tight"); plt.close(fig)
# data table
rows = []
for plat, d in data.items():
    for end, col in [("telomere", "rank_tel"), ("centromere", "rank_cen")]:
        for r in range(1, NK + 1):
            v = d.loc[d[col] == r, "pct"].dropna().values
            if len(v) < MIN_N: continue
            rows.append({"platform": plat, "array_end": end, "copy_from_end": r,
                         "n_haplotypes": len(v), "median": round(np.median(v), 1),
                         "mean": round(v.mean(), 1)})
pd.DataFrame(rows).to_csv(OUT / "data/array_end_methylation_ONT_vs_HiFi.tsv", sep="\t", index=False)
print(f"Saved: {out}")
for plat in ["ONT", "HiFi"]:
    print(f"  {plat}: {data[plat][['s','hp']].drop_duplicates().shape[0]} haplotypes, "
          f"interior median {ref[plat]:.0f}%")
