#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 48_gbgc_resolved.py — Quantifies GC-biased gene conversion (gBGC) as a gradient
# of the GC-gaining variant fraction across within-array carrier count.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
48_gbgc_resolved.py

GC-biased gene conversion (gBGC) quantified as a gradient across carrier count k
(= number of distinct carrier copies of a variant within its array), instead of a
singletons-vs-tracts split.

For each (haplotype, pos, ref, alt) variant we record k (distinct carrier copies)
and its strong<->weak class (S=G/C, W=A/T). Per carrier-count bin the GC-gaining
fraction (W->S among weak<->strong changes) and the composition of the four S/W
categories are computed with Wilson confidence intervals.

Outputs:
  figures/48_gbgc_resolved.pdf
  exports/gbgc_by_carrier_count.csv
"""

import os
import sqlite3
from math import sqrt
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DB = Path(os.environ.get("FIVES_DB", "5S_rDNA.db"))
FIG = Path(os.environ.get("FIVES_OUT", "output")) / "figures"
EXP = Path(os.environ.get("FIVES_DATA", "data")) / "exports"; EXP.mkdir(exist_ok=True)
S = set("CG")
CATS = ["W→S", "S→S", "W→W", "S→W"]
CAT_COLOR = {"W→S": "#2ca02c", "S→S": "#1f77b4", "W→W": "#999999", "S→W": "#d62728"}


def cat(ref, alt):
    return ("S" if ref in S else "W") + "→" + ("S" if alt in S else "W")


def wilson(x, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    p = x/n; d = 1 + z*z/n
    c = (p + z*z/(2*n))/d
    h = z*sqrt(p*(1-p)/n + z*z/(4*n*n))/d
    return c-h, c+h


def main():
    con = sqlite3.connect(DB)
    v = pd.read_sql_query("""
        SELECT v.ref, v.alt, COUNT(DISTINCT v.copy_id) AS k
        FROM variant v JOIN copy c USING(copy_id)
        JOIN haplotype h USING(haplotype_id) JOIN assembly a USING(assembly_id)
        WHERE c.border_note='interior' AND a.cohort IN ('HPRC_Year1','HPRC_Release2')
          AND v.alignment_source='consensus_t2t' AND v.var_type='snp'
          AND length(v.ref)=1 AND length(v.alt)=1
          AND v.ref IN ('A','C','G','T') AND v.alt IN ('A','C','G','T')
          AND v.masked=0
        GROUP BY h.haplotype_id, v.consensus_pos, v.ref, v.alt
    """, con)
    con.close()
    v["cat"] = [cat(r, a) for r, a in zip(v.ref, v.alt)]

    # carrier-count bins — finer at the high-k end
    edges = [1,2,3,4,5,6,7,8,11,15,21,31,46,66,96,141,10**9]
    labels = ["1","2","3","4","5","6","7","8-10","11-14","15-20","21-30",
              "31-45","46-65","66-95","96-140","141+"]
    v["kbin"] = pd.cut(v["k"], bins=edges, right=False, labels=labels)

    rows = []
    for kb in labels:
        d = v[v["kbin"] == kb]
        if len(d) == 0:
            continue
        cc = {c: int((d["cat"] == c).sum()) for c in CATS}
        ws, sw = cc["W→S"], cc["S→W"]
        tot = len(d)
        frac_gc = ws/(ws+sw) if (ws+sw) else np.nan      # GC-gaining among W<->S
        lo, hi = wilson(ws, ws+sw)
        row = {"k_bin": kb, "mean_k": round(d["k"].mean(), 2), "n_variants": tot,
               "n_W2S": ws, "n_S2W": sw, "n_S2S": cc["S→S"], "n_W2W": cc["W→W"],
               "ratio_W2S_over_S2W": round(ws/sw, 3) if sw else np.nan,
               "frac_GC_gaining": round(frac_gc, 4),
               "gc_ci_low": round(lo, 4), "gc_ci_high": round(hi, 4)}
        # per-category proportion + Wilson CI (for Panel B error bars)
        for c in CATS:
            tag = c.replace('→', '2')
            p = cc[c]/tot
            clo, chi = wilson(cc[c], tot)
            row[f"frac_{tag}"] = round(p, 4)
            row[f"frac_{tag}_lo"] = round(clo, 4)
            row[f"frac_{tag}_hi"] = round(chi, 4)
        rows.append(row)
    tab = pd.DataFrame(rows)
    tab.to_csv(EXP/"gbgc_by_carrier_count.csv", index=False)
    print(tab[["k_bin","mean_k","n_variants","n_W2S","n_S2W","ratio_W2S_over_S2W",
               "frac_GC_gaining"]].to_string(index=False))
    print(f"\nExported -> {EXP/'gbgc_by_carrier_count.csv'}")

    # ── figure ─────────────────────────────────────────────────────────────────
    x = tab["mean_k"].values
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))

    # A: GC-gaining fraction vs carrier count (gBGC gradient)
    ax = axes[0]
    y = tab["frac_GC_gaining"].values*100
    yerr = np.vstack([y - tab["gc_ci_low"].values*100, tab["gc_ci_high"].values*100 - y])
    ax.errorbar(x, y, yerr=yerr, fmt="o-", color="#2ca02c", lw=2, ms=6, capsize=3,
                ecolor="#2ca02c", markeredgecolor="black")
    base = y[0]
    ax.axhline(base, color="gray", ls="--", lw=1, label=f"new mutations (k=1): {base:.0f}%")
    ax.set_xscale("log")
    ax.set_xlabel("Carrier copies per variant  (k)")
    ax.set_ylabel("% of weak↔strong changes that are GC-gaining (W→S)")
    ax.set_title("A   GC-gaining fraction (W→S) among weak↔strong changes\nby carrier count",
                 fontweight="bold", loc="left", fontsize=11)
    ax.legend(fontsize=9, loc="upper left")
    ax.grid(alpha=0.3, which="both")

    # B: composition of all four S/W categories vs carrier count (Wilson 95% CI)
    ax2 = axes[1]
    for c in CATS:
        tag = c.replace('→', '2')
        p = tab[f"frac_{tag}"].values*100
        ylo = p - tab[f"frac_{tag}_lo"].values*100
        yhi = tab[f"frac_{tag}_hi"].values*100 - p
        ax2.errorbar(x, p, yerr=np.vstack([ylo, yhi]), fmt="o-", color=CAT_COLOR[c],
                     lw=2, ms=4.5, capsize=2, elinewidth=0.9, label=c)
    ax2.set_xscale("log")
    ax2.set_xlabel("Carrier copies per variant  (k)")
    ax2.set_ylabel("% of variants")
    ax2.set_title("B   Composition of the four S/W categories by carrier count",
                  fontweight="bold", loc="left", fontsize=11)
    ax2.legend(fontsize=9, title="G/C change", ncol=2)
    ax2.grid(alpha=0.3, which="both")

    fig.suptitle("GC-biased gene conversion in the 5S array — resolved by carrier count",
                 fontsize=12.5, fontweight="bold")
    fig.tight_layout(rect=[0,0,1,0.96])
    out = FIG/"48_gbgc_resolved.pdf"
    fig.savefig(out, bbox_inches="tight", dpi=160); plt.close(fig)
    print(f"Saved {out.name}")


if __name__ == "__main__":
    main()
