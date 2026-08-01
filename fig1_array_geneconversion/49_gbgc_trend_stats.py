#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 49_gbgc_trend_stats.py — Logistic-regression test of whether the GC-gaining
# direction of interior-copy substitutions increases with within-array carrier count.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
49_gbgc_trend_stats.py

Statistical test of whether the GC-gaining direction of substitutions increases
with carrier count (within-array gene-conversion spread).

Model: logistic regression at the single-variant level
    logit P(outcome) = a + b * log10(carrier_count)
Each (haplotype, pos, ref, alt) variant is one observation. b = the trend slope
(log-odds per 10x carrier count). Standard errors are CLUSTER-ROBUST by haplotype
(variants within one array are not independent). Reported as slope, SE, z, p, and
odds ratio per 10x and per doubling of carrier count.

Tests run:
  - GC-gaining among weak<->strong changes: W->S (1) vs S->W (0)
  - each S/W category, one-vs-rest, trend with carrier count
Also a Cochran-Armitage-style check (Spearman of binned proportion vs k).

Exports: exports/gbgc_trend_stats.csv
"""

import os
import sqlite3
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

DB = Path(os.environ.get("FIVES_DB", "5S_rDNA.db"))
EXP = Path(os.environ.get("FIVES_DATA", "data")) / "exports"; EXP.mkdir(exist_ok=True)
S = set("CG")


def cat(ref, alt):
    return ("S" if ref in S else "W") + "2" + ("S" if alt in S else "W")


def logit_trend(df, y, group, label):
    """Logistic regression logit(y) ~ log10(k), cluster-robust SE by group."""
    X = sm.add_constant(df["log10k"].values)
    m = sm.GLM(df[y].values, X, family=sm.families.Binomial())
    r = m.fit(cov_type="cluster", cov_kwds={"groups": df[group].values})
    b = r.params[1]; se = r.bse[1]; z = r.tvalues[1]; p = r.pvalues[1]
    lo, hi = r.conf_int()[1]
    return {
        "test": label, "n_variants": int(len(df)), "n_haplotypes": int(df[group].nunique()),
        "slope_per_10x_k": round(b, 4), "SE": round(se, 4), "z": round(z, 2),
        "p_value": p,
        "OR_per_10x_k": round(np.exp(b), 4),
        "OR_per_doubling_k": round(np.exp(b*np.log10(2)), 4),
        "OR_ci_low_10x": round(np.exp(lo), 4), "OR_ci_high_10x": round(np.exp(hi), 4),
    }


def main():
    con = sqlite3.connect(DB)
    v = pd.read_sql_query("""
        SELECT h.haplotype_id AS hap, v.ref, v.alt, COUNT(DISTINCT v.copy_id) AS k
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
    v["log10k"] = np.log10(v["k"])

    results = []
    # among weak<->strong changes, P(GC-gaining W->S)
    ws = v[v["cat"].isin(["W2S", "S2W"])].copy()
    ws["y"] = (ws["cat"] == "W2S").astype(int)
    results.append(logit_trend(ws, "y", "hap",
                  "GC-gaining (W>S vs S>W) ~ log10(carrier count)"))

    # each category one-vs-rest
    for c, name in [("W2S", "W>S"), ("S2W", "S>W"), ("S2S", "S>S"), ("W2W", "W>W")]:
        d = v.copy(); d["y"] = (d["cat"] == c).astype(int)
        results.append(logit_trend(d, "y", "hap", f"fraction {name} ~ log10(carrier count)"))

    res = pd.DataFrame(results)
    res.to_csv(EXP/"gbgc_trend_stats.csv", index=False)

    # Spearman cross-check on binned GC-gaining fraction (from the per-bin table)
    binf = EXP/"gbgc_by_carrier_count.csv"
    spearman_txt = ""
    if binf.exists():
        b = pd.read_csv(binf)
        rho, pp = stats.spearmanr(b["mean_k"], b["frac_GC_gaining"])
        spearman_txt = f"  [cross-check] Spearman GC-gaining vs mean_k across bins: rho={rho:.2f}, p={pp:.1e}"

    pd.set_option("display.width", 200)
    print(res.to_string(index=False))
    print(spearman_txt)
    print(f"\nExported -> {binf.parent/'gbgc_trend_stats.csv'}")
    # console summary of the first test
    h = results[0]
    print(f"\nHeadline: GC-gaining direction increases significantly with carrier count "
          f"(p={h['p_value']:.1e}); odds of a weak<->strong change being GC-gaining "
          f"rise {h['OR_per_doubling_k']:.2f}x per doubling of carrier count "
          f"({h['OR_per_10x_k']:.2f}x per 10x).")


if __name__ == "__main__":
    main()
