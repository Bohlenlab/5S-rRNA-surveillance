#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 47_clustering_single_panels.py — Produces two single-panel figures summarizing
# within-array variant clustering: the per-haplotype clustering index and the
# fraction of multi-copy variants whose carriers form a single contiguous tract.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
47_clustering_single_panels.py

Two single-panel figures summarizing within-array variant clustering:

  47a  Per-haplotype clustering index
       index = (observed adjacencies - expected) / (max - expected)
       0 = carriers as scattered as random; 1 = one contiguous tract.

  47b  Per haplotype, % of multi-copy variants whose carriers
       form a single contiguous tract, observed vs expected-if-independent
       (analytical null P(1 run) = (n-k+1)/C(n,k)).

Outputs: figures/47a_clustering_index_per_haplotype.pdf
         figures/47b_single_tract_fraction.pdf
"""

import os
import sqlite3
from math import comb
from collections import defaultdict
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

EXP = Path(os.environ.get("FIVES_DATA", "data")) / "exports"
EXP.mkdir(exist_ok=True)

DB = Path(os.environ.get("FIVES_DB", "5S_rDNA.db"))
FIG = Path(os.environ.get("FIVES_OUT", "output")) / "figures"


def load():
    con = sqlite3.connect(DB)
    cp = con.execute("""
        SELECT h.haplotype_id, c.copy_number
        FROM copy c JOIN haplotype h USING(haplotype_id) JOIN assembly a USING(assembly_id)
        WHERE c.border_note='interior' AND a.cohort IN ('HPRC_Year1','HPRC_Release2')
    """).fetchall()
    var = con.execute("""
        SELECT h.haplotype_id, c.copy_number, v.consensus_pos, v.alt
        FROM variant v JOIN copy c USING(copy_id)
        JOIN haplotype h USING(haplotype_id) JOIN assembly a USING(assembly_id)
        WHERE c.border_note='interior' AND a.cohort IN ('HPRC_Year1','HPRC_Release2')
          AND v.alignment_source='consensus_t2t' AND v.var_type='snp'
    """).fetchall()
    meta = {hid: (s, hl) for hid, s, hl in con.execute("""
        SELECT h.haplotype_id, a.sample_id, h.hap_label
        FROM haplotype h JOIN assembly a USING(assembly_id)
        WHERE a.cohort IN ('HPRC_Year1','HPRC_Release2')""").fetchall()}
    con.close()
    return cp, var, meta


def main():
    cp, var, meta = load()
    # rank map: copy_number -> 1..n per haplotype
    bycopy = defaultdict(list)
    for hid, cn in cp:
        bycopy[hid].append(cn)
    rank = {}; nmap = {}
    for hid, cns in bycopy.items():
        order = sorted(cns); nmap[hid] = len(order)
        rank[hid] = {c: i+1 for i, c in enumerate(order)}
    # variants -> carrier ranks per (hid,pos,alt)
    carr = defaultdict(lambda: defaultdict(set))
    for hid, cn, pos, alt in var:
        if cn in rank.get(hid, {}):
            carr[hid][(pos, alt)].add(rank[hid][cn])

    # original analysis: all multi-copy variants (2 <= k < n); singletons & fully
    # fixed (k=n) excluded because they admit no adjacency comparison.
    ci_rows = []
    obs_frac = []; null_frac = []
    pooled_obs = pooled_n = 0; pooled_null = 0.0
    for hid, vd in carr.items():
        n = nmap[hid]
        sumA = sumE = sumM = 0.0
        o = nn = 0; nu = 0.0
        for (pos, alt), ranks in vd.items():
            rs = sorted(ranks); k = len(rs)
            if k < 2 or k >= n:
                continue
            runs = 1 + sum(rs[i]-rs[i-1] > 1 for i in range(1, k))
            A_obs = k - runs
            A_exp = k*(k-1)/n
            A_max = k-1
            sumA += A_obs; sumE += A_exp; sumM += A_max
            nn += 1; o += (runs == 1); nu += (n-k+1)/comb(n, k)
        if sumM > sumE and nn > 0:
            s, hl = meta.get(hid, ("", ""))
            ci_rows.append({"sample_id": s, "hap_label": hl, "haplotype_id": hid,
                            "n_copies": n, "n_multicopy_variant_sites": nn,
                            "adjacencies_observed": round(sumA, 4),
                            "adjacencies_expected": round(sumE, 4),
                            "adjacencies_max": int(sumM),
                            "clustering_index": round((sumA - sumE)/(sumM - sumE), 5)})
            obs_frac.append(o/nn*100); null_frac.append(nu/nn*100)
            pooled_obs += o; pooled_n += nn; pooled_null += nu
    ci_df = pd.DataFrame(ci_rows).sort_values("clustering_index", ascending=False)
    ci_df.to_csv(EXP/"clustering_index_per_haplotype.csv", index=False)
    ci_hap = ci_df["clustering_index"].values
    obs_frac = np.array(obs_frac); null_frac = np.array(null_frac)
    w_p = stats.wilcoxon(ci_hap, alternative="greater")[1]
    po, pn = pooled_obs/pooled_n*100, pooled_null/pooled_n*100
    print(f"Clustering index: median {np.median(ci_hap):.3f}, n={len(ci_hap)} haps, "
          f"{(ci_hap>0).mean()*100:.0f}% >0, Wilcoxon p={w_p:.1e}")
    print(f"Exported -> {EXP/'clustering_index_per_haplotype.csv'}")
    print(f"Single-tract fraction: observed {po:.0f}% vs null {pn:.1f}% ({po/pn:.0f}x), n={pooled_n:,} variants")

    # ── 47a: clustering index ─────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(7.5, 6))
    ax.hist(ci_hap, bins=30, color="#2c7fb8", edgecolor="white", linewidth=0.4)
    ax.axvline(0, color="black", ls="--", lw=1.5, label="independent mutations (0)")
    ax.axvline(np.median(ci_hap), color="#c0392b", lw=2.5,
               label=f"median = {np.median(ci_hap):.2f}")
    ax.set_xlabel("Variant clustering index per haplotype\n"
                  "(0 = mutations independent   ·   1 = single contiguous tract)", fontsize=11)
    ax.set_ylabel("Number of haplotypes", fontsize=11)
    ax.set_title("Per-haplotype within-array variant clustering index",
                 fontsize=12, fontweight="bold")
    ax.text(0.97, 0.78, f"{len(ci_hap)} haplotypes\n{(ci_hap>0).mean()*100:.0f}% clustered (>0)\n"
            f"Wilcoxon p = {w_p:.0e}",
            transform=ax.transAxes, ha="right", va="top", fontsize=10,
            bbox=dict(boxstyle="round", fc="white", ec="gray", alpha=0.9))
    ax.legend(fontsize=9.5, loc="upper left")
    fig.tight_layout()
    fig.savefig(FIG/"47a_clustering_index_per_haplotype.pdf", bbox_inches="tight", dpi=160)
    plt.close(fig)

    # ── 47b: single-tract fraction observed vs null ───────────────────────────
    fig, ax = plt.subplots(figsize=(7.5, 6))
    bins = np.linspace(0, max(obs_frac.max(), 1), 30)
    ax.hist(null_frac, bins=bins, color="#bdbdbd", edgecolor="white", linewidth=0.3,
            alpha=0.85, label=f"if mutations independent (mean {pn:.0f}%)")
    ax.hist(obs_frac, bins=bins, color="#2166ac", edgecolor="white", linewidth=0.3,
            alpha=0.7, label=f"observed (mean {po:.0f}%)")
    ax.axvline(po, color="#2166ac", lw=2); ax.axvline(pn, color="#525252", lw=2, ls="--")
    ax.set_xlabel("% of a haplotype's multi-copy variants\nwhose copies form ONE contiguous tract", fontsize=11)
    ax.set_ylabel("Number of haplotypes", fontsize=11)
    ax.set_title("Multi-copy variants whose carriers form one contiguous tract\n"
                 f"observed {po:.0f}% vs null {pn:.1f}%", fontsize=12, fontweight="bold")
    ax.text(0.97, 0.6, f"{len(obs_frac)} haplotypes\npooled: {po:.0f}% obs\nvs {pn:.1f}% null\n"
            f"= {po/pn:.0f}× enrichment", transform=ax.transAxes, ha="right", va="top",
            fontsize=10, bbox=dict(boxstyle="round", fc="white", ec="gray", alpha=0.9))
    ax.legend(fontsize=9.5, loc="upper right")
    fig.tight_layout()
    fig.savefig(FIG/"47b_single_tract_fraction.pdf", bbox_inches="tight", dpi=160)
    plt.close(fig)
    print("Saved 47a + 47b")


if __name__ == "__main__":
    main()
