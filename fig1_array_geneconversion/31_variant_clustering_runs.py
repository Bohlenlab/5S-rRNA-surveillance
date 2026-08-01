#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 31_variant_clustering_runs.py
# Within-array variant clustering by an exact runs test.
#
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the terms in the LICENSE file at the repository root.
# -----------------------------------------------------------------------------
"""
Within-array variant clustering by an exact runs test.

For a variant carried by k of an array's n linearly ordered interior copies, the
carriers are modelled under the null as a uniformly random k-subset of the n
positions. The clustering statistic is the number of maximal contiguous carrier
runs R (fewer runs = more clustered); adjacencies A = k - R.

Exact null (Wald-Wolfowitz runs of one type):
    P(R = r) = C(k-1, r-1) * C(n-k+1, r) / C(n, k)
    E[R]     = k (n - k + 1) / n
yielding, per site, an exact one-sided p-value P(R <= R_obs), valid for any k.

Per-site values are aggregated to:
  (A) the exact per-site p-value distribution and the fraction below 0.05,
  (B) adjacency fold-enrichment (obs/exp) versus carrier count k,
  (C) a per-haplotype clustering index tested across haplotypes (Wilcoxon),
  (D) the clustering index by repeat-unit region and carrier-count class.

Input : 5S_rDNA.db (tables copy, variant, haplotype, assembly).
Output: <FIVES_OUT>/05_gene_conversion/31_variant_clustering_runs.pdf

Paths are read from environment variables (see repository README):
    FIVES_DB   path to 5S_rDNA.db
    FIVES_OUT  output directory
"""

import os
import sqlite3
from math import comb
from functools import lru_cache
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
from scipy import stats

DB_PATH = Path(os.environ.get("FIVES_DB", "5S_rDNA.db"))
FIG_DIR = Path(os.environ.get("FIVES_OUT", "output")) / "05_gene_conversion"

GENE_START, GENE_END = 630, 748
REGION_ORDER = ["nts_pre", "gene", "nts_post"]
REGION_LABELS = {"nts_pre": "NTS-pre", "gene": "5S gene", "nts_post": "NTS-post"}
REGION_COLORS = {"nts_pre": "#4575b4", "gene": "#ff9900", "nts_post": "#d73027"}


# ── exact runs-distribution statistics, memoized per (k, n) ───────────────────

@lru_cache(maxsize=None)
def runs_stats(k, n):
    """Return (cum_le, E, Var) for the number of success-runs R of k carriers
    among n linearly ordered positions. cum_le[r] = P(R <= r)."""
    m = n - k                      # non-carriers
    rmax = min(k, m + 1)
    denom = comb(n, k)
    pmf = np.zeros(rmax + 1)        # index by r (1..rmax)
    for r in range(1, rmax + 1):
        pmf[r] = comb(k - 1, r - 1) * comb(m + 1, r) / denom
    s = pmf.sum()
    if s > 0:
        pmf /= s
    rs = np.arange(rmax + 1)
    E = float((rs * pmf).sum())
    Var = float((rs**2 * pmf).sum() - E**2)
    cum = np.cumsum(pmf)
    return cum, E, Var


def n_runs(sorted_ranks):
    """Number of maximal contiguous runs among 1-based integer ranks."""
    R = 1
    for i in range(1, len(sorted_ranks)):
        if sorted_ranks[i] - sorted_ranks[i-1] > 1:
            R += 1
    return R


# ── load ──────────────────────────────────────────────────────────────────────

def load():
    con = sqlite3.connect(DB_PATH)
    # interior copies per haplotype -> ordered rank map
    cp = pd.read_sql_query("""
        SELECT h.haplotype_id, c.copy_number
        FROM copy c JOIN haplotype h USING(haplotype_id)
                    JOIN assembly a USING(assembly_id)
        WHERE c.border_note='interior'
          AND a.cohort IN ('HPRC_Year1','HPRC_Release2')
    """, con)
    var = pd.read_sql_query("""
        SELECT h.haplotype_id, c.copy_number, v.consensus_pos, v.alt, v.region
        FROM variant v
        JOIN copy c USING(copy_id)
        JOIN haplotype h USING(haplotype_id)
        JOIN assembly a USING(assembly_id)
        WHERE c.border_note='interior'
          AND a.cohort IN ('HPRC_Year1','HPRC_Release2')
          AND v.alignment_source='consensus_t2t' AND v.var_type='snp'
          AND length(v.ref)=1 AND length(v.alt)=1
          AND v.ref IN ('A','C','G','T') AND v.alt IN ('A','C','G','T')
    """, con)
    con.close()
    return cp, var


def build_sites(cp, var):
    """One record per (haplotype, pos, alt) with k>=2 carriers among interior copies."""
    # rank map: copy_number -> 1-based rank within each haplotype's interior copies
    rankmap = {}
    nmap = {}
    for hid, g in cp.groupby("haplotype_id"):
        order = sorted(g["copy_number"].tolist())
        rankmap[hid] = {c: i + 1 for i, c in enumerate(order)}
        nmap[hid] = len(order)

    recs = []
    for (hid, pos, alt), g in var.groupby(["haplotype_id", "consensus_pos", "alt"]):
        rm = rankmap.get(hid)
        if rm is None:
            continue
        ranks = sorted({rm[c] for c in g["copy_number"] if c in rm})
        k = len(ranks)
        n = nmap[hid]
        if k < 2 or k >= n:          # need >=2 carriers and room for structure
            continue
        R = n_runs(ranks)
        cum, E, Var = runs_stats(k, n)
        p = float(cum[R])                      # P(R <= R_obs)
        A_obs = k - R                          # observed adjacencies
        A_exp = k * (k - 1) / n                # expected adjacencies = k - E[R]
        A_max = k - 1
        ci = (A_obs - A_exp) / (A_max - A_exp) if A_max > A_exp else np.nan
        recs.append((hid, pos, g["region"].iloc[0], k, n, R, E, Var,
                     p, A_obs, A_exp, A_max, ci))
    cols = ["haplotype_id", "pos", "region", "k", "n", "R", "E_R", "Var_R",
            "p", "A_obs", "A_exp", "A_max", "ci"]
    return pd.DataFrame(recs, columns=cols)


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    print("Loading ...", flush=True)
    cp, var = load()
    sites = build_sites(cp, var)
    n_sites = len(sites)
    n_haps = sites["haplotype_id"].nunique()
    print(f"{n_sites:,} multi-copy variant sites (k>=2) across {n_haps} haplotypes")

    # pooled clustering index (0 = random, 1 = contiguous)
    def ci_pool(d):
        num = d["A_obs"].sum() - d["A_exp"].sum()
        den = d["A_max"].sum() - d["A_exp"].sum()
        return num / den if den > 0 else np.nan

    # headline aggregate
    tot_obs = sites["A_obs"].sum()
    tot_exp = sites["A_exp"].sum()
    fold = tot_obs / tot_exp
    ci_overall = ci_pool(sites)
    frac_sig = (sites["p"] < 0.05).mean()
    print(f"\nOverall clustering index (pooled): {ci_overall:.3f}  (0=random, 1=contiguous)")
    print(f"Overall adjacency fold-enrichment: {fold:.2f}x "
          f"(obs {tot_obs:.0f} vs exp {tot_exp:.0f})")
    print(f"Fraction of sites with exact p<0.05: {frac_sig*100:.1f}% (null 5%)")

    # ── clustering index by carrier count k ───────────────────────────────────
    kbin_edges = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 21, 31, 51, 101, 10**9]
    kbin_labels = ["2","3","4","5","6","7","8","9","10","11-15","16-20",
                   "21-30","31-50","51-100","100+"]
    sites["kbin"] = pd.cut(sites["k"], bins=kbin_edges, right=False,
                           labels=kbin_labels)
    rows = []
    for kb, d in sites.groupby("kbin", observed=True):
        rows.append({"kbin": kb, "ci": ci_pool(d),
                     "fold": d["A_obs"].sum()/d["A_exp"].sum() if d["A_exp"].sum()>0 else np.nan,
                     "n": len(d), "frac_sig": (d["p"] < 0.05).mean()})
    perk = pd.DataFrame(rows)
    print("\nClustering index by carrier count k:")
    for _, r in perk.iterrows():
        print(f"  k={str(r['kbin']):7s}: CI={r['ci']:.3f}  fold={r['fold']:.2f}  "
              f"%sig={r['frac_sig']*100:.0f}  n={int(r['n'])}")

    # ── per-haplotype clustering index + Wilcoxon ─────────────────────────────
    rows = []
    for hid, d in sites.groupby("haplotype_id"):
        rows.append({"haplotype_id": hid, "ci_hap": ci_pool(d), "n_sites": len(d)})
    perhap = pd.DataFrame(rows).dropna(subset=["ci_hap"])
    w_stat, w_p = stats.wilcoxon(perhap["ci_hap"], alternative="greater")
    print(f"\nPer-haplotype clustering index: median {perhap['ci_hap'].median():.3f}; "
          f"Wilcoxon vs 0 p={w_p:.2e} ({len(perhap)} haplotypes)")

    # ── clustering index by region and carrier-count class ────────────────────
    sites["kclass"] = np.where(sites["k"] >= 8, "high-k (≥8)", "low-k (2–4)")
    sites.loc[(sites["k"] >= 5) & (sites["k"] < 8), "kclass"] = "mid-k (5–7)"
    perreg = pd.DataFrame([
        {"region": rg, "ci": ci_pool(d), "fold": d["A_obs"].sum()/d["A_exp"].sum(),
         "frac_sig": (d["p"] < 0.05).mean(), "n": len(d)}
        for rg, d in sites.groupby("region")
    ]).set_index("region").reindex(REGION_ORDER).reset_index()
    print("\nBy region (pooled clustering index / fold):")
    for _, r in perreg.iterrows():
        print(f"  {REGION_LABELS[r['region']]:9s}: CI={r['ci']:.3f}  fold={r['fold']:.2f}  "
              f"%sig={r['frac_sig']*100:.0f}  n={int(r['n'])}")

    # region x kclass clustering-index grid
    reg_k = {}
    for rg in REGION_ORDER:
        reg_k[rg] = {}
        for kc in ["low-k (2–4)", "mid-k (5–7)", "high-k (≥8)"]:
            d = sites[(sites["region"] == rg) & (sites["kclass"] == kc)]
            reg_k[rg][kc] = (ci_pool(d) if len(d) >= 10 else np.nan, len(d))
    print("\nRegion × carrier-count clustering index (CI):")
    for rg in REGION_ORDER:
        cells = "  ".join(f"{kc.split()[0]}={reg_k[rg][kc][0]:.2f}(n={reg_k[rg][kc][1]})"
                          for kc in ["low-k (2–4)", "mid-k (5–7)", "high-k (≥8)"])
        print(f"  {REGION_LABELS[rg]:9s}: {cells}")

    # ── figure ─────────────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(16, 12))
    gs = gridspec.GridSpec(2, 2, hspace=0.33, wspace=0.27)

    # A: per-site exact p-value histogram
    ax = fig.add_subplot(gs[0, 0])
    ax.hist(sites["p"], bins=np.linspace(0, 1, 41), color="#4292c6",
            edgecolor="white", linewidth=0.3)
    ax.axhline(n_sites / 40, color="black", ls="--", lw=1.2,
               label="uniform (no clustering)")
    ax.set_xlabel("Exact runs-test p-value  P(R ≤ R_obs)")
    ax.set_ylabel("Number of variant sites")
    ax.set_title("A   Per-site clustering p-values (exact, any k)",
                 fontweight="bold", loc="left")
    ax.text(0.97, 0.95,
            f"{frac_sig*100:.0f}% of sites p<0.05\n(null: 5%)  →  {frac_sig/0.05:.1f}× excess",
            transform=ax.transAxes, ha="right", va="top", fontsize=8.5,
            bbox=dict(boxstyle="round", fc="white", ec="gray", alpha=0.9))
    ax.legend(fontsize=8, loc="upper center")

    # B: clustering index vs carrier count k
    ax2 = fig.add_subplot(gs[0, 1])
    xx = np.arange(len(perk))
    ax2.axhline(0.0, color="black", ls="--", lw=1.2, label="null (random placement)")
    ax2.axhline(1.0, color="gray", ls=":", lw=1.0, label="fully contiguous tract")
    ax2.plot(xx, perk["ci"], "o-", color="#c0392b", lw=2, ms=7)
    for xi, (_, r) in zip(xx, perk.iterrows()):
        ax2.annotate(f"n={int(r['n'])}", (xi, r["ci"]), fontsize=6,
                     textcoords="offset points", xytext=(0, 8), ha="center",
                     rotation=90, color="dimgray")
    ax2.set_xticks(xx)
    ax2.set_xticklabels(perk["kbin"], rotation=45, ha="right", fontsize=8)
    ax2.set_xlabel("Carrier copies per variant within one array (k)")
    ax2.set_ylabel("Clustering index  (0 = random, 1 = contiguous)")
    ax2.set_title("B   Clustering index by carrier count",
                  fontweight="bold", loc="left", fontsize=10)
    ax2.legend(fontsize=8, loc="upper left")
    ax2.set_ylim(-0.1, 1.05)

    # C: per-haplotype clustering index distribution
    ax3 = fig.add_subplot(gs[1, 0])
    ax3.hist(perhap["ci_hap"], bins=30, color="#41ab5d", edgecolor="white",
             linewidth=0.3)
    med = perhap["ci_hap"].median()
    ax3.axvline(0, color="black", ls="--", lw=1.2, label="random (CI=0)")
    ax3.axvline(med, color="#c0392b", lw=2, label=f"median = {med:.2f}")
    ax3.set_xlabel("Per-haplotype clustering index\n(0 = random, 1 = fully contiguous)")
    ax3.set_ylabel("Number of haplotypes")
    ax3.set_title("C   Clustering per haplotype (independent units)",
                  fontweight="bold", loc="left")
    ax3.text(0.03, 0.95, f"Wilcoxon vs 0:\np = {w_p:.1e}\nn = {len(perhap)} haplotypes",
             transform=ax3.transAxes, va="top", fontsize=9,
             bbox=dict(boxstyle="round", fc="white", ec="gray", alpha=0.9))
    ax3.legend(fontsize=8, loc="upper right")

    # D: clustering index by region and carrier-count class
    ax4 = fig.add_subplot(gs[1, 1])
    kclasses = ["low-k (2–4)", "mid-k (5–7)", "high-k (≥8)"]
    kc_colors = ["#9ecae1", "#fdae6b", "#c0392b"]
    xr = np.arange(len(REGION_ORDER)); w = 0.26
    ax4.axhline(0.0, color="black", ls="--", lw=1.0)
    for j, kc in enumerate(kclasses):
        vals = [reg_k[rg][kc][0] for rg in REGION_ORDER]
        ns   = [reg_k[rg][kc][1] for rg in REGION_ORDER]
        ax4.bar(xr + (j-1)*w, vals, w, color=kc_colors[j], label=kc,
                edgecolor="black", linewidth=0.5)
        for xi, v, nn in zip(xr + (j-1)*w, vals, ns):
            if not np.isnan(v):
                ax4.text(xi, v + (0.02 if v >= 0 else -0.05),
                         f"n={nn}", ha="center",
                         va="bottom" if v >= 0 else "top", fontsize=6, rotation=90)
    ax4.set_xticks(xr)
    ax4.set_xticklabels([REGION_LABELS[r] for r in REGION_ORDER])
    ax4.set_ylabel("Clustering index (0 = random, 1 = contiguous)")
    ax4.set_title("D   Clustering index by region × carrier count",
                  fontweight="bold", loc="left", fontsize=10)
    ax4.legend(fontsize=8, title="carrier count")

    fig.suptitle(
        f"5S rDNA within-array variant clustering — exact runs test  "
        f"({n_sites:,} sites, {n_haps} haplotypes)",
        fontsize=12, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    out = FIG_DIR / "31_variant_clustering_runs.pdf"
    fig.savefig(out, bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"\nSaved {out}")


if __name__ == "__main__":
    main()
