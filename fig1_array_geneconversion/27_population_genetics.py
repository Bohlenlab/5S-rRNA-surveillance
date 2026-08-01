#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 27_population_genetics.py — computes 5S rDNA population-genetic summaries across HPRC superpopulations (copy number, folded site-frequency spectrum, per-site Hudson Fst, Tajima's D, and variant sharing) and writes the figures and tables.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
27_population_genetics.py

Population genetics of 5S rDNA array variation across HPRC superpopulations,
using the consensus_t2t variant pool (ref = population-major allele).

Unit of analysis: the HAPLOTYPE. Copies within a haplotype are treated as
non-independent, so a variant's population frequency is the fraction of
HAPLOTYPES carrying it in >=1 copy.

Figures (figures/):
  27a_copy_number_by_population.pdf  — array length distribution & differentiation
  27b_folded_sfs.pdf                 — folded site-frequency spectrum by region
  27c_fst.pdf                        — per-site Fst between superpopulations
  27d_selection_tests.pdf            — Tajima's D by region; private/shared variants

Tables:
  tables/27_fst_per_site.tsv
  tables/27_population_summary.tsv
"""

import os
import sqlite3
from itertools import combinations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
from scipy import stats

HPRC    = Path(os.environ.get("FIVES_OUT", "output"))
DB_PATH = Path(os.environ.get("FIVES_DB", "5S_rDNA.db"))
FIG_DIR = HPRC / "figures" / "06_population_genetics"
TAB_DIR = HPRC / "tables"
FIG_DIR.mkdir(exist_ok=True); TAB_DIR.mkdir(exist_ok=True)

SUPERPOPS = ["AFR", "EUR", "SAS", "EAS", "AMR"]
POP_COLORS = {"AFR": "#E41A1C", "EUR": "#4DAF4A", "SAS": "#984EA3",
              "EAS": "#377EB8", "AMR": "#FF7F00"}
REGION_ORDER = ["nts_pre", "gene", "nts_post"]
REGION_COLORS = {"nts_pre": "#4575b4", "gene": "#ff9900", "nts_post": "#d73027"}
REGION_LABELS = {"nts_pre": "NTS-pre", "gene": "5S gene", "nts_post": "NTS-post"}
REGION_LEN = {"nts_pre": 629, "gene": 119, "nts_post": 1420}


# ── data ──────────────────────────────────────────────────────────────────────

def load():
    con = sqlite3.connect(DB_PATH)

    # per-haplotype metadata + interior copy count
    haps = pd.read_sql_query("""
        SELECT h.haplotype_id, a.sample_id, a.superpopulation AS superpop,
               a.population AS pop,
               SUM(CASE WHEN c.border_note='interior' THEN 1 ELSE 0 END) AS n_interior
        FROM haplotype h
        JOIN assembly a USING(assembly_id)
        JOIN copy c USING(haplotype_id)
        WHERE a.cohort IN ('HPRC_Year1','HPRC_Release2')
        GROUP BY h.haplotype_id
    """, con)

    # all interior SNPs (consensus_t2t), including masked
    snps = pd.read_sql_query("""
        SELECT h.haplotype_id, a.superpopulation AS superpop,
               v.consensus_pos, v.alt, v.region
        FROM variant v
        JOIN copy c USING(copy_id)
        JOIN haplotype h USING(haplotype_id)
        JOIN assembly a USING(assembly_id)
        WHERE c.border_note='interior'
          AND a.cohort IN ('HPRC_Year1','HPRC_Release2')
          AND v.alignment_source='consensus_t2t' AND v.var_type='snp'
          AND length(v.alt)=1 AND v.alt IN ('A','C','G','T')
    """, con)
    con.close()

    haps = haps[haps["superpop"].isin(SUPERPOPS)].copy()
    snps = snps[snps["superpop"].isin(SUPERPOPS)].copy()
    return haps, snps


# ── haplotype-level allele frequencies ────────────────────────────────────────

def build_site_table(haps, snps):
    """
    Per variant site (pos, alt): number of haplotypes carrying it, overall and
    per superpopulation. Returns a DataFrame.
    """
    n_hap_total = haps["haplotype_id"].nunique()
    n_by_pop = haps.groupby("superpop")["haplotype_id"].nunique().to_dict()

    # haplotypes carrying each (pos, alt)
    car = (snps.groupby(["consensus_pos", "alt", "region"])
           .agg(hap_set=("haplotype_id", lambda x: frozenset(x)))
           .reset_index())
    car["n_total"] = car["hap_set"].apply(len)
    car["freq"] = car["n_total"] / n_hap_total

    # per-pop carrier counts
    hap_pop = haps.set_index("haplotype_id")["superpop"].to_dict()
    for pop in SUPERPOPS:
        car[f"n_{pop}"] = car["hap_set"].apply(
            lambda s: sum(1 for h in s if hap_pop.get(h) == pop))
        car[f"f_{pop}"] = car[f"n_{pop}"] / n_by_pop[pop]
    return car, n_hap_total, n_by_pop


# ── Hudson's Fst per site ─────────────────────────────────────────────────────

def hudson_fst(p1, n1, p2, n2):
    """Hudson's Fst estimator for one site between two populations.
    p = allele frequency, n = number of haplotypes sampled."""
    if n1 < 2 or n2 < 2:
        return np.nan
    # numerator = (p1-p2)^2 - p1(1-p1)/(n1-1) - p2(1-p2)/(n2-1)
    num = (p1 - p2)**2 - p1*(1-p1)/(n1-1) - p2*(1-p2)/(n2-1)
    den = p1*(1-p2) + p2*(1-p1)
    if den == 0:
        return np.nan
    return num / den


# ── Tajima's D ────────────────────────────────────────────────────────────────

def tajimas_d(allele_counts, n):
    """
    Tajima's D from a list of derived-allele counts (per segregating site)
    in a sample of n haplotypes.
    """
    S = len(allele_counts)
    if S == 0 or n < 4:
        return np.nan
    a1 = sum(1.0/i for i in range(1, n))
    a2 = sum(1.0/i**2 for i in range(1, n))
    b1 = (n+1)/(3*(n-1))
    b2 = 2*(n**2 + n + 3)/(9*n*(n-1))
    c1 = b1 - 1/a1
    c2 = b2 - (n+2)/(a1*n) + a2/a1**2
    e1 = c1/a1
    e2 = c2/(a1**2 + a2)
    theta_w = S / a1
    # pi = mean pairwise differences = sum over sites of 2*p*(1-p)*n/(n-1)
    pi = 0.0
    for c in allele_counts:
        p = c / n
        pi += 2 * p * (1-p) * n/(n-1)
    var = e1*S + e2*S*(S-1)
    if var <= 0:
        return np.nan
    return (pi - theta_w) / np.sqrt(var)


# ── Figure A: copy number ─────────────────────────────────────────────────────

def fig_copy_number(haps):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    ax = axes[0]
    data = [haps[haps["superpop"]==p]["n_interior"].values for p in SUPERPOPS]
    parts = ax.violinplot(data, positions=range(len(SUPERPOPS)),
                          showmedians=True, showextrema=False, widths=0.8)
    for b, p in zip(parts["bodies"], SUPERPOPS):
        b.set_facecolor(POP_COLORS[p]); b.set_alpha(0.65)
    parts["cmedians"].set_color("black")
    # jittered points
    for i, p in enumerate(SUPERPOPS):
        y = haps[haps["superpop"]==p]["n_interior"].values
        x = np.random.normal(i, 0.06, len(y))
        ax.scatter(x, y, s=8, c="black", alpha=0.35, zorder=3)
    ax.set_xticks(range(len(SUPERPOPS)))
    ax.set_xticklabels([f"{p}\n(n={len(haps[haps.superpop==p])})" for p in SUPERPOPS])
    ax.set_ylabel("Interior copies per haplotype")
    ax.set_title("A   Array copy number by superpopulation", fontweight="bold", loc="left")

    # Kruskal-Wallis
    H, pval = stats.kruskal(*data)
    ax.text(0.02, 0.96, f"Kruskal-Wallis H={H:.1f}, p={pval:.3f}",
            transform=ax.transAxes, va="top", fontsize=9,
            bbox=dict(boxstyle="round", fc="white", ec="gray", alpha=0.9))

    # right: ECDF
    ax2 = axes[1]
    for p in SUPERPOPS:
        y = np.sort(haps[haps["superpop"]==p]["n_interior"].values)
        ax2.step(y, np.arange(1, len(y)+1)/len(y), where="post",
                 color=POP_COLORS[p], lw=2, label=p)
    ax2.set_xlabel("Interior copies per haplotype")
    ax2.set_ylabel("Cumulative fraction")
    ax2.set_title("B   Copy-number ECDF", fontweight="bold", loc="left")
    ax2.legend(fontsize=9, title="Superpop")

    med_overall = haps["n_interior"].median()
    fig.suptitle(f"5S rDNA Array Copy Number — {len(haps)} haplotypes "
                 f"(median {med_overall:.0f} copies)", fontsize=13, fontweight="bold")
    fig.tight_layout()
    out = FIG_DIR / "27a_copy_number_by_population.pdf"
    fig.savefig(out, bbox_inches="tight"); plt.close(fig)
    print(f"  Saved {out.name}  (KW p={pval:.3g})")
    return {"kw_H": H, "kw_p": pval}


# ── Figure B: folded SFS ──────────────────────────────────────────────────────

def fig_folded_sfs(site, n_hap_total):
    """Folded SFS: minor-allele-frequency spectrum. ref is already the major
    allele globally, but per-site the carrier set can exceed 50%, so fold."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # minor allele count = min(carriers, total - carriers)
    site = site.copy()
    site["mac"] = np.minimum(site["n_total"], n_hap_total - site["n_total"])

    # Panel A: overall folded SFS
    ax = axes[0]
    bins = np.arange(0.5, min(site["mac"].max(), n_hap_total//2) + 1.5)
    ax.hist(site["mac"], bins=bins, color="#555", edgecolor="white", linewidth=0.2)
    ax.set_yscale("log")
    ax.set_xlabel("Minor allele count (haplotypes)")
    ax.set_ylabel("Number of sites")
    ax.set_title("A   Folded SFS (all regions)", fontweight="bold", loc="left")
    sing = (site["mac"]==1).sum()
    ax.text(0.55, 0.9, f"Sites: {len(site):,}\nSingletons: {sing:,} "
            f"({sing/len(site)*100:.0f}%)", transform=ax.transAxes, va="top",
            fontsize=9, bbox=dict(boxstyle="round", fc="white", ec="gray", alpha=0.9))

    # Panel B: folded SFS by region (normalized)
    ax2 = axes[1]
    maxmac = 20
    x = np.arange(1, maxmac+1)
    for region in REGION_ORDER:
        s = site[site["region"]==region]["mac"]
        if len(s)==0: continue
        h = np.array([(s==k).sum() for k in x], float)
        h = h/h.sum()
        ax2.plot(x, h, "o-", color=REGION_COLORS[region], ms=4, lw=1.5,
                 label=f"{REGION_LABELS[region]} (n={len(s):,})")
    ax2.set_yscale("log")
    ax2.set_xlabel("Minor allele count (haplotypes)")
    ax2.set_ylabel("Fraction of sites")
    ax2.set_title("B   Folded SFS by region", fontweight="bold", loc="left")
    ax2.legend(fontsize=8)

    # Panel C: proportion rare (MAC<=2) vs common by region
    ax3 = axes[2]
    fracs = []
    for region in REGION_ORDER:
        s = site[site["region"]==region]["mac"]
        rare = (s<=2).mean() if len(s) else 0
        common = (s >= 0.05*n_hap_total).mean() if len(s) else 0
        fracs.append((rare, common))
    x = np.arange(len(REGION_ORDER))
    rares = [f[0] for f in fracs]; commons = [f[1] for f in fracs]
    ax3.bar(x-0.2, rares, 0.4, color="#9970ab", label="Rare (MAC≤2)")
    ax3.bar(x+0.2, commons, 0.4, color="#5aae61", label="Common (≥5%)")
    ax3.set_xticks(x); ax3.set_xticklabels([REGION_LABELS[r] for r in REGION_ORDER])
    ax3.set_ylabel("Fraction of sites")
    ax3.set_title("C   Rare vs common by region", fontweight="bold", loc="left")
    ax3.legend(fontsize=9)

    fig.suptitle("5S rDNA Folded Site-Frequency Spectrum", fontsize=13, fontweight="bold")
    fig.tight_layout()
    out = FIG_DIR / "27b_folded_sfs.pdf"
    fig.savefig(out, bbox_inches="tight"); plt.close(fig)
    print(f"  Saved {out.name}")


# ── Figure C: Fst ─────────────────────────────────────────────────────────────

def fig_fst(site, n_by_pop):
    # Compute per-site mean pairwise Fst across the 10 superpop pairs
    pairs = list(combinations(SUPERPOPS, 2))
    fst_rows = []
    for _, row in site.iterrows():
        vals = []
        for a, b in pairs:
            f = hudson_fst(row[f"f_{a}"], n_by_pop[a], row[f"f_{b}"], n_by_pop[b])
            if not np.isnan(f):
                vals.append(max(f, 0.0))  # clamp small negatives to 0
        if vals:
            fst_rows.append({"pos": row["consensus_pos"], "alt": row["alt"],
                             "region": row["region"], "mean_fst": np.mean(vals),
                             "freq": row["freq"]})
    fst = pd.DataFrame(fst_rows)
    fst.to_csv(TAB_DIR / "27_fst_per_site.tsv", sep="\t", index=False)

    fig = plt.figure(figsize=(16, 10))
    gs = gridspec.GridSpec(2, 2, hspace=0.35, wspace=0.3)

    # A: Fst distribution
    ax = fig.add_subplot(gs[0,0])
    ax.hist(fst["mean_fst"], bins=50, color="#666", edgecolor="white", linewidth=0.2)
    ax.axvline(fst["mean_fst"].mean(), color="red", ls="--",
               label=f"mean={fst['mean_fst'].mean():.3f}")
    ax.set_xlabel("Mean pairwise Fst (Hudson)")
    ax.set_ylabel("Number of sites"); ax.set_yscale("log")
    ax.set_title("A   Per-site Fst distribution", fontweight="bold", loc="left")
    ax.legend(fontsize=9)

    # B: Fst by region (boxplot)
    ax2 = fig.add_subplot(gs[0,1])
    data = [fst[fst["region"]==r]["mean_fst"].values for r in REGION_ORDER]
    bp = ax2.boxplot(data, positions=range(len(REGION_ORDER)), widths=0.6,
                     patch_artist=True, showfliers=False)
    for patch, r in zip(bp["boxes"], REGION_ORDER):
        patch.set_facecolor(REGION_COLORS[r]); patch.set_alpha(0.6)
    ax2.set_xticks(range(len(REGION_ORDER)))
    ax2.set_xticklabels([REGION_LABELS[r] for r in REGION_ORDER])
    ax2.set_ylabel("Mean pairwise Fst")
    ax2.set_title("B   Fst by region", fontweight="bold", loc="left")

    # C: Fst vs frequency (Manhattan-style by position)
    ax3 = fig.add_subplot(gs[1,:])
    for r in REGION_ORDER:
        s = fst[fst["region"]==r]
        ax3.scatter(s["pos"], s["mean_fst"], s=12, c=REGION_COLORS[r],
                    alpha=0.5, label=REGION_LABELS[r])
    # highlight top differentiated sites
    top = fst.nlargest(10, "mean_fst")
    ax3.scatter(top["pos"], top["mean_fst"], s=45, facecolors="none",
                edgecolors="black", linewidths=1.0, zorder=5)
    for _, t in top.iterrows():
        ax3.annotate(f"{int(t['pos'])}", (t["pos"], t["mean_fst"]),
                     fontsize=7, ha="center", va="bottom")
    ax3.axvspan(630, 748, color="#ffe7b3", alpha=0.4, zorder=0)
    ax3.set_xlabel("Position in repeat unit")
    ax3.set_ylabel("Mean pairwise Fst")
    ax3.set_title("C   Fst landscape along the repeat unit (top 10 labeled)",
                  fontweight="bold", loc="left")
    ax3.legend(fontsize=8)
    ax3.set_xlim(0, 2168)

    fig.suptitle("5S rDNA Population Differentiation (Fst across 5 superpopulations)",
                 fontsize=13, fontweight="bold")
    out = FIG_DIR / "27c_fst.pdf"
    fig.savefig(out, bbox_inches="tight"); plt.close(fig)
    print(f"  Saved {out.name}  (mean Fst={fst['mean_fst'].mean():.4f})")
    return fst


# ── Figure D: selection tests + private variants ──────────────────────────────

def fig_selection(site, snps, haps, n_hap_total, n_by_pop):
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # Tajima's D by region (using haplotype-level allele counts)
    ax = axes[0]
    tajd = {}
    for region in REGION_ORDER:
        s = site[site["region"]==region]
        # minor allele counts
        mac = np.minimum(s["n_total"], n_hap_total - s["n_total"]).values
        tajd[region] = tajimas_d(list(mac), n_hap_total)
    bars = ax.bar(range(len(REGION_ORDER)), [tajd[r] for r in REGION_ORDER],
                  color=[REGION_COLORS[r] for r in REGION_ORDER],
                  edgecolor="black", linewidth=0.8)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xticks(range(len(REGION_ORDER)))
    ax.set_xticklabels([REGION_LABELS[r] for r in REGION_ORDER])
    ax.set_ylabel("Tajima's D")
    ax.set_title("A   Tajima's D by region", fontweight="bold", loc="left")
    for i, r in enumerate(REGION_ORDER):
        ax.text(i, tajd[r], f"{tajd[r]:.2f}", ha="center",
                va="bottom" if tajd[r]>=0 else "top", fontsize=9)

    # Private vs shared variants across superpopulations
    ax2 = axes[1]
    pop_presence = site[[f"n_{p}" for p in SUPERPOPS]].gt(0).sum(axis=1)
    counts = pop_presence.value_counts().sort_index()
    ax2.bar(counts.index, counts.values, color="#4393c3",
            edgecolor="black", linewidth=0.6)
    ax2.set_xlabel("Number of superpopulations carrying variant")
    ax2.set_ylabel("Number of sites")
    ax2.set_title("B   Variant sharing across superpopulations",
                  fontweight="bold", loc="left")
    priv = (pop_presence==1).sum()
    allp = (pop_presence==5).sum()
    ax2.text(0.5, 0.9, f"Private (1 pop): {priv:,} ({priv/len(site)*100:.0f}%)\n"
             f"Shared (all 5): {allp:,} ({allp/len(site)*100:.0f}%)",
             transform=ax2.transAxes, va="top", fontsize=9,
             bbox=dict(boxstyle="round", fc="white", ec="gray", alpha=0.9))

    # Private variant count by superpopulation (normalized by sample size)
    ax3 = axes[2]
    priv_by_pop = {}
    private_sites = site[pop_presence==1]
    for p in SUPERPOPS:
        priv_by_pop[p] = (private_sites[f"n_{p}"]>0).sum()
    # normalize by n haplotypes
    vals = [priv_by_pop[p]/n_by_pop[p] for p in SUPERPOPS]
    ax3.bar(range(len(SUPERPOPS)), vals,
            color=[POP_COLORS[p] for p in SUPERPOPS], edgecolor="black", linewidth=0.6)
    ax3.set_xticks(range(len(SUPERPOPS)))
    ax3.set_xticklabels(SUPERPOPS)
    ax3.set_ylabel("Private variants per haplotype")
    ax3.set_title("C   Population-private variants (size-normalized)",
                  fontweight="bold", loc="left")
    for i, p in enumerate(SUPERPOPS):
        ax3.text(i, vals[i], f"{priv_by_pop[p]}", ha="center", va="bottom", fontsize=8)

    fig.suptitle("5S rDNA Selection & Population Structure", fontsize=13, fontweight="bold")
    fig.tight_layout()
    out = FIG_DIR / "27d_selection_tests.pdf"
    fig.savefig(out, bbox_inches="tight"); plt.close(fig)
    print(f"  Saved {out.name}")
    print(f"    Tajima's D: " + ", ".join(f"{r}={tajd[r]:.2f}" for r in REGION_ORDER))
    print(f"    Private variants: {priv:,}/{len(site):,} ({priv/len(site)*100:.0f}%); "
          f"shared by all 5: {allp:,}")
    return tajd


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    print("Loading data...")
    haps, snps = load()
    print(f"  {haps['haplotype_id'].nunique()} haplotypes, {len(snps):,} SNP calls")
    print("  Superpop sizes: " +
          ", ".join(f"{p}={haps[haps.superpop==p]['haplotype_id'].nunique()}" for p in SUPERPOPS))

    site, n_hap_total, n_by_pop = build_site_table(haps, snps)
    print(f"  {len(site):,} unique variant sites\n")

    # population summary table
    summ = haps.groupby("superpop").agg(
        n_haplotypes=("haplotype_id","nunique"),
        median_copies=("n_interior","median"),
        mean_copies=("n_interior","mean"),
        sd_copies=("n_interior","std")).reset_index()
    summ.to_csv(TAB_DIR / "27_population_summary.tsv", sep="\t", index=False)

    print("--- A: copy number ---");      fig_copy_number(haps)
    print("--- B: folded SFS ---");        fig_folded_sfs(site, n_hap_total)
    print("--- C: Fst ---");               fst = fig_fst(site, n_by_pop)
    print("--- D: selection tests ---");   fig_selection(site, snps, haps, n_hap_total, n_by_pop)

    # report top differentiated sites
    print("\nTop 10 most population-differentiated sites:")
    print(fst.nlargest(10, "mean_fst")[["pos","alt","region","mean_fst","freq"]].to_string(index=False))
    print(f"\nFigures + tables saved.")


if __name__ == "__main__":
    main()
