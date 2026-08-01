#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 39_substitution_normalized_cpg.py — Computes the base-content-normalized
# substitution spectrum and the CpG frequency / CpG-mutation profile along the
# 5S rDNA repeat unit.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
39_substitution_normalized_cpg.py

(1) Substitution spectrum normalized to nucleotide content — for the full copy
    and for each sub-feature (NTS-pre / gene / NTS-post). Each pyrimidine-collapsed
    class rate = variants / eligible reference bases / haplotype:
        C>A,C>G,C>T  per (C+G) ;  T>A,T>C,T>G  per (T+A)

(2) CpG frequency, CpG-mutation frequency, and their ratio (per-CpG mutability)
    along the 2168 bp repeat unit. CpG mutation = C>T at the C, or G>A at the G,
    of a CpG (5mC-deamination signature), counted at the haplotype level
    (distinct haplotypes mutated).

Reference = population consensus (consensus_reference.json), matching the
polarized consensus_t2t variants.

Outputs:
  39a_substitution_spectrum_normalized.pdf
  39b_cpg_along_copy.pdf
"""

import json, os, sqlite3
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd

HPRC    = Path(os.environ.get("FIVES_OUT", "output"))
DB_PATH = Path(os.environ.get("FIVES_DB", "5S_rDNA.db"))
FIG_DIR = HPRC / "figures" / "05_gene_conversion"
CONS    = Path(os.environ.get("FIVES_DATA", "data")) / "consensus_reference.json"

T2T_LEN = 2168
GENE_START, GENE_END = 630, 748
REGIONS = {"nts_pre": (1, 629), "gene": (630, 748), "nts_post": (749, 2168)}
REGION_ORDER = ["nts_pre", "gene", "nts_post"]
REGION_LABELS = {"nts_pre": "NTS-pre", "gene": "5S gene", "nts_post": "NTS-post"}
REGION_COLORS = {"nts_pre": "#4575b4", "gene": "#ff9900", "nts_post": "#d73027"}

MUT6 = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
MUT6_COLORS = ["#1f77b4", "#aec7e8", "#d62728", "#ff9896", "#2ca02c", "#98df8a"]
COMP = str.maketrans("ACGT", "TGCA")
PYR = set("CT")


def pyr_class(ref, alt):
    if ref not in "ACGT" or alt not in "ACGT":
        return None
    if ref not in PYR:
        ref = ref.translate(COMP); alt = alt.translate(COMP)
    return f"{ref}>{alt}"


def load():
    con = sqlite3.connect(DB_PATH)
    n_haps = con.execute("""
        SELECT COUNT(DISTINCT h.haplotype_id) FROM haplotype h
        JOIN assembly a USING(assembly_id)
        WHERE a.cohort IN ('HPRC_Year1','HPRC_Release2')""").fetchone()[0]
    # haplotype-level variant incidence per (pos, ref, alt)
    v = pd.read_sql_query("""
        SELECT v.consensus_pos AS pos, v.ref, v.alt, v.region,
               COUNT(DISTINCT v.copy_id) AS n_copies,
               COUNT(DISTINCT h.haplotype_id) AS n_haps
        FROM variant v JOIN copy c USING(copy_id)
        JOIN haplotype h USING(haplotype_id) JOIN assembly a USING(assembly_id)
        WHERE c.border_note='interior' AND a.cohort IN ('HPRC_Year1','HPRC_Release2')
          AND v.alignment_source='consensus_t2t' AND v.var_type='snp'
          AND length(v.ref)=1 AND length(v.alt)=1
          AND v.ref IN ('A','C','G','T') AND v.alt IN ('A','C','G','T')
        GROUP BY v.consensus_pos, v.ref, v.alt, v.region
    """, con)
    con.close()
    cons = json.load(open(CONS))["consensus"].upper()
    return v, cons, n_haps


# ── (1) normalized substitution spectrum ──────────────────────────────────────

def base_counts(seq, lo, hi):
    sub = seq[lo-1:hi]
    return {b: sub.count(b) for b in "ACGT"}


def normalized_spectrum(v, cons, n_haps):
    v = v.copy()
    v["cls"] = [pyr_class(r, a) for r, a in zip(v["ref"], v["alt"])]
    v = v[v["cls"].notna()]
    # eligible-base denominators per region and full
    def denom(bc, cls):
        return (bc["C"] + bc["G"]) if cls.startswith("C") else (bc["T"] + bc["A"])
    scopes = {"full": (1, T2T_LEN)} | REGIONS
    rates = {}            # scope -> {class: rate per eligible base per hap}
    for scope, (lo, hi) in scopes.items():
        bc = base_counts(cons, lo, hi)
        if scope == "full":
            sub = v
        else:
            sub = v[v["region"] == scope]
        r = {}
        for cls in MUT6:
            n = sub.loc[sub["cls"] == cls, "n_haps"].sum()   # hap-level incidence
            r[cls] = n / max(denom(bc, cls), 1) / n_haps
        rates[scope] = r
    return rates, {s: base_counts(cons, *rg) for s, rg in scopes.items()}


def fig_spectrum(rates, bcs, n_haps):
    fig = plt.figure(figsize=(16, 11))
    gs = gridspec.GridSpec(2, 2, hspace=0.34, wspace=0.27)

    # A: full-copy normalized spectrum
    ax = fig.add_subplot(gs[0, 0])
    vals = [rates["full"][c] for c in MUT6]
    ax.bar(range(6), vals, color=MUT6_COLORS, edgecolor="black", linewidth=0.5)
    ax.set_xticks(range(6)); ax.set_xticklabels(MUT6)
    ax.set_ylabel("Substitutions per eligible base per haplotype")
    ax.set_title("A   Full-copy spectrum, normalized to base content",
                 fontweight="bold", loc="left")

    # B: per-region normalized spectrum
    ax2 = fig.add_subplot(gs[0, 1])
    x = np.arange(6); w = 0.26
    for j, rg in enumerate(REGION_ORDER):
        ax2.bar(x + (j-1)*w, [rates[rg][c] for c in MUT6], w,
                color=REGION_COLORS[rg], edgecolor="black", linewidth=0.4,
                label=REGION_LABELS[rg])
    ax2.set_xticks(x); ax2.set_xticklabels(MUT6)
    ax2.set_ylabel("Substitutions per eligible base per haplotype")
    ax2.set_title("B   Per-region spectrum, each normalized to its own base content",
                  fontweight="bold", loc="left", fontsize=10)
    ax2.legend(fontsize=8)

    # C: base composition per region (why normalization matters)
    ax3 = fig.add_subplot(gs[1, 0])
    xr = np.arange(len(REGION_ORDER)); w = 0.2
    for j, b in enumerate("ACGT"):
        fr = [bcs[rg][b]/sum(bcs[rg].values()) for rg in REGION_ORDER]
        ax3.bar(xr + (j-1.5)*w, fr, w, label=b, edgecolor="white", linewidth=0.3)
    ax3.set_xticks(xr); ax3.set_xticklabels([REGION_LABELS[r] for r in REGION_ORDER])
    ax3.set_ylabel("Base fraction in consensus")
    ax3.set_title("C   Base composition by region (GC differs → normalize)",
                  fontweight="bold", loc="left", fontsize=10)
    ax3.legend(fontsize=8, ncol=4)

    # D: overall substitution rate per eligible base by region + CpG share
    ax4 = fig.add_subplot(gs[1, 1])
    tot = [sum(rates[rg].values()) for rg in REGION_ORDER]
    ct  = [rates[rg]["C>T"] for rg in REGION_ORDER]
    ax4.bar(xr, tot, 0.5, color=[REGION_COLORS[r] for r in REGION_ORDER],
            edgecolor="black", linewidth=0.5, label="all classes")
    ax4.bar(xr, ct, 0.5, color="#7a0000", edgecolor="black", linewidth=0.5,
            alpha=0.85, label="C>T component")
    ax4.set_xticks(xr); ax4.set_xticklabels([REGION_LABELS[r] for r in REGION_ORDER])
    ax4.set_ylabel("Substitutions per eligible base per haplotype")
    ax4.set_title("D   Total normalized substitution rate by region",
                  fontweight="bold", loc="left", fontsize=10)
    ax4.legend(fontsize=8)

    fig.suptitle(f"5S rDNA substitution spectrum normalized to nucleotide content "
                 f"({n_haps} haplotypes)", fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out = FIG_DIR / "39a_substitution_spectrum_normalized.pdf"
    fig.savefig(out, bbox_inches="tight", dpi=150); plt.close(fig)
    print(f"Saved {out.name}")
    # print table
    print("\nNormalized substitution rate (x1e3 per eligible base per hap):")
    print(f"  {'class':5s} " + "  ".join(f"{REGION_LABELS[r]:>9s}" for r in ['nts_pre','gene','nts_post']) + f"  {'FULL':>9s}")
    for c in MUT6:
        print(f"  {c:5s} " + "  ".join(f"{rates[r][c]*1e3:9.3f}" for r in ['nts_pre','gene','nts_post']) + f"  {rates['full'][c]*1e3:9.3f}")


# ── (2) CpG along copy length ─────────────────────────────────────────────────

def cpg_profile(v, cons, n_haps, BIN=40):
    # CpG dinucleotides: 1-based C-position p where cons[p-1]=='C', cons[p]=='G'
    cpg_C = [p for p in range(1, T2T_LEN) if cons[p-1] == "C" and cons[p] == "G"]
    cpg_set_C = set(cpg_C)                       # C positions of CpGs
    cpg_set_G = set(p+1 for p in cpg_C)          # G positions of CpGs

    # haplotype-level mutation incidence at each position for deamination alleles
    vd = v.set_index(["pos", "ref", "alt"])["n_haps"]
    def nhap(pos, ref, alt):
        try:    return int(vd.loc[(pos, ref, alt)].sum() if hasattr(vd.loc[(pos,ref,alt)],'sum') else vd.loc[(pos,ref,alt)])
        except KeyError: return 0

    # per CpG site (indexed by C position): deamination incidence (C>T + G>A)
    cpg_mut = {}
    for p in cpg_C:
        cpg_mut[p] = nhap(p, "C", "T") + nhap(p+1, "G", "A")

    # non-CpG C>T incidence per C/G position (for comparison)
    noncpg_pos, noncpg_mut = [], []
    for p in range(1, T2T_LEN+1):
        b = cons[p-1]
        if b == "C" and p not in cpg_set_C:
            noncpg_pos.append(p); noncpg_mut.append(nhap(p, "C", "T"))
        elif b == "G" and p not in cpg_set_G:
            noncpg_pos.append(p); noncpg_mut.append(nhap(p, "G", "A"))

    # window aggregation
    edges = np.arange(0, T2T_LEN + BIN, BIN)
    centers = edges[:-1] + BIN/2
    cpg_freq = np.zeros(len(centers))      # CpG sites per window
    cpg_mfreq = np.zeros(len(centers))     # hap-incidence per window / n_haps
    ncpg_sites = np.zeros(len(centers))
    ncpg_mfreq = np.zeros(len(centers))
    for p in cpg_C:
        b = min(int(p // BIN), len(centers)-1)
        cpg_freq[b] += 1
        cpg_mfreq[b] += cpg_mut[p] / n_haps
    for p, m in zip(noncpg_pos, noncpg_mut):
        b = min(int(p // BIN), len(centers)-1)
        ncpg_sites[b] += 1
        ncpg_mfreq[b] += m / n_haps
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(cpg_freq > 0, cpg_mfreq / cpg_freq, np.nan)
        ncpg_ratio = np.where(ncpg_sites > 0, ncpg_mfreq / ncpg_sites, np.nan)
    # genome-wide CpG C>T enrichment
    tot_cpg_mut = sum(cpg_mut.values())          # C>T + G>A at CpG (both cytosines)
    n_cpg_cytosines = 2 * len(cpg_C)
    tot_ncpg_mut = sum(noncpg_mut)
    n_ncpg_cyt = len(noncpg_pos)
    rate_cpg = tot_cpg_mut / max(n_cpg_cytosines, 1) / n_haps
    rate_ncpg = tot_ncpg_mut / max(n_ncpg_cyt, 1) / n_haps
    enrich = rate_cpg / rate_ncpg if rate_ncpg > 0 else float("nan")
    return dict(centers=centers, cpg_freq=cpg_freq, cpg_mfreq=cpg_mfreq,
                ratio=ratio, ncpg_ratio=ncpg_ratio, BIN=BIN,
                n_cpg=len(cpg_C), enrich=enrich,
                rate_cpg=rate_cpg, rate_ncpg=rate_ncpg)


def fig_cpg(prof, n_haps):
    c = prof["centers"]; BIN = prof["BIN"]
    fig = plt.figure(figsize=(15, 12))
    gs = gridspec.GridSpec(3, 1, hspace=0.18)

    def shade(ax):
        ax.axvspan(GENE_START, GENE_END, color="#aec6cf", alpha=0.4, zorder=0)
        for x, lbl in [(GENE_START/2,"NTS-pre"),((GENE_START+GENE_END)/2,"gene"),
                       ((GENE_END+T2T_LEN)/2,"NTS-post")]:
            ax.text(x, ax.get_ylim()[1], lbl, fontsize=7, ha="center", va="bottom",
                    color="dimgrey")

    ax = fig.add_subplot(gs[0]);
    ax.bar(c, prof["cpg_freq"], width=BIN*0.9, color="#2c7fb8")
    ax.set_xlim(0, T2T_LEN); ax.set_ylabel(f"CpG sites / {BIN} bp")
    ax.set_xticklabels([]); shade(ax)
    ax.set_title(f"A   CpG frequency along the repeat unit  (total {prof['n_cpg']} CpGs)",
                 fontweight="bold", loc="left")

    ax2 = fig.add_subplot(gs[1])
    ax2.bar(c, prof["cpg_mfreq"], width=BIN*0.9, color="#d62728")
    ax2.set_xlim(0, T2T_LEN); ax2.set_ylabel(f"CpG mutations / {BIN} bp\n(hap-incidence / hap)")
    ax2.set_xticklabels([]); shade(ax2)
    ax2.set_title("B   CpG-mutation frequency (C>T / G>A at CpG, deamination)",
                  fontweight="bold", loc="left")

    ax3 = fig.add_subplot(gs[2])
    ax3.plot(c, prof["ratio"], "-o", color="#7a0000", ms=3, lw=1.5,
             label="per-CpG mutability (CpG mut / CpG site)")
    ax3.plot(c, prof["ncpg_ratio"], "-", color="#999", lw=1.2,
             label="non-CpG C>T per site (comparison)")
    ax3.set_xlim(0, T2T_LEN); ax3.set_ylabel("Mutated haplotype fraction per site")
    ax3.set_xlabel("Position in 2168 bp repeat unit (bp)")
    shade(ax3)
    ax3.set_title("C   Ratio = per-CpG mutability along the copy",
                  fontweight="bold", loc="left")
    ax3.legend(fontsize=8, loc="upper right")

    fig.suptitle(f"CpG frequency, CpG-mutation frequency, and per-CpG mutability "
                 f"along the 5S unit ({n_haps} haplotypes)",
                 fontsize=13, fontweight="bold")
    out = FIG_DIR / "39b_cpg_along_copy.pdf"
    fig.savefig(out, bbox_inches="tight", dpi=150); plt.close(fig)
    print(f"Saved {out.name}")
    # region summary
    print(f"\nMean per-CpG mutability by region:")
    for rg, (lo, hi) in REGIONS.items():
        mask = (c >= lo) & (c <= hi)
        print(f"  {REGION_LABELS[rg]:9s}: {np.nanmean(prof['ratio'][mask]):.3f} "
              f"(non-CpG {np.nanmean(prof['ncpg_ratio'][mask]):.3f})")


def main():
    v, cons, n_haps = load()
    print(f"Loaded {len(v):,} (pos,ref,alt,region) variant groups, {n_haps} haplotypes")
    rates, bcs = normalized_spectrum(v, cons, n_haps)
    fig_spectrum(rates, bcs, n_haps)
    prof = cpg_profile(v, cons, n_haps)
    print(f"\nGenome-wide CpG C>T enrichment: "
          f"{prof['enrich']:.2f}x  (CpG {prof['rate_cpg']*1e3:.3f} vs "
          f"non-CpG {prof['rate_ncpg']*1e3:.3f} per cytosine per hap)")
    fig_cpg(prof, n_haps)


if __name__ == "__main__":
    main()
