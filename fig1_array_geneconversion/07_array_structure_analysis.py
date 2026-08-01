#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 07_array_structure_analysis.py — summarizes 5S rDNA array structure (copy number, repeat-unit length, identity, and per-position variant frequency) across HPRC haplotypes and writes the figures and summary tables.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
HPRC 5S array structure analysis — 47 samples × 2 haplotypes = 94 haplotypes.

Loads all assembly databases, merges with population metadata, and summarizes
5S rDNA array organization.

Outputs:
  figures/07_array_structure.pdf   — multi-panel figure
  figures/07_variant_heatmap.pdf   — position-level variant heatmap
  tables/07_array_summary.tsv      — per-haplotype summary table
  tables/07_variant_catalog.tsv    — pan-population variant catalog
"""

import os
import re
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
from scipy import stats

HPRC  = Path(os.environ.get("FIVES_OUT", "output"))
DB    = Path(os.environ.get("FIVES_DATA", "data")) / "databases"
BLAST = Path(os.environ.get("FIVES_DATA", "data")) / "blast"
INV   = Path(os.environ.get("FIVES_DATA", "data")) / "data_inventory.tsv"
FIG   = HPRC / "figures"
TAB   = HPRC / "tables"
FIG.mkdir(exist_ok=True)
TAB.mkdir(exist_ok=True)

# Region lengths from T2T consensus (2168 bp total)
# NTS-pre: 0–628, 5S gene: 629–747, NTS-post: 748–2167
NTS_PRE_LEN  = 629
GENE_LEN     = 119
NTS_POST_LEN = 1419
UNIT_LEN     = NTS_PRE_LEN + GENE_LEN + NTS_POST_LEN   # 2167
GENE_START   = NTS_PRE_LEN                               # 629
GENE_END     = NTS_PRE_LEN + GENE_LEN                   # 748

# gene_variants positions are 0-based in the samtools faidx -i extracted sequence.
# For minus-strand arrays: pos maps directly to consensus position.
# For plus-strand arrays: consensus_pos = UNIT_LEN - 1 - pos  (sequence is RC'd).
# Strand is read from BLAST files.

BLAST_COLS = ["qseqid","sseqid","pident","length","mismatch","gapopen",
              "qstart","qend","sstart","send","evalue","bitscore","qlen","slen","sstrand"]

SUPERPOPS = ["AFR", "AMR", "EAS", "SAS"]  # HPRC Year 1 has no EUR

# Sub-region definitions for the repeat unit.
# NTS-pre fractions (0–1 of NTS_PRE_LEN), gene as literal, NTS-post fractions (0–1 of NTS_POST_LEN).
REGIONS = [
    # name                  start_pos                                  end_pos                                     color
    ("CA repeats",          0,                                         int(0.22 * NTS_PRE_LEN),                    "#74add1"),
    ("Spacer prom.",        int(0.22 * NTS_PRE_LEN),                  int(0.99 * NTS_PRE_LEN),                    "#4575b4"),
    ("UPE",                 int(0.99 * NTS_PRE_LEN),                  NTS_PRE_LEN,                                "#313695"),
    ("5S gene",             GENE_START,                                GENE_END,                                   "#ff9900"),
    ("Terminator",          GENE_END,                                  GENE_END + int(0.03 * NTS_POST_LEN),        "#d7301f"),
    ("Alu-like",            GENE_END + int(0.03 * NTS_POST_LEN),      GENE_END + int(0.34 * NTS_POST_LEN),        "#fc8d59"),
    ("Spacer",              GENE_END + int(0.34 * NTS_POST_LEN),      GENE_END + int(0.78 * NTS_POST_LEN),        "#fee090"),
    ("CTTCAA/TC-rich",      GENE_END + int(0.78 * NTS_POST_LEN),      UNIT_LEN,                                   "#abdda4"),
]


# ── load data ─────────────────────────────────────────────────────────────────

def get_array_strands():
    """
    Read BLAST files to determine the predominant array strand for each haplotype.
    Returns dict: {(sample_id, haplotype): 'plus' or 'minus'}.

    gene_variants positions are 0-based in the samtools faidx -i extracted sequence.
    For minus-strand arrays the -i RC gives same orientation as the consensus.
    For plus-strand arrays the -i RC reverses orientation, so positions must be
    mapped as: consensus_pos = UNIT_LEN - 1 - gene_seqs_pos.
    """
    strands = {}
    for f in sorted(BLAST.glob("*_blast.txt")):
        name = f.stem.replace("_blast", "")
        parts = name.rsplit("_", 1)
        if len(parts) != 2:
            continue
        sample, hap = parts
        try:
            df = pd.read_csv(f, sep="\t", names=BLAST_COLS, usecols=["pident","length","sseqid","sstart","send","sstrand"])
        except Exception:
            continue
        hits = df[(df["length"] >= 115) & (df["pident"] >= 95)]
        if hits.empty:
            continue
        chrom = hits["sseqid"].value_counts().idxmax()
        chr1  = hits[hits["sseqid"] == chrom]
        main_strand = chr1["sstrand"].value_counts().idxmax()
        strands[(sample, hap)] = main_strand
    return strands


def load_all():
    inv = pd.read_csv(INV, sep="\t")
    pop = inv[["sample_id", "population", "superpopulation"]].drop_duplicates()

    frames = []
    for f in sorted(DB.glob("*.tsv")):
        df = pd.read_csv(f, sep="\t")
        frames.append(df)
    df = pd.concat(frames, ignore_index=True)
    df = df.merge(pop, on="sample_id", how="left")

    # Attach array strand from BLAST files
    strands = get_array_strands()
    df["array_strand"] = df.apply(
        lambda r: strands.get((r["sample_id"], r["haplotype"]), "minus"), axis=1)
    return df


def parse_variants(variant_str):
    """Return list of (pos, ref, alt) from variant string like '4:T>A; 100:C>G'."""
    if not isinstance(variant_str, str) or not variant_str.strip():
        return []
    variants = []
    for token in variant_str.split(";"):
        token = token.strip()
        m = re.match(r"(\d+):([A-Z-]+)>([A-Z-]+)", token)
        if m:
            variants.append((int(m.group(1)), m.group(2), m.group(3)))
    return variants


# ── analysis ──────────────────────────────────────────────────────────────────

def compute_hap_summary(df):
    """Per-haplotype: copy number, mean unit length, mean gene identity, SNV counts."""
    g = df.groupby(["sample_id", "haplotype", "superpopulation", "population"])
    summary = g.agg(
        n_copies        = ("copy_id",           "count"),
        mean_unit_bp    = ("unit_length_bp",     "mean"),
        median_unit_bp  = ("unit_length_bp",     "median"),
        sd_unit_bp      = ("unit_length_bp",     "std"),
        mean_identity   = ("gene_pct_identity",  "mean"),
        min_identity    = ("gene_pct_identity",  "min"),
        n_identical     = ("category",           lambda x: (x == "identical").sum()),
        mean_snv_gene   = ("n_snv_gene",         "mean"),
        mean_snv_pre    = ("n_snv_nts_pre",      "mean"),
        mean_snv_post   = ("n_snv_nts_post",     "mean"),
    ).reset_index()
    summary["diploid_copies"] = summary.groupby("sample_id")["n_copies"].transform("sum")
    return summary


def build_variant_catalog(df, region):
    """
    For a given region ('gene', 'nts_pre', 'nts_post'), collect all variants
    across all copies. Returns DataFrame with columns:
      pos, ref, alt, n_copies, n_haplotypes, n_samples, freq_copies, freq_haplotypes
    """
    col = f"{region}_variants"
    total_copies = len(df)
    total_haps   = df[["sample_id", "haplotype"]].drop_duplicates().shape[0]

    records = defaultdict(lambda: {"copies": set(), "haps": set(), "samples": set()})
    for _, row in df.iterrows():
        for pos, ref, alt in parse_variants(row[col]):
            key = (pos, ref, alt)
            records[key]["copies"].add((row["sample_id"], row["haplotype"], row["copy_id"]))
            records[key]["haps"].add((row["sample_id"], row["haplotype"]))
            records[key]["samples"].add(row["sample_id"])

    rows = []
    for (pos, ref, alt), d in records.items():
        rows.append({
            "region": region,
            "pos": pos, "ref": ref, "alt": alt,
            "n_copies":      len(d["copies"]),
            "n_haplotypes":  len(d["haps"]),
            "n_samples":     len(d["samples"]),
            "freq_copies":   len(d["copies"]) / total_copies,
            "freq_haplotypes": len(d["haps"]) / total_haps,
        })
    return pd.DataFrame(rows).sort_values("pos")


def position_frequency_vector(df, region, length):
    """Return array[pos] = fraction of copies carrying any variant at that position."""
    col = f"{region}_variants"
    total = len(df)
    counts = np.zeros(length + 1)
    for vstr in df[col].dropna():
        for pos, _, _ in parse_variants(vstr):
            if 1 <= pos <= length:
                counts[pos] += 1
    return counts[1:] / total   # 0-indexed, pos 1..length


def position_frequency_vector_gene(df):
    """
    Return array[cons_pos] = fraction of copies with a variant at that consensus position,
    using gene_variants which covers the full 2168 bp repeat unit.

    Positions in gene_variants are 0-based in the samtools faidx -i extracted sequence.
    Orientation depends on the BLAST hit strand stored in df['array_strand']:
      minus: consensus_pos = gene_seqs_pos  (same orientation as T2T consensus)
      plus:  consensus_pos = UNIT_LEN - 1 - gene_seqs_pos  (RC'd → reversed)
    """
    total  = len(df)
    counts = np.zeros(UNIT_LEN)
    for _, row in df.iterrows():
        is_plus = (row.get("array_strand", "minus") == "plus")
        for pos, _, _ in parse_variants(row["gene_variants"]):
            if is_plus:
                cpos = UNIT_LEN - 1 - pos
            else:
                cpos = pos
            if 0 <= cpos < UNIT_LEN:
                counts[cpos] += 1
    return counts / total


# ── figures ───────────────────────────────────────────────────────────────────

POP_COLORS = {
    "AFR": "#E41A1C", "AMR": "#FF7F00", "EAS": "#377EB8",
    "EUR": "#4DAF4A", "SAS": "#984EA3",
}

def fig_copy_number(summary):
    """Panel: copy number distribution by superpopulation."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # left: violin per superpop
    ax = axes[0]
    pops = [p for p in SUPERPOPS if p in summary["superpopulation"].values]
    data = [summary[summary["superpopulation"] == p]["n_copies"].values for p in pops]
    vp = ax.violinplot(data, positions=range(len(pops)), showmedians=True, showextrema=True)
    for i, (body, pop) in enumerate(zip(vp["bodies"], pops)):
        body.set_facecolor(POP_COLORS[pop])
        body.set_alpha(0.7)
    ax.set_xticks(range(len(pops)))
    ax.set_xticklabels(pops)
    ax.set_ylabel("Copies per haplotype")
    ax.set_title("Copy number by superpopulation")

    # right: histogram all haplotypes
    ax = axes[1]
    ax.hist(summary["n_copies"], bins=30, color="#555", edgecolor="white", linewidth=0.5)
    ax.axvline(summary["n_copies"].median(), color="red", linestyle="--", label=f"Median {summary['n_copies'].median():.0f}")
    ax.set_xlabel("Copies per haplotype")
    ax.set_ylabel("Haplotypes")
    ax.set_title("Copy number distribution (n=94 haplotypes)")
    ax.legend()

    fig.tight_layout()
    return fig


def fig_unit_length(df):
    """Panel: repeat unit length distribution."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    ax = axes[0]
    ax.hist(df["unit_length_bp"], bins=60, color="#2166ac", edgecolor="white", linewidth=0.3)
    ax.axvline(2168, color="red", linestyle="--", label="T2T consensus (2168 bp)")
    ax.set_xlabel("Repeat unit length (bp)")
    ax.set_ylabel("Copies")
    ax.set_title(f"Unit length distribution (n={len(df):,} copies)")
    ax.legend()

    # unit length vs copy number
    ax = axes[1]
    hap = df.groupby(["sample_id", "haplotype"]).agg(
        n_copies=("copy_id", "count"), mean_unit=("unit_length_bp", "mean")).reset_index()
    for pop, grp in hap.merge(df[["sample_id", "haplotype", "superpopulation"]].drop_duplicates(),
                               on=["sample_id", "haplotype"]).groupby("superpopulation"):
        ax.scatter(grp["n_copies"], grp["mean_unit"], c=POP_COLORS.get(pop, "gray"),
                   label=pop, alpha=0.7, s=40)
    r, p = stats.pearsonr(hap["n_copies"], hap["mean_unit"])
    ax.set_xlabel("Copies per haplotype")
    ax.set_ylabel("Mean unit length (bp)")
    ax.set_title(f"Unit length vs copy number (r={r:.2f}, p={p:.2e})")
    ax.legend(title="Superpop", fontsize=8)

    fig.tight_layout()
    return fig


def fig_identity(df):
    """
    Repeat unit pairwise identity to T2T consensus.
    Note: gene_pct_identity is the BLAST pident of the full ~2167 bp repeat unit
    alignment, not just the 119 bp gene — the column name is a legacy artefact.
    """
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.hist(df["gene_pct_identity"], bins=50, color="#d6604d", edgecolor="white", linewidth=0.3)
    ax.axvline(100, color="k", linestyle="--", alpha=0.5, label="100% (identical unit)")
    ax.axvline(99,  color="gray", linestyle=":", alpha=0.7, label="99% threshold")
    ax.set_xlabel("Repeat unit % identity to T2T consensus (full ~2167 bp unit)")
    ax.set_ylabel("Copies")
    ax.set_title(f"Repeat unit identity to T2T consensus  (n={len(df):,} copies, 94 haplotypes)")
    pct_99 = (df["gene_pct_identity"] >= 99).mean() * 100
    pct_id = (df["gene_pct_identity"] == 100).mean() * 100
    ax.text(0.02, 0.95, f"≥99%: {pct_99:.1f}%\n=100%: {pct_id:.1f}%",
            transform=ax.transAxes, va="top", fontsize=9)
    ax.legend()
    fig.tight_layout()
    return fig


def fig_variant_position_profile(df):
    """
    Two-panel figure:
      Top: per-position variant frequency (dark red bars) on coloured sub-region background.
      Bottom: repeat unit structure strip with labelled sub-region blocks.
    """
    freq = position_frequency_vector_gene(df)   # shape (UNIT_LEN,), consensus coords
    x    = np.arange(UNIT_LEN)

    fig, (ax_freq, ax_map) = plt.subplots(
        2, 1, figsize=(16, 5),
        gridspec_kw={"height_ratios": [5, 1]},
    )
    fig.subplots_adjust(hspace=0.08)

    # ── top panel: coloured backgrounds + dark-red frequency bars ─────────────
    for name, lo, hi, color in REGIONS:
        ax_freq.axvspan(lo, hi, color=color, alpha=0.25, linewidth=0)

    ax_freq.bar(x, freq * 100, width=1, color="#8b1a1a", alpha=0.85)

    # Boundary lines at each region edge
    seen = set()
    for _, lo, hi, _ in REGIONS:
        for xpos in (lo, hi):
            if xpos not in seen and 0 < xpos < UNIT_LEN:
                ax_freq.axvline(xpos, color="gray", linewidth=0.5, linestyle=":")
                seen.add(xpos)

    # Sub-region name labels at top of panel
    ymax = freq.max() * 100
    ax_freq.set_ylim(0, ymax * 1.25)
    for name, lo, hi, color in REGIONS:
        mid = (lo + hi) / 2
        ax_freq.text(mid, ax_freq.get_ylim()[1] * 0.97, name,
                     ha="center", va="top", fontsize=7.5,
                     color=color if name != "Spacer" else "#a08020",
                     fontweight="bold" if name == "5S gene" else "normal")

    ax_freq.set_xlim(0, UNIT_LEN)
    ax_freq.set_ylabel("% copies\nwith variant", fontsize=9)
    ax_freq.set_xticklabels([])
    n_plus  = (df["array_strand"] == "plus").sum()
    n_minus = (df["array_strand"] == "minus").sum()
    ax_freq.set_title(
        f"Variant frequency across full 5S rDNA repeat unit  "
        f"(n={len(df):,} copies, 94 haplotypes, 47 individuals)",
        fontsize=11,
    )

    # ── bottom panel: repeat unit structure strip ─────────────────────────────
    ax_map.set_xlim(0, UNIT_LEN)
    ax_map.set_ylim(0, 1)
    ax_map.axis("off")

    for name, lo, hi, color in REGIONS:
        rect = plt.Rectangle((lo, 0.1), hi - lo, 0.8,
                              facecolor=color, edgecolor="white", linewidth=0.5)
        ax_map.add_patch(rect)
        mid = (lo + hi) / 2
        bp  = hi - lo
        # Skip label if region too narrow to fit text
        if bp > 30:
            ax_map.text(mid, 0.5, f"{bp}", ha="center", va="center",
                        fontsize=7, color="white" if name in ("Spacer prom.", "UPE", "Terminator") else "black")
        ax_map.text(mid, -0.05, name, ha="center", va="top", fontsize=6.5, color="#333")

    ax_map.set_xlabel("Position within repeat unit (bp)", fontsize=9,
                      labelpad=14)
    # Restore x-ticks via invisible axis
    ax_twin = ax_map.twiny()
    ax_twin.set_xlim(0, UNIT_LEN)
    ax_twin.xaxis.set_ticks_position("bottom")
    ax_twin.xaxis.set_label_position("bottom")
    ax_twin.tick_params(axis="x", labelsize=8, pad=2)

    return fig


def fig_snv_burden(df):
    """Scatter: per-copy SNV burden in gene vs NTS, coloured by identity."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    ax = axes[0]
    sc = ax.scatter(df["n_snv_nts_pre"] + df["n_snv_nts_post"], df["n_snv_gene"],
                    c=df["gene_pct_identity"], cmap="RdYlGn", alpha=0.3, s=8, vmin=95, vmax=100)
    plt.colorbar(sc, ax=ax, label="Gene % identity")
    ax.set_xlabel("NTS SNVs (pre + post)")
    ax.set_ylabel("Gene SNVs")
    ax.set_title("Per-copy SNV burden: gene vs NTS")
    r, p = stats.pearsonr(df["n_snv_nts_pre"] + df["n_snv_nts_post"], df["n_snv_gene"])
    ax.text(0.05, 0.95, f"r={r:.2f}", transform=ax.transAxes, va="top")

    # right: SNV burden by border position
    ax = axes[1]
    has_border = df["border_note"].isin(["5-prime_array_border", "3-prime_array_border"])
    ax.boxplot(
        [df[~has_border]["n_snv_gene"].values, df[has_border]["n_snv_gene"].values],
        labels=["Interior copies", "Border copies"],
        notch=False, patch_artist=True,
        boxprops=dict(facecolor="#aec7e8")
    )
    ax.set_ylabel("Gene SNVs per copy")
    ax.set_title("SNV burden: interior vs border copies")
    t, p2 = stats.mannwhitneyu(
        df[~has_border]["n_snv_gene"], df[has_border]["n_snv_gene"], alternative="two-sided")
    ax.text(0.5, 0.95, f"MWU p={p2:.2e}", transform=ax.transAxes, ha="center", va="top")

    fig.tight_layout()
    return fig


def fig_top_variants(catalog_gene, catalog_pre, catalog_post, top_n=20):
    """Bar chart of most frequent variants in each region."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    for ax, cat, title, color in [
        (axes[0], catalog_gene, "Gene",     "#d6604d"),
        (axes[1], catalog_pre,  "NTS-pre",  "#2166ac"),
        (axes[2], catalog_post, "NTS-post", "#1a9850"),
    ]:
        top = cat.nlargest(top_n, "freq_copies")
        labels = top.apply(lambda r: f"{r['pos']}:{r['ref']}>{r['alt']}", axis=1)
        ax.barh(range(len(top)), top["freq_copies"] * 100, color=color, alpha=0.8)
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels(labels, fontsize=7)
        ax.invert_yaxis()
        ax.set_xlabel("% copies")
        ax.set_title(f"Top {top_n} variants — {title}")
    fig.suptitle("Most frequent variants across 94 haplotypes", y=1.02)
    fig.tight_layout()
    return fig


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    print("Loading data...")
    df = load_all()
    print(f"  {len(df):,} copies across {df['sample_id'].nunique()} samples, "
          f"{df[['sample_id','haplotype']].drop_duplicates().shape[0]} haplotypes")

    # ── summary table ─────────────────────────────────────────────────────────
    summary = compute_hap_summary(df)
    summary.to_csv(TAB / "07_array_summary.tsv", sep="\t", index=False)
    print(f"  Haplotype summary → {TAB / '07_array_summary.tsv'}")

    # ── variant catalogs ──────────────────────────────────────────────────────
    print("Building variant catalogs...")
    cat_gene = build_variant_catalog(df, "gene")
    cat_pre  = build_variant_catalog(df, "nts_pre")
    cat_post = build_variant_catalog(df, "nts_post")
    catalog  = pd.concat([cat_gene, cat_pre, cat_post], ignore_index=True)
    catalog.to_csv(TAB / "07_variant_catalog.tsv", sep="\t", index=False)
    print(f"  Variant catalog → {TAB / '07_variant_catalog.tsv'}")
    print(f"    gene: {len(cat_gene)} unique variants")
    print(f"    NTS-pre: {len(cat_pre)} unique variants")
    print(f"    NTS-post: {len(cat_post)} unique variants")

    # ── key statistics ────────────────────────────────────────────────────────
    print("\n── Array structure summary ──")
    print(f"  Copy number range:  {summary['n_copies'].min()}–{summary['n_copies'].max()} per haplotype")
    print(f"  Diploid copy range: {summary['diploid_copies'].min()}–{summary['diploid_copies'].max()}")
    print(f"  Median haploid CN:  {summary['n_copies'].median():.0f}")
    print(f"  Unit length:        {df['unit_length_bp'].min()}–{df['unit_length_bp'].max()} bp "
          f"(median {df['unit_length_bp'].median():.0f})")
    print(f"  Gene identity:      {df['gene_pct_identity'].min():.2f}–100% "
          f"(mean {df['gene_pct_identity'].mean():.3f}%)")
    pct_id = (df["gene_pct_identity"] == 100).mean() * 100
    pct_hs = (df["gene_pct_identity"] >= 99).mean() * 100
    print(f"  Identical copies:   {pct_id:.1f}%")
    print(f"  ≥99% identity:      {pct_hs:.1f}%")

    # Strand breakdown
    n_plus  = (df["array_strand"] == "plus").sum()
    n_minus = (df["array_strand"] == "minus").sum()
    print(f"\n  Array strand: {n_minus:,} copies minus-strand, {n_plus:,} copies plus-strand")
    hap_strands = df.groupby(["sample_id","haplotype"])["array_strand"].first()
    print(f"  Haplotypes: {(hap_strands=='minus').sum()} minus, {(hap_strands=='plus').sum()} plus")

    # Hotspot positions from strand-normalised gene_variants (consensus coordinates)
    freq_unit = position_frequency_vector_gene(df)
    hotspots = list(np.where(freq_unit > 0.05)[0])
    pre_hot  = [p for p in hotspots if p < GENE_START]
    gene_hot = [p for p in hotspots if GENE_START <= p < GENE_END]
    post_hot = [p for p in hotspots if p >= GENE_END]
    print(f"\n  Hotspots >5% (strand-normalised consensus coords):")
    print(f"    NTS-pre (0–{GENE_START-1}): {pre_hot}")
    print(f"    Gene ({GENE_START}–{GENE_END-1}): {gene_hot}")
    print(f"    NTS-post ({GENE_END}–{UNIT_LEN-1}): {post_hot[:10]}{'...' if len(post_hot)>10 else ''}")

    top_gene = cat_gene.nlargest(5, "freq_copies")[["pos","ref","alt","freq_copies","n_samples"]]
    print(f"\n  Top 5 gene_variants entries (full-unit, strand-mixed positions):\n{top_gene.to_string(index=False)}")

    # Population CN differences
    print("\n  Copy number by superpopulation:")
    for pop, grp in summary.groupby("superpopulation"):
        print(f"    {pop}: median {grp['n_copies'].median():.0f}  "
              f"(range {grp['n_copies'].min()}–{grp['n_copies'].max()})")

    # Kruskal-Wallis test on CN across superpops
    groups = [summary[summary["superpopulation"] == p]["n_copies"].values
              for p in SUPERPOPS if p in summary["superpopulation"].values]
    h, p_kw = stats.kruskal(*groups)
    print(f"  Kruskal-Wallis CN ~ superpopulation: H={h:.2f}, p={p_kw:.3e}")

    # Border copies
    border = df["border_note"].isin(["5-prime_array_border", "3-prime_array_border"])
    print(f"\n  Border copies: {border.sum()} ({border.mean()*100:.1f}%)")
    print(f"  Mean gene SNVs — interior: {df[~border]['n_snv_gene'].mean():.2f}, "
          f"border: {df[border]['n_snv_gene'].mean():.2f}")

    # Shared variants (present in >10% of samples)
    widespread_gene = cat_gene[cat_gene["n_samples"] >= int(0.1 * df["sample_id"].nunique())]
    print(f"\n  Gene variants in ≥10% of samples: {len(widespread_gene)}")
    private_gene = cat_gene[cat_gene["n_samples"] == 1]
    print(f"  Sample-private gene variants: {len(private_gene)}")

    # ── figures ───────────────────────────────────────────────────────────────
    print("\nGenerating figures...")

    fig = fig_copy_number(summary)
    fig.savefig(FIG / "07_copy_number.pdf", bbox_inches="tight")
    plt.close(fig)

    fig = fig_unit_length(df)
    fig.savefig(FIG / "07_unit_length.pdf", bbox_inches="tight")
    plt.close(fig)

    fig = fig_identity(df)
    fig.savefig(FIG / "07_gene_identity.pdf", bbox_inches="tight")
    plt.close(fig)

    fig = fig_variant_position_profile(df)
    fig.savefig(FIG / "07_variant_profile.pdf", bbox_inches="tight")
    plt.close(fig)

    fig = fig_snv_burden(df)
    fig.savefig(FIG / "07_snv_burden.pdf", bbox_inches="tight")
    plt.close(fig)

    fig = fig_top_variants(cat_gene, cat_pre, cat_post)
    fig.savefig(FIG / "07_top_variants.pdf", bbox_inches="tight")
    plt.close(fig)

    print("Done. Figures saved to", FIG)


if __name__ == "__main__":
    main()
