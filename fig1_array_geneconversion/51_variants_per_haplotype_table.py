#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 51_variants_per_haplotype_table.py — Builds a per-haplotype table of interior
# copy count and unique SNV/indel counts split by 5S gene vs non-transcribed spacer.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
51_variants_per_haplotype_table.py

Per-haplotype variant summary table.

One row per haplotype: interior copy count and number of unique SNVs
(distinct position+alt allele within that haplotype, consensus_t2t pool),
split into 5S gene vs non-transcribed-spacer.

Outputs: tables/q_variants_per_haplotype.tsv + .xlsx
"""

import os
import sqlite3
from pathlib import Path
import pandas as pd

DB = Path(os.environ.get("FIVES_DB", "5S_rDNA.db"))
TABLES = Path(os.environ.get("FIVES_OUT", "output")) / "tables"
GENE = (630, 749)
COHORTS = ("HPRC_Year1", "HPRC_Release2")


def main():
    con = sqlite3.connect(DB)
    # metadata + interior copy count per haplotype
    meta = pd.read_sql_query(f"""
        SELECT a.sample_id, a.cohort, a.superpopulation, a.population,
               h.hap_label, h.haplotype_id,
               (SELECT COUNT(*) FROM copy c
                WHERE c.haplotype_id = h.haplotype_id AND c.border_note='interior')
               AS n_interior_copies
        FROM haplotype h JOIN assembly a USING(assembly_id)
        WHERE a.cohort IN {COHORTS}
    """, con)
    # distinct (pos,alt) SNVs per haplotype (unmasked)
    var = pd.read_sql_query(f"""
        SELECT DISTINCT h.haplotype_id, v.consensus_pos AS pos, v.alt
        FROM variant v JOIN copy c USING(copy_id)
        JOIN haplotype h USING(haplotype_id) JOIN assembly a USING(assembly_id)
        WHERE c.border_note='interior' AND a.cohort IN {COHORTS}
          AND v.alignment_source='consensus_t2t' AND v.var_type='snp'
          AND length(v.ref)=1 AND length(v.alt)=1
          AND v.ref IN ('A','C','G','T') AND v.alt IN ('A','C','G','T')
          AND v.masked=0
    """, con)
    BOUNDARY_EXCL = set(range(1, 90)) | set(range(2034, 2169))

    # indels per haplotype — count per (pos, var_type)
    # Separates copy-variable (n_carrying < n_interior) from haplotype-fixed (all copies)
    indel_raw = pd.read_sql_query(f"""
        SELECT h.haplotype_id, v.consensus_pos AS pos, v.var_type, v.region,
               COUNT(DISTINCT c.copy_id) AS n_carrying
        FROM variant v JOIN copy c USING(copy_id)
        JOIN haplotype h USING(haplotype_id) JOIN assembly a USING(assembly_id)
        WHERE c.border_note='interior' AND a.cohort IN {COHORTS}
          AND v.alignment_source='consensus_t2t' AND v.var_type IN ('ins','del')
        GROUP BY h.haplotype_id, v.consensus_pos, v.var_type
    """, con)
    con.close()

    var["is_gene"] = var["pos"].between(*GENE)
    g = var.groupby("haplotype_id")
    counts = pd.DataFrame({
        "n_unique_snvs": g.size(),
        "n_unique_gene_snvs": g["is_gene"].sum(),
    }).reset_index()
    counts["n_unique_nts_snvs"] = counts["n_unique_snvs"] - counts["n_unique_gene_snvs"]

    n_int = meta[["haplotype_id", "n_interior_copies"]].drop_duplicates()
    indel_raw = indel_raw.merge(n_int, on="haplotype_id")
    indel_raw = indel_raw[~indel_raw["pos"].isin(BOUNDARY_EXCL)]

    # copy-variable: present in ≥1 but <all copies (includes singletons/doubletons)
    indel = indel_raw[indel_raw["n_carrying"] < indel_raw["n_interior_copies"]].copy()
    # haplotype-fixed: present in every interior copy
    fixed = indel_raw[indel_raw["n_carrying"] == indel_raw["n_interior_copies"]].copy()
    indel["is_gene"] = indel["pos"].between(*GENE)
    gi = indel.groupby("haplotype_id")
    indel_counts = pd.DataFrame({
        "n_unique_indels": gi.size(),
        "n_unique_gene_indels": gi["is_gene"].sum(),
    }).reset_index()
    indel_counts["n_unique_nts_indels"] = (
        indel_counts["n_unique_indels"] - indel_counts["n_unique_gene_indels"])

    fixed["is_gene"] = fixed["pos"].between(*GENE)
    gf = fixed.groupby("haplotype_id")
    fixed_counts = pd.DataFrame({
        "n_hap_fixed_indels":      gf.size(),
        "n_hap_fixed_gene_indels": gf["is_gene"].sum(),
    }).reset_index()
    fixed_counts["n_hap_fixed_nts_indels"] = (
        fixed_counts["n_hap_fixed_indels"] - fixed_counts["n_hap_fixed_gene_indels"])

    df = meta.merge(counts, on="haplotype_id", how="left")
    df = df.merge(indel_counts, on="haplotype_id", how="left")
    df = df.merge(fixed_counts, on="haplotype_id", how="left")
    fill_zero = ["n_unique_snvs", "n_unique_gene_snvs", "n_unique_nts_snvs",
                 "n_unique_indels", "n_unique_gene_indels", "n_unique_nts_indels",
                 "n_hap_fixed_indels", "n_hap_fixed_gene_indels", "n_hap_fixed_nts_indels"]
    df = df.fillna({c: 0 for c in fill_zero})
    for c in fill_zero:
        df[c] = df[c].astype(int)
    df["n_unique_variants"] = df["n_unique_snvs"] + df["n_unique_indels"]
    df["variants_per_interior_copy"] = (df["n_unique_variants"] /
                                        df["n_interior_copies"].replace(0, pd.NA)).round(3)
    df["snvs_per_interior_copy"] = (df["n_unique_snvs"] /
                                    df["n_interior_copies"].replace(0, pd.NA)).round(3)

    cols = ["sample_id", "cohort", "superpopulation", "population", "hap_label",
            "n_interior_copies",
            "n_unique_snvs", "n_unique_gene_snvs", "n_unique_nts_snvs",
            "n_unique_indels", "n_unique_gene_indels", "n_unique_nts_indels",
            "n_hap_fixed_indels", "n_hap_fixed_gene_indels", "n_hap_fixed_nts_indels",
            "n_unique_variants", "variants_per_interior_copy", "snvs_per_interior_copy"]
    df = df[cols].sort_values(["sample_id", "hap_label"]).reset_index(drop=True)

    TABLES.mkdir(exist_ok=True)
    df.to_csv(TABLES/"q_variants_per_haplotype.tsv", sep="\t", index=False)
    df.to_excel(TABLES/"q_variants_per_haplotype.xlsx", index=False)
    print(f"{len(df)} haplotypes from {df['sample_id'].nunique()} individuals")
    print(df.head(3).to_string(index=False))
    print(f"\nmean SNVs: {df['n_unique_snvs'].mean():.1f}  "
          f"copy-variable indels: {df['n_unique_indels'].mean():.1f}  "
          f"hap-fixed indels: {df['n_hap_fixed_indels'].mean():.1f}")
    print(f"Wrote tables/q_variants_per_haplotype.tsv + .xlsx")

    # ── Separate table: haplotype-fixed indel positions ───────────────────────
    if not fixed.empty:
        fixed_detail = fixed.merge(
            meta[["haplotype_id","sample_id","cohort","superpopulation","population","hap_label"]],
            on="haplotype_id")
        fixed_detail = fixed_detail[
            ["sample_id","cohort","superpopulation","population","hap_label",
             "pos","var_type","region","n_carrying","n_interior_copies"]
        ].sort_values(["sample_id","hap_label","pos"]).reset_index(drop=True)
        fixed_detail["copy_fraction"] = (
            fixed_detail["n_carrying"] / fixed_detail["n_interior_copies"]).round(3)
        fixed_detail.to_csv(TABLES/"q_haplotype_fixed_indels.tsv", sep="\t", index=False)
        print(f"\nWrote q_haplotype_fixed_indels.tsv: {len(fixed_detail)} rows "
              f"across {fixed_detail['sample_id'].nunique()} samples")

    # ── Separate table: copy-variable indel positions ─────────────────────────
    if not indel.empty:
        indel_detail = indel.merge(
            meta[["haplotype_id","sample_id","cohort","superpopulation","population","hap_label"]],
            on="haplotype_id")
        indel_detail = indel_detail[
            ["sample_id","cohort","superpopulation","population","hap_label",
             "pos","var_type","region","n_carrying","n_interior_copies"]
        ].sort_values(["sample_id","hap_label","pos"]).reset_index(drop=True)
        indel_detail["copy_fraction"] = (
            indel_detail["n_carrying"] / indel_detail["n_interior_copies"]).round(3)
        indel_detail.to_csv(TABLES/"q_copy_variable_indels.tsv", sep="\t", index=False)
        print(f"\nWrote q_copy_variable_indels.tsv: {len(indel_detail)} rows "
              f"across {indel_detail['sample_id'].nunique()} samples")


if __name__ == "__main__":
    main()
