#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 52_variant_catalog_table.py — Builds a pan-population catalog of unique 5S
# variants with per-variant copy and donor counts and frequencies.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
52_variant_catalog_table.py

Pan-population variant catalog: one row per unique variant position/type,
with copy and donor counts across all HPRC_Year1 + HPRC_Release2 haplotypes.

SNVs:    distinct (pos, alt), masked=0, interior copies only.
Indels:  distinct (pos, var_type), copy-variable only — positions where ALL
         copies of a haplotype carry the indel (haplotype-fixed) are excluded.

Output: tables/q_variant_catalog.tsv
"""

import os
import sqlite3
from pathlib import Path
import pandas as pd

DB     = Path(os.environ.get("FIVES_DB", "5S_rDNA.db"))
TABLES = Path(os.environ.get("FIVES_OUT", "output")) / "tables"
GENE   = (630, 749)
COHORTS = ("HPRC_Year1", "HPRC_Release2")
BOUNDARY_EXCL = set(range(1, 90)) | set(range(2034, 2169))


def main():
    con = sqlite3.connect(DB)

    # ── total interior copies (denominator) ───────────────────────────────────
    total_copies = pd.read_sql_query(f"""
        SELECT COUNT(*) AS n FROM copy c
        JOIN haplotype h USING(haplotype_id)
        JOIN assembly a USING(assembly_id)
        WHERE c.border_note='interior' AND a.cohort IN {COHORTS}
    """, con).iloc[0, 0]

    total_donors = pd.read_sql_query(f"""
        SELECT COUNT(DISTINCT a.sample_id) AS n FROM assembly a
        WHERE a.cohort IN {COHORTS}
    """, con).iloc[0, 0]

    # ── SNVs: per-copy occurrences ────────────────────────────────────────────
    snv = pd.read_sql_query(f"""
        SELECT a.sample_id, v.consensus_pos AS pos, v.ref, v.alt,
               v.region, COUNT(DISTINCT c.copy_id) AS n_copies_hap
        FROM variant v JOIN copy c USING(copy_id)
        JOIN haplotype h USING(haplotype_id) JOIN assembly a USING(assembly_id)
        WHERE c.border_note='interior' AND a.cohort IN {COHORTS}
          AND v.alignment_source='consensus_t2t' AND v.var_type='snp'
          AND length(v.ref)=1 AND length(v.alt)=1
          AND v.ref IN ('A','C','G','T') AND v.alt IN ('A','C','G','T')
          AND v.masked=0
        GROUP BY a.sample_id, v.consensus_pos, v.ref, v.alt, v.region
    """, con)

    snv = snv[~snv["pos"].isin(BOUNDARY_EXCL)]
    snv["var_type"] = "snp"
    snv_cat = (snv.groupby(["region", "pos", "ref", "alt", "var_type"])
               .agg(n_copies=("n_copies_hap", "sum"),
                    n_donors=("sample_id", "nunique"))
               .reset_index())

    # ── Indels: copy-variable only (exclude haplotype-fixed) ─────────────────
    # For each (haplotype, pos, var_type), get n_carrying and n_interior_copies.
    # Keep only rows where n_carrying < n_interior_copies.
    indel_raw = pd.read_sql_query(f"""
        SELECT a.sample_id, h.haplotype_id, v.consensus_pos AS pos, v.var_type,
               v.ref, v.region,
               COUNT(DISTINCT c.copy_id) AS n_carrying,
               (SELECT COUNT(*) FROM copy c2
                WHERE c2.haplotype_id = h.haplotype_id
                  AND c2.border_note = 'interior') AS n_interior
        FROM variant v JOIN copy c USING(copy_id)
        JOIN haplotype h USING(haplotype_id) JOIN assembly a USING(assembly_id)
        WHERE c.border_note='interior' AND a.cohort IN {COHORTS}
          AND v.alignment_source='consensus_t2t' AND v.var_type IN ('ins','del')
        GROUP BY a.sample_id, h.haplotype_id, v.consensus_pos, v.var_type, v.ref, v.region
    """, con)
    con.close()

    indel_raw = indel_raw[~indel_raw["pos"].isin(BOUNDARY_EXCL)]
    # copy-variable only
    indel_var = indel_raw[indel_raw["n_carrying"] < indel_raw["n_interior"]].copy()

    # alt = "-" for indels (symbolic, not stored per-base in catalog)
    indel_var["alt"] = "-"
    indel_cat = (indel_var.groupby(["region", "pos", "ref", "alt", "var_type"])
                 .agg(n_copies=("n_carrying", "sum"),
                      n_donors=("sample_id", "nunique"))
                 .reset_index())

    # ── combine ───────────────────────────────────────────────────────────────
    cat = pd.concat([snv_cat, indel_cat], ignore_index=True)
    cat["is_gene_variant"] = cat["pos"].between(*GENE).astype(int)
    cat["global_copy_freq"] = cat["n_copies"] / total_copies
    cat["donor_freq"]       = cat["n_donors"] / total_donors
    cat = cat.sort_values("n_copies", ascending=False).reset_index(drop=True)
    cols = ["region", "pos", "ref", "alt", "is_gene_variant", "var_type",
            "n_copies", "n_donors", "global_copy_freq", "donor_freq"]
    cat = cat[cols]

    TABLES.mkdir(exist_ok=True)
    cat.to_csv(TABLES / "q_variant_catalog.tsv", sep="\t", index=False)

    snvs  = (cat["var_type"] == "snp").sum()
    indels = (cat["var_type"] != "snp").sum()
    print(f"{len(cat)} variants ({snvs} SNVs, {indels} indels) across {total_donors} donors")
    print(f"Total interior copies (denominator): {total_copies:,}")
    print(f"Top 5 variants:")
    print(cat.head(5).to_string(index=False))
    print(f"\nWrote tables/q_variant_catalog.tsv")


if __name__ == "__main__":
    main()
