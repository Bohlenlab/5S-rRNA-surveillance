#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# import_legacy_assemblies.py — import CHM13 and HG002 copies/variants/read-variants into the database.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
Import CHM13 and HG002 data into 5S_rDNA.db.

These pre-date the HPRC pipeline and use a slightly different database format:
  - gene_variants: 0-based positions within the 119 bp 5S gene only
    → consensus_pos = 630 + pos_0based   (vs HPRC: pos + 1)
  - nts_pre_variants: 0-based within NTS-pre → consensus_pos = pos + 1  (same)
  - nts_post_variants: 0-based within NTS-post → consensus_pos = 749 + pos  (same)
  - No n_snv_5s_gene column (gene_variants IS the gene; n_snv_5s_gene = n_snv_gene)

Assemblies added:
  CHM13      — haploid, one haplotype (databases/5S_array_database.tsv)
  HG002_GIAB — diploid GIAB assembly: MATERNAL=hap1, PATERNAL=hap2
               (distinct from "HG002" in DB which is the HPRC Year 1 assembly)

Read variants added (read_variant table):
  HG002_GIAB illumina — shortread_pipeline/tsv/BGIseq_150bp_t2tcons.tsv (bcftools query)

Note: CHM13_mRNAseq_150bp_t2tcons.tsv is in raw samtools pileup format (not bcftools
query format) so it is not imported here; it requires a separate parser.

Usage:
  python3 import_legacy_assemblies.py
"""

import os
import re
import sqlite3
from pathlib import Path

import pandas as pd

T2T  = Path(os.environ.get("FIVES_DATA", "data"))
DB   = Path(os.environ.get("FIVES_DB", "5S_rDNA.db"))

GENE_START_1     = 630   # 1-based inclusive
NTS_POST_START_1 = 749   # 1-based start of NTS-post

# Old-format database files
CHM13_DB  = T2T / "databases/5S_array_database.tsv"
HG002_MAT = T2T / "databases/5S_array_database_HG002_MATERNAL.tsv"
HG002_PAT = T2T / "databases/5S_array_database_HG002_PATERNAL.tsv"

# Read-variant source files
HG002_SR_TSV = T2T / "HG002/shortread_pipeline/tsv/BGIseq_150bp_t2tcons.tsv"
CHM13_RNA_TSV = T2T / "CHM13/shortread_pipeline/tsv/CHM13_mRNAseq_150bp_t2tcons.tsv"


# ── helpers ───────────────────────────────────────────────────────────────────

def region_of(pos1: int) -> str:
    if pos1 < GENE_START_1:
        return "nts_pre"
    if pos1 <= 748:
        return "gene"
    return "nts_post"


def parse_variants_str(s: str) -> list:
    """Parse "pos:ref>alt; ..." → [(0-based pos, ref, alt), ...]."""
    if not s or str(s).strip().lower() in ("none", ".", ""):
        return []
    out = []
    for token in str(s).split(";"):
        token = token.strip()
        m = re.match(r"(\d+):([A-Za-z]+)>([A-Za-z]+)", token)
        if m:
            out.append((int(m.group(1)), m.group(2), m.group(3)))
    return out


def load_mpileup_tsv(path: Path) -> pd.DataFrame:
    """Load raw bcftools query TSV (no header): pos, ref, alt, dp, ad_str."""
    rows = []
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        pos, ref, alts_str, dp, ad_str = line.split("\t")
        pos, dp = int(pos), int(dp)
        ads = list(map(int, ad_str.split(",")))
        for i, alt in enumerate(alts_str.split(",")):
            if alt in (".", "<*>"):
                continue
            ad  = ads[i + 1] if i + 1 < len(ads) else 0
            vaf = ad / dp if dp > 0 else 0.0
            rows.append({"pos": pos, "ref": ref, "alt": alt,
                         "dp": dp, "ad": ad, "vaf": vaf})
    return pd.DataFrame(rows) if rows else pd.DataFrame(
        columns=["pos", "ref", "alt", "dp", "ad", "vaf"])


def safe_int(v):
    try:
        f = float(v)
        return None if pd.isna(f) else int(f)
    except (TypeError, ValueError):
        return None


def safe_float(v):
    try:
        f = float(v)
        return None if pd.isna(f) else f
    except (TypeError, ValueError):
        return None


def safe_str(v, null_vals=("nan", "none", "", ".")):
    s = str(v).strip() if v is not None else ""
    return None if s.lower() in null_vals else s


# ── import one old-format database ───────────────────────────────────────────

def import_old_db(con, df, assembly_id, haplotype_id, n_cop_ref):
    """
    Import copies and variants from an old-format assembly database.

    Old format differences vs HPRC:
      - gene_variants: 0-based pos in 119 bp gene → consensus_pos = 630 + pos
      - nts_pre/post: same as HPRC
      - No n_snv_5s_gene: use n_snv_gene (they are equivalent in old format)
    """
    n_cop = n_var = 0

    for _, r in df.iterrows():
        copy_num = int(r["copy_id"])

        # Local coordinates: present in HG002, absent in CHM13 (genomic only)
        usl = safe_int(r.get("unit_start_local") or r.get("unit_start_genomic"))
        uel = safe_int(r.get("unit_end_local")   or r.get("unit_end_genomic"))
        glo = safe_int(r.get("gene_lo_local")    or r.get("gene_lo_genomic"))
        ghi = safe_int(r.get("gene_hi_local")    or r.get("gene_hi_genomic"))

        cur = con.execute("""
            INSERT OR IGNORE INTO copy
              (haplotype_id, copy_number,
               unit_start_local, unit_end_local, unit_length_bp,
               spacing_to_next_bp,
               gene_lo_local, gene_hi_local,
               gene_pct_identity, gene_mismatches, gene_gaps,
               n_snv_gene, n_snv_5s_gene, n_snv_nts_pre, n_snv_nts_post,
               category, border_note)
            VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
        """, (
            haplotype_id, copy_num,
            usl, uel,
            safe_int(r.get("unit_length_bp")),
            safe_float(r.get("spacing_to_next_bp")),
            glo, ghi,
            safe_float(r.get("gene_pct_identity")),
            safe_int(r.get("gene_mismatches")),
            safe_int(r.get("gene_gaps")),
            safe_int(r.get("n_snv_gene")),
            safe_int(r.get("n_snv_gene")),     # n_snv_5s_gene = n_snv_gene in old format
            safe_int(r.get("n_snv_nts_pre")),
            safe_int(r.get("n_snv_nts_post")),
            safe_str(r.get("category")),
            safe_str(str(r.get("border_note", "")).split("|")[0]),  # CHM13 has extra notes
        ))
        copy_db_id = cur.lastrowid
        n_cop += 1

        variant_rows = []

        # gene_unit variants: old format = 0-based pos in 119 bp gene
        # consensus_pos = GENE_START_1 + pos_0based = 630 + pos
        for pos0, ref, alt in parse_variants_str(str(r.get("gene_variants", ""))):
            cp1 = GENE_START_1 + pos0
            variant_rows.append((copy_db_id, "gene_unit", cp1, ref, alt, "gene"))

        # nts_pre_aln: 0-based in NTS-pre → consensus_pos = pos + 1
        for pos0, ref, alt in parse_variants_str(str(r.get("nts_pre_variants", ""))):
            cp1 = pos0 + 1
            variant_rows.append((copy_db_id, "nts_pre_aln", cp1, ref, alt, "nts_pre"))

        # nts_post_aln: 0-based in NTS-post → consensus_pos = 749 + pos
        for pos0, ref, alt in parse_variants_str(str(r.get("nts_post_variants", ""))):
            cp1 = NTS_POST_START_1 + pos0
            variant_rows.append((copy_db_id, "nts_post_aln", cp1, ref, alt, "nts_post"))

        if variant_rows:
            con.executemany("""
                INSERT INTO variant
                  (copy_id, alignment_source, consensus_pos, ref, alt, region)
                VALUES (?,?,?,?,?,?)
            """, variant_rows)
            n_var += len(variant_rows)

    return n_cop, n_var


# ── import read variants ──────────────────────────────────────────────────────

def import_read_variants(con, assembly_id, modality, rows_iter):
    """Insert read_variant rows; returns count inserted."""
    n = 0
    for r in rows_iter:
        con.execute("""
            INSERT INTO read_variant
              (assembly_id, modality, consensus_pos, ref, alt,
               depth, alt_depth, vaf, region)
            VALUES (?,?,?,?,?,?,?,?,?)
        """, r)
        n += 1
    return n


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    if not DB.exists():
        print(f"ERROR: {DB} not found — run build_database.py first")
        return

    con = sqlite3.connect(DB)
    con.execute("PRAGMA foreign_keys = ON")

    print(f"Adding legacy assemblies to {DB.name}\n")

    # ── CHM13 ─────────────────────────────────────────────────────────────────
    print("── CHM13 ──")
    con.execute("""
        INSERT OR IGNORE INTO assembly
          (sample_id, cohort, sex, has_hifi, has_illumina, has_rnaseq)
        VALUES ('CHM13','CHM13','F',0,0,1)
    """)
    con.commit()
    chm13_asm = con.execute(
        "SELECT assembly_id FROM assembly WHERE sample_id='CHM13'"
    ).fetchone()[0]

    df_chm13 = pd.read_csv(CHM13_DB, sep="\t")
    n_copies_chm13 = len(df_chm13)

    all_coords = (df_chm13["unit_start_genomic"].dropna().tolist() +
                  df_chm13["unit_end_genomic"].dropna().tolist())

    cur = con.execute("""
        INSERT OR IGNORE INTO haplotype
          (assembly_id, array_id, hap_label, array_chrom, strand,
           n_copies, array_start_local, array_end_local)
        VALUES (?,1,'hap1','chr1','minus',?,?,?)
    """, (chm13_asm, n_copies_chm13,
          int(min(all_coords)) if all_coords else None,
          int(max(all_coords)) if all_coords else None))
    chm13_hap = cur.lastrowid

    n_cop, n_var = import_old_db(con, df_chm13, chm13_asm, chm13_hap, n_copies_chm13)
    con.commit()
    print(f"  haplotype: 1   copies: {n_cop}   variants: {n_var:,}")

    # CHM13 RNA-seq: TSV is raw samtools pileup format, not bcftools query format.
    # Skipped — requires a separate parser.
    print(f"  read_variant (rnaseq): skipped (non-standard TSV format)")

    # ── HG002 ─────────────────────────────────────────────────────────────────
    print("\n── HG002_GIAB ──")
    # Use "HG002_GIAB" to distinguish from "HG002" (HPRC Year 1 assembly already in DB)
    con.execute("""
        INSERT OR IGNORE INTO assembly
          (sample_id, cohort, sex, has_hifi, has_illumina)
        VALUES ('HG002_GIAB','HG002_GIAB','M',1,1)
    """)
    con.commit()
    hg002_asm = con.execute(
        "SELECT assembly_id FROM assembly WHERE sample_id='HG002_GIAB'"
    ).fetchone()[0]

    for hap_label, db_path, array_chrom in [
        ("hap1", HG002_MAT, "chr1_MATERNAL"),
        ("hap2", HG002_PAT, "chr1_PATERNAL"),
    ]:
        df = pd.read_csv(db_path, sep="\t")
        n_c = len(df)
        all_c = (df["unit_start_local"].dropna().tolist() +
                 df["unit_end_local"].dropna().tolist())
        cur = con.execute("""
            INSERT OR IGNORE INTO haplotype
              (assembly_id, array_id, hap_label, array_chrom, strand,
               n_copies, array_start_local, array_end_local)
            VALUES (?,1,?,?,'minus',?,?,?)
        """, (hg002_asm, hap_label, array_chrom, n_c,
              int(min(all_c)) if all_c else None,
              int(max(all_c)) if all_c else None))
        hap_id = cur.lastrowid
        n_cop, n_var = import_old_db(con, df, hg002_asm, hap_id, n_c)
        con.commit()
        print(f"  {hap_label} ({array_chrom}): copies={n_cop}  variants={n_var:,}")

    # HG002_GIAB illumina read variants (BGIseq 150bp, T2T consensus aligned)
    if HG002_SR_TSV.exists():
        sr_df = load_mpileup_tsv(HG002_SR_TSV)
        rows = []
        for _, r in sr_df.iterrows():
            pos = int(r["pos"])
            rows.append((hg002_asm, "illumina", pos,
                         str(r["ref"]), str(r["alt"]),
                         safe_int(r["dp"]), safe_int(r["ad"]), safe_float(r["vaf"]),
                         region_of(pos)))
        n_rv = import_read_variants(con, hg002_asm, "illumina", rows)
        con.commit()
        print(f"  read_variant (illumina/BGIseq_t2tcons): {n_rv} rows")

    # ── final summary ─────────────────────────────────────────────────────────
    print("\n── Database summary ─────────────────────────────────────────────")
    for tbl in ["assembly", "haplotype", "copy", "variant", "read_variant"]:
        n = con.execute(f"SELECT COUNT(*) FROM [{tbl}]").fetchone()[0]
        print(f"  {tbl:20s}: {n:>8,} rows")

    print("\n  Assemblies:")
    for row in con.execute(
        "SELECT sample_id, cohort, sex FROM assembly ORDER BY cohort, sample_id"
    ).fetchall():
        print(f"    {row[0]:<12} cohort={row[1]}  sex={row[2]}")

    con.close()
    print(f"\nDone.")


if __name__ == "__main__":
    main()
