#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# import_r2_assemblies.py — import HPRC Release 2 assembly-analysis results into the database.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
Import HPRC Release 2 assembly analysis results into 5S_rDNA.db.

Reads databases/{sample_id}_{hap}.tsv and blast/{sample_id}_{hap}_blast.txt
and adds them to the existing database.

Skips samples already in the assembly table.
Idempotent: safe to re-run after adding more completed samples.

Usage:
  python3 import_r2_assemblies.py
  python3 import_r2_assemblies.py --dry-run   # count what would be imported
"""

import os
import re
import sqlite3
import sys
from pathlib import Path

import pandas as pd

T2T  = Path(os.environ.get("FIVES_DATA", "data"))
HPRC = T2T
DB   = Path(os.environ.get("FIVES_DB", "5S_rDNA.db"))
INV  = HPRC / "release2_inventory.tsv"

GENE_START_1     = 630
NTS_POST_START_1 = 749

BLAST_COLS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
              "qstart", "qend", "sstart", "send", "evalue", "bitscore",
              "qlen", "slen", "sstrand"]


def region_of(pos1: int) -> str:
    if pos1 < 630:
        return "nts_pre"
    if pos1 <= 748:
        return "gene"
    return "nts_post"


def parse_variants_str(s: str) -> list:
    if not s or s.strip().lower() in ("none", ".", ""):
        return []
    out = []
    for token in s.split(";"):
        token = token.strip()
        if not token:
            continue
        m = re.match(r"(\d+):([A-Za-z]+)>([A-Za-z]+)", token)
        if m:
            out.append((int(m.group(1)), m.group(2), m.group(3)))
    return out


def get_strand(blast_file: Path) -> str:
    if not blast_file.exists():
        return "unknown"
    try:
        df = pd.read_csv(blast_file, sep="\t", names=BLAST_COLS, usecols=["sstrand"])
        if not df.empty:
            return str(df["sstrand"].value_counts().idxmax())
    except Exception:
        pass
    return "unknown"


def safe_int(v):
    try:
        f = float(v)
        return int(f) if not pd.isna(f) else None
    except (TypeError, ValueError):
        return None


def safe_float(v):
    try:
        f = float(v)
        return f if not pd.isna(f) else None
    except (TypeError, ValueError):
        return None


def safe_str(v, null_vals=("nan", "none", "", ".")):
    s = str(v).strip() if v is not None else ""
    return None if s.lower() in null_vals else s


def main():
    dry_run = "--dry-run" in sys.argv

    if not DB.exists():
        print(f"ERROR: database not found at {DB}")
        sys.exit(1)

    inv = pd.read_csv(INV, sep="\t")
    print(f"R2 inventory: {len(inv)} samples, {len(inv)*2} haplotypes")

    con = sqlite3.connect(DB)
    con.execute("PRAGMA foreign_keys = ON")
    con.execute("PRAGMA journal_mode = WAL")

    # Already-imported sample IDs
    existing = set(r[0] for r in con.execute("SELECT sample_id FROM assembly").fetchall())
    print(f"Already in DB: {len(existing)} assemblies")

    n_asm = n_hap = n_cop = n_var = 0
    skipped = []

    for _, inv_row in inv.iterrows():
        sid = safe_str(inv_row["sample_id"])
        if sid is None:
            continue

        if sid in existing:
            skipped.append(sid)
            continue

        # Check which haplotypes have completed TSVs
        hap_pairs = [
            (safe_str(inv_row["hap1_label"]), safe_str(inv_row["assembly_hap1_s3"])),
            (safe_str(inv_row["hap2_label"]), safe_str(inv_row["assembly_hap2_s3"])),
        ]

        available = []
        for hap_label, s3_path in hap_pairs:
            if hap_label is None:
                continue
            db_file = HPRC / "databases" / f"{sid}_{hap_label}.tsv"
            if db_file.exists() and db_file.stat().st_size > 100:
                available.append((hap_label, s3_path, db_file))

        if not available:
            continue   # nothing to import for this sample yet

        if dry_run:
            print(f"  [dry-run] would import {sid}: {[h for h, _, _ in available]}")
            n_hap += len(available)
            continue

        # Insert assembly row
        h1_s3 = safe_str(inv_row.get("assembly_hap1_s3"))
        h2_s3 = safe_str(inv_row.get("assembly_hap2_s3"))

        cur = con.execute("""
            INSERT OR IGNORE INTO assembly
              (sample_id, cohort, population, superpopulation,
               assembly_hap1_s3, assembly_hap2_s3,
               has_hifi, has_illumina, has_methylation, has_rnaseq)
            VALUES (?,?,?,?,?,?,0,0,0,0)
        """, (sid, "HPRC_Release2", None, None, h1_s3, h2_s3))
        con.commit()

        asm_row = con.execute(
            "SELECT assembly_id FROM assembly WHERE sample_id=?", (sid,)
        ).fetchone()
        if asm_row is None:
            print(f"  ERROR: could not insert assembly for {sid} — skipping")
            continue
        assembly_id = asm_row[0]
        n_asm += 1

        for hap_label, s3_path, db_file in available:
            df = pd.read_csv(db_file, sep="\t")
            if df.empty:
                print(f"  WARNING: empty TSV for {sid} {hap_label}")
                continue

            blast_file = HPRC / "blast" / f"{sid}_{hap_label}_blast.txt"
            strand = get_strand(blast_file)

            array_chrom = safe_str(df["array_chrom"].iloc[0]) or "unknown"
            n_copies = len(df)

            all_coords = (
                df["unit_start_local"].dropna().tolist() +
                df["unit_end_local"].dropna().tolist()
            )
            arr_lo = int(min(all_coords)) if all_coords else None
            arr_hi = int(max(all_coords)) if all_coords else None

            cur2 = con.execute("""
                INSERT INTO haplotype
                  (assembly_id, array_id, hap_label, array_chrom, strand,
                   n_copies, array_start_local, array_end_local)
                VALUES (?,1,?,?,?,?,?,?)
            """, (assembly_id, hap_label, array_chrom, strand, n_copies, arr_lo, arr_hi))
            haplotype_id = cur2.lastrowid
            n_hap += 1

            for _, r in df.iterrows():
                copy_num = int(r["copy_id"])

                cur3 = con.execute("""
                    INSERT INTO copy
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
                    safe_int(r.get("unit_start_local")),
                    safe_int(r.get("unit_end_local")),
                    safe_int(r.get("unit_length_bp")),
                    safe_float(r.get("spacing_to_next_bp")),
                    safe_int(r.get("gene_lo_local")),
                    safe_int(r.get("gene_hi_local")),
                    safe_float(r.get("gene_pct_identity")),
                    safe_int(r.get("gene_mismatches")),
                    safe_int(r.get("gene_gaps")),
                    safe_int(r.get("n_snv_gene")),
                    safe_int(r.get("n_snv_5s_gene")),
                    safe_int(r.get("n_snv_nts_pre")),
                    safe_int(r.get("n_snv_nts_post")),
                    safe_str(r.get("category")),
                    safe_str(r.get("border_note")),
                ))
                copy_db_id = cur3.lastrowid
                n_cop += 1

                variant_rows = []

                # gene_unit: full 2168 bp unit alignment, 0-based pos → +1 for 1-based
                for pos0, ref, alt in parse_variants_str(str(r.get("gene_variants", ""))):
                    cp1 = pos0 + 1
                    variant_rows.append((copy_db_id, "gene_unit", cp1, ref, alt, region_of(cp1)))

                # nts_pre_aln: NTS-pre only, 0-based → +1 (range 1–629)
                for pos0, ref, alt in parse_variants_str(str(r.get("nts_pre_variants", ""))):
                    cp1 = pos0 + 1
                    variant_rows.append((copy_db_id, "nts_pre_aln", cp1, ref, alt, "nts_pre"))

                # nts_post_aln: NTS-post only, 0-based → +749 (range 749–2168)
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

        con.commit()
        print(f"  {sid}: {len(available)} haplotype(s) imported")

    if dry_run:
        print(f"\n[dry-run] would import {n_hap} haplotypes ({len(inv) - len(skipped) - len(skipped)} new samples)")
        con.close()
        return

    # ── summary ───────────────────────────────────────────────────────────────
    print(f"\n── Import complete ─────────────────────────────────────────────")
    print(f"  New assemblies:  {n_asm}")
    print(f"  New haplotypes:  {n_hap}")
    print(f"  New copies:      {n_cop:,}")
    print(f"  New variants:    {n_var:,}")
    if skipped:
        print(f"  Skipped (already in DB): {len(skipped)}")

    print(f"\n── Database totals ─────────────────────────────────────────────")
    for tbl in ["assembly", "haplotype", "copy", "variant", "read_variant"]:
        n = con.execute(f"SELECT COUNT(*) FROM [{tbl}]").fetchone()[0]
        print(f"  {tbl:20s}: {n:>8,} rows")

    print(f"\n  Cohort breakdown:")
    for cohort, cnt in con.execute(
        "SELECT cohort, COUNT(*) FROM assembly GROUP BY cohort ORDER BY cohort"
    ).fetchall():
        print(f"    {cohort:20s}: {cnt} samples")

    con.close()
    print(f"\nSaved: {DB}  ({DB.stat().st_size/1024/1024:.1f} MB)")


if __name__ == "__main__":
    main()
