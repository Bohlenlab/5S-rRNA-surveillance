#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# build_database.py — build the 5S rDNA SQLite database from assembly-analysis outputs.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
Build the 5S rDNA SQLite database from HPRC assembly analysis outputs.

Tables (hierarchical):
  array_reference  — one row: the T2T 5S rDNA consensus unit
  assembly         — one row per individual (sample_id)
  haplotype        — one row per (sample_id, hap_label)
  copy             — one row per repeat copy within a haplotype
  variant          — per-copy variants from three MAFFT alignments
  read_variant     — per-sample variants from read-based pipelines (LR / SR)

Variant positions are always 1-based relative to the T2T consensus unit
(length 2168 bp, gene at positions 630–748).

Sources imported:
  databases/*_hap*.tsv       — per-haplotype copy/variant tables
  blast/*_blast.txt          — strand determination
  bam/*/*_hifi_summary.tsv   — HiFi read-variant summaries
  data_inventory.tsv         — sample metadata

Paths are read from environment variables:
  FIVES_DB    path to 5S_rDNA.db
  FIVES_DATA  input derived-data directory (databases/, blast/, bam/, inventory)
  FIVES_REFS  reference directory (consensus FASTA)

Usage:
  python3 build_database.py
"""

import os
import re
import sqlite3
from pathlib import Path

import pandas as pd

T2T  = Path(os.environ.get("FIVES_DATA", "data"))
HPRC = T2T
CONS = Path(os.environ.get("FIVES_REFS", "refs")) / "5S_t2t_consensus.fa"
INV  = HPRC / "data_inventory.tsv"
DB   = Path(os.environ.get("FIVES_DB", "5S_rDNA.db"))

GENE_START_1     = 630   # 1-based inclusive
GENE_END_1       = 748   # 1-based inclusive
NTS_POST_START_1 = 749   # 1-based start of NTS-post

BLAST_COLS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
              "qstart", "qend", "sstart", "send", "evalue", "bitscore",
              "qlen", "slen", "sstrand"]


# ── helpers ───────────────────────────────────────────────────────────────────

def region_of(pos1: int) -> str:
    if pos1 < GENE_START_1:
        return "nts_pre"
    if pos1 <= GENE_END_1:
        return "gene"
    return "nts_post"


def parse_variants_str(s: str) -> list:
    """Parse "pos:ref>alt; ..." → list of (0-based int pos, str ref, str alt)."""
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
    """Return 'plus' or 'minus' from the majority sstrand in a BLAST output file."""
    if not blast_file.exists():
        return "unknown"
    try:
        df = pd.read_csv(blast_file, sep="\t", names=BLAST_COLS, usecols=["sstrand"])
        if not df.empty:
            return str(df["sstrand"].value_counts().idxmax())
    except Exception:
        pass
    return "unknown"


def read_consensus_seq() -> str:
    lines = CONS.read_text().splitlines()
    return "".join(l.strip() for l in lines if not l.startswith(">"))


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


# ── schema ────────────────────────────────────────────────────────────────────

SCHEMA = """
PRAGMA journal_mode = WAL;
PRAGMA foreign_keys = ON;

CREATE TABLE IF NOT EXISTS array_reference (
    array_id       INTEGER PRIMARY KEY,
    name           TEXT NOT NULL,          -- "5S_rDNA"
    ref_label      TEXT NOT NULL,          -- FASTA header of the reference unit
    sequence       TEXT NOT NULL,          -- full T2T consensus sequence (2168 bp)
    length_bp      INTEGER NOT NULL,
    nts_pre_start  INTEGER NOT NULL,       -- 1 (1-based, inclusive)
    nts_pre_end    INTEGER NOT NULL,       -- 629
    gene_start     INTEGER NOT NULL,       -- 630
    gene_end       INTEGER NOT NULL,       -- 748
    nts_post_start INTEGER NOT NULL,       -- 749
    nts_post_end   INTEGER NOT NULL        -- 2168
);

CREATE TABLE IF NOT EXISTS assembly (
    assembly_id      INTEGER PRIMARY KEY AUTOINCREMENT,
    sample_id        TEXT NOT NULL UNIQUE,
    cohort           TEXT NOT NULL,        -- "HPRC_Year1", "CHM13", "HG002_GIAB"
    population       TEXT,                 -- e.g. "CHS"
    superpopulation  TEXT,                 -- e.g. "EAS"
    sex              TEXT,                 -- "M", "F", NULL if unknown
    age              INTEGER,              -- NULL if unknown
    assembly_hap1_s3 TEXT,
    assembly_hap2_s3 TEXT,
    has_hifi         INTEGER NOT NULL DEFAULT 0,
    has_illumina     INTEGER NOT NULL DEFAULT 0,
    has_methylation  INTEGER NOT NULL DEFAULT 0,
    has_rnaseq       INTEGER NOT NULL DEFAULT 0
);

CREATE TABLE IF NOT EXISTS haplotype (
    haplotype_id      INTEGER PRIMARY KEY AUTOINCREMENT,
    assembly_id       INTEGER NOT NULL REFERENCES assembly(assembly_id),
    array_id          INTEGER NOT NULL REFERENCES array_reference(array_id),
    hap_label         TEXT NOT NULL,       -- "hap1" or "hap2"
    array_chrom       TEXT NOT NULL,       -- contig / chromosome name in the assembly
    strand            TEXT NOT NULL,       -- "plus" or "minus" (from BLAST sstrand)
    n_copies          INTEGER NOT NULL,
    array_start_local INTEGER,             -- leftmost coord of array in assembly (1-based)
    array_end_local   INTEGER,             -- rightmost coord
    UNIQUE(assembly_id, hap_label)
);

CREATE TABLE IF NOT EXISTS copy (
    copy_id            INTEGER PRIMARY KEY AUTOINCREMENT,
    haplotype_id       INTEGER NOT NULL REFERENCES haplotype(haplotype_id),
    copy_number        INTEGER NOT NULL,   -- 1-based, 5'→3' order (strand-normalised)
    unit_start_local   INTEGER,            -- start of full repeat unit in assembly (1-based)
    unit_end_local     INTEGER,            -- end
    unit_length_bp     INTEGER,
    spacing_to_next_bp REAL,               -- NULL for the last copy in array
    gene_lo_local      INTEGER,            -- 5S gene start in assembly coords
    gene_hi_local      INTEGER,            -- 5S gene end
    gene_pct_identity  REAL,
    gene_mismatches    INTEGER,
    gene_gaps          INTEGER,
    n_snv_gene         INTEGER,            -- variants in full-unit MAFFT alignment
    n_snv_5s_gene      INTEGER,            -- subset in 5S gene region (pos 630–748)
    n_snv_nts_pre      INTEGER,            -- variants in NTS-pre MAFFT alignment
    n_snv_nts_post     INTEGER,            -- variants in NTS-post MAFFT alignment
    category           TEXT,               -- "identical" (no gene variants) | "highly_similar"
    border_note        TEXT,               -- "interior" | "5-prime_array_border" | "3-prime_array_border"
    UNIQUE(haplotype_id, copy_number)
);

-- Assembly-derived variants (from MAFFT alignments).
-- alignment_source disambiguates the three separate alignments:
--   "gene_unit"  : full 2168 bp repeat unit alignment (consensus_pos = 0-based + 1)
--   "nts_pre_aln": NTS-pre only alignment             (consensus_pos = 0-based + 1)
--   "nts_post_aln": NTS-post only alignment            (consensus_pos = 749 + 0-based)
-- ref/alt are relative to the MAFFT majority-vote consensus of that alignment,
-- which closely approximates the T2T reference sequence.
CREATE TABLE IF NOT EXISTS variant (
    variant_id       INTEGER PRIMARY KEY AUTOINCREMENT,
    copy_id          INTEGER NOT NULL REFERENCES copy(copy_id),
    alignment_source TEXT NOT NULL,
    consensus_pos    INTEGER NOT NULL,     -- 1-based position in T2T consensus unit
    ref              TEXT NOT NULL,
    alt              TEXT NOT NULL,
    region           TEXT NOT NULL        -- "nts_pre" | "gene" | "nts_post"
);

-- Read-derived variants (LR HiFi, Illumina SR, RNA-seq).
-- consensus_pos is the bcftools 1-based mpileup position (already 1-based, no conversion).
CREATE TABLE IF NOT EXISTS read_variant (
    read_variant_id INTEGER PRIMARY KEY AUTOINCREMENT,
    assembly_id     INTEGER NOT NULL REFERENCES assembly(assembly_id),
    modality        TEXT NOT NULL,         -- "hifi" | "illumina" | "rnaseq"
    consensus_pos   INTEGER NOT NULL,      -- 1-based T2T consensus position
    ref             TEXT NOT NULL,
    alt             TEXT NOT NULL,
    depth           INTEGER,
    alt_depth       INTEGER,
    vaf             REAL,
    region          TEXT                   -- "nts_pre" | "gene" | "nts_post"
);

CREATE INDEX IF NOT EXISTS idx_variant_copy   ON variant(copy_id);
CREATE INDEX IF NOT EXISTS idx_variant_pos    ON variant(consensus_pos);
CREATE INDEX IF NOT EXISTS idx_variant_source ON variant(alignment_source);
CREATE INDEX IF NOT EXISTS idx_read_var_asm   ON read_variant(assembly_id);
CREATE INDEX IF NOT EXISTS idx_read_var_pos   ON read_variant(consensus_pos);
CREATE INDEX IF NOT EXISTS idx_copy_hap       ON copy(haplotype_id);
CREATE INDEX IF NOT EXISTS idx_hap_asm        ON haplotype(assembly_id);
"""


# ── build ─────────────────────────────────────────────────────────────────────

def build_db():
    if DB.exists():
        DB.unlink()
        print(f"Removed existing {DB.name}")

    con = sqlite3.connect(DB)
    con.executescript(SCHEMA)
    print(f"Schema created: {DB}")

    # ── 1. array_reference ───────────────────────────────────────────────────
    seq = read_consensus_seq()
    con.execute("""
        INSERT INTO array_reference
          (array_id, name, ref_label, sequence, length_bp,
           nts_pre_start, nts_pre_end, gene_start, gene_end,
           nts_post_start, nts_post_end)
        VALUES (1,'5S_rDNA','5S_rDNA_consensus_CHM13',?,?,1,629,630,748,749,?)
    """, (seq, len(seq), len(seq)))
    print(f"  array_reference: 1 row  (consensus length = {len(seq)} bp)")

    # ── 2. assembly ──────────────────────────────────────────────────────────
    inv = pd.read_csv(INV, sep="\t")
    for _, row in inv.iterrows():
        con.execute("""
            INSERT OR IGNORE INTO assembly
              (sample_id, cohort, population, superpopulation,
               assembly_hap1_s3, assembly_hap2_s3,
               has_hifi, has_illumina, has_methylation, has_rnaseq)
            VALUES (?,?,?,?,?,?,?,?,?,?)
        """, (
            safe_str(row["sample_id"]),
            "HPRC_Year1",
            safe_str(row.get("population")),
            safe_str(row.get("superpopulation")),
            safe_str(row.get("assembly_hap1_s3")),
            safe_str(row.get("assembly_hap2_s3")),
            int(bool(row.get("has_hifi",      False))),
            int(bool(row.get("has_illumina",  False))),
            int(bool(row.get("has_hifi_meth", False))),
            int(bool(row.get("has_rnaseq",    False))),
        ))
    con.commit()
    n_asm = con.execute("SELECT COUNT(*) FROM assembly").fetchone()[0]
    print(f"  assembly:        {n_asm} rows")

    # ── 3+4. haplotype + copy + variant ──────────────────────────────────────
    db_files = sorted((HPRC / "databases").glob("*.tsv"))
    n_hap = n_cop = n_var = 0

    for db_file in db_files:
        stem = db_file.stem   # e.g. "HG00438_hap1"
        parts = stem.rsplit("_", 1)
        if len(parts) != 2:
            continue
        sample_id, hap_label = parts

        df = pd.read_csv(db_file, sep="\t")
        if df.empty:
            print(f"  WARNING: empty database file {db_file.name}")
            continue

        asm_row = con.execute(
            "SELECT assembly_id FROM assembly WHERE sample_id=?", (sample_id,)
        ).fetchone()
        if asm_row is None:
            print(f"  WARNING: {sample_id} not in assembly table — skipping")
            continue
        assembly_id = asm_row[0]

        # True strand from BLAST file
        blast_file = HPRC / "blast" / f"{sample_id}_{hap_label}_blast.txt"
        strand = get_strand(blast_file)

        array_chrom = safe_str(df["array_chrom"].iloc[0]) or "unknown"
        n_copies = len(df)

        # Array span in assembly coordinates
        all_coords = (
            df["unit_start_local"].dropna().tolist() +
            df["unit_end_local"].dropna().tolist()
        )
        arr_lo = int(min(all_coords)) if all_coords else None
        arr_hi = int(max(all_coords)) if all_coords else None

        cur = con.execute("""
            INSERT INTO haplotype
              (assembly_id, array_id, hap_label, array_chrom, strand,
               n_copies, array_start_local, array_end_local)
            VALUES (?,1,?,?,?,?,?,?)
        """, (assembly_id, hap_label, array_chrom, strand, n_copies, arr_lo, arr_hi))
        haplotype_id = cur.lastrowid
        n_hap += 1

        for _, r in df.iterrows():
            copy_num = int(r["copy_id"])

            cur2 = con.execute("""
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
            copy_db_id = cur2.lastrowid
            n_cop += 1

            variant_rows = []

            # gene_unit: full 2168 bp unit alignment, 0-based pos → +1 for 1-based
            for pos0, ref, alt in parse_variants_str(str(r.get("gene_variants", ""))):
                cp1 = pos0 + 1
                variant_rows.append((copy_db_id, "gene_unit", cp1, ref, alt, region_of(cp1)))

            # nts_pre_aln: NTS-pre only alignment, 0-based pos → +1 for 1-based (range 1–629)
            for pos0, ref, alt in parse_variants_str(str(r.get("nts_pre_variants", ""))):
                cp1 = pos0 + 1
                variant_rows.append((copy_db_id, "nts_pre_aln", cp1, ref, alt, "nts_pre"))

            # nts_post_aln: NTS-post only alignment, 0-based pos → +749 for 1-based (range 749–2168)
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
    print(f"  haplotype:       {n_hap} rows")
    print(f"  copy:            {n_cop:,} rows")
    print(f"  variant:         {n_var:,} rows")

    # ── 5. read_variant (LR HiFi summaries, locally available) ───────────────
    n_rv = 0
    for lr_file in sorted(HPRC.glob("bam/*/*_hifi_summary.tsv")):
        sample_id = lr_file.parent.name
        asm_row = con.execute(
            "SELECT assembly_id FROM assembly WHERE sample_id=?", (sample_id,)
        ).fetchone()
        if asm_row is None:
            print(f"  WARNING: {sample_id} not in assembly (read_variant skipped)")
            continue
        assembly_id = asm_row[0]

        try:
            lr = pd.read_csv(lr_file, sep="\t")
        except Exception as e:
            print(f"  WARNING: could not read {lr_file.name}: {e}")
            continue

        rows = []
        for _, r in lr.iterrows():
            pos = safe_int(r.get("pos"))   # bcftools pos is already 1-based
            if pos is None:
                continue
            rows.append((
                assembly_id, "hifi", pos,
                safe_str(r.get("ref")), safe_str(r.get("alt")),
                safe_int(r.get("dp")),
                safe_int(r.get("ad")),
                safe_float(r.get("vaf")),
                safe_str(r.get("region")),
            ))
        if rows:
            con.executemany("""
                INSERT INTO read_variant
                  (assembly_id, modality, consensus_pos, ref, alt,
                   depth, alt_depth, vaf, region)
                VALUES (?,?,?,?,?,?,?,?,?)
            """, rows)
            n_rv += len(rows)
            print(f"    {sample_id}: {len(rows)} HiFi read variants")

    con.commit()
    print(f"  read_variant:    {n_rv} rows")

    # ── summary ───────────────────────────────────────────────────────────────
    print("\n── Database summary ─────────────────────────────────────────────")
    for tbl in ["array_reference", "assembly", "haplotype", "copy", "variant", "read_variant"]:
        n = con.execute(f"SELECT COUNT(*) FROM [{tbl}]").fetchone()[0]
        print(f"  {tbl:20s}: {n:>8,} rows")

    # Spot-check: copies per haplotype
    stats = con.execute("""
        SELECT MIN(n_copies), MAX(n_copies), AVG(n_copies)
        FROM haplotype
    """).fetchone()
    print(f"\n  Copies/haplotype: min={stats[0]}  max={stats[1]}  mean={stats[2]:.1f}")

    # Variants per alignment source
    print("\n  Variants by alignment_source:")
    for src, cnt in con.execute("""
        SELECT alignment_source, COUNT(*) FROM variant GROUP BY alignment_source ORDER BY 1
    """).fetchall():
        print(f"    {src:20s}: {cnt:>8,}")

    # Variants by region
    print("\n  Variants by region:")
    for rgn, cnt in con.execute("""
        SELECT region, COUNT(*) FROM variant GROUP BY region ORDER BY 1
    """).fetchall():
        print(f"    {rgn:20s}: {cnt:>8,}")

    con.close()
    print(f"\nSaved: {DB}  ({DB.stat().st_size/1024/1024:.1f} MB)")


if __name__ == "__main__":
    build_db()
