#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 76_import_ukbb_population.py — import UK Biobank 5S rDNA variant calls into the
# T2T coordinate system and store the full per-variant carrier VAF distribution.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
76_import_ukbb_population.py

Import UK Biobank (n=490,075) 5S rDNA variant calls into the T2T coordinate
system and store the full per-variant VAF distribution for threshold-free
downstream analysis.

Source:
  WGS-Variant-Identification/Trios_3.0/results_500k/carriers_AD1_with_dp.tsv
  (98.6 M rows; all per-sample calls at AD ≥ 1, UKBB consensus positions 450–950)

Coordinate transform:
  T2T_pos  = 1417 − UKBB_pos
  T2T_ref  = complement(UKBB_ref)     ← gene is on minus strand in both consensuses
  T2T_alt  = complement(UKBB_alt)
  Region (T2T consensus_pos):
    nts_post : > 748       (UKBB 450–668)
    gene     : 630–748     (UKBB 669–787)
    nts_pre  : < 630       (UKBB 788–950)

The 501-bp covered window (UKBB 450–950 ↔ T2T 467–967) is 100% identical
between the two consensus sequences (BLAST-verified; 0 mismatches, 0 gaps).

Design: per-variant VAF distribution
  The core column `vaf_array` is a sorted float32 numpy array (one value per
  carrier at AD ≥ 1), stored as a raw binary BLOB.  This allows applying any
  VAF threshold t post-hoc without re-importing:

    vafs = np.frombuffer(blob, dtype=np.float32)
    n_carriers_at_t = int(len(vafs) - np.searchsorted(vafs, t))

  Two integer convenience columns are also stored for fast SQL queries:
    n_carriers_ad1  – len(vaf_array), i.e. carriers at the baseline AD ≥ 1
    n_carriers_ad5  – the UKBB-validated germline threshold (trio heritability ≥ 80%)

  Mean and median VAF (computed from vaf_array) are stored for quick scanning.

Creates two tables in 5S_rDNA.db:
  ukbb_population_variants  – per-variant frequency data
  ukbb_depth_profile         – per-position depth stats in T2T coords

Idempotent: drops and recreates both tables on each run.

Memory note: accumulates ~394 MB of float32 data during processing; peak
RSS ~700 MB.  Runtime: ~3–5 min on Mac (≈20 chunks × 5 M rows each).
"""

import os
import sqlite3
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

# ── paths ─────────────────────────────────────────────────────────────────────

T2T    = Path(os.environ.get("FIVES_DATA", "data"))
UKBB   = Path(os.environ.get("FIVES_DATA", "data")) / "Trios_3.0/results_500k"
DB     = Path(os.environ.get("FIVES_DB", "5S_rDNA.db"))

CARRIERS_FILE   = UKBB / "carriers_AD1_with_dp.tsv"
DEPTH_FILE      = UKBB / "dp_by_position.tsv"

N_TOTAL_SAMPLES = 490_075
CHUNKSIZE       = 5_000_000

# ── coordinate helpers ────────────────────────────────────────────────────────

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

def region_of(t2t_pos: int) -> str:
    if t2t_pos < 630:  return "nts_pre"
    if t2t_pos <= 748: return "gene"
    return "nts_post"


# ── pass 1: accumulate per-variant VAF arrays ─────────────────────────────────
# Each variant accumulates a list of small float32 arrays (one per chunk).
# At the end we concatenate + sort once per variant.

print(f"Reading {CARRIERS_FILE.name}  "
      f"({CARRIERS_FILE.stat().st_size / 1e9:.1f} GB) …", flush=True)

# vaf_chunks[(t2t_pos, t2t_ref, t2t_alt)] = [array_chunk1, array_chunk2, ...]
vaf_chunks: dict[tuple, list[np.ndarray]] = defaultdict(list)

# AD ≥ 5 counts accumulated the cheap way (no extra memory vs what we already store)
ad5_counts: dict[tuple, int] = defaultdict(int)

n_chunks = 0
n_rows   = 0

for chunk in pd.read_csv(
        CARRIERS_FILE, sep="\t", chunksize=CHUNKSIZE,
        usecols=["POS", "REF", "ALT", "AD", "VAF"],
        dtype={"POS": int, "REF": str, "ALT": str, "AD": int, "VAF": float}):

    # Filter to covered window (guard against stray positions)
    chunk = chunk[(chunk["POS"] >= 450) & (chunk["POS"] <= 950)].copy()
    if chunk.empty:
        n_chunks += 1
        continue

    # Vectorised coordinate transform
    chunk["t2t_pos"] = 1417 - chunk["POS"]
    chunk["t2t_ref"] = chunk["REF"].str.translate(_COMP)
    chunk["t2t_alt"] = chunk["ALT"].str.translate(_COMP)

    # Accumulate VAF arrays per variant (all rows are AD ≥ 1 already)
    grp = chunk.groupby(["t2t_pos", "t2t_ref", "t2t_alt"], sort=False)
    for key, sub in grp:
        vaf_chunks[key].append(sub["VAF"].values.astype(np.float32))

    # AD ≥ 5 convenience count
    above5 = chunk[chunk["AD"] >= 5].groupby(
        ["t2t_pos", "t2t_ref", "t2t_alt"], sort=False).size()
    for key, cnt in above5.items():
        ad5_counts[key] += int(cnt)

    n_rows   += len(chunk)
    n_chunks += 1
    if n_chunks % 5 == 0:
        print(f"  chunk {n_chunks:3d}  rows: {n_rows:>12,}  "
              f"variants: {len(vaf_chunks):,}", flush=True)

print(f"Done reading. {n_rows:,} rows → {len(vaf_chunks):,} unique variants.", flush=True)


# ── pass 2: sort arrays + build DB rows ──────────────────────────────────────

print("Sorting VAF arrays …", flush=True)

rows = []
for (t2t_pos, t2t_ref, t2t_alt), arrays in vaf_chunks.items():
    vafs = np.sort(np.concatenate(arrays))   # sorted float32
    n_ad1 = len(vafs)
    n_ad5 = int(ad5_counts.get((t2t_pos, t2t_ref, t2t_alt), 0))
    mean_vaf  = float(vafs.mean())
    median_vaf = float(np.median(vafs))
    rows.append((
        int(t2t_pos), t2t_ref, t2t_alt,
        1417 - int(t2t_pos),        # ukbb_pos
        region_of(int(t2t_pos)),
        n_ad1, n_ad5,
        N_TOTAL_SAMPLES,
        mean_vaf, median_vaf,
        vafs.tobytes(),             # BLOB
    ))

print(f"{len(rows):,} variants ready for import.", flush=True)


# ── depth profile ─────────────────────────────────────────────────────────────

print("Reading depth profile …", flush=True)
dp = pd.read_csv(DEPTH_FILE, sep="\t")
dp_rows = [
    (int(1417 - r.POS), int(r.POS), int(r.N_SAMPLES),
     float(r.MEAN_DP), float(r.MEDIAN_DP),
     float(r.P5_DP), float(r.P25_DP), float(r.P75_DP), float(r.P95_DP),
     region_of(int(1417 - r.POS)))
    for _, r in dp.iterrows()
]


# ── write to SQLite ───────────────────────────────────────────────────────────

print(f"Writing to {DB} …", flush=True)
con = sqlite3.connect(DB)

con.executescript("""
DROP TABLE IF EXISTS ukbb_population_variants;
CREATE TABLE ukbb_population_variants (
    ukbb_pop_id      INTEGER PRIMARY KEY AUTOINCREMENT,
    t2t_pos          INTEGER NOT NULL,
    t2t_ref          TEXT    NOT NULL,
    t2t_alt          TEXT    NOT NULL,
    ukbb_pos         INTEGER NOT NULL,
    region           TEXT,
    n_carriers_ad1   INTEGER DEFAULT 0,
    n_carriers_ad5   INTEGER DEFAULT 0,
    n_total_samples  INTEGER DEFAULT 490075,
    mean_vaf         REAL,
    median_vaf       REAL,
    vaf_array        BLOB NOT NULL,
    UNIQUE(t2t_pos, t2t_ref, t2t_alt)
);
CREATE INDEX IF NOT EXISTS ix_ukbb_pop_pos
    ON ukbb_population_variants(t2t_pos);

DROP TABLE IF EXISTS ukbb_depth_profile;
CREATE TABLE ukbb_depth_profile (
    t2t_pos   INTEGER PRIMARY KEY,
    ukbb_pos  INTEGER NOT NULL,
    n_samples INTEGER,
    mean_dp   REAL,
    median_dp REAL,
    p5_dp     REAL,
    p25_dp    REAL,
    p75_dp    REAL,
    p95_dp    REAL,
    region    TEXT
);
""")

con.executemany("""
    INSERT INTO ukbb_population_variants
        (t2t_pos, t2t_ref, t2t_alt, ukbb_pos, region,
         n_carriers_ad1, n_carriers_ad5,
         n_total_samples, mean_vaf, median_vaf, vaf_array)
    VALUES (?,?,?,?,?,?,?,?,?,?,?)
""", rows)

con.executemany("""
    INSERT INTO ukbb_depth_profile
        (t2t_pos, ukbb_pos, n_samples, mean_dp, median_dp,
         p5_dp, p25_dp, p75_dp, p95_dp, region)
    VALUES (?,?,?,?,?,?,?,?,?,?)
""", dp_rows)

con.commit()

# ── summary ───────────────────────────────────────────────────────────────────

print("\n── import summary ───────────────────────────────────────────────")
for region, n in con.execute("""
        SELECT region, COUNT(*) FROM ukbb_population_variants
        GROUP BY region ORDER BY region"""):
    print(f"  {region or 'NULL':10s}: {n:,} variants")

gene_common = con.execute("""
    SELECT COUNT(*) FROM ukbb_population_variants
    WHERE region='gene' AND n_carriers_ad5 >= 100""").fetchone()[0]
print(f"  Gene-body variants with ≥100 UKBB carriers at AD≥5: {gene_common:,}")

blob_size = con.execute("""
    SELECT SUM(LENGTH(vaf_array)) FROM ukbb_population_variants
""").fetchone()[0]
print(f"  Total BLOB storage: {blob_size / 1e6:.1f} MB")

print("\nUsage example:")
print("  import numpy as np, sqlite3")
print("  vafs = np.frombuffer(blob, dtype=np.float32)")
print("  n_carriers_at_vaf_003 = int(len(vafs) - np.searchsorted(vafs, 0.003))")

con.close()
print("Done.")
