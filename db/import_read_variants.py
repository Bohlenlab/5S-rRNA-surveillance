#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# import_read_variants.py — load HiFi and Illumina read-level variants into the read_variant table.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
Import HiFi (LR) and Illumina (SR) read-level variants into 5S_rDNA.db.

Sources:
  HPRC/bam/{sample}/{sample}_hifi_summary.tsv  → modality='hifi'
  HPRC/bam/{sample}/{sample}_sr_vs_lr.tsv      → modality='illumina'

Idempotent: skips (assembly_id, modality) pairs already in read_variant.

Usage:
  python3 import_read_variants.py
  python3 import_read_variants.py --dry-run
"""

import os
import sqlite3
import sys
from pathlib import Path

import pandas as pd

T2T  = Path(os.environ.get("FIVES_DATA", "data"))
HPRC = T2T
DB   = Path(os.environ.get("FIVES_DB", "5S_rDNA.db"))
BAM  = HPRC / "bam"

dry_run = "--dry-run" in sys.argv

con = sqlite3.connect(DB)
con.execute("PRAGMA foreign_keys = ON")
con.execute("PRAGMA journal_mode = WAL")

# Build sample_id → assembly_id map (prefer HPRC_Year1 over GIAB duplicate)
asm_map = {}
for sid, aid, cohort in con.execute(
    "SELECT sample_id, assembly_id, cohort FROM assembly ORDER BY cohort"
).fetchall():
    if sid not in asm_map or cohort == "HPRC_Year1":
        asm_map[sid] = aid

# Already-imported (assembly_id, modality) pairs
done = set(con.execute(
    "SELECT assembly_id, modality FROM read_variant GROUP BY assembly_id, modality"
).fetchall())

n_hifi = n_sr = n_skip = 0

for sample_dir in sorted(BAM.iterdir()):
    if not sample_dir.is_dir():
        continue
    sid = sample_dir.name

    if sid not in asm_map:
        continue
    aid = asm_map[sid]

    for tsv, modality in [
        (sample_dir / f"{sid}_hifi_summary.tsv",  "hifi"),
        (sample_dir / f"{sid}_sr_vs_lr.tsv",      "illumina"),
    ]:
        if not tsv.exists() or tsv.stat().st_size < 50:
            continue

        if (aid, modality) in done:
            n_skip += 1
            continue

        df = pd.read_csv(tsv, sep="\t")
        if df.empty:
            continue

        rows = []
        for _, r in df.iterrows():
            rows.append((
                aid,
                modality,
                int(r["pos"]),
                str(r["ref"]),
                str(r["alt"]),
                int(r["dp"])   if pd.notna(r.get("dp"))  else None,
                int(r["ad"])   if pd.notna(r.get("ad"))  else None,
                float(r["vaf"]) if pd.notna(r.get("vaf")) else None,
                str(r["region"]) if pd.notna(r.get("region")) else None,
            ))

        if dry_run:
            print(f"  [dry-run] {sid} {modality}: {len(rows)} variants")
            if modality == "hifi":
                n_hifi += len(rows)
            else:
                n_sr += len(rows)
            continue

        con.executemany("""
            INSERT INTO read_variant
              (assembly_id, modality, consensus_pos, ref, alt,
               depth, alt_depth, vaf, region)
            VALUES (?,?,?,?,?,?,?,?,?)
        """, rows)
        con.commit()
        print(f"  {sid} {modality}: {len(rows)} variants imported")
        if modality == "hifi":
            n_hifi += len(rows)
        else:
            n_sr += len(rows)

print(f"\n── Import complete ──────────────────────────────────────")
print(f"  HiFi variants:      {n_hifi:>6,}")
print(f"  Illumina variants:  {n_sr:>6,}")
print(f"  Skipped (done):     {n_skip:>6,}")

for mod, cnt, n_samp in con.execute(
    "SELECT modality, COUNT(*), COUNT(DISTINCT assembly_id) FROM read_variant GROUP BY modality"
).fetchall():
    print(f"  DB {mod:10s}: {cnt:>7,} rows  ({n_samp} samples)")

con.close()
