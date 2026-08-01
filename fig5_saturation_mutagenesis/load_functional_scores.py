#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# load_functional_scores.py — Load saturation-mutagenesis RNA-expression and 60S-
# incorporation scores into the functional_annotation table of 5S_rDNA.db.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
Load saturation mutagenesis functional data into 5S_rDNA.db.

RNA_expression.csv  — "Ref" = WT/reference base; values for each alt base
                       in non-Ref columns (2 replicates). WT (Ref) col empty.
60S_incorporation.csv — "Ref" = WT/reference base; values for each alt base
                         in non-Ref columns (2 replicates). WT (Ref) col empty.

gene_pos (1-based within 5S gene)  →  consensus_pos = 629 + gene_pos
"""

import os
import sqlite3
from pathlib import Path
import numpy as np
import pandas as pd

BASE   = Path(os.environ.get("FIVES_DATA", "data"))
DB     = Path(os.environ.get("FIVES_DB", "5S_rDNA.db"))
RNA_F  = BASE / "RNA_expression.csv"
S60_F  = BASE / "60S_incorporation.csv"

GENE_OFFSET = 629  # gene_pos 1 → consensus_pos 630

# ─── parse RNA_expression.csv ─────────────────────────────────────────────────
# 10 columns: Ref, position, A_r1, A_r2, C_r1, C_r2, G_r1, G_r2, T_r1, T_r2
# Same format as 60S: Ref = WT base, values for non-Ref alts only
rna_cols = ['ref_base', 'gene_pos',
            'A_r1','A_r2',
            'C_r1','C_r2',
            'G_r1','G_r2',
            'T_r1','T_r2']
rna_raw  = pd.read_csv(RNA_F, header=0, names=rna_cols, encoding='utf-8-sig')
rna_raw  = rna_raw[pd.to_numeric(rna_raw['gene_pos'], errors='coerce').notna()].copy()
rna_raw['gene_pos'] = rna_raw['gene_pos'].astype(int)

ref_at_pos_rna = {}  # gene_pos → WT base (from RNA Ref col)
rna_dict = {}        # (gene_pos, alt_base) → (r1, r2, mean)
for _, row in rna_raw.iterrows():
    ref = row['ref_base']
    pos = row['gene_pos']
    if ref not in ('A','C','G','T'):
        continue
    ref_at_pos_rna[pos] = ref
    for alt in ('A','C','G','T'):
        if alt == ref:
            continue
        r1, r2 = row[f'{alt}_r1'], row[f'{alt}_r2']
        vals = [v for v in (r1, r2) if pd.notna(v)]
        if not vals:
            continue
        mean = float(np.mean(vals))
        rna_dict[(pos, alt)] = (
            float(r1) if pd.notna(r1) else None,
            float(r2) if pd.notna(r2) else None,
            mean
        )

print(f"RNA: {len(rna_dict)} (pos, alt) entries; WT bases at {len(ref_at_pos_rna)} positions")

# ─── parse 60S_incorporation.csv ─────────────────────────────────────────────
# 10 columns: Ref, Position, A_r1, A_r2, C_r1, C_r2, G_r1, G_r2, T_r1, T_r2
s60_cols = ['ref_base', 'gene_pos',
            'A_r1','A_r2',
            'C_r1','C_r2',
            'G_r1','G_r2',
            'T_r1','T_r2']
s60_raw  = pd.read_csv(S60_F, header=0, names=s60_cols)
s60_raw  = s60_raw[pd.to_numeric(s60_raw['gene_pos'], errors='coerce').notna()].copy()
s60_raw['gene_pos'] = s60_raw['gene_pos'].astype(int)

ref_at_pos = {}          # gene_pos → WT base (from 60S Ref col)
s60_dict   = {}          # (gene_pos, alt_base) → (r1, r2, mean)

for _, row in s60_raw.iterrows():
    ref = row['ref_base']
    pos = row['gene_pos']
    if ref not in ('A','C','G','T'):
        continue
    ref_at_pos[pos] = ref
    for alt in ('A','C','G','T'):
        if alt == ref:
            continue
        r1, r2 = row[f'{alt}_r1'], row[f'{alt}_r2']
        vals = [v for v in (r1, r2) if pd.notna(v)]
        if not vals:
            continue
        mean = float(np.mean(vals))
        s60_dict[(pos, alt)] = (
            float(r1) if pd.notna(r1) else None,
            float(r2) if pd.notna(r2) else None,
            mean
        )

print(f"60S: {len(s60_dict)} (pos, alt) entries; WT bases at {len(ref_at_pos)} positions")

# Merge ref_at_pos from both sources (should agree; 60S takes precedence)
ref_at_pos = {**ref_at_pos_rna, **ref_at_pos}
mismatches = [(p, ref_at_pos_rna[p], ref_at_pos[p])
              for p in ref_at_pos_rna if p in ref_at_pos and ref_at_pos_rna[p] != ref_at_pos[p]]
if mismatches:
    print(f"WARNING: {len(mismatches)} ref-base mismatches between RNA and 60S:")
    for p, r, s in mismatches[:5]:
        print(f"  pos {p}: RNA says {r}, 60S says {s}")

# ─── build merged records ─────────────────────────────────────────────────────
all_pos_alts = set(rna_dict.keys()) | set(s60_dict.keys())
rows = []
for (pos, alt) in sorted(all_pos_alts):
    ref = ref_at_pos.get(pos, None)
    if ref == alt:
        continue
    cpos = GENE_OFFSET + pos
    rna = rna_dict.get((pos, alt))
    s60 = s60_dict.get((pos, alt))
    rows.append((
        pos, cpos, ref, alt,
        rna[0] if rna else None,
        rna[1] if rna else None,
        None,                        # rna_expr_rep3 not present in new format
        rna[2] if rna else None,     # mean
        s60[0] if s60 else None,
        s60[1] if s60 else None,
        s60[2] if s60 else None,
    ))

print(f"Total (pos, alt) pairs to insert: {len(rows)}")

# ─── insert into SQLite ───────────────────────────────────────────────────────
con = sqlite3.connect(DB)
con.execute("DROP TABLE IF EXISTS functional_annotation")
con.execute("""
CREATE TABLE functional_annotation (
    gene_pos        INTEGER NOT NULL,
    consensus_pos   INTEGER NOT NULL,
    ref_base        TEXT,
    alt_base        TEXT NOT NULL,
    rna_expr_rep1   REAL,
    rna_expr_rep2   REAL,
    rna_expr_rep3   REAL,
    rna_expr_mean   REAL,
    incorp_60s_rep1 REAL,
    incorp_60s_rep2 REAL,
    incorp_60s_mean REAL,
    PRIMARY KEY (gene_pos, alt_base)
)
""")
con.executemany("""
INSERT INTO functional_annotation VALUES (?,?,?,?,?,?,?,?,?,?,?)
""", rows)
con.commit()

n = con.execute("SELECT COUNT(*) FROM functional_annotation").fetchone()[0]
print(f"Inserted {n} rows into functional_annotation")
n_rna = con.execute("SELECT COUNT(*) FROM functional_annotation WHERE rna_expr_mean IS NOT NULL").fetchone()[0]
n_60s = con.execute("SELECT COUNT(*) FROM functional_annotation WHERE incorp_60s_mean IS NOT NULL").fetchone()[0]
n_both = con.execute("SELECT COUNT(*) FROM functional_annotation WHERE rna_expr_mean IS NOT NULL AND incorp_60s_mean IS NOT NULL").fetchone()[0]
print(f"  RNA only: {n_rna - n_both}  |  60S only: {n_60s - n_both}  |  Both: {n_both}")

con.close()
print("Done.")
