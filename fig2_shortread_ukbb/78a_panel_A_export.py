#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 78a_panel_A_export.py — export the UKBB 5S carrier-landscape panel (carriers at
# VAF >= 0.30% vs consensus position, with a read-depth overlay) at 300 DPI.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""Export the UKBB 5S carrier-landscape panel at 300 DPI."""

import os
import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

T2T    = Path(os.environ.get("FIVES_DATA", "data"))
DB     = Path(os.environ.get("FIVES_DB", "5S_rDNA.db"))
OUTDIR = Path(os.environ.get("FIVES_OUT", "output")) / "02_variant_calling_qc"

WIN_LO, WIN_HI = 467, 967
GENE_LO, GENE_HI = 630, 748
N_TOTAL = 490_075
THRESH_MAIN = 0.003

GENE_COLOR = "#ffe14d"
REG_COLORS = {"gene": "#e67e22", "nts_pre": "#888888", "nts_post": "#888888"}
REG_ORDER  = ["gene", "nts_pre", "nts_post"]

print("Loading data …", flush=True)
con = sqlite3.connect(DB)

raw = con.execute("""
    SELECT t2t_pos, t2t_ref, t2t_alt, region,
           n_carriers_ad1, n_carriers_ad5, mean_vaf, median_vaf, vaf_array
    FROM ukbb_population_variants ORDER BY t2t_pos
""").fetchall()

depth_df = pd.DataFrame(con.execute("""
    SELECT t2t_pos, median_dp, p5_dp, p25_dp, p75_dp, p95_dp
    FROM ukbb_depth_profile ORDER BY t2t_pos
""").fetchall(), columns=["t2t_pos","median_dp","p5_dp","p25_dp","p75_dp","p95_dp"])

tp_rows = con.execute("""
    SELECT DISTINCT consensus_pos, ref, alt FROM (
        SELECT v.consensus_pos, v.ref, v.alt
        FROM variant v
        JOIN copy c ON v.copy_id = c.copy_id
        JOIN haplotype h ON c.haplotype_id = h.haplotype_id
        JOIN assembly a ON h.assembly_id = a.assembly_id
        WHERE v.alignment_source = 'gene_unit_t2t'
          AND a.cohort = 'HPRC_Year1'
          AND v.consensus_pos BETWEEN ? AND ?
        UNION
        SELECT rv.consensus_pos, rv.ref, rv.alt
        FROM read_variant rv
        JOIN assembly a ON rv.assembly_id = a.assembly_id
        WHERE rv.modality = 'hifi'
          AND rv.vaf IS NOT NULL
          AND a.cohort = 'HPRC_Year1'
          AND rv.consensus_pos BETWEEN ? AND ?
    )
""", (WIN_LO, WIN_HI, WIN_LO, WIN_HI)).fetchall()
tp_set = {(int(r[0]), r[1], r[2]) for r in tp_rows}
con.close()

records = []
for pos, ref, alt, reg, n1, n5, mv, mdv, blob in raw:
    vafs = np.frombuffer(blob, dtype=np.float32)
    nc   = int(len(vafs) - np.searchsorted(vafs, THRESH_MAIN))
    records.append(dict(t2t_pos=int(pos), ref=ref, alt=alt,
                        region=reg or "unknown", n_main=nc))

df = pd.DataFrame(records)
df["is_tp"] = [(int(r["t2t_pos"]), r["ref"], r["alt"]) in tp_set for r in records]
df_called   = df[df["n_main"] >= 1].copy()
print(f"  {len(df_called):,} variants called at 0.30%  |  {df_called['is_tp'].sum():,} TP")

# ── data table into the figure folder ─────────────────────────────────────────
DATADIR = OUTDIR / "data"; DATADIR.mkdir(parents=True, exist_ok=True)
df_called[["t2t_pos","ref","alt","region","n_main","is_tp"]].rename(columns={
    "n_main":"carriers_vaf_ge_0p3pct","is_tp":"assembly_or_hifi_confirmed"}).to_csv(
    DATADIR/"Fig_78a_carrier_landscape.tsv", sep="\t", index=False)
print("data table ->", DATADIR/"Fig_78a_carrier_landscape.tsv")

# ── Panel A figure (11.5 cm × 4 cm) ──────────────────────────────────────────
CM = 1 / 2.54
fig, ax_A = plt.subplots(figsize=(11.5 * CM, 4 * CM))

for reg in REG_ORDER:
    sub    = df_called[(df_called["region"] == reg) & ~df_called["is_tp"]]
    sub_tp = df_called[(df_called["region"] == reg) &  df_called["is_tp"]]
    for subset, edge, lw, zo in [(sub, "none", 0, 3), (sub_tp, "#c0392b", 0.35, 4)]:
        sizes = np.clip(np.log10(subset["n_main"].values + 1) * 2.5, 0.8, 17)
        ax_A.scatter(subset["t2t_pos"], subset["n_main"],
                     s=sizes, c=REG_COLORS[reg], alpha=0.70, zorder=zo,
                     edgecolors=edge, linewidths=lw,
                     label=f"{reg} (n={len(sub)+len(sub_tp):,})" if edge == "none" else None)

ax_A.scatter([], [], s=6, c="grey", edgecolors="#c0392b", linewidths=0.35,
             label=f"assembly ∪ HiFi confirmed (n={df_called['is_tp'].sum():,})")

ax_A.axvspan(GENE_LO, GENE_HI, color=GENE_COLOR, alpha=0.35, zorder=1)
ax_A.set_yscale("log")
ax_A.set_xlim(WIN_LO - 5, WIN_HI + 5)
ax_A.set_xlabel("T2T consensus position", fontsize=7)
ax_A.set_ylabel("Carriers at VAF ≥ 0.30%", fontsize=7)
ax_A.tick_params(labelsize=6)
ax_A.legend(fontsize=5.5, title="Region", title_fontsize=5.5, loc="upper right",
            handletextpad=0.3, borderpad=0.4, labelspacing=0.25)
ax_A.grid(True, lw=0.3, alpha=0.4)

ax_A2 = ax_A.twinx()
ax_A2.fill_between(depth_df["t2t_pos"], depth_df["p25_dp"], depth_df["p75_dp"],
                   alpha=0.12, color="grey")
ax_A2.plot(depth_df["t2t_pos"], depth_df["median_dp"],
           color="grey", lw=0.8, alpha=0.5)
ax_A2.set_ylabel("Read depth", color="grey", fontsize=6)
ax_A2.tick_params(axis="y", labelcolor="grey", labelsize=5.5)
ax_A2.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v/1e3:.1f}k"))

fig.tight_layout(pad=0.5)
out = OUTDIR / "78a_panel_A.png"
fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(str(out).replace(".png", ".pdf"), bbox_inches="tight")
plt.close()
print(f"Saved: {out}  (11.5 cm × 4 cm, 300 DPI)")
