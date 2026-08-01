#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 26_array_maps.py — draws per-haplotype 5S rDNA array maps plotting SNV and indel positions along the 2168 bp repeat unit against copy number.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
26_array_maps.py

Per-haplotype 5S array maps from the consensus_t2t variant data.

Each map plots copy number (y) against position in the 2168 bp repeat unit (x):
  - SNVs              → dark dots (colored by region)
  - Indels (ins/del)  → yellow stars
  - SNVs masked near indels → light grey dots (flagged as possible artifact)
  - Region bands (NTS-pre / gene / NTS-post) shaded in background

Usage:
  python3 26_array_maps.py --sample HG00438            # both haplotypes
  python3 26_array_maps.py --sample HG00438 --hap hap1
  python3 26_array_maps.py --all                       # every HPRC haplotype
  python3 26_array_maps.py --all --jobs 6
"""

import os
import sys, sqlite3
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

DB_PATH = Path(os.environ.get("FIVES_DB", "5S_rDNA.db"))
OUT_DIR = Path(os.environ.get("FIVES_OUT", "output")) / "figures" / "08_array_maps" / "per_haplotype_maps"
OUT_DIR.mkdir(parents=True, exist_ok=True)

T2T_LEN = 2168
GENE_START, GENE_END = 630, 748

REGION_BANDS = [
    ("NTS-pre",   1,         629,  "#dce6f1"),
    ("5S gene",   630,       748,  "#ffe7b3"),
    ("NTS-post",  749,       2168, "#f5dcdc"),
]
REGION_DOT = {"nts_pre": "#2c5d8f", "gene": "#cc7a00", "nts_post": "#a02020"}


def parse_args():
    a = sys.argv[1:]
    def opt(flag):
        if flag in a:
            i = a.index(flag)
            if i + 1 < len(a): return a[i+1]
        return None
    return {
        "sample": opt("--sample"),
        "hap":    opt("--hap"),
        "all":    "--all" in a,
        "jobs":   int(opt("--jobs") or 4),
    }


def get_haplotype_list(sample=None, hap=None, all_hprc=False):
    con = sqlite3.connect(DB_PATH)
    where = "a.cohort IN ('HPRC_Year1','HPRC_Release2','CHM13','HG002_GIAB')"
    params = []
    if sample:
        where += " AND a.sample_id = ?"; params.append(sample)
    if hap:
        where += " AND h.hap_label = ?"; params.append(hap)
    rows = con.execute(f"""
        SELECT a.sample_id, a.cohort, h.hap_label, h.haplotype_id, h.n_copies
        FROM assembly a JOIN haplotype h USING(assembly_id)
        WHERE {where}
        ORDER BY a.sample_id, h.hap_label
    """, params).fetchall()
    con.close()
    return rows


def make_array_map(task):
    sample_id, cohort, hap_label, hap_id, n_copies = task

    con = sqlite3.connect(DB_PATH)
    # All consensus_t2t variants for this haplotype (interior + border)
    rows = con.execute("""
        SELECT c.copy_number, c.border_note, v.consensus_pos, v.region,
               v.var_type, v.masked
        FROM variant v JOIN copy c USING(copy_id)
        WHERE c.haplotype_id = ? AND v.alignment_source = 'consensus_t2t'
    """, (hap_id,)).fetchall()
    # copy count / numbering
    copy_rows = con.execute("""
        SELECT copy_number, border_note FROM copy WHERE haplotype_id = ?
        ORDER BY copy_number
    """, (hap_id,)).fetchall()
    con.close()

    if not copy_rows:
        return f"  {sample_id} {hap_label}: no copies"

    max_copy = max(cn for cn, _ in copy_rows)

    # Partition variant points
    snp_x, snp_y, snp_c = [], [], []
    msk_x, msk_y = [], []
    # Collapse consecutive same-type indel columns into single events per copy.
    # A multi-bp indel spans several alignment columns → one star at its start.
    indel_by_copy = {}   # copy_number → list of (pos, vtype)
    for cn, border, pos, region, vtype, masked in rows:
        if vtype == "snp":
            if masked:
                msk_x.append(pos); msk_y.append(cn)
            else:
                snp_x.append(pos); snp_y.append(cn)
                snp_c.append(REGION_DOT.get(region, "#444"))
        else:  # ins or del
            indel_by_copy.setdefault(cn, []).append((pos, vtype))

    ind_x, ind_y = [], []
    for cn, plist in indel_by_copy.items():
        plist.sort()
        prev_pos, prev_vt = None, None
        for pos, vt in plist:
            # new event if type changes or gap > 1 bp from previous column
            if prev_pos is None or vt != prev_vt or pos - prev_pos > 1:
                ind_x.append(pos); ind_y.append(cn)
            prev_pos, prev_vt = pos, vt

    # ── figure ────────────────────────────────────────────────────────────────
    height = max(4, min(0.10 * max_copy + 1.5, 26))
    fig, ax = plt.subplots(figsize=(13, height))

    # region bands
    for name, lo, hi, col in REGION_BANDS:
        ax.axvspan(lo, hi, color=col, zorder=0)
        ax.text((lo + hi) / 2, max_copy + max_copy*0.02 + 0.6, name,
                ha="center", va="bottom", fontsize=8, color="#333")

    # gene boundary lines
    ax.axvline(GENE_START, color="#cc7a00", lw=0.6, ls=":", alpha=0.7, zorder=1)
    ax.axvline(GENE_END,   color="#cc7a00", lw=0.6, ls=":", alpha=0.7, zorder=1)

    # masked SNVs (light grey, behind)
    if msk_x:
        ax.scatter(msk_x, msk_y, s=7, c="#bbbbbb", marker="o",
                   linewidths=0, zorder=2, label="SNV near indel (masked)")
    # SNVs (region-colored dots)
    if snp_x:
        ax.scatter(snp_x, snp_y, s=9, c=snp_c, marker="o",
                   linewidths=0, alpha=0.85, zorder=3)
    # Indels (yellow stars)
    if ind_x:
        ax.scatter(ind_x, ind_y, s=55, c="#ffd500", marker="*",
                   edgecolors="#7a6200", linewidths=0.4, zorder=4,
                   label="indel (ins/del)")

    ax.set_xlim(0, T2T_LEN)
    ax.set_ylim(0.5, max_copy + max_copy*0.05 + 1.5)
    ax.invert_yaxis()  # copy 1 (5' end) at top
    ax.set_xlabel("Position in 2168 bp repeat unit (consensus coordinates)")
    ax.set_ylabel("Copy number (5′ → 3′)")
    ax.set_title(f"{sample_id} · {hap_label} · {cohort}  —  {max_copy} copies "
                 f"(vs population consensus)", fontsize=11, fontweight="bold")

    # legend (region dots + indel + masked)
    handles = [
        mpatches.Patch(color=REGION_DOT["nts_pre"],  label="SNV NTS-pre"),
        mpatches.Patch(color=REGION_DOT["gene"],     label="SNV gene"),
        mpatches.Patch(color=REGION_DOT["nts_post"], label="SNV NTS-post"),
        plt.Line2D([], [], marker="*", color="#ffd500", markeredgecolor="#7a6200",
                   linestyle="", markersize=11, label="indel"),
        plt.Line2D([], [], marker="o", color="#bbbbbb", linestyle="",
                   markersize=6, label="masked SNV (near indel)"),
    ]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.005, 1.0),
              fontsize=8, frameon=True, borderaxespad=0)

    fig.tight_layout()
    out = OUT_DIR / f"{sample_id}_{hap_label}_array_map.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return f"  {sample_id} {hap_label}: {max_copy} copies, {len(snp_x)} SNVs, {len(ind_x)} indels → {out.name}"


def main():
    args = parse_args()
    if not (args["sample"] or args["all"]):
        print("Specify --sample SAMPLE [--hap hapX] or --all")
        return

    haps = get_haplotype_list(args["sample"], args["hap"], args["all"])
    print(f"Generating {len(haps)} array map(s)...")

    if args["jobs"] == 1 or len(haps) <= 2:
        for h in haps:
            print(make_array_map(h))
    else:
        with ProcessPoolExecutor(max_workers=args["jobs"]) as pool:
            for res in as_completed([pool.submit(make_array_map, h) for h in haps]):
                print(res.result())

    print(f"\nSaved to {OUT_DIR}/")


if __name__ == "__main__":
    main()
