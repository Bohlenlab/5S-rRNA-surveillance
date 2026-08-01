#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 93b_ukbb_assembly_consistency_all.py — cross-validate UK Biobank short-read 5S
# variant calls against all HPRC Release 2 assembly haplotypes, using mean VAF.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
93b_ukbb_assembly_consistency_all.py

Four-panel cross-validation using:
  - All HPRC Release 2 haplotypes (372 haps, all superpopulations)
  - Mean VAF for both assembly and UKBB (instead of median)

Panels:
  A. VAF-weighted substitution spectrum (all assembly haps vs UKBB)
  B. W→S enrichment by UKBB mean within-carrier VAF bin
  C. Assembly mean within-haplotype VAF vs UKBB mean within-carrier VAF (log-log)
  D. W→S by UKBB population carrier count

Coverage: UKBB pos 467–967 (501 bp).
"""

import os
import sqlite3
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, fisher_exact

# ── Paths ─────────────────────────────────────────────────────────────────────
T2T_DIR       = Path(os.environ.get("FIVES_DATA", "data"))
DB            = Path(os.environ.get("FIVES_DB", "5S_rDNA.db"))
RES_CSV       = Path(os.environ.get("FIVES_DATA", "data")) / "per_variant_results.csv"
VAF_STATS_CSV = Path(os.environ.get("FIVES_DATA", "data")) / "82_vaf_stats.csv"
OUT_DIR       = Path(os.environ.get("FIVES_OUT", "output")) / "02_variant_calling_qc/93_consistency"
OUT_PDF       = OUT_DIR / "93b_ukbb_assembly_consistency_all.pdf"
OUT_PNG       = OUT_DIR / "93b_ukbb_assembly_consistency_all.png"

UKBB_POS_MIN = 467
UKBB_POS_MAX = 967
REG_COLORS   = {"gene": "#e41a1c", "nts_pre": "#4daf4a", "nts_post": "#377eb8"}
COMP         = str.maketrans("ACGTacgt", "TGCAtgca")

def classify_sub(ref, alt):
    r, a = ref.upper(), alt.upper()
    if r in "GA":
        r = r.translate(COMP)
        a = a.translate(COMP)
    ws_r = "W" if r in "AT" else "S"
    ws_a = "W" if a in "AT" else "S"
    return f"{r}>{a}", f"{ws_r}>{ws_a}"

# ── Load UKBB ─────────────────────────────────────────────────────────────────
print("Loading UKBB results …")
res = pd.read_csv(RES_CSV)
vaf_stats = pd.read_csv(VAF_STATS_CSV)

ukbb_vars = (res[["t2t_pos","t2t_ref","t2t_alt","n_carriers","region"]]
             .drop_duplicates(["t2t_pos","t2t_ref","t2t_alt"])
             .copy())
ukbb_vars = ukbb_vars.merge(
    vaf_stats[["t2t_pos","t2t_alt","mean_vaf","median_vaf"]],
    on=["t2t_pos","t2t_alt"], how="left")
ukbb_vars[["sub_cat","ws_class"]] = pd.DataFrame(
    ukbb_vars.apply(lambda r: classify_sub(r.t2t_ref, r.t2t_alt), axis=1).tolist(),
    index=ukbb_vars.index)
print(f"  {len(ukbb_vars):,} unique UKBB variants")

# ── Load ALL assembly haplotypes ──────────────────────────────────────────────
print("Loading ALL HPRC assembly data …")
con = sqlite3.connect(DB)

seq    = con.execute("SELECT sequence FROM array_reference").fetchone()[0]
ref_at = {pos+1: base for pos, base in enumerate(seq)}

all_hap_ids = [r[0] for r in con.execute("""
    SELECT h.haplotype_id FROM haplotype h
    JOIN assembly a ON a.assembly_id = h.assembly_id
    WHERE a.cohort = 'HPRC_Release2'
""").fetchall()]
n_all_haps = len(all_hap_ids)
print(f"  {n_all_haps} haplotypes (all superpopulations)")

ph = ",".join("?" * n_all_haps)
asm_rows = con.execute(f"""
    SELECT pos, alt, region, n_carriers, n_copies, vaf
    FROM hap_site_freq
    WHERE haplotype_id IN ({ph})
      AND pos BETWEEN {UKBB_POS_MIN} AND {UKBB_POS_MAX}
""", all_hap_ids).fetchall()
con.close()

asm_df = pd.DataFrame(asm_rows, columns=["pos","alt","region","n_carriers_asm","n_copies","vaf_asm"])
asm_df["ref"] = asm_df["pos"].map(ref_at)
asm_df[["sub_cat","ws_class"]] = pd.DataFrame(
    asm_df.apply(lambda r: classify_sub(r["ref"], r["alt"]), axis=1).tolist(),
    index=asm_df.index)

asm_agg = (asm_df.groupby(["pos","ref","alt","region","sub_cat","ws_class"])
           .agg(
               n_haps_with_var = ("vaf_asm","count"),
               mean_vaf_asm    = ("vaf_asm","mean"),
               median_vaf_asm  = ("vaf_asm","median"),
               max_copies_asm  = ("n_carriers_asm","max"),
           )
           .reset_index())
asm_agg["hap_freq"] = asm_agg["n_haps_with_var"] / n_all_haps
print(f"  {len(asm_agg):,} unique assembly variants in 467–967 window")

# ── Match assembly → UKBB using mean_vaf ─────────────────────────────────────
print("Matching …")
ukbb_lookup = ukbb_vars.set_index(["t2t_pos","t2t_alt"])[["n_carriers","mean_vaf"]].to_dict("index")

asm_agg["ukbb_n_carriers"] = asm_agg.apply(
    lambda r: ukbb_lookup.get((r["pos"], r["alt"]), {}).get("n_carriers", 0), axis=1)
asm_agg["ukbb_mean_vaf"]   = asm_agg.apply(
    lambda r: ukbb_lookup.get((r["pos"], r["alt"]), {}).get("mean_vaf", np.nan), axis=1)
asm_agg["ukbb_detected"]   = asm_agg["ukbb_n_carriers"] > 0

n_asm = len(asm_agg)
n_det = asm_agg["ukbb_detected"].sum()
print(f"  {n_det}/{n_asm} assembly variants detected in UKBB ({100*n_det/n_asm:.1f}%)")

# ── Panel A: VAF-weighted substitution spectrum ───────────────────────────────
CATS = ["C>A","C>G","C>T","T>A","T>C","T>G"]

ukbb_vaf_sum = ukbb_vars.groupby("sub_cat")["mean_vaf"].sum()
asm_vaf_sum  = asm_df.groupby("sub_cat")["vaf_asm"].sum()

ukbb_frac = {c: ukbb_vaf_sum.get(c,0)/ukbb_vaf_sum.sum() for c in CATS}
asm_frac  = {c: asm_vaf_sum.get(c,0)/asm_vaf_sum.sum()   for c in CATS}

spec_df = pd.DataFrame({
    "category":    CATS,
    "ukbb_vaf_sum": [ukbb_vaf_sum.get(c,0) for c in CATS],
    "ukbb_frac":   [ukbb_frac[c] for c in CATS],
    "asm_vaf_sum": [asm_vaf_sum.get(c,0) for c in CATS],
    "asm_frac":    [asm_frac[c] for c in CATS],
})

# ── Panel B: W→S by UKBB mean within-carrier VAF bin ─────────────────────────
vaf_bins   = [0, 0.003, 0.005, 0.01, 0.02, 0.05, 1.0]
vaf_labels = ["<0.3%","0.3–0.5%","0.5–1%","1–2%","2–5%",">5%"]
ukbb_v = ukbb_vars.dropna(subset=["mean_vaf"]).copy()
ukbb_v["vaf_bin"] = pd.cut(ukbb_v["mean_vaf"], bins=vaf_bins, labels=vaf_labels)

ws_vaf = (ukbb_v.groupby("vaf_bin", observed=True)
          .apply(lambda g: pd.Series({
              "n_variants": len(g),
              "frac_ws": (g["ws_class"]=="W>S").mean(),
              "frac_sw": (g["ws_class"]=="S>W").mean(),
              "frac_ww": (g["ws_class"]=="W>W").mean(),
              "frac_ss": (g["ws_class"]=="S>S").mean(),
          }), include_groups=False)
          .reset_index())

# Fisher: W→S at >5% vs rest
ws_h = ws_vaf[ws_vaf["vaf_bin"]==">5%"].iloc[0]
n_h  = int(ws_h["n_variants"]); ws_hn = int(round(ws_h["frac_ws"]*n_h))
ws_r_df = ws_vaf[ws_vaf["vaf_bin"]!=">5%"]
n_r = int(ws_r_df["n_variants"].sum()); ws_rn = int(round((ws_r_df["frac_ws"]*ws_r_df["n_variants"]).sum()))
_, pval_ws_high = fisher_exact([[ws_hn, n_h-ws_hn],[ws_rn, n_r-ws_rn]])
print(f"  Fisher W→S >5% mean VAF vs rest: {ws_hn}/{n_h} vs {ws_rn}/{n_r}, p={pval_ws_high:.4f}")

# ── Panel C: assembly mean VAF vs UKBB mean VAF ───────────────────────────────
matched_vaf = asm_agg[asm_agg["ukbb_detected"]].dropna(subset=["ukbb_mean_vaf"]).copy()
rho_vaf, pval_vaf = spearmanr(matched_vaf["mean_vaf_asm"], matched_vaf["ukbb_mean_vaf"])
print(f"  Spearman ρ (asm mean VAF vs ukbb mean VAF) = {rho_vaf:.3f}, p={pval_vaf:.2e}, n={len(matched_vaf)}")

# ── Panel D: W→S by UKBB carrier count ───────────────────────────────────────
carr_bins   = [0, 100, 500, 2000, 10000, 500000]
carr_labels = ["<100","100–500","500–2k","2k–10k",">10k"]
ukbb_vars["carrier_bin"] = pd.cut(ukbb_vars["n_carriers"], bins=carr_bins, labels=carr_labels)

ws_carr = (ukbb_vars.groupby("carrier_bin", observed=True)
           .apply(lambda g: pd.Series({
               "n_variants": len(g),
               "frac_ws": (g["ws_class"]=="W>S").mean(),
               "frac_sw": (g["ws_class"]=="S>W").mean(),
               "frac_ww": (g["ws_class"]=="W>W").mean(),
               "frac_ss": (g["ws_class"]=="S>S").mean(),
           }), include_groups=False)
           .reset_index())

# ── Save tables ───────────────────────────────────────────────────────────────
spec_df.to_csv(OUT_DIR / "93b_substitution_spectrum_all.tsv",     sep="\t", index=False)
ws_vaf.to_csv(OUT_DIR / "93b_ws_by_vaf_bin_all.tsv",              sep="\t", index=False)
ws_carr.to_csv(OUT_DIR / "93b_ws_by_carrier_bin_all.tsv",          sep="\t", index=False)
matched_vaf.to_csv(OUT_DIR / "93b_assembly_ukbb_vaf_matched_all.tsv", sep="\t", index=False)
print("Tables saved.")

# ── Figure ────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(13, 10))
ax_A, ax_B, ax_C, ax_D = axes.flat

CAT_COLORS = {
    "C>A": "#56B4E9", "C>G": "#E69F00", "C>T": "#D55E00",
    "T>A": "#009E73", "T>C": "#CC79A7", "T>G": "#0072B2"
}

# ── Panel A ───────────────────────────────────────────────────────────────────
x, width = np.arange(len(CATS)), 0.35
ax_A.bar(x - width/2, [ukbb_frac[c] for c in CATS], width,
         color=[CAT_COLORS[c] for c in CATS], edgecolor="white", linewidth=0.5,
         label="UKBB (mean_vaf weighted)")
ax_A.bar(x + width/2, [asm_frac[c] for c in CATS], width,
         color=[CAT_COLORS[c] for c in CATS], edgecolor="black", linewidth=0.8,
         alpha=0.5, label=f"All assembly (vaf weighted; n={n_all_haps} haps)")
ax_A.set_xticks(x)
ax_A.set_xticklabels(CATS, fontsize=9)
ax_A.set_ylabel("Fraction of total copy burden\n(sum of VAFs, normalized)", fontsize=8)
ax_A.set_title("A   Substitution spectrum: VAF-weighted burden (UKBB vs all assemblies)",
               fontsize=9, loc="left", fontweight="bold")
ax_A.legend(fontsize=7.5)
for cat in ["T>C","T>G"]:
    ax_A.axvspan(CATS.index(cat) - 0.5, CATS.index(cat) + 0.5, color="#ffffcc", alpha=0.5, zorder=0)
ax_A.text(0.98, 0.97, "W→S (gBGC)", transform=ax_A.transAxes,
          fontsize=7, ha="right", va="top", color="#aaaa00", style="italic")

# ── Panel B ───────────────────────────────────────────────────────────────────
xs_b = np.arange(len(ws_vaf))
ax_B.bar(xs_b, ws_vaf["frac_ws"], color="#e6b800", label="W→S (gBGC direction)")
ax_B.bar(xs_b, ws_vaf["frac_ww"], bottom=ws_vaf["frac_ws"],
         color="#aaaaaa", alpha=0.7, label="W→W")
ax_B.bar(xs_b, ws_vaf["frac_ss"], bottom=ws_vaf["frac_ws"]+ws_vaf["frac_ww"],
         color="#444444", alpha=0.7, label="S→S")
ax_B.bar(xs_b, ws_vaf["frac_sw"], bottom=ws_vaf["frac_ws"]+ws_vaf["frac_ww"]+ws_vaf["frac_ss"],
         color="#5599cc", alpha=0.7, label="S→W")
ax_B.set_xticks(xs_b)
ax_B.set_xticklabels(ws_vaf["vaf_bin"].astype(str), fontsize=7.5, rotation=20)
ax_B.set_xlabel("UKBB mean within-carrier VAF bin\n(proxy: copies with variant / total array copies)", fontsize=8)
ax_B.set_ylabel("Fraction of variants", fontsize=9)
ax_B.set_title("B   W→S enrichment by within-carrier copy frequency",
               fontsize=9, loc="left", fontweight="bold")
ax_B.legend(fontsize=7, loc="upper right")
for i, row in ws_vaf.iterrows():
    ax_B.text(i, 0.02, f"n={int(row.n_variants)}", ha="center", va="bottom",
              fontsize=6.5, color="white", fontweight="bold")
ax_B.annotate(f"W→S={ws_hn/n_h*100:.0f}%\np={pval_ws_high:.4f}",
              xy=(len(ws_vaf)-1, ws_vaf.iloc[-1]["frac_ws"]),
              xytext=(-35, 18), textcoords="offset points",
              fontsize=7, color="#996600", fontweight="bold",
              arrowprops=dict(arrowstyle="-", color="#ccaa00", lw=0.8))

# ── Panel C ───────────────────────────────────────────────────────────────────
for reg, col in REG_COLORS.items():
    sub = matched_vaf[matched_vaf["region"]==reg]
    ax_C.scatter(sub["mean_vaf_asm"]*100, sub["ukbb_mean_vaf"]*100,
                 s=10, color=col, alpha=0.55, linewidths=0,
                 label=reg.replace("_"," "))

ax_C.set_xscale("log")
ax_C.set_yscale("log")
lo = min(matched_vaf["mean_vaf_asm"].min(), matched_vaf["ukbb_mean_vaf"].min()) * 100 * 0.8
hi = max(matched_vaf["mean_vaf_asm"].max(), matched_vaf["ukbb_mean_vaf"].max()) * 100 * 1.3
ax_C.plot([lo, hi], [lo, hi], color="#aaa", lw=0.7, ls="--", zorder=0)
ax_C.set_xlabel("All assemblies: mean within-haplotype VAF (%)\n(copies with variant / total copies)", fontsize=8)
ax_C.set_ylabel("UKBB: mean within-carrier VAF (%)\n(copies with variant / total copies)", fontsize=8)
ax_C.set_title("C   Within-carrier copy frequency: all assemblies vs UKBB",
               fontsize=9, loc="left", fontweight="bold")
handles = [mpatches.Patch(color=c, label=r.replace("_"," ")) for r, c in REG_COLORS.items()]
ax_C.legend(handles=handles, fontsize=8, loc="upper left")
ax_C.text(0.97, 0.05,
          f"Spearman ρ = {rho_vaf:.2f}\np = {pval_vaf:.1e}\nn = {len(matched_vaf)}",
          transform=ax_C.transAxes, fontsize=8.5, ha="right", va="bottom",
          bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#cccccc", alpha=0.85))

# Annotate two major common W→S variants
for pos, alt, label in [(569,"G","pos569 A>G (W→S)"), (823,"G","pos823 A>G (W→S)")]:
    row = matched_vaf[(matched_vaf["pos"]==pos) & (matched_vaf["alt"]==alt)]
    if len(row):
        r = row.iloc[0]
        ax_C.scatter([r["mean_vaf_asm"]*100], [r["ukbb_mean_vaf"]*100],
                     s=50, color="black", zorder=10, linewidths=0)
        ax_C.annotate(label, xy=(r["mean_vaf_asm"]*100, r["ukbb_mean_vaf"]*100),
                      xytext=(5,3), textcoords="offset points", fontsize=6.5, color="#333")

# ── Panel D ───────────────────────────────────────────────────────────────────
xs_d = np.arange(len(ws_carr))
ax_D.bar(xs_d, ws_carr["frac_ws"], color="#e6b800", label="W→S (gBGC direction)")
ax_D.bar(xs_d, ws_carr["frac_ww"], bottom=ws_carr["frac_ws"],
         color="#aaaaaa", alpha=0.7, label="W→W")
ax_D.bar(xs_d, ws_carr["frac_ss"], bottom=ws_carr["frac_ws"]+ws_carr["frac_ww"],
         color="#444444", alpha=0.7, label="S→S")
ax_D.bar(xs_d, ws_carr["frac_sw"], bottom=ws_carr["frac_ws"]+ws_carr["frac_ww"]+ws_carr["frac_ss"],
         color="#5599cc", alpha=0.7, label="S→W (CpG deamination)")
ax_D.set_xticks(xs_d)
ax_D.set_xticklabels(ws_carr["carrier_bin"].astype(str), fontsize=8, rotation=20)
ax_D.set_xlabel("UKBB carrier count bin\n(allele frequency in population)", fontsize=8)
ax_D.set_ylabel("Fraction of variants", fontsize=9)
ax_D.set_title("D   W→S enrichment by population allele frequency",
               fontsize=9, loc="left", fontweight="bold")
ax_D.legend(fontsize=7, loc="upper right")
for i, row in ws_carr.iterrows():
    ax_D.text(i, 0.02, f"n={int(row.n_variants)}", ha="center", va="bottom",
              fontsize=6.5, color="white", fontweight="bold")

# ── Finalise ─────────────────────────────────────────────────────────────────
fig.suptitle(
    "UKBB 5S rDNA — cross-validation with ALL HPRC assemblies (mean VAF)\n"
    f"(HPRC Release 2: {n_all_haps} haplotypes, all superpopulations; "
    f"UKBB: 430k participants; pos 467–967)",
    fontsize=10, fontweight="bold", y=1.01)
fig.tight_layout()
fig.savefig(OUT_PDF, bbox_inches="tight")
fig.savefig(OUT_PNG, dpi=150, bbox_inches="tight")
plt.close()
print(f"\nSaved: {OUT_PDF}")
print(f"Saved: {OUT_PNG}")

print("\n=== VAF-weighted spectrum (all assemblies) ===")
print(spec_df.to_string(index=False))
print(f"\nW→S assembly copy burden: {asm_frac['T>C']+asm_frac['T>G']:.1%}")
print(f"W→S UKBB copy burden:     {ukbb_frac['T>C']+ukbb_frac['T>G']:.1%}")
