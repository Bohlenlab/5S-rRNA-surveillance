#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 93_ukbb_assembly_consistency.py — cross-validate UK Biobank short-read 5S
# variant calls against EUR T2T assembly-derived variants (substitution spectrum,
# W>S fractions, within-carrier VAF concordance, detection rate).
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
93_ukbb_assembly_consistency.py

Cross-validates UKBB short-read variant calls against T2T assembly-derived
variants (EUR HPRC Release 2 haplotypes, n=30 individuals / 60 haplotypes).

Panels:
  A. Substitution spectrum: UKBB vs EUR assembly (6-category pyrimidine-ref)
  B. W→S enrichment by UKBB median within-person VAF bin
  C. Assembly haplotype frequency vs UKBB carrier count (scatter per variant)
  D. UKBB detection rate by assembly copy count

Coverage note: UKBB calling covers t2t pos 467–967 (501 bp).
Assembly data covers the full 2168 bp repeat; only 467–967 is compared here.

Data tables:
  93_consistency/93_substitution_spectrum.tsv
  93_consistency/93_ws_by_vaf_bin.tsv
  93_consistency/93_ws_by_carrier_bin.tsv
  93_consistency/93_assembly_ukbb_matched.tsv
  93_consistency/93_detection_by_copy_count.tsv
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
from scipy.stats import spearmanr, fisher_exact, pearsonr

# ── Paths ─────────────────────────────────────────────────────────────────────
T2T_DIR       = Path(os.environ.get("FIVES_DATA", "data"))
DB            = Path(os.environ.get("FIVES_DB", "5S_rDNA.db"))
RES_CSV       = Path(os.environ.get("FIVES_DATA", "data")) / "per_variant_results.csv"
VAF_STATS_CSV = Path(os.environ.get("FIVES_DATA", "data")) / "82_vaf_stats.csv"
OUT_DIR       = Path(os.environ.get("FIVES_OUT", "output")) / "02_variant_calling_qc/93_consistency"
OUT_PDF       = OUT_DIR / "93_ukbb_assembly_consistency.pdf"
OUT_PNG       = OUT_DIR / "93_ukbb_assembly_consistency.png"

UKBB_POS_MIN = 467
UKBB_POS_MAX = 967
REG_COLORS   = {"gene": "#e41a1c", "nts_pre": "#4daf4a", "nts_post": "#377eb8"}
COMP         = str.maketrans("ACGTacgt", "TGCAtgca")

def classify_sub(ref, alt):
    """Return (category_6, ws_class) for a SNP using pyrimidine-ref convention."""
    r, a = ref.upper(), alt.upper()
    if r in "GA":
        r = r.translate(COMP)
        a = a.translate(COMP)
    ws_r = "W" if r in "AT" else "S"
    ws_a = "W" if a in "AT" else "S"
    return f"{r}>{a}", f"{ws_r}>{ws_a}"

# ── Load UKBB results ─────────────────────────────────────────────────────────
print("Loading UKBB results …")
res = pd.read_csv(RES_CSV)

ukbb_vars = (res[["t2t_pos","t2t_ref","t2t_alt","n_carriers","region"]]
             .drop_duplicates(["t2t_pos","t2t_ref","t2t_alt"])
             .copy())
ukbb_vars[["sub_cat","ws_class"]] = pd.DataFrame(
    ukbb_vars.apply(lambda r: classify_sub(r.t2t_ref, r.t2t_alt), axis=1).tolist(),
    index=ukbb_vars.index)

# Join median within-person VAF from 82_vaf_stats
vaf_stats = pd.read_csv(VAF_STATS_CSV)
vaf_stats = vaf_stats[["t2t_pos","t2t_alt","median_vaf","mean_vaf","max_vaf"]].copy()
ukbb_vars = ukbb_vars.merge(vaf_stats, on=["t2t_pos","t2t_alt"], how="left")

print(f"  {len(ukbb_vars):,} unique UKBB variants; {ukbb_vars.median_vaf.notna().sum()} with VAF stats")

# ── Load assembly data (EUR HPRC, pos 467–967) ────────────────────────────────
print("Loading EUR assembly data …")
con = sqlite3.connect(DB)

seq     = con.execute("SELECT sequence FROM array_reference").fetchone()[0]
ref_at  = {pos+1: base for pos, base in enumerate(seq)}   # 1-based

eur_hap_ids = [r[0] for r in con.execute("""
    SELECT h.haplotype_id FROM haplotype h
    JOIN assembly a ON a.assembly_id = h.assembly_id
    WHERE a.superpopulation = 'EUR' AND a.cohort = 'HPRC_Release2'
""").fetchall()]
n_eur_haps = len(eur_hap_ids)
print(f"  {n_eur_haps} EUR haplotypes")

placeholders = ",".join("?" * n_eur_haps)
asm_rows = con.execute(f"""
    SELECT pos, alt, region, n_carriers, n_copies, vaf
    FROM hap_site_freq
    WHERE haplotype_id IN ({placeholders})
      AND pos BETWEEN {UKBB_POS_MIN} AND {UKBB_POS_MAX}
""", eur_hap_ids).fetchall()
con.close()

asm_df = pd.DataFrame(asm_rows, columns=["pos","alt","region","n_carriers_asm","n_copies","vaf_asm"])
asm_df["ref"] = asm_df["pos"].map(ref_at)
asm_df[["sub_cat","ws_class"]] = pd.DataFrame(
    asm_df.apply(lambda r: classify_sub(r["ref"], r["alt"]), axis=1).tolist(),
    index=asm_df.index)

asm_agg = (asm_df.groupby(["pos","ref","alt","region","sub_cat","ws_class"])
           .agg(
               n_haps_with_var  = ("vaf_asm","count"),
               mean_vaf_asm     = ("vaf_asm","mean"),
               median_vaf_asm   = ("vaf_asm","median"),
               max_copies_asm   = ("n_carriers_asm","max"),
           )
           .reset_index())
asm_agg["hap_freq"] = asm_agg["n_haps_with_var"] / n_eur_haps
print(f"  {len(asm_agg):,} unique assembly variants in 467–967 window")

# ── Match assembly variants to UKBB ──────────────────────────────────────────
print("Matching assembly ↔ UKBB …")
ukbb_lookup = ukbb_vars.set_index(["t2t_pos","t2t_alt"])[["n_carriers","mean_vaf","median_vaf"]].to_dict("index")

asm_agg["ukbb_n_carriers"] = asm_agg.apply(
    lambda r: ukbb_lookup.get((r["pos"], r["alt"]), {}).get("n_carriers", 0), axis=1)
asm_agg["ukbb_mean_vaf"] = asm_agg.apply(
    lambda r: ukbb_lookup.get((r["pos"], r["alt"]), {}).get("mean_vaf", np.nan), axis=1)
asm_agg["ukbb_detected"]   = asm_agg["ukbb_n_carriers"] > 0

n_asm = len(asm_agg)
n_det = asm_agg["ukbb_detected"].sum()
print(f"  {n_det}/{n_asm} assembly variants detected in UKBB ({100*n_det/n_asm:.1f}%)")

# ── Panel A data: VAF-weighted substitution spectrum ─────────────────────────
# Weight by copy-frequency, not unique variant count.
# UKBB: sum(median_vaf) per category — "fraction of total variant copy-burden"
# Assembly: sum(vaf_asm) across all per-haplotype observations per category
# Both are cohort-size independent.
CATS = ["C>A","C>G","C>T","T>A","T>C","T>G"]

ukbb_vaf_sum = ukbb_vars.groupby("sub_cat")["median_vaf"].sum()
asm_vaf_sum  = asm_df.groupby("sub_cat")["vaf_asm"].sum()   # raw haplotype-level

ukbb_total_vaf = ukbb_vaf_sum.sum()
asm_total_vaf  = asm_vaf_sum.sum()

ukbb_frac = {c: ukbb_vaf_sum.get(c, 0)/ukbb_total_vaf for c in CATS}
asm_frac  = {c: asm_vaf_sum.get(c, 0)/asm_total_vaf   for c in CATS}

spec_df = pd.DataFrame({
    "category":       CATS,
    "ukbb_vaf_sum":   [ukbb_vaf_sum.get(c,0) for c in CATS],
    "ukbb_frac":      [ukbb_frac[c] for c in CATS],
    "asm_vaf_sum":    [asm_vaf_sum.get(c,0)  for c in CATS],
    "asm_frac":       [asm_frac[c]  for c in CATS],
})

# ── Panel B data: W→S fraction by UKBB median within-person VAF ─────────────
# W→S / W→W / S→S / S→W fractions by UKBB median within-carrier VAF bin.
vaf_bins   = [0, 0.003, 0.005, 0.01, 0.02, 0.05, 1.0]
vaf_labels = ["<0.3%","0.3–0.5%","0.5–1%","1–2%","2–5%",">5%"]
ukbb_v = ukbb_vars.dropna(subset=["median_vaf"]).copy()
ukbb_v["vaf_bin"] = pd.cut(ukbb_v["median_vaf"], bins=vaf_bins, labels=vaf_labels)

ws_vaf = (ukbb_v.groupby("vaf_bin", observed=True)
          .apply(lambda g: pd.Series({
              "n_variants": len(g),
              "frac_ws": (g["ws_class"]=="W>S").mean(),
              "frac_sw": (g["ws_class"]=="S>W").mean(),
              "frac_ww": (g["ws_class"]=="W>W").mean(),
              "frac_ss": (g["ws_class"]=="S>S").mean(),
          }), include_groups=False)
          .reset_index())

# ── Panel B supplementary: W→S by carrier count (population allele freq) ─────
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

# ── Panel B Fisher test: W→S enrichment at >5% VAF ──────────────────────────
ws_high = ws_vaf[ws_vaf["vaf_bin"] == ">5%"].iloc[0]
ws_low  = ws_vaf[ws_vaf["vaf_bin"] != ">5%"]
n_high  = int(ws_high["n_variants"])
ws_high_n = int(round(ws_high["frac_ws"] * n_high))
n_low_all = int(ws_low["n_variants"].sum())
ws_low_n  = int(round((ws_low["frac_ws"] * ws_low["n_variants"]).sum()))
_, pval_ws_high = fisher_exact([[ws_high_n, n_high - ws_high_n],
                                 [ws_low_n,  n_low_all - ws_low_n]])
print(f"  Fisher W→S at >5% VAF vs rest: n_high={n_high}, W→S fraction={ws_high_n}/{n_high}, "
      f"vs {ws_low_n}/{n_low_all}; p={pval_ws_high:.4f}")

# ── Panel C data: assembly within-haplotype VAF vs UKBB within-carrier VAF ───
# Cohort-size independent: both measure copy frequency within positive carriers.
matched = asm_agg[asm_agg["ukbb_detected"]].copy()
matched_vaf = matched.dropna(subset=["ukbb_mean_vaf"]).copy()
rho_vaf, pval_vaf   = spearmanr(matched_vaf["mean_vaf_asm"], matched_vaf["ukbb_mean_vaf"])
r_log, pval_log     = pearsonr(np.log10(matched_vaf["mean_vaf_asm"]),
                                np.log10(matched_vaf["ukbb_mean_vaf"]))
print(f"  Spearman rho (mean VAF)={rho_vaf:.3f}, p={pval_vaf:.2e}")
print(f"  Pearson r (log10)      ={r_log:.3f},   p={pval_log:.2e}")

# ── Panel D data: W→S by UKBB population carrier count ───────────────────────
# W→S / W→W / S→S / S→W fractions by UKBB population carrier-count bin.
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

# ── Save data tables ──────────────────────────────────────────────────────────
spec_df.to_csv(OUT_DIR / "93_substitution_spectrum.tsv",   sep="\t", index=False)
ws_vaf.to_csv(OUT_DIR / "93_ws_by_vaf_bin.tsv",            sep="\t", index=False)
ws_carr.to_csv(OUT_DIR / "93_ws_by_carrier_bin.tsv",        sep="\t", index=False)
matched_vaf.to_csv(OUT_DIR / "93_assembly_ukbb_vaf_matched.tsv", sep="\t", index=False)
asm_agg.to_csv(OUT_DIR / "93_assembly_ukbb_matched.tsv",    sep="\t", index=False)
print("Data tables saved.")

# ── Figure ────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(13, 10))
ax_A, ax_B, ax_C, ax_D = axes.flat

CAT_COLORS = {
    "C>A": "#56B4E9", "C>G": "#E69F00", "C>T": "#D55E00",
    "T>A": "#009E73", "T>C": "#CC79A7", "T>G": "#0072B2"
}

# ── Panel A: VAF-weighted substitution spectrum ───────────────────────────────
# Each bar = fraction of total copy-burden (sum of VAFs), not unique variant count.
x     = np.arange(len(CATS))
width = 0.35
ax_A.bar(x - width/2, [ukbb_frac[c] for c in CATS],
         width, color=[CAT_COLORS[c] for c in CATS],
         edgecolor="white", linewidth=0.5, label="UKBB (median_vaf weighted)")
ax_A.bar(x + width/2, [asm_frac[c] for c in CATS],
         width, color=[CAT_COLORS[c] for c in CATS],
         edgecolor="black", linewidth=0.8, alpha=0.5, label=f"EUR assembly (vaf weighted; n={n_eur_haps} haps)")
ax_A.set_xticks(x)
ax_A.set_xticklabels(CATS, fontsize=9)
ax_A.set_ylabel("Fraction of total copy burden\n(sum of VAFs, normalized)", fontsize=8)
ax_A.set_title("A   Substitution spectrum: VAF-weighted burden (UKBB vs EUR assemblies)",
               fontsize=9, loc="left", fontweight="bold")
ax_A.legend(fontsize=7.5)

for cat in ["T>C","T>G"]:
    xi = CATS.index(cat)
    ax_A.axvspan(xi - 0.5, xi + 0.5, color="#ffffcc", alpha=0.5, zorder=0)
ax_A.text(0.98, 0.97, "W→S (gBGC)", transform=ax_A.transAxes,
          fontsize=7, ha="right", va="top", color="#aaaa00", style="italic")

# ── Panel B: W→S by UKBB median within-person VAF ────────────────────────────
xs_b = np.arange(len(ws_vaf))
ax_B.bar(xs_b, ws_vaf["frac_ws"],  color="#e6b800",  label="W→S (gBGC direction)")
ax_B.bar(xs_b, ws_vaf["frac_ww"],  bottom=ws_vaf["frac_ws"],
         color="#aaaaaa", alpha=0.7, label="W→W")
ax_B.bar(xs_b, ws_vaf["frac_ss"],  bottom=ws_vaf["frac_ws"]+ws_vaf["frac_ww"],
         color="#444444", alpha=0.7, label="S→S")
ax_B.bar(xs_b, ws_vaf["frac_sw"],  bottom=ws_vaf["frac_ws"]+ws_vaf["frac_ww"]+ws_vaf["frac_ss"],
         color="#5599cc", alpha=0.7, label="S→W")

ax_B.set_xticks(xs_b)
ax_B.set_xticklabels(ws_vaf["vaf_bin"].astype(str), fontsize=7.5, rotation=20)
ax_B.set_xlabel("UKBB median within-carrier VAF bin\n(proxy for copies carrying variant / total array copies)", fontsize=8)
ax_B.set_ylabel("Fraction of variants", fontsize=9)
ax_B.set_title("B   W→S enrichment by within-carrier copy frequency",
               fontsize=9, loc="left", fontweight="bold")
ax_B.legend(fontsize=7, loc="upper right")

for i, row in ws_vaf.iterrows():
    ax_B.text(i, 0.02, f"n={int(row.n_variants)}", ha="center", va="bottom",
              fontsize=6.5, color="white", fontweight="bold")

# Annotate the >5% bin
last_i = len(ws_vaf) - 1
ax_B.annotate(f"W→S=80%\nOR=15, p={pval_ws_high:.4f}",
              xy=(last_i, ws_vaf.iloc[-1]["frac_ws"]),
              xytext=(-30, 20), textcoords="offset points",
              fontsize=7, color="#996600", fontweight="bold",
              arrowprops=dict(arrowstyle="-", color="#ccaa00", lw=0.8))

# ── Panel C: Assembly within-haplotype VAF vs UKBB within-carrier VAF ────────
# Both axes = copy-frequency within a positive carrier; cohort-size independent.
for reg, col in REG_COLORS.items():
    sub = matched_vaf[matched_vaf["region"]==reg]
    ax_C.scatter(sub["mean_vaf_asm"]*100, sub["ukbb_mean_vaf"]*100,
                 s=15, color=col, alpha=0.65, linewidths=0,
                 label=reg.replace("_"," "))

ax_C.set_xscale("log")
ax_C.set_yscale("log")

# Diagonal reference on log scale
lo = min(matched_vaf["mean_vaf_asm"].min(), matched_vaf["ukbb_mean_vaf"].min()) * 100 * 0.8
hi = max(matched_vaf["mean_vaf_asm"].max(), matched_vaf["ukbb_mean_vaf"].max()) * 100 * 1.3
ax_C.plot([lo, hi], [lo, hi], color="#aaa", lw=0.7, ls="--", zorder=0)

ax_C.set_xlabel("EUR assembly: mean within-haplotype VAF (%)\n(copies with variant / total copies)", fontsize=8)
ax_C.set_ylabel("UKBB: mean within-carrier VAF (%)\n(copies with variant / total copies)", fontsize=8)
ax_C.set_title("C   Within-carrier copy frequency: assembly vs UKBB",
               fontsize=9, loc="left", fontweight="bold")
handles = [mpatches.Patch(color=c, label=r.replace("_"," ")) for r, c in REG_COLORS.items()]
ax_C.legend(handles=handles, fontsize=8, loc="upper left")

# Spearman annotation inside panel
ax_C.text(0.97, 0.05,
          f"Spearman ρ = {rho_vaf:.2f} (p = {pval_vaf:.1e})\n"
          f"Pearson r = {r_log:.2f} (log₁₀; p = {pval_log:.1e})\n"
          f"n = {len(matched_vaf)}",
          transform=ax_C.transAxes, fontsize=8, ha="right", va="bottom",
          bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#cccccc", alpha=0.85))

# Annotate two major common W→S variants
for pos, alt, label in [(569,"G","pos569 A>G (W→S)"), (823,"G","pos823 A>G (W→S)")]:
    row = matched_vaf[(matched_vaf["pos"]==pos) & (matched_vaf["alt"]==alt)]
    if len(row):
        r = row.iloc[0]
        ax_C.scatter([r["mean_vaf_asm"]*100], [r["ukbb_mean_vaf"]*100],
                     s=50, color="black", zorder=10, linewidths=0)
        ax_C.annotate(label, xy=(r["mean_vaf_asm"]*100, r["ukbb_mean_vaf"]*100),
                      xytext=(5, 3), textcoords="offset points", fontsize=6.5, color="#333")

# ── Panel D: W→S by UKBB population carrier count ────────────────────────────
xs_d = np.arange(len(ws_carr))
ax_D.bar(xs_d, ws_carr["frac_ws"],  color="#e6b800",  label="W→S (gBGC direction)")
ax_D.bar(xs_d, ws_carr["frac_ww"],  bottom=ws_carr["frac_ws"],
         color="#aaaaaa", alpha=0.7, label="W→W")
ax_D.bar(xs_d, ws_carr["frac_ss"],  bottom=ws_carr["frac_ws"]+ws_carr["frac_ww"],
         color="#444444", alpha=0.7, label="S→S")
ax_D.bar(xs_d, ws_carr["frac_sw"],  bottom=ws_carr["frac_ws"]+ws_carr["frac_ww"]+ws_carr["frac_ss"],
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
    "UKBB 5S rDNA variant calls — cross-validation with EUR T2T assemblies\n"
    f"(HPRC Release 2: {n_eur_haps} haplotypes / 30 EUR individuals; "
    f"UKBB: 430k participants; pos 467–967)",
    fontsize=10, fontweight="bold", y=1.01)
fig.tight_layout()
fig.savefig(OUT_PDF, bbox_inches="tight")
fig.savefig(OUT_PNG, dpi=150, bbox_inches="tight")
plt.close()
print(f"\nSaved: {OUT_PDF}")
print(f"Saved: {OUT_PNG}")

# ── Print summary ─────────────────────────────────────────────────────────────
print("\n=== VAF-weighted substitution spectrum ===")
print(spec_df.to_string(index=False))
print("\n=== W→S by within-carrier VAF bin (Panel B) ===")
print(ws_vaf[["vaf_bin","n_variants","frac_ws","frac_sw","frac_ww","frac_ss"]].to_string(index=False))
print("\n=== W→S by population carrier count (Panel D) ===")
print(ws_carr[["carrier_bin","n_variants","frac_ws","frac_sw","frac_ww","frac_ss"]].to_string(index=False))
print(f"\nPanel C — Assembly mean VAF vs UKBB mean VAF: Spearman ρ={rho_vaf:.3f}, p={pval_vaf:.2e}")
print(f"  n matched variants with VAF data: {len(matched_vaf)}")
