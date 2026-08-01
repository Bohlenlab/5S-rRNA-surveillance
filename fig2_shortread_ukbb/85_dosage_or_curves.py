#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 85_dosage_or_curves.py — empirical dose-response odds-ratio curves relating
# 5S variant copy-number bins to disease odds for FDR-significant associations.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
85_dosage_or_curves.py

For the top FDR-significant variant-phenotype associations, show the empirical
dose-response relationship between variant copy-number and disease odds.

Method:
  - Carriers binned by est_copies = VAF × 160 (5 log-spaced bins)
  - Non-carriers (cohort members with no call for that variant) = reference
  - OR and 95% CI computed from 2×2 tables (Haldane-corrected)
  - Covariates NOT adjusted here — this is purely descriptive/visual
"""

import os
import sqlite3
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd

# ── Paths ──────────────────────────────────────────────────────────────────────
_SRV = Path(os.environ.get("FIVES_DATA", "data"))
ASSOC_FILE   = _SRV / "results_500k/association/association_input.tsv"
COHORT_DB    = Path(os.environ.get("FIVES_ICD10_DB", "cohort_icd10.db"))
RESULTS_CSV  = _SRV / "81_results/per_variant_results.csv"
OUT_DIR      = Path(os.environ.get("FIVES_OUT", "output")) / "81_results"
FIG_DIR      = Path(os.environ.get("FIVES_OUT", "output")) / "02_variant_calling_qc"

VAF_THRESH    = 0.003
EST_MULT      = 160
TOP_N         = 9       # number of hits to plot when PLOT_ALL=False
PLOT_ALL      = True     # plot every FDR<0.05 association (one panel per variant×phenotype)
MIN_BIN_N     = 5       # minimum carriers per bin to show point

REG_COLORS = {"gene": "#e41a1c", "nts_pre": "#4daf4a", "nts_post": "#377eb8"}

# ICD-10 three-character category titles (WHO ICD-10; verify against UKBB coding19 if needed)
ICD10_DESC = {
    "D09": "Carcinoma in situ of other and unspecified sites",
    "M96": "Postprocedural musculoskeletal disorders, NEC",
    "F42": "Obsessive-compulsive disorder",
    "X50": "Overexertion and strenuous or repetitive movements",
    "M34": "Systemic sclerosis [scleroderma]",
    "I36": "Nonrheumatic tricuspid valve disorders",
    "Y46": "Antiepileptic and antiparkinsonism drugs, adverse effect",
    "B37": "Candidiasis",
    "K28": "Gastrojejunal ulcer",
    "C73": "Malignant neoplasm of thyroid gland",
    "Y59": "Other/unspecified vaccines & biologicals, adverse effect",
    "N31": "Neuromuscular dysfunction of bladder, NEC",
    "Z46": "Encounter for fitting/adjustment of other devices",
    "M80": "Osteoporosis with current pathological fracture",
    "B02": "Zoster [herpes zoster]",
    "F03": "Unspecified dementia",
    "I07": "Rheumatic tricuspid valve diseases",
    "I86": "Varicose veins of other sites",
    "J03": "Acute tonsillitis",
    "J12": "Viral pneumonia, not elsewhere classified",
    "J15": "Bacterial pneumonia, not elsewhere classified",
    "N60": "Benign mammary dysplasia",
    "O67": "Labour/delivery complicated by intrapartum haemorrhage, NEC",
    "Q54": "Hypospadias",
    "R50": "Fever of other and unknown origin",
    "R68": "Other general symptoms and signs",
    "R80": "Isolated proteinuria",
    "T42": "Poisoning by antiepileptic/sedative-hypnotic/antiparkinsonism drugs",
    "T90": "Sequelae of injuries of head",
    "X61": "Intentional self-poisoning by antiepileptic/sedative/psychotropic drugs",
    "Y47": "Sedative-hypnotic and antianxiety drugs, adverse effect",
    "Z04": "Encounter for examination/observation for other reasons",
    "Z81": "Family history of mental and behavioural disorders",
}

# 5S gene-body functional scores (deep-mutational expr/incorp; keyed by T2T position)
FUNC5S = {
    719:"expr=1.05 incorp=1.24 [tolerated]",  648:"expr=0.39 incorp=0.81 [tolerated]",
    694:"expr=0.00 incorp=NA [EXPRESSION-NULL]", 715:"expr=0.17 incorp=2.40 [low-expr]",
    700:"expr=0.00 incorp=NA [EXPRESSION-NULL]", 704:"expr=4.41 incorp=1.39 [tolerated]",
    737:"expr=1.06 incorp=1.72 [tolerated]",  686:"expr=0.00 incorp=NA [EXPRESSION-NULL]",
    745:"expr=1.09 incorp=0.95 [tolerated]",  643:"expr=0.53 incorp=2.31 [tolerated]",
    664:"expr=0.86 incorp=0.38 [low-incorp]", 711:"expr=0.37 incorp=1.27 [tolerated]",
    716:"expr=0.23 incorp=2.17 [low-expr]",   717:"expr=0.38 incorp=1.96 [tolerated]",
}

# Log-spaced copy-number bins (edges in VAF units)
# 0.003 → 0.006 → 0.012 → 0.025 → 0.05 → 0.25
BIN_EDGES_VAF  = [0.003, 0.006, 0.012, 0.025, 0.05, 0.25]
BIN_LABELS     = ["0.5–1", "1–2", "2–4", "4–8", ">8"]   # est. copies

# ── Helpers ────────────────────────────────────────────────────────────────────

def table_or_ci(a, b, c, d):
    """2×2 OR and 95% CI with Haldane correction (add 0.5 to each cell if any zero)."""
    if min(a, b, c, d) == 0:
        a, b, c, d = a+0.5, b+0.5, c+0.5, d+0.5
    or_ = (a * d) / (b * c)
    se  = np.sqrt(1/a + 1/b + 1/c + 1/d)
    lo  = np.exp(np.log(or_) - 1.96 * se)
    hi  = np.exp(np.log(or_) + 1.96 * se)
    return or_, lo, hi

def icd_match(icd_str, code):
    """True if code appears anywhere in the pipe-separated ICD string."""
    if not isinstance(icd_str, str):
        return False
    return any(entry.strip().startswith(code) for entry in icd_str.split('|'))

# ── 1. Load top hits ───────────────────────────────────────────────────────────
print("Loading per-variant results …", flush=True)
res = pd.read_csv(RESULTS_CSV)
res = res[res['fdr'] < 0.05].copy()
# ALL block-level FDR<0.05 associations (one panel per variant x phenotype; a variant with
# >=2 significant phenotypes appears once per phenotype). Set PLOT_ALL=False to revert to top-N.
res_block = res[res['phenotype_level'] == 'block']
if PLOT_ALL:
    top = (res_block[res_block['n_carrier_cases'] >= 5]
           .sort_values('pval').reset_index(drop=True))
else:
    top = (res_block.sort_values('pval').drop_duplicates(subset=['variant_id'])
           .head(TOP_N * 2))
    top = top[top['n_carrier_cases'] >= 5].head(TOP_N)
print(f"  {len(top)} hits selected", flush=True)
print(top[['t2t_pos','t2t_ref','t2t_alt','region','phenotype',
           'n_carriers','pval','fdr','or_']].to_string(index=False))

# ── 2. Load cohort ─────────────────────────────────────────────────────────────
print("\nLoading cohort …", flush=True)
con = sqlite3.connect(str(COHORT_DB))
cohort = pd.read_sql("""
    SELECT participant_id, icd10_codes
    FROM carriers
    WHERE pc1 IS NOT NULL AND age_recruitment IS NOT NULL
      AND phenotypic_sex IS NOT NULL
      AND (sex_chr_aneuploidy IS NULL OR sex_chr_aneuploidy != 'Yes')
      AND icd10_codes IS NOT NULL AND LENGTH(icd10_codes) > 0
""", con)
con.close()
cohort['participant_id'] = cohort['participant_id'].astype(str)
pid_set = set(cohort['participant_id'])
pid_to_icd = dict(zip(cohort['participant_id'], cohort['icd10_codes']))
print(f"  {len(cohort):,} participants", flush=True)

# ── 3. Load association_input.tsv (once, keep needed variant_ids) ─────────────
needed_vids = set(top['variant_id'])
print(f"\nStreaming association_input.tsv for {len(needed_vids)} variants …", flush=True)

carrier_rows = []
for chunk in pd.read_csv(ASSOC_FILE, sep='\t',
                         dtype={'SAMPLE_ID': str, 'VAF': float},
                         chunksize=500_000):
    sub = chunk[chunk['variant_id'].isin(needed_vids) &
                (chunk['VAF'] >= VAF_THRESH) &
                chunk['SAMPLE_ID'].isin(pid_set)]
    if len(sub):
        carrier_rows.append(sub[['SAMPLE_ID', 'variant_id', 'VAF']])

carriers = pd.concat(carrier_rows, ignore_index=True)
carriers['est_copies'] = carriers['VAF'] * EST_MULT
print(f"  {len(carriers):,} carrier rows loaded", flush=True)

# ── 4. Build dose-response per hit ────────────────────────────────────────────
print("\nComputing dose-response curves …", flush=True)

results_dr = []   # (hit_row, bin_label, or_, lo, hi, n_carriers, n_cases)

for _, hit in top.iterrows():
    vid     = hit['variant_id']
    pheno   = hit['phenotype']
    region  = hit['region']

    car = carriers[carriers['variant_id'] == vid].copy()
    carrier_pids = set(car['SAMPLE_ID'])
    noncarrier_pids = pid_set - carrier_pids

    # Non-carrier baseline
    nc_cases    = sum(1 for p in noncarrier_pids if icd_match(pid_to_icd.get(p,''), pheno))
    nc_controls = len(noncarrier_pids) - nc_cases

    # Bin carriers
    car['bin'] = pd.cut(car['est_copies'],
                        bins=[e * EST_MULT for e in BIN_EDGES_VAF],
                        labels=BIN_LABELS, right=False)

    for lbl in BIN_LABELS:
        grp = car[car['bin'] == lbl]
        if len(grp) < MIN_BIN_N:
            continue
        n_cases    = sum(1 for p in grp['SAMPLE_ID'] if icd_match(pid_to_icd.get(p,''), pheno))
        n_controls = len(grp) - n_cases
        or_, lo, hi = table_or_ci(n_cases, n_controls, nc_cases, nc_controls)
        results_dr.append({
            'variant_id': vid,
            't2t_pos':    hit['t2t_pos'],
            't2t_ref':    hit['t2t_ref'],
            't2t_alt':    hit['t2t_alt'],
            'region':     region,
            'phenotype':  pheno,
            'phenotype_desc': hit['phenotype_desc'],
            'pval':       hit['pval'],
            'fdr':        hit['fdr'],
            'or_logistic': hit['or_'],
            'bin':        lbl,
            'n_bin':      len(grp),
            'n_cases':    n_cases,
            'or_':        or_,
            'ci_lo':      lo,
            'ci_hi':      hi,
        })
    print(f"  {vid} × {pheno}: {len(car):,} carriers", flush=True)

dr_df = pd.DataFrame(results_dr)
dr_df.to_csv(OUT_DIR / "85_dosage_or.csv", index=False)
print(f"Saved {len(dr_df)} rows → {OUT_DIR}/85_dosage_or.csv", flush=True)

# ── 5. Figure ──────────────────────────────────────────────────────────────────
print("\nGenerating figure …", flush=True)

n_hits = len(top)
ncols  = 3
nrows  = int(np.ceil(n_hits / ncols))

fig, axes = plt.subplots(nrows, ncols, figsize=(14, 4.5 * nrows),
                         sharey=False, sharex=False)
axes = np.array(axes).flatten()

bin_x = np.arange(len(BIN_LABELS))  # x positions for bins

for ax_i, (_, hit) in enumerate(top.iterrows()):
    ax   = axes[ax_i]
    vid  = hit['variant_id']
    pheno = hit['phenotype']
    sub  = dr_df[(dr_df['variant_id'] == vid) & (dr_df['phenotype'] == pheno)]

    # binary gene (red) vs NTS (blue) marking — applied to points AND title
    region_class = 'GENE' if hit['region'] == 'gene' else 'NTS'
    col = '#e41a1c' if region_class == 'GENE' else '#377eb8'

    # Reference line at OR=1
    ax.axhline(1.0, color='gray', lw=0.8, ls='--')

    # Plot each bin
    xs, ors, los, his, ns = [], [], [], [], []
    for lbl in BIN_LABELS:
        row = sub[sub['bin'] == lbl]
        if row.empty:
            continue
        row = row.iloc[0]
        xi  = BIN_LABELS.index(lbl)
        xs.append(xi)
        ors.append(row['or_'])
        los.append(row['or_'] - row['ci_lo'])
        his.append(row['ci_hi'] - row['or_'])
        ns.append(int(row['n_bin']))

    if xs:
        ax.errorbar(xs, ors,
                    yerr=[los, his],
                    fmt='o', color=col, markersize=7,
                    lw=1.8, capsize=4, capthick=1.5, zorder=5)
        # n labels
        for xi, or_v, n in zip(xs, ors, ns):
            ax.text(xi, or_v, f' n={n}', fontsize=6.5, va='center', color='#333333')

    # Logistic regression OR as horizontal dashed reference
    ax.axhline(hit['or_'], color=col, lw=1.0, ls=':', alpha=0.6,
               label=f"Logistic OR={hit['or_']:.2f}")

    ax.set_xticks(range(len(BIN_LABELS)))
    ax.set_xticklabels([f'{l}\ncopies' for l in BIN_LABELS], fontsize=8)
    ax.set_ylabel('Odds ratio vs non-carriers', fontsize=8)
    ax.set_yscale('log')
    ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())

    p_str  = f"p={hit['pval']:.1e}"
    fdr_str = f"FDR={hit['fdr']:.3f}"
    title = (f"T2T {int(hit['t2t_pos'])}{hit['t2t_ref']}>{hit['t2t_alt']}  "
             f"[{region_class}: {hit['region']}]\n"
             f"ICD10 {pheno}: {ICD10_DESC.get(str(pheno).strip(), str(hit['phenotype_desc']))[:55]}\n"
             f"n carriers={int(hit['n_carriers']):,}  {p_str}  {fdr_str}")
    ax.set_title(title, fontsize=7.5, loc='left', fontweight='bold', color=col)
    if hit['region'] == 'gene' and int(hit['t2t_pos']) in FUNC5S:   # 5S functional consequence
        ax.text(0.03, 0.03, "5S "+FUNC5S[int(hit['t2t_pos'])], transform=ax.transAxes,
                fontsize=6, va='bottom', ha='left', color='#5a2d82', fontweight='bold')
    ax.legend(fontsize=7, loc='upper left')

# Hide unused axes
for ax in axes[n_hits:]:
    ax.set_visible(False)

fig.suptitle(
    'UKBB 5S rDNA — Variant dosage vs disease odds\n'
    '(empirical ORs per copy-number bin vs non-carriers; unadjusted)',
    fontsize=11, fontweight='bold', y=1.01)
fig.tight_layout()

out_pdf = FIG_DIR / "85_dosage_or_curves.pdf"
FIG_DIR.mkdir(parents=True, exist_ok=True)
fig.savefig(out_pdf, bbox_inches='tight')
fig.savefig(str(out_pdf).replace('.pdf', '.png'), dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved: {out_pdf}")
print("Done.")
