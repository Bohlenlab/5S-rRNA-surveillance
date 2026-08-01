#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 81_ukbb_phenotype_associations.py — ICD10 phenotype association analysis for
# UK Biobank 5S rDNA variants (logistic burden tests, Manhattan/volcano/forest plots).
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
81_ukbb_phenotype_associations.py

ICD10 phenotype association analysis for UKBB 5S rDNA variants.

Features:
  - Calibrated VAF ≥ 0.30% threshold (vs raw AD-based tiers)
  - T2T coordinate system throughout
  - T2T functional annotations (incorp_60s_mean, rna_expr_mean) from 5S_rDNA.db
  - Region-stratified burden: gene / nts_pre / nts_post
  - Incorp-defective burden (incorp_60s_mean < 0.5)
  - Figures: Manhattan along T2T locus, volcano, forest plot, functional table

Input files:
  association_input.tsv  – per-sample carrier calls (UKBB coordinates, 15.6M rows)
  cohort ICD10 database  – cohort table with ICD10, age, sex, PC1–10 (430k participants)
  5S_rDNA.db             – functional_annotation, ukbb_population_variants

Coordinate transform (UKBB → T2T):
  t2t_pos = 1417 − ukbb_pos
  t2t_ref = complement(ukbb_ref)
  t2t_alt = complement(ukbb_alt)

Runtime: ~15 min on 10 cores. Results cached; re-run regenerates figures only.
"""

import sqlite3
import sys
import time
import threading
import warnings
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy import stats as scipy_stats
from scipy.optimize import minimize as scipy_minimize
from joblib import Parallel, delayed

warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', message='overflow encountered')
warnings.filterwarnings('ignore', message='divide by zero')
warnings.filterwarnings('ignore', message='invalid value')

import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

# ── Paths ──────────────────────────────────────────────────────────────────────

T2T_DIR    = Path(os.environ.get("FIVES_DATA", "data"))
DB_5S      = Path(os.environ.get("FIVES_DB", "5S_rDNA.db"))
TRNA_DIR   = Path(os.environ.get("FIVES_DATA", "data"))
ASSOC_FILE = TRNA_DIR / "5S_rRNA_Analysis" / "association_input.tsv"
COHORT_DB  = Path(os.environ.get("FIVES_ICD10_DB", "cohort_icd10.db"))
OUT_DIR    = Path(os.environ.get("FIVES_OUT", "output")) / "81_results"
FIG_DIR    = Path(os.environ.get("FIVES_OUT", "output")) / "02_variant_calling_qc"

# ── Configuration ──────────────────────────────────────────────────────────────

VAF_THRESHOLD    = 0.003   # calibrated 0.30% threshold
EST_COPIES_MULT  = 160     # VAF × 160 = estimated copies (80 copies/haplotype × 2)
MIN_CARRIERS     = 50      # minimum carriers per variant to run regression
MIN_CARRIER_CASES = 5      # minimum carrier-cases per phenotype×variant
CHI2_PRESCREEN_P  = 0.10   # skip IRLS if 2×2 chi-sq p > this (very loose gate)
MIN_CASES_BLOCK  = 100     # minimum total ICD10-block cases to include phenotype
MIN_CASES_CHAPTER = 50
INCORP_DEFECT_THRESH = 0.5 # incorp_60s_mean < 0.5 → defective
N_JOBS           = 10

# Continuous model thresholds
CONT_THRESH_BFGS = 20
CONT_THRESH_IRLS = 200

# ICD10 chapters
ICD10_CHAPTERS = {
    'A': 'I – Infectious', 'B': 'I – Infectious',
    'C': 'II – Neoplasms', 'D': 'III – Blood/Immune',
    'E': 'IV – Endocrine', 'F': 'V – Mental',
    'G': 'VI – Nervous', 'H': 'VII/VIII – Eye/Ear',
    'I': 'IX – Circulatory', 'J': 'X – Respiratory',
    'K': 'XI – Digestive', 'L': 'XII – Skin',
    'M': 'XIII – Musculoskeletal', 'N': 'XIV – Genitourinary',
    'O': 'XV – Pregnancy', 'Q': 'XVII – Congenital',
    'R': 'XVIII – Symptoms',
}
SEX_SPECIFIC = {'O'}

# Region colours
REG_COLORS = {"gene": "#e41a1c", "nts_pre": "#4daf4a", "nts_post": "#377eb8"}
REG_ORDER  = ["nts_pre", "gene", "nts_post"]

# ── Coordinate helpers ─────────────────────────────────────────────────────────

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

def to_t2t(ukbb_pos, ref, alt):
    return 1417 - ukbb_pos, ref.translate(_COMP), alt.translate(_COMP)

def region_of(t2t_pos):
    if t2t_pos < 630:  return "nts_pre"
    if t2t_pos <= 748: return "gene"
    return "nts_post"

# ── Logging ────────────────────────────────────────────────────────────────────

_t0 = time.time()
_log_lock = threading.Lock()

def log(msg):
    elapsed = time.time() - _t0
    m, s = divmod(int(elapsed), 60)
    line = f"[{m:02d}:{s:02d}] {msg}"
    with _log_lock:
        print(line, flush=True)

# ── Statistics ─────────────────────────────────────────────────────────────────

def _run_logistic(y, X, use_bfgs=False):
    n, p = X.shape
    beta = np.zeros(p)
    converged = False
    method = 'logistic'

    for _ in range(100 if use_bfgs else 50):
        eta = np.clip(X @ beta, -500, 500)
        mu  = 1.0 / (1.0 + np.exp(-eta))
        W   = np.maximum(mu * (1.0 - mu), 1e-10)
        XtW  = X.T * W
        XtWX = XtW @ X
        adj  = eta + (y - mu) / W
        try:
            beta_new = np.linalg.solve(XtWX, XtW @ adj)
        except np.linalg.LinAlgError:
            break
        if np.max(np.abs(beta_new - beta)) < 1e-7:
            beta = beta_new; converged = True; break
        beta = beta_new

    if not converged and use_bfgs:
        def neg_ll(b):
            mu = 1.0 / (1.0 + np.exp(-np.clip(X @ b, -500, 500)))
            return -np.sum(y * np.log(mu + 1e-15) + (1-y) * np.log(1-mu + 1e-15))
        def grad(b):
            mu = 1.0 / (1.0 + np.exp(-np.clip(X @ b, -500, 500)))
            return X.T @ (mu - y)
        try:
            res = scipy_minimize(neg_ll, beta, jac=grad, method='BFGS',
                                 options={'maxiter': 200, 'disp': False})
            beta = res.x; converged = res.success; method = 'logistic_bfgs'
        except Exception:
            pass

    try:
        eta  = np.clip(X @ beta, -500, 500)
        mu   = 1.0 / (1.0 + np.exp(-eta))
        W    = np.maximum(mu * (1.0 - mu), 1e-10)
        cov  = np.linalg.inv((X.T * W) @ X)
        se   = np.sqrt(np.diag(cov))
    except np.linalg.LinAlgError:
        return np.nan, np.nan, np.nan, False, 'failed'

    b, s = beta[1], se[1]
    if not (np.isfinite(b) and np.isfinite(s) and s > 0):
        return np.nan, np.nan, np.nan, False, 'failed:se'
    return b, s, 2.0 * float(scipy_stats.norm.sf(abs(b/s))), converged, method


def _or_ci(beta, se):
    if not (np.isfinite(beta) and np.isfinite(se) and se > 0):
        return np.nan, np.nan, np.nan
    return np.exp(beta), np.exp(beta - 1.96*se), np.exp(beta + 1.96*se)


def _chi2_prescreen(carrier_mask, case_vec):
    """Fast 2×2 chi-squared p-value. Returns 1.0 if any expected count < 1."""
    a = int(( carrier_mask &  case_vec).sum())
    b = int(( carrier_mask & ~case_vec).sum())
    c = int((~carrier_mask &  case_vec).sum())
    d = int((~carrier_mask & ~case_vec).sum())
    n = a + b + c + d
    if n == 0: return 1.0
    r1, r2 = a + b, c + d
    c1, c2 = a + c, b + d
    e_a = r1 * c1 / n
    if min(e_a, r1*c2/n, r2*c1/n, r2*c2/n) < 1.0:
        return 1.0
    chi2 = ((a-e_a)**2/e_a + (b-r1*c2/n)**2/(r1*c2/n) +
            (c-r2*c1/n)**2/(r2*c1/n) + (d-r2*c2/n)**2/(r2*c2/n))
    return float(scipy_stats.chi2.sf(chi2, df=1))


def bh_fdr(pvalues, m_total=None):
    """BH-FDR. m_total sets the denominator (total tests including skipped ones).
    When a prescreen is used, pass the full test-space size so skipped tests
    count as p=1.0 — otherwise FDR is artificially liberal."""
    arr = np.asarray(pvalues, dtype=float)
    q   = np.full(len(arr), np.nan)
    idx = np.where(np.isfinite(arr))[0]
    if not len(idx): return q
    p = arr[idx]
    m = m_total if m_total is not None else len(idx)
    order = np.argsort(p)
    q_s = p[order] * m / (np.arange(1, len(p)+1))
    for i in range(len(p)-2, -1, -1): q_s[i] = min(q_s[i], q_s[i+1])
    q_s = np.minimum(q_s, 1.0)
    q[idx] = q_s[np.argsort(order)]
    return q

# ── Data loading ───────────────────────────────────────────────────────────────

def load_functional_annotations():
    """Load incorp/expr from 5S_rDNA.db as DataFrame (for vectorised merge)."""
    log("Loading T2T functional annotations …")
    con = sqlite3.connect(DB_5S)
    func_df = pd.read_sql("""
        SELECT consensus_pos AS t2t_pos, alt_base AS t2t_alt,
               incorp_60s_mean AS incorp, rna_expr_mean AS expr
        FROM functional_annotation
    """, con)
    con.close()
    func_df['t2t_pos'] = func_df['t2t_pos'].astype(int)
    func_df['t2t_alt'] = func_df['t2t_alt'].astype(str)
    log(f"  {len(func_df):,} position×alt annotations loaded")
    return func_df


def load_association_data(func_df):
    """
    Load association_input.tsv, filter VAF ≥ 0.003, map to T2T coordinates.
    Returns:
      variant_data: {variant_id: (sample_ids, est_copies)}
      variant_meta: {variant_id: {t2t_pos, t2t_ref, t2t_alt, ukbb_pos, region, incorp, expr, ...}}
      burden_df: per-sample burden scores (indexed by SAMPLE_ID)
    """
    log(f"Loading {ASSOC_FILE.name}  "
        f"({ASSOC_FILE.stat().st_size/1e9:.1f} GB) …")
    t = time.time()
    df = pd.read_csv(ASSOC_FILE, sep='\t',
                     dtype={'SAMPLE_ID': str, 'VAF': float, 'est_copies': float})
    log(f"  Loaded {len(df):,} rows in {time.time()-t:.1f}s")

    # Apply VAF threshold
    n_before = len(df)
    df = df[df['VAF'] >= VAF_THRESHOLD].copy()
    log(f"  After VAF ≥ {VAF_THRESHOLD:.3f}: {len(df):,} rows "
        f"({len(df)/n_before:.1%} kept)")

    # Recompute est_copies from VAF (overrides pre-computed column)
    df['est_copies'] = (df['VAF'] * EST_COPIES_MULT).astype(np.float32)

    # Map to T2T coordinates (vectorised)
    df['t2t_pos'] = (1417 - df['POS']).astype(int)
    df['t2t_ref'] = df['REF'].str.translate(_COMP)
    df['t2t_alt'] = df['ALT'].str.translate(_COMP)
    df['region']  = df['t2t_pos'].map(region_of)

    # Join functional annotations (vectorised merge)
    df = df.merge(func_df[['t2t_pos', 't2t_alt', 'incorp', 'expr']],
                  on=['t2t_pos', 't2t_alt'], how='left')

    # Burden scores per sample
    df['incorp_defect']  = (df['incorp'] < INCORP_DEFECT_THRESH) & df['incorp'].notna()
    df['ec_incorp_def']  = np.where(df['incorp_defect'],  df['est_copies'], 0.0)
    df['ec_gene']        = np.where(df['region'] == 'gene', df['est_copies'], 0.0)
    df['ec_nts_pre']     = np.where(df['region'] == 'nts_pre', df['est_copies'], 0.0)
    df['ec_nts_post']    = np.where(df['region'] == 'nts_post', df['est_copies'], 0.0)

    burden = df.groupby('SAMPLE_ID').agg(
        burden_total      = ('est_copies',   'sum'),
        burden_gene       = ('ec_gene',      'sum'),
        burden_nts_pre    = ('ec_nts_pre',   'sum'),
        burden_nts_post   = ('ec_nts_post',  'sum'),
        burden_incorp_def = ('ec_incorp_def','sum'),
    ).reset_index()

    # Per-variant lookup
    log("  Building per-variant lookup …")
    variant_data = {}
    variant_meta = {}
    for vid, grp in df.groupby('variant_id'):
        t2t_pos = int(grp['t2t_pos'].iloc[0])
        t2t_ref = str(grp['t2t_ref'].iloc[0])
        t2t_alt = str(grp['t2t_alt'].iloc[0])
        incorp  = float(grp['incorp'].iloc[0]) if pd.notna(grp['incorp'].iloc[0]) else np.nan
        expr    = float(grp['expr'].iloc[0])   if pd.notna(grp['expr'].iloc[0])   else np.nan
        variant_data[vid] = (
            grp['SAMPLE_ID'].values,
            grp['est_copies'].values.astype(np.float32),
        )
        variant_meta[vid] = {
            'ukbb_pos': int(grp['POS'].iloc[0]),
            't2t_pos':  t2t_pos,
            't2t_ref':  t2t_ref,
            't2t_alt':  t2t_alt,
            'region':   region_of(t2t_pos),
            'n_carriers': len(grp),
            'incorp':   incorp,
            'expr':     expr,
            'incorp_defect': np.isfinite(incorp) and incorp < INCORP_DEFECT_THRESH,
        }

    del df
    log(f"  {len(variant_data):,} unique variants after VAF filter; "
        f"{burden.shape[0]:,} samples with any call")
    return variant_data, variant_meta, burden


def load_cohort():
    log("Loading cohort …")
    con = sqlite3.connect(str(COHORT_DB))
    cohort = pd.read_sql("""
        SELECT participant_id, icd10_codes, age_recruitment, phenotypic_sex,
               pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9,pc10
        FROM carriers
        WHERE pc1 IS NOT NULL AND age_recruitment IS NOT NULL
          AND phenotypic_sex IS NOT NULL
          AND (sex_chr_aneuploidy IS NULL OR sex_chr_aneuploidy != 'Yes')
          AND icd10_codes IS NOT NULL AND LENGTH(icd10_codes) > 0
    """, con)
    con.close()
    log(f"  {len(cohort):,} participants")
    cohort['sex_binary'] = (cohort['phenotypic_sex'] == 'Male').astype(float)
    for col in [f'pc{i}' for i in range(1, 11)]:
        mu, sd = cohort[col].mean(), cohort[col].std()
        cohort[col] = (cohort[col] - mu) / max(sd, 1e-10)
    return cohort


def build_phenotype_vectors(cohort):
    log("Building ICD10 phenotype vectors …")
    t = time.time()
    n = len(cohort)
    pids = cohort['participant_id'].values
    pid_to_idx = {str(pid): i for i, pid in enumerate(pids)}

    cov_cols = ['age_recruitment', 'sex_binary'] + [f'pc{i}' for i in range(1, 11)]
    covariates = np.column_stack([np.ones(n), cohort[cov_cols].values.astype(float)])

    # Parse ICD10 codes
    participant_codes = {}
    for row in cohort.itertuples(index=False):
        pid = str(row.participant_id)
        codes = set()
        icd_str = row.icd10_codes
        if icd_str and isinstance(icd_str, str):
            for entry in icd_str.split('|'):
                entry = entry.strip()
                if entry:
                    code = entry.split(' ', 1)[0].strip()
                    if code:
                        codes.add(code)
        participant_codes[pid] = codes

    all_codes = set().union(*participant_codes.values())

    def _build_level(level, min_cases):
        if level == 'chapter':
            label_fn = lambda c: c[0].upper() if c else ''
            desc_map = ICD10_CHAPTERS
            valid = {label_fn(c) for c in all_codes
                     if label_fn(c) in ICD10_CHAPTERS and label_fn(c) not in SEX_SPECIFIC}
        else:
            label_fn = lambda c: c[:3].upper() if len(c) >= 3 else ''
            valid = {label_fn(c) for c in all_codes
                     if len(c) >= 3 and label_fn(c) not in
                     # skip pregnancy blocks
                     {'O00','O01','O02','O03','O04','O05','O06','O07','O08','O09',
                      'O10','O11','O12','O13','O14','O15','O16','O20','O21','O22'}}
            desc_map = {}

        label_list = sorted(valid)
        label_idx  = {lb: i for i, lb in enumerate(label_list)}
        vecs = np.zeros((len(label_list), n), dtype=bool)

        for pid, codes in participant_codes.items():
            i = pid_to_idx.get(pid)
            if i is None: continue
            for code in codes:
                lb = label_fn(code)
                li = label_idx.get(lb)
                if li is not None:
                    vecs[li, i] = True

        phenos = {}
        for lb, li in label_idx.items():
            nc = int(vecs[li].sum())
            if nc >= min_cases:
                desc = desc_map.get(lb, lb)
                phenos[lb] = (vecs[li], desc, nc)
        return phenos

    chapter_phenos = _build_level('chapter', MIN_CASES_CHAPTER)
    block_phenos   = _build_level('block',   MIN_CASES_BLOCK)
    log(f"  {len(chapter_phenos)} chapter | {len(block_phenos)} block phenotypes "
        f"({time.time()-t:.1f}s)")
    return pid_to_idx, covariates, chapter_phenos, block_phenos

# ── Per-variant worker ─────────────────────────────────────────────────────────

def process_variant(vid, carrier_sids, carrier_copies,
                    pid_to_idx, n_cohort, covariates,
                    chapter_phenos, block_phenos, meta):
    est_copies = np.zeros(n_cohort, dtype=np.float32)
    idxs = np.array([pid_to_idx.get(str(s), -1) for s in carrier_sids])
    valid = idxs >= 0
    est_copies[idxs[valid]] = carrier_copies[valid]

    n_carriers   = int((est_copies > 0).sum())
    if n_carriers < MIN_CARRIERS:
        return []

    carrier_mask = est_copies > 0
    results = []

    for level, pheno_dict in [('chapter', chapter_phenos), ('block', block_phenos)]:
        for pheno_label, (case_vec, pheno_desc, n_total_cases) in pheno_dict.items():
            n_exp_cases = int((carrier_mask & case_vec).sum())
            if n_exp_cases < MIN_CARRIER_CASES:
                continue
            if n_total_cases < CONT_THRESH_BFGS:
                continue

            # Fast pre-screen: skip IRLS if no nominal 2×2 signal
            if _chi2_prescreen(carrier_mask, case_vec) > CHI2_PRESCREEN_P:
                continue

            use_bfgs = (n_total_cases < CONT_THRESH_IRLS)
            y     = case_vec.astype(float)
            ec    = est_copies.astype(float)
            X     = np.column_stack([np.ones(n_cohort), ec, covariates[:, 1:]])
            b, s, p, conv, meth = _run_logistic(y, X, use_bfgs=use_bfgs)
            or_, lo, hi = _or_ci(b, s)

            # Z-scored (per-SD-of-carriers)
            sd_c = float(ec[carrier_mask].std()) if carrier_mask.sum() >= 2 else 1.0
            sd_c = max(sd_c, 1e-10)

            results.append({
                'variant_id': vid,
                'ukbb_pos':   meta['ukbb_pos'],
                't2t_pos':    meta['t2t_pos'],
                't2t_ref':    meta['t2t_ref'],
                't2t_alt':    meta['t2t_alt'],
                'region':     meta['region'],
                'incorp':     meta['incorp'],
                'expr':       meta['expr'],
                'incorp_defect': meta['incorp_defect'],
                'n_carriers': n_carriers,
                'phenotype':  pheno_label,
                'phenotype_desc': pheno_desc,
                'phenotype_level': level,
                'n_total_cases': n_total_cases,
                'n_carrier_cases': n_exp_cases,
                'prev_carrier':    n_exp_cases / max(n_carriers, 1),
                'prev_noncarrier': int((~carrier_mask & case_vec).sum()) / max(n_cohort - n_carriers, 1),
                'beta':  b,  'se':    s,   'pval':   p,
                'or_':   or_, 'ci_lo': lo,  'ci_hi':  hi,
                'sd_carriers': sd_c,
                'or_z':  np.exp(b * sd_c) if np.isfinite(b) else np.nan,
                'method': meth, 'converged': conv,
            })

    return results


def run_per_variant_analysis(variant_data, variant_meta,
                              pid_to_idx, n_cohort, covariates,
                              chapter_phenos, block_phenos):
    vids = [v for v in variant_data
            if variant_meta[v]['n_carriers'] >= MIN_CARRIERS]
    log(f"\nPer-variant regression: {len(vids):,} variants × "
        f"{len(chapter_phenos)+len(block_phenos)} phenotypes, n_jobs={N_JOBS}")

    batch = max(1, len(vids) // 20)
    all_results = []
    t0 = time.time()
    for start in range(0, len(vids), batch):
        chunk = vids[start: start+batch]
        res_lists = Parallel(n_jobs=N_JOBS, backend='loky')(
            delayed(process_variant)(
                v, variant_data[v][0], variant_data[v][1],
                pid_to_idx, n_cohort, covariates,
                chapter_phenos, block_phenos, variant_meta[v],
            )
            for v in chunk
        )
        for rl in res_lists: all_results.extend(rl)
        done = min(start+batch, len(vids))
        elapsed = time.time()-t0
        eta = (elapsed/done*(len(vids)-done)) if done > 0 else 0
        log(f"  {done}/{len(vids)} variants ({done/len(vids):.0%}) "
            f"| {elapsed/60:.1f}m elapsed | ETA {eta/60:.1f}m "
            f"| associations so far: {len(all_results):,}")

    log(f"  Per-variant done: {len(all_results):,} associations tested")
    return all_results


def run_burden_analysis(burden_df, pid_to_idx, n_cohort, covariates,
                        chapter_phenos, block_phenos):
    log("\nBurden analysis …")
    score_cols = ['burden_total', 'burden_gene', 'burden_nts_pre',
                  'burden_nts_post', 'burden_incorp_def']
    score_labels = {
        'burden_total':      'All variants',
        'burden_gene':       'Gene-body only',
        'burden_nts_pre':    'NTS-pre only',
        'burden_nts_post':   'NTS-post only',
        'burden_incorp_def': 'Incorp-defective (<0.5)',
    }

    # Build z-scored burden vectors over cohort (vectorised)
    burden_vecs = {}
    sids   = burden_df['SAMPLE_ID'].astype(str).values
    idxs   = np.array([pid_to_idx.get(s, -1) for s in sids])
    valid  = idxs >= 0
    for sc in score_cols:
        vec = np.zeros(n_cohort, dtype=float)
        vals = burden_df[sc].values.astype(float)
        vec[idxs[valid]] = vals[valid]
        mu, sd = vec.mean(), vec.std()
        burden_vecs[sc] = (vec - mu) / max(sd, 1e-10)

    # Pre-build base design matrix (intercept + covariates); swap col 1 per score
    X_base = np.column_stack([np.ones(n_cohort), np.zeros(n_cohort), covariates[:, 1:]])

    results = []
    for sc, bvec in burden_vecs.items():
        X_base[:, 1] = bvec          # swap burden column in-place
        X = X_base                   # reuse same allocation
        for level, pheno_dict in [('chapter', chapter_phenos), ('block', block_phenos)]:
            for pheno_label, (case_vec, pheno_desc, n_total_cases) in pheno_dict.items():
                if n_total_cases < CONT_THRESH_BFGS: continue
                y = case_vec.astype(float)
                use_bfgs = (n_total_cases < CONT_THRESH_IRLS)
                b, s, p, conv, meth = _run_logistic(y, X, use_bfgs=use_bfgs)
                or_, lo, hi = _or_ci(b, s)
                results.append({
                    'score': sc,
                    'score_label': score_labels[sc],
                    'phenotype': pheno_label,
                    'phenotype_desc': pheno_desc,
                    'phenotype_level': level,
                    'n_total_cases': n_total_cases,
                    'beta': b, 'se': s, 'pval': p,
                    'or_': or_, 'ci_lo': lo, 'ci_hi': hi,
                    'method': meth, 'converged': conv,
                })

    log(f"  Burden done: {len(results):,} associations")
    return results

# ── Figures ────────────────────────────────────────────────────────────────────

def make_figures(var_df, burden_df, out_pdf):
    log("Generating figures …")
    import matplotlib.backends.backend_pdf as pdf_backend
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Prepare block-level per-variant data
    bl = var_df[var_df['phenotype_level'] == 'block'].copy()
    bl = bl[bl['pval'].notna() & bl['or_'].notna()]
    bl['logp']   = -np.log10(bl['pval'].clip(1e-300))
    bl['log2or'] = np.log2(bl['or_'].clip(1e-4, 1e4))

    # ── Figure layout: 2×3 top + 2×2 bottom ──────────────────────────────────
    fig = plt.figure(figsize=(18, 22))
    gs  = fig.add_gridspec(3, 3, hspace=0.45, wspace=0.38,
                           height_ratios=[1.2, 1.2, 1.0])

    # ── Panel A: Manhattan along T2T locus ────────────────────────────────────
    ax_A = fig.add_subplot(gs[0, :])
    fdr_thresh = -np.log10(0.05)

    for reg in REG_ORDER:
        sub = bl[bl['region'] == reg]
        ax_A.scatter(sub['t2t_pos'], sub['logp'],
                     s=5, c=REG_COLORS[reg], alpha=0.35, zorder=3, linewidths=0,
                     label=f"{reg} (n={len(sub):,})")

    # FDR-significant highlights
    sig = bl[bl['fdr'] < 0.05]
    ax_A.scatter(sig['t2t_pos'], sig['logp'],
                 s=20, c='black', zorder=5, linewidths=0, label='FDR<0.05')

    # Gene-body shading
    ax_A.axvspan(630, 748, color='lightgray', alpha=0.3, zorder=1)
    ax_A.axhline(fdr_thresh, color='gray', lw=0.8, ls='--')
    ax_A.set_xlabel('T2T consensus position', fontsize=10)
    ax_A.set_ylabel('-log₁₀(p)', fontsize=10)
    ax_A.set_title('A  Manhattan: UKBB 5S rDNA variant × ICD10 phenotype associations '
                   f'(VAF ≥ 0.30%, block-level)',
                   fontsize=10, loc='left', fontweight='bold')
    ax_A.set_xlim(460, 970)
    ax_A.legend(fontsize=7, ncol=4, loc='upper right')

    # Annotate top-5 per region
    for reg in REG_ORDER:
        sub_sig = bl[(bl['region'] == reg) & (bl['fdr'] < 0.05)]
        for _, row in sub_sig.nsmallest(5, 'pval').iterrows():
            ax_A.annotate(f"{row['phenotype']}",
                          (row['t2t_pos'], row['logp']),
                          fontsize=5, ha='center', va='bottom',
                          xytext=(0, 2), textcoords='offset points',
                          color=REG_COLORS[reg])

    # ── Panel B: Volcano ───────────────────────────────────────────────────────
    ax_B = fig.add_subplot(gs[1, 0])
    for reg in REG_ORDER:
        sub = bl[bl['region'] == reg]
        ax_B.scatter(sub['log2or'], sub['logp'],
                     s=4, c=REG_COLORS[reg], alpha=0.4, zorder=3, linewidths=0)
    ax_B.scatter(sig['log2or'], sig['logp'],
                 s=15, facecolors='none', edgecolors='black', lw=0.6, zorder=5)
    ax_B.axhline(fdr_thresh, color='gray', lw=0.8, ls='--')
    ax_B.axvline(0, color='gray', lw=0.5)
    ax_B.set_xlabel('log₂(OR)', fontsize=9)
    ax_B.set_ylabel('-log₁₀(p)', fontsize=9)
    ax_B.set_title('B  Volcano (block-level, per-variant)',
                   fontsize=9, loc='left', fontweight='bold')
    patches = [mpatches.Patch(color=REG_COLORS[r], label=r) for r in REG_ORDER]
    ax_B.legend(handles=patches, fontsize=7)

    # ── Panel C: Volcano — gene body only with functional colour ─────────────
    ax_C = fig.add_subplot(gs[1, 1])
    gene_bl = bl[bl['region'] == 'gene'].copy()
    has_incorp = gene_bl['incorp'].notna()
    def_mask   = has_incorp & (gene_bl['incorp'] < INCORP_DEFECT_THRESH)
    norm_mask  = has_incorp & (gene_bl['incorp'] >= INCORP_DEFECT_THRESH)
    no_mask    = ~has_incorp

    for mask, col, label in [
        (no_mask,   '#aaaaaa', 'No incorp data'),
        (norm_mask, '#2166ac', f'incorp ≥ {INCORP_DEFECT_THRESH}'),
        (def_mask,  '#d6604d', f'incorp < {INCORP_DEFECT_THRESH} (defective)'),
    ]:
        sub = gene_bl[mask]
        ax_C.scatter(sub['log2or'], sub['logp'],
                     s=8, c=col, alpha=0.6, zorder=3, linewidths=0, label=label)

    ax_C.axhline(fdr_thresh, color='gray', lw=0.8, ls='--')
    ax_C.axvline(0, color='gray', lw=0.5)
    ax_C.set_xlabel('log₂(OR)', fontsize=9)
    ax_C.set_ylabel('-log₁₀(p)', fontsize=9)
    ax_C.set_title('C  Gene-body variants — functional annotation',
                   fontsize=9, loc='left', fontweight='bold')
    ax_C.legend(fontsize=7)

    # ── Panel D: Forest plot — top hits ───────────────────────────────────────
    ax_D = fig.add_subplot(gs[1, 2])
    top_hits = bl[bl['n_carrier_cases'] >= 20].nsmallest(20, 'pval').copy()
    top_hits = top_hits.drop_duplicates(subset=['t2t_pos', 'phenotype']).head(15)

    if not top_hits.empty:
        top_hits = top_hits.reset_index(drop=True)
        ys = np.arange(len(top_hits))
        for i, row in top_hits.iterrows():
            col = REG_COLORS.get(row['region'], 'gray')
            ax_D.errorbar(row['log2or'], ys[i],
                          xerr=[[row['log2or'] - np.log2(row['ci_lo'])],
                                [np.log2(row['ci_hi']) - row['log2or']]],
                          fmt='o', color=col, markersize=4, lw=1.2, capsize=2)
        ax_D.axvline(0, color='gray', lw=0.8, ls='--')
        labels = [f"{r['phenotype']} | T2T {r['t2t_pos']}{r['t2t_ref']}>{r['t2t_alt']}"
                  for _, r in top_hits.iterrows()]
        ax_D.set_yticks(ys); ax_D.set_yticklabels(labels, fontsize=6)
        ax_D.set_xlabel('log₂(OR)', fontsize=9)
        ax_D.set_title('D  Forest — top hits (≥20 carrier-cases)',
                       fontsize=9, loc='left', fontweight='bold')
    else:
        ax_D.text(0.5, 0.5, 'No hits with ≥20 carrier-cases',
                  ha='center', va='center', transform=ax_D.transAxes, fontsize=9)
        ax_D.set_title('D  Forest', fontsize=9, loc='left', fontweight='bold')

    # ── Panel E: Burden analysis — top block hits ──────────────────────────────
    ax_E = fig.add_subplot(gs[2, :2])
    if not burden_df.empty and 'phenotype_level' in burden_df.columns:
        brd_block = burden_df[(burden_df['phenotype_level'] == 'block') &
                              burden_df['fdr'].notna()].copy()
        if not brd_block.empty:
            brd_block['logp'] = -np.log10(brd_block['pval'].clip(1e-300))
    else:
        brd_block = pd.DataFrame()

    score_colors = {
        'burden_total':      '#666666',
        'burden_gene':       REG_COLORS['gene'],
        'burden_nts_pre':    REG_COLORS['nts_pre'],
        'burden_nts_post':   REG_COLORS['nts_post'],
        'burden_incorp_def': '#ff7f00',
    }

    _score_labels = {
        'burden_total': 'All variants', 'burden_gene': 'Gene-body',
        'burden_nts_pre': 'NTS-pre', 'burden_nts_post': 'NTS-post',
        'burden_incorp_def': 'Incorp-defective',
    }
    brd_plot = brd_block.nsmallest(100, 'pval') if not brd_block.empty else brd_block
    for sc, col in score_colors.items():
        sub = brd_plot[brd_plot['score'] == sc] if not brd_plot.empty else pd.DataFrame()
        if sub.empty: continue
        ax_E.scatter(np.log2(sub['or_'].clip(1e-4, 1e4)), sub['logp'],
                     s=8, c=col, alpha=0.5, linewidths=0,
                     label=_score_labels.get(sc, sc))

    ax_E.axhline(fdr_thresh, color='gray', lw=0.8, ls='--')
    ax_E.axvline(0, color='gray', lw=0.5)
    ax_E.set_xlabel('log₂(OR) of z-scored burden', fontsize=9)
    ax_E.set_ylabel('-log₁₀(p)', fontsize=9)
    ax_E.set_title('E  Burden analysis — z-scored burden score × ICD10 block',
                   fontsize=9, loc='left', fontweight='bold')
    ax_E.legend(fontsize=7, ncol=2)

    # ── Panel F: Summary table ─────────────────────────────────────────────────
    ax_F = fig.add_subplot(gs[2, 2])
    ax_F.axis('off')

    n_var_sig   = int((var_df['fdr'] < 0.05).sum())
    n_brd_sig   = int((burden_df['fdr'] < 0.05).sum()) if 'fdr' in burden_df else 0
    n_var_tests = len(var_df)
    n_brd_tests = len(burden_df)

    # Region breakdown of FDR hits
    reg_counts = var_df[var_df['fdr'] < 0.05].groupby('region').size().to_dict()
    gene_sig  = reg_counts.get('gene', 0)
    pre_sig   = reg_counts.get('nts_pre', 0)
    post_sig  = reg_counts.get('nts_post', 0)

    # Top variants (unique positions) by min pval
    top5 = var_df.groupby(['t2t_pos', 't2t_ref', 't2t_alt'])['pval'].min().nsmallest(5)

    lines = [
        f"Per-variant analysis",
        f"  Tested: {n_var_tests:,} variant × phenotype",
        f"  FDR < 0.05: {n_var_sig:,} associations",
        f"    gene: {gene_sig}  nts_pre: {pre_sig}  nts_post: {post_sig}",
        f"",
        f"Burden analysis",
        f"  Tested: {n_brd_tests:,} score × phenotype",
        f"  FDR < 0.05: {n_brd_sig:,}",
        f"",
        f"Top variants (by min p across phenotypes):",
    ]
    for (pos, ref, alt), p in top5.items():
        lines.append(f"  T2T {pos} {ref}>{alt}  p={p:.2e}")

    ax_F.text(0.05, 0.95, '\n'.join(lines), transform=ax_F.transAxes,
              fontsize=8, va='top', fontfamily='monospace',
              bbox=dict(boxstyle='round', fc='#f8f8f8', ec='#cccccc', lw=0.8))
    ax_F.set_title('F  Summary', fontsize=9, loc='left', fontweight='bold')

    fig.suptitle(
        f'UKBB 5S rDNA Variant Phenotype Associations  '
        f'(n≈430k, VAF ≥ 0.30%, BH-FDR corrected)',
        fontsize=11, y=0.995, fontweight='bold'
    )

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_pdf, bbox_inches='tight')
    fig.savefig(str(out_pdf).replace('.pdf', '.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    log(f"  Saved: {out_pdf}")


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--full', action='store_true',
                        help='Disable chi-sq pre-screen (run all variant×phenotype pairs)')
    parser.add_argument('--server', action='store_true',
                        help='Use alternate (server) paths')
    parser.add_argument('--njobs', type=int, default=None,
                        help='Override N_JOBS (default: use config value)')
    args = parser.parse_args()

    global OUT_DIR, FIG_DIR, ASSOC_FILE, DB_5S, COHORT_DB, CHI2_PRESCREEN_P, N_JOBS

    if args.server:
        _srv = Path(os.environ.get("FIVES_DATA", "data"))
        ASSOC_FILE = _srv / "results_500k/association/association_input.tsv"
        DB_5S      = Path(os.environ.get("FIVES_DB", "5S_rDNA.db"))
        COHORT_DB  = Path(os.environ.get("FIVES_ICD10_DB", "cohort_icd10.db"))
        OUT_DIR    = Path(os.environ.get("FIVES_OUT", "output")) / "81_results"
        FIG_DIR    = Path(os.environ.get("FIVES_OUT", "output")) / "02_variant_calling_qc"
        N_JOBS     = 70

    if args.full:
        OUT_DIR = OUT_DIR.parent / "81_results_full"
        CHI2_PRESCREEN_P = 1.0

    if args.njobs is not None:
        N_JOBS = args.njobs

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    mode = "FULL (no prescreen)" if args.full else f"PRESCREENED (chi-sq p<{CHI2_PRESCREEN_P})"
    log("=" * 65)
    log(f"81_ukbb_phenotype_associations.py  [{mode}]")
    log(f"VAF threshold: {VAF_THRESHOLD:.3f} | n_jobs: {N_JOBS}")
    log(f"Output: {OUT_DIR}")
    log("=" * 65)

    # ── 1. Load data ───────────────────────────────────────────────────────────
    func              = load_functional_annotations()
    variant_data, variant_meta, burden_sample_df = load_association_data(func)
    cohort            = load_cohort()
    pid_to_idx, covariates, chapter_phenos, block_phenos = build_phenotype_vectors(cohort)
    n_cohort          = len(cohort)
    del cohort

    # ── 2. Per-variant regression ──────────────────────────────────────────────
    cache_var = OUT_DIR / "per_variant_results.csv"

    if cache_var.exists():
        log(f"\nLoading cached per-variant results from {cache_var} …")
        var_df = pd.read_csv(cache_var)
        log(f"  {len(var_df):,} rows loaded; FDR<0.05: {(var_df['fdr']<0.05).sum():,}")
    else:
        var_results = run_per_variant_analysis(
            variant_data, variant_meta, pid_to_idx, n_cohort,
            covariates, chapter_phenos, block_phenos)

        var_df = pd.DataFrame(var_results)
        if not var_df.empty:
            # m_total = full test space (variants tested × all phenotypes),
            # so skipped tests (prescreen p>0.10) count as p=1.0 in BH
            n_vids  = sum(1 for v in variant_meta
                          if variant_meta[v]['n_carriers'] >= MIN_CARRIERS)
            m_total = n_vids * (len(chapter_phenos) + len(block_phenos))
            log(f"  FDR denominator: {m_total:,} (full test space, "
                f"{len(var_df):,} tests ran)")
            var_df['fdr']  = bh_fdr(var_df['pval'].values, m_total=m_total)
            var_df['bonf'] = np.minimum(var_df['pval'].values * m_total, 1.0)
            var_df['sig_fdr']  = var_df['fdr']  < 0.05
            var_df['sig_bonf'] = var_df['bonf'] < 0.05
            var_df.to_csv(cache_var, index=False)
            log(f"  Saved {len(var_df):,} rows → {cache_var}")
            log(f"  FDR<0.05: {var_df['sig_fdr'].sum():,}  "
                f"Bonf<0.05: {var_df['sig_bonf'].sum():,}")
        else:
            log("  No per-variant results.")
            var_df = pd.DataFrame()

    # Save strong hits
    if not var_df.empty:
        strong = var_df[(var_df['fdr'] < 0.05) & (var_df['n_carrier_cases'] >= 50)]
        strong.to_csv(OUT_DIR / "strong_hits_50carriercases.csv", index=False)
        log(f"  Strong hits (FDR<0.05, ≥50 carrier-cases): {len(strong):,}")
        gene_hits = var_df[(var_df['fdr'] < 0.05) & (var_df['region'] == 'gene')]
        gene_hits.to_csv(OUT_DIR / "gene_body_hits.csv", index=False)
        log(f"  Gene-body FDR hits: {len(gene_hits):,}")

    # ── 3. Burden analysis ─────────────────────────────────────────────────────
    cache_brd = OUT_DIR / "burden_results.csv"

    if cache_brd.exists():
        log(f"\nLoading cached burden results from {cache_brd} …")
        brd_df = pd.read_csv(cache_brd)
        log(f"  {len(brd_df):,} rows; FDR<0.05: {(brd_df['fdr']<0.05).sum():,}")
    else:
        brd_results = run_burden_analysis(
            burden_sample_df, pid_to_idx, n_cohort, covariates,
            chapter_phenos, block_phenos)
        brd_df = pd.DataFrame(brd_results)
        if not brd_df.empty:
            brd_df['fdr']  = bh_fdr(brd_df['pval'].values)
            brd_df['bonf'] = np.minimum(brd_df['pval'].values *
                                         brd_df['pval'].notna().sum(), 1.0)
            brd_df.to_csv(cache_brd, index=False)
            log(f"  Saved {len(brd_df):,} rows → {cache_brd}")
            log(f"  FDR<0.05: {(brd_df['fdr']<0.05).sum():,}")
        else:
            brd_df = pd.DataFrame()

    # ── 4. Figures ─────────────────────────────────────────────────────────────
    suffix  = "_full" if args.full else ""
    out_pdf = FIG_DIR / f"81_ukbb_phenotype_associations{suffix}.pdf"
    if not var_df.empty:
        make_figures(var_df, brd_df if not brd_df.empty else pd.DataFrame(),
                     out_pdf)
    else:
        log("No results to plot.")

    elapsed = time.time() - _t0
    log(f"\nDone. Total: {elapsed/60:.1f} min")
    log(f"Results: {OUT_DIR}")
    log(f"Figure:  {out_pdf}")


if __name__ == "__main__":
    main()
