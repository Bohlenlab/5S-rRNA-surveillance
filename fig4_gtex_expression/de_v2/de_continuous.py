# -----------------------------------------------------------------------------
# de_continuous.py — per-tissue DESeq2 DE against a continuous per-donor 5S-variant
# dosage predictor, combined across tissues downstream.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""Whole-cohort continuous 5S-variant DE (per-tissue, continuous dosage predictor).

Unlike de_pertissue.py's binary carrier-vs-non-carrier contrasts, here the 5S predictor is a
CONTINUOUS per-donor dosage metric (from donor_metrics.tsv, standardized; non-carriers = 0), so a
single 'metric' coefficient measures the expression gradient with dosage across ALL donors. Fit per
tissue with the standard covariates + library batch; results are combined across tissues by de_meta.py.

  DE_METRIC=dna_burden_z  -> cumulative mutated 5S copies (genetic dose)
  DE_METRIC=rna_excess_z  -> variant-5S expression above non-carrier background
  (SHORT also maps log_*/cn_z/mut_z variants; scope tag uses the short name.)

Design = ~ SMRIN + SMTSISCH + gPC1-5 + SEX + DTHHRDY [+ batch] + metric ; Wald test on 'metric'.
The contrast is built as a unit vector selecting the 'metric' column of the fitted design matrix.

Env: DE_METRIC selects the predictor; DE_BATCH (default seqbatch) = explicit library-batch term;
     DE_ANCESTRY=WHITE restricts to RACE==3 donors (scope gets a 'W'); DE_PERMUTE=<seed> shuffles the
     dosage label among donors for a permutation null (tag gets _perm<seed>); DE_MIN_TISSUE_N (default
     300) minimum samples per tissue; DE_RESUME skips existing outputs; DE_CPUS.
Output de_v2/out/de_<scope>_<tissue>_cont.tsv ; scope = cont[W]_{dna|rna}[_perm<seed>].
Usage: DE_METRIC=... [DE_ANCESTRY=WHITE] [DE_PERMUTE=N] de_continuous.py [--shard i/k] [--tissue T]
"""
import sys, os, re, warnings, numpy as np, pandas as pd
warnings.simplefilter("ignore")
import de_common as C
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

MIN_CASE = 20   # minimum samples in a tissue to attempt a fit
N_CPUS = int(os.environ.get("DE_CPUS", "4"))
OUT = f"{C.ROOT}/de_v2/out"; os.makedirs(OUT, exist_ok=True)
COV = ["SMRIN", "SMTSISCH", "gPC1", "gPC2", "gPC3", "gPC4", "gPC5"]  # continuous covariates
METRIC = os.environ.get("DE_METRIC", "dna_burden_z")  # which donor_metrics.tsv column is the predictor
# donor_metrics.tsv column -> short token used in the output scope tag (de_<scope>_<tissue>_cont.tsv)
SHORT = {"dna_burden_z": "dna", "rna_excess_z": "rna", "log_dna_burden_z": "dnalog", "log_rna_excess_z": "rnalog",
         "cn_z": "cn", "mut_z": "mut",
         "dna_lowIE_z": "lowIE", "dna_highIE_z": "highIE"}  # DNA copy-load of low- vs high-internal-promoter (I/E) variant classes


def run_tissue(A, mmap, tissue, tag_prefix):
    """Fit the continuous-dosage DE for ONE tissue and write de_<tag>.tsv. mmap: donor -> dosage value."""
    tag = f"{tag_prefix}_{re.sub(r'[^A-Za-z0-9]+','_',tissue)}_cont"
    if os.environ.get("DE_RESUME") and os.path.exists(f"{OUT}/de_{tag}.tsv"):
        return  # DE_RESUME: skip tissues already written
    sub = A[A.obs.SMTSD == tissue].copy()
    sub.obs["m"] = sub.obs.donor.map(mmap)          # attach each sample's donor dosage metric
    sub = sub[sub.obs.m.notna()].copy()             # require a dosage value (donors present in the metric table)
    need = COV + ["SEX", "DTHHRDY"]
    ok = sub.obs[need].notna().all(1) & (sub.obs.SEX != "<NA>")  # complete-covariate samples only
    sub = sub[ok.values].copy()
    if sub.n_obs < MIN_CASE:
        return
    cnt = sub.X.toarray().astype(int) if not isinstance(sub.X, np.ndarray) else sub.X.astype(int)
    kf = C.gene_filter(cnt, 10, 0.5); cnt = cnt[:, kf]; gn = sub.var_names[kf]  # keep genes >=10 in >=50% samples
    ob = sub.obs.copy()
    # assemble the pydeseq2 metadata frame (numeric covariates + categorical SEX/Hardy + optional batch + metric)
    md = pd.DataFrame(index=ob.index)
    for c in COV:
        md[c] = pd.to_numeric(ob[c], errors="coerce").astype(float)
    md["SEX"] = ob.SEX.astype(str).astype("category")
    md["DTHHRDY"] = ob.DTHHRDY.astype(str).astype("category")
    # explicit library-batch term (default seqbatch = SMGEBTCH); merge rare levels (<5 samples) into
    # 'other' so the term stays estimable, and only include it if >=2 levels survive.
    bterm = ""
    bc = os.environ.get("DE_BATCH", "seqbatch")
    if bc:
        b = ob[bc].astype(str).fillna("NA"); vc = b.value_counts()
        b = b.where(~b.isin(vc[vc < 5].index), "other")
        if b.nunique() >= 2:
            md["batch"] = pd.Categorical(b); bterm = " + batch"
    md["metric"] = ob.m.astype(float)   # already standardized in donor_metrics.tsv
    cdf = pd.DataFrame(cnt, index=ob.index, columns=gn)  # samples x genes raw counts
    design = f"~ SMRIN + SMTSISCH + gPC1 + gPC2 + gPC3 + gPC4 + gPC5 + SEX + DTHHRDY{bterm} + metric"
    try:
        dds = DeseqDataSet(counts=cdf, metadata=md, design=design, quiet=True, n_cpus=N_CPUS)
        dds.deseq2()
        # continuous predictor has no factor levels -> Wald-test it via a contrast vector that picks
        # out the single 'metric' column of the fitted design matrix (log2FC per 1 SD of dosage).
        cols = list(dds.obsm["design_matrix"].columns)
        vec = np.array([1.0 if c == "metric" else 0.0 for c in cols])
        st = DeseqStats(dds, contrast=vec, quiet=True); st.summary()
        r = st.results_df.copy()  # per-gene: baseMean, log2FoldChange, lfcSE, stat, pvalue, padj
    except Exception as e:
        print(f"  FAIL {tag_prefix} {tissue}: {type(e).__name__}: {str(e)[:100]}"); return
    # provenance columns for de_meta.py; nB/k_sv unused here (single continuous term, no reference arm/SVs)
    r["scope"] = tag_prefix; r["tissue"] = tissue; r["contrast"] = "cont"
    r["nA"] = sub.n_obs; r["nB"] = 0; r["k_sv"] = 0; r.index.name = "ensg"
    r.to_csv(f"{OUT}/de_{tag}.tsv", sep="\t")
    print(f"  {tissue:32s} n={sub.n_obs} genes={len(r)} padj<0.1={(r.padj<0.1).sum()}", flush=True)


def main():
    # load raw counts + covariates and attach per-sample covariates (RIN/ischemic/PCs/sex/Hardy/batch)
    meta = C.load_meta(); A = C.load_counts()
    A.obs = A.obs.join(meta.drop(columns=["donor"]), how="left"); A.obs["donor"] = A.obs["donor"].astype(str)
    W = ""
    # DE_ANCESTRY=WHITE: restrict to the largest homogeneous ancestry (RACE==3) as an ancestry-confound check
    if os.environ.get("DE_ANCESTRY", "ALL") == "WHITE":
        SP = pd.read_csv(C.SUBJPHENO, sep="\t", low_memory=False).set_index("SUBJID")
        white = set(SP.index[SP.RACE == 3]); A = A[A.obs.donor.isin(white)].copy(); W = "W"
        print(f"[ANCESTRY=WHITE] {A.n_obs} samples / {A.obs.donor.nunique()} donors", flush=True)
    # per-donor dosage predictor: pick the DE_METRIC column (non-carriers already encoded as 0)
    dm = pd.read_csv(f"{C.ROOT}/de_v2/donor_metrics.tsv", sep="\t").set_index("donor")
    mmap = dm[METRIC].to_dict()
    tag_prefix = f"cont{W}_{SHORT[METRIC]}"
    if os.environ.get("DE_PERMUTE"):   # null: shuffle the dosage label among donors (preserves distribution + covariates)
        # break the donor<->dosage link (values permuted, donor set + covariates untouched) to calibrate
        # the false-positive rate / genomic inflation; output tag gets a _perm<seed> suffix.
        seed = int(os.environ["DE_PERMUTE"]); rng = np.random.default_rng(seed)
        donors = list(mmap.keys()); vals = np.array([mmap[d] for d in donors]); rng.shuffle(vals)
        mmap = dict(zip(donors, vals)); tag_prefix += f"_perm{seed}"
        print(f"[PERMUTE seed={seed}] shuffled metric among {len(donors)} donors", flush=True)
    # tissues largest-first; keep only those with >=DE_MIN_TISSUE_N samples (continuous fit needs power)
    tissues = sorted(A.obs.SMTSD.dropna().unique(), key=lambda t: -(A.obs.SMTSD == t).sum())
    mn = int(os.environ.get("DE_MIN_TISSUE_N", "300"))
    tissues = [t for t in tissues if (A.obs.SMTSD == t).sum() >= mn]
    if "--tissue" in sys.argv:  # run a single named tissue
        tissues = [sys.argv[sys.argv.index("--tissue") + 1]]
    if "--shard" in sys.argv:   # "i/k": this process takes every k-th tissue starting at i (parallelism)
        i, k = map(int, sys.argv[sys.argv.index("--shard") + 1].split("/")); tissues = tissues[i::k]
    print(f"[{tag_prefix}] {len(tissues)} tissues, metric={METRIC}", flush=True)
    for t in tissues:
        run_tissue(A, mmap, t, tag_prefix)


if __name__ == "__main__":
    main()
