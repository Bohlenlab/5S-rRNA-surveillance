# -----------------------------------------------------------------------------
# de_dose.py — within-carrier RNA-VAF dose model: per-tissue DESeq2 effect of how
# much a carried 5S variant is expressed (beta_RNA), per variant and for the aggregate.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""Within-carrier RNA-VAF dose model: per-gene effect of variant expression level (beta_RNA).

The marginal DNA carrier model confounds whether a donor carries variant V with how much of V
they express. This script conditions on carriage: for each variant V (and the GEN aggregate) it
restricts to CARRIERS of V and regresses expression on the continuous RNA-VAF dose of V. Because
every donor in the fit carries V, the 'dose' coefficient measures the effect of how much V is
expressed alone.

MODEL (per tissue, negative-binomial GLM via pydeseq2; same covariates + library batch as the DNA model):
  expression ~ SMRIN + SMTSISCH + gPC1-5 + SEX + DTHHRDY + batch + dose
  dose = standardized log10 of V's RNA-VAF; Wald test on the single 'dose' coefficient.

INPUTS:  counts + covariates via de_common (H5AD, GTEx metadata, genotype PCs); per-donor per-variant
         VAF table de_common.RNAVAF (donor_variant_rnavaf.tsv) supplies each carrier's VAF of V.
OUTPUT:  de_v2/out/de_dose[W]_{V|AGG}_<tissue>_dose.tsv  (one row per gene; scope=dose[W]_{V|AGG},
         contrast='dose'). de_meta.py inverse-variance meta-analyzes these across tissues.
         The DE_DOSE_SRC=wgs_vaf variant writes the parallel 'dnadose*' scope instead.
Usage: [DE_ANCESTRY=WHITE] [DE_DOSE_SRC=rna_vaf|wgs_vaf] de_dose.py [--shard i/k] [--min-carrier 25] [--tissue T]
"""
import sys, os, re, warnings, numpy as np, pandas as pd
warnings.simplefilter("ignore")
import de_common as C
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

MIN_CASE = 10          # min carrier-samples in a tissue to attempt a fit
MIN_CARRIER = 25       # min carriers of V (across tissues) for V to get its own dose model
N_CPUS = int(os.environ.get("DE_CPUS", "4"))
OUT = f"{C.ROOT}/de_v2/out"; os.makedirs(OUT, exist_ok=True)
COV = ["SMRIN", "SMTSISCH", "gPC1", "gPC2", "gPC3", "gPC4", "gPC5"]  # continuous covariates in the design


def run_tissue_dose(A, dosemap, tissue, tag_prefix):
    """Fit the within-carrier dose model for one variant (tag_prefix) in one tissue and write its TSV.
    dosemap: donor -> that donor's VAF of the variant (carriers only; non-carriers are absent/NaN)."""
    tag = f"{tag_prefix}_{re.sub(r'[^A-Za-z0-9]+','_',tissue)}_dose"
    if os.environ.get("DE_RESUME") and os.path.exists(f"{OUT}/de_{tag}.tsv"):
        return
    sub = A[A.obs.SMTSD == tissue].copy()
    sub.obs["dose_raw"] = sub.obs.donor.map(dosemap)
    sub = sub[sub.obs.dose_raw.notna()].copy()   # CARRIER RESTRICTION: keep only samples whose donor carries V
    need = COV + ["SEX", "DTHHRDY"]
    ok = sub.obs[need].notna().all(1) & (sub.obs.SEX != "<NA>")   # drop samples missing any covariate
    sub = sub[ok.values].copy()
    if sub.n_obs < MIN_CASE or sub.obs.dose_raw.nunique() < 4:   # too few carriers / too little dose spread to fit
        return
    cnt = sub.X.toarray().astype(int) if not isinstance(sub.X, np.ndarray) else sub.X.astype(int)
    kf = C.gene_filter(cnt, 10, 0.5); cnt = cnt[:, kf]; gn = sub.var_names[kf]   # keep adequately-expressed genes
    ob = sub.obs.copy()
    # dose predictor = standardized log10 VAF. log10 compresses the heavy VAF tail so the coefficient
    # is per-decade (fold-change) rather than per-VAF-unit; +1e-4 pseudocount avoids log10(0); then
    # z-score so the Wald effect is comparable across variants/tissues before meta-analysis.
    dose = np.log10(ob.dose_raw.values.astype(float) + 1e-4)
    ob["dose"] = (dose - dose.mean()) / (dose.std() + 1e-9)
    md = pd.DataFrame(index=ob.index)
    for c in COV:
        md[c] = pd.to_numeric(ob[c], errors="coerce").astype(float)
    md["SEX"] = ob.SEX.astype(str).astype("category")
    md["DTHHRDY"] = ob.DTHHRDY.astype(str).astype("category")
    # explicit library batch (SMGEBTCH, same correction as the DNA model); levels with <5 samples are
    # collapsed to 'other' so DESeq2 does not choke on near-singleton dummies. Skipped if <2 usable levels.
    bterm = ""
    bc = os.environ.get("DE_BATCH", "seqbatch")
    if bc:
        b = ob[bc].astype(str).fillna("NA"); vc = b.value_counts()
        b = b.where(~b.isin(vc[vc < 5].index), "other")
        if b.nunique() >= 2:
            md["batch"] = pd.Categorical(b); bterm = " + batch"
    md["dose"] = ob.dose.astype(float)
    cdf = pd.DataFrame(cnt, index=ob.index, columns=gn)
    design = f"~ SMRIN + SMTSISCH + gPC1 + gPC2 + gPC3 + gPC4 + gPC5 + SEX + DTHHRDY{bterm} + dose"
    try:
        dds = DeseqDataSet(counts=cdf, metadata=md, design=design, quiet=True, n_cpus=N_CPUS)
        dds.deseq2()
        # Wald test on 'dose' alone: build a contrast vector that is 1 on the dose column of the fitted
        # design matrix and 0 everywhere else, so results_df reports the dose log2FC / p adjusted for all covariates+batch.
        cols = list(dds.obsm["design_matrix"].columns)
        vec = np.array([1.0 if c == "dose" else 0.0 for c in cols])
        st = DeseqStats(dds, contrast=vec, quiet=True); st.summary()
        r = st.results_df.copy()
    except Exception as e:
        print(f"  FAIL {tag_prefix} {tissue}: {type(e).__name__}: {str(e)[:100]}"); return
    # annotate rows with the meta-analysis keys de_meta groups on (scope, tissue, contrast) and the
    # sample count (nA); nB/k_sv are 0 here (no second group, no RUVr SVs) but kept for a uniform schema.
    r["scope"] = tag_prefix; r["tissue"] = tissue; r["contrast"] = "dose"
    r["nA"] = sub.n_obs; r["nB"] = 0; r["k_sv"] = 0; r.index.name = "ensg"
    tag = f"{tag_prefix}_{re.sub(r'[^A-Za-z0-9]+','_',tissue)}_dose"
    r.to_csv(f"{OUT}/de_{tag}.tsv", sep="\t")
    print(f"  {tissue:32s} dose: n={sub.n_obs} genes={len(r)} padj<0.1={(r.padj<0.1).sum()}", flush=True)


def main():
    mc = int(sys.argv[sys.argv.index("--min-carrier") + 1]) if "--min-carrier" in sys.argv else MIN_CARRIER
    only = sys.argv[sys.argv.index("--tissue") + 1] if "--tissue" in sys.argv else None
    meta = C.load_meta(); A = C.load_counts()
    # attach per-sample covariates/batch to the count matrix; obs.donor is the GTEX-xxxx key for VAF mapping
    A.obs = A.obs.join(meta.drop(columns=["donor"]), how="left"); A.obs["donor"] = A.obs["donor"].astype(str)
    W = ""
    if os.environ.get("DE_ANCESTRY", "ALL") == "WHITE":  # ancestry-controlled within-carrier dose
        SP = pd.read_csv(C.SUBJPHENO, sep="\t", low_memory=False).set_index("SUBJID")
        white = set(SP.index[SP.RACE == 3]); A = A[A.obs.donor.isin(white)].copy(); W = "W"
        print(f"[ANCESTRY=WHITE] dose restricted to {A.n_obs} samples / {A.obs.donor.nunique()} donors", flush=True)
    d = pd.read_csv(C.RNAVAF, sep="\t")
    # DE_DOSE_SRC: rna_vaf (expressed fraction, default) or wgs_vaf (intragenomic DNA copy-fraction = genetic dose)
    SRC = os.environ.get("DE_DOSE_SRC", "rna_vaf")
    BASE = "dose" if SRC == "rna_vaf" else "dnadose"   # scope stem drives the output-file / meta grouping name
    # build one dose job per variant carried by >=mc donors, ordered most-carried first (best-powered first).
    cc = d.groupby("variant").donor.nunique(); vs = sorted(cc[cc >= mc].index, key=lambda v: -cc[v])
    g = d[d.variant.isin(C.GEN)]
    agg = g.groupby("donor")[SRC].max()   # GEN aggregate dose = a donor's max VAF over the functional GEN set
    # jobs = (scope_tag, donor->VAF map); AGG first, then each single variant. set_index('donor') makes the per-tissue .map lookup a carrier-only dose map.
    jobs = [(f"{BASE}{W}_AGG", agg)] + [(f"{BASE}{W}_{v}", d[d.variant == v].set_index("donor")[SRC]) for v in vs]
    if "--shard" in sys.argv:   # optional round-robin split of the variant jobs across k parallel workers
        i, k = (int(x) for x in sys.argv[sys.argv.index("--shard") + 1].split("/"))
        jobs = [j for n, j in enumerate(jobs) if n % k == i]
    tissues = A.obs.SMTSD.value_counts()
    big = list(tissues[tissues >= int(os.environ.get("DE_MIN_TISSUE_N", "300"))].index)   # tissues with enough samples to fit
    print(f"=== dose model: {len(jobs)} variants (this shard), {len(big)} tissues ===", flush=True)
    for tag_prefix, dosemap in jobs:
        print(f"-- {tag_prefix} (n_carriers={dosemap.notna().sum()}) --", flush=True)
        for t in ([only] if only else big):
            run_tissue_dose(A, dosemap, t, tag_prefix)


if __name__ == "__main__":
    main()
