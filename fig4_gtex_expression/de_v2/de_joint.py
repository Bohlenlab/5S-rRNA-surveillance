# -----------------------------------------------------------------------------
# de_joint.py — joint dosage model: per-tissue DESeq2 fitting the DNA mutation-load
# and RNA-expression 5S dosage axes together as partial coefficients.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""Joint dosage model: fit the DNA mutation-load and RNA-expression dosage axes simultaneously.

Marginal models fit each dosage axis alone and cannot separate them because they are correlated.
This script fits BOTH continuous, whole-cohort axes simultaneously so each reported effect is a
PARTIAL coefficient (that axis adjusted for the other):
  mut         = dna_burden_z  -- genetic dose (mutant 5S copies carried), from donor_metrics.tsv
  rna         = rna_excess_z  -- expression dose (variant 5S made above background)
Non-carriers enter with dosage 0 (whole-cohort, not carrier-restricted). The two donor-level axes
are near-orthogonal, so the joint fit is well-conditioned.

MODEL (per tissue, pydeseq2 negative-binomial GLM):
  expression ~ SMRIN + SMTSISCH + gPC1-5 + SEX + DTHHRDY + batch + mut + rna
  A separate Wald contrast is extracted for 'mut' and for 'rna' from the one shared fit.
INPUT:   de_v2/donor_metrics.tsv (mut_z, rna_excess_z) + counts/covariates via de_common.
OUTPUT:  de_v2/out/de_joint[W]_{mut,rna}_<tissue>_cont.tsv (scopes joint[W]_mut, joint[W]_rna,
         contrast='cont'); de_meta.py meta-analyzes each scope across tissues.
Usage: [DE_ANCESTRY=WHITE] de_joint.py [--shard i/k]"""
import sys, os, re, warnings, numpy as np, pandas as pd
warnings.simplefilter("ignore")
import de_common as C
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

MIN_CASE = 20          # min samples in a tissue to attempt a fit
N_CPUS = int(os.environ.get("DE_CPUS", "4"))
OUT = f"{C.ROOT}/de_v2/out"; os.makedirs(OUT, exist_ok=True)
COV = ["SMRIN", "SMTSISCH", "gPC1", "gPC2", "gPC3", "gPC4", "gPC5"]  # continuous covariates in the design


def run_tissue(A, mmap, rmap, tissue, W):
    """Fit the joint mut+rna model for one tissue and write both partial-coefficient TSVs.
    mmap/rmap: donor -> mut_z / rna_excess_z (whole cohort; non-carriers already 0 in donor_metrics)."""
    # tag template with a {} slot filled by 'mut'/'rna' for the two output files
    tag = f"joint{W}_{{}}_{re.sub(r'[^A-Za-z0-9]+','_',tissue)}_cont"
    if os.environ.get("DE_RESUME") and os.path.exists(f"{OUT}/de_{tag.format('mut')}.tsv") and os.path.exists(f"{OUT}/de_{tag.format('rna')}.tsv"):
        return
    sub = A[A.obs.SMTSD == tissue].copy()
    sub.obs["mmet"] = sub.obs.donor.map(mmap); sub.obs["rmet"] = sub.obs.donor.map(rmap)   # attach both dosage axes per donor
    sub = sub[sub.obs.mmet.notna() & sub.obs.rmet.notna()].copy()   # keep donors that have both metrics
    ok = sub.obs[COV + ["SEX", "DTHHRDY"]].notna().all(1) & (sub.obs.SEX != "<NA>")   # drop samples missing any covariate
    sub = sub[ok.values].copy()
    if sub.n_obs < MIN_CASE:
        return
    cnt = sub.X.toarray().astype(int) if not isinstance(sub.X, np.ndarray) else sub.X.astype(int)
    kf = C.gene_filter(cnt, 10, 0.5); cnt = cnt[:, kf]; gn = sub.var_names[kf]   # keep adequately-expressed genes
    ob = sub.obs.copy()
    md = pd.DataFrame(index=ob.index)
    for c in COV:
        md[c] = pd.to_numeric(ob[c], errors="coerce").astype(float)
    md["SEX"] = ob.SEX.astype(str).astype("category"); md["DTHHRDY"] = ob.DTHHRDY.astype(str).astype("category")
    # explicit library batch (SMGEBTCH); levels with <5 samples collapsed to 'other'; skipped if <2 levels
    bterm = ""
    bc = os.environ.get("DE_BATCH", "seqbatch")
    if bc:
        b = ob[bc].astype(str).fillna("NA"); vc = b.value_counts(); b = b.where(~b.isin(vc[vc < 5].index), "other")
        if b.nunique() >= 2:
            md["batch"] = pd.Categorical(b); bterm = " + batch"
    md["mut"] = ob.mmet.astype(float); md["rna"] = ob.rmet.astype(float)   # both already z-scored in donor_metrics.tsv
    cdf = pd.DataFrame(cnt, index=ob.index, columns=gn)
    design = f"~ SMRIN + SMTSISCH + gPC1 + gPC2 + gPC3 + gPC4 + gPC5 + SEX + DTHHRDY{bterm} + mut + rna"
    try:
        dds = DeseqDataSet(counts=cdf, metadata=md, design=design, quiet=True, n_cpus=N_CPUS)
        dds.deseq2()   # one shared fit with both dosage terms; partials extracted below
        cols = list(dds.obsm["design_matrix"].columns)
        # extract each PARTIAL coefficient from the single joint fit: a contrast vector that is 1 on the
        # target column ('mut' or 'rna') and 0 elsewhere isolates that axis adjusted for the other + covariates.
        for coef in ["mut", "rna"]:
            vec = np.array([1.0 if c == coef else 0.0 for c in cols])
            st = DeseqStats(dds, contrast=vec, quiet=True); st.summary(); r = st.results_df.copy()
            # tag with meta keys: scope encodes which partial (joint[W]_mut / joint[W]_rna); nB/k_sv=0 (uniform schema)
            r["scope"] = f"joint{W}_{coef}"; r["tissue"] = tissue; r["contrast"] = "cont"; r["nA"] = sub.n_obs; r["nB"] = 0; r["k_sv"] = 0; r.index.name = "ensg"
            r.to_csv(f"{OUT}/de_{tag.format(coef)}.tsv", sep="\t")
        print(f"  {tissue:30s} n={sub.n_obs} genes={len(gn)}", flush=True)
    except Exception as e:
        print(f"  FAIL {tissue}: {type(e).__name__}: {str(e)[:90]}")


def main():
    meta = C.load_meta(); A = C.load_counts()
    # attach per-sample covariates/batch; obs.donor is the GTEX-xxxx key that maps to the dosage metrics
    A.obs = A.obs.join(meta.drop(columns=["donor"]), how="left"); A.obs["donor"] = A.obs["donor"].astype(str)
    W = ""
    if os.environ.get("DE_ANCESTRY", "ALL") == "WHITE":   # within-White (RACE==3) ancestry sensitivity; tags scope with 'W'
        SP = pd.read_csv(C.SUBJPHENO, sep="\t", low_memory=False).set_index("SUBJID")
        white = set(SP.index[SP.RACE == 3]); A = A[A.obs.donor.isin(white)].copy(); W = "W"
    # the two whole-cohort dosage axes, precomputed + z-scored by build_metrics.py (non-carriers = 0)
    dm = pd.read_csv(f"{C.ROOT}/de_v2/donor_metrics.tsv", sep="\t").set_index("donor")
    mmap = dm["mut_z"].to_dict(); rmap = dm["rna_excess_z"].to_dict()
    tissues = sorted(A.obs.SMTSD.dropna().unique(), key=lambda t: -(A.obs.SMTSD == t).sum())   # largest tissue first
    mn = int(os.environ.get("DE_MIN_TISSUE_N", "300")); tissues = [t for t in tissues if (A.obs.SMTSD == t).sum() >= mn]
    if "--shard" in sys.argv:   # optional strided split of tissues across k parallel workers
        i, k = map(int, sys.argv[sys.argv.index("--shard") + 1].split("/")); tissues = tissues[i::k]
    print(f"[joint{W}] {len(tissues)} tissues, mut+rna partial coefficients", flush=True)
    for t in tissues:
        run_tissue(A, mmap, rmap, t, W)


if __name__ == "__main__":
    main()
