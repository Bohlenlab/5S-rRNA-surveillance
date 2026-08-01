# -----------------------------------------------------------------------------
# de_pertissue.py — per-tissue DESeq2 negative-binomial GLM for 5S-variant donor
# contrasts (carrier/expresser/tertile/per-variant), with covariate + SV correction.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""Per-tissue DESeq2 for binary 5S-variant contrasts (carrier-vs-non-carrier etc.).

Fits a negative-binomial GLM (pydeseq2) separately per tissue, comparing gene expression between
two donor groups defined by 5S rRNA variant status, adjusting for the standard covariate set plus a
batch/surrogate-variable correction. One output TSV is written per (test, tissue) and later combined
across tissues by de_meta.py.

Usage:
  de_pertissue.py <contrast> [--scan] [--tissue "Whole Blood"]   # aggregate contrast
  de_pertissue.py pervariant [--scan] [--min-carrier 25]          # per-variant DNA carrier-vs-noncarrier
  contrast in {dna, rna, 3group, tertile, pervariant, aggregate}
    - dna/rna/3group/tertile : binary/multi-level group columns from groups_donor.tsv (see CONTRASTS)
    - aggregate              : runs rna, 3group, tertile, dna in one pass
    - pervariant             : one carrier-vs-non-carrier test per single GEN variant (>= --min-carrier)
  --tissue T : run only tissue T (default = every tissue with enough donors per arm)

Design (per tissue) =
  ~ SMRIN + SMTSISCH + gPC1..gPC5 + SEX + DTHHRDY [+ batch] [+ SV1..SVk] + group
Batch/latent correction is one of (see the per-test block below):
  - explicit library batch (DE_BATCH=seqbatch|seqbin|nabatch), the current/preferred mode; and/or
  - surrogate variables (DE_SV_MODE): 'null' (default) = residual PCs after the COVARIATES-only design
    (batch-capturing, kept even where batch correlates with group); 'ortho' = classic RUVr (residual
    after the FULL design incl. group -> group-orthogonal); 'none' = skip when a batch term is present.
The Wald test extracts the group A-vs-B coefficient (B = reference).

Env: DE_ANCESTRY=WHITE restricts to RACE==3 donors; DE_BATCH / DE_SV_MODE select correction mode;
     DE_MIN_TISSUE_N caps to well-powered tissues; DE_RESUME skips already-written outputs; DE_CPUS.
Writes de_v2/out/de_<tag>.tsv (tag = <scope>_<tissue>_<Alvl>vs<Blvl>), one per test/tissue.
"""
import sys, os, re, warnings, numpy as np, pandas as pd, anndata as ad
warnings.simplefilter("ignore")
import de_common as C
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# Aggregate contrasts: each maps a groups_donor.tsv column ('col') to the ordered list of pairwise
# tests to run on it, with 'ref' = the baseline level. 'tests' entries are (A_level, B_level) where
# B is the DESeq2 reference. dna = DNA carrier vs non-carrier; rna = RNA expresser vs non-expresser;
# 3group = expresser/silent/NC three-way; tertile = high/low dosage tertiles vs NC.
CONTRASTS = {
    "dna":     {"col": "g_dna",     "ref": "NC",      "tests": [("carrier", "NC")]},
    "rna":     {"col": "g_rna",     "ref": "nonexpr", "tests": [("expr", "nonexpr")]},
    "3group":  {"col": "g_3group",  "ref": "NC",      "tests": [("expr", "silent"), ("expr", "NC"), ("silent", "NC")]},
    "tertile": {"col": "g_tertile", "ref": "NC",      "tests": [("high", "NC"), ("low", "NC"), ("high", "low")]},
}
MIN_CASE = 10    # min donors in the smaller arm for a (tissue, test) to be fit
MIN_CARRIER = 25  # per-variant: min carrier-donors for a variant to be tested
N_CPUS = int(os.environ.get("DE_CPUS", "4"))
OUT = f"{C.ROOT}/de_v2/out"; os.makedirs(OUT, exist_ok=True)


def load_all():
    """Load donor group labels + raw counts, attach per-sample covariates, optional ancestry filter.
    Returns (AnnData A with covariate columns in .obs, labels frame lab from groups_donor.tsv)."""
    lab = pd.read_csv(f"{C.ROOT}/de_v2/groups_donor.tsv", sep="\t")
    meta = C.load_meta()
    A = C.load_counts()
    A.obs = A.obs.join(meta.drop(columns=["donor"]), how="left")  # add RIN/ischemic/PCs/sex/Hardy/batch per SAMPID
    A.obs["donor"] = A.obs["donor"].astype(str)
    # DE_ANCESTRY=WHITE -> stratify to the largest homogeneous ancestry group (RACE==3) to remove
    # population-structure confounding from the DNA carrier-vs-noncarrier contrast.
    if os.environ.get("DE_ANCESTRY", "ALL") == "WHITE":
        SP = pd.read_csv(C.SUBJPHENO, sep="\t", low_memory=False).set_index("SUBJID")
        white = set(SP.index[SP.RACE == 3])
        A = A[A.obs.donor.isin(white)].copy()
        print(f"[ANCESTRY=WHITE] restricted to {A.n_obs} samples / {A.obs.donor.nunique()} donors", flush=True)
    return A, lab


def pervariant_labels(min_carrier=MIN_CARRIER):
    """donor->{'carrier','ref'} maps, one per variant with >=min_carrier carriers.
    ref = every donor NOT carrying that variant (standard single-variant reference)."""
    d = pd.read_csv(C.RNAVAF, sep="\t")
    cc = d.groupby("variant").donor.nunique()
    vs = sorted(cc[cc >= min_carrier].index, key=lambda v: -cc[v])
    car = {v: set(d[d.variant == v].donor) for v in vs}
    return vs, car


def run_label(A, donor_grp, tests, tag_prefix, only=None):
    """Generic per-tissue DESeq2 over a donor->group label. tests=[(A_lvl,B_lvl)] (B=reference).
    Selects which tissues to fit (enough donors per arm), then delegates each to _run_tissue_label."""
    levels = sorted({x for t in tests for x in t})
    A.obs["group"] = A.obs.donor.map(donor_grp)
    # tissue x group-level count table -> keep tissues where at least one test has >=MIN_CASE per arm
    tab = A.obs.groupby("SMTSD")["group"].value_counts().unstack(fill_value=0)
    tab = tab[[c for c in levels if c in tab.columns]]
    qual = [t for t in tab.index if any(min(tab.loc[t].get(a, 0), tab.loc[t].get(b, 0)) >= MIN_CASE for a, b in tests)]
    minn = int(os.environ.get("DE_MIN_TISSUE_N", "0"))  # cap to well-powered tissues (meta is large-tissue dominated)
    if minn:
        tn = A.obs.SMTSD.value_counts()
        qual = [t for t in qual if tn.get(t, 0) >= minn]
    tissues = [only] if only else qual
    out = []
    for tissue in tissues:
        out += _run_tissue_label(A, donor_grp, tests, tissue, tag_prefix)
    return out


def _run_tissue_label(A, donor_grp, tests, tissue, tag_prefix):
    """Fit every requested contrast for ONE tissue; write de_<tag>.tsv per test and return the frames."""
    sub = A[A.obs.SMTSD == tissue].copy()
    sub.obs["group"] = sub.obs.donor.map(donor_grp)
    levels = sorted({x for t in tests for x in t})
    sub = sub[sub.obs.group.isin(levels)].copy()  # keep only samples in one of the contrast levels
    # drop samples missing any model covariate (DESeq2 needs a complete design matrix)
    need = ["SMRIN", "SMTSISCH", "gPC1", "gPC2", "gPC3", "gPC4", "gPC5", "SEX", "DTHHRDY"]
    ok = sub.obs[need].notna().all(1) & (sub.obs.SEX != "<NA>")
    sub = sub[ok.values].copy()
    if sub.n_obs < 40:  # too few samples for a stable per-tissue fit
        return []
    # one sample per donor per tissue is typical; keep all samples (rare dup tissue aliquots)
    counts = sub.X.toarray().astype(int) if not isinstance(sub.X, np.ndarray) else sub.X.astype(int)
    # expression filter: keep genes with >=10 counts in >=50% of samples (drops low/undetected genes)
    keep = C.gene_filter(counts, min_count=10, min_frac=0.5)
    counts = counts[:, keep]; genes = sub.var_names[keep]
    obs = sub.obs.copy()
    results = []
    for A_lvl, B_lvl in tests:
        nA = (obs.group == A_lvl).sum(); nB = (obs.group == B_lvl).sum()
        if min(nA, nB) < MIN_CASE:  # not enough donors in one arm of this particular test
            continue
        tag = f"{tag_prefix}_{re.sub(r'[^A-Za-z0-9]+','_',tissue)}_{A_lvl}vs{B_lvl}"
        if os.environ.get("DE_RESUME") and os.path.exists(f"{OUT}/de_{tag}.tsv"):
            continue  # DE_RESUME: skip tests whose output already exists
        # restrict to the two levels of this test and re-apply the gene filter on that subset
        m = obs.group.isin([A_lvl, B_lvl]).values
        cnt = counts[m]; ob = obs[m].copy()
        kf = C.gene_filter(cnt, min_count=10, min_frac=0.5)
        cnt = cnt[:, kf]; gn = genes[kf]
        ob["group"] = pd.Categorical(ob.group, categories=[B_lvl, A_lvl])  # B = reference
        # --- surrogate variables ---
        # DE_SV_MODE: 'null' (default, GTEx/PEER-style) estimates SVs from the COVARIATES-ONLY residual,
        #   so they capture extraction-batch structure even where it is correlated with group, and are
        #   symmetric across contrasts. 'ortho' = classic RUVr (residual incl. group -> group-orthogonal).
        gd = pd.get_dummies(ob.group, drop_first=True)
        sv_mode = os.environ.get("DE_SV_MODE", "null")
        # optional explicit batch covariate (DE_BATCH=seqbin|seqbatch|nabatch); merge rare within-tissue
        # levels (<MINB samples) into 'other' so the term stays estimable and not group-nested.
        batch_col = os.environ.get("DE_BATCH", "")
        bterm = ""
        if batch_col:
            MINB = int(os.environ.get("DE_BATCH_MINSAMP", "5"))
            b = ob[batch_col].astype(str).fillna("NA")
            vc = b.value_counts(); small = vc[vc < MINB].index
            b = b.where(~b.isin(small), "other")
            if b.nunique() >= 2:
                ob["batch"] = pd.Categorical(b)
                bterm = " + batch"
        # SVs skipped only when explicitly disabled AND we have working correction (explicit batch
        # succeeded or no batch requested). If a batch was requested but degenerated (<2 levels),
        # fall back to null-SVs so the tissue is still batch-corrected.
        if sv_mode == "none" and (not batch_col or bterm):
            k = 0
        else:
            gd_for_sv = gd if sv_mode == "ortho" else gd.iloc[:, :0]
            D, _ = C.build_design(ob, gd_for_sv)
            kmax = int(min(15, max(1, cnt.shape[0] // 25)))
            SV, k = C.compute_ruvr_svs(cnt, D, k_max=kmax)
            for j in range(k):
                ob[f"SV{j+1}"] = SV[:, j]
        # --- DESeq2 ---
        # Build the pydeseq2 metadata frame: numeric covariates, categorical SEX/Hardy, optional batch,
        # any estimated SVs, and the group column (last term -> the coefficient of interest).
        md = pd.DataFrame(index=ob.index)
        for c in ["SMRIN", "SMTSISCH", "gPC1", "gPC2", "gPC3", "gPC4", "gPC5"]:
            md[c] = pd.to_numeric(ob[c], errors="coerce").astype(float)
        md["SEX"] = ob.SEX.astype(str).astype("category")
        md["DTHHRDY"] = ob.DTHHRDY.astype(str).astype("category")
        if bterm:
            md["batch"] = ob["batch"]
        for j in range(k):
            md[f"SV{j+1}"] = ob[f"SV{j+1}"].astype(float)
        md["group"] = ob.group
        cdf = pd.DataFrame(cnt, index=ob.index, columns=gn)  # samples x genes raw counts
        svterms = "".join(f" + SV{j+1}" for j in range(k))
        design = f"~ SMRIN + SMTSISCH + gPC1 + gPC2 + gPC3 + gPC4 + gPC5 + SEX + DTHHRDY{bterm}{svterms} + group"
        try:
            # fit the NB GLM, then Wald-test the group A-vs-B contrast (log2FC of A relative to B)
            dds = DeseqDataSet(counts=cdf, metadata=md, design=design, quiet=True, n_cpus=N_CPUS)
            dds.deseq2()
            st = DeseqStats(dds, contrast=["group", A_lvl, B_lvl], quiet=True)
            st.summary()
            r = st.results_df.copy()  # per-gene: baseMean, log2FoldChange, lfcSE, stat, pvalue, padj
        except Exception as e:
            print(f"  FAIL {tissue} {A_lvl}vs{B_lvl}: {type(e).__name__}: {str(e)[:120]}")
            continue
        # annotate provenance columns that de_meta.py groups/keys on, then write one TSV per test
        r["scope"] = tag_prefix; r["tissue"] = tissue; r["contrast"] = f"{A_lvl}vs{B_lvl}"
        r["nA"] = nA; r["nB"] = nB; r["k_sv"] = k
        r.index.name = "ensg"
        tag = f"{tag_prefix}_{re.sub(r'[^A-Za-z0-9]+','_',tissue)}_{A_lvl}vs{B_lvl}"
        r.to_csv(f"{OUT}/de_{tag}.tsv", sep="\t")
        nsig = (r.padj < 0.1).sum()
        print(f"  {tissue:32s} {A_lvl}vs{B_lvl}: nA={nA} nB={nB} k={k} genes={len(r)} padj<0.1={nsig}")
        results.append(r)
    return results


def main():
    contrast = sys.argv[1]  # first positional arg: contrast name (see CONTRASTS / 'aggregate' / 'pervariant')
    only = sys.argv[sys.argv.index("--tissue") + 1] if "--tissue" in sys.argv else None
    A, lab = load_all()
    if contrast == "aggregate":
        # run all four aggregate contrasts in one process, reusing the loaded counts/covariates
        for ct in ["rna", "3group", "tertile", "dna"]:
            spec = CONTRASTS[ct]
            donor_grp = lab.set_index("donor")[spec["col"]]  # donor -> level from that groups_donor column
            print(f"=== contrast={ct} ===", flush=True)
            run_label(A, donor_grp, spec["tests"], ct, only=only)
        return
    if contrast == "pervariant":
        # per-variant DNA carrier-vs-non-carrier: one binary test for each sufficiently common GEN variant
        mc = int(sys.argv[sys.argv.index("--min-carrier") + 1]) if "--min-carrier" in sys.argv else MIN_CARRIER
        vs, car = pervariant_labels(mc)
        if "--shard" in sys.argv:  # "i/k": process variants where index % k == i (parallel processes)
            i, k = (int(x) for x in sys.argv[sys.argv.index("--shard") + 1].split("/"))
            vs = [v for n, v in enumerate(vs) if n % k == i]
        all_donors = set(A.obs.donor)
        pfx = "varW_" if os.environ.get("DE_ANCESTRY", "ALL") == "WHITE" else "var_"  # separate White outputs
        print(f"=== pervariant DNA ({pfx}): {len(vs)} variants (this shard) with >= {mc} carriers ===")
        for v in vs:
            cset = car[v]
            # carrier = donor carries variant v; ref = every other donor (single-variant reference)
            donor_grp = pd.Series({d: ("carrier" if d in cset else "ref") for d in all_donors})
            print(f"-- {v} (carriers={len(cset & all_donors)}) --")
            run_label(A, donor_grp, [("carrier", "ref")], f"{pfx}{v}", only=only)
        return
    # single named aggregate contrast (dna / rna / 3group / tertile)
    spec = CONTRASTS[contrast]
    donor_grp = lab.set_index("donor")[spec["col"]]
    print(f"=== contrast={contrast} ===")
    run_label(A, donor_grp, spec["tests"], contrast, only=only)


if __name__ == "__main__":
    main()
