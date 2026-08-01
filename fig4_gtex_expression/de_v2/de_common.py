# -----------------------------------------------------------------------------
# de_common.py — shared layer for the GTEx bulk RNA-seq differential-expression
# pipeline over 5S-variant donor groups: paths, covariates, and surrogate variables.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""de_v2 shared layer — GTEx bulk RNA-seq DE pipeline for 5S-variant donor groups.

Covariate design:
  design = ~ gPC1 + gPC2 + gPC3 + SEX + SMRIN + DTHHRDY + SMTSISCH + (k RUVr SVs) + group
  - extraction batch (SMNABTCH) contributes most hidden structure but has many near-singleton
    levels, so it cannot be modelled as explicit dummies and is absorbed by RUVr surrogate variables.
  - genotype PCs adjust for ancestry, since low-frequency 5S variant carriage is correlated with ancestry.
RUVr SVs = top PCs of residuals after the FULL design (covariates+group), which makes them
group-orthogonal by construction and symmetric across contrasts.
"""
import os, numpy as np, pandas as pd, anndata as ad

# Paths are env-overridable (DE_* variables).
ROOT = os.environ.get("DE_ROOT", os.environ.get("FIVES_DATA", "data"))
H5AD = os.environ.get("DE_H5AD", os.environ.get("FIVES_GTEX_COUNTS", "GTEx_v11_bulk_gene_counts.h5ad"))
RNAVAF = os.environ.get("DE_RNAVAF", f"{ROOT}/results/eqtl/extreme/donor_variant_rnavaf.tsv")
SAMPLEATTR = os.environ.get("DE_SAMPLEATTR", f"{ROOT}/metadata/GTEx_Analysis_2025-08-22_v11_Annotations_SampleAttributesDS.txt")
SUBJPHENO = os.environ.get("DE_SUBJPHENO", f"{ROOT}/eqtl_inputs/SubjectPhenotypesDS.txt")
PCS = os.environ.get("DE_PCS", f"{ROOT}/results/eqtl_inputs/genotype_20PCs.eigenvec.txt")

# gene-region functional 5S variant set
GEN = ["687G", "701C", "725C", "726T", "730G", "733T", "734T", "743G"]

# measured covariates kept in every design (5 genotype PCs adjust for ancestry, since 5S variant
# carriage is correlated with ancestry; RUVr SVs are group-orthogonal and do not control ancestry).
CONT_COV = ["SMRIN", "SMTSISCH", "gPC1", "gPC2", "gPC3", "gPC4", "gPC5"]
CAT_COV = ["SEX", "DTHHRDY"]


def load_meta():
    """Per-SAMPID covariate frame: tissue + RIN + ischemic + sex + Hardy death class + genotype PC1-3."""
    SA = pd.read_csv(SAMPLEATTR, sep="\t", low_memory=False).set_index("SAMPID")[["SMTSD", "SMRIN", "SMTSISCH"]]
    donmap = pd.Series(["-".join(s.split("-")[:2]) for s in SA.index], index=SA.index)
    SP = pd.read_csv(SUBJPHENO, sep="\t", low_memory=False).set_index("SUBJID")
    SA["SEX"] = donmap.map(SP["SEX"]).astype("Int64").astype(str)
    SA["DTHHRDY"] = donmap.map(SP["DTHHRDY"]).fillna(-1).astype(int).astype(str)  # Hardy scale; -1 = unknown
    pc = pd.read_csv(PCS, sep="\t")
    pc["donor"] = pc.IID.str.extract(r'(GTEX-[A-Z0-9]+)')
    pc = pc.drop_duplicates("donor").set_index("donor")
    for i in (1, 2, 3, 4, 5):
        SA[f"gPC{i}"] = donmap.map(pc[f"PC{i}"])
    SA["donor"] = donmap
    # explicit batch annotations (for DE_BATCH mode): the fine library batch + a coarse temporal
    # clustering of it = sequencing date binned into NBINS quantiles (groups LCSETs sequenced together).
    sa_full = pd.read_csv(SAMPLEATTR, sep="\t", low_memory=False).set_index("SAMPID")
    SA["seqbatch"] = sa_full["SMGEBTCH"].reindex(SA.index)        # library/sequencing batch (LCSET)
    SA["nabatch"] = sa_full["SMNABTCH"].reindex(SA.index)         # nucleic-acid isolation batch
    sd = pd.to_datetime(sa_full["SMGEBTCHD"].reindex(SA.index), errors="coerce")
    nb = int(os.environ.get("DE_BATCH_NBINS", "24"))
    ok = sd.notna()
    SA["seqbin"] = np.nan
    SA.loc[ok, "seqbin"] = pd.qcut(sd[ok].astype("int64"), nb, labels=False, duplicates="drop").astype(float)
    return SA


def load_counts():
    """Raw-count AnnData, version-stripped Ensembl, duplicate genes dropped. (DESeq2 wants raw ints.)"""
    A = ad.read_h5ad(H5AD)
    A.var["ensg"] = [v.split(".")[0] for v in A.var_names]
    A = A[:, ~A.var.ensg.duplicated()].copy()
    A.var_names = A.var.ensg.values
    return A


def gene_filter(counts, min_count=10, min_frac=0.5):
    """Keep genes with >=min_count in >=min_frac of samples (within the tissue subset)."""
    X = counts if isinstance(counts, np.ndarray) else counts.toarray()
    keep = ((X >= min_count).mean(0) >= min_frac)
    return keep


def compute_ruvr_svs(raw_counts, design_mat, k_max=20, seed=0):
    """RUVr surrogate variables: regress log-CPM on the FULL design (covariates+group), SVD the
    residuals, return top-k PCs. k chosen by parallel analysis (permuted-eigenvalue null).
    raw_counts: samples x genes (ints). design_mat: samples x p (float, incl. intercept+group)."""
    X = raw_counts.astype(np.float64)
    lib = X.sum(1, keepdims=True); lib[lib == 0] = 1
    Y = np.log1p(X / lib * 1e4)
    Y = Y - Y.mean(0)
    D = design_mat
    # residuals after projecting out the full design
    beta, *_ = np.linalg.lstsq(D, Y, rcond=None)
    R = Y - D @ beta
    # SVD of residuals
    U, S, Vt = np.linalg.svd(R, full_matrices=False)
    ev = S ** 2
    # parallel analysis: permute genes within each column to break structure, get null eigenvalue spectrum
    rng = np.random.default_rng(seed)
    n_perm = 5
    null = np.zeros((n_perm, len(ev)))
    for p in range(n_perm):
        Rp = np.column_stack([rng.permutation(R[:, j]) for j in rng.choice(R.shape[1], min(2000, R.shape[1]), replace=False)])
        # scale null spectrum to same total residual variance per-feature
        sp = np.linalg.svd(Rp - Rp.mean(0), compute_uv=False) ** 2
        sp = sp * (ev.sum() / (R.shape[1] / Rp.shape[1]) / sp.sum())
        null[p, :len(sp)] = sp[:len(ev)]
    thr = null.mean(0)
    k = int(min(k_max, max(1, np.sum(ev[:k_max] > thr[:k_max]))))
    SV = U[:, :k] * S[:k]
    SV = (SV - SV.mean(0)) / (SV.std(0) + 1e-9)
    return SV, k


def build_design(obs, group_dummies, cont=CONT_COV, cat=CAT_COV):
    """Assemble a numeric design matrix (intercept + covariates + group dummies) for SV estimation.
    Missing continuous -> median impute; categorical -> one-hot drop-first. Returns (D, colnames)."""
    cols, names = [np.ones(len(obs))], ["intercept"]
    for c in cont:
        x = pd.to_numeric(obs[c], errors="coerce").values.astype(float)
        x = np.where(np.isnan(x), np.nanmedian(x), x)
        cols.append((x - x.mean()) / (x.std() + 1e-9)); names.append(c)
    for c in cat:
        d = pd.get_dummies(obs[c].astype(str), prefix=c, drop_first=True)
        for col in d.columns:
            cols.append(d[col].values.astype(float)); names.append(col)
    for col in group_dummies.columns:
        cols.append(group_dummies[col].values.astype(float)); names.append(col)
    return np.column_stack(cols), names
