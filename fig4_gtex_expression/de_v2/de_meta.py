# -----------------------------------------------------------------------------
# de_meta.py — inverse-variance cross-tissue meta-analysis and directional
# consistency of the per-tissue DESeq2 results, one row per gene per contrast.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""Cross-tissue meta-analysis + overlap-count of the per-tissue DESeq2 results.

Downstream stage of the de_v2 5S-variant DE pipeline: collapses the per-tissue
DESeq2 tables into one row per gene per (scope, contrast) by combining each gene's
log2FC across tissues under an inverse-variance fixed-effect model.

Usage:
    python de_meta.py            # no CLI args; scans all of de_v2/out/de_*.tsv

Inputs  (upstream, from de_pertissue.py / de_continuous.py / de_dose.py / de_joint*.py):
    de_v2/out/de_*.tsv           one per (scope, contrast, tissue); the (scope, contrast)
                                 grouping key is read from the file *content* (scope/contrast
                                 columns), not the filename. Columns consumed here:
                                 ensg, log2FoldChange, lfcSE, padj, scope, contrast.

Per gene, within each (scope, contrast) group, computes:
  - inverse-variance (fixed-effect) meta of log2FC across tissues -> beta_meta, se_meta, z, pval, padj
  - Cochran Q / I^2 heterogeneity across tissues
  - n_sig_same_dir: overlap-count = #tissues with padj<0.1 AND same sign as meta beta
                    (threshold-based "robust across tissues" readout; biased toward big effects)
  - frac_consistent / consistent_dir / sign_p / sign_padj: a non-thresholded directional
    (sign) consistency test, magnitude-independent, crediting small coordinated effects.

Outputs:
    de_v2/meta/meta_<scope>_<contrast>.tsv   one row per gene, columns:
        ensg, sum_w, sum_wl, n_tissue, beta_meta, se_meta, z, pval, Q, df, I2,
        n_sig_same_dir, frac_consistent, consistent_dir, sign_p, padj, sign_padj
    de_v2/meta/SUMMARY.tsv                   one row per (scope, contrast): tissue/gene counts,
        n_meta_sig (padj<0.1), n_robust_3tiss (n_sig_same_dir>=3).

CAVEAT: the global genomic inflation (lambda ~= 3-4) in the meta z/p is a structural consequence
of a fixed-effect meta over tissues that share donors (a permutation reproduces it), so padj here
is anti-conservative; use a permutation-calibrated test (DE_PERMUTE) for inference.
"""
import os, glob, numpy as np, pandas as pd
from scipy import stats
import de_common as C

OUT = f"{C.ROOT}/de_v2/out"; MOUT = f"{C.ROOT}/de_v2/meta"; os.makedirs(MOUT, exist_ok=True)


def meta_group(files):
    """Inverse-variance fixed-effect meta over one (scope, contrast) group of per-tissue files.

    Stacks all tissues' per-gene log2FC + lfcSE into one long frame, then aggregates by gene.
    Returns one row per gene, sorted by meta p-value.
    """
    # Stack every tissue's table; each row is one (gene, tissue) DESeq2 estimate.
    frames = []
    for f in files:
        d = pd.read_csv(f, sep="\t", usecols=["ensg", "log2FoldChange", "lfcSE", "padj"])
        frames.append(d.rename(columns={"log2FoldChange": "lfc", "lfcSE": "se"}))
    big = pd.concat(frames, ignore_index=True)
    # Drop rows DESeq2 could not estimate (NA lfc/se, or se<=0 -> weight undefined).
    big = big[np.isfinite(big.lfc) & np.isfinite(big.se) & (big.se > 0)].copy()
    # Inverse-variance weights: each tissue's estimate is weighted by 1/se^2 (precision),
    # so tissues that estimated the effect more precisely dominate the pooled beta.
    big["w"] = 1.0 / big.se.values ** 2
    big["wl"] = big.w * big.lfc                              # weight * effect, for the weighted sum
    g = big.groupby("ensg")
    m = pd.DataFrame({"sum_w": g.w.sum(), "sum_wl": g.wl.sum(), "n_tissue": g.lfc.size()})
    m["beta_meta"] = m.sum_wl / m.sum_w                      # pooled effect = Sum(w*lfc) / Sum(w)
    m["se_meta"] = np.sqrt(1.0 / m.sum_w)                    # pooled SE = 1/sqrt(Sum(w))
    m["z"] = m.beta_meta / m.se_meta                         # Wald z of the pooled effect
    m["pval"] = 2 * stats.norm.sf(np.abs(m.z))               # two-sided normal tail p (see lambda caveat above)
    # Heterogeneity: Cochran Q = Sum(w * (lfc - beta_meta)^2), the weighted dispersion of the
    # per-tissue estimates about the pooled effect. Join beta_meta back onto each row to vectorize.
    big = big.join(m["beta_meta"], on="ensg")
    big["q"] = big.w * (big.lfc - big.beta_meta) ** 2
    m["Q"] = big.groupby("ensg").q.sum()
    m["df"] = m.n_tissue - 1                                 # Q ~ chi2_{df} under homogeneity
    # I^2 = fraction of total variation due to between-tissue heterogeneity (0 if Q<=df), as a percent.
    m["I2"] = np.where((m.Q > 0) & (m.df > 0), np.maximum(0.0, (m.Q - m.df) / m.Q) * 100, 0.0)
    # overlap-count: tissues padj<0.1 with sign == meta sign (threshold-based; biased to big effects)
    big["sig_same"] = (big.padj < 0.1) & (np.sign(big.lfc) == np.sign(big.beta_meta))
    m["n_sig_same_dir"] = big.groupby("ensg").sig_same.sum().astype(int)
    # NON-THRESHOLDED directional consistency (magnitude-independent; credits small coordinated effects):
    # per gene, count +/- tissues, majority fraction, and a two-sided binomial sign-test p (H0: 50/50).
    npos = big.assign(pos=(big.lfc > 0).astype(int)).groupby("ensg").pos.sum()
    ntot = m.n_tissue.astype(int)
    nmaj = np.maximum(npos, ntot - npos)
    m["frac_consistent"] = nmaj / ntot                                  # 0.5 (random) .. 1.0 (all one way)
    m["consistent_dir"] = np.where(npos >= ntot - npos, 1, -1)
    # Two-sided binomial sign-test p: P(majority count >= nmaj | H0 = 50/50 up/down). sf(k-1)=P(>=k).
    m["sign_p"] = np.clip(2.0 * stats.binom.sf(nmaj.values - 1, ntot.values, 0.5), 0.0, 1.0)
    # Benjamini-Hochberg FDR, applied independently to the meta z-test and to the sign test.
    def bh(pv):
        pv = pv.astype(float); out = np.full(len(pv), np.nan)
        ok = np.isfinite(pv) & (pv >= 0) & (pv <= 1)                 # only adjust valid p-values
        if ok.any(): out[ok] = stats.false_discovery_control(np.clip(pv[ok], 0.0, 1.0))
        return out
    m["padj"] = bh(m.pval.values)                           # NB: anti-conservative (lambda~3-4 inflation)
    m["sign_padj"] = bh(m.sign_p.values)
    return m.reset_index().sort_values("pval")


def main():
    files = glob.glob(f"{OUT}/de_*.tsv")
    if not files:
        print("no per-tissue results yet"); return
    # Bucket the per-tissue files by (scope, contrast). The key is read from the file's first data
    # row (scope/contrast columns), which is robust to scope/contrast/tissue names that contain '_'
    # and would be ambiguous to parse from the filename.
    keyed = {}
    for f in files:
        h = pd.read_csv(f, sep="\t", usecols=["scope", "contrast"], nrows=1)
        keyed.setdefault((h.scope[0], h.contrast[0]), []).append(f)
    summary = []
    for (scope, contrast), fs in sorted(keyed.items()):
        m = meta_group(fs)                                  # meta over the tissues in this group
        tag = f"{scope}_{contrast}"
        m.to_csv(f"{MOUT}/meta_{tag}.tsv", sep="\t", index=False)
        n_meta = int((m.padj < 0.1).sum())                  # genes passing the meta FDR
        n_robust = int((m.n_sig_same_dir >= 3).sum())       # genes significant+concordant in >=3 tissues
        summary.append((scope, contrast, len(fs), len(m), n_meta, n_robust))
        print(f"{tag:40s} tissues={len(fs):3d} genes={len(m):6d} meta_padj<0.1={n_meta:5d} robust(>=3 tiss)={n_robust:4d}")
    pd.DataFrame(summary, columns=["scope", "contrast", "n_tissue", "n_gene",
                                   "n_meta_sig", "n_robust_3tiss"]).to_csv(f"{MOUT}/SUMMARY.tsv", sep="\t", index=False)
    print(f"\nwrote {MOUT}/SUMMARY.tsv")


if __name__ == "__main__":
    main()
