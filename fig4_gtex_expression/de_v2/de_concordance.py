# -----------------------------------------------------------------------------
# de_concordance.py — DNA-signature x RNA-dose concordance test: Spearman /
# sign-concordance / permutation p between the carrier and RNA-dose meta effects.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""DNA-signature x RNA-dose concordance test for the de_v2 5S-variant DE pipeline.

Logic: if 5S variant V has an effect, the genes perturbed by merely CARRYING V (a DNA/genotype
contrast) should be perturbed in the SAME DIRECTION when carriers EXPRESS MORE of V (an RNA-dose
contrast). The DNA contrast is used as the anchor; the within-carrier RNA-dose contrast is the
independent comparison.

For each variant V (and the aggregate):
    beta_DNA = carrier-vs-noncarrier meta effect per gene   (anchor "signature")
    beta_RNA = within-carrier RNA-VAF dose meta effect per gene
We measure agreement by Spearman(beta_DNA, beta_RNA) both genome-wide and restricted to the top
DNA-signature genes (largest |z_DNA|), a sign-concordance fraction on that top set, and a
permutation p-value (shuffle RNA gene labels). Positive and significant indicates directional agreement.

Usage:
    python de_concordance.py     # no CLI args; scans de_v2/meta/meta_dose_*.tsv (run de_meta.py first)

Inputs  (all produced by de_meta.py, one row per gene):
    de_v2/meta/meta_dose_<scope>_dose.tsv    RNA-dose meta (beta_RNA); the scan target
    de_v2/meta/meta_dna_carriervsNC.tsv      DNA anchor for the 'dose_AGG' aggregate scope
    de_v2/meta/meta_var_<V>_carriervsref.tsv DNA anchor for a per-variant scope 'dose_<V>'
    (each read for columns beta_meta, z, n_tissue)

Output:
    de_v2/meta/CONCORDANCE.tsv   one row per variant, sorted by perm_p_top; columns:
        variant, n_genes, n_top, rho_all, rho_top, sign_concord_top, perm_p_top
"""
import os, re, glob, numpy as np, pandas as pd
from scipy import stats
import de_common as C
MDIR = f"{C.ROOT}/de_v2/meta"


def load_beta(path):
    """Load a meta table indexed by gene, keeping the effect (beta_meta), z, and tissue count."""
    d = pd.read_csv(path, sep="\t").set_index("ensg")
    return d[["beta_meta", "z", "n_tissue"]]


def dna_path_for(scope):
    """Map an RNA-dose scope to its matching DNA-anchor meta file, or None if that file is absent.

    'dose_AGG' (the aggregate across variants) pairs with the aggregate carrier-vs-NC anchor;
    a per-variant scope 'dose_<V>' pairs with that variant's carrier-vs-ref anchor.
    """
    if scope == "dose_AGG":
        p = f"{MDIR}/meta_dna_carriervsNC.tsv"
    else:
        v = scope.replace("dose_", "")                      # recover the variant id V from 'dose_<V>'
        p = f"{MDIR}/meta_var_{v}_carriervsref.tsv"
    return p if os.path.exists(p) else None


def concordance(bd, br, top=300, seed=0):
    """Correlate the DNA-anchor (bd) and RNA-dose (br) per-gene effects; return None if too sparse.

    bd/br are load_beta frames indexed by gene. `top` = size of the DNA-signature gene subset;
    `seed` seeds the permutation null. Returns a dict of summary stats (or None if <100 shared genes).
    """
    # Inner-join the two effect vectors on gene; keep only genes with a finite effect in both.
    j = bd.join(br, lsuffix="_dna", rsuffix="_rna", how="inner").dropna(subset=["beta_meta_dna", "beta_meta_rna"])
    if len(j) < 100:                                        # not enough shared genes to say anything
        return None
    x = j.beta_meta_dna.values; y = j.beta_meta_rna.values
    rho_all = stats.spearmanr(x, y).statistic              # genome-wide rank correlation of the two effects
    # Restrict to the top DNA-signature genes (largest |z_dna|): where the DNA effect is strongest,
    # the RNA-dose effect should track it most cleanly.
    t = j.reindex(j.z_dna.abs().sort_values(ascending=False).index).head(top)
    rho_top = stats.spearmanr(t.beta_meta_dna, t.beta_meta_rna).statistic
    # Sign-concordance: fraction of top genes where DNA and RNA effects point the same way (0.5 = chance).
    sign_top = float(np.mean(np.sign(t.beta_meta_dna) == np.sign(t.beta_meta_rna)))
    # Permutation null: break the DNA<->RNA gene pairing by shuffling RNA labels, recompute rho 2000x.
    # p = fraction of |null rho| at least as extreme as observed |rho_top| (+1 smoothing on num & denom).
    rng = np.random.default_rng(seed)
    null = np.array([stats.spearmanr(t.beta_meta_dna.values, rng.permutation(t.beta_meta_rna.values)).statistic
                     for _ in range(2000)])
    p_top = (np.sum(np.abs(null) >= abs(rho_top)) + 1) / (len(null) + 1)
    return dict(n_genes=len(j), n_top=len(t), rho_all=rho_all, rho_top=rho_top,
                sign_concord_top=sign_top, perm_p_top=p_top)


def main():
    rows = []
    # Each RNA-dose meta file is meta_<scope>_dose.tsv with scope starting 'dose_'; iterate them.
    for f in sorted(glob.glob(f"{MDIR}/meta_dose_*.tsv")):
        m = re.match(r'^meta_(dose_.+)_dose\.tsv$', os.path.basename(f))  # capture scope = dose_<V> or dose_AGG
        if not m:
            continue
        scope = m.group(1)
        dnap = dna_path_for(scope)                          # locate the matching DNA anchor
        if not dnap:
            print(f"{scope}: no DNA anchor -> skip"); continue
        res = concordance(load_beta(dnap), load_beta(f))    # DNA anchor vs RNA-dose for this variant
        if res is None:
            print(f"{scope}: too few shared genes"); continue
        res["variant"] = scope.replace("dose_", "")
        rows.append(res)
        print(f"{res['variant']:8s} genes={res['n_genes']:5d} rho_all={res['rho_all']:+.3f} "
              f"rho_top{res['n_top']}={res['rho_top']:+.3f} sign={res['sign_concord_top']:.2f} "
              f"perm_p={res['perm_p_top']:.4f}")
    if rows:
        out = pd.DataFrame(rows)[["variant", "n_genes", "n_top", "rho_all", "rho_top", "sign_concord_top", "perm_p_top"]]
        out = out.sort_values("perm_p_top")
        out.to_csv(f"{MDIR}/CONCORDANCE.tsv", sep="\t", index=False)
        print(f"\nwrote {MDIR}/CONCORDANCE.tsv")


if __name__ == "__main__":
    main()
