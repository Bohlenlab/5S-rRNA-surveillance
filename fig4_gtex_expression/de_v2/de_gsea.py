# -----------------------------------------------------------------------------
# de_gsea.py — pre-ranked GSEA (gseapy, local Hallmark + Reactome gmts) on the
# meta or per-tissue DESeq2 statistics.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""GSEA (gseapy prerank, local Hallmark + Reactome gmts) on meta or per-tissue DESeq2 stats.

Downstream stage of the de_v2 5S-variant DE pipeline: turns a per-gene DE statistic into a
ranked gene list and runs pre-ranked GSEA against two MSigDB-derived gene-set libraries
(Hallmark 2020, Reactome 2022) shipped as local .gmt files.

Genes are ranked by their signed DE statistic (meta z, or the per-tissue Wald `stat`), highest
first (most up-regulated). ENSG ids are mapped to HGNC symbols from the gencode v49 GTF -- this
matters because the .gmt gene sets are keyed by symbol, not ENSG.

Usage:
    de_gsea.py meta                 # GSEA on every de_v2/meta/meta_*.tsv (ranked by meta z)
    de_gsea.py tissue <glob>        # GSEA on per-tissue de_v2/out/<glob> files (ranked by Wald stat)

Inputs:
    de_v2/meta/meta_*.tsv           ('meta' mode)   from de_meta.py; needs columns ensg, z
    de_v2/out/de_*.tsv              ('tissue' mode)  from de_pertissue.py etc.; needs ensg, stat
    gencode v49 GTF (DE_GTF)        for the ENSG->symbol map (cached to ensg2symbol.tsv)
    Hallmark / Reactome .gmt        (DE_GMT_HALLMARK / DE_GMT_REACTOME)

Outputs:
    de_v2/gsea/gsea_<tag>_<lib>.tsv one per input file x library (lib in {hallmark, reactome});
                                    <tag> is the input basename minus the meta_/de_ prefix and .tsv.
                                    Columns are gseapy prerank's res2d (Term, ES, NES, NOM p-val,
                                    FDR q-val, ...). Also prints the top-3 terms by FDR per run.
"""
import os, sys, glob, gzip, re, numpy as np, pandas as pd, gseapy as gp
from pathlib import Path
import de_common as C

GTF = os.environ.get("DE_GTF", str(Path(os.environ.get("FIVES_REFS", "refs")) / "gencode.v49.annotation.gtf.gz"))
GMT = {"hallmark": os.environ.get("DE_GMT_HALLMARK", str(Path(os.environ.get("FIVES_REFS", "refs")) / "Enrichr.MSigDB_Hallmark_2020.gmt")),
       "reactome": os.environ.get("DE_GMT_REACTOME", str(Path(os.environ.get("FIVES_REFS", "refs")) / "Enrichr.Reactome_2022.gmt"))}
GOUT = f"{C.ROOT}/de_v2/gsea"; os.makedirs(GOUT, exist_ok=True)
SYMCACHE = f"{C.ROOT}/de_v2/ensg2symbol.tsv"


def ensg2symbol():
    """Return {ENSG -> gene symbol}, built once from the gencode GTF and cached to SYMCACHE.

    ENSG ids are captured WITHOUT their version suffix ([^."]+ stops at the '.'), so they match
    the un-versioned ensg column in the DE tables.
    """
    if os.path.exists(SYMCACHE):                            # reuse the cached map if already built
        s = pd.read_csv(SYMCACHE, sep="\t")
        return dict(zip(s.ensg, s.symbol))
    m = {}
    with gzip.open(GTF, "rt") as fh:
        for ln in fh:
            if ln.startswith("#") or "\tgene\t" not in ln:  # only 'gene' feature lines carry the id<->name pair
                continue
            gid = re.search(r'gene_id "([^."]+)', ln)        # ENSG, up to (excluding) the version dot
            gn = re.search(r'gene_name "([^"]+)"', ln)       # HGNC symbol (what the .gmt gene sets use)
            if gid and gn:
                m[gid.group(1)] = gn.group(1)
    pd.DataFrame({"ensg": list(m), "symbol": list(m.values())}).to_csv(SYMCACHE, sep="\t", index=False)
    return m


def run_prerank(rnk, tag):
    """Run gseapy pre-ranked GSEA for each gene-set library and write/report the results.

    `rnk` is a 2-column (symbol, score) frame already sorted descending; prerank walks this
    ranking and scores each gene set's enrichment at the top vs bottom, with a 1000-permutation
    null for the NES/FDR. min_size/max_size gate out gene sets too small or large to test.
    """
    for lib, gmt in GMT.items():
        try:
            pre = gp.prerank(rnk=rnk, gene_sets=gmt, min_size=10, max_size=500,
                             permutation_num=1000, seed=0, threads=4, no_plot=True, verbose=False)
            res = pre.res2d.copy()
            res.to_csv(f"{GOUT}/gsea_{tag}_{lib}.tsv", sep="\t", index=False)   # full enrichment table
            # Print the 3 most significant terms (smallest FDR) as a console summary.
            res["FDR q-val"] = pd.to_numeric(res["FDR q-val"], errors="coerce")
            top = res.sort_values("FDR q-val").head(3)
            for _, t in top.iterrows():
                print(f"    [{lib}] {t['Term'][:48]:48s} NES={float(t['NES']):+.2f} FDR={float(t['FDR q-val']):.3f}")
        except Exception as e:
            # A library can fail (e.g. no gene set meets min/max size after symbol mapping); keep going.
            print(f"    [{lib}] FAIL {type(e).__name__}: {str(e)[:90]}")


def main():
    mode = sys.argv[1]                                      # "meta" or "tissue"
    sym = ensg2symbol()
    if mode == "meta":
        # Cross-tissue meta tables: rank genes by the fixed-effect meta z-score.
        files = sorted(glob.glob(f"{C.ROOT}/de_v2/meta/meta_*.tsv"))
        for f in files:
            tag = re.sub(r'^meta_|\.tsv$', '', os.path.basename(f))         # drop meta_ prefix + .tsv
            d = pd.read_csv(f, sep="\t")
            d["symbol"] = d.ensg.map(sym)
            # Need a symbol and a rank score; collapse duplicate symbols (keep first) so the ranking is 1:1.
            d = d.dropna(subset=["symbol", "z"]).drop_duplicates("symbol")
            rnk = d[["symbol", "z"]].sort_values("z", ascending=False).reset_index(drop=True)
            print(f"{tag} (n={len(rnk)} genes):")
            run_prerank(rnk, tag)
    else:
        # Per-tissue mode: rank by the DESeq2 Wald `stat` of each tissue file matched by sys.argv[2].
        for f in sorted(glob.glob(f"{C.ROOT}/de_v2/out/{sys.argv[2]}")):
            tag = re.sub(r'^de_|\.tsv$', '', os.path.basename(f))           # drop de_ prefix + .tsv
            d = pd.read_csv(f, sep="\t")
            d["symbol"] = d.ensg.map(sym)
            d = d.dropna(subset=["symbol", "stat"]).drop_duplicates("symbol")
            rnk = d[["symbol", "stat"]].sort_values("stat", ascending=False).reset_index(drop=True)
            print(f"{tag} (n={len(rnk)} genes):")
            run_prerank(rnk, tag)


if __name__ == "__main__":
    main()
