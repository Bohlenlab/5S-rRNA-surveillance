# Reusable methods

The transferable computational ideas from the study, factored out and documented so they can be
reused on other loci. Each maps back to the figure code (see `../CODE_MANIFEST.md`).

`collapsed_array_caller_demo.py` is a self-contained, dependency-free demonstration of the
collapsed-array low-VAF caller (building block 1). It runs in under a second and needs no external
data or packages.

1. **Collapsed-array short-read variant caller** — map short reads to a single consensus repeat
   unit so all array copies pool into one pileup; call at low VAF; rescue apparent false positives
   with orthogonal long reads; set the threshold by F1. (from `5S_setup` + `74/75`, `61–73`)
2. **Multi-contig per-copy array reconstruction & 5′→3′ ordering** — recover copies split across
   assembly contigs; strand-aware extraction; flank-anchored ordering. (`57–60b`)
3. **Runs-based within-array clustering index** — Wald–Wolfowitz exact null for k carriers among n
   ordered copies; per-haplotype index of gene-conversion tract structure. (`31`, `47`)
4. **Content-adjusted regional site-frequency spectrum** — per-context (CpG-aware) expected rates
   to compare rare/common variant density across sub-regions. (`27`)
5. **Edge-distance methylation/accessibility profiling on collapsed arrays** — competitive/anchored
   alignment + signed distance-into-array binning; generalizes across 5S, 45S, RNU2. (`40`,
   `45S_methylation`, `fiberseq_5S`)
6. **Fraction-seq → per-variant expression & incorporation scoring** — count per-position variant
   frequencies across gradient fractions; subtract a GFP-control error model. (`process_rep_mutation_ratios.py`)
7. **Variant-dose → per-tissue DESeq2 + cross-tissue meta + GSEA** — continuous dosage predictor,
   covariate-adjusted NB-GLM per tissue, inverse-variance meta, prerank GSEA. (`de_v2`)
8. **p53-stratified genomic-dose → RP-repression interaction** — expanded p53 classification and
   dose × p53-status interaction on the cytosolic-RP module. (`cancer_5S`)
