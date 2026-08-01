# Figure 6 + S6 — cancer p53 axis and cross-cohort selection

Cancer cohorts controlled (dbGaP phs000178 TCGA, phs001287 CPTAC; project #44042).
Cancer DE ran on a compute server; final panels assembled from `cancer_5S/surveillance_v2/`.

| Panel | Script | Where |
|---|---|---|
| 6A–C function-stratified trans-effects (IE) | `GTEx/de_v2/make_figs3to6.py`, `make_fig13_IE_*.py` | [local]/[server] |
| 6D–L p53-stratified DE + GSEA | `cancer_5S/surveillance_v2/` + `scripts/{03_slice_call,54–77}` | [server] run → [local] figs |
| 6Q somatic 5S gains × p53 | `cancer_5S/scripts/{100–111}` | [server] |
| 6M–P cross-cohort selection (HPRC/UKBB/GTEx) | `97_incorporation_depletion_three_cohort.py`, `98_binary_classifications_two_ways.py`, `2E_fair_perdonor_vaf.py`, `de_v2/make_figs3to6.py` | [local] |
| S6A methylation by variant consequence | `24_gene_methylation_by_functional_consequence.py` (+`24b–g`) | [local] |
| robustness (substitution-covariate) | `97c_substitution_covariate_robustness.py` | [local] |
| MDM2–p53 structure asset | `MDM2/1RV1.cif` | [local] |

Cancer 5S genotyping streams region slices from the GDC slicing API (`03_slice_call.py`) → same
bwa/bcftools consensus-unit recipe as GTEx. p53 status = TP53 mutation + MDM2/MDM4 amp + CDKN2A/ARF
del (+ HPV⁺ for CPTAC). `surveillance_v2/00_CHECKLIST.md` records which DE model fed the figure.
