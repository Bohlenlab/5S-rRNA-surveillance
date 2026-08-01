# Figure 4 + S4 — GTEx variant 5S rRNA expression and trans-effects

GTEx is controlled (dbGaP phs000424). The differential-expression pipeline is `GTEx/de_v2/`,
run on a compute server (`meta/SUMMARY.tsv`).

| Panel | Script |
|---|---|
| 4B RNA-seq depth profile | `fig4_coverage_pooldata.py` + `fig4_coverage_panel.py` |
| 4C per-variant carrier vs non-carrier | `fig4_panel4_pooldata.py` + `fig4_build_panels.py` |
| 4D rank-skew volcano | `fig4_build_panels.py` (from `wgs_rna_rank_skew.py`) |
| 4E expressed-variant distribution | `fig4_p3_detection.py` |
| 4F prevalence of expressed variants | `genotyping_expression/prevalence_strict.py` |
| 4G DNA vs RNA VAF | `fig4_build_panels.py` |
| 4H trans-effect volcano | `ex05j_p10_4x4_biogenesis.py` (DE from `ex05e_fullmodel_gsea.py`) |
| 4I–L dose-response / GSEA / meta | `ex05g_dnarna_vaf_panel.py`, `de_v2/{de_common,de_meta,de_gsea}.py` |
| genotyping (DNA/RNA VAF) | `GTEx/scripts/{10–14,20,23}` |

Panel 4F reports prevalence using a strict, background-exceeding call (`prevalence_strict.py`),
which gives ~14%.

DE model: per-tissue negative-binomial (pydeseq2) ~ RIN + ischemic time + genotype PCs + sex +
Hardy + batch + dosage; cross-tissue inverse-variance meta; GSEA (gseapy prerank). See building block 7.
