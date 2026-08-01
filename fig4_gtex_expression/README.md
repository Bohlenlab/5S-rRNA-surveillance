# Figure 4 + S4 — GTEx variant 5S rRNA expression and trans-effects

GTEx is controlled (dbGaP phs000424). Canonical DE pipeline = **`GTEx/de_v2/`** (PI-confirmed);
earlier iterations (`figures/de_expression`, `scripts/eq*`, `ex01–04`, `vm_archive`) are superseded.
The DE run itself executed on trr237 (`~/de_v2/`, SLURM logs + `meta/SUMMARY.tsv`).

| Panel | Script |
|---|---|
| 4B RNA-seq depth profile | `fig4_coverage_pooldata.py` + `fig4_coverage_panel.py` |
| 4C per-variant carrier vs non-carrier | `fig4_panel4_pooldata.py` + `fig4_build_panels.py` |
| 4D rank-skew volcano | `fig4_build_panels.py` (from `wgs_rna_rank_skew.py`) |
| 4E expressed-variant distribution | `fig4_p3_detection.py` |
| 4F prevalence (AUC q99) | `fig4_panel4_build.py` — **DECISION D2: q99 37% vs strict 14%** |
| 4G DNA vs RNA VAF | `fig4_build_panels.py` |
| 4H trans-effect volcano | `ex05j_p10_4x4_biogenesis.py` (DE from `ex05e_fullmodel_gsea.py`) |
| 4I–L dose-response / GSEA / meta | `ex05g_dnarna_vaf_panel.py`, `de_v2/{de_common,de_meta,de_gsea}.py` |
| genotyping (DNA/RNA VAF) | `GTEx/scripts/{10–14,20,23}` |

DE model: per-tissue negative-binomial (pydeseq2) ~ RIN + ischemic time + genotype PCs + sex +
Hardy + batch + dosage; cross-tissue inverse-variance meta; GSEA (gseapy prerank). See building block 7.
