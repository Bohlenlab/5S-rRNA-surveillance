# Figure 2 + S2 — short-read 5S caller and UK Biobank

**UK Biobank individual data is controlled (Application 98772)** and not included here.

## Short-read caller (public / HPRC benchmark)
| Panel | Script | Where |
|---|---|---|
| 2A,2B SR vs assembly VAF; HiFi rescue | `75_platform_comparison.py` | [local] |
| 2C calls-by-group vs VAF | `74_summary_figure.py` | [local] |
| S2A,B pseudogene specificity (wgsim/bwa) | `figS_pseudogene.py` + `WGS-Variant-Identification/` sim reads | [local] |
| S2C HPRC SR coverage | `13c_hprc_sr_coverage_S2C.py` | [local] |
| calibration backbone | `61–73`, `fig2_rescue_*.py` | [local] |

## UK Biobank
| Stage / panel | Script | Where |
|---|---|---|
| Read extraction from full CRAMs | runs on the UK Biobank Research Analysis Platform; described in Methods, not included here | [RAP] |
| extraction→align→variant→trio-QC | `5S_setup/{0,1,2,4,6}_*` | [server] |
| 2D carrier landscape | `78a_panel_A_export.py` | [local] |
| 2F VAF concordance vs assembly | `93_ukbb_assembly_consistency.py` | [local] |
| 2G twin/parent-child sharing; S2D Falconer | `94_heritability_vaf.py` | [local] |
| 2H gene-conversion by VAF | `93b_ukbb_assembly_consistency_all.py` | [local] |
| 2I variants per carrier | `76_import_ukbb_population.py` (carriers table) | [server] |
| PheWAS: ICD-10 dosage logistic regressions | `81_ukbb_phenotype_associations.py` (vs `cohort_icd10.db`) | [server] |
| 2J ICD-10 positional enrichment | `95_icd10_positional_enrichment.py` | [server] |
| 2K,2L dosage-OR curves (OCD/SSc) | `85_dosage_or_curves.py` | [server] |
| variant landscape / selection setup | `78_/79_ukbb_variant_landscape*.py`, `80_ukbb_purifying_selection.py` | [server] |

The upstream UK Biobank read-extraction step runs on the UK Biobank Research Analysis Platform and
is described in the Methods; it is not included in this repository.
