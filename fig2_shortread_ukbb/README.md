# Figure 2 + S2 — short-read 5S caller and UK Biobank

Split across machines. **UK Biobank individual data is controlled (Application 98772)** and not here.

## Short-read caller (public / HPRC benchmark)
| Panel | Script | Where |
|---|---|---|
| 2A,2B SR vs assembly VAF; HiFi rescue | `75_platform_comparison.py` | Mac |
| 2C calls-by-group vs VAF (sens 83/prec 69) | `74_summary_figure.py` | Mac |
| S2A,B pseudogene specificity (wgsim/bwa) | `figS_pseudogene.py` + `WGS-Variant-Identification/` sim reads | Mac |
| S2C HPRC SR coverage | `13c_hprc_sr_coverage_S2C.py` | Mac |
| calibration backbone | `61–73`, `fig2_rescue_*.py` | Mac |

## UK Biobank
| Stage / panel | Script | Where |
|---|---|---|
| **Read extraction from full CRAMs** | **DNAnexus RAP applet — TO RETRIEVE or describe in Methods** | RAP |
| extraction→align→variant→trio-QC | `5S_setup/{0,1,2,4,6}_*` | trr237 (authoritative) |
| 2D carrier landscape | `78a_panel_A_export.py` | Mac |
| 2F VAF concordance vs assembly | `93_ukbb_assembly_consistency.py` | Mac |
| 2G twin/parent-child sharing; S2D Falconer | `94_heritability_vaf.py` | Mac |
| 2H gene-conversion by VAF | `93b_ukbb_assembly_consistency_all.py` | Mac |
| 2I variants per carrier | `76_import_ukbb_population.py` (carriers table) | immuno2 |
| **PheWAS: ICD-10 dosage logistic regressions** | `81_ukbb_phenotype_associations.py` (vs `cohort_icd10.db`) | immuno2 (authoritative) |
| 2J ICD-10 positional enrichment | `95_icd10_positional_enrichment.py` | immuno2 |
| 2K,2L dosage-OR curves (OCD/SSc) | `85_dosage_or_curves.py` | immuno2 |
| variant landscape / selection setup | `78_/79_ukbb_variant_landscape*.py`, `80_ukbb_purifying_selection.py` | immuno2 |

The one missing code step is the RAP-side read extraction (decision D1).
