# Figure 1 + S1 — 5S array structure and gene-conversion-shaped variation

All local (`HPRC/scripts`). Reads the assembly variant catalog in `5S_rDNA.db`
(filter `array_member=1`; SNP analyses `var_type='snp'`). Public inputs (HPRC/CPC assemblies).

| Panel | Script |
|---|---|
| 1B per-copy variant map | `26_array_maps.py` |
| 1C / S1A copy number | `07_array_structure_analysis.py` |
| 1D, 1I variants per hap / gene | `07_...`, `51_variants_per_haplotype_table.py` |
| 1E copies-vs-donors per variant | `07_...`, `52_variant_catalog_table.py` |
| 1F substitution spectrum (normalized) | `39_substitution_normalized_cpg.py` |
| 1G / S1B clustering index | `31_variant_clustering_runs.py`, `47_clustering_single_panels.py` |
| 1H tract-size spectrum (gBGC) | `48_gbgc_resolved.py`, `49_gbgc_trend_stats.py` |
| 1J, 1K tracts by region | `50_geneconv_gene_alu_nts.py` |
| 1K content-adjusted regional SFS | `27_population_genetics.py` |
| Table S2 repeat-unit annotation | `44_repeat_unit_annotation.py` |

DB construction (shared): `build_database.py`, `12_consensus_rederive.py`, multi-contig `57–60b`.
See `../db/`. Core methods generalized in `../building_blocks/`.
