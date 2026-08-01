# Database — schema and build/query code

The analyses read a central SQLite store, `5S_rDNA.db`. **The populated database is not
shipped** — it contains derived controlled-access genotypes (UK Biobank, GTEx, TCGA/CPTAC). Its
shareable, public-assembly-derived content is published as the paper's Supplementary Tables S1–S11.

Build and query scripts present here:
- `build_database.py`, `import_legacy_assemblies.py`, `import_r2_assemblies.py`,
  `import_read_variants.py` — construction and import.
- `12_consensus_rederive.py`, `59b_rederive_consensus_migrated.py`, `59c_rederive_gene_unit_t2t.py`
  — variant polarization vs the population consensus.
- `57_extract_gfa_haplotypes.py`, `58_multicontig_extract.py`, `59_migrate_multicontig.py`,
  `60_multicontig_methylation.py`, `60b_load_multicontig_methylation.py` — multi-contig array
  correction (recovers copies on non-dominant contigs; orientation fix; methylation load).
- `SCHEMA.md` — schema reference.

Key tables: `assembly`, `haplotype`, `copy` (filter `array_member=1`), `variant`
(`alignment_source='consensus_t2t'` primary; SNP analyses `var_type='snp'`), `copy_methylation`,
`read_variant`, `hap_site_freq`. Coordinates are 1-based in the 2168 bp unit (gene 630–748).
