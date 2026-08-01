# Database — schema and build/query code

The analyses read a central SQLite store, `5S_rDNA.db` (~960 MB). **The populated database is not
shipped** — it contains derived controlled-access genotypes (UK Biobank, GTEx, TCGA/CPTAC). Its
shareable, public-assembly-derived content is published as the paper's Supplementary Tables S1–S11.

What belongs here (to be copied during curation):
- `build_database.py`, `import_legacy_assemblies.py`, `import_r2_assemblies.py`, `import_read_variants.py` — construction.
- `12_consensus_rederive.py`, `09/10_*rederive_t2t*` — variant polarization vs population consensus.
- `57–60b` — multi-contig array correction (recovers copies on non-dominant contigs; orientation fix).
- `T2T/5S_rDNA_database_documentation.md` — schema reference.

Key tables: `assembly`, `haplotype`, `copy` (filter `array_member=1`), `variant`
(`alignment_source='consensus_t2t'` primary; SNP analyses `var_type='snp'`), `copy_methylation`,
`read_variant`, `hap_site_freq`. Coordinates are 1-based in the 2168 bp unit (gene 630–748).
