# Reference inputs (public)

Small public reference assets the scripts consume. To be copied during curation:

- `5S_population_consensus_unit.fa` — the 2168 bp population-consensus repeat unit (calling/mapping target).
- `repeat_unit_features.tsv` — feature map (5S gene 630–748, antisense AluY 787–1066, ICR Box A/IE/Box C,
  Pol III terminator, microsatellites).
- `consensus_reference.json` — consensus + masking metadata from `12_consensus_rederive.py`.
- Gene sets for GSEA (MSigDB Hallmark + Reactome `.gmt`, custom cyto-/mito-RP sets) — from `GTEx/de_v2/`.

Large public references (GRCh38, CHM13v2.0) are obtained from their sources (see top-level README),
not vendored here.
