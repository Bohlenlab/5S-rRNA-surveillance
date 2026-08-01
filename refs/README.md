# Reference inputs (public)

Small public reference assets the scripts consume:

- `5S_consensus_unit.fa` — the 2168 bp population-consensus repeat unit (calling/mapping target).
- `repeat_unit_features.tsv` — feature map (5S gene 630–748, antisense AluY 787–1066, ICR Box A/IE/Box C,
  Pol III terminator, microsatellites).
- `consensus_reference.json` — consensus + masking metadata from `12_consensus_rederive.py`.

Gene sets for GSEA (MSigDB Hallmark + Reactome `.gmt`) are not redistributed here; they are obtained
from MSigDB.

Large public references (GRCh38, CHM13v2.0) are obtained from their sources (see top-level README),
not vendored here.
