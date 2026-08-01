# Code for: Surveillance and selection of 5S ribosomal RNA genes in the human genome

Sengl, Bagaric, Conil, Seeleuthner, Müller, Cobat, Bohlen.

This repository contains the analysis code for the study. It is organized by manuscript figure,
plus a `building_blocks/` folder highlighting the reusable methods. **`CODE_MANIFEST.md` is the
index**: it maps every figure panel to the exact script that produced it.

## What is and isn't here

- **Code, yes.** The scripts that implement the analyses and the novel methods.
- **Data, no.** All human-subjects data is controlled-access and is *not* redistributed here. The
  repository references datasets by accession; obtain them through the paths below.
- **Infrastructure, no.** Compute ran on a Linux compute server; server paths and job
  scripts are not part of the release. Captured software environments are in `env/`.
- **Not every panel is reproducible end-to-end from raw data.** Some final panels were assembled
  in GraphPad Prism / Adobe Illustrator from the tables the scripts emit; those steps are noted in
  each figure's README.

## Data availability (accessions)

| Dataset | Accession / access |
|---|---|
| Sequencing generated here | NCBI GEO **GSE339543** |
| T2T-CHM13v2.0 | GenBank GCA_009914755.4 |
| GIAB HG002 | NCBI BioProject PRJNA200694 |
| HPRC (Year 1 + Release 2) assemblies + ONT/HiFi | PRJNA730823 |
| CPC Phase I / T2T-YAO / CN1 | POG portal; GWH / PRJCA016397 |
| Fiber-seq (FIRE) | SRA SRR16356599 (+ broadly-consented reprocessing) |
| UK Biobank | Application 98772 (controlled) |
| GTEx v9 WGS / v11 RNA-seq | dbGaP phs000424 (project #43672) |
| TCGA | dbGaP phs000178 (project #44042) |
| CPTAC | dbGaP phs001287 (project #44042) |

## Layout

```
env/                        captured conda environments (see env/README.md)
refs/                       consensus repeat unit, feature map, gene sets (public inputs)
db/                         5S_rDNA.db schema + build/query scripts (DB itself not shipped)
fig1_array_geneconversion/  Fig 1 + S1
fig2_shortread_ukbb/        Fig 2 + S2  (short-read caller + UK Biobank)
fig3_methylation_fiberseq/  Fig 3 + S3  (CpG methylation + Fiber-seq accessibility + 45S/RNU2)
fig4_gtex_expression/       Fig 4 + S4  (GTEx variant expression + DESeq2/GSEA)
fig5_saturation_mutagenesis/Fig 5 + S5  (SUNi library, fraction-seq, RNAfold)
fig6_cancer_selection/      Fig 6 + S6  (CPTAC/TCGA p53 + cross-cohort selection)
building_blocks/            the reusable methods, generalized
CODE_MANIFEST.md            panel -> script index (start here)
```

## Software environment

Analyses used the core analysis environment `5s_pipeline` (provided as two captured snapshots
with minor tool-version differences), `fiberseq` for the Fiber-seq tools, and `bioinfo_pipeline`
for the saturation-mutagenesis / functional-assay analysis. The exact, versioned exports are in
`env/`. Core tools: Python 3.11 (NumPy, pandas, SciPy, statsmodels, pysam, pydeseq2, gseapy),
samtools/bcftools, bwa, minimap2, BLAST+, MAFFT, modkit; ViennaRNA + bowtie2/cutadapt for the
functional assay.

## Running the scripts

Scripts read paths from environment variables, with in-repo defaults:

| Variable | Meaning | Default |
|---|---|---|
| `FIVES_DB` | path to `5S_rDNA.db` (not distributed) | `5S_rDNA.db` |
| `FIVES_OUT` | output directory for figures/tables | `output` |
| `FIVES_DATA` | derived-data inputs (per-donor tables, pileups, exports) | `data` |
| `FIVES_REFS` | public reference inputs (`refs/`) | `refs` |
| `FIVES_BIN` | directory of external tool binaries, if not on `PATH` | (uses `PATH`) |

Example: `FIVES_DB=/path/to/5S_rDNA.db python fig1_array_geneconversion/31_variant_clustering_runs.py`
