# Software environments

Exact, versioned conda exports.

| File | Used for |
|---|---|
| `5s_pipeline.yml` | core analysis environment (alignment, variant calling, methylation, DE) |
| `5s_pipeline.alt.yml` | second captured snapshot of the core environment, with minor tool-version differences |
| `fiberseq.yml` | Fiber-seq tools (fibertools-rs, winnowmap, meryl) |
| `bioinfo_pipeline.yml` | saturation-mutagenesis / functional-assay analysis (cutadapt, Trimmomatic, bowtie2, samtools, ViennaRNA) |

The core `5s_pipeline` environment is provided as two captured snapshots (`5s_pipeline.yml` and
`5s_pipeline.alt.yml`) that carry the same set of tools with minor version differences (e.g.
samtools and modkit). The Python analysis libraries (pydeseq2, gseapy, statsmodels, scikit-learn)
match across both.

To recreate an environment: `conda env create -f 5s_pipeline.yml`.
