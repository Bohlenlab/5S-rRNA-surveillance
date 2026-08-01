# Software environments

Exact, versioned conda exports captured 2026-08-01 from the two compute servers.

| File | Server | Used for |
|---|---|---|
| `5s_pipeline.immuno2.yml` | immuno2 | cancer GDC slicing + DE, HPRC HiFi/ONT streaming, 45S methylation (half 1), consensus rederivation |
| `5s_pipeline.trr237.yml` | trr237 ("server2") | UK Biobank 5S extraction+calling, GTEx DESeq2 run (`de_v2`), 45S methylation (half 2) |
| `fiberseq.trr237.yml` | trr237 | Fiber-seq (fibertools-rs, winnowmap, meryl) |

**The `5s_pipeline` env is not identical between servers** — the same env name carries different
tool versions (e.g. immuno2 samtools 1.23 / modkit 0.6.1 vs trr237 samtools 1.6 / modkit 0.4.1).
Both are archived; `CODE_MANIFEST.md` states which server (hence which snapshot) ran each analysis.
Python analysis libraries (pydeseq2 0.5.4, gseapy, statsmodels, scikit-learn) match across both.

To recreate an environment: `conda env create -f 5s_pipeline.trr237.yml`.

TODO: capture the wet-lab functional-assay env (`bioinfo_pipeline`: cutadapt, Trimmomatic,
bowtie2, samtools, ViennaRNA) from the local Mac — used for Figure 5.
