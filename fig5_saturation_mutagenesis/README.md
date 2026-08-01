# Figure 5 + S5 — saturation mutagenesis of the 5S rRNA gene

Functional assay generated here (raw reads: GEO GSE339543). Definitive experiment folder:
`20251006 5S rRNA Sequencing experiment 1/`.

| Panel | Script |
|---|---|
| 5C–G fraction-seq → expression + 60S-incorporation scores | `5S/` + `GFP/` stages (cutadapt demux → Trimmomatic → bowtie2 → `4_count/psysam-script*.py` → `5_analyse/summarize script.py`) → `process_rep_mutation_ratios.py` (GFP error subtraction) |
| 5 region panels (E/G) | `GTEx/de_v2/make_fig5_region.py` |
| S5D expression by ICR region | fraction-seq analysis output |
| S5F RNAfold ΔΔG | `RNA-fold/regenerate_sense_folding.py` + `make_folding_figure.py` |
| load scores → DB | `T2T/scripts/load_functional_annotation.py` |

Environment: `bioinfo_pipeline` (cutadapt, Trimmomatic, bowtie2, samtools) + ViennaRNA (see `../env/`).

Cautions:
- Use the terminal-trimmed counting branch (`4_count_trim/psysam-script-noterminal.py`, `--snv-end-trim 3`).
- RNAfold: the scripts fold the sense (transcript) strand.
- Final score aggregation (Table S8) and structure schematics (ChimeraX, `VIsualizing 5S rRNA/`)
  involve manual steps — documented, not fully code-reproducible.
