# 5S rDNA paper — code archive manifest (DRAFT)

Sengl et al., "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
Scope: share **code + conceptual building blocks**, not data and not infrastructure.
Controlled data (UKBB, GTEx phs000424, TCGA phs000178, CPTAC phs001287, CPC) is never included —
the repo references accessions + access paths only. Panel-for-panel raw-data reproduction is NOT a goal;
the goal is that a competent bioinformatician can read each method and adapt it.

Legend: ★ = canonical script to archive · [Mac] local · [imm] immuno2 · [trr] trr237/server2 · [RAP] DNAnexus.

---

## 1. Proposed repo layout

```
5S_paper_code/
  README.md                     # overview, data-access statement, how to read this repo
  CODE_MANIFEST.md              # this file (panel -> script provenance)
  env/                          # captured conda envs (per server), + a curated env.yml
    5s_pipeline.immuno2.yml
    5s_pipeline.trr237.yml
    fiberseq.trr237.yml
  refs/                         # consensus unit FASTA, feature map, gene sets (public only)
  db/                           # schema + build/query scripts (NOT the populated DB)
  fig1_array_geneconversion/
  fig2_shortread_ukbb/
  fig3_methylation_fiberseq/
  fig4_gtex_expression/
  fig5_saturation_mutagenesis/
  fig6_cancer_selection/
  building_blocks/              # the ~8 reusable methods, lightly generalized + documented
```

Each `figN_*/` gets a short README that (a) lists the panels it makes, (b) names the input
tables/DB it reads, (c) states the data-access requirement. Hardcoded paths are lifted to a
`config` block at the top of each script (or a small `paths.py`).

---

## 2. Panel -> final script (pinned via `HPRC/manuscript_tables` READMEs)

NOTE: `manuscript_tables` uses an EARLIER figure-numbering than the V8 manuscript. Mapping is by
SCRIPT, not panel letter. Where they differ, V8 letters are given in ( ).

### Figure 1 — array structure & gene conversion   (all [Mac] `HPRC/scripts`)
| V8 panel | Analysis | ★ Script |
|---|---|---|
| 1B | per-copy variant map (HG002) | `26_array_maps.py` |
| 1C/S1A | copy number | `07_array_structure_analysis.py` |
| 1D,1I | variants per hap / per gene | `07_...` + `51_variants_per_haplotype_table.py` |
| 1E | copies-vs-donors per variant | `07_...` / catalog `52_variant_catalog_table.py` |
| 1F | substitution spectrum (normalized) | `39_substitution_normalized_cpg.py` |
| 1G/S1B | within-array clustering index | `31_variant_clustering_runs.py` (+ `47_clustering_single_panels.py`) |
| 1H | tract-size substitution spectrum (gBGC) | `48_gbgc_resolved.py` (+ `49_gbgc_trend_stats.py`) |
| 1J,1K | large-tract rate & formation by region | `50_geneconv_gene_alu_nts.py` |
| 1K SFS | content-adjusted regional SFS | `27_population_genetics.py` |
| S2/Table S2 | repeat-unit annotation | `44_repeat_unit_annotation.py` |
| DB build | consensus + catalog | `build_database.py`, `12_consensus_rederive.py`, multi-contig `57–60b` |

### Figure 2 — short-read caller + UKBB   (SPLIT across [Mac]/[trr]/[imm]/[RAP])
| V8 panel | Analysis | ★ Script | Where |
|---|---|---|---|
| 2A,2B | SR vs assembly VAF; HiFi rescue | `75_platform_comparison.py` | [Mac] |
| 2C | calls-by-group vs VAF (sens 83/prec 69) | `74_summary_figure.py` | [Mac] |
| S2A,B | pseudogene specificity (sim reads) | `figS_pseudogene.py` + `WGS-Variant-Identification/` sim | [Mac] |
| S2C | HPRC SR coverage | `13c_hprc_sr_coverage_S2C.py` | [Mac] |
| pipeline | UKBB read extraction (full CRAMs) | **RAP applet — TO RETRIEVE/DESCRIBE** | [RAP] |
| pipeline | extraction→align→variant→trio-QC | `5S_setup/{0,1,2,4,6}_*` | ★[trr] |
| 2D | carrier landscape | `78a_panel_A_export.py` | [Mac] |
| 2F(=t.2E) | UKBB vs assembly VAF concordance | `93_ukbb_assembly_consistency.py` | [Mac]/[imm] |
| 2G(=t.2F) | twin/parent-child sharing | `94_heritability_vaf.py` | [Mac] |
| S2D | Falconer midparent heritability | `94_...` (falconer mode) | [Mac] |
| 2H(=t.2G) | gene-conversion by UKBB VAF | `93b_ukbb_assembly_consistency_all.py` | [Mac] |
| 2I(=t.2H) | variants per UKBB carrier | carriers table (`76_import_ukbb_population.py`) | [imm] |
| 2J(=t.2I) | ICD-10 positional enrichment | `95_icd10_positional_enrichment.py` | [imm] |
| 2K,2L(=t.2J) | dosage-OR curves (OCD/SSc) | `85_dosage_or_curves.py` | [imm] |
| PheWAS core | ICD-10 dosage logistic regressions | `81_ukbb_phenotype_associations.py` | ★[imm] (vs `cohort_icd10.db`) |
| S2C landscape | UKBB variant landscape | `78_/79_ukbb_variant_landscape*.py` | [imm] |
| selection setup | UKBB purifying selection | `80_ukbb_purifying_selection.py` | [imm] |
SR-caller calibration backbone (not all panelled): `61–73` [Mac] (VAF sweep, filter tests, rescue nulls `fig2_rescue_*`).

### Figure 3 — methylation & accessibility   (SPLIT)
| V8 panel | Analysis | ★ Script | Where |
|---|---|---|---|
| 3A/S3A,B | array methylation vs position (ONT/HiFi) | `40_methylation_full215*.py` | [Mac]+[imm/trr] stream |
| 3B/S3C | per-copy methylation distribution | `40_...` / `41_methylation_overview*.py` | [Mac] |
| 3C,D(class)/S3F,G | within-copy positional by class | `40_methylation_full215.py` (★ authoritative) | [Mac] |
| 3E(=dosage) | dosage compensation (set-point) | `40_...` + `36_lowmeth_copies_donor_stability.py` | [Mac] |
| 3F,G,H | whole-copy meth by variant count | `60_variant_hypo_proportion*.py` (+`24b–g`) | [Mac] |
| 3(edge)/S3D,E | array-edge / border methylation | `62_/64_border|array_end...` , `65_/66_` | [Mac] |
| 3E,F(fiberseq)/S3H,I | Fiber-seq m6A accessibility | `fiberseq_5S/scripts/build_spanning_figure.py`, `within_copy_profile_pub.py` | ★[Mac]+[trr] stream |
| S3K,L | 45S NOR edge methylation | `45S_methylation/scripts/{00–04}` | ★[Mac]+[imm+trr] (both halves) |
| S3M | RNU2 own-assembly methylation | `rDNA_dosage_control/09.../RNU2/compute_rnu2_ownasm.py` | [Mac] |
| S3N | edge-biased gene-body variation | `nascent_edge_analysis/edge_analysis.py`,`edge_expression.py` | [Mac] |
| DB meth table | build `copy_methylation` | `37_/38_export|load_copy_methylation*.py`, `60_multicontig_methylation.py` | [Mac]+[imm] |

### Figure 4 — GTEx expression   ([Mac] `GTEx/` + ★[trr] `de_v2` run)
| V8 panel | Analysis | ★ Script |
|---|---|---|
| 4B | RNA-seq depth profile | `fig4_coverage_pooldata.py`+`fig4_coverage_panel.py` |
| 4C | per-variant carrier vs non-carrier | `fig4_panel4_pooldata.py`+`fig4_build_panels.py` |
| 4D | rank-skew volcano | `fig4_build_panels.py` (from `wgs_rna_rank_skew.py`) |
| 4E | expressed-variant distribution | `fig4_p3_detection.py` |
| 4F | prevalence (AUC q99) | `fig4_panel4_build.py`  *(DECISION: q99 37% vs strict 14% — pick one)* |
| 4G | DNA vs RNA VAF | `fig4_build_panels.py` |
| 4H | trans-effect volcano | `ex05j_p10_4x4_biogenesis.py` (DE from `ex05e_fullmodel_gsea.py`) |
| 4I,J,K,L | dose-response + GSEA + meta | `ex05g_dnarna_vaf_panel.py`, `de_v2/{de_common,de_meta,de_gsea}.py` |
| genotyping | GTEx WGS/RNA 5S VAF | `GTEx/scripts/{10–14,20,23}` |
DE executed run + logs: ★[trr] `~/de_v2/` (`meta/SUMMARY.tsv`). Superseded: `figures/de_expression`, `eq*`, `ex0[1-4]`, `vm_archive`.

### Figure 5 — saturation mutagenesis   ([Mac] dated folders + `de_v2`)
| V8 panel | Analysis | ★ Script |
|---|---|---|
| 5C–G | fraction-seq → expr/incorporation scores | `20251006 5S rRNA Sequencing experiment 1/` (`5S/`+`GFP/` stages → `process_rep_mutation_ratios.py`) |
| 5(region E/G) | region expr/incorp panels | `GTEx/de_v2/make_fig5_region.py` |
| S5F | RNAfold ΔΔG | `RNA-fold/regenerate_sense_folding.py` (+`make_folding_figure.py`) — use `*_CORRECTED.*`, NOT `Folding energies.csv` |
| load | functional scores → DB | `T2T/scripts/load_functional_annotation.py` |
| S5D | expression by ICR region | part of fraction-seq analysis |
Structure schematics (5D/F, S5E): `VIsualizing 5S rRNA/` (ChimeraX 8BGU) — not code-reproducible (documented).

### Figure 6 — cancer p53 + cross-cohort selection   (SPLIT)
| V8 panel | Analysis | ★ Script | Where |
|---|---|---|---|
| 6A–C | function-stratified trans-effects (IE) | `de_v2/make_figs3to6.py` (+`make_fig13_IE_*`) | [Mac]/[trr] |
| 6D–L | cancer p53-stratified DE, GSEA | `cancer_5S/surveillance_v2/` + `scripts/{03_slice_call,54–77}` | ★[imm] run → [Mac] figs |
| 6Q | somatic 5S gains × p53 | `cancer_5S/scripts/{100–111}` | [imm] |
| 6M–P | cross-cohort selection (HPRC/UKBB/GTEx) | `97_incorporation_depletion_three_cohort.py`, `98_binary_classifications_two_ways.py`, `2E_fair_perdonor_vaf.py`, `de_v2/make_figs3to6.py` | [Mac] |
| S6A | methylation by variant consequence | `24_gene_methylation_by_functional_consequence.py` (+`24b–g`) | [Mac] |
| robustness | substitution-covariate control | `97c_substitution_covariate_robustness.py` | [Mac] |
| support | MDM2–p53 structure | `MDM2/1RV1.cif` (asset) | [Mac] |

---

## 3. Reusable building blocks (the shareable "conceptual" core)
1. **Collapsed-array short-read caller** — map to a single consensus repeat unit, low-VAF calling, HiFi rescue, F1-calibrated threshold (`5S_setup` + `74/75/61–73`). *Explicitly pitched in the Discussion as a template for other collapsed arrays.*
2. **Multi-contig per-copy array reconstruction & ordering** (`57–60b`).
3. **Runs-based within-array clustering index** (Wald–Wolfowitz) (`31`/`47`).
4. **Content-adjusted regional SFS** (`27`).
5. **Edge-distance methylation/accessibility profiling on collapsed arrays** (competitive anchoring) (`40`, `45S_methylation`, `fiberseq_5S`).
6. **Fraction-seq → per-variant expression/incorporation scoring** (`process_rep_mutation_ratios.py`).
7. **Variant-dose → per-tissue DESeq2 + cross-tissue meta + GSEA** (`de_v2`).
8. **p53-stratified genomic-dose → RP-repression interaction** (`cancer_5S`).

---

## 4. Environments (captured 2026-08-01)
- `5s_pipeline` differs by server — archive BOTH snapshots and note which analysis used which:
  - immuno2: samtools 1.23, bcftools 1.21, bwa 0.7.19, minimap2 2.30, modkit 0.6.1, pydeseq2 0.5.4, gseapy — used for cancer GDC, HPRC streaming, 45S half1.
  - trr237: samtools 1.6, bcftools 1.9, bwa 0.7.18, minimap2 2.28, modkit 0.4.1, pydeseq2 0.5.4, gseapy 1.3.0 — used for UKBB extraction, GTEx de_v2, 45S half2, fiberseq.
- `fiberseq` (trr237): fibertools-rs 0.9.0, winnowmap 2.03.
- Wet-lab (Fig 5): `bioinfo_pipeline` (cutadapt/Trimmomatic/bowtie2/samtools) + ViennaRNA — env not yet captured (local Mac).

---

## 5. Excluded (present in tree, NOT in paper)
Evo_5S (cross-species); all Hi-C (`HiC/`, `hic_analysis`, `scripts/hic`); most `rDNA_dosage_control` (08 tRNA-chr6, 01 chromatin/H3K27me3, 03/07 planning); REH (confirmed not used); pilots (`20250603`, `20250915`, `T2T/README` CHM13/HG002 pilot, GTEx `de_expression`/`eq*`/`vm_archive`, cancer v1 `figures/`+`results_*`, Trios/Trios_2.0, `Folding energies.csv`).

---

## 6. Open decisions / gaps
- **[D1] UKBB RAP extraction code** — retrieve the DNAnexus applet/script, or describe in Methods only? (It's the one missing code step.)
- **[D2] Fig 4F prevalence number** — 37% (AUC q99) vs 14% (strict); the Figure-4 README flags "pick one before submission."
- **[D3] DB** — ship schema+build scripts only (agreed); confirm no populated DB shipped.
- **[D4] Path handling** — lift hardcoded `/home/j.bohlen`, `/home/bohlen_lab`, `/Users/bohlen` (237/68/336 files) to per-script config blocks during curation.
- **[D5] HPRC/scripts README** only indexes 01–56; 57–99 (UKBB, SR calibration, selection) undocumented — extend during curation.
- **[D6] Wet-lab env** (`bioinfo_pipeline`) still to capture from the Mac.
- **[D7] Manual steps** (Prism/Excel/Affinity) behind a few panels — document as "assembled in GraphPad/Illustrator" in the figure READMEs.
