# 5S_rDNA.db — Database Documentation

**File:** `5S_rDNA.db` (SQLite)

---

## 1. Purpose

A hierarchical relational database of 5S rDNA tandem repeat arrays and their per-copy sequence variants across human T2T-quality assemblies. Positions are numbered on the 2168 bp 5S rDNA repeat unit. The primary variant reference is the **population-majority consensus** of that unit (the `consensus_t2t` polarization): at each position the reference allele is the base carried by most copies across all haplotypes. This population consensus provides the coordinate scaffold and differs from the raw CHM13 5S unit at 8 positions.

Current content (as of 2026-05-31):

| Table | Rows |
|---|---|
| `array_reference` | 1 |
| `assembly` | 50 |
| `haplotype` | 97 |
| `copy` | 8,090 |
| `variant` | 133,750 |
| `read_variant` | 1,738 |

---

## 2. Data Sources

### 2.1 HPRC Year 1 (`cohort = "HPRC_Year1"`)

47 individuals from the Human Pangenome Reference Consortium Year 1 release (Liao et al. 2023, *Nature*). All assemblies are haplotype-resolved diploid (hap1/hap2) produced with Hifiasm from PacBio HiFi reads. Assemblies are available at `s3://human-pangenomics/`.

- **Assembly source:** `s3://human-pangenomics/working/HPRC/{SAMPLE}/assemblies/year1_f1_assembly_v2_genbank/`
- **HiFi reads:** `s3://human-pangenomics/working/HPRC/{SAMPLE}/raw_data/PacBio_HiFi/`  
  Three samples use HPRC_PLUS bucket: HG002, HG005, HG00733
- **Illumina reads:** `s3://human-pangenomics/working/HPRC/{SAMPLE}/raw_data/Illumina/child/{SAMPLE}.final.cram`  
  GRCh38-aligned CRAMs; chr1 5S region extracted before re-alignment to T2T consensus
- **Metadata inventory:** `HPRC/data_inventory.tsv`

Populations covered: AFR, AMR, EAS, EUR, SAS (population codes follow 1000 Genomes Project).

### 2.2 CHM13 (`cohort = "CHM13"`)

The CHM13 T2T reference genome (v2.0, Nurk et al. 2022, *Science*). Haploid hydatidiform mole cell line. The 5S array sits on chr1 at approximately chr1:227,700,000–228,700,000.

- **Assembly:** `chm13v2.0.fa.gz` (local)
- **Array database:** `databases/5S_array_database.tsv`  
  123 copies on chr1, all minus strand, BLAST-filtered (≥95% identity, ≥115 bp)

### 2.3 HG002_GIAB (`cohort = "HG002_GIAB"`)

The GIAB benchmark individual HG002 (Ashkenazi Jewish male). This entry uses the **hg002v1** assembly (older diploid assembly from the HG002 genome paper), distinct from the HPRC Year 1 HG002 assembly (sample_id `"HG002"`) which is also in the database.

- **Assembly:** hg002v1 (MATERNAL / PATERNAL haplotypes)
- **Array database (MATERNAL = hap1):** `databases/5S_array_database_HG002_MATERNAL.tsv` — 104 copies
- **Array database (PATERNAL = hap2):** `databases/5S_array_database_HG002_PATERNAL.tsv` — 69 copies
- **SR read variants:** BGIseq 150 bp, aligned to T2T consensus (bcftools mpileup, ~2,750× depth)

---

## 3. Processing Pipelines

### 3.1 Assembly Analysis

Runs on per-haplotype assembly FASTAs. Produces per-copy databases that populate the `haplotype`, `copy`, and `variant` tables.

**Steps:**
1. **BLAST** — `blastn` with the T2T consensus as query against the assembly FASTA  
   Filter: ≥95% identity (`-perc_identity 90`), ≥115 bp alignment length  
   Format: 15 columns including `sstrand` (the source of the `strand` field in `haplotype`)
2. **Sequence extraction** — `samtools faidx` to extract, per copy:
   - Full ~2168 bp unit (gene lo−1420 to gene hi+629)
   - 5S gene region (119 bp, from BLAST hit coordinates)
   - NTS-pre: 629 bp immediately 3′ of the gene in assembly coordinates (= 5′ NTS in sense-strand orientation for minus-strand arrays)
   - NTS-post: 1419 bp immediately 5′ of the gene in assembly coordinates
3. **MAFFT alignment** — `mafft --auto --quiet --thread 4` run separately for three sequence sets: full unit, NTS-pre, NTS-post
4. **Variant calling** — Majority-vote consensus across all copies in the alignment; each copy's deviations from that consensus are recorded as `pos:ref>alt` strings  
   `pos` is 0-based count of non-gap characters in the alignment row up to that column

### 3.2 HiFi LR Pipeline

Runs on PacBio HiFi reads. Produces `read_variant` rows with `modality = "hifi"`.

**Steps:**
1. Stream HiFi reads from S3 (BAM or FASTQ) or SRA
2. Extract 5S-aligning reads: `minimap2 -ax map-hifi` against T2T consensus, MAPQ ≥ 20
3. Chop reads into 2168 bp chunks (MIN_CHUNK = 800 bp) to normalize coverage across the full unit
4. Re-align chunks: `minimap2 -ax map-hifi --secondary=no | samtools view -F 0x904 | samtools sort`
5. `bcftools mpileup -d 10000000 -q 0 -Q 10 | bcftools query -f '%POS\t%REF\t%ALT[\t%DP\t%AD]\n'`
6. Filter: AD ≥ 3, VAF ≥ 2%

Positions from bcftools are 1-based — stored directly with no coordinate conversion.

### 3.3 Illumina SR Pipeline

Runs on Illumina reads. Produces `read_variant` rows with `modality = "illumina"`.

**Steps:**
1. Stream the GRCh38-aligned CRAM from S3, extract reads overlapping the 5S region (`chr1:228,400,000–229,300,000`)
2. Convert to interleaved FASTQ: `samtools collate | samtools fastq`
3. QC: `fastp --interleaved_in -A -Q -L` (adapter trimming, no quality/length filter)
4. Align to T2T consensus: `bwa mem -B 6 -O 8 -L 5,5 -T 30` (clipping-penalised for tandem repeats)
5. Filter: `samtools view -q 30 -F 0x904`
6. `bcftools mpileup -d 10000000 -q 30 -Q 30 | bcftools query`

---

## 4. Schema

### Relationships

```
array_reference (1)
       │
       └──< haplotype (97)
                │
assembly (50) ──┤
                │
                └──< copy (8,090)
                           │
                           └──< variant (133,750)

assembly (50) ──< read_variant (1,738)
```

---

### 4.1 `array_reference`

One row. Stores the reference repeat unit used as the alignment target for all analyses.

| Column | Type | Value / Description |
|---|---|---|
| `array_id` | INTEGER PK | Always 1 |
| `name` | TEXT | `"5S_rDNA"` |
| `ref_label` | TEXT | Contig name of the consensus reference FASTA used in the alignment commands |
| `sequence` | TEXT | Full 2168 bp T2T consensus nucleotide sequence |
| `length_bp` | INTEGER | 2168 |
| `nts_pre_start` | INTEGER | 1 (1-based, inclusive) |
| `nts_pre_end` | INTEGER | 629 |
| `gene_start` | INTEGER | 630 — first base of the 5S rRNA gene |
| `gene_end` | INTEGER | 748 — last base of the 5S rRNA gene |
| `nts_post_start` | INTEGER | 749 |
| `nts_post_end` | INTEGER | 2168 |

The 5S rRNA gene occupies positions 630–748 (119 bp). Positions 1–629 are the Non-Transcribed Spacer upstream of the gene (NTS-pre in sense-strand orientation), and 749–2168 are the NTS downstream (NTS-post).

---

### 4.2 `assembly`

One row per individual.

| Column | Type | Description |
|---|---|---|
| `assembly_id` | INTEGER PK | Auto-increment |
| `sample_id` | TEXT UNIQUE | Identifier (e.g. `"HG00438"`, `"CHM13"`) |
| `cohort` | TEXT | Dataset of origin (see §5.1) |
| `population` | TEXT | 1000 Genomes population code (e.g. `"CHS"`, `"YRI"`). NULL for non-HPRC samples |
| `superpopulation` | TEXT | Five-letter superpopulation (e.g. `"EAS"`, `"AFR"`). NULL for non-HPRC |
| `sex` | TEXT | `"M"` / `"F"` / NULL. Currently populated only for CHM13 (F) and HG002_GIAB (M) |
| `age` | INTEGER | NULL (not yet populated) |
| `assembly_hap1_s3` | TEXT | S3 URI of hap1 assembly FASTA. NULL for non-HPRC |
| `assembly_hap2_s3` | TEXT | S3 URI of hap2 assembly FASTA. NULL for non-HPRC |
| `has_hifi` | INTEGER | 1 if HiFi reads were available and processed |
| `has_illumina` | INTEGER | 1 if Illumina reads were available and processed |
| `has_methylation` | INTEGER | 1 if HiFi kinetics (MM/ML tags) were available for methylation analysis |
| `has_rnaseq` | INTEGER | 1 if RNA-seq data is available |

---

### 4.3 `haplotype`

One row per (individual, haplotype) pair. For diploid assemblies there are two rows per individual; for haploid assemblies (CHM13) one row.

| Column | Type | Description |
|---|---|---|
| `haplotype_id` | INTEGER PK | Auto-increment |
| `assembly_id` | INTEGER FK | References `assembly` |
| `array_id` | INTEGER FK | References `array_reference` (always 1) |
| `hap_label` | TEXT | `"hap1"` or `"hap2"`. For HPRC: hap1 = paternal, hap2 = maternal (per hifiasm convention, but not guaranteed) |
| `array_chrom` | TEXT | Contig / chromosome name in the source assembly (e.g. `"HG00438#1#JAHBCB010000053.1"`, `"chr1"`) |
| `strand` | TEXT | `"plus"` or `"minus"` — majority `sstrand` from BLAST output file. All 97 haplotypes are minus-strand |
| `n_copies` | INTEGER | Number of BLAST-filtered repeat copies (≥95% identity, ≥115 bp) |
| `array_start_local` | INTEGER | Leftmost assembly coordinate across all copies (1-based) |
| `array_end_local` | INTEGER | Rightmost assembly coordinate |

UNIQUE constraint on `(assembly_id, hap_label)`.

**Copy count range:** 28–152, mean 83 (across 97 haplotypes).

---

### 4.4 `copy`

One row per numbered repeat copy within a haplotype. Copies are numbered 1-based in 5′→3′ order after strand normalisation (i.e., copy 1 is the most 5′ copy in the sense direction of the haplotype).

| Column | Type | Description |
|---|---|---|
| `copy_id` | INTEGER PK | Auto-increment (internal DB key, not the original copy number) |
| `haplotype_id` | INTEGER FK | References `haplotype` |
| `copy_number` | INTEGER | 1-based, 5′→3′ in sense-strand orientation |
| `unit_start_local` | INTEGER | Start of the full repeat unit in assembly coordinates (1-based) |
| `unit_end_local` | INTEGER | End |
| `unit_length_bp` | INTEGER | Length of the full unit (typically ~2168 bp; can vary especially at array borders) |
| `spacing_to_next_bp` | REAL | Distance to the next copy's unit start. NULL for the last copy in the array |
| `gene_lo_local` | INTEGER | 5S gene start in assembly coordinates (from BLAST `min(sstart, send)`) |
| `gene_hi_local` | INTEGER | 5S gene end |
| `gene_pct_identity` | REAL | BLAST percent identity of this copy to the T2T consensus (range ~95–100%) |
| `gene_mismatches` | INTEGER | BLAST mismatch count |
| `gene_gaps` | INTEGER | BLAST gap-open count |
| `n_snv_gene` | INTEGER | Number of variants called in the full-unit MAFFT alignment (covers all three sub-regions) |
| `n_snv_5s_gene` | INTEGER | Subset of `n_snv_gene` that fall in the 5S gene region (consensus positions 630–748). This is the primary measure of gene sequence divergence |
| `n_snv_nts_pre` | INTEGER | Variants in the separate NTS-pre MAFFT alignment |
| `n_snv_nts_post` | INTEGER | Variants in the separate NTS-post MAFFT alignment |
| `category` | TEXT | `"identical"` if `n_snv_5s_gene = 0` (no 5S gene variants); `"highly_similar"` otherwise. For legacy HPRC databases this was originally BLAST-pident-based and was corrected to the gene-variant definition |
| `border_note` | TEXT | Array position annotation (see §5.4) |

UNIQUE constraint on `(haplotype_id, copy_number)`.

---

### 4.5 `variant`

One row per variant per copy, derived from MAFFT multiple-sequence alignments. These are **assembly-derived** variants — differences between a specific copy and the MAFFT majority-vote consensus of all copies in the haplotype.

| Column | Type | Description |
|---|---|---|
| `variant_id` | INTEGER PK | Auto-increment |
| `copy_id` | INTEGER FK | References `copy` |
| `alignment_source` | TEXT | Which MAFFT alignment produced this variant (see §5.2) |
| `consensus_pos` | INTEGER | **1-based position in the T2T consensus unit** (1–2168) |
| `ref` | TEXT | Reference base(s) per the MAFFT alignment consensus. For the full-unit alignment this closely approximates the CHM13 5S unit; minor deviations are possible at hypervariable NTS positions |
| `alt` | TEXT | Alternate base(s) observed in this copy |
| `region` | TEXT | `"nts_pre"`, `"gene"`, or `"nts_post"` (see §5.3) |

**Important:** The `ref` is the MAFFT majority-vote consensus across all copies in the haplotype, not necessarily the population consensus (or the CHM13 5S unit). At NTS positions with high copy-to-copy variation, these may differ.

**Coordinate conversion applied during import:**

| `alignment_source` | 0-based source pos | 1-based consensus_pos stored |
|---|---|---|
| `gene_unit` (HPRC) | Position in full 2168 bp unit alignment | `pos_0based + 1` |
| `gene_unit` (legacy: CHM13, HG002_GIAB) | Position in 119 bp gene alignment only | `630 + pos_0based` |
| `nts_pre_aln` | Position in NTS-pre alignment (629 bp) | `pos_0based + 1` |
| `nts_post_aln` | Position in NTS-post alignment (1419 bp) | `749 + pos_0based` |

---

### 4.6 `read_variant`

One row per variant detected from sequencing reads via bcftools mpileup. Not copy-assignable (reads align to the consensus and cannot be phased to individual copies). Links to `assembly` rather than `haplotype`.

| Column | Type | Description |
|---|---|---|
| `read_variant_id` | INTEGER PK | Auto-increment |
| `assembly_id` | INTEGER FK | References `assembly` |
| `modality` | TEXT | `"hifi"`, `"illumina"`, or `"rnaseq"` (see §5.5) |
| `consensus_pos` | INTEGER | **1-based position in the T2T consensus unit**, directly from bcftools. No coordinate conversion applied |
| `ref` | TEXT | Reference base at this position per the T2T consensus |
| `alt` | TEXT | Alt base observed in reads |
| `depth` | INTEGER | Total read depth at this position |
| `alt_depth` | INTEGER | Number of reads supporting the alt allele |
| `vaf` | REAL | `alt_depth / depth` |
| `region` | TEXT | `"nts_pre"`, `"gene"`, or `"nts_post"` |

**Current coverage:**

| Modality | Rows | Samples |
|---|---|---|
| `hifi` | 149 | HG01891, HG02257 |
| `illumina` | 1,589 | HG002_GIAB (BGIseq 150 bp) |

Additional HPRC HiFi and Illumina `read_variant` rows are added as the read-based pipelines complete.

---

## 5. Enumerated Values and Labels

### 5.1 `assembly.cohort`

| Value | Description |
|---|---|
| `"HPRC_Year1"` | 47 HPRC Year 1 individuals (Liao et al. 2023) |
| `"CHM13"` | CHM13 T2T reference (Nurk et al. 2022) |
| `"HG002_GIAB"` | HG002 hg002v1 GIAB assembly (distinct from HPRC HG002) |

### 5.2 `variant.alignment_source`

Each copy's variants are stored under a labelled `alignment_source`. The **primary** source used throughout the study is `consensus_t2t` (variants polarized against the population-majority consensus).

| Value | Reference | Note |
|---|---|---|
| `"consensus_t2t"` | Population-majority consensus of the 2168 bp unit | **Primary.** Reference allele = the base most copies carry across all haplotypes; used for the population-genetic and gene-conversion analyses. |
| `"gene_unit_t2t"` | Raw CHM13 T2T 5S unit | Variants relative to the CHM13 unit (differs from the population consensus at 8 positions). |
| `"gene_unit"` | Per-haplotype MAFFT consensus | Full ~2168 bp repeat-unit alignment; variants relative to each haplotype's own majority-vote consensus. |
| `"nts_pre_aln"` | Per-haplotype MAFFT, NTS-pre only (pos 1–629) | Higher NTS-pre sensitivity than `gene_unit`. |
| `"nts_post_aln"` | Per-haplotype MAFFT, NTS-post only (pos 749–2168) | Higher NTS-post sensitivity. |

The region-specific MAFFT alignments (`nts_pre_aln`, `nts_post_aln`) are more sensitive for NTS variant discovery because the full-unit alignment is dominated by the conserved gene region.

### 5.3 `variant.region` / `read_variant.region`

Assigned from `consensus_pos`:

| Value | Position range | Description |
|---|---|---|
| `"nts_pre"` | 1–629 | Non-Transcribed Spacer upstream of the gene (5′ NTS in sense orientation) |
| `"gene"` | 630–748 | 5S rRNA gene (119 bp) |
| `"nts_post"` | 749–2168 | Non-Transcribed Spacer downstream of the gene (3′ NTS) |

### 5.4 `copy.border_note`

| Value | Description |
|---|---|
| `"interior"` | A copy not at the array boundary. NTS coordinates reflect genuine inter-copy spacer sequence |
| `"5-prime_array_border"` | The outermost copy at the 5′ end of the array. Its NTS-pre window extends into flanking non-repetitive genomic sequence, so NTS-pre variant counts are inflated (~109 variants vs ~3 for interior copies) |
| `"3-prime_array_border"` | The outermost copy at the 3′ end. NTS-post window extends into flanking sequence (~289 variants vs ~3 for interior copies) |

Border copies should be excluded from analyses of NTS variant rates. The edge inflation drops to interior levels by copy 3.

### 5.5 `read_variant.modality`

| Value | Read type | Alignment | Typical depth |
|---|---|---|---|
| `"hifi"` | PacBio HiFi CCS (>99.9% accuracy) | minimap2 `map-hifi` after 2168 bp chunking | ~30–60× per copy × copy number |
| `"illumina"` | Illumina paired-end (150–250 bp) | bwa mem | ~30× WGS × copy number |
| `"rnaseq"` | Illumina RNA-seq | bwa mem to gene-only reference | Variable (5S rRNA is highly expressed) |

### 5.6 `haplotype.strand`

All 97 haplotypes currently have `strand = "minus"`, meaning the 5S array is on the minus strand of the assembly contig. The extracted sequences are reverse-complemented by `samtools faidx -i` before alignment, so all stored coordinates and variants are in sense-strand (5′→3′ NTS-pre → gene → NTS-post) orientation.

### 5.7 `copy.category`

| Value | Definition | Count |
|---|---|---|
| `"identical"` | `n_snv_5s_gene = 0` — no variants in the 119 bp 5S gene relative to the MAFFT consensus | 6,964 (86%) |
| `"highly_similar"` | `n_snv_5s_gene > 0` — at least one variant in the 5S gene | 1,126 (14%) |

Note: these categories reflect variation relative to the per-haplotype MAFFT consensus, not the population consensus. A copy that is "identical" may still differ from the population consensus if the consensus differs.

---

## 6. Coordinate System

All positions in the database are **1-based, relative to the T2T consensus repeat unit (2168 bp)**.

```
Position:   1          629|630         748|749                       2168
            ├──────────────┤├──────────────┤├──────────────────────────┤
Region:         NTS-pre         5S gene             NTS-post
                (629 bp)        (119 bp)             (1420 bp)

Sub-regions (approximate positions):
  NTS-pre:
    CA repeats        1–138
    Spacer/promoter   139–622
    UPE               623–629

  5S gene:
    Entire gene       630–748

  NTS-post:
    Terminator        749–790
    Alu-like repeat   791–1230
    Spacer            1231–1854
    CTTCAA/TC-rich    1855–2168
```

---

## 7. Known Limitations and Caveats

**MAFFT consensus vs the population consensus.** Variants in the `variant` table are relative to the MAFFT majority-vote consensus of each haplotype's copies, not the population consensus directly. At highly variable NTS positions, the MAFFT consensus may differ from the population consensus by 1–2 bases, making `ref` in the variant table unreliable as an absolute reference. For gene-region variants this is not a concern (the gene is conserved enough that the MAFFT consensus matches the population consensus at all positions). For polarized, absolute-reference variants use the `consensus_t2t` source.

**Border copies.** The first and last copy of each array have inflated NTS variant counts because the extraction window extends into non-repetitive flanking sequence. Exclude `border_note IN ('5-prime_array_border', '3-prime_array_border')` from NTS variant rate analyses.

**Read variants are diploid.** `read_variant` entries represent the diploid mixture of both haplotypes (reads cannot be assigned to specific copies or haplotypes). VAF therefore reflects the frequency of an alt allele across all copies in both haplotypes combined.

**Strand is uniformly minus.** All 97 haplotypes have minus-strand arrays. If a plus-strand array is ever encountered, the coordinate convention would remain identical (samtools faidx -i reverse-complements before alignment), but `copy_number` ordering would need verification.

**HG002_GIAB vs HG002 (HPRC).** Both sample_ids refer to the same biological individual but different assemblies (hg002v1 GIAB vs HPRC Year 1 hifiasm). Copy numbers differ (GIAB: 104 MAT, 69 PAT; HPRC: check haplotype table). Do not merge these rows without accounting for the assembly difference.

**Legacy coordinate encoding for gene_unit variants.** For CHM13 and HG002_GIAB, `gene_unit` variants have `alignment_source = "gene_unit"` but were generated from a gene-only (119 bp) alignment, not a full-unit alignment as in HPRC. The `consensus_pos` values are correctly converted to T2T coordinates in both cases, but the number of variants per copy may not be directly comparable across cohorts for NTS positions in `gene_unit`.

**Incomplete read_variant table.** As of the current build, only a subset of assemblies have HiFi and Illumina read variants. The table is repopulated as the read-based pipelines complete for the remaining samples.

---

## 8. Example Queries

```sql
-- Copy number per haplotype for all HPRC samples
SELECT a.sample_id, a.superpopulation, h.hap_label, h.n_copies
FROM haplotype h JOIN assembly a USING(assembly_id)
WHERE a.cohort = 'HPRC_Year1'
ORDER BY h.n_copies DESC;

-- Gene variants at a specific position across all haplotypes
SELECT a.sample_id, h.hap_label, c.copy_number, v.ref, v.alt
FROM variant v
JOIN copy c USING(copy_id)
JOIN haplotype h USING(haplotype_id)
JOIN assembly a USING(assembly_id)
WHERE v.consensus_pos = 714 AND v.region = 'gene'
ORDER BY a.sample_id, c.copy_number;

-- Fraction of copies with 5S gene variants, by sample
SELECT a.sample_id,
       COUNT(*) AS total_copies,
       SUM(CASE WHEN c.n_snv_5s_gene > 0 THEN 1 ELSE 0 END) AS copies_with_gene_variants,
       ROUND(100.0 * SUM(CASE WHEN c.n_snv_5s_gene > 0 THEN 1 ELSE 0 END) / COUNT(*), 1) AS pct
FROM copy c
JOIN haplotype h USING(haplotype_id)
JOIN assembly a USING(assembly_id)
WHERE c.border_note = 'interior'
GROUP BY a.sample_id
ORDER BY pct DESC;

-- HiFi read variants in the 5S gene for one sample
SELECT rv.consensus_pos, rv.ref, rv.alt, rv.alt_depth, rv.vaf
FROM read_variant rv JOIN assembly a USING(assembly_id)
WHERE a.sample_id = 'HG01891' AND rv.modality = 'hifi' AND rv.region = 'gene'
ORDER BY rv.consensus_pos;

-- Variants shared between CHM13 and any HPRC haplotype (exact position + alt)
SELECT v1.consensus_pos, v1.ref, v1.alt,
       a2.sample_id AS hprc_sample, h2.hap_label, COUNT(c2.copy_id) AS n_copies
FROM variant v1
JOIN copy c1 USING(copy_id)
JOIN haplotype h1 USING(haplotype_id)
JOIN assembly a1 USING(assembly_id)
JOIN variant v2 ON v2.consensus_pos = v1.consensus_pos AND v2.alt = v1.alt
JOIN copy c2 ON v2.copy_id = c2.copy_id
JOIN haplotype h2 ON c2.haplotype_id = h2.haplotype_id
JOIN assembly a2 ON h2.assembly_id = a2.assembly_id
WHERE a1.sample_id = 'CHM13' AND a2.cohort = 'HPRC_Year1'
GROUP BY v1.consensus_pos, v1.alt, a2.sample_id, h2.hap_label
ORDER BY v1.consensus_pos;
```
