#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# 2_align-sort-analyze.sh — align filtered 5S reads to the consensus repeat with
# bwa, sort, and call per-position variants with bcftools.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
# =============================================================================
# 2_align-sort-analyze.sh — Align 5S reads to consensus, sort, call variants.
#
# Defaults derive from this script's location (scripts/ inside the project).
# Override any path by exporting the corresponding env var before running.
#
# Resume mode: export RESUME=1 to skip samples whose .tsv already exists.
# =============================================================================
set -euo pipefail

_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
_PIPELINE_DIR="$(dirname "$_SCRIPT_DIR")"

# ---- Paths (env-var overridable) ----
FASTQ_DIR="${FASTQ_DIR:-$_PIPELINE_DIR/results/filtered_fastq_5S}"
CONS_INDEX="${CONS_INDEX:-}"   # must be set via env var on each machine
ALIGN_DIR="${ALIGN_DIR:-$_PIPELINE_DIR/results/alignments}"
VCF_DIR="${VCF_DIR:-$_PIPELINE_DIR/results/vcf}"
TSV_DIR="${TSV_DIR:-$_PIPELINE_DIR/results/tsv}"

if [[ -z "$CONS_INDEX" ]]; then
    echo "[ERROR] CONS_INDEX environment variable is not set."
    echo "        Point it to the 5S_repeats_consensus.fa file."
    exit 1
fi

mkdir -p "$ALIGN_DIR" "$TSV_DIR"

SUMMARY_FILE="$ALIGN_DIR/alignment_summary.tsv"

# ---- Auto-detect cores (Linux: nproc, macOS: sysctl) ----
if command -v nproc &>/dev/null; then
    CORES="$(nproc)"
else
    CORES="$(sysctl -n hw.ncpu)"
fi
RESERVE=${RESERVE_CORES:-1}
USE=$(( CORES - RESERVE ))
BWA_T=${BWA_T_OVERRIDE:-2}
SAM_T=${SAM_T_OVERRIDE:-2}
BCF_T=${BCF_T_OVERRIDE:-2}
THREADS_PER_JOB=$(( BWA_T + SAM_T + BCF_T ))
JOBS=${JOBS_OVERRIDE:-$(( USE / THREADS_PER_JOB ))}
(( JOBS < 1 )) && JOBS=1

# ---- Alignment parameters (env-var overridable for indel-aware mode) ----
# Defaults are tuned for SNV calling. For indel calling, use:
#   BWA_B=4 BWA_O=6 BWA_MINSCORE=20 MAPQ_FILTER=0
BWA_B=${BWA_B:-6}
BWA_O=${BWA_O:-8}
BWA_MINSCORE=${BWA_MINSCORE:-30}
MAPQ_FILTER=${MAPQ_FILTER:-30}    # set to 0 to disable

echo "Cores: $CORES (reserving $RESERVE) → using ~${USE}"
echo "JOBS=$JOBS, BWA_T=$BWA_T, SAM_T=$SAM_T, BCF_T=$BCF_T"
echo "BWA params: -B $BWA_B -O $BWA_O -T $BWA_MINSCORE  MAPQ_FILTER=$MAPQ_FILTER"
echo "Input dir:  $FASTQ_DIR"
echo "Consensus:  $CONS_INDEX"
echo "Align dir:  $ALIGN_DIR"
echo "TSV dir:    $TSV_DIR"

# ---- Dependency checks ----
command -v parallel >/dev/null || { echo "parallel not found. Run: bash setup_conda_env.sh"; exit 1; }
command -v bwa      >/dev/null || { echo "bwa not found.      Run: bash setup_conda_env.sh"; exit 1; }
command -v samtools >/dev/null || { echo "samtools not found. Run: bash setup_conda_env.sh"; exit 1; }
command -v bcftools >/dev/null || { echo "bcftools not found. Run: bash setup_conda_env.sh"; exit 1; }

# ---- Aligner selection: prefer bwa-mem2 (2-3× faster on AVX2/AVX512 servers) ----
if command -v bwa-mem2 &>/dev/null; then
    BWA_CMD="bwa-mem2"
    # bwa-mem2 requires its own index format (.bwt.2bit.64); generate if absent
    if [[ ! -f "${CONS_INDEX}.bwt.2bit.64" ]]; then
        echo "[INFO] Generating bwa-mem2 index for $CONS_INDEX ..."
        bwa-mem2 index "$CONS_INDEX"
    fi
    echo "[INFO] Aligner: bwa-mem2 ($(bwa-mem2 version 2>/dev/null | tail -1))"
else
    BWA_CMD="bwa"
    echo "[INFO] Aligner: bwa mem (bwa-mem2 not found — install for speedup)"
fi

# ---- Per-sample worker ----
process_one() {
    set -euo pipefail
    fq="$1"
    b=$(basename "$fq")
    base="${b%.filtered.fq.gz}"
    base="${base%.fastq.gz}"

    log="$ALIGN_DIR/$base.bwa.log"
    mapped_txt="$ALIGN_DIR/$base.mapped.txt"
    bam="$ALIGN_DIR/$base.sorted.bam"
    tsv="$TSV_DIR/$base.tsv"

    # Resume: skip if final TSV already exists
    if [[ "${RESUME:-0}" == "1" && -f "$tsv" ]]; then
        return
    fi

    # Align → filter → sort → write BAM + index → variant call → TSV
    SAM_MAPQ_ARGS=()
    (( MAPQ_FILTER > 0 )) && SAM_MAPQ_ARGS=(-q "$MAPQ_FILTER")
    $BWA_CMD mem -t "$BWA_T" -B "$BWA_B" -O "$BWA_O" -L 5,5 -T "$BWA_MINSCORE" \
        "$CONS_INDEX" "$fq" 2>"$log" \
    | samtools view -b -@ "$SAM_T" "${SAM_MAPQ_ARGS[@]}" -F 0x904 \
    | samtools sort -@ "$SAM_T" -o "$bam"

    samtools index -@ "$SAM_T" "$bam"
    samtools view -c -F 4 "$bam" > "$mapped_txt"

    bcftools mpileup -Ou \
        -d 1000000 \
        -q 30 \
        -Q 30 \
        -C 50 \
        -a FORMAT/AD,FORMAT/DP \
        --fasta-ref "$CONS_INDEX" "$bam" \
    | bcftools query -f '%POS\t%REF\t%ALT[\t%DP\t%AD]\n' > "$tsv"
}
export -f process_one
export FASTQ_DIR CONS_INDEX ALIGN_DIR TSV_DIR BWA_T SAM_T BCF_T RESUME BWA_CMD BWA_B BWA_O BWA_MINSCORE MAPQ_FILTER

# ---- Build job list, pre-filtering already-done samples in resume mode ----
find "$FASTQ_DIR" -type f \( -name '*.filtered.fq.gz' -o -name '*.fastq.gz' \) -print0 \
| if [[ "${RESUME:-0}" == "1" ]]; then
    while IFS= read -r -d '' fq; do
        b=$(basename "$fq")
        base="${b%.filtered.fq.gz}"; base="${base%.fastq.gz}"
        [[ -f "$TSV_DIR/$base.tsv" ]] || printf '%s\0' "$fq"
    done
else
    cat
fi \
| parallel -0 --jobs "$JOBS" --halt soon,fail=1% --bar process_one {}

# ---- Build alignment summary ----
echo -e "Sample\tMapped_Reads" > "$SUMMARY_FILE"
find "$ALIGN_DIR" -name "*.mapped.txt" -print0 \
| xargs -0 awk 'FNR==1{split(FILENAME,a,"/"); n=a[length(a)]; sub(/\.mapped\.txt$/,"",n); print n"\t"$0}' \
>> "$SUMMARY_FILE"
echo "Summary written to $SUMMARY_FILE"
