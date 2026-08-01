#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# 1_filter.sh — quality-filter 5S rDNA FASTQ files with fastp.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
# =============================================================================
# 1_filter.sh — Quality-filter 5S FASTQ files with fastp.
#
# Defaults derive from this script's location (scripts/ inside the project).
# Override any path by exporting the corresponding env var before running.
#
# Resume mode: export RESUME=1 to skip samples whose .filtered.fq.gz exists.
# =============================================================================
set -euo pipefail

_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
_PIPELINE_DIR="$(dirname "$_SCRIPT_DIR")"

# ---- Paths (env-var overridable) ----
FASTQ_DIR="${FASTQ_DIR:-$_PIPELINE_DIR/results/fastq_5S}"
FILTERED_DIR="${FILTERED_DIR:-$_PIPELINE_DIR/results/filtered_fastq_5S}"
REPORT="${REPORT:-$_PIPELINE_DIR/results/filter_report.tsv}"

mkdir -p "$FILTERED_DIR"

# ---- Auto-detect cores (Linux: nproc, macOS: sysctl) ----
if command -v nproc &>/dev/null; then
    CORES="$(nproc)"
else
    CORES="$(sysctl -n hw.ncpu)"
fi
RESERVE=${RESERVE_CORES:-1}
USE=$(( CORES - RESERVE ))
FASTP_T=${FASTP_T_OVERRIDE:-2}
JOBS=${JOBS_OVERRIDE:-$(( USE / FASTP_T ))}
if (( JOBS < 1 )); then JOBS=1; fi

echo "Detected ${CORES} cores; reserving ${RESERVE}. Using ~${USE} cores:"
echo "JOBS=${JOBS}, FASTP_T=${FASTP_T}"
echo "Input dir:  $FASTQ_DIR"
echo "Output dir: $FILTERED_DIR"
echo "Report:     $REPORT"

# ---- Dependency checks ----
command -v fastp    >/dev/null || { echo "fastp not found.    Run: bash setup_conda_env.sh"; exit 1; }
command -v jq       >/dev/null || { echo "jq not found.       Run: bash setup_conda_env.sh"; exit 1; }
command -v parallel >/dev/null || { echo "parallel not found. Run: bash setup_conda_env.sh"; exit 1; }

# ---- Run fastp in parallel ----
find "$FASTQ_DIR" -type f \( -name '*.fq.gz' -o -name '*.fastq.gz' \) -print0 \
| parallel -0 --jobs "$JOBS" --halt now,fail=1 --bar '
    set -euo pipefail
    fq="{}"
    b=$(basename "$fq")
    base="${b%.fastq.gz}"
    base="${base%.fq.gz}"

    out="'"$FILTERED_DIR"'/${base}.filtered.fq.gz"
    json="'"$FILTERED_DIR"'/${base}.fastp.json"

    # Resume: skip if filtered output already exists
    if [[ "'"${RESUME:-0}"'" == "1" && -f "$out" ]]; then
        exit 0
    fi

    fastp \
      -i "$fq" \
      -o "$out" \
      -w '"$FASTP_T"' \
      -A \
      -3 \
      -W 4 \
      -M 30 \
      -l 75 \
      -n 0 \
      -e 30 \
      --json "$json" \
      --html /dev/null \
      >/dev/null
'

# ---- Build report after all jobs ----
echo -e "Sample\tBefore\tAfter\tTrimmed\tRemoved" > "$REPORT"
find "$FILTERED_DIR" -name "*.fastp.json" -print0 \
| xargs -0 jq -r '[input_filename,
    (.summary.before_filtering.total_reads|tostring),
    (.summary.after_filtering.total_reads|tostring),
    (.cutting_result.trimmed_reads|tostring),
    (.filtering_result.reads_failed_filter|tostring)]
    | join("\t")' \
| sed 's|.*/||; s|\.fastp\.json\t|\t|' \
>> "$REPORT"

echo "Filter report written to $REPORT"
