#!/usr/bin/env bash
if [[ -z "${BASH_VERSION:-}" ]]; then
  exec /usr/bin/env bash "$0" "$@"
fi
set -Eeuo pipefail
IFS=$'\n\t'
ulimit -n 16384 || true

################################################################################
# functional_annotation_v2.sh
#
# Modular RNA-seq + functional annotation + descriptive expression plotting
#
# Main outputs:
#   - SRA downloads
#   - FASTQ files
#   - fastp-trimmed FASTQ files and QC reports
#   - MultiQC report
#   - HISAT2 sorted/indexed BAM files
#   - Rsubread::featureCounts gene-level count matrix
#   - eggNOG-mapper protein annotation
#   - InterProScan protein/domain annotation
#   - DIAMOND/Swiss-Prot curated protein-similarity annotation
#   - KofamScan KEGG Orthology annotation
#   - optional dbCAN CAZyme annotation
#   - merged/descriptive RNA-seq plots and tables from one sample per organ
#
# Important:
#   This script can generate descriptive "volcano-like" fold-change plots from one
#   sample per organ, but these are NOT formal differential-expression volcano
#   plots because there are no biological replicates per organ.
#
# Usage:
#   bash functional_annotation_v2.sh [threads_per_job] [parallel_jobs]
#
# Example:
#   bash functional_annotation_v2.sh 8 4 \
#     2>&1 | tee -a functional_annotation_v2.$(date +%F).log
################################################################################

# ==============================================================================
# 0) Configuration
# ==============================================================================

_ts() {
  date +"%Y-%m-%d %H:%M:%S"
}

INFO() {
  printf '[%s] [INFO]  %s\n' "$(_ts)" "$*"
}

WARN() {
  printf '[%s] [WARN]  %s\n' "$(_ts)" "$*" >&2
}

ERROR() {
  printf '[%s] [ERROR] %s\n' "$(_ts)" "$*" >&2
  exit 1
}

log() {
  INFO "$*"
}

die() {
  ERROR "$*"
}

is_true() {
  [[ "${1^^}" == "TRUE" || "${1}" == "1" || "${1^^}" == "YES" ]]
}

check_file() {
  [[ -s "$1" ]] || die "Missing or empty file: $1"
}

check_dir() {
  [[ -d "$1" ]] || die "Missing directory: $1"
}

check_executable() {
  [[ -x "$1" ]] || die "Required executable not found or not executable: $1"
}

require_tool() {
  local label="$1"
  local tool_path="$2"
  [[ -n "${tool_path}" ]] || die "Required tool '${label}' is unset"
  if [[ "${tool_path}" == */* ]]; then
    [[ -x "${tool_path}" ]] || die "Required tool '${label}' not found or not executable at: ${tool_path}"
  else
    command -v "${tool_path}" >/dev/null 2>&1 || die "Required command '${label}' not found in PATH"
  fi
  log "Found tool: ${label} at ${tool_path}"
}

bam_index_exists() {
  local bam="$1"
  [[ -s "${bam}.bai" || -s "${bam}.csi" ]]
}

bam_is_valid() {
  local bam="$1"
  [[ -s "${bam}" ]] || return 1
  bam_index_exists "${bam}" || return 1
  [[ -x "${SAMTOOLS_BIN:-}" ]] || return 1
  "${SAMTOOLS_BIN}" quickcheck "${bam}" >/dev/null 2>&1
}

count_outputs_valid() {
  [[ -s "${COUNTS_RAW}" && -s "${COUNTS_MATRIX}" && -s "${COUNTS_RAW}.summary" ]] || return 1
  awk 'NR == 1 { exit ($1 == "gene_id" ? 0 : 1) }' "${COUNTS_MATRIX}"
}

dbcan_dir_outputs_valid() {
  local output_dir="$1"
  [[ -d "${output_dir}" ]] || return 1
  find "${output_dir}" -type f ! -name ".dbcan.done" -print -quit | grep -q .
}

dbcan_outputs_valid() {
  dbcan_dir_outputs_valid "${DBCAN_OUT_DIR}"
}

gzip_is_valid() {
  [[ -s "$1" ]] || return 1
  gzip -t "$1" >/dev/null 2>&1
}

json_is_valid() {
  [[ -s "$1" ]] || return 1
  if command -v jq >/dev/null 2>&1; then
    jq -e . "$1" >/dev/null 2>&1
  else
    awk '
      BEGIN { first = ""; last = "" }
      {
        if (first == "") {
          line = $0
          sub(/^[[:space:]]+/, "", line)
          if (line != "") first = substr(line, 1, 1)
        }
        line = $0
        sub(/[[:space:]]+$/, "", line)
        if (line != "") last = substr(line, length(line), 1)
      }
      END { exit !((first == "{" && last == "}") || (first == "[" && last == "]")) }
    ' "$1"
  fi
}

file_identity() {
  local path="$1"
  if [[ -e "${path}" ]]; then
    stat -Lc '%n|%s|%Y' "${path}"
  else
    printf '%s|missing\n' "${path}"
  fi
}

directory_identity() {
  local path="$1"
  if [[ -d "${path}" ]]; then
    stat -Lc '%n|%Y' "${path}"
  else
    printf '%s|missing\n' "${path}"
  fi
}

make_signature() {
  printf '%s\0' "$@" | sha256sum | awk '{print $1}'
}

manifest_matches() {
  local manifest="$1"
  local expected_signature="$2"
  shift 2
  [[ -s "${manifest}" ]] || return 1
  [[ "$(awk -F '\t' '$1 == "signature" { print $2; exit }' "${manifest}")" == "${expected_signature}" ]] || return 1
  local output
  for output in "$@"; do
    [[ -s "${output}" ]] || return 1
  done
}

write_manifest() {
  local manifest="$1"
  local signature="$2"
  local stage="$3"
  shift 3
  local tmp_manifest="${manifest}.tmp.$$"
  {
    printf 'signature\t%s\n' "${signature}"
    printf 'stage\t%s\n' "${stage}"
    printf 'completed_at\t%s\n' "$(_ts)"
    printf 'output\t%s\n' "$@"
  } > "${tmp_manifest}"
  mv -f "${tmp_manifest}" "${manifest}"
}

stage_manifest_path() {
  local stage="$1"
  local key="${2:-global}"
  printf '%s/%s.%s.manifest.tsv\n' "${MANIFEST_DIR}" "${stage}" "${key}"
}

mark_environment_checked() {
  ENVIRONMENT_CHECK_CACHE["$1"]=1
}

environment_was_checked() {
  [[ "${ENVIRONMENT_CHECK_CACHE[$1]:-0}" == "1" ]]
}

run_version_cmd() {
  local label="$1"
  shift
  log "Version check: ${label}"
  "$@" >&2 || WARN "Version check failed for ${label}"
}

cleanup_workdir_tmp_artifacts() {
  local stale_tmp
  local marker=".functional_annotation_v2.tmp"

  if ! is_true "${RUN_CLEANUP_TMP}"; then
    log "Final temporary-file cleanup disabled: RUN_CLEANUP_TMP=${RUN_CLEANUP_TMP}"
    return 0
  fi

  [[ -n "${SCRATCH_ROOT:-}" && -d "${SCRATCH_ROOT}" && "${SCRATCH_ROOT}" != "/" ]] \
    || die "Refusing temporary-file cleanup because SCRATCH_ROOT is invalid: ${SCRATCH_ROOT:-unset}"

  log "Cleaning script-owned temporary directories in scratch root: ${SCRATCH_ROOT}"

  while IFS= read -r -d '' stale_tmp; do
    [[ -f "${stale_tmp}/${marker}" ]] || continue
    printf '%s\n' "${stale_tmp}" >&2
    rm -rf "${stale_tmp}"
  done < <(
    find "${SCRATCH_ROOT}" \
      -mindepth 1 \
      -maxdepth 1 \
      -type d \
      -name "functional_annotation_v2.tmp.*" \
      ! -path "${tmpdir}" \
      -print0
  )
}

cleanup_current_tmpdir() {
  if [[ -z "${tmpdir:-}" || ! -d "${tmpdir}" || "${tmpdir}" == "/" ]]; then
    return 0
  fi
  if is_true "${RUN_CLEANUP_TMP}"; then
    INFO "Cleaning up tmpdir ${tmpdir}"
    rm -rf "${tmpdir}"
  else
    INFO "Preserving tmpdir because RUN_CLEANUP_TMP=${RUN_CLEANUP_TMP}: ${tmpdir}"
  fi
}

run_rscript() {
  if is_true "${USE_R_MAMBA_ENV}"; then
    check_executable "${MAMBA_BIN}"
    "${MAMBA_BIN}" run -n "${R_ENV}" Rscript "$@"
  else
    check_executable "${RSCRIPT_BIN}"
    "${RSCRIPT_BIN}" "$@"
  fi
}

download_with_curl() {
  local url="$1"
  local target="$2"

  curl \
    --fail \
    --location \
    --continue-at - \
    --retry 5 \
    --retry-all-errors \
    --connect-timeout 30 \
    --output "${target}" \
    "${url}"
}

start_time=${START_TIME:-$(date +%s)}

# ----------------------------- Runtime arguments ------------------------------

NTHREADS="${1:-${NTHREADS:-8}}"
NJOBS="${2:-${NJOBS:-4}}"          # 0 = auto-size from detected cores
PARALLEL_LOAD="${PARALLEL_LOAD:-80%}"
MEMFREE_DEFAULT="${MEMFREE_DEFAULT:-4G}"
PARALLEL_PROGRESS="${PARALLEL_PROGRESS:-FALSE}"
PARALLEL_PROGRESS_OPTS=()
if is_true "${PARALLEL_PROGRESS}" && [[ -t 2 ]]; then
  PARALLEL_PROGRESS_OPTS+=(--bar)
fi
REQUESTED_NJOBS="${NJOBS}"

[[ "${NTHREADS}" =~ ^[0-9]+$ ]] || { printf "NTHREADS must be an integer: %s\n" "${NTHREADS}" >&2; exit 1; }
[[ "${NJOBS}" =~ ^[0-9]+$ ]] || { printf "NJOBS must be an integer: %s\n" "${NJOBS}" >&2; exit 1; }
(( NTHREADS >= 1 )) || { printf "NTHREADS must be >= 1\n" >&2; exit 1; }

NPROC_VAL="$(nproc 2>/dev/null || printf "1")"
if (( NJOBS == 0 )); then
  NJOBS=$((NPROC_VAL / NTHREADS))
  (( NJOBS < 1 )) && NJOBS=1
fi
CPU_JOB_CAP=$((NPROC_VAL / NTHREADS))
(( CPU_JOB_CAP < 1 )) && CPU_JOB_CAP=1
if (( NJOBS > CPU_JOB_CAP )); then
  WARN "Reducing parallel jobs from ${NJOBS} to ${CPU_JOB_CAP} to keep NTHREADS * NJOBS within ${NPROC_VAL} detected CPUs"
  NJOBS="${CPU_JOB_CAP}"
fi
TOTAL_CPUS=$((NTHREADS * NJOBS))
(( TOTAL_CPUS > NPROC_VAL )) && TOTAL_CPUS="${NPROC_VAL}"
(( TOTAL_CPUS < 1 )) && TOTAL_CPUS=1

# ----------------------------- Stage switches ---------------------------------

RUN_MODE="${RUN_MODE:-annotation_only}"  # all, annotation_only, plots_only, preflight, custom
PREFLIGHT_ONLY="FALSE"

STAGE_SWITCHES=(
  RUN_CREATE_METADATA RUN_DOWNLOAD RUN_FASTQ_DUMP RUN_COMPRESS_FASTQ
  RUN_FASTP RUN_MULTIQC_AFTER_FASTP RUN_HISAT2_INDEX RUN_HISAT2_MAPPING
  RUN_FEATURECOUNTS RUN_EGGNOG_MAPPER RUN_INTERPROSCAN
  RUN_DIAMOND_ANNOTATION RUN_KOFAMSCAN RUN_DBCAN RUN_CREATE_R_SCRIPT
  RUN_RNASEQ_PLOTS_TABLES RUN_MULTIQC_FINAL
  RUN_CREATE_GENOME_GRAPHICS_R_SCRIPT RUN_GENOME_ANNOTATION_GRAPHICS
)

for stage_switch in "${STAGE_SWITCHES[@]}"; do
  printf -v "DEFAULT_${stage_switch}" '%s' "FALSE"
done

enable_default_stages() {
  local stage_switch
  for stage_switch in "$@"; do
    printf -v "DEFAULT_${stage_switch}" '%s' "TRUE"
  done
}

case "${RUN_MODE}" in
  all)
    enable_default_stages "${STAGE_SWITCHES[@]}"
    ;;
  annotation_only)
    enable_default_stages RUN_EGGNOG_MAPPER RUN_INTERPROSCAN RUN_DIAMOND_ANNOTATION RUN_KOFAMSCAN RUN_DBCAN
    ;;
  plots_only)
    enable_default_stages RUN_CREATE_R_SCRIPT RUN_RNASEQ_PLOTS_TABLES RUN_MULTIQC_FINAL \
      RUN_CREATE_GENOME_GRAPHICS_R_SCRIPT RUN_GENOME_ANNOTATION_GRAPHICS
    ;;
  preflight)
    enable_default_stages "${STAGE_SWITCHES[@]}"
    PREFLIGHT_ONLY="TRUE"
    ;;
  startup_check|custom)
    ;;
  *)
    die "Invalid RUN_MODE='${RUN_MODE}'. Use all, annotation_only, plots_only, startup_check, preflight, or custom."
    ;;
esac

for stage_switch in "${STAGE_SWITCHES[@]}"; do
  default_switch="DEFAULT_${stage_switch}"
  if [[ -z "${!stage_switch:-}" ]]; then
    printf -v "${stage_switch}" '%s' "${!default_switch}"
  fi
done

# Downloads and database replacement stay opt-in for every profile.
EGGNOG_DOWNLOAD_DATA="${EGGNOG_DOWNLOAD_DATA:-FALSE}"
EGGNOG_FORCE_DOWNLOAD="${EGGNOG_FORCE_DOWNLOAD:-FALSE}"
RUN_INTERPROSCAN_DB_DOWNLOAD="${RUN_INTERPROSCAN_DB_DOWNLOAD:-FALSE}"
INTERPROSCAN_FORCE_DB_UPDATE="${INTERPROSCAN_FORCE_DB_UPDATE:-FALSE}"
RUN_SWISSPROT_DB_DOWNLOAD="${RUN_SWISSPROT_DB_DOWNLOAD:-FALSE}"
SWISSPROT_FORCE_DB_UPDATE="${SWISSPROT_FORCE_DB_UPDATE:-FALSE}"
RUN_DIAMOND_DB_PREP="${RUN_DIAMOND_DB_PREP:-FALSE}"
DIAMOND_FORCE_DB_REBUILD="${DIAMOND_FORCE_DB_REBUILD:-FALSE}"
RUN_KOFAM_DB_DOWNLOAD="${RUN_KOFAM_DB_DOWNLOAD:-FALSE}"
RUN_KOFAM_DB_PREP="${RUN_KOFAM_DB_PREP:-FALSE}"
KOFAM_FORCE_DB_UPDATE="${KOFAM_FORCE_DB_UPDATE:-FALSE}"
RUN_DBCAN_DB_PREP="${RUN_DBCAN_DB_PREP:-FALSE}"
DBCAN_FORCE_DB_UPDATE="${DBCAN_FORCE_DB_UPDATE:-FALSE}"

RUN_TIDYPLOTS="${RUN_TIDYPLOTS:-TRUE}"
RUN_CLEANUP_TMP="${RUN_CLEANUP_TMP:-TRUE}"

# Formal volcano plots should remain FALSE with one sample per organ.
# The R script will instead generate descriptive volcano-like fold-change plots.
RUN_FORMAL_DESEQ2_VOLCANO="${RUN_FORMAL_DESEQ2_VOLCANO:-FALSE}"
RUN_DESCRIPTIVE_VOLCANO="${RUN_DESCRIPTIVE_VOLCANO:-TRUE}"

# ----------------------------- Project directories ----------------------------

PROJECT_DIR="${PROJECT_DIR:-${WORKDIR:-${PWD}}}"
workdir="${PROJECT_DIR}"
[[ -n "${workdir}" && "${workdir}" != "/" ]] || die "Refusing invalid PROJECT_DIR: ${workdir:-unset}"
if ! is_true "${PREFLIGHT_ONLY}"; then
  mkdir -p "${workdir}" || die "Could not create PROJECT_DIR: ${workdir}"
fi
if is_true "${PREFLIGHT_ONLY}"; then
  DEFAULT_SCRATCH_ROOT="${TMPDIR:-/tmp}"
else
  DEFAULT_SCRATCH_ROOT="${TMPDIR:-${workdir}}"
fi
SCRATCH_ROOT="${SCRATCH_ROOT:-${DEFAULT_SCRATCH_ROOT}}"
[[ -n "${SCRATCH_ROOT}" && "${SCRATCH_ROOT}" != "/" ]] || die "Refusing invalid SCRATCH_ROOT: ${SCRATCH_ROOT:-unset}"
mkdir -p "${SCRATCH_ROOT}" || die "Could not create SCRATCH_ROOT: ${SCRATCH_ROOT}"
tmpdir=$(mktemp -d -p "${SCRATCH_ROOT}" functional_annotation_v2.tmp.XXXXXXXXXX)
touch "${tmpdir}/.functional_annotation_v2.tmp"
trap cleanup_current_tmpdir EXIT

DIR_00_LOGS="${DIR_00_LOGS:-${PROJECT_DIR}/00_logs}"
DIR_01_METADATA="${DIR_01_METADATA:-${PROJECT_DIR}/01_metadata}"
DIR_02_SRA="${DIR_02_SRA:-${PROJECT_DIR}/02_sra}"
DIR_03_FASTQ="${DIR_03_FASTQ:-${PROJECT_DIR}/03_fastq}"

DIR_04_FASTP="${DIR_04_FASTP:-${PROJECT_DIR}/04_fastp}"
DIR_04_TRIMMED="${DIR_04_TRIMMED:-${DIR_04_FASTP}/trimmed_fastq}"
DIR_04_FASTP_JSON="${DIR_04_FASTP_JSON:-${DIR_04_FASTP}/json}"
DIR_04_FASTP_HTML="${DIR_04_FASTP_HTML:-${DIR_04_FASTP}/html}"
DIR_04_FASTP_FAILED="${DIR_04_FASTP_FAILED:-${DIR_04_FASTP}/failed_reads}"
DIR_04_FASTP_UNPAIRED="${DIR_04_FASTP_UNPAIRED:-${DIR_04_FASTP}/unpaired_reads}"

DIR_05_MULTIQC="${DIR_05_MULTIQC:-${PROJECT_DIR}/05_multiqc}"

DIR_06_HISAT2="${DIR_06_HISAT2:-${PROJECT_DIR}/06_hisat2}"
DIR_06_INDEX="${DIR_06_INDEX:-${DIR_06_HISAT2}/index}"
DIR_06_BAM="${DIR_06_BAM:-${DIR_06_HISAT2}/bam}"
DIR_06_SUMMARY="${DIR_06_SUMMARY:-${DIR_06_HISAT2}/summary}"

DIR_07_COUNTS="${DIR_07_COUNTS:-${PROJECT_DIR}/07_featurecounts}"
DIR_08_EGGNOG="${DIR_08_EGGNOG:-${PROJECT_DIR}/08_eggnog_mapper}"
DIR_09_INTERPRO="${DIR_09_INTERPRO:-${PROJECT_DIR}/09_interproscan}"
DIR_10_DIAMOND="${DIR_10_DIAMOND:-${PROJECT_DIR}/10_diamond_annotation}"
DIR_11_KOFAM="${DIR_11_KOFAM:-${PROJECT_DIR}/11_kofamscan}"
DIR_12_DBCAN="${DIR_12_DBCAN:-${PROJECT_DIR}/12_dbcan}"
DIR_13_CONSENSUS="${DIR_13_CONSENSUS:-${PROJECT_DIR}/13_consensus_annotation}"
DIR_13_SCRIPTS="${DIR_13_SCRIPTS:-${DIR_13_CONSENSUS}/scripts}"
DIR_13_TABLES="${DIR_13_TABLES:-${DIR_13_CONSENSUS}/tables}"
DIR_13_PLOTS="${DIR_13_PLOTS:-${DIR_13_CONSENSUS}/plots}"
DIR_13_ANNOTATION="${DIR_13_ANNOTATION:-${DIR_13_CONSENSUS}/annotation}"
MANIFEST_DIR="${MANIFEST_DIR:-${DIR_00_LOGS}/manifests}"

if ! is_true "${PREFLIGHT_ONLY}"; then
  mkdir -p \
    "${DIR_00_LOGS}" \
  "${DIR_01_METADATA}" \
  "${DIR_02_SRA}" \
  "${DIR_03_FASTQ}" \
  "${DIR_04_FASTP}" \
  "${DIR_04_TRIMMED}" \
  "${DIR_04_FASTP_JSON}" \
  "${DIR_04_FASTP_HTML}" \
  "${DIR_04_FASTP_FAILED}" \
  "${DIR_04_FASTP_UNPAIRED}" \
  "${DIR_05_MULTIQC}" \
  "${DIR_06_HISAT2}" \
  "${DIR_06_INDEX}" \
  "${DIR_06_BAM}" \
  "${DIR_06_SUMMARY}" \
  "${DIR_07_COUNTS}" \
  "${DIR_08_EGGNOG}" \
  "${DIR_09_INTERPRO}" \
  "${DIR_10_DIAMOND}" \
  "${DIR_11_KOFAM}" \
  "${DIR_12_DBCAN}" \
  "${DIR_13_CONSENSUS}" \
  "${DIR_13_SCRIPTS}" \
  "${DIR_13_TABLES}" \
  "${DIR_13_PLOTS}" \
    "${DIR_13_ANNOTATION}" \
    "${MANIFEST_DIR}"
fi

declare -A ENVIRONMENT_CHECK_CACHE=()

LOG_TO_FILE="${LOG_TO_FILE:-$([[ "${PREFLIGHT_ONLY}" == "TRUE" ]] && printf '0' || printf '1')}"
LOG_FILE="${LOG_FILE:-${DIR_00_LOGS}/runtime_$(date +%Y%m%d_%H%M%S).log}"
if [[ "${LOG_TO_FILE}" == "1" ]]; then
  exec > >(tee -ia "${LOG_FILE}") 2> >(tee -ia "${LOG_FILE}" >&2)
  printf "[%s] [INFO]  Logging to %s\n" "$(date '+%F %T')" "${LOG_FILE}"
fi

# ----------------------------- SRA metadata -----------------------------------

SRA_METADATA="${DIR_01_METADATA}/sra_metadata.tsv"
SAMPLE_METADATA="${DIR_01_METADATA}/sample_metadata.tsv"

if ! is_true "${PREFLIGHT_ONLY}" && [[ ! -s "${SRA_METADATA}" || "${RUN_CREATE_METADATA}" == "TRUE" ]]; then
  cat > "${SRA_METADATA}" <<'EOF'
accession	organ
SRR29142956	Testes
SRR29142957	Muscle
SRR29142958	Liver
SRR29142959	Gills
EOF
fi

if ! is_true "${PREFLIGHT_ONLY}" && [[ ! -s "${SAMPLE_METADATA}" || "${RUN_CREATE_METADATA}" == "TRUE" ]]; then
  awk '
    BEGIN { FS = OFS = "\t" }
    NR == 1 { print "sample", "accession", "organ"; next }
    NR > 1 { print $1 "." $2, $1, $2 }
  ' "${SRA_METADATA}" > "${SAMPLE_METADATA}"
fi

SRA_RUNS=()
if [[ -s "${SRA_METADATA}" ]]; then
  mapfile -t SRA_RUNS < <(awk 'NR > 1 {print $1}' "${SRA_METADATA}")
elif ! is_true "${RUN_CREATE_METADATA}" && ( is_true "${RUN_DOWNLOAD}" || is_true "${RUN_FASTQ_DUMP}" || is_true "${RUN_COMPRESS_FASTQ}" \
  || is_true "${RUN_FASTP}" || is_true "${RUN_HISAT2_MAPPING}" ); then
  die "Missing SRA metadata required by enabled sample stages: ${SRA_METADATA}"
fi

# ----------------------------- EGAPx outputs ----------------------------------

EGAPX_OUT_DIR="${EGAPX_OUT_DIR:-/data/bferrette/chondrichthyes/batoids/narcine/assembly/brasiliensis/hic_assembly/egapx/out}"

EGAPX_GENOME="${EGAPX_GENOME:-${EGAPX_OUT_DIR}/complete.genomic.fna}"
EGAPX_GFF="${EGAPX_GFF:-${EGAPX_OUT_DIR}/complete.genomic.gff}"
EGAPX_GTF="${EGAPX_GTF:-${EGAPX_OUT_DIR}/complete.genomic.gtf}"
EGAPX_CDS="${EGAPX_CDS:-${EGAPX_OUT_DIR}/complete.cds.fna}"
EGAPX_TRANSCRIPTS="${EGAPX_TRANSCRIPTS:-${EGAPX_OUT_DIR}/complete.transcripts.fna}"
EGAPX_PROTEINS="${EGAPX_PROTEINS:-${EGAPX_OUT_DIR}/complete.proteins.faa}"

# ----------------------------- Tool paths -------------------------------------

MINIFORGE_ROOT="${MINIFORGE_ROOT:-/data/bferrette/miniforge3}"

HISAT2_BIN="${HISAT2_BIN:-${HISAT2:-/home/bferrette/software/hisat2/hisat2}}"
HISAT2_BUILD_BIN="${HISAT2_BUILD_BIN:-${HISAT2_BUILD:-/home/bferrette/software/hisat2/hisat2-build}}"
SAMTOOLS_BIN="${SAMTOOLS_BIN:-${SAMTOOLS:-/home/bferrette/software/samtools-1.23.1/samtools}}"

PREFETCH_BIN="${PREFETCH_BIN:-${PREFETCH:-/home/bferrette/software/sratoolkit.3.4.1-ubuntu64/bin/prefetch}}"
FASTERQ_DUMP_BIN="${FASTERQ_DUMP_BIN:-${FASTERQ_DUMP:-/home/bferrette/software/sratoolkit.3.4.1-ubuntu64/bin/fasterq-dump}}"

FASTP_BIN="${FASTP_BIN:-${FASTP:-/home/bferrette/software/fastp/fastp}}"

PIGZ_BIN="${PIGZ_BIN:-${PIGZ:-$(command -v pigz || printf '%s/bin/pigz' "${MINIFORGE_ROOT}")}}"
PARALLEL_BIN="${PARALLEL_BIN:-${PARALLEL_CMD:-$(command -v parallel || printf '%s/bin/parallel' "${MINIFORGE_ROOT}")}}"
MAMBA_BIN="${MAMBA_BIN:-${MAMBA:-$(command -v mamba || printf '%s/bin/mamba' "${MINIFORGE_ROOT}")}}"
MULTIQC_BIN="${MULTIQC_BIN:-${MULTIQC:-/data/bferrette/miniforge3/bin/multiqc}}"
RSCRIPT_BIN="${RSCRIPT_BIN:-${RSCRIPT:-${MINIFORGE_ROOT}/bin/Rscript}}"

# Several command-name environment variables are parsed by the tools themselves
# or inherited by GNU Parallel children. Keep executable paths in *_BIN variables.
unset FASTP FASTERQ_DUMP HISAT2 HISAT2_BUILD MAMBA MULTIQC PARALLEL PARALLEL_CMD PIGZ PREFETCH RSCRIPT SAMTOOLS
export SHELL=/bin/bash

# R plotting/summary environment. Recommended: use the dedicated mamba env.
USE_R_MAMBA_ENV="${USE_R_MAMBA_ENV:-TRUE}"
R_ENV="${R_ENV:-Renv}"
R_REQUIRED_PACKAGES=(
  tidyverse
  optparse
  data.table
  gt
  gtsummary
  scales
  ggrepel
  svglite
  ragg
  Cairo
  Rsubread
  RColorBrewer
  GenomicRanges
  Rsamtools
  Biostrings
  karyoploteR
)

SINGULARITY_BIN="${SINGULARITY_BIN:-${SINGULARITY:-$(command -v singularity || command -v apptainer || true)}}"
unset SINGULARITY

# ----------------------------- SRA / FASTQ options ----------------------------

PREFETCH_MAX_SIZE="${PREFETCH_MAX_SIZE:-500G}"
FASTERQ_MEM="${FASTERQ_MEM:-16G}"

# ----------------------------- fastp options ----------------------------------

FASTP_QUAL="${FASTP_QUAL:-20}"
FASTP_LENGTH_REQUIRED="${FASTP_LENGTH_REQUIRED:-36}"
FASTP_N_BASE_LIMIT="${FASTP_N_BASE_LIMIT:-5}"
FASTP_UNQUALIFIED_PERCENT_LIMIT="${FASTP_UNQUALIFIED_PERCENT_LIMIT:-40}"
FASTP_COMPLEXITY_FILTER="${FASTP_COMPLEXITY_FILTER:-TRUE}"
FASTP_COMPLEXITY_THRESHOLD="${FASTP_COMPLEXITY_THRESHOLD:-30}"
FASTP_OVERREPRESENTATION="${FASTP_OVERREPRESENTATION:-TRUE}"
FASTP_OVERREPRESENTATION_SAMPLING="${FASTP_OVERREPRESENTATION_SAMPLING:-20}"
FASTP_POLYG="${FASTP_POLYG:-TRUE}"
FASTP_POLYG_MIN_LEN="${FASTP_POLYG_MIN_LEN:-10}"
FASTP_POLYX="${FASTP_POLYX:-TRUE}"
FASTP_POLYX_MIN_LEN="${FASTP_POLYX_MIN_LEN:-10}"
FASTP_CORRECTION="${FASTP_CORRECTION:-TRUE}"
FASTP_COMPRESSION="${FASTP_COMPRESSION:-6}"
FASTP_SAVE_FAILED="${FASTP_SAVE_FAILED:-TRUE}"
FASTP_SAVE_UNPAIRED="${FASTP_SAVE_UNPAIRED:-TRUE}"
FASTP_DONT_OVERWRITE="${FASTP_DONT_OVERWRITE:-FALSE}"

# ----------------------------- HISAT2 options ---------------------------------

HISAT2_INDEX_PREFIX="${HISAT2_INDEX_PREFIX:-${DIR_06_INDEX}/egapx_reference_hisat2}"

RNA_STRANDNESS="${RNA_STRANDNESS:-}"       # "", RF, or FR
KNOWN_SPLICESITES="${KNOWN_SPLICESITES:-}" # optional; generated externally if desired

HISAT2_PRESET="${HISAT2_PRESET:---dta}"
HISAT2_EXTRA_OPTS="${HISAT2_EXTRA_OPTS:---new-summary --time --mm}"
MAPPING_INPUT_SOURCE="${MAPPING_INPUT_SOURCE:-auto}" # auto, trimmed, raw

SAMTOOLS_SORT_MEM="${SAMTOOLS_SORT_MEM:-4G}"
if (( NTHREADS == 1 )); then
  HISAT2_THREADS="${HISAT2_THREADS:-1}"
  SAMTOOLS_SORT_THREADS="${SAMTOOLS_SORT_THREADS:-1}"
else
  HISAT2_THREADS="${HISAT2_THREADS:-$((NTHREADS * 3 / 4))}"
  (( HISAT2_THREADS < 1 )) && HISAT2_THREADS=1
  SAMTOOLS_SORT_THREADS="${SAMTOOLS_SORT_THREADS:-$((NTHREADS - HISAT2_THREADS))}"
  (( SAMTOOLS_SORT_THREADS < 1 )) && SAMTOOLS_SORT_THREADS=1
fi
SAMTOOLS_AUX_THREADS="${SAMTOOLS_AUX_THREADS:-$((SAMTOOLS_SORT_THREADS < 2 ? SAMTOOLS_SORT_THREADS : 2))}"
for resource_value in HISAT2_THREADS SAMTOOLS_SORT_THREADS SAMTOOLS_AUX_THREADS; do
  if [[ ! "${!resource_value}" =~ ^[0-9]+$ ]] || (( ${!resource_value} < 1 )); then
    die "${resource_value} must be a positive integer: ${!resource_value}"
  fi
done
MAPPING_CPU_PER_JOB=$((HISAT2_THREADS + SAMTOOLS_SORT_THREADS))
MAPPING_CPU_JOB_CAP=$((NPROC_VAL / MAPPING_CPU_PER_JOB))
(( MAPPING_CPU_JOB_CAP < 1 )) && MAPPING_CPU_JOB_CAP=1
MAPPING_JOBS="${MAPPING_JOBS:-${NJOBS}}"
if [[ ! "${MAPPING_JOBS}" =~ ^[0-9]+$ ]] || (( MAPPING_JOBS < 1 )); then
  die "MAPPING_JOBS must be a positive integer: ${MAPPING_JOBS}"
fi
if (( MAPPING_JOBS > MAPPING_CPU_JOB_CAP )); then
  WARN "Reducing mapping jobs from ${MAPPING_JOBS} to ${MAPPING_CPU_JOB_CAP}; each mapping pipeline can use ${MAPPING_CPU_PER_JOB} concurrent CPU threads"
  MAPPING_JOBS="${MAPPING_CPU_JOB_CAP}"
fi
MEM_AVAILABLE_BYTES="$(awk '/^MemAvailable:/ { printf "%.0f", $2 * 1024; exit }' /proc/meminfo 2>/dev/null || printf '0')"
SAMTOOLS_SORT_MEM_BYTES="$(numfmt --from=iec "${SAMTOOLS_SORT_MEM^^}" 2>/dev/null || printf '0')"
MAPPING_FIXED_MEMORY_BYTES="$((4 * 1024 * 1024 * 1024))"
if (( MEM_AVAILABLE_BYTES > 0 && SAMTOOLS_SORT_MEM_BYTES > 0 )); then
  MAPPING_ESTIMATED_MEMORY_BYTES=$((SAMTOOLS_SORT_MEM_BYTES * SAMTOOLS_SORT_THREADS + MAPPING_FIXED_MEMORY_BYTES))
  MAPPING_MEMORY_JOB_CAP=$((MEM_AVAILABLE_BYTES * 85 / 100 / MAPPING_ESTIMATED_MEMORY_BYTES))
  (( MAPPING_MEMORY_JOB_CAP < 1 )) && MAPPING_MEMORY_JOB_CAP=1
  if (( MAPPING_JOBS > MAPPING_MEMORY_JOB_CAP )); then
    WARN "Reducing mapping jobs from ${MAPPING_JOBS} to ${MAPPING_MEMORY_JOB_CAP} based on samtools sort memory and current MemAvailable"
    MAPPING_JOBS="${MAPPING_MEMORY_JOB_CAP}"
  fi
fi

RG_PLATFORM="${RG_PLATFORM:-ILLUMINA}"
RG_CENTER="${RG_CENTER:-NA}"

# ----------------------------- featureCounts options --------------------------

FEATURECOUNTS_FEATURE_TYPE="${FEATURECOUNTS_FEATURE_TYPE:-exon}"
FEATURECOUNTS_GROUP_ATTRIBUTE="${FEATURECOUNTS_GROUP_ATTRIBUTE:-gene_id}"
FEATURECOUNTS_STRANDNESS="${FEATURECOUNTS_STRANDNESS:-0}" # 0 unstranded, 1 stranded, 2 reverse stranded
FEATURECOUNTS_IS_PAIRED="${FEATURECOUNTS_IS_PAIRED:-TRUE}"
FEATURECOUNTS_COUNT_MULTIMAPPERS="${FEATURECOUNTS_COUNT_MULTIMAPPERS:-FALSE}"
FEATURECOUNTS_ALLOW_MULTI_OVERLAP="${FEATURECOUNTS_ALLOW_MULTI_OVERLAP:-FALSE}"

COUNTS_RAW="${DIR_07_COUNTS}/gene_counts.Rsubread.featureCounts.tsv"
COUNTS_MATRIX="${DIR_07_COUNTS}/gene_count_matrix.tsv"
FEATURECOUNTS_R_SCRIPT="${DIR_13_SCRIPTS}/run_featurecounts_rsubread.R"

if is_true "${RUN_RNASEQ_PLOTS_TABLES}" && [[ ! -s "${COUNTS_MATRIX}" ]] && ! is_true "${RUN_FEATURECOUNTS}"; then
  WARN "RUN_RNASEQ_PLOTS_TABLES=TRUE but count matrix is missing: ${COUNTS_MATRIX}"
  WARN "Provide an existing count matrix or set RUN_FEATURECOUNTS=TRUE with valid BAM inputs."
fi

# ----------------------------- eggNOG-mapper options --------------------------

EGGNOG_ENV="${EGGNOG_ENV:-EGGNOGMAPPERenv}"
EGGNOG_DATA_DIR="${EGGNOG_DATA_DIR:-/data/bferrette/databases/eggnog_mapper}"

EGGNOG_PREFIX="${EGGNOG_PREFIX:-narcine_brasiliensis_eggnog}"
EGGNOG_TMP_DIR="${tmpdir}/eggnog_tmp"
EGGNOG_SCRATCH_DIR="${tmpdir}/eggnog_scratch"

EGGNOG_METHOD="${EGGNOG_METHOD:-diamond}"
EGGNOG_ITYPE="${EGGNOG_ITYPE:-proteins}"
EGGNOG_TAX_SCOPE="${EGGNOG_TAX_SCOPE:-eukaryota}"
EGGNOG_GO_EVIDENCE="${EGGNOG_GO_EVIDENCE:-non-electronic}"
EGGNOG_SENSMODE="${EGGNOG_SENSMODE:-ultra-sensitive}"
EGGNOG_PFAM_REALIGN="${EGGNOG_PFAM_REALIGN:-none}"

EGGNOG_OVERRIDE="${EGGNOG_OVERRIDE:-FALSE}"
EGGNOG_EXCEL="${EGGNOG_EXCEL:-TRUE}"

mkdir -p "${EGGNOG_TMP_DIR}" "${EGGNOG_SCRATCH_DIR}"

# ----------------------------- InterProScan options ---------------------------

INTERPROSCAN_SIF="${INTERPROSCAN_SIF:-/data/bferrette/containers/interproscan.sif}"
INTERPROSCAN_DATA_DIR="${INTERPROSCAN_DATA_DIR:-/data/bferrette/databases/interproscan-5.78-109.0/data}"
INTERPROSCAN_DB_DIR="${INTERPROSCAN_DB_DIR:-/data/bferrette/databases}"
INTERPROSCAN_DATA_VERSION="${INTERPROSCAN_DATA_VERSION:-5.78-109.0}"
INTERPROSCAN_DATA_URL="${INTERPROSCAN_DATA_URL:-http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/${INTERPROSCAN_DATA_VERSION}/alt/interproscan-data-${INTERPROSCAN_DATA_VERSION}.tar.gz}"
INTERPROSCAN_DATA_MD5_URL="${INTERPROSCAN_DATA_MD5_URL:-${INTERPROSCAN_DATA_URL}.md5}"
INTERPROSCAN_DATA_ARCHIVE="${INTERPROSCAN_DATA_ARCHIVE:-${INTERPROSCAN_DB_DIR}/interproscan-data-${INTERPROSCAN_DATA_VERSION}.tar.gz}"
INTERPROSCAN_DATA_MD5="${INTERPROSCAN_DATA_MD5:-${INTERPROSCAN_DATA_ARCHIVE}.md5}"

INTERPROSCAN_INPUT="${INTERPROSCAN_INPUT:-${EGAPX_PROTEINS}}"
INTERPROSCAN_PREFIX="${INTERPROSCAN_PREFIX:-narcine_brasiliensis_interproscan}"
INTERPROSCAN_TMP_DIR="${tmpdir}/interproscan_tmp"

INTERPROSCAN_CONTAINER_BIN="${INTERPROSCAN_CONTAINER_BIN:-/opt/interproscan/interproscan.sh}"
INTERPROSCAN_CONTAINER_DATA="${INTERPROSCAN_CONTAINER_DATA:-/opt/interproscan/data}"
INTERPROSCAN_CONTAINER_INPUT="/input/$(basename "${INTERPROSCAN_INPUT}")"
INTERPROSCAN_CONTAINER_OUTPUT="${INTERPROSCAN_CONTAINER_OUTPUT:-/output}"
INTERPROSCAN_CONTAINER_TEMP="${INTERPROSCAN_CONTAINER_TEMP:-/temp}"

INTERPROSCAN_FORMATS="${INTERPROSCAN_FORMATS:-TSV,GFF3,XML}"
INTERPROSCAN_SEQTYPE="${INTERPROSCAN_SEQTYPE:-p}"
INTERPROSCAN_APPLICATIONS="${INTERPROSCAN_APPLICATIONS:-}"

INTERPROSCAN_DISABLE_PRECALC="${INTERPROSCAN_DISABLE_PRECALC:-TRUE}"
INTERPROSCAN_GO_TERMS="${INTERPROSCAN_GO_TERMS:-TRUE}"
INTERPROSCAN_PATHWAYS="${INTERPROSCAN_PATHWAYS:-TRUE}"
INTERPROSCAN_IPRLOOKUP="${INTERPROSCAN_IPRLOOKUP:-TRUE}"
INTERPROSCAN_HIGHMEM="${INTERPROSCAN_HIGHMEM:-FALSE}"
INTERPROSCAN_VERBOSE="${INTERPROSCAN_VERBOSE:-FALSE}"
INTERPROSCAN_OVERRIDE="${INTERPROSCAN_OVERRIDE:-FALSE}"

mkdir -p "${INTERPROSCAN_TMP_DIR}"

# ----------------------------- DIAMOND / Swiss-Prot options --------------------

DIAMOND_ENV="${DIAMOND_ENV:-DIAMONDenv}"

DIAMOND_DB_DIR="${DIAMOND_DB_DIR:-/data/bferrette/databases/swiss_prot}"
SWISSPROT_URL="${SWISSPROT_URL:-https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz}"
SWISSPROT_FASTA_GZ="${SWISSPROT_FASTA_GZ:-${DIAMOND_DB_DIR}/uniprot_sprot.fasta.gz}"
SWISSPROT_FASTA="${SWISSPROT_FASTA:-${DIAMOND_DB_DIR}/uniprot_sprot.fasta}"
SWISSPROT_MANIFEST="${SWISSPROT_MANIFEST:-${DIAMOND_DB_DIR}/uniprot_sprot.manifest.tsv}"
DIAMOND_DB_FASTA="${DIAMOND_DB_FASTA:-${SWISSPROT_FASTA}}"
DIAMOND_DB_PREFIX="${DIAMOND_DB_PREFIX:-${DIAMOND_DB_DIR}/uniprot_sprot}"
DIAMOND_DB_FILE="${DIAMOND_DB_FILE:-${DIAMOND_DB_PREFIX}.dmnd}"

DIAMOND_PREFIX="${DIAMOND_PREFIX:-narcine_brasiliensis_swissprot}"
DIAMOND_OUT="${DIAMOND_OUT:-${DIR_10_DIAMOND}/${DIAMOND_PREFIX}.diamond.tsv}"

DIAMOND_EVALUE="${DIAMOND_EVALUE:-1e-5}"
DIAMOND_MAX_TARGET_SEQS="${DIAMOND_MAX_TARGET_SEQS:-5}"
DIAMOND_MIN_IDENTITY="${DIAMOND_MIN_IDENTITY:-30}"
DIAMOND_QUERY_COVER="${DIAMOND_QUERY_COVER:-50}"
DIAMOND_SUBJECT_COVER="${DIAMOND_SUBJECT_COVER:-50}"
DIAMOND_SENSITIVITY="${DIAMOND_SENSITIVITY:-ultra-sensitive}"
DIAMOND_COMP_BASED_STATS="${DIAMOND_COMP_BASED_STATS:-3}"
DIAMOND_GAPPED_FILTER_EVALUE="${DIAMOND_GAPPED_FILTER_EVALUE:-1.0}"
DIAMOND_EXT_CHUNK_SIZE="${DIAMOND_EXT_CHUNK_SIZE:-400}"
DIAMOND_BIN="${DIAMOND_BIN:-64}"
DIAMOND_EXT="${DIAMOND_EXT:-}"
DIAMOND_BLOCK_SIZE="${DIAMOND_BLOCK_SIZE:-4}"
DIAMOND_INDEX_CHUNKS="${DIAMOND_INDEX_CHUNKS:-1}"
DIAMOND_TMP_DIR="${DIAMOND_TMP_DIR:-${tmpdir}/diamond_tmp}"

mkdir -p "${DIAMOND_TMP_DIR}"

# ----------------------------- KofamScan options -------------------------------

KOFAM_ENV="${KOFAM_ENV:-KOFAMSCANenv}"

KOFAM_DB_DIR="${KOFAM_DB_DIR:-/data/bferrette/databases/kofam}"
KOFAM_KO_LIST_URL="${KOFAM_KO_LIST_URL:-https://www.genome.jp/ftp/db/kofam/ko_list.gz}"
KOFAM_PROFILE_ARCHIVE_URL="${KOFAM_PROFILE_ARCHIVE_URL:-https://www.genome.jp/ftp/db/kofam/profiles.tar.gz}"
KOFAM_PROFILE_ARCHIVE="${KOFAM_PROFILE_ARCHIVE:-${KOFAM_DB_DIR}/profiles.tar.gz}"
KOFAM_KO_LIST_GZ="${KOFAM_KO_LIST_GZ:-${KOFAM_DB_DIR}/ko_list.gz}"

KOFAM_PROFILE_DIR="${KOFAM_PROFILE_DIR:-${KOFAM_DB_DIR}/profiles}"
KOFAM_KO_LIST="${KOFAM_KO_LIST:-${KOFAM_DB_DIR}/ko_list}"

KOFAM_PREFIX="${KOFAM_PREFIX:-narcine_brasiliensis_kofamscan}"
KOFAM_OUT_DETAIL="${KOFAM_OUT_DETAIL:-${DIR_11_KOFAM}/${KOFAM_PREFIX}.detail.tsv}"
KOFAM_OUT_MAPPER="${KOFAM_OUT_MAPPER:-${DIR_11_KOFAM}/${KOFAM_PREFIX}.mapper.tsv}"

KOFAM_TMP_DIR="${tmpdir}/kofamscan_tmp"
KOFAM_REANNOTATION_CACHE_DIR="${KOFAM_REANNOTATION_CACHE_DIR:-${DIR_11_KOFAM}/.kofamscan_reannotation}"

mkdir -p "${KOFAM_TMP_DIR}"

# ----------------------------- dbCAN options ----------------------------------

DBCAN_ENV="${DBCAN_ENV:-DBCANenv}"

DBCAN_DB_DIR="${DBCAN_DB_DIR:-/data/bferrette/databases/dbcan}"
DBCAN_PREFIX="${DBCAN_PREFIX:-narcine_brasiliensis_dbcan}"
DBCAN_OUT_DIR="${DBCAN_OUT_DIR:-${DIR_12_DBCAN}/${DBCAN_PREFIX}}"

DBCAN_MODE="${DBCAN_MODE:-protein}"
DBCAN_METHODS="${DBCAN_METHODS:-diamond,hmm,dbCANsub}"

DBCAN_LOG_FILE="${DBCAN_LOG_FILE:-${DIR_12_DBCAN}/${DBCAN_PREFIX}.run_dbcan.log}"
DBCAN_LOG_LEVEL="${DBCAN_LOG_LEVEL:-INFO}"

DBCAN_DIAMOND_EVALUE="${DBCAN_DIAMOND_EVALUE:-1e-102}"
DBCAN_HMM_EVALUE="${DBCAN_HMM_EVALUE:-1e-15}"
DBCAN_HMM_COVERAGE="${DBCAN_HMM_COVERAGE:-0.35}"
DBCAN_DBSUB_EVALUE="${DBCAN_DBSUB_EVALUE:-1e-15}"
DBCAN_DBSUB_COVERAGE="${DBCAN_DBSUB_COVERAGE:-0.35}"

DBCAN_LARGE="${DBCAN_LARGE:-FALSE}"
DBCAN_MEMORY_MONITORING="${DBCAN_MEMORY_MONITORING:-TRUE}"
DBCAN_MAX_RETRIES="${DBCAN_MAX_RETRIES:-3}"
DBCAN_MEMORY_SAFETY_FACTOR="${DBCAN_MEMORY_SAFETY_FACTOR:-0.5}"
DBCAN_MAX_MEMORY_USAGE="${DBCAN_MAX_MEMORY_USAGE:-0.8}"
DBCAN_CSV_BUFFER_SIZE="${DBCAN_CSV_BUFFER_SIZE:-5000}"

DBCAN_RUN_DEEPTMHMM="${DBCAN_RUN_DEEPTMHMM:-FALSE}"
DBCAN_DEEPTMHMM_DIR="${DBCAN_DEEPTMHMM_DIR:-}"
DBCAN_DEEPTMHMM_PYTHON="${DBCAN_DEEPTMHMM_PYTHON:-python}"

DBCAN_RUN_SIGNALP="${DBCAN_RUN_SIGNALP:-FALSE}"
DBCAN_SIGNALP_ORG="${DBCAN_SIGNALP_ORG:-euk}"

if ! is_true "${PREFLIGHT_ONLY}"; then
  mkdir -p "${DBCAN_OUT_DIR}"
fi

# ----------------------------- R plotting options -----------------------------

R_SCRIPT="${DIR_13_SCRIPTS}/rnaseq_single_sample_plots.R"
GENOME_GRAPHICS_R_SCRIPT="${DIR_13_SCRIPTS}/egapx_genome_annotation_graphics_karyoploter.R"

TOP_N_EXPRESSED="${TOP_N_EXPRESSED:-30}"
PSEUDOCOUNT="${PSEUDOCOUNT:-1}"
VOLCANO_LFC_THRESHOLD="${VOLCANO_LFC_THRESHOLD:-1}"
DESCRIPTIVE_VOLCANO_DENOMINATOR="${DESCRIPTIVE_VOLCANO_DENOMINATOR:-Muscle}"

PLOT_DPI="${PLOT_DPI:-600}"
PLOT_WIDTH="${PLOT_WIDTH:-10}"
PLOT_HEIGHT="${PLOT_HEIGHT:-7}"
PLOT_FORMATS="$(printf '%s' "${PLOT_FORMATS:-png,pdf,svg}" | tr '[:upper:]' '[:lower:]' | tr -d '[:space:]')"

# ----------------------------- EGAPx genome-graphics options ------------------

# All graphics outputs from the genome-annotation plotting script are written
# under DIR_13_PLOTS by using an output prefix inside that directory.
GENOME_PLOT_TOP_N="${GENOME_PLOT_TOP_N:-15}"
GENOME_GRAPHICS_PREFIX="${GENOME_GRAPHICS_PREFIX:-egapx_annotation_top${GENOME_PLOT_TOP_N}_functional}"
GENOME_GRAPHICS_OUT_PREFIX="${GENOME_GRAPHICS_OUT_PREFIX:-${DIR_13_PLOTS}/${GENOME_GRAPHICS_PREFIX}}"
GENOME_PLOT_WINDOW="${GENOME_PLOT_WINDOW:-500000}"
GENOME_PLOT_STEP="${GENOME_PLOT_STEP:-0}"
GENOME_CDS_METRIC="${GENOME_CDS_METRIC:-count}"
GENOME_LABEL_STYLE="${GENOME_LABEL_STYLE:-chr}"
GENOME_OFFICIAL_COUNT_SOURCE="${GENOME_OFFICIAL_COUNT_SOURCE:-egapx}"
GENOME_PLOT_COUNT_SOURCE="${GENOME_PLOT_COUNT_SOURCE:-gff}"
GENOME_PLOT_COUNT_COMPARISON="${GENOME_PLOT_COUNT_COMPARISON:-TRUE}"
GENOME_METHOD_SUMMARY_PLOTS="${GENOME_METHOD_SUMMARY_PLOTS:-TRUE}"
GENOME_METHOD_KARYOPLOTS="${GENOME_METHOD_KARYOPLOTS:-TRUE}"
GENOME_CREATE_ID_MAP="${GENOME_CREATE_ID_MAP:-TRUE}"
GENOME_FORCE_ID_MAP="${GENOME_FORCE_ID_MAP:-FALSE}"
GENOME_GC_RIBBON_K="${GENOME_GC_RIBBON_K:-7}"
GENOME_TOP_MARKER_N="${GENOME_TOP_MARKER_N:-40}"
GENOME_PLOT_FUNCTIONAL_COVERAGE="${GENOME_PLOT_FUNCTIONAL_COVERAGE:-TRUE}"
GENOME_PLOT_FUNCTIONAL_SUBTRACKS="${GENOME_PLOT_FUNCTIONAL_SUBTRACKS:-TRUE}"
GENOME_DPI="${GENOME_DPI:-600}"
GENOME_PLOT_FORMATS="$(printf '%s' "${GENOME_PLOT_FORMATS:-png,pdf,svg}" | tr '[:upper:]' '[:lower:]' | tr -d '[:space:]')"

validate_plot_formats() {
  local formats="$1"
  local format
  local parsed=()
  IFS=',' read -r -a parsed <<< "${formats}"
  ((${#parsed[@]})) || die "Plot format list must not be empty"
  for format in "${parsed[@]}"; do
    case "${format}" in
      png|pdf|svg)
        ;;
      *)
        die "Invalid plot format '${format}' in '${formats}'. Use a comma-separated subset of png,pdf,svg."
        ;;
    esac
  done
}
validate_plot_formats "${PLOT_FORMATS}"
validate_plot_formats "${GENOME_PLOT_FORMATS}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
R_SOURCE_DIR="${R_SOURCE_DIR:-${SCRIPT_DIR}/functional_annotation_v2_R}"
FEATURECOUNTS_R_SOURCE="${FEATURECOUNTS_R_SOURCE:-${R_SOURCE_DIR}/run_featurecounts_rsubread.R}"
RNASEQ_R_SOURCE="${RNASEQ_R_SOURCE:-${R_SOURCE_DIR}/rnaseq_single_sample_plots.R}"
GENOME_GRAPHICS_R_SOURCE="${GENOME_GRAPHICS_R_SOURCE:-${R_SOURCE_DIR}/egapx_genome_annotation_graphics_karyoploter.R}"

DEFAULT_ANNOTATION_PARALLEL_JOBS=1
if (( TOTAL_CPUS >= 16 && MEM_AVAILABLE_BYTES >= 128 * 1024 * 1024 * 1024 )); then
  DEFAULT_ANNOTATION_PARALLEL_JOBS=2
fi
ANNOTATION_PARALLEL_JOBS="${ANNOTATION_PARALLEL_JOBS:-${DEFAULT_ANNOTATION_PARALLEL_JOBS}}"
if [[ ! "${ANNOTATION_PARALLEL_JOBS}" =~ ^[0-9]+$ ]] || (( ANNOTATION_PARALLEL_JOBS < 1 )); then
  die "ANNOTATION_PARALLEL_JOBS must be a positive integer: ${ANNOTATION_PARALLEL_JOBS}"
fi
if (( ANNOTATION_PARALLEL_JOBS > TOTAL_CPUS )); then
  ANNOTATION_PARALLEL_JOBS="${TOTAL_CPUS}"
fi
ANNOTATION_THREADS_PER_TASK="${ANNOTATION_THREADS_PER_TASK:-$((TOTAL_CPUS / ANNOTATION_PARALLEL_JOBS))}"
if [[ ! "${ANNOTATION_THREADS_PER_TASK}" =~ ^[0-9]+$ ]] || (( ANNOTATION_THREADS_PER_TASK < 1 )); then
  die "ANNOTATION_THREADS_PER_TASK must be a positive integer: ${ANNOTATION_THREADS_PER_TASK}"
fi
(( ANNOTATION_THREADS_PER_TASK < 1 )) && ANNOTATION_THREADS_PER_TASK=1
if (( ANNOTATION_PARALLEL_JOBS * ANNOTATION_THREADS_PER_TASK > TOTAL_CPUS )); then
  ANNOTATION_THREADS_PER_TASK=$((TOTAL_CPUS / ANNOTATION_PARALLEL_JOBS))
  (( ANNOTATION_THREADS_PER_TASK < 1 )) && ANNOTATION_THREADS_PER_TASK=1
fi

# Optional manual size overrides. If unset, dimensions are computed from
# GENOME_PLOT_TOP_N so plots grow with the number of contigs/scaffolds/chromosomes.
GENOME_PLOT_WIDTH="${GENOME_PLOT_WIDTH:-}"
GENOME_PLOT_HEIGHT="${GENOME_PLOT_HEIGHT:-}"
GENOME_LINEAR_WIDTH="${GENOME_LINEAR_WIDTH:-}"
GENOME_LINEAR_HEIGHT="${GENOME_LINEAR_HEIGHT:-}"
GENOME_KARYO_WIDTH="${GENOME_KARYO_WIDTH:-}"
GENOME_KARYO_HEIGHT="${GENOME_KARYO_HEIGHT:-}"
GENOME_KARYO_BOTTOM_OMA="${GENOME_KARYO_BOTTOM_OMA:-9}"
GENOME_KARYO_LEGEND_INSET="${GENOME_KARYO_LEGEND_INSET:-0.045}"
GENOME_KARYO_BASE_CEX="${GENOME_KARYO_BASE_CEX:-0.55}"
GENOME_KARYO_LABEL_CEX="${GENOME_KARYO_LABEL_CEX:-0.68}"
GENOME_KARYO_LEGEND_CEX="${GENOME_KARYO_LEGEND_CEX:-0.80}"

GENOME_ID_MAP="${GENOME_ID_MAP:-${DIR_13_ANNOTATION}/egapx_gene_transcript_protein_map.tsv}"
EGAPX_FEATURE_COUNTS="${EGAPX_FEATURE_COUNTS:-${EGAPX_OUT_DIR}/stats/feature_counts.txt}"

# ==============================================================================
# 1) Environment check functions
# ==============================================================================

check_r_environment() {
  environment_was_checked "r:${R_ENV}:${USE_R_MAMBA_ENV}" && return 0
  log "Checking R environment"
  local required_packages=("${R_REQUIRED_PACKAGES[@]}")

  if is_true "${USE_R_MAMBA_ENV}"; then
    check_executable "${MAMBA_BIN}"
    log "R execution mode: mamba environment '${R_ENV}'"
    "${MAMBA_BIN}" run -n "${R_ENV}" Rscript --version >&2
  else
    check_executable "${RSCRIPT_BIN}"
    log "R execution mode: direct Rscript '${RSCRIPT_BIN}'"
    "${RSCRIPT_BIN}" --version >&2
  fi

  if is_true "${RUN_TIDYPLOTS}"; then
    required_packages+=(tidyplots)
  fi

  local pkg_expr
  local pkg_vec
  printf -v pkg_vec '"%s",' "${required_packages[@]}"
  pkg_vec="${pkg_vec%,}"
  pkg_expr="pkgs <- c(${pkg_vec}); missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]; if (length(missing)) { cat('Missing R packages:', paste(missing, collapse=', '), '\n', file=stderr()); quit(status=1) } else { cat('All required R packages are available:', paste(pkgs, collapse=', '), '\n') }"
  run_rscript -e "${pkg_expr}"
  mark_environment_checked "r:${R_ENV}:${USE_R_MAMBA_ENV}"
}

check_mamba_env_command() {
  local env_name="$1"
  local cmd_name="$2"

  environment_was_checked "mamba:${env_name}:${cmd_name}" && return 0
  check_executable "${MAMBA_BIN}"

  "${MAMBA_BIN}" run -n "${env_name}" which "${cmd_name}" >/dev/null 2>&1 \
    || die "Command '${cmd_name}' was not found inside mamba environment '${env_name}'"

  log "Found ${cmd_name} in mamba environment: ${env_name}"
  mark_environment_checked "mamba:${env_name}:${cmd_name}"
}

check_diamond_environment() {
  check_mamba_env_command "${DIAMOND_ENV}" diamond
}

check_kofamscan_environment() {
  check_mamba_env_command "${KOFAM_ENV}" exec_annotation
}

check_dbcan_environment() {
  check_mamba_env_command "${DBCAN_ENV}" run_dbcan
}

check_dbcan_annotation_command() {
  environment_was_checked "dbcan:${DBCAN_ENV}:CAZyme_annotation" && return 0
  check_dbcan_environment

  "${MAMBA_BIN}" run -n "${DBCAN_ENV}" run_dbcan CAZyme_annotation --help >/dev/null 2>&1 \
    || die "Your dbCAN installation does not expose 'run_dbcan CAZyme_annotation'"
  mark_environment_checked "dbcan:${DBCAN_ENV}:CAZyme_annotation"
}

validate_dbcan_options() {
  local dbcan_methods_array=()
  local method

  case "${DBCAN_MODE}" in
    prok|meta|protein)
      ;;
    *)
      die "Invalid DBCAN_MODE='${DBCAN_MODE}'. Use prok, meta, or protein."
      ;;
  esac

  case "${DBCAN_LOG_LEVEL}" in
    DEBUG|INFO|WARNING|ERROR|CRITICAL)
      ;;
    *)
      die "Invalid DBCAN_LOG_LEVEL='${DBCAN_LOG_LEVEL}'. Use DEBUG, INFO, WARNING, ERROR, or CRITICAL."
      ;;
  esac

  case "${DBCAN_SIGNALP_ORG}" in
    other|euk)
      ;;
    *)
      die "Invalid DBCAN_SIGNALP_ORG='${DBCAN_SIGNALP_ORG}'. Use other or euk."
      ;;
  esac

  IFS=',' read -r -a dbcan_methods_array <<< "${DBCAN_METHODS}"
  for method in "${dbcan_methods_array[@]}"; do
    case "${method}" in
      diamond|hmm|dbCANsub)
        ;;
      *)
        die "Invalid dbCAN method '${method}' in DBCAN_METHODS='${DBCAN_METHODS}'. Use comma-separated choices from diamond,hmm,dbCANsub."
        ;;
    esac
  done
}

required_dbcan_database_files() {
  local dbcan_methods_array=()
  local method

  IFS=',' read -r -a dbcan_methods_array <<< "${DBCAN_METHODS}"
  for method in "${dbcan_methods_array[@]}"; do
    case "${method}" in
      diamond)
        printf '%s\n' "${DBCAN_DB_DIR}/CAZy.dmnd"
        ;;
      hmm)
        printf '%s\n' "${DBCAN_DB_DIR}/dbCAN.hmm"
        ;;
      dbCANsub)
        printf '%s\n' "${DBCAN_DB_DIR}/dbCAN-sub.hmm"
        ;;
    esac
  done
}

get_organ_from_accession() {
  local acc="$1"
  awk -v acc="${acc}" 'BEGIN{FS=OFS="\t"} NR>1 && $1==acc{print $2; exit}' "${SRA_METADATA}"
}

create_sample_metadata() {
  check_file "${SRA_METADATA}"

  awk '
    BEGIN { FS = OFS = "\t" }
    NR == 1 { print "sample", "accession", "organ"; next }
    NR > 1 { print $1 "." $2, $1, $2 }
  ' "${SRA_METADATA}" > "${SAMPLE_METADATA}"

  log "Created sample metadata: ${SAMPLE_METADATA}"
}

find_sra_path() {
  local acc="$1"

  if [[ -s "${DIR_02_SRA}/${acc}/${acc}.sra" ]]; then
    printf "%s\n" "${DIR_02_SRA}/${acc}/${acc}.sra"
  elif [[ -s "${DIR_02_SRA}/${acc}.sra" ]]; then
    printf "%s\n" "${DIR_02_SRA}/${acc}.sra"
  elif [[ -d "${DIR_02_SRA}/${acc}" ]]; then
    printf "%s\n" "${DIR_02_SRA}/${acc}"
  else
    die "Could not locate downloaded SRA file for ${acc} in ${DIR_02_SRA}"
  fi
}

fastq_pair_is_valid() {
  local acc="$1"
  local extension="$2"
  local r1="${DIR_03_FASTQ}/${acc}_1.fastq${extension}"
  local r2="${DIR_03_FASTQ}/${acc}_2.fastq${extension}"
  if [[ "${extension}" == ".gz" ]]; then
    gzip_is_valid "${r1}" && gzip_is_valid "${r2}"
  else
    [[ -s "${r1}" && -s "${r2}" ]]
  fi
}

# ==============================================================================
# 2) SRA, FASTQ, fastp, and mapping functions
# ==============================================================================

download_one_sra() {
  local acc="$1"
  local organ
  organ="$(get_organ_from_accession "${acc}")"
  [[ -n "${organ}" ]] || die "Could not find organ metadata for ${acc}"

  log "Downloading SRA accession: ${acc} (${organ})"

  "${PREFETCH_BIN}" \
    --max-size "${PREFETCH_MAX_SIZE}" \
    --progress \
    --resume yes \
    --output-directory "${DIR_02_SRA}" \
    --log-level info \
    "${acc}"

  log "Finished downloading: ${acc} (${organ})"
}

convert_one_sra() {
  local acc="$1"
  local organ
  local sra_path
  local tmp_acc
  local force_opts=()

  organ="$(get_organ_from_accession "${acc}")"
  [[ -n "${organ}" ]] || die "Could not find organ metadata for ${acc}"

  if fastq_pair_is_valid "${acc}" ".gz"; then
    log "Skipping FASTQ conversion; valid compressed mates already exist for ${acc}"
    return 0
  fi

  if fastq_pair_is_valid "${acc}" ""; then
    log "Skipping FASTQ conversion; uncompressed mates already exist for ${acc}"
    return 0
  fi

  sra_path="$(find_sra_path "${acc}")"
  tmp_acc="${tmpdir}/fasterq_${acc}"
  mkdir -p "${tmp_acc}"

  if is_true "${FASTERQ_FORCE:-FALSE}"; then
    force_opts+=(--force)
  elif [[ -e "${DIR_03_FASTQ}/${acc}_1.fastq" || -e "${DIR_03_FASTQ}/${acc}_2.fastq" ]]; then
    die "Incomplete FASTQ pair exists for ${acc}; set FASTERQ_FORCE=TRUE to replace it"
  fi

  log "Converting ${acc} (${organ}) to FASTQ"

  "${FASTERQ_DUMP_BIN}" "${sra_path}" \
    --outdir "${DIR_03_FASTQ}" \
    --temp "${tmp_acc}" \
    --threads "${NTHREADS}" \
    --mem "${FASTERQ_MEM}" \
    --skip-technical \
    --split-files \
    --split-3 \
    --progress \
    "${force_opts[@]}"

  fastq_pair_is_valid "${acc}" "" \
    || die "fasterq-dump did not create a complete paired FASTQ set for ${acc}"

  log "Finished FASTQ conversion: ${acc} (${organ})"
}

convert_and_compress_one_sra() {
  local acc="$1"
  local fq

  convert_one_sra "${acc}"
  if fastq_pair_is_valid "${acc}" ".gz"; then
    return 0
  fi
  for fq in "${DIR_03_FASTQ}/${acc}_1.fastq" "${DIR_03_FASTQ}/${acc}_2.fastq"; do
    check_file "${fq}"
    compress_one_fastq "${fq}"
  done
  fastq_pair_is_valid "${acc}" ".gz" \
    || die "Compression did not create a valid paired FASTQ set for ${acc}"
}

compress_one_fastq() {
  local fq="$1"

  if gzip_is_valid "${fq}.gz"; then
    log "Skipping compression; already exists: ${fq}.gz"
    return 0
  fi

  if [[ -e "${fq}.gz" ]]; then
    die "Existing compressed FASTQ is invalid and will not be overwritten automatically: ${fq}.gz"
  fi

  log "Compressing FASTQ: ${fq}"
  env -u PIGZ "${PIGZ_BIN}" -p "${NTHREADS}" --best "${fq}"
  gzip_is_valid "${fq}.gz" || die "pigz created an invalid gzip file: ${fq}.gz"
}

run_fastp_one_sample() {
  local acc="$1"
  local organ
  local r1
  local r2
  local out1
  local out2
  local json
  local html
  local failed_out
  local unpaired1
  local unpaired2
  local complexity_opts=()
  local overrep_opts=()
  local polyg_opts=()
  local polyx_opts=()
  local correction_opts=()
  local failed_opts=()
  local unpaired_opts=()
  local overwrite_opts=()
  local manifest
  local signature
  local run_dir
  local run_out1
  local run_out2
  local run_json
  local run_html

  organ="$(get_organ_from_accession "${acc}")"
  [[ -n "${organ}" ]] || die "Could not find organ metadata for ${acc}"

  r1="${DIR_03_FASTQ}/${acc}_1.fastq.gz"
  r2="${DIR_03_FASTQ}/${acc}_2.fastq.gz"

  out1="${DIR_04_TRIMMED}/${acc}.${organ}_1.trimmed.fastq.gz"
  out2="${DIR_04_TRIMMED}/${acc}.${organ}_2.trimmed.fastq.gz"
  json="${DIR_04_FASTP_JSON}/${acc}.${organ}.fastp.json"
  html="${DIR_04_FASTP_HTML}/${acc}.${organ}.fastp.html"
  failed_out="${DIR_04_FASTP_FAILED}/${acc}.${organ}.failed.fastq.gz"
  unpaired1="${DIR_04_FASTP_UNPAIRED}/${acc}.${organ}_unpaired1.fastq.gz"
  unpaired2="${DIR_04_FASTP_UNPAIRED}/${acc}.${organ}_unpaired2.fastq.gz"
  run_dir="${tmpdir}/fastp_${acc}.$$"
  run_out1="${run_dir}/$(basename "${out1}")"
  run_out2="${run_dir}/$(basename "${out2}")"
  run_json="${run_dir}/$(basename "${json}")"
  run_html="${run_dir}/$(basename "${html}")"

  check_file "${r1}"
  check_file "${r2}"
  gzip_is_valid "${r1}" || die "Invalid gzip input for fastp: ${r1}"
  gzip_is_valid "${r2}" || die "Invalid gzip input for fastp: ${r2}"

  manifest="$(stage_manifest_path fastp "${acc}")"
  signature="$(make_signature \
    "$(file_identity "${r1}")" "$(file_identity "${r2}")" "$(file_identity "${FASTP_BIN}")" \
    "${FASTP_QUAL}" "${FASTP_LENGTH_REQUIRED}" "${FASTP_N_BASE_LIMIT}" \
    "${FASTP_UNQUALIFIED_PERCENT_LIMIT}" "${FASTP_COMPLEXITY_FILTER}" "${FASTP_COMPLEXITY_THRESHOLD}" \
    "${FASTP_OVERREPRESENTATION}" "${FASTP_OVERREPRESENTATION_SAMPLING}" \
    "${FASTP_POLYG}" "${FASTP_POLYG_MIN_LEN}" "${FASTP_POLYX}" "${FASTP_POLYX_MIN_LEN}" \
    "${FASTP_CORRECTION}" "${FASTP_SAVE_FAILED}" "${FASTP_SAVE_UNPAIRED}" \
    "${FASTP_COMPRESSION}" "${NTHREADS}")"

  if manifest_matches "${manifest}" "${signature}" "${out1}" "${out2}" "${json}" "${html}" \
    && gzip_is_valid "${out1}" && gzip_is_valid "${out2}" && json_is_valid "${json}"; then
    log "Skipping fastp; validated outputs and matching manifest exist for ${acc} (${organ})"
    return 0
  fi

  if is_true "${FASTP_COMPLEXITY_FILTER}"; then
    complexity_opts+=(
      --low_complexity_filter
      --complexity_threshold "${FASTP_COMPLEXITY_THRESHOLD}"
    )
  fi

  if is_true "${FASTP_OVERREPRESENTATION}"; then
    overrep_opts+=(
      --overrepresentation_analysis
      --overrepresentation_sampling "${FASTP_OVERREPRESENTATION_SAMPLING}"
    )
  fi

  if is_true "${FASTP_POLYG}"; then
    polyg_opts+=(--trim_poly_g --poly_g_min_len "${FASTP_POLYG_MIN_LEN}")
  fi

  if is_true "${FASTP_POLYX}"; then
    polyx_opts+=(--trim_poly_x --poly_x_min_len "${FASTP_POLYX_MIN_LEN}")
  fi

  if is_true "${FASTP_CORRECTION}"; then
    correction_opts+=(--correction)
  fi

  if is_true "${FASTP_SAVE_FAILED}"; then
    failed_opts+=(--failed_out "${failed_out}")
  fi

  if is_true "${FASTP_SAVE_UNPAIRED}"; then
    unpaired_opts+=(--unpaired1 "${unpaired1}" --unpaired2 "${unpaired2}")
  fi

  if is_true "${FASTP_DONT_OVERWRITE}"; then
    overwrite_opts+=(--dont_overwrite)
  fi

  log "Running fastp for ${acc} (${organ})"
  log "fastp outputs: ${out1}, ${out2}, ${json}, ${html}"
  rm -rf "${run_dir}"
  mkdir -p "${run_dir}"

  "${FASTP_BIN}" \
    --in1 "${r1}" \
    --in2 "${r2}" \
    --out1 "${run_out1}" \
    --out2 "${run_out2}" \
    --detect_adapter_for_pe \
    --cut_tail \
    --cut_tail_window_size 4 \
    --cut_tail_mean_quality "${FASTP_QUAL}" \
    --qualified_quality_phred "${FASTP_QUAL}" \
    --unqualified_percent_limit "${FASTP_UNQUALIFIED_PERCENT_LIMIT}" \
    --n_base_limit "${FASTP_N_BASE_LIMIT}" \
    --length_required "${FASTP_LENGTH_REQUIRED}" \
    --compression "${FASTP_COMPRESSION}" \
    --thread "${NTHREADS}" \
    --json "${run_json}" \
    --html "${run_html}" \
    --report_title "${acc}.${organ}" \
    "${complexity_opts[@]}" \
    "${overrep_opts[@]}" \
    "${polyg_opts[@]}" \
    "${polyx_opts[@]}" \
    "${correction_opts[@]}" \
    "${failed_opts[@]}" \
    "${unpaired_opts[@]}" \
    "${overwrite_opts[@]}"

  gzip_is_valid "${run_out1}" || die "fastp created an invalid gzip file: ${run_out1}"
  gzip_is_valid "${run_out2}" || die "fastp created an invalid gzip file: ${run_out2}"
  json_is_valid "${run_json}" || die "fastp created invalid JSON: ${run_json}"
  check_file "${run_html}"
  mv -f "${run_out1}" "${out1}"
  mv -f "${run_out2}" "${out2}"
  mv -f "${run_json}" "${json}"
  mv -f "${run_html}" "${html}"
  write_manifest "${manifest}" "${signature}" fastp "${out1}" "${out2}" "${json}" "${html}"
  log "Finished fastp for ${acc} (${organ})"
}

run_multiqc_report() {
  local label="$1"
  local outdir="${DIR_05_MULTIQC}/${label}"
  local input_dirs=()
  local verbose_opts=()

  mkdir -p "${outdir}"

  log "Running MultiQC report: ${label}"

  case "${label}" in
    after_fastp)
      input_dirs=("${DIR_04_FASTP_JSON}" "${DIR_04_FASTP_HTML}")
      ;;
    final)
      input_dirs=("${DIR_04_FASTP_JSON}" "${DIR_04_FASTP_HTML}" "${DIR_06_SUMMARY}" "${DIR_07_COUNTS}")
      ;;
    *)
      die "Unknown MultiQC report label: ${label}"
      ;;
  esac

  is_true "${MULTIQC_VERBOSE:-FALSE}" && verbose_opts+=(--verbose)

  "${MULTIQC_BIN}" \
    "${input_dirs[@]}" \
    --outdir "${outdir}" \
    --filename "multiqc_${label}" \
    --interactive \
    --force \
    "${verbose_opts[@]}"

  log "Finished MultiQC report: ${label}"
}

check_hisat2_index() {
  local idx="$1"
  local suffix
  local part
  if [[ -s "${idx}.1.ht2" ]]; then
    suffix="ht2"
  elif [[ -s "${idx}.1.ht2l" ]]; then
    suffix="ht2l"
  else
    return 1
  fi
  for part in 1 2 3 4 5 6 7 8; do
    [[ -s "${idx}.${part}.${suffix}" ]] || return 1
  done
}

build_hisat2_index() {
  local manifest
  local signature
  local tmp_index_prefix="${tmpdir}/hisat2_index.$$.tmp"
  local suffix
  local part
  local expected_outputs=()
  local run_expected_outputs=()
  local run_output_dir="${tmpdir}/interproscan_output.$$"

  manifest="$(stage_manifest_path hisat2_index global)"
  signature="$(make_signature "$(file_identity "${EGAPX_GENOME}")" \
    "$(file_identity "${HISAT2_BUILD_BIN}")" "${NTHREADS}")"

  if check_hisat2_index "${HISAT2_INDEX_PREFIX}" \
    && [[ "$(awk -F '\t' '$1 == "signature" { print $2; exit }' "${manifest}" 2>/dev/null || true)" == "${signature}" ]]; then
    log "Validated HISAT2 index and matching manifest already exist: ${HISAT2_INDEX_PREFIX}"
    return 0
  fi

  check_file "${EGAPX_GENOME}"
  check_executable "${HISAT2_BUILD_BIN}"

  log "Building HISAT2 index from EGAPx genome"
  rm -f "${tmp_index_prefix}".*.ht2 "${tmp_index_prefix}".*.ht2l

  "${HISAT2_BUILD_BIN}" \
    -p "${NTHREADS}" \
    "${EGAPX_GENOME}" \
    "${tmp_index_prefix}"

  check_hisat2_index "${tmp_index_prefix}" || die "HISAT2 created an incomplete index: ${tmp_index_prefix}"
  if [[ -s "${tmp_index_prefix}.1.ht2" ]]; then
    suffix="ht2"
  else
    suffix="ht2l"
  fi
  mkdir -p "$(dirname "${HISAT2_INDEX_PREFIX}")"
  rm -f "${HISAT2_INDEX_PREFIX}".*.ht2 "${HISAT2_INDEX_PREFIX}".*.ht2l
  for part in 1 2 3 4 5 6 7 8; do
    mv -f "${tmp_index_prefix}.${part}.${suffix}" "${HISAT2_INDEX_PREFIX}.${part}.${suffix}"
    expected_outputs+=("${HISAT2_INDEX_PREFIX}.${part}.${suffix}")
  done
  write_manifest "${manifest}" "${signature}" hisat2_index "${expected_outputs[@]}"

  log "Finished HISAT2 index"
}

map_one_sample_hisat2() {
  local acc="$1"
  local organ
  local r1
  local r2
  local out_bam
  local tmp_bam
  local tmp_prefix
  local summary_file
  local metrics_file
  local stderr_log
  local rg_id
  local rg_sm
  local rg_lb
  local hisat2_opts=()
  local preset_opts=()
  local extra_opts=()
  local manifest
  local signature

  organ="$(get_organ_from_accession "${acc}")"
  [[ -n "${organ}" ]] || die "Could not find organ metadata for ${acc}"

  case "${MAPPING_INPUT_SOURCE}" in
    auto)
      if gzip_is_valid "${DIR_04_TRIMMED}/${acc}.${organ}_1.trimmed.fastq.gz" \
        && gzip_is_valid "${DIR_04_TRIMMED}/${acc}.${organ}_2.trimmed.fastq.gz"; then
        r1="${DIR_04_TRIMMED}/${acc}.${organ}_1.trimmed.fastq.gz"
        r2="${DIR_04_TRIMMED}/${acc}.${organ}_2.trimmed.fastq.gz"
      else
        r1="${DIR_03_FASTQ}/${acc}_1.fastq.gz"
        r2="${DIR_03_FASTQ}/${acc}_2.fastq.gz"
      fi
      ;;
    trimmed)
      r1="${DIR_04_TRIMMED}/${acc}.${organ}_1.trimmed.fastq.gz"
      r2="${DIR_04_TRIMMED}/${acc}.${organ}_2.trimmed.fastq.gz"
      ;;
    raw)
      r1="${DIR_03_FASTQ}/${acc}_1.fastq.gz"
      r2="${DIR_03_FASTQ}/${acc}_2.fastq.gz"
      ;;
    *)
      die "Invalid MAPPING_INPUT_SOURCE='${MAPPING_INPUT_SOURCE}'. Use auto, trimmed, or raw."
      ;;
  esac

  out_bam="${DIR_06_BAM}/${acc}.${organ}.hisat2.sorted.bam"
  tmp_bam="${tmpdir}/${acc}.${organ}.hisat2.sorted.bam"
  tmp_prefix="${tmpdir}/${acc}.${organ}.samtools_sort"

  summary_file="${DIR_06_SUMMARY}/${acc}.${organ}.hisat2.summary.txt"
  metrics_file="${DIR_06_SUMMARY}/${acc}.${organ}.hisat2.metrics.txt"
  stderr_log="${DIR_06_SUMMARY}/${acc}.${organ}.hisat2.stderr.log"

  rg_id="${acc}"
  rg_sm="${acc}.${organ}"
  rg_lb="${organ}_RNAseq"

  check_file "${r1}"
  check_file "${r2}"
  gzip_is_valid "${r1}" || die "Invalid gzip input for HISAT2: ${r1}"
  gzip_is_valid "${r2}" || die "Invalid gzip input for HISAT2: ${r2}"

  manifest="$(stage_manifest_path hisat2_mapping "${acc}")"
  signature="$(make_signature \
    "$(file_identity "${r1}")" "$(file_identity "${r2}")" "$(file_identity "${EGAPX_GENOME}")" \
    "$(file_identity "${HISAT2_INDEX_PREFIX}.1.ht2")" "$(file_identity "${HISAT2_INDEX_PREFIX}.1.ht2l")" \
    "${HISAT2_PRESET}" "${HISAT2_EXTRA_OPTS}" "${KNOWN_SPLICESITES}" "${RNA_STRANDNESS}" \
    "${MAPPING_INPUT_SOURCE}" "${RG_PLATFORM}" "${RG_CENTER}" \
    "${HISAT2_THREADS}" "${SAMTOOLS_SORT_THREADS}" "${SAMTOOLS_SORT_MEM}")"

  if manifest_matches "${manifest}" "${signature}" "${out_bam}" \
    "${DIR_06_SUMMARY}/${acc}.${organ}.flagstat.txt" \
    "${DIR_06_SUMMARY}/${acc}.${organ}.samtools.stats.txt" \
    && bam_is_valid "${out_bam}"; then
    log "Skipping HISAT2 mapping; validated BAM and matching manifest exist for ${acc} (${organ})"
    return 0
  fi

  if [[ -e "${out_bam}" || -e "${out_bam}.bai" || -e "${out_bam}.csi" ]]; then
    WARN "Existing HISAT2 BAM or index is missing/invalid and will be replaced for ${acc} (${organ})"
    rm -f "${out_bam}" "${out_bam}.bai" "${out_bam}.csi"
  fi

  rm -f "${tmp_bam}" "${tmp_bam}.bai" "${tmp_bam}.csi"

  if [[ -n "${HISAT2_PRESET}" ]]; then
    IFS=' ' read -r -a preset_opts <<< "${HISAT2_PRESET}"
    hisat2_opts+=("${preset_opts[@]}")
  fi

  if [[ -n "${HISAT2_EXTRA_OPTS}" ]]; then
    IFS=' ' read -r -a extra_opts <<< "${HISAT2_EXTRA_OPTS}"
    hisat2_opts+=("${extra_opts[@]}")
  fi

  if [[ -n "${KNOWN_SPLICESITES}" ]]; then
    hisat2_opts+=(--known-splicesite-infile "${KNOWN_SPLICESITES}")
  fi

  if [[ -n "${RNA_STRANDNESS}" ]]; then
    hisat2_opts+=(--rna-strandness "${RNA_STRANDNESS}")
  fi

  log "Mapping ${acc} (${organ}) with HISAT2"

  "${HISAT2_BIN}" \
    -x "${HISAT2_INDEX_PREFIX}" \
    -1 "${r1}" \
    -2 "${r2}" \
    -p "${HISAT2_THREADS}" \
    "${hisat2_opts[@]}" \
    --summary-file "${summary_file}" \
    --met-file "${metrics_file}" \
    --rg-id "${rg_id}" \
    --rg "SM:${rg_sm}" \
    --rg "PL:${RG_PLATFORM}" \
    --rg "LB:${rg_lb}" \
    --rg "CN:${RG_CENTER}" \
    2> "${stderr_log}" \
  | "${SAMTOOLS_BIN}" sort \
      --threads "${SAMTOOLS_SORT_THREADS}" \
      --reference "${EGAPX_GENOME}" \
      -m "${SAMTOOLS_SORT_MEM}" \
      -l 9 \
      -O BAM \
      -T "${tmp_prefix}" \
      -o "${tmp_bam}" \
      --write-index \
      -

  "${SAMTOOLS_BIN}" quickcheck "${tmp_bam}" \
    || die "HISAT2/samtools produced an invalid BAM for ${acc} (${organ}): ${tmp_bam}"

  mv -f "${tmp_bam}" "${out_bam}"
  rm -f "${out_bam}.bai" "${out_bam}.csi"

  if [[ -s "${tmp_bam}.bai" ]]; then
    mv -f "${tmp_bam}.bai" "${out_bam}.bai"
  elif [[ -s "${tmp_bam}.csi" ]]; then
    mv -f "${tmp_bam}.csi" "${out_bam}.csi"
  else
    "${SAMTOOLS_BIN}" index --threads "${SAMTOOLS_AUX_THREADS}" "${out_bam}"
  fi

  "${SAMTOOLS_BIN}" flagstat \
    --threads "${SAMTOOLS_AUX_THREADS}" \
    "${out_bam}" \
    > "${DIR_06_SUMMARY}/${acc}.${organ}.flagstat.txt"

  "${SAMTOOLS_BIN}" stats \
    --threads "${SAMTOOLS_AUX_THREADS}" \
    --reference "${EGAPX_GENOME}" \
    "${out_bam}" \
    > "${DIR_06_SUMMARY}/${acc}.${organ}.samtools.stats.txt"

  write_manifest "${manifest}" "${signature}" hisat2_mapping \
    "${out_bam}" "${DIR_06_SUMMARY}/${acc}.${organ}.flagstat.txt" \
    "${DIR_06_SUMMARY}/${acc}.${organ}.samtools.stats.txt"
  log "Finished HISAT2 mapping: ${acc} (${organ})"
}

# ==============================================================================
# 3) Counting and annotation functions
# ==============================================================================

write_featurecounts_r_script() {
  check_file "${FEATURECOUNTS_R_SOURCE}"
  mkdir -p "$(dirname "${FEATURECOUNTS_R_SCRIPT}")"
  cp -f "${FEATURECOUNTS_R_SOURCE}" "${FEATURECOUNTS_R_SCRIPT}"
  chmod +x "${FEATURECOUNTS_R_SCRIPT}"
  log "Installed R script: ${FEATURECOUNTS_R_SCRIPT}"
}

run_featurecounts() {
  local manifest
  local signature
  local bam
  local bam_identities=()
  local run_counts="${tmpdir}/featurecounts.counts.$$"
  local run_matrix="${tmpdir}/featurecounts.matrix.$$"
  local run_summary="${tmpdir}/featurecounts.summary.$$"
  check_file "${EGAPX_GTF}"
  check_executable "${MAMBA_BIN}"

  mapfile -t BAM_FILES < <(
    find "${DIR_06_BAM}" \
      -maxdepth 1 \
      -type f \
      -name "*.hisat2.sorted.bam" \
      | sort
  )

  if [[ "${#BAM_FILES[@]}" -eq 0 ]]; then
    die "No sorted BAM files found in ${DIR_06_BAM}"
  fi

  for bam in "${BAM_FILES[@]}"; do
    bam_identities+=("$(file_identity "${bam}")")
  done

  manifest="$(stage_manifest_path featurecounts global)"
  signature="$(make_signature \
    "$(file_identity "${EGAPX_GTF}")" \
    "${bam_identities[@]}" \
    "${FEATURECOUNTS_FEATURE_TYPE}" "${FEATURECOUNTS_GROUP_ATTRIBUTE}" \
    "${FEATURECOUNTS_STRANDNESS}" "${FEATURECOUNTS_IS_PAIRED}" \
    "${FEATURECOUNTS_ALLOW_MULTI_OVERLAP}" "${FEATURECOUNTS_COUNT_MULTIMAPPERS}")"

  if manifest_matches "${manifest}" "${signature}" "${COUNTS_RAW}" "${COUNTS_MATRIX}" "${COUNTS_RAW}.summary" \
    && count_outputs_valid; then
    log "Skipping Rsubread::featureCounts; valid outputs already exist"
    log "Counts table: ${COUNTS_RAW}"
    log "Count matrix: ${COUNTS_MATRIX}"
    return 0
  fi

  if [[ -e "${COUNTS_RAW}" || -e "${COUNTS_MATRIX}" || -e "${COUNTS_RAW}.summary" ]]; then
    WARN "Existing featureCounts outputs are incomplete or invalid and will be replaced"
    rm -f "${COUNTS_RAW}" "${COUNTS_MATRIX}" "${COUNTS_RAW}.summary"
  fi

  write_featurecounts_r_script
  rm -f "${run_counts}" "${run_matrix}" "${run_summary}"

  log "Running Rsubread::featureCounts through ${R_ENV}"
  log "BAM directory: ${DIR_06_BAM}"
  log "GTF: ${EGAPX_GTF}"
  log "Threads: ${TOTAL_CPUS}"
  log "featureCounts strandness: ${FEATURECOUNTS_STRANDNESS}"
  log "featureCounts feature type: ${FEATURECOUNTS_FEATURE_TYPE}"
  log "featureCounts group attribute: ${FEATURECOUNTS_GROUP_ATTRIBUTE}"

  run_rscript "${FEATURECOUNTS_R_SCRIPT}" \
    --bam_dir "${DIR_06_BAM}" \
    --gtf "${EGAPX_GTF}" \
    --out_counts "${run_counts}" \
    --out_matrix "${run_matrix}" \
    --out_summary "${run_summary}" \
    --threads "${TOTAL_CPUS}" \
    --feature_type "${FEATURECOUNTS_FEATURE_TYPE}" \
    --attr_type "${FEATURECOUNTS_GROUP_ATTRIBUTE}" \
    --strandness "${FEATURECOUNTS_STRANDNESS}" \
    --paired "${FEATURECOUNTS_IS_PAIRED}" \
    --count_read_pairs "TRUE" \
    --allow_multi_overlap "${FEATURECOUNTS_ALLOW_MULTI_OVERLAP}" \
    --count_multimappers "${FEATURECOUNTS_COUNT_MULTIMAPPERS}"

  check_file "${run_counts}"
  check_file "${run_matrix}"
  check_file "${run_summary}"
  mv -f "${run_counts}" "${COUNTS_RAW}"
  mv -f "${run_matrix}" "${COUNTS_MATRIX}"
  mv -f "${run_summary}" "${COUNTS_RAW}.summary"

  write_manifest "${manifest}" "${signature}" featurecounts \
    "${COUNTS_RAW}" "${COUNTS_MATRIX}" "${COUNTS_RAW}.summary"
  log "Finished Rsubread::featureCounts"
}

ensure_eggnog_data() {
  local required_files=()
  local missing_files=()
  local force_flag=()
  local file

  if is_true "${EGGNOG_DOWNLOAD_DATA}"; then
    mkdir -p "${EGGNOG_DATA_DIR}"
  else
    check_dir "${EGGNOG_DATA_DIR}"
  fi

  case "${EGGNOG_METHOD}" in
    diamond)
      required_files+=("${EGGNOG_DATA_DIR}/eggnog_proteins.dmnd")
      ;;
    hmmer|no_search|cache|none)
      ;;
    *)
      WARN "No eggNOG database preflight rule for EGGNOG_METHOD=${EGGNOG_METHOD}"
      ;;
  esac

  if is_true "${EGGNOG_PFAM_REALIGN}" || [[ "${EGGNOG_PFAM_REALIGN}" != "none" ]]; then
    required_files+=("${EGGNOG_DATA_DIR}/Pfam-A.hmm.h3f")
  fi

  for file in "${required_files[@]}"; do
    [[ -s "${file}" ]] || missing_files+=("${file}")
  done

  if [[ "${#missing_files[@]}" -eq 0 ]] && ! is_true "${EGGNOG_FORCE_DOWNLOAD}"; then
    log "eggNOG data preflight passed: ${EGGNOG_DATA_DIR}"
    return 0
  fi

  if ! is_true "${EGGNOG_DOWNLOAD_DATA}"; then
    if [[ "${#missing_files[@]}" -gt 0 ]]; then
      printf 'Missing eggNOG database file(s):\n' >&2
      printf '  %s\n' "${missing_files[@]}" >&2
    fi
    die "Install/update eggNOG data only when intended: EGGNOG_DOWNLOAD_DATA=TRUE ${MAMBA_BIN} run -n ${EGGNOG_ENV} download_eggnog_data.py --data_dir ${EGGNOG_DATA_DIR} -y"
  fi

  is_true "${EGGNOG_FORCE_DOWNLOAD}" && force_flag+=(-f)

  log "Downloading eggNOG data into ${EGGNOG_DATA_DIR}"
  WARN "This may replace eggNOG files only because EGGNOG_DOWNLOAD_DATA=TRUE was set"
  "${MAMBA_BIN}" run -n "${EGGNOG_ENV}" download_eggnog_data.py \
    --data_dir "${EGGNOG_DATA_DIR}" \
    -y \
    "${force_flag[@]}"

  for file in "${required_files[@]}"; do
    check_file "${file}"
  done
}

run_eggnog_mapper() {
  local annotation_file
  local override_flag=()
  local excel_flag=()
  local manifest
  local signature
  local run_output_dir="${tmpdir}/eggnog_output.$$"
  local produced

  annotation_file="${DIR_08_EGGNOG}/${EGGNOG_PREFIX}.emapper.annotations"

  check_file "${EGAPX_PROTEINS}"
  check_executable "${MAMBA_BIN}"
  check_dir "${EGGNOG_DATA_DIR}"
  ensure_eggnog_data

  manifest="$(stage_manifest_path eggnog global)"
  signature="$(make_signature \
    "$(file_identity "${EGAPX_PROTEINS}")" "$(directory_identity "${EGGNOG_DATA_DIR}")" \
    "${EGGNOG_ENV}" "${EGGNOG_METHOD}" "${EGGNOG_ITYPE}" "${EGGNOG_TAX_SCOPE}" \
    "${EGGNOG_GO_EVIDENCE}" "${EGGNOG_SENSMODE}" "${EGGNOG_PFAM_REALIGN}" \
    "${EGGNOG_EXCEL}" "${ANNOTATION_THREADS_PER_TASK}")"

  if manifest_matches "${manifest}" "${signature}" "${annotation_file}" && ! is_true "${EGGNOG_OVERRIDE}"; then
    log "Skipping eggNOG-mapper; validated output and matching manifest exist: ${annotation_file}"
    return 0
  fi

  if is_true "${EGGNOG_OVERRIDE}"; then
    override_flag+=(--override)
  fi

  if is_true "${EGGNOG_EXCEL}"; then
    excel_flag+=(--excel)
  fi

  log "Running eggNOG-mapper"
  log "CPUs used by eggNOG-mapper: ${ANNOTATION_THREADS_PER_TASK}"
  rm -rf "${run_output_dir}"
  mkdir -p "${run_output_dir}"

  "${MAMBA_BIN}" run -n "${EGGNOG_ENV}" emapper.py \
    -i "${EGAPX_PROTEINS}" \
    --itype "${EGGNOG_ITYPE}" \
    -m "${EGGNOG_METHOD}" \
    --data_dir "${EGGNOG_DATA_DIR}" \
    --cpu "${ANNOTATION_THREADS_PER_TASK}" \
    --output "${EGGNOG_PREFIX}" \
    --output_dir "${run_output_dir}" \
    --scratch_dir "${EGGNOG_SCRATCH_DIR}" \
    --temp_dir "${EGGNOG_TMP_DIR}" \
    --tax_scope "${EGGNOG_TAX_SCOPE}" \
    --go_evidence "${EGGNOG_GO_EVIDENCE}" \
    --sensmode "${EGGNOG_SENSMODE}" \
    --pfam_realign "${EGGNOG_PFAM_REALIGN}" \
    --report_orthologs \
    --md5 \
    "${excel_flag[@]}" \
    "${override_flag[@]}"

  check_file "${run_output_dir}/${EGGNOG_PREFIX}.emapper.annotations"
  find "${DIR_08_EGGNOG}" -maxdepth 1 -type f -name "${EGGNOG_PREFIX}*" -delete
  while IFS= read -r -d '' produced; do
    mv -f "${produced}" "${DIR_08_EGGNOG}/$(basename "${produced}")"
  done < <(find "${run_output_dir}" -maxdepth 1 -type f -name "${EGGNOG_PREFIX}*" -print0)
  check_file "${annotation_file}"
  write_manifest "${manifest}" "${signature}" eggnog "${annotation_file}"
  log "Finished eggNOG-mapper"
}

run_interproscan() {
  local output_base
  local interpro_opts=()
  local expected_outputs=()
  local format
  local manifest
  local signature

  output_base="${DIR_09_INTERPRO}/${INTERPROSCAN_PREFIX}"

  check_file "${INTERPROSCAN_INPUT}"
  check_file "${INTERPROSCAN_SIF}"
  check_interproscan_database

  [[ -n "${SINGULARITY_BIN}" ]] || die "Neither singularity nor apptainer was found in PATH"

  IFS=',' read -r -a interpro_formats <<< "${INTERPROSCAN_FORMATS}"
  for format in "${interpro_formats[@]}"; do
    expected_outputs+=("${output_base}.${format,,}")
    run_expected_outputs+=("${run_output_dir}/${INTERPROSCAN_PREFIX}.${format,,}")
  done

  manifest="$(stage_manifest_path interproscan global)"
  signature="$(make_signature \
    "$(file_identity "${INTERPROSCAN_INPUT}")" "$(file_identity "${INTERPROSCAN_SIF}")" \
    "$(directory_identity "${INTERPROSCAN_DATA_DIR}")" "${INTERPROSCAN_FORMATS}" \
    "${INTERPROSCAN_SEQTYPE}" "${INTERPROSCAN_APPLICATIONS}" "${INTERPROSCAN_DISABLE_PRECALC}" \
    "${INTERPROSCAN_GO_TERMS}" "${INTERPROSCAN_PATHWAYS}" "${INTERPROSCAN_IPRLOOKUP}" \
    "${INTERPROSCAN_HIGHMEM}" "${ANNOTATION_THREADS_PER_TASK}")"

  if manifest_matches "${manifest}" "${signature}" "${expected_outputs[@]}" && ! is_true "${INTERPROSCAN_OVERRIDE}"; then
    log "Skipping InterProScan; validated outputs and matching manifest exist"
    return 0
  fi

  rm -rf "${run_output_dir}"
  mkdir -p "${run_output_dir}"

  interpro_opts+=(
    --input "${INTERPROSCAN_CONTAINER_INPUT}"
    --seqtype "${INTERPROSCAN_SEQTYPE}"
    --formats "${INTERPROSCAN_FORMATS}"
    --output-file-base "${INTERPROSCAN_CONTAINER_OUTPUT}/${INTERPROSCAN_PREFIX}"
    --tempdir "${INTERPROSCAN_CONTAINER_TEMP}"
    --cpu "${ANNOTATION_THREADS_PER_TASK}"
  )

  is_true "${INTERPROSCAN_DISABLE_PRECALC}" && interpro_opts+=(--disable-precalc)
  is_true "${INTERPROSCAN_GO_TERMS}" && interpro_opts+=(--goterms)
  is_true "${INTERPROSCAN_PATHWAYS}" && interpro_opts+=(--pathways)
  is_true "${INTERPROSCAN_IPRLOOKUP}" && interpro_opts+=(--iprlookup)
  is_true "${INTERPROSCAN_HIGHMEM}" && interpro_opts+=(--highmem)
  is_true "${INTERPROSCAN_VERBOSE}" && interpro_opts+=(--verbose)

  if [[ -n "${INTERPROSCAN_APPLICATIONS}" ]]; then
    interpro_opts+=(--applications "${INTERPROSCAN_APPLICATIONS}")
  fi

  log "Running InterProScan"
  log "CPUs used by InterProScan: ${ANNOTATION_THREADS_PER_TASK}"

  "${SINGULARITY_BIN}" exec \
    -B "${INTERPROSCAN_DATA_DIR}:${INTERPROSCAN_CONTAINER_DATA}" \
    -B "$(dirname "${INTERPROSCAN_INPUT}"):/input" \
    -B "${run_output_dir}:${INTERPROSCAN_CONTAINER_OUTPUT}" \
    -B "${INTERPROSCAN_TMP_DIR}:${INTERPROSCAN_CONTAINER_TEMP}" \
    "${INTERPROSCAN_SIF}" \
    "${INTERPROSCAN_CONTAINER_BIN}" \
    "${interpro_opts[@]}"

  for format in "${run_expected_outputs[@]}"; do
    check_file "${format}"
  done
  rm -f "${expected_outputs[@]}"
  for format in "${run_expected_outputs[@]}"; do
    mv -f "${format}" "${DIR_09_INTERPRO}/$(basename "${format}")"
  done
  write_manifest "${manifest}" "${signature}" interproscan "${expected_outputs[@]}"
  log "Finished InterProScan"
}

download_archive_once() {
  local url="$1"
  local target="$2"
  local force_update="$3"
  local backup_file

  if [[ -s "${target}" ]] && ! is_true "${force_update}"; then
    log "Archive already exists; skipping download: ${target}"
    return 0
  fi

  mkdir -p "$(dirname "${target}")"

  if [[ -s "${target}" ]] && is_true "${force_update}"; then
    backup_file="${target}.bak.$(date +%Y%m%d_%H%M%S)"
    mv -f "${target}" "${backup_file}"
    log "Backed up existing archive before download: ${backup_file}"
  fi

  log "Downloading archive: ${url}"
  download_with_curl "${url}" "${target}"
  check_file "${target}"
}

resume_archive_download() {
  local url="$1"
  local target="$2"

  mkdir -p "$(dirname "${target}")"
  log "Resuming archive download: ${url}"
  download_with_curl "${url}" "${target}"
  check_file "${target}"
}

check_gzip_archive() {
  local archive="$1"

  check_file "${archive}"
  require_tool gzip gzip
  gzip -t "${archive}"
}

write_swissprot_manifest() {
  {
    printf 'key\tvalue\n'
    printf 'url\t%s\n' "${SWISSPROT_URL}"
    printf 'archive\t%s\n' "${SWISSPROT_FASTA_GZ}"
    printf 'fasta\t%s\n' "${DIAMOND_DB_FASTA}"
    printf 'diamond_db\t%s\n' "${DIAMOND_DB_FILE}"
    printf 'download_enabled\t%s\n' "${RUN_SWISSPROT_DB_DOWNLOAD}"
    printf 'force_update\t%s\n' "${SWISSPROT_FORCE_DB_UPDATE}"
    printf 'created_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  } > "${SWISSPROT_MANIFEST}"
}

prepare_swissprot_fasta() {
  local backup_file
  local tmp_fasta

  if [[ -s "${DIAMOND_DB_FASTA}" ]] && ! is_true "${SWISSPROT_FORCE_DB_UPDATE}"; then
    log "Swiss-Prot FASTA already exists: ${DIAMOND_DB_FASTA}"
    return 0
  fi

  if [[ ! -s "${SWISSPROT_FASTA_GZ}" ]]; then
    if is_true "${RUN_SWISSPROT_DB_DOWNLOAD}"; then
      require_tool curl curl
      download_archive_once "${SWISSPROT_URL}" "${SWISSPROT_FASTA_GZ}" "${SWISSPROT_FORCE_DB_UPDATE}"
    else
      die "Missing Swiss-Prot FASTA archive: ${SWISSPROT_FASTA_GZ}. Download only when intended: RUN_SWISSPROT_DB_DOWNLOAD=TRUE"
    fi
  elif is_true "${RUN_SWISSPROT_DB_DOWNLOAD}" && is_true "${SWISSPROT_FORCE_DB_UPDATE}"; then
    require_tool curl curl
    download_archive_once "${SWISSPROT_URL}" "${SWISSPROT_FASTA_GZ}" "${SWISSPROT_FORCE_DB_UPDATE}"
  fi

  if ! check_gzip_archive "${SWISSPROT_FASTA_GZ}"; then
    if is_true "${RUN_SWISSPROT_DB_DOWNLOAD}"; then
      WARN "Swiss-Prot archive gzip validation failed; attempting resumable download"
      resume_archive_download "${SWISSPROT_URL}" "${SWISSPROT_FASTA_GZ}"
      check_gzip_archive "${SWISSPROT_FASTA_GZ}" \
        || die "Swiss-Prot archive gzip validation still fails after resume. Set SWISSPROT_FORCE_DB_UPDATE=TRUE to replace the archive."
    else
      die "Swiss-Prot archive gzip validation failed. Redownload only when intended: RUN_SWISSPROT_DB_DOWNLOAD=TRUE"
    fi
  fi

  mkdir -p "$(dirname "${DIAMOND_DB_FASTA}")"
  tmp_fasta="${tmpdir}/$(basename "${DIAMOND_DB_FASTA}").$$"
  rm -f "${tmp_fasta}"

  log "Extracting Swiss-Prot FASTA from archive"
  gzip -dc "${SWISSPROT_FASTA_GZ}" > "${tmp_fasta}"
  check_file "${tmp_fasta}"

  if [[ -s "${DIAMOND_DB_FASTA}" ]]; then
    backup_file="${DIAMOND_DB_FASTA}.bak.$(date +%Y%m%d_%H%M%S)"
    mv -f "${DIAMOND_DB_FASTA}" "${backup_file}"
    log "Backed up existing Swiss-Prot FASTA before replacement: ${backup_file}"
  fi

  mv -f "${tmp_fasta}" "${DIAMOND_DB_FASTA}"
  write_swissprot_manifest
  log "Swiss-Prot FASTA ready: ${DIAMOND_DB_FASTA}"
}

check_interproscan_archive_checksum() {
  check_file "${INTERPROSCAN_DATA_ARCHIVE}"
  check_file "${INTERPROSCAN_DATA_MD5}"
  require_tool md5sum md5sum

  (
    cd "$(dirname "${INTERPROSCAN_DATA_ARCHIVE}")"
    md5sum -c "$(basename "${INTERPROSCAN_DATA_MD5}")"
  )
}

prepare_interproscan_database() {
  local backup_stamp
  local data_parent
  local extracted_data_dir
  local tmp_extract

  backup_stamp="$(date +%Y%m%d_%H%M%S)"
  data_parent="$(dirname "${INTERPROSCAN_DATA_DIR}")"

  if [[ -d "${INTERPROSCAN_DATA_DIR}" ]] && ! is_true "${INTERPROSCAN_FORCE_DB_UPDATE}"; then
    log "InterProScan data directory already exists: ${INTERPROSCAN_DATA_DIR}"
    return 0
  fi

  if [[ ! -s "${INTERPROSCAN_DATA_ARCHIVE}" || ! -s "${INTERPROSCAN_DATA_MD5}" ]]; then
    if is_true "${RUN_INTERPROSCAN_DB_DOWNLOAD}"; then
      require_tool curl curl
      download_archive_once "${INTERPROSCAN_DATA_URL}" "${INTERPROSCAN_DATA_ARCHIVE}" "${INTERPROSCAN_FORCE_DB_UPDATE}"
      download_archive_once "${INTERPROSCAN_DATA_MD5_URL}" "${INTERPROSCAN_DATA_MD5}" "${INTERPROSCAN_FORCE_DB_UPDATE}"
    else
      die "Missing InterProScan data archive or md5. Download only when intended: RUN_INTERPROSCAN_DB_DOWNLOAD=TRUE"
    fi
  fi

  if ! check_interproscan_archive_checksum; then
    if is_true "${RUN_INTERPROSCAN_DB_DOWNLOAD}"; then
      WARN "InterProScan archive checksum failed; attempting resumable download"
      resume_archive_download "${INTERPROSCAN_DATA_URL}" "${INTERPROSCAN_DATA_ARCHIVE}"
      check_interproscan_archive_checksum \
        || die "InterProScan archive checksum still fails after resume. Set INTERPROSCAN_FORCE_DB_UPDATE=TRUE to replace the archive."
    else
      die "InterProScan archive checksum failed. Redownload only when intended: RUN_INTERPROSCAN_DB_DOWNLOAD=TRUE"
    fi
  fi

  if [[ -d "${INTERPROSCAN_DATA_DIR}" ]] && is_true "${INTERPROSCAN_FORCE_DB_UPDATE}"; then
    mv -f "${INTERPROSCAN_DATA_DIR}" "${INTERPROSCAN_DATA_DIR}.bak.${backup_stamp}"
    log "Backed up existing InterProScan data directory before replacement"
  fi

  tmp_extract="${tmpdir}/interproscan_data_extract.$$"
  rm -rf "${tmp_extract}"
  mkdir -p "${tmp_extract}"

  log "Extracting InterProScan data archive"
  tar -xzf "${INTERPROSCAN_DATA_ARCHIVE}" -C "${tmp_extract}"

  extracted_data_dir="$(find "${tmp_extract}" -type d -name data -print -quit)"
  [[ -n "${extracted_data_dir}" ]] || die "InterProScan data archive did not contain a data directory: ${INTERPROSCAN_DATA_ARCHIVE}"

  mkdir -p "${data_parent}"
  mv -f "${extracted_data_dir}" "${INTERPROSCAN_DATA_DIR}"
  rm -rf "${tmp_extract}"

  check_dir "${INTERPROSCAN_DATA_DIR}"
  log "InterProScan data directory ready: ${INTERPROSCAN_DATA_DIR}"
}

check_interproscan_database() {
  if [[ -d "${INTERPROSCAN_DATA_DIR}" ]] && ! is_true "${INTERPROSCAN_FORCE_DB_UPDATE}"; then
    log "InterProScan data preflight passed: ${INTERPROSCAN_DATA_DIR}"
    return 0
  fi

  prepare_interproscan_database
}

prepare_diamond_database() {
  local tmp_prefix
  local tmp_dmnd
  local backup_file

  tmp_prefix="${tmpdir}/$(basename "${DIAMOND_DB_PREFIX}").tmp.$$"
  tmp_dmnd="${tmp_prefix}.dmnd"

  check_diamond_environment
  mkdir -p "${DIAMOND_DB_DIR}"

  if [[ ! -s "${DIAMOND_DB_FASTA}" ]] || is_true "${SWISSPROT_FORCE_DB_UPDATE}" || is_true "${RUN_SWISSPROT_DB_DOWNLOAD}"; then
    prepare_swissprot_fasta
  fi

  check_file "${DIAMOND_DB_FASTA}"

  if [[ -s "${DIAMOND_DB_FILE}" ]] && ! is_true "${DIAMOND_FORCE_DB_REBUILD}"; then
    log "DIAMOND database already exists: ${DIAMOND_DB_FILE}"
    return 0
  fi

  if [[ -s "${DIAMOND_DB_FILE}" ]] && is_true "${DIAMOND_FORCE_DB_REBUILD}"; then
    WARN "DIAMOND_FORCE_DB_REBUILD=TRUE; existing database will be backed up before replacement"
  fi

  rm -f "${tmp_dmnd}"

  log "Creating DIAMOND database"
  log "Input FASTA: ${DIAMOND_DB_FASTA}"
  log "Temporary DB prefix: ${tmp_prefix}"
  log "Final DB file: ${DIAMOND_DB_FILE}"

  "${MAMBA_BIN}" run -n "${DIAMOND_ENV}" diamond makedb \
    --in "${DIAMOND_DB_FASTA}" \
    --db "${tmp_prefix}" \
    --threads "${TOTAL_CPUS}"

  check_file "${tmp_dmnd}"

  if [[ -s "${DIAMOND_DB_FILE}" ]]; then
    backup_file="${DIAMOND_DB_FILE}.bak.$(date +%Y%m%d_%H%M%S)"
    mv -f "${DIAMOND_DB_FILE}" "${backup_file}"
    log "Backed up existing DIAMOND database: ${backup_file}"
  fi

  mv -f "${tmp_dmnd}" "${DIAMOND_DB_FILE}"
  log "DIAMOND database ready: ${DIAMOND_DB_FILE}"
}

prepare_kofam_database() {
  local backup_stamp
  local n_profiles
  local tmp_ko_list
  local tmp_extract

  if is_true "${RUN_KOFAM_DB_DOWNLOAD}"; then
    require_tool curl curl
    download_archive_once "${KOFAM_KO_LIST_URL}" "${KOFAM_KO_LIST_GZ}" "${KOFAM_FORCE_DB_UPDATE}"
    download_archive_once "${KOFAM_PROFILE_ARCHIVE_URL}" "${KOFAM_PROFILE_ARCHIVE}" "${KOFAM_FORCE_DB_UPDATE}"
  elif [[ ! -s "${KOFAM_KO_LIST_GZ}" || ! -s "${KOFAM_PROFILE_ARCHIVE}" ]]; then
    die "Missing KofamScan archive(s). Download only when intended: RUN_KOFAM_DB_DOWNLOAD=TRUE"
  fi

  check_file "${KOFAM_KO_LIST_GZ}"
  check_file "${KOFAM_PROFILE_ARCHIVE}"

  backup_stamp="$(date +%Y%m%d_%H%M%S)"

  if [[ -d "${KOFAM_PROFILE_DIR}" ]]; then
    n_profiles="$(find "${KOFAM_PROFILE_DIR}" -type f -name "*.hmm" | wc -l)"
  else
    n_profiles=0
  fi

  if [[ -s "${KOFAM_KO_LIST}" && "${n_profiles}" -ge 1000 ]] && ! is_true "${KOFAM_FORCE_DB_UPDATE}"; then
    log "KofamScan database already prepared: ${KOFAM_DB_DIR}"
    return 0
  fi

  if is_true "${KOFAM_FORCE_DB_UPDATE}"; then
    WARN "KOFAM_FORCE_DB_UPDATE=TRUE; existing prepared KofamScan files will be backed up before replacement"
    [[ -s "${KOFAM_KO_LIST}" ]] && mv -f "${KOFAM_KO_LIST}" "${KOFAM_KO_LIST}.bak.${backup_stamp}"
    [[ -d "${KOFAM_PROFILE_DIR}" ]] && mv -f "${KOFAM_PROFILE_DIR}" "${KOFAM_PROFILE_DIR}.bak.${backup_stamp}"
  elif [[ -d "${KOFAM_PROFILE_DIR}" && "${n_profiles}" -lt 1000 ]]; then
    die "KofamScan profiles directory exists but is incomplete. Set KOFAM_FORCE_DB_UPDATE=TRUE to replace it: ${KOFAM_PROFILE_DIR}"
  fi

  if [[ -s "${KOFAM_KO_LIST}" ]]; then
    log "KofamScan ko_list already exists: ${KOFAM_KO_LIST}"
  else
    log "Preparing KofamScan ko_list"
    tmp_ko_list="${tmpdir}/$(basename "${KOFAM_KO_LIST}").$$"
    rm -f "${tmp_ko_list}"
    gunzip -c "${KOFAM_KO_LIST_GZ}" > "${tmp_ko_list}"
    check_file "${tmp_ko_list}"
    mv -f "${tmp_ko_list}" "${KOFAM_KO_LIST}"
  fi

  if [[ -d "${KOFAM_PROFILE_DIR}" ]]; then
    log "KofamScan profiles directory already exists: ${KOFAM_PROFILE_DIR}"
  else
    log "Extracting KofamScan profiles"
    tmp_extract="${tmpdir}/kofam_extract.$$"
    rm -rf "${tmp_extract}"
    mkdir -p "${tmp_extract}"
    tar -xzf "${KOFAM_PROFILE_ARCHIVE}" -C "${tmp_extract}"

    if [[ -d "${tmp_extract}/profiles" ]]; then
      mv -f "${tmp_extract}/profiles" "${KOFAM_PROFILE_DIR}"
    else
      mkdir -p "${KOFAM_PROFILE_DIR}"
      find "${tmp_extract}" -type f -name "*.hmm" -exec mv -t "${KOFAM_PROFILE_DIR}" {} +
    fi

    rm -rf "${tmp_extract}"
  fi

  n_profiles="$(find "${KOFAM_PROFILE_DIR}" -type f -name "*.hmm" | wc -l)"
  if [[ "${n_profiles}" -lt 1000 ]]; then
    die "KofamScan profile extraction appears incomplete: ${KOFAM_PROFILE_DIR}"
  fi

  log "KofamScan database ready: ${KOFAM_DB_DIR}"
}

prepare_dbcan_database() {
  local dbcan_required_files=()
  local missing_files=()
  local file

  check_dbcan_environment

  "${MAMBA_BIN}" run -n "${DBCAN_ENV}" run_dbcan database --help >/dev/null 2>&1 \
    || die "Your dbCAN installation does not expose 'run_dbcan database'"

  mkdir -p "${DBCAN_DB_DIR}"
  validate_dbcan_options

  if find "${DBCAN_DB_DIR}" -name "*.part" | grep -q .; then
    if is_true "${DBCAN_FORCE_DB_UPDATE}"; then
      WARN "Removing partial dbCAN database files from ${DBCAN_DB_DIR}"
      find "${DBCAN_DB_DIR}" -name "*.part" -delete
    else
      die "dbCAN database contains partial download files. Set DBCAN_FORCE_DB_UPDATE=TRUE to remove and redownload: ${DBCAN_DB_DIR}"
    fi
  fi

  mapfile -t dbcan_required_files < <(required_dbcan_database_files)
  for file in "${dbcan_required_files[@]}"; do
    [[ -s "${file}" ]] || missing_files+=("${file}")
  done

  if [[ "${#missing_files[@]}" -eq 0 ]] && ! is_true "${DBCAN_FORCE_DB_UPDATE}"; then
    log "dbCAN required database files already exist: ${DBCAN_DB_DIR}"
    return 0
  fi

  if [[ "${#missing_files[@]}" -gt 0 ]]; then
    printf 'Missing dbCAN database file(s):\n' >&2
    printf '  %s\n' "${missing_files[@]}" >&2
  fi

  log "Downloading/creating dbCAN database with AWS S3"
  log "dbCAN database directory: ${DBCAN_DB_DIR}"

  "${MAMBA_BIN}" run -n "${DBCAN_ENV}" run_dbcan database \
    --db_dir "${DBCAN_DB_DIR}" \
    --aws_s3

  if find "${DBCAN_DB_DIR}" -name "*.part" | grep -q .; then
    die "dbCAN database download left partial files in ${DBCAN_DB_DIR}"
  fi

  missing_files=()
  for file in "${dbcan_required_files[@]}"; do
    [[ -s "${file}" ]] || missing_files+=("${file}")
  done

  [[ "${#missing_files[@]}" -eq 0 ]] || {
    printf 'Missing dbCAN database file(s):\n' >&2
    printf '  %s\n' "${missing_files[@]}" >&2
    die "dbCAN database command did not create all files required by DBCAN_METHODS='${DBCAN_METHODS}'"
  }

  log "dbCAN database ready: ${DBCAN_DB_DIR}"
}

check_diamond_database() {
  if [[ ! -s "${DIAMOND_DB_FASTA}" ]] || is_true "${SWISSPROT_FORCE_DB_UPDATE}" || is_true "${RUN_SWISSPROT_DB_DOWNLOAD}"; then
    prepare_swissprot_fasta
  fi

  if is_true "${SWISSPROT_FORCE_DB_UPDATE}" && [[ -s "${DIAMOND_DB_FILE}" ]] && ! is_true "${DIAMOND_FORCE_DB_REBUILD}"; then
    WARN "Swiss-Prot FASTA was updated but existing DIAMOND database was not rebuilt. Set DIAMOND_FORCE_DB_REBUILD=TRUE to rebuild ${DIAMOND_DB_FILE}."
  fi

  if [[ -s "${DIAMOND_DB_FILE}" ]] && is_true "${RUN_DIAMOND_DB_PREP}" && is_true "${DIAMOND_FORCE_DB_REBUILD}"; then
    prepare_diamond_database
    return 0
  fi

  if [[ ! -s "${DIAMOND_DB_FILE}" ]]; then
    if is_true "${RUN_DIAMOND_DB_PREP}"; then
      prepare_diamond_database
    else
      die "Missing DIAMOND database: ${DIAMOND_DB_FILE}. Build it only when intended: RUN_DIAMOND_DB_PREP=TRUE DIAMOND_DB_FASTA=${DIAMOND_DB_FASTA}"
    fi
  fi
}

check_kofam_database() {
  local n_profiles

  if [[ -s "${KOFAM_KO_LIST}" && -d "${KOFAM_PROFILE_DIR}" ]] \
    && is_true "${RUN_KOFAM_DB_PREP}" \
    && is_true "${KOFAM_FORCE_DB_UPDATE}"; then
    prepare_kofam_database
    return 0
  fi

  if [[ ! -s "${KOFAM_KO_LIST}" || ! -d "${KOFAM_PROFILE_DIR}" ]]; then
    if is_true "${RUN_KOFAM_DB_PREP}"; then
      prepare_kofam_database
    else
      die "Missing KofamScan prepared files: ${KOFAM_KO_LIST}, ${KOFAM_PROFILE_DIR}. Extract local archives only when intended: RUN_KOFAM_DB_PREP=TRUE"
    fi
  fi

  n_profiles="$(find "${KOFAM_PROFILE_DIR}" -type f -name "*.hmm" | wc -l)"
  if [[ "${n_profiles}" -lt 1000 ]]; then
    die "KofamScan profile directory appears incomplete: ${KOFAM_PROFILE_DIR}"
  fi
}

check_dbcan_database() {
  local dbcan_required_files=()
  local missing_files=()
  local file

  if [[ -d "${DBCAN_DB_DIR}" ]]; then
    :
  elif is_true "${RUN_DBCAN_DB_PREP}"; then
    mkdir -p "${DBCAN_DB_DIR}"
  else
    die "Missing dbCAN database directory: ${DBCAN_DB_DIR}. Download only when intended: RUN_DBCAN_DB_PREP=TRUE"
  fi
  validate_dbcan_options

  if find "${DBCAN_DB_DIR}" -name "*.part" | grep -q .; then
    if is_true "${RUN_DBCAN_DB_PREP}" && is_true "${DBCAN_FORCE_DB_UPDATE}"; then
      prepare_dbcan_database
    else
      die "dbCAN database contains partial download files: ${DBCAN_DB_DIR}. Redownload only when intended: RUN_DBCAN_DB_PREP=TRUE DBCAN_FORCE_DB_UPDATE=TRUE"
    fi
  fi

  mapfile -t dbcan_required_files < <(required_dbcan_database_files)

  for file in "${dbcan_required_files[@]}"; do
    [[ -s "${file}" ]] || missing_files+=("${file}")
  done

  if [[ "${#missing_files[@]}" -gt 0 ]]; then
    if is_true "${RUN_DBCAN_DB_PREP}"; then
      prepare_dbcan_database
      missing_files=()
      for file in "${dbcan_required_files[@]}"; do
        [[ -s "${file}" ]] || missing_files+=("${file}")
      done
      [[ "${#missing_files[@]}" -eq 0 ]] || {
        printf 'Missing dbCAN database file(s):\n' >&2
        printf '  %s\n' "${missing_files[@]}" >&2
        die "dbCAN database preparation did not create all files required by DBCAN_METHODS='${DBCAN_METHODS}'"
      }
    else
      printf 'Missing dbCAN database file(s):\n' >&2
      printf '  %s\n' "${missing_files[@]}" >&2
      die "dbCAN database is incomplete for DBCAN_METHODS='${DBCAN_METHODS}'. Download/update only when intended: RUN_DBCAN_DB_PREP=TRUE"
    fi
  fi

  log "dbCAN database preflight passed: ${DBCAN_DB_DIR}"
}

run_diamond_annotation() {
  local diamond_sensitivity="${DIAMOND_SENSITIVITY#--}"
  local sensitivity_opts=()
  local diamond_advanced_opts=()
  local manifest
  local signature
  local run_output="${tmpdir}/${DIAMOND_PREFIX}.diamond.tsv.$$"

  check_file "${EGAPX_PROTEINS}"
  check_diamond_environment
  check_diamond_database

  manifest="$(stage_manifest_path diamond global)"
  signature="$(make_signature \
    "$(file_identity "${EGAPX_PROTEINS}")" "$(file_identity "${DIAMOND_DB_FILE}")" \
    "${DIAMOND_ENV}" "${DIAMOND_EVALUE}" "${DIAMOND_MAX_TARGET_SEQS}" \
    "${DIAMOND_MIN_IDENTITY}" "${DIAMOND_QUERY_COVER}" "${DIAMOND_SUBJECT_COVER}" \
    "${DIAMOND_SENSITIVITY}" "${DIAMOND_COMP_BASED_STATS}" "${DIAMOND_GAPPED_FILTER_EVALUE}" \
    "${DIAMOND_EXT_CHUNK_SIZE}" "${DIAMOND_BIN}" "${DIAMOND_EXT}" \
    "${DIAMOND_BLOCK_SIZE}" "${DIAMOND_INDEX_CHUNKS}" "${ANNOTATION_THREADS_PER_TASK}")"

  if manifest_matches "${manifest}" "${signature}" "${DIAMOND_OUT}"; then
    log "Skipping DIAMOND annotation; validated output and matching manifest exist: ${DIAMOND_OUT}"
    return 0
  fi

  log "Running DIAMOND blastp against Swiss-Prot"
  log "Query proteins: ${EGAPX_PROTEINS}"
  log "Database: ${DIAMOND_DB_FILE}"
  log "Output: ${DIAMOND_OUT}"

  case "${diamond_sensitivity}" in
    default|"")
      ;;
    fast|mid-sensitive|sensitive|more-sensitive|very-sensitive|ultra-sensitive)
      sensitivity_opts+=("--${diamond_sensitivity}")
      ;;
    *)
      die "Invalid DIAMOND_SENSITIVITY='${DIAMOND_SENSITIVITY}'. Use default, fast, mid-sensitive, sensitive, more-sensitive, very-sensitive, or ultra-sensitive."
      ;;
  esac

  if [[ -n "${DIAMOND_COMP_BASED_STATS}" ]]; then
    case "${DIAMOND_COMP_BASED_STATS}" in
      0|1|2|3|4)
        diamond_advanced_opts+=(--comp-based-stats "${DIAMOND_COMP_BASED_STATS}")
        ;;
      *)
        die "Invalid DIAMOND_COMP_BASED_STATS='${DIAMOND_COMP_BASED_STATS}'. Use 0, 1, 2, 3, or 4."
        ;;
    esac
  fi

  if [[ -n "${DIAMOND_GAPPED_FILTER_EVALUE}" ]]; then
    diamond_advanced_opts+=(--gapped-filter-evalue "${DIAMOND_GAPPED_FILTER_EVALUE}")
  fi

  if [[ -n "${DIAMOND_EXT_CHUNK_SIZE}" ]]; then
    [[ "${DIAMOND_EXT_CHUNK_SIZE}" =~ ^[0-9]+$ ]] \
      || die "Invalid DIAMOND_EXT_CHUNK_SIZE='${DIAMOND_EXT_CHUNK_SIZE}'. Use a positive integer."
    diamond_advanced_opts+=(--ext-chunk-size "${DIAMOND_EXT_CHUNK_SIZE}")
  fi

  if [[ -n "${DIAMOND_BIN}" ]]; then
    [[ "${DIAMOND_BIN}" =~ ^[0-9]+$ ]] \
      || die "Invalid DIAMOND_BIN='${DIAMOND_BIN}'. Use a positive integer."
    diamond_advanced_opts+=(--bin "${DIAMOND_BIN}")
  fi

  if [[ -n "${DIAMOND_EXT}" ]]; then
    case "${DIAMOND_EXT}" in
      banded-fast|banded-slow|full)
        diamond_advanced_opts+=(--ext "${DIAMOND_EXT}")
        ;;
      *)
        die "Invalid DIAMOND_EXT='${DIAMOND_EXT}'. Use banded-fast, banded-slow, or full."
        ;;
    esac
  fi

  rm -f "${run_output}"
  "${MAMBA_BIN}" run -n "${DIAMOND_ENV}" diamond blastp \
    --query "${EGAPX_PROTEINS}" \
    --db "${DIAMOND_DB_PREFIX}" \
    --out "${run_output}" \
    --outfmt 6 \
      qseqid qlen sseqid slen pident length mismatch gapopen \
      qstart qend sstart send evalue bitscore qcovhsp scovhsp stitle \
    --header simple \
    --evalue "${DIAMOND_EVALUE}" \
    --max-target-seqs "${DIAMOND_MAX_TARGET_SEQS}" \
    --id "${DIAMOND_MIN_IDENTITY}" \
    --query-cover "${DIAMOND_QUERY_COVER}" \
    --subject-cover "${DIAMOND_SUBJECT_COVER}" \
    "${sensitivity_opts[@]}" \
    "${diamond_advanced_opts[@]}" \
    --threads "${ANNOTATION_THREADS_PER_TASK}" \
    --block-size "${DIAMOND_BLOCK_SIZE}" \
    --index-chunks "${DIAMOND_INDEX_CHUNKS}" \
    --tmpdir "${DIAMOND_TMP_DIR}"

  check_file "${run_output}"
  mv -f "${run_output}" "${DIAMOND_OUT}"
  write_manifest "${manifest}" "${signature}" diamond "${DIAMOND_OUT}"
  log "Finished DIAMOND annotation"
}

run_kofamscan() {
  local manifest
  local signature
  local rerun_detail="FALSE"

  check_file "${EGAPX_PROTEINS}"
  check_kofamscan_environment
  check_kofam_database

  mkdir -p "${KOFAM_REANNOTATION_CACHE_DIR}"

  manifest="$(stage_manifest_path kofamscan global)"
  signature="$(make_signature \
    "$(file_identity "${EGAPX_PROTEINS}")" "$(file_identity "${KOFAM_KO_LIST}")" \
    "$(directory_identity "${KOFAM_PROFILE_DIR}")" "${KOFAM_ENV}" \
    "${ANNOTATION_THREADS_PER_TASK}")"

  if manifest_matches "${manifest}" "${signature}" "${KOFAM_OUT_DETAIL}" "${KOFAM_OUT_MAPPER}"; then
    log "Skipping KofamScan; validated detail and mapper outputs with matching manifest exist"
    return 0
  fi

  if [[ "$(awk -F '\t' '$1 == "signature" { print $2; exit }' "${manifest}" 2>/dev/null || true)" != "${signature}" ]] \
    || [[ ! -s "${KOFAM_OUT_DETAIL}" || ! -d "${KOFAM_REANNOTATION_CACHE_DIR}" ]] \
    || ! find "${KOFAM_REANNOTATION_CACHE_DIR}" -type f -print -quit | grep -q .; then
    rerun_detail="TRUE"
  fi

  if is_true "${rerun_detail}"; then
    log "Running KofamScan detail-tsv annotation"
    rm -rf "${KOFAM_REANNOTATION_CACHE_DIR}"
    mkdir -p "${KOFAM_REANNOTATION_CACHE_DIR}"
    rm -f "${KOFAM_OUT_DETAIL}" "${KOFAM_OUT_MAPPER}"

    "${MAMBA_BIN}" run -n "${KOFAM_ENV}" exec_annotation \
      -o "${KOFAM_OUT_DETAIL}" \
      -p "${KOFAM_PROFILE_DIR}" \
      -k "${KOFAM_KO_LIST}" \
      --cpu "${ANNOTATION_THREADS_PER_TASK}" \
      --tmp-dir "${KOFAM_REANNOTATION_CACHE_DIR}" \
      -f detail-tsv \
      "${EGAPX_PROTEINS}"
  fi

  check_file "${KOFAM_OUT_DETAIL}"

  log "Generating KofamScan mapper output from cached HMM results"
  rm -f "${KOFAM_OUT_MAPPER}"

  "${MAMBA_BIN}" run -n "${KOFAM_ENV}" exec_annotation \
    -o "${KOFAM_OUT_MAPPER}" \
    -p "${KOFAM_PROFILE_DIR}" \
    -k "${KOFAM_KO_LIST}" \
    --cpu "${ANNOTATION_THREADS_PER_TASK}" \
    --tmp-dir "${KOFAM_REANNOTATION_CACHE_DIR}" \
    --reannotation \
    -f mapper \
    "${EGAPX_PROTEINS}"

  check_file "${KOFAM_OUT_MAPPER}"
  write_manifest "${manifest}" "${signature}" kofamscan "${KOFAM_OUT_DETAIL}" "${KOFAM_OUT_MAPPER}"
  log "Finished KofamScan"
}

run_dbcan_annotation() {
  local done_file="${DBCAN_OUT_DIR}/.dbcan.done"
  local large_flag=()
  local memory_flag=()
  local deeptmhmm_flag=()
  local signalp_flag=()
  local manifest
  local signature
  local run_out_dir="${tmpdir}/dbcan.${DBCAN_PREFIX}.$$"
  local stale_out_dir

  check_file "${EGAPX_PROTEINS}"
  check_dbcan_annotation_command
  validate_dbcan_options
  check_dbcan_database

  manifest="$(stage_manifest_path dbcan global)"
  signature="$(make_signature \
    "$(file_identity "${EGAPX_PROTEINS}")" "$(directory_identity "${DBCAN_DB_DIR}")" \
    "${DBCAN_ENV}" "${DBCAN_MODE}" "${DBCAN_METHODS}" "${DBCAN_DIAMOND_EVALUE}" \
    "${DBCAN_HMM_EVALUE}" "${DBCAN_HMM_COVERAGE}" "${DBCAN_DBSUB_EVALUE}" \
    "${DBCAN_DBSUB_COVERAGE}" "${DBCAN_LARGE}" "${DBCAN_RUN_DEEPTMHMM}" \
    "${DBCAN_RUN_SIGNALP}" "${DBCAN_SIGNALP_ORG}" "${DBCAN_MEMORY_MONITORING}" \
    "${DBCAN_MAX_RETRIES}" "${DBCAN_MEMORY_SAFETY_FACTOR}" "${DBCAN_MAX_MEMORY_USAGE}" \
    "${DBCAN_CSV_BUFFER_SIZE}" "${ANNOTATION_THREADS_PER_TASK}")"

  if [[ -f "${done_file}" ]] && dbcan_outputs_valid \
    && [[ "$(awk -F '\t' '$1 == "signature" { print $2; exit }' "${manifest}" 2>/dev/null || true)" == "${signature}" ]]; then
    log "Skipping dbCAN; validated outputs and matching manifest exist"
    return 0
  fi

  if [[ -f "${done_file}" ]]; then
    WARN "Ignoring stale dbCAN completion marker because no dbCAN outputs were detected: ${done_file}"
    rm -f "${done_file}"
  fi

  if is_true "${DBCAN_LARGE}"; then
    large_flag+=(--large)
  else
    large_flag+=(--no-large)
  fi

  if is_true "${DBCAN_MEMORY_MONITORING}"; then
    memory_flag+=(--enable_memory_monitoring)
  else
    memory_flag+=(--no-enable_memory_monitoring)
  fi

  if is_true "${DBCAN_RUN_DEEPTMHMM}"; then
    [[ -n "${DBCAN_DEEPTMHMM_DIR}" ]] || die "DBCAN_RUN_DEEPTMHMM=TRUE but DBCAN_DEEPTMHMM_DIR is unset"
    check_dir "${DBCAN_DEEPTMHMM_DIR}"
    deeptmhmm_flag+=(
      --run_deeptmhmm
      --deeptmhmm_dir "${DBCAN_DEEPTMHMM_DIR}"
      --deeptmhmm_python "${DBCAN_DEEPTMHMM_PYTHON}"
    )
  else
    deeptmhmm_flag+=(--no-run_deeptmhmm)
  fi

  if is_true "${DBCAN_RUN_SIGNALP}"; then
    signalp_flag+=(
      --run_signalp
      --signalp_org "${DBCAN_SIGNALP_ORG}"
    )
  else
    signalp_flag+=(--no-run_signalp)
  fi

  log "Running dbCAN CAZyme annotation"
  log "Input proteins: ${EGAPX_PROTEINS}"
  log "dbCAN database: ${DBCAN_DB_DIR}"
  log "Output directory: ${DBCAN_OUT_DIR}"
  log "dbCAN methods: ${DBCAN_METHODS}"
  rm -rf "${run_out_dir}"
  mkdir -p "${run_out_dir}"

  "${MAMBA_BIN}" run -n "${DBCAN_ENV}" run_dbcan CAZyme_annotation \
    --log-file "${DBCAN_LOG_FILE}" \
    --log-level "${DBCAN_LOG_LEVEL}" \
    --mode "${DBCAN_MODE}" \
    --output_dir "${run_out_dir}" \
    --input_raw_data "${EGAPX_PROTEINS}" \
    --db_dir "${DBCAN_DB_DIR}" \
    --methods "${DBCAN_METHODS}" \
    --threads "${ANNOTATION_THREADS_PER_TASK}" \
    --e_value_threshold "${DBCAN_DIAMOND_EVALUE}" \
    --coverage_threshold_dbcan "${DBCAN_HMM_COVERAGE}" \
    --e_value_threshold_dbcan "${DBCAN_HMM_EVALUE}" \
    --coverage_threshold_dbsub "${DBCAN_DBSUB_COVERAGE}" \
    --e_value_threshold_dbsub "${DBCAN_DBSUB_EVALUE}" \
    --max_retries "${DBCAN_MAX_RETRIES}" \
    --memory_safety_factor "${DBCAN_MEMORY_SAFETY_FACTOR}" \
    --max_memory_usage "${DBCAN_MAX_MEMORY_USAGE}" \
    --csv_buffer_size "${DBCAN_CSV_BUFFER_SIZE}" \
    "${large_flag[@]}" \
    "${memory_flag[@]}" \
    "${deeptmhmm_flag[@]}" \
    "${signalp_flag[@]}"

  dbcan_dir_outputs_valid "${run_out_dir}" \
    || die "dbCAN completed without detectable output files in ${run_out_dir}"

  touch "${run_out_dir}/.dbcan.done"
  if dbcan_outputs_valid; then
    stale_out_dir="${DBCAN_OUT_DIR}.stale.$(date +%Y%m%d_%H%M%S)"
    mv -f "${DBCAN_OUT_DIR}" "${stale_out_dir}"
    WARN "Preserved superseded dbCAN output directory: ${stale_out_dir}"
  else
    rm -rf "${DBCAN_OUT_DIR}"
  fi
  mv -f "${run_out_dir}" "${DBCAN_OUT_DIR}"
  write_manifest "${manifest}" "${signature}" dbcan "${done_file}"
  log "Finished dbCAN"
}

# ==============================================================================
# 4) R plotting script
# ==============================================================================

create_rnaseq_r_script() {
  check_file "${RNASEQ_R_SOURCE}"
  mkdir -p "$(dirname "${R_SCRIPT}")"
  cp -f "${RNASEQ_R_SOURCE}" "${R_SCRIPT}"
  chmod +x "${R_SCRIPT}"
  log "Installed R script: ${R_SCRIPT}"
}


create_genome_annotation_graphics_r_script() {
  check_file "${GENOME_GRAPHICS_R_SOURCE}"
  mkdir -p "$(dirname "${GENOME_GRAPHICS_R_SCRIPT}")"
  cp -f "${GENOME_GRAPHICS_R_SOURCE}" "${GENOME_GRAPHICS_R_SCRIPT}"
  chmod +x "${GENOME_GRAPHICS_R_SCRIPT}"
  log "Installed R script: ${GENOME_GRAPHICS_R_SCRIPT}"
}

calculate_genome_plot_dimensions() {
  local n="${GENOME_PLOT_TOP_N}"

  [[ "${n}" =~ ^[0-9]+$ ]] || die "GENOME_PLOT_TOP_N must be an integer: ${n}"
  (( n >= 1 )) || die "GENOME_PLOT_TOP_N must be >= 1"

  # These formulas reproduce the previous good default for 15 contigs
  # (karyoploteR 36 x 32 in) and scale up/down with the number of plotted sequences.
  if [[ -z "${GENOME_PLOT_WIDTH}" ]]; then
    GENOME_PLOT_WIDTH="$(awk -v n="${n}" 'BEGIN{x=12+n*0.4; if(x<18)x=18; printf "%.1f", x}')"
  fi
  if [[ -z "${GENOME_PLOT_HEIGHT}" ]]; then
    GENOME_PLOT_HEIGHT="$(awk -v n="${n}" 'BEGIN{x=8+n*0.27; if(x<12)x=12; printf "%.1f", x}')"
  fi
  if [[ -z "${GENOME_LINEAR_WIDTH}" ]]; then
    GENOME_LINEAR_WIDTH="$(awk -v n="${n}" 'BEGIN{x=14+n*0.3; if(x<18)x=18; printf "%.1f", x}')"
  fi
  if [[ -z "${GENOME_LINEAR_HEIGHT}" ]]; then
    GENOME_LINEAR_HEIGHT="$(awk -v n="${n}" 'BEGIN{x=8+n*0.4; if(x<14)x=14; printf "%.1f", x}')"
  fi
  if [[ -z "${GENOME_KARYO_WIDTH}" ]]; then
    GENOME_KARYO_WIDTH="$(awk -v n="${n}" 'BEGIN{x=18+n*1.2; if(x<24)x=24; printf "%.1f", x}')"
  fi
  if [[ -z "${GENOME_KARYO_HEIGHT}" ]]; then
    GENOME_KARYO_HEIGHT="$(awk -v n="${n}" 'BEGIN{x=14+n*1.2; if(x<20)x=20; printf "%.1f", x}')"
  fi

  log "Genome graphics dimensions for top_n=${n}: circular=${GENOME_PLOT_WIDTH}x${GENOME_PLOT_HEIGHT}, linear=${GENOME_LINEAR_WIDTH}x${GENOME_LINEAR_HEIGHT}, karyoploteR=${GENOME_KARYO_WIDTH}x${GENOME_KARYO_HEIGHT}"
}

run_genome_annotation_graphics() {
  local eggnog_file="${DIR_08_EGGNOG}/${EGGNOG_PREFIX}.emapper.annotations"
  local interpro_file="${DIR_09_INTERPRO}/${INTERPROSCAN_PREFIX}.tsv"
  local diamond_file="${DIAMOND_OUT}"
  local kofam_file="${KOFAM_OUT_DETAIL}"
  local dbcan_dir="${DBCAN_OUT_DIR}"
  local feature_counts_file="${EGAPX_FEATURE_COUNTS}"
  local id_map_file="${GENOME_ID_MAP}"
  local format_flags=()
  local genome_format
  local manifest
  local signature
  local primary_format="${GENOME_PLOT_FORMATS%%,*}"
  local primary_output="${GENOME_GRAPHICS_OUT_PREFIX}.${primary_format}"

  check_r_environment
  check_file "${GENOME_GRAPHICS_R_SCRIPT}"
  check_file "${EGAPX_GFF}"
  check_file "${EGAPX_GENOME}"
  mkdir -p "${DIR_13_PLOTS}" "${DIR_13_TABLES}" "${DIR_13_ANNOTATION}"

  [[ -s "${feature_counts_file}" ]] || feature_counts_file=""
  [[ -s "${eggnog_file}" ]] || eggnog_file=""
  [[ -s "${interpro_file}" ]] || interpro_file=""
  [[ -s "${diamond_file}" ]] || diamond_file=""
  [[ -s "${kofam_file}" ]] || kofam_file=""
  [[ -d "${dbcan_dir}" ]] || dbcan_dir=""

  calculate_genome_plot_dimensions

  for genome_format in png pdf svg; do
    if [[ ",${GENOME_PLOT_FORMATS,,}," != *",${genome_format},"* ]]; then
      format_flags+=("--no_${genome_format}")
    fi
  done

  manifest="$(stage_manifest_path genome_graphics global)"
  signature="$(make_signature \
    "$(file_identity "${EGAPX_GFF}")" "$(file_identity "${EGAPX_GENOME}")" \
    "$(file_identity "${GENOME_GRAPHICS_R_SOURCE}")" "$(file_identity "${eggnog_file}")" \
    "$(file_identity "${interpro_file}")" "$(file_identity "${diamond_file}")" \
    "$(file_identity "${kofam_file}")" "$(directory_identity "${dbcan_dir}")" \
    "${GENOME_PLOT_TOP_N}" "${GENOME_PLOT_WINDOW}" "${GENOME_CDS_METRIC}" \
    "${GENOME_LABEL_STYLE}" "${GENOME_OFFICIAL_COUNT_SOURCE}" "${GENOME_PLOT_COUNT_SOURCE}" \
    "${GENOME_PLOT_COUNT_COMPARISON}" "${GENOME_PLOT_FUNCTIONAL_COVERAGE}" \
    "${GENOME_PLOT_FUNCTIONAL_SUBTRACKS}" "${GENOME_METHOD_SUMMARY_PLOTS}" \
    "${GENOME_METHOD_KARYOPLOTS}" "${GENOME_CREATE_ID_MAP}" "${GENOME_FORCE_ID_MAP}" \
    "${GENOME_PLOT_FORMATS}" "${GENOME_DPI}")"
  if manifest_matches "${manifest}" "${signature}" "${primary_output}"; then
    log "Skipping EGAPx genome graphics; matching manifest and primary output exist"
    return 0
  fi

  log "Running EGAPx genome annotation graphics"
  log "Genome graphics output prefix: ${GENOME_GRAPHICS_OUT_PREFIX}"
  log "All genome-annotation plots will be written under: ${DIR_13_PLOTS}"

  run_rscript "${GENOME_GRAPHICS_R_SCRIPT}" \
    --gff "${EGAPX_GFF}" \
    --fasta "${EGAPX_GENOME}" \
    --feature_counts "${feature_counts_file}" \
    --id_map "${id_map_file}" \
    --id_map_out "${id_map_file}" \
    --eggnog "${eggnog_file}" \
    --interpro "${interpro_file}" \
    --diamond "${diamond_file}" \
    --kofam "${kofam_file}" \
    --dbcan_dir "${dbcan_dir}" \
    --out "${GENOME_GRAPHICS_OUT_PREFIX}" \
    --top_n "${GENOME_PLOT_TOP_N}" \
    --window "${GENOME_PLOT_WINDOW}" \
    --step "${GENOME_PLOT_STEP}" \
    --cds_metric "${GENOME_CDS_METRIC}" \
    --label_style "${GENOME_LABEL_STYLE}" \
    --official_count_source "${GENOME_OFFICIAL_COUNT_SOURCE}" \
    --plot_count_source "${GENOME_PLOT_COUNT_SOURCE}" \
    --plot_count_comparison "${GENOME_PLOT_COUNT_COMPARISON}" \
    --width "${GENOME_PLOT_WIDTH}" \
    --height "${GENOME_PLOT_HEIGHT}" \
    --linear_width "${GENOME_LINEAR_WIDTH}" \
    --linear_height "${GENOME_LINEAR_HEIGHT}" \
    --karyo_width "${GENOME_KARYO_WIDTH}" \
    --karyo_height "${GENOME_KARYO_HEIGHT}" \
    --karyo_bottom_oma "${GENOME_KARYO_BOTTOM_OMA}" \
    --karyo_legend_inset "${GENOME_KARYO_LEGEND_INSET}" \
    --karyo_base_cex "${GENOME_KARYO_BASE_CEX}" \
    --karyo_label_cex "${GENOME_KARYO_LABEL_CEX}" \
    --karyo_legend_cex "${GENOME_KARYO_LEGEND_CEX}" \
    --gc_ribbon_k "${GENOME_GC_RIBBON_K}" \
    --top_marker_n "${GENOME_TOP_MARKER_N}" \
    --plot_functional_coverage "${GENOME_PLOT_FUNCTIONAL_COVERAGE}" \
    --plot_functional_subtracks "${GENOME_PLOT_FUNCTIONAL_SUBTRACKS}" \
    --method_summary_plots "${GENOME_METHOD_SUMMARY_PLOTS}" \
    --method_karyoplots "${GENOME_METHOD_KARYOPLOTS}" \
    --create_id_map "${GENOME_CREATE_ID_MAP}" \
    --force_id_map "${GENOME_FORCE_ID_MAP}" \
    --dpi "${GENOME_DPI}" \
    "${format_flags[@]}"

  # Keep the generated ID map in the annotation directory if the R script wrote it
  # beside the output prefix and the requested ID map was missing before execution.
  if [[ ! -s "${id_map_file}" && -s "${GENOME_GRAPHICS_OUT_PREFIX}.egapx_gene_transcript_protein_map.tsv" ]]; then
    cp -f "${GENOME_GRAPHICS_OUT_PREFIX}.egapx_gene_transcript_protein_map.tsv" "${id_map_file}"
    log "Copied generated ID map to: ${id_map_file}"
  fi

  # Move non-plot TSV/HTML outputs generated with the plot prefix into the table or annotation directories,
  # while leaving all PNG/PDF/SVG files under DIR_13_PLOTS as requested.
  find "${DIR_13_PLOTS}" -maxdepth 1 -type f -name "${GENOME_GRAPHICS_PREFIX}*.tsv" -print0 \
    | while IFS= read -r -d '' f; do
        case "$(basename "$f")" in
          *egapx_gene_transcript_protein_map.tsv|*protein_to_gene_id_map*|*protein_to_gene_id_lookup*)
            mv -f "$f" "${DIR_13_ANNOTATION}/$(basename "$f")"
            ;;
          *)
            mv -f "$f" "${DIR_13_TABLES}/$(basename "$f")"
            ;;
        esac
      done

  find "${DIR_13_PLOTS}" -maxdepth 1 -type f -name "${GENOME_GRAPHICS_PREFIX}*.html" -print0 \
    | while IFS= read -r -d '' f; do
        mv -f "$f" "${DIR_13_TABLES}/$(basename "$f")"
      done

  write_manifest "${manifest}" "${signature}" genome_graphics "${primary_output}"
  log "Finished EGAPx genome annotation graphics"
}

run_consensus_annotation_outputs() {
  check_file "${COUNTS_MATRIX}"
  check_file "${SAMPLE_METADATA}"
  check_r_environment
  check_file "${R_SCRIPT}"

  local eggnog_file="${DIR_08_EGGNOG}/${EGGNOG_PREFIX}.emapper.annotations"
  local interpro_file="${DIR_09_INTERPRO}/${INTERPROSCAN_PREFIX}.tsv"
  local diamond_file="${DIAMOND_OUT}"
  local kofam_file="${KOFAM_OUT_DETAIL}"
  local dbcan_dir="${DBCAN_OUT_DIR}"
  local manifest
  local signature
  local primary_format="${PLOT_FORMATS%%,*}"
  local primary_output="${DIR_13_PLOTS}/library_size_by_organ.${primary_format}"

  if [[ ! -s "${eggnog_file}" ]]; then
    eggnog_file=""
  fi

  if [[ ! -s "${interpro_file}" ]]; then
    interpro_file=""
  fi

  if [[ ! -s "${diamond_file}" ]]; then
    diamond_file=""
  fi

  if [[ ! -s "${kofam_file}" ]]; then
    kofam_file=""
  fi

  if [[ ! -d "${dbcan_dir}" ]]; then
    dbcan_dir=""
  fi

  manifest="$(stage_manifest_path consensus_plots global)"
  signature="$(make_signature \
    "$(file_identity "${COUNTS_MATRIX}")" "$(file_identity "${SAMPLE_METADATA}")" \
    "$(file_identity "${RNASEQ_R_SOURCE}")" "$(file_identity "${eggnog_file}")" \
    "$(file_identity "${interpro_file}")" "$(file_identity "${diamond_file}")" \
    "$(file_identity "${kofam_file}")" "$(directory_identity "${dbcan_dir}")" \
    "${TOP_N_EXPRESSED}" "${PSEUDOCOUNT}" "${VOLCANO_LFC_THRESHOLD}" \
    "${DESCRIPTIVE_VOLCANO_DENOMINATOR}" "${RUN_DESCRIPTIVE_VOLCANO}" "${RUN_TIDYPLOTS}" \
    "${PLOT_FORMATS}" "${PLOT_DPI}" "${PLOT_WIDTH}" "${PLOT_HEIGHT}")"
  if manifest_matches "${manifest}" "${signature}" "${primary_output}"; then
    log "Skipping consensus plots and tables; matching manifest and primary output exist"
    return 0
  fi

  log "Running consensus annotation plots and tables"

  run_rscript "${R_SCRIPT}" \
    --counts "${COUNTS_MATRIX}" \
    --metadata "${SAMPLE_METADATA}" \
    --eggnog "${eggnog_file}" \
    --interpro "${interpro_file}" \
    --diamond "${diamond_file}" \
    --kofam "${kofam_file}" \
    --dbcan_dir "${dbcan_dir}" \
    --hisat2_summary_dir "${DIR_06_SUMMARY}" \
    --outdir "${DIR_13_CONSENSUS}" \
    --plotdir "${DIR_13_PLOTS}" \
    --tabledir "${DIR_13_TABLES}" \
    --annotation_dir "${DIR_13_ANNOTATION}" \
    --run_descriptive_volcano "${RUN_DESCRIPTIVE_VOLCANO}" \
    --run_tidyplots "${RUN_TIDYPLOTS}" \
    --denominator "${DESCRIPTIVE_VOLCANO_DENOMINATOR}" \
    --top_n "${TOP_N_EXPRESSED}" \
    --pseudocount "${PSEUDOCOUNT}" \
    --lfc_threshold "${VOLCANO_LFC_THRESHOLD}" \
    --dpi "${PLOT_DPI}" \
    --formats "${PLOT_FORMATS}" \
    --width "${PLOT_WIDTH}" \
    --height "${PLOT_HEIGHT}"

  write_manifest "${manifest}" "${signature}" consensus_plots "${primary_output}"
  log "Finished consensus annotation plots and tables"
}

run_annotation_task() {
  local task="$1"
  case "${task}" in
    eggnog)
      log "Step 05 [eggNOG]: started"
      run_eggnog_mapper
      log "Step 05 [eggNOG]: completed"
      ;;
    interproscan)
      log "Step 06 [interproscan]: started"
      run_interproscan
      log "Step 06 [interproscan]: completed"
      ;;
    diamond)
      log "Step 07 [diamond_swissprot]: started"
      run_diamond_annotation
      log "Step 07 [diamond_swissprot]: completed"
      ;;
    kofamscan)
      log "Step 08 [kofamscan]: started"
      run_kofamscan
      log "Step 08 [kofamscan]: completed"
      ;;
    dbcan)
      log "Step 09 [dbcan]: started"
      run_dbcan_annotation
      log "Step 09 [dbcan]: completed"
      ;;
    *)
      die "Unknown annotation task: ${task}"
      ;;
  esac
}

run_enabled_annotation_tasks() {
  local tasks=()
  local task
  local active=0
  local failed=0

  is_true "${RUN_EGGNOG_MAPPER}" && tasks+=(eggnog)
  is_true "${RUN_INTERPROSCAN}" && tasks+=(interproscan)
  is_true "${RUN_DIAMOND_ANNOTATION}" && tasks+=(diamond)
  is_true "${RUN_KOFAMSCAN}" && tasks+=(kofamscan)
  is_true "${RUN_DBCAN}" && tasks+=(dbcan)
  ((${#tasks[@]})) || return 0

  if (( ANNOTATION_PARALLEL_JOBS == 1 )); then
    for task in "${tasks[@]}"; do
      run_annotation_task "${task}"
    done
    return 0
  fi

  log "Running up to ${ANNOTATION_PARALLEL_JOBS} annotation methods concurrently with ${ANNOTATION_THREADS_PER_TASK} threads per method"
  for task in "${tasks[@]}"; do
    while (( active >= ANNOTATION_PARALLEL_JOBS )); do
      if ! wait -n; then
        failed=1
      fi
      active=$((active - 1))
    done
    run_annotation_task "${task}" &
    active=$((active + 1))
  done

  while (( active > 0 )); do
    if ! wait -n; then
      failed=1
    fi
    active=$((active - 1))
  done

  (( failed == 0 )) || die "One or more concurrent annotation methods failed"
}

# ==============================================================================
# 5) Export variables and functions for GNU Parallel
# ==============================================================================

export SRA_METADATA
export DIR_02_SRA DIR_03_FASTQ DIR_04_TRIMMED DIR_04_FASTP_JSON DIR_04_FASTP_HTML
export DIR_04_FASTP_FAILED DIR_04_FASTP_UNPAIRED
export DIR_06_BAM DIR_06_SUMMARY tmpdir MANIFEST_DIR
export NTHREADS TOTAL_CPUS PREFETCH_MAX_SIZE FASTERQ_MEM
export HISAT2_THREADS SAMTOOLS_SORT_THREADS SAMTOOLS_AUX_THREADS
export FASTP_BIN FASTP_QUAL FASTP_LENGTH_REQUIRED FASTP_N_BASE_LIMIT FASTP_UNQUALIFIED_PERCENT_LIMIT
export FASTP_COMPLEXITY_FILTER FASTP_COMPLEXITY_THRESHOLD FASTP_OVERREPRESENTATION FASTP_OVERREPRESENTATION_SAMPLING
export FASTP_POLYG FASTP_POLYG_MIN_LEN FASTP_POLYX FASTP_POLYX_MIN_LEN FASTP_CORRECTION FASTP_COMPRESSION
export FASTP_SAVE_FAILED FASTP_SAVE_UNPAIRED FASTP_DONT_OVERWRITE
export PREFETCH_BIN FASTERQ_DUMP_BIN PIGZ_BIN HISAT2_BIN SAMTOOLS_BIN
export RUN_FASTP HISAT2_INDEX_PREFIX HISAT2_PRESET HISAT2_EXTRA_OPTS EGAPX_GENOME
export KNOWN_SPLICESITES RNA_STRANDNESS RG_PLATFORM RG_CENTER SAMTOOLS_SORT_MEM
export PARALLEL_LOAD MEMFREE_DEFAULT

export -f _ts INFO WARN ERROR log die is_true check_file check_dir check_executable require_tool
export -f bam_index_exists bam_is_valid gzip_is_valid json_is_valid file_identity directory_identity
export -f make_signature manifest_matches write_manifest stage_manifest_path run_version_cmd
export -f get_organ_from_accession find_sra_path
export -f fastq_pair_is_valid download_one_sra convert_one_sra convert_and_compress_one_sra
export -f compress_one_fastq run_fastp_one_sample map_one_sample_hisat2

# ==============================================================================
# 6) Preflight checks
# ==============================================================================

log "Starting RNA-seq functional annotation pipeline"

if is_true "${RUN_DOWNLOAD}" || is_true "${RUN_FASTQ_DUMP}" || is_true "${RUN_COMPRESS_FASTQ}" \
  || is_true "${RUN_FASTP}" || is_true "${RUN_HISAT2_MAPPING}"; then
  require_tool parallel "${PARALLEL_BIN}"
fi

is_true "${RUN_DOWNLOAD}" && require_tool prefetch "${PREFETCH_BIN}"
is_true "${RUN_FASTQ_DUMP}" && require_tool fasterq-dump "${FASTERQ_DUMP_BIN}"
is_true "${RUN_COMPRESS_FASTQ}" && require_tool pigz "${PIGZ_BIN}"
is_true "${RUN_FASTP}" && require_tool fastp "${FASTP_BIN}"
if is_true "${RUN_MULTIQC_AFTER_FASTP}" || is_true "${RUN_MULTIQC_FINAL}"; then
  require_tool multiqc "${MULTIQC_BIN}"
fi

if is_true "${RUN_HISAT2_INDEX}" || is_true "${RUN_HISAT2_MAPPING}"; then
  check_file "${EGAPX_GENOME}"
  require_tool hisat2 "${HISAT2_BIN}"
  require_tool samtools "${SAMTOOLS_BIN}"
fi

is_true "${RUN_HISAT2_INDEX}" && require_tool hisat2-build "${HISAT2_BUILD_BIN}"

if [[ -n "${KNOWN_SPLICESITES}" ]]; then
  check_file "${KNOWN_SPLICESITES}"
fi

if is_true "${RUN_FEATURECOUNTS}"; then
  check_file "${EGAPX_GTF}"
  check_file "${FEATURECOUNTS_R_SOURCE}"
  check_r_environment
fi

if is_true "${RUN_EGGNOG_MAPPER}"; then
  check_file "${EGAPX_PROTEINS}"
  require_tool mamba "${MAMBA_BIN}"
  if [[ -d "${EGGNOG_DATA_DIR}" ]]; then
    check_dir "${EGGNOG_DATA_DIR}"
  elif is_true "${EGGNOG_DOWNLOAD_DATA}"; then
    mkdir -p "${EGGNOG_DATA_DIR}"
  else
    check_dir "${EGGNOG_DATA_DIR}"
  fi
fi

if is_true "${RUN_INTERPROSCAN}" || is_true "${RUN_INTERPROSCAN_DB_DOWNLOAD}"; then
  if is_true "${PREFLIGHT_ONLY}"; then
    check_dir "${INTERPROSCAN_DATA_DIR}"
  else
    check_interproscan_database
  fi

  if is_true "${RUN_INTERPROSCAN}"; then
    check_file "${EGAPX_PROTEINS}"
    check_file "${INTERPROSCAN_SIF}"
    [[ -n "${SINGULARITY_BIN}" ]] || die "Neither singularity nor apptainer was found in PATH"
  fi
fi

if is_true "${RUN_DIAMOND_ANNOTATION}" || is_true "${RUN_DIAMOND_DB_PREP}" || is_true "${RUN_SWISSPROT_DB_DOWNLOAD}"; then
  if is_true "${RUN_DIAMOND_ANNOTATION}"; then
    check_file "${EGAPX_PROTEINS}"
  fi
  require_tool mamba "${MAMBA_BIN}"
  check_diamond_environment
  if is_true "${PREFLIGHT_ONLY}"; then
    check_file "${DIAMOND_DB_FILE}"
  else
    check_diamond_database
  fi
fi

if is_true "${RUN_KOFAMSCAN}" || is_true "${RUN_KOFAM_DB_PREP}" || is_true "${RUN_KOFAM_DB_DOWNLOAD}"; then
  if is_true "${RUN_KOFAMSCAN}"; then
    check_file "${EGAPX_PROTEINS}"
    require_tool mamba "${MAMBA_BIN}"
    check_kofamscan_environment
  fi
  check_kofam_database
fi

if is_true "${RUN_DBCAN}" || is_true "${RUN_DBCAN_DB_PREP}"; then
  check_file "${EGAPX_PROTEINS}"
  require_tool mamba "${MAMBA_BIN}"
  check_dbcan_environment
  validate_dbcan_options
  check_dbcan_database

  if is_true "${RUN_DBCAN}"; then
    check_dbcan_annotation_command
  fi
fi

if is_true "${RUN_RNASEQ_PLOTS_TABLES}"; then
  if ! is_true "${RUN_FEATURECOUNTS}"; then
    check_file "${COUNTS_MATRIX}"
  fi
  if ! is_true "${RUN_CREATE_METADATA}"; then
    check_file "${SAMPLE_METADATA}"
  fi
  check_file "${RNASEQ_R_SOURCE}"
  if ! is_true "${RUN_CREATE_R_SCRIPT}"; then
    check_file "${R_SCRIPT}"
  fi
  check_r_environment
fi

if is_true "${RUN_GENOME_ANNOTATION_GRAPHICS}"; then
  check_file "${EGAPX_GFF}"
  check_file "${EGAPX_GENOME}"
  check_file "${GENOME_GRAPHICS_R_SOURCE}"
  if ! is_true "${RUN_CREATE_GENOME_GRAPHICS_R_SCRIPT}"; then
    check_file "${GENOME_GRAPHICS_R_SCRIPT}"
  fi
  check_r_environment
fi

log "System configuration:"
uname -a >&2

log "Tool versions:"
if [[ -x "${PARALLEL_BIN}" ]]; then
  "${PARALLEL_BIN}" --version | head -n 2 >&2 || WARN "Version check failed for parallel"
fi
if [[ -x "${PREFETCH_BIN}" ]]; then
  run_version_cmd prefetch "${PREFETCH_BIN}" --version
fi
if [[ -x "${FASTERQ_DUMP_BIN}" ]]; then
  run_version_cmd fasterq-dump "${FASTERQ_DUMP_BIN}" --version
fi
if [[ -x "${FASTP_BIN}" ]]; then
  run_version_cmd fastp "${FASTP_BIN}" --version
fi
if [[ -x "${MULTIQC_BIN}" ]]; then
  run_version_cmd multiqc "${MULTIQC_BIN}" --version
fi
if [[ -x "${HISAT2_BIN}" ]]; then
  "${HISAT2_BIN}" --version | head -n 2 >&2 || WARN "Version check failed for hisat2"
fi
if [[ -x "${SAMTOOLS_BIN}" ]]; then
  "${SAMTOOLS_BIN}" --version | head -n 2 >&2 || WARN "Version check failed for samtools"
fi
if is_true "${RUN_FEATURECOUNTS}"; then
  run_rscript -e 'cat("Rsubread version:", as.character(packageVersion("Rsubread")), "\n")' >&2 || WARN "Version check failed for Rsubread"
fi
if [[ -x "${MAMBA_BIN}" ]]; then
  run_version_cmd mamba "${MAMBA_BIN}" --version
fi
if is_true "${RUN_DIAMOND_ANNOTATION}" || is_true "${RUN_DIAMOND_DB_PREP}"; then
  "${MAMBA_BIN}" run -n "${DIAMOND_ENV}" diamond version >&2 \
    || WARN "Version check failed for DIAMOND"
fi
if is_true "${RUN_KOFAMSCAN}" || is_true "${RUN_KOFAM_DB_PREP}"; then
  "${MAMBA_BIN}" run -n "${KOFAM_ENV}" exec_annotation --help | head -n 8 >&2 \
    || WARN "Version check failed for KofamScan"
fi
if is_true "${RUN_DBCAN}" || is_true "${RUN_DBCAN_DB_PREP}"; then
  "${MAMBA_BIN}" run -n "${DBCAN_ENV}" run_dbcan --help | head -n 12 >&2 \
    || WARN "Version check failed for dbCAN"
fi
if [[ -n "${SINGULARITY_BIN}" ]]; then
  run_version_cmd singularity "${SINGULARITY_BIN}" --version
fi
if is_true "${USE_R_MAMBA_ENV}"; then
  if [[ -x "${MAMBA_BIN}" ]]; then
    "${MAMBA_BIN}" run -n "${R_ENV}" Rscript --version >&2 || WARN "Version check failed for Rscript in ${R_ENV}"
  fi
else
  if [[ -x "${RSCRIPT_BIN}" ]]; then
    run_version_cmd Rscript "${RSCRIPT_BIN}" --version
  fi
fi

log "Project directory: ${PROJECT_DIR}"
log "Temporary directory: ${tmpdir}"
log "Run mode: ${RUN_MODE}"
log "Threads per job: ${NTHREADS}"
log "Parallel jobs: ${NJOBS}"
log "Requested parallel jobs: ${REQUESTED_NJOBS}"
log "Mapping jobs: ${MAPPING_JOBS}"
log "Mapping threads per job: HISAT2=${HISAT2_THREADS}, samtools sort=${SAMTOOLS_SORT_THREADS}"
log "Total CPUs for single-process annotation/counting steps: ${TOTAL_CPUS}"
log "Annotation concurrency: jobs=${ANNOTATION_PARALLEL_JOBS}, threads per method=${ANNOTATION_THREADS_PER_TASK}"
log "Detected cores: ${NPROC_VAL}"
log "GNU Parallel load limit: ${PARALLEL_LOAD}"
log "GNU Parallel memfree limit: ${MEMFREE_DEFAULT}"
log "Runtime log file: ${LOG_FILE}"
log "Plot export formats: ${PLOT_FORMATS}; raster DPI: ${PLOT_DPI}"
log "R mamba environment enabled: ${USE_R_MAMBA_ENV}"
log "R environment name: ${R_ENV}"
log "tidyplots plotting enabled: ${RUN_TIDYPLOTS}"
log "featureCounts enabled: ${RUN_FEATURECOUNTS}"

if is_true "${PREFLIGHT_ONLY}"; then
  log "Preflight checks completed successfully; no pipeline stages were executed"
  exit 0
fi

# ==============================================================================
# 7) Pipeline execution
# ==============================================================================

if is_true "${RUN_CREATE_METADATA}" || is_true "${RUN_DOWNLOAD}" || is_true "${RUN_FASTQ_DUMP}" || is_true "${RUN_COMPRESS_FASTQ}"; then
  log "Step 01 [sra2fastq]: started"

  if is_true "${RUN_CREATE_METADATA}"; then
    log "sra2fastq: creating metadata tables"
    create_sample_metadata
  fi

  if is_true "${RUN_DOWNLOAD}"; then
    log "sra2fastq: downloading SRA files"
    "${PARALLEL_BIN}" --jobs "${NJOBS}" --line-buffer --load "${PARALLEL_LOAD}" --memfree "${MEMFREE_DEFAULT}" --tmpdir "${tmpdir}" \
      download_one_sra {} ::: "${SRA_RUNS[@]}"
  fi

  if is_true "${RUN_FASTQ_DUMP}" && is_true "${RUN_COMPRESS_FASTQ}"; then
    log "sra2fastq: converting and immediately compressing each accession"
    "${PARALLEL_BIN}" --jobs "${NJOBS}" "${PARALLEL_PROGRESS_OPTS[@]}" --line-buffer --load "${PARALLEL_LOAD}" --memfree "${MEMFREE_DEFAULT}" --tmpdir "${tmpdir}" \
      convert_and_compress_one_sra {} ::: "${SRA_RUNS[@]}"
  elif is_true "${RUN_FASTQ_DUMP}"; then
    log "sra2fastq: converting SRA to FASTQ"
    "${PARALLEL_BIN}" --jobs "${NJOBS}" "${PARALLEL_PROGRESS_OPTS[@]}" --line-buffer --load "${PARALLEL_LOAD}" --memfree "${MEMFREE_DEFAULT}" --tmpdir "${tmpdir}" \
      convert_one_sra {} ::: "${SRA_RUNS[@]}"
  fi

  if is_true "${RUN_COMPRESS_FASTQ}" && ! is_true "${RUN_FASTQ_DUMP}"; then
    log "sra2fastq: compressing FASTQ files"
    mapfile -t FASTQ_FILES < <(find "${DIR_03_FASTQ}" -maxdepth 1 -type f -name "*.fastq" | sort)

    if [[ "${#FASTQ_FILES[@]}" -eq 0 ]]; then
      log "sra2fastq: no uncompressed FASTQ files found"
    else
      "${PARALLEL_BIN}" --jobs "${NJOBS}" "${PARALLEL_PROGRESS_OPTS[@]}" --line-buffer --load "${PARALLEL_LOAD}" --memfree "${MEMFREE_DEFAULT}" --tmpdir "${tmpdir}" \
        compress_one_fastq {} ::: "${FASTQ_FILES[@]}"
    fi
  fi

  log "Step 01 [sra2fastq]: completed"
fi

if is_true "${RUN_FASTP}" || is_true "${RUN_MULTIQC_AFTER_FASTP}"; then
  log "Step 02 [fastp]: started"

  if is_true "${RUN_FASTP}"; then
    log "fastp: running QC and trimming"
    "${PARALLEL_BIN}" --jobs "${NJOBS}" "${PARALLEL_PROGRESS_OPTS[@]}" --line-buffer --load "${PARALLEL_LOAD}" --memfree "${MEMFREE_DEFAULT}" --tmpdir "${tmpdir}" \
      run_fastp_one_sample {} ::: "${SRA_RUNS[@]}"
  fi

  if is_true "${RUN_MULTIQC_AFTER_FASTP}"; then
    log "fastp: running MultiQC after fastp"
    run_multiqc_report "after_fastp"
  fi

  log "Step 02 [fastp]: completed"
fi

if is_true "${RUN_HISAT2_INDEX}" || is_true "${RUN_HISAT2_MAPPING}"; then
  log "Step 03 [hisat2]: started"

  if is_true "${RUN_HISAT2_INDEX}"; then
    log "hisat2: building/checking HISAT2 index"
    build_hisat2_index
  fi

  if is_true "${RUN_HISAT2_MAPPING}"; then
    log "hisat2: mapping RNA-seq reads"

    if ! check_hisat2_index "${HISAT2_INDEX_PREFIX}"; then
      die "HISAT2 index not found: ${HISAT2_INDEX_PREFIX}. Activate RUN_HISAT2_INDEX=\"TRUE\"."
    fi

    "${PARALLEL_BIN}" --jobs "${MAPPING_JOBS}" "${PARALLEL_PROGRESS_OPTS[@]}" --line-buffer --load "${PARALLEL_LOAD}" --memfree "${MEMFREE_DEFAULT}" --tmpdir "${tmpdir}" \
      map_one_sample_hisat2 {} ::: "${SRA_RUNS[@]}"
  fi

  log "Step 03 [hisat2]: completed"
fi

if is_true "${RUN_FEATURECOUNTS}"; then
  log "Step 04 [featurecounts]: started"
  log "featurecounts: running featureCounts"
  run_featurecounts
  log "Step 04 [featurecounts]: completed"
fi

run_enabled_annotation_tasks

if is_true "${RUN_CREATE_R_SCRIPT}" || is_true "${RUN_RNASEQ_PLOTS_TABLES}" || is_true "${RUN_CREATE_GENOME_GRAPHICS_R_SCRIPT}" || is_true "${RUN_GENOME_ANNOTATION_GRAPHICS}" || is_true "${RUN_MULTIQC_FINAL}"; then
  log "Step 10 [consensus_annotation]: started"

  if is_true "${RUN_CREATE_R_SCRIPT}"; then
    log "consensus_annotation: creating R plotting script"
    create_rnaseq_r_script
  fi

  if is_true "${RUN_RNASEQ_PLOTS_TABLES}"; then
    log "consensus_annotation: running descriptive RNA-seq plots and consensus tables"
    run_consensus_annotation_outputs
  fi

  if is_true "${RUN_CREATE_GENOME_GRAPHICS_R_SCRIPT}"; then
    log "consensus_annotation: creating EGAPx genome-annotation graphics R script"
    create_genome_annotation_graphics_r_script
  fi

  if is_true "${RUN_GENOME_ANNOTATION_GRAPHICS}"; then
    log "consensus_annotation: running EGAPx genome-annotation graphics"
    run_genome_annotation_graphics
  fi

  if is_true "${RUN_MULTIQC_FINAL}"; then
    log "consensus_annotation: running final MultiQC report"
    run_multiqc_report "final"
  fi

  log "Step 10 [consensus_annotation]: completed"
fi

# ==============================================================================
# 8) Final summary
# ==============================================================================

log "Final output summary"

log "Metadata files:"
find "${DIR_01_METADATA}" -maxdepth 1 -type f -print | sort >&2 || true

log "Raw FASTQ files:"
find "${DIR_03_FASTQ}" -maxdepth 1 -type f \( -name "*.fastq.gz" -o -name "*.fastq" \) -print | sort >&2 || true

log "Trimmed FASTQ files:"
find "${DIR_04_TRIMMED}" -maxdepth 1 -type f -name "*.fastq.gz" -print | sort >&2 || true

log "BAM files:"
find "${DIR_06_BAM}" -maxdepth 1 -type f \( -name "*.bam" -o -name "*.bam.bai" -o -name "*.bam.csi" \) -print | sort >&2 || true

log "featureCounts outputs:"
find "${DIR_07_COUNTS}" -maxdepth 1 -type f -print | sort >&2 || true

log "eggNOG outputs:"
find "${DIR_08_EGGNOG}" -maxdepth 1 -type f -name "${EGGNOG_PREFIX}*" -print | sort >&2 || true

log "InterProScan outputs:"
find "${DIR_09_INTERPRO}" -maxdepth 1 -type f -name "${INTERPROSCAN_PREFIX}*" -print | sort >&2 || true

log "DIAMOND/Swiss-Prot outputs:"
find "${DIR_10_DIAMOND}" -maxdepth 1 -type f -name "${DIAMOND_PREFIX}*" -print | sort >&2 || true

log "KofamScan outputs:"
find "${DIR_11_KOFAM}" -maxdepth 1 -type f -name "${KOFAM_PREFIX}*" -print | sort >&2 || true

log "dbCAN outputs:"
find "${DIR_12_DBCAN}" -maxdepth 3 -type f -print | sort >&2 || true

log "Consensus annotation outputs:"
find "${DIR_13_CONSENSUS}" -maxdepth 2 -type f -print | sort >&2 || true

log "Consensus annotation plots:"
find "${DIR_13_PLOTS}" -maxdepth 1 -type f \( -name "*.png" -o -name "*.pdf" -o -name "*.svg" \) -print | sort >&2 || true

log "Consensus annotation tables:"
find "${DIR_13_TABLES}" -maxdepth 1 -type f -print | sort >&2 || true

end_time=$(date +%s)
runtime=$((end_time - start_time))

cleanup_workdir_tmp_artifacts

printf "Pipeline completed successfully.\n" >&2
printf "End time: %s\n" "$(date)" >&2
printf "Total runtime: %02d:%02d:%02d (%d seconds)\n" \
  "$((runtime / 3600))" \
  "$(((runtime % 3600) / 60))" \
  "$((runtime % 60))" \
  "${runtime}" >&2
