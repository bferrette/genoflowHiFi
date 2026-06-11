#!/usr/bin/env bash

# Re-exec with bash if not already running under bash
if [ -z "${BASH_VERSION:-}" ]; then
  exec /usr/bin/env bash "$0" "$@"
fi

# Usage:
# bash genoflowHiFi_v0.1.0.sh 2>&1 | tee -a genoflowHiFi_v0.1.0.$(date +%Y%m%d_%H%M%S).log

# =============================== GLOBAL SHELL MODE ===============================
set -Eeuo pipefail
shopt -s extglob
umask 0022
export LC_ALL=C
ulimit -c 0   # production default (disable core dumps); set to 'unlimited' when debugging

# Display header with description and requirements
cat <<'EOF'
##############################################################################################################                                                                                                            #
#                            🧬 genoflowHiFi v0.1.0+ 🧬                                            #
#                                                                                                            #
#  Purpose                                                                                                   #
#    End-to-end vertebrate reference-genome assembly from PacBio long-read data. The workflow performs       #
#    primary assembly, optional decontamination, iterative consensus polishing, robust scaffold selection,   #
#    final evaluation, repeat discovery/annotation, and consolidated reporting.                              #
#                                                                                                            #
#  Main Workflow                                                                                             #
#                                                                                                            #
#    1) 🧩 Hifiasm                                                                                           #
#       Primary haplotype-aware contig assembly from PacBio HiFi reads. Produces the initial draft genome    #
#       assembly used as the foundation for all downstream refinement.                                       #
#                                                                                                            #
#    2) 🧼 FCS-Decon                                                                                         #
#       Optional NCBI FCS-Adaptor and FCS-GX screening to remove sequencing adaptors, foreign contaminants,  #
#       and taxonomically inconsistent sequence before polishing.                                            #
#                                                                                                            #
#    3) 🛠️  DeepVariant Pre-Polish                                                                           #
#       Iterative consensus polishing before scaffolding. HiFi reads are mapped back to the draft assembly,  #
#       high-confidence homozygous corrections are called, filtered, and applied across multiple rounds.     #
#                                                                                                            #
#    4) 🔍 Inspector                                                                                         #
#       Reference-free assembly evaluation and correction using long-read evidence. Detects and corrects     #
#       local assembly errors before scaffolding.                                                            #
#                                                                                                            #
#    5) 🧱 Robust Scaffolding Selection                                                                      #
#       LongStitch parameter-grid search with GenomeScope2 genome profiling, QUAST structural metrics,       #
#       optional yak screening, Merqury k-mer validation, robust multi-metric ranking, and diagnostic plots. #
#                                                                                                            #
#    6) 💎 Final Polish                                                                                      #
#       Post-scaffold DeepVariant polishing followed by optional final FCS-GX cleanup. Produces the final    #
#       polished and decontaminated genome assembly.                                                         #
#                                                                                                            #
#    7) 📊 Final Evaluation                                                                                  #
#       Comprehensive assembly assessment using QUAST, BUSCO v6, BlobToolKit, Qualimap, and MultiQC.         #
#       Generates structural, completeness, taxonomic, coverage, and mapping-quality summaries.              #
#                                                                                                            #
#    8) 🧬 Repeat Annotation                                                                                 #
#       Repeat discovery, curation, masking, and annotation using RepeatModeler2, TEtrimmer, RepeatMasker,   #
#       and Earl Grey annotation-only reporting.                                                             #
#                                                                                                            #
#    9) 📦 Final Summary                                                                                     #
#       Consolidated run report containing selected metrics, software provenance, resource usage, logs,      #
#       final deliverables, and important intermediate artifacts.                                            #
#                                                                                                            #
#  Design Principles                                                                                         #
#    ✅ Project-specific inputs, resources, and toggles are concentrated in Section 1.                       #
#    ✅ Software paths, conda environments, and containers are concentrated in Section 2.                    #
#    ✅ Checkpoints validate expected outputs before a step is considered complete.                          #
#    ✅ Existing valid outputs are reused to support safe restart after interruption.                        #
#    ✅ nthreads controls per-tool threading; njobs controls parallel grid/candidate execution.              #
#    ✅ Large intermediate files are removed only after successful checkpoint validation.                    #
#                                                                                                            #
#  Key Tool Groups                                                                                           #
#    🧩 Assembly / polishing     : hifiasm, DeepVariant, Inspector                                           #
#    🧱 Scaffolding selection    : LongStitch, GenomeScope2, QUAST, yak, Merqury                             #
#    📊 Assembly evaluation      : QUAST, BUSCO v6, BlobToolKit, Qualimap, MultiQC                           #
#    🧬 Repeat analysis          : RepeatModeler2, TEtrimmer, RepeatMasker, Earl Grey                        #
#    ⚙️  Infrastructure          : minimap2, samtools, bcftools, GATK, GNU Parallel, Apptainer/Singularity   #
#                                                                                                            #
#  Primary Deliverable                                                                                       #
#    📄 06.final_asm/<Scientific_Name>.fasta                                                                 #
#                                                                                                            #
#  Recommended Inputs                                                                                        #
#    🧬 PacBio HiFi reads       : required for primary assembly and polishing                                #
#    🧱 HiFi FASTA or CLR reads : optional for LongStitch scaffolding                                        #
#    🧪 Illumina reads          : optional but preferred for GenomeScope2, yak, and Merqury k-mer validation #
#                                                                                                            #
##############################################################################################################
EOF

# ===================================================================================================================
# 1) USER/PROJECT OVERRIDES — put ALL per-project knobs right here
#    (inputs, species/taxonomy, resources). Everything below should rarely change.
# ===================================================================================================================

# ---- Project & basic paths -------------------------------------------------------
working_dir=${working_dir:-"$PWD"}

# ---- Inputs (edit per project) ---------------------------------------------------
# PacBio HiFi reads (FASTQ/FASTA or gz)
pacbio_hifi=${pacbio_hifi:-"/full/path/to/rawdata/pacbio/hifi.fq.gz"}
# Optional: pre-converted HiFi FASTA for scaffolding (skip if not set)
pacbio_hifi_fasta=${pacbio_hifi_fasta:-"/full/path/to/rawdata/pacbio/hifi.fasta"}
# PacBio CLR reads for LongStitch scaffolding (skip if not set)
pacbio_clr=${pacbio_clr:-}

# ---- Species / Taxonomy (edit per project) --------------------------------------
TAXID=${TAXID:-13604}              # NCBI Taxonomy ID for species
SPECIES_NAME=${SPECIES_NAME:-""}    # explicit override, e.g. "Makaira nigricans"
SPECIES_HINT=${SPECIES_HINT:-""}    # fallback label if taxdump unavailable

# ---- Resources (edit per project) -----------------------------------------------
nthreads=${nthreads:-16}            # threads per heavy task
njobs=${njobs:-4}                   # parallel jobs where supported

# ---- Step toggles (edit per project) --------------------------------------------
run_hifiasm=${run_hifiasm:-"TRUE"}
run_fcs_early=${run_fcs_early:-${run_fcs:-"TRUE"}}
run_deepvariant=${run_deepvariant:-"TRUE"}
run_inspector=${run_inspector:-"TRUE"}
run_scaffolding=${run_scaffolding:-"TRUE"}
run_final_asm=${run_final_asm:-"TRUE"}
run_stats=${run_stats:-"TRUE"}
run_repeats=${run_repeats:-"TRUE"}

# ===================================================================================================================
# 2) SOFTWARE / ENVIRONMENTS / CONTAINERS
# ===================================================================================================================

# ---- Conda bootstrap -------------------------------------------------------------
CONDA_SH=${CONDA_SH:-""}
if [[ -z "$CONDA_SH" || ! -s "$CONDA_SH" ]]; then
  candidate="$HOME/anaconda3/etc/profile.d/conda.sh"
  [[ -s "$candidate" ]] && CONDA_SH="$candidate"
fi

if [[ -s "$CONDA_SH" ]]; then
  # shellcheck source=/dev/null
  source "$CONDA_SH"
else
  echo "[ERROR] Cannot source conda from $HOME anaconda3 paths" >&2
  exit 1
fi

MAMBA_SH=${MAMBA_SH:-""}
if [[ -z "$MAMBA_SH" || ! -s "$MAMBA_SH" ]]; then
  candidate="$HOME/anaconda3/etc/profile.d/mamba.sh"
  [[ -s "$candidate" ]] && MAMBA_SH="$candidate"
fi
if [[ -n "$MAMBA_SH" && -s "$MAMBA_SH" ]]; then
  # shellcheck source=/dev/null
  source "$MAMBA_SH"
fi

# ---- Environment runner ----------------------------------------------------------
ENV_RUNNER=${ENV_RUNNER:-mamba}                     # mamba | conda

# ---- Conda env names -------------------------------------------------------------
ENV_HIFIASM=${ENV_HIFIASM:-"HIFIASMenv"}
ENV_INSPECTOR=${ENV_INSPECTOR:-"INSPECTORenv"}
ENV_LONGSTITCH=${ENV_LONGSTITCH:-"LONGSTITCHenv"}
ENV_MERQURY=${ENV_MERQURY:-"MERQURYenv"}
ENV_GENOMESCOPE2=${ENV_GENOMESCOPE2:-"GENOMESCOPE2env"}
ENV_R=${ENV_R:-"Renv"}
ENV_QUALIMAP=${ENV_QUALIMAP:-"QUALIMAPenv"}
ENV_QUAST=${ENV_QUAST:-"QUASTenv"}
ENV_BUSCO=${ENV_BUSCO:-"BUSCOenv"}
ENV_TETRIMMER=${ENV_TETRIMMER:-"TETRIMMERenv"}
ENV_EARLGREY=${ENV_EARLGREY:-"EARLGREYenv"}

# ---- Tool binaries (compiled installs) ------------------------------------------
bcftools_bin=${bcftools_bin:-"$HOME/software/bcftools-1.23.1/bcftools"}
bedtools_bin=${bedtools_bin:-"$HOME/software/bedtools2/bin/bedtools"}
blastn_bin=${blastn_bin:-"$HOME/software/ncbi-blast-2.17.0+/bin/blastn"}
minimap2_bin=${minimap2_bin:-"$HOME/software/minimap2/minimap2"}
samtools_bin=${samtools_bin:-"$HOME/software/samtools-1.23.1/samtools"}

# ---- Containers (Singularity/Apptainer) -----------------------------------------
sif_blobtoolkit=${sif_blobtoolkit:-"$HOME/software/blobtoolkit.sif"}
sif_deepvariant=${sif_deepvariant:-"$HOME/software/deepvariant.sif"}
sif_fcs_adaptor=${sif_fcs_adaptor:-${ADAPTOR_SIF:-"$HOME/software/fcs-adaptor.sif"}}
sif_fcs_gx=${sif_fcs_gx:-${GX_SIF:-"$HOME/software/fcs-gx.sif"}}
sif_gatk4=${sif_gatk4:-"$HOME/software/gatk4.sif"}
sif_repeatmasker=${sif_repeatmasker:-${REPEATMASKER_SIF:-"$HOME/software/repeatmasker.sif"}}

# ---- NCBI FCS (central knobs) ----------------------------------------------------
ADAPTOR_RUNNER=${ADAPTOR_RUNNER:-"$HOME/software/fcs/dist/run_fcsadaptor.sh"}
GX_ROOT=${GX_ROOT:-"/full/path/to/ncbi/refseq/fcs-gx/gxdb"}
GX_DB=${GX_DB:-"$GX_ROOT/gxdb"}
GX_MANIFEST=${GX_MANIFEST:-"https://ncbi-fcs-gx.s3.amazonaws.com/gxdb/latest/all.manifest"}
GX_SIF=${GX_SIF:-"${sif_fcs_gx}"}
FCS_TINI_SUBREAPER=${FCS_TINI_SUBREAPER:-"TRUE"}
FCS_ADAPTOR_MIN_SEQ_LEN=${FCS_ADAPTOR_MIN_SEQ_LEN:-10}

# ---- External DBs / References ---------------------------------------------------
BLAST_DB_ROOT=${BLAST_DB_ROOT:-"/full/path/to/ncbi/refseq/core_nt/core_nt"}
BUSCO_LINEAGE=${BUSCO_LINEAGE:-"actinopterygii_odb12"}

# ---- Step 5 enhanced scaffolding and assembly selection --------------------------
GENOME_SIZE=${GENOME_SIZE:-""}              # optional; avoid hard-coded false precision
GENOME_SIZE_EST=${GENOME_SIZE_EST:-""}

LONGSTITCH_K_VALUES=${LONGSTITCH_K_VALUES:-"21 24 27 32 35 40 45 51"}
LONGSTITCH_W_VALUES=${LONGSTITCH_W_VALUES:-"50 100 150 250 350 500 750 1000"}
LONGSTITCH_ROUNDS=${LONGSTITCH_ROUNDS:-5}
LONGSTITCH_GAP_FILL="True"                  # fixed: maximize contiguity before HiFi polishing
LONGSTITCH_MODE=${LONGSTITCH_MODE:-"ntLink-arks"}  # tigmint-ntLink-arks | ntLink-arks

QUAST_MIN_CONTIG=${QUAST_MIN_CONTIG:-1}
QUAST_KMER_SIZE=${QUAST_KMER_SIZE:-51}
QUAST_USE_READS=${QUAST_USE_READS:-"FALSE"}
QUAST_SPLIT_SCAFFOLDS=${QUAST_SPLIT_SCAFFOLDS:-"FALSE"}

GS_K=${GS_K:-31}
GS_PLOIDY=${GS_PLOIDY:-2}
GS_NAME_PREFIX=${GS_NAME_PREFIX:-""}
GS_MAX_KMERCOV=${GS_MAX_KMERCOV:-""}
GS_LAMBDA=${GS_LAMBDA:-""}
GS_TRANSFORM_EXP=${GS_TRANSFORM_EXP:-1}
GS_TYPICAL_ERROR=${GS_TYPICAL_ERROR:-""}
GS_JSON_REPORT=${GS_JSON_REPORT:-"TRUE"}
GS_FITTED_HIST=${GS_FITTED_HIST:-"TRUE"}
GS_VERBOSE=${GS_VERBOSE:-"TRUE"}
GS_NO_UNIQUE_SEQUENCE=${GS_NO_UNIQUE_SEQUENCE:-"FALSE"}
GS_TOPOLOGY=${GS_TOPOLOGY:-""}
GS_INITIAL_REPETITIVENESS=${GS_INITIAL_REPETITIVENESS:-""}
GS_INITIAL_HETEROZYGOSITIES=${GS_INITIAL_HETEROZYGOSITIES:-""}
GS_TESTING=${GS_TESTING:-"FALSE"}
GS_TRUE_PARAMS=${GS_TRUE_PARAMS:-""}
GS_TRACE_FLAG=${GS_TRACE_FLAG:-"FALSE"}
GS_NUM_ROUNDS=${GS_NUM_ROUNDS:-""}
GS_START_SHIFT=${GS_START_SHIFT:-""}
GENOMESCOPE2_BIN=${GENOMESCOPE2_BIN:-""}

MERQURY_K=${MERQURY_K:-31}
MERQURY_TOP_N=${MERQURY_TOP_N:-10}
MERQURY=${MERQURY:-""}                       # optional path to Merqury repo if merqury.sh is not on PATH

YAK_K=${YAK_K:-31}
YAK_BLOOM_BITS=${YAK_BLOOM_BITS:-37}
STEP5_SAVE_SVG=${STEP5_SAVE_SVG:-"TRUE"}
STEP5_PLOT_DPI=${STEP5_PLOT_DPI:-600}

# ---- Repeat annotation knobs (edit per project) ---------------------------------
# Empty REPEATS_* path/name overrides are resolved at Step 8 after final assembly paths exist.
REPEATS_PROJECT_ROOT=${REPEATS_PROJECT_ROOT:-""}   # default: ${working_dir}/08.repeats
REPEATS_GENOME=${REPEATS_GENOME:-""}               # default: final species FASTA from Step 6
REPEATS_SPECIES=${REPEATS_SPECIES:-"1489904"}      # default: resolved species_name
REPEATS_THREADS=${REPEATS_THREADS:-$(( nthreads * njobs ))} # threads for RepeatModeler/TEtrimmer/RepeatMasker/Earl Grey

AUTO_RESUME=${AUTO_RESUME:-"TRUE"}
AUTO_ARCHIVE_UNRECOVERABLE_RM_DIRS=${AUTO_ARCHIVE_UNRECOVERABLE_RM_DIRS:-"TRUE"}

KNOWN_LIBRARY_OVERRIDE=${KNOWN_LIBRARY_OVERRIDE:-""}
REPEATMODELER_LIBRARY_OVERRIDE=${REPEATMODELER_LIBRARY_OVERRIDE:-""}
TETRIMMER_LIBRARY_OVERRIDE=${TETRIMMER_LIBRARY_OVERRIDE:-""}
FINAL_LIBRARY_OVERRIDE=${FINAL_LIBRARY_OVERRIDE:-""}

DFAM_FAMDB_DIR=${DFAM_FAMDB_DIR:-""}               # default: ${REPEATS_PROJECT_ROOT}/Libraries/famdb
DFAM_TAXON=${DFAM_TAXON:-${REPEATS_SPECIES}}
DFAM_FLAGS=${DFAM_FLAGS:---curated}
REPBASE_EMBL=${REPBASE_EMBL:-""}                   # default: ${REPEATS_PROJECT_ROOT}/Libraries/RepBase-20240529.embl
ARTEFACTS_FASTA=${ARTEFACTS_FASTA:-""}

USE_LTRSTRUCT=${USE_LTRSTRUCT:-"yes"}
REPEATMODELER_RECOVER=${REPEATMODELER_RECOVER:-"TRUE"}
REPEATMODELER_RECOVER_DIR=${REPEATMODELER_RECOVER_DIR:-""}
REPEATMODELER_SANITIZE_FASTA_HEADERS=${REPEATMODELER_SANITIZE_FASTA_HEADERS:-"TRUE"}
RM_ADDL_ROUNDS=${RM_ADDL_ROUNDS:-0}

TETRIMMER_MODE=${TETRIMMER_MODE:-"curated"}              # fast | curated | debug
TETRIMMER_PRESET=${TETRIMMER_PRESET:-"divergent"}        # conserved | divergent | empty
TETRIMMER_CONTINUE=${TETRIMMER_CONTINUE:-"TRUE"}
TETRIMMER_CONTINUE_MODE=${TETRIMMER_CONTINUE_MODE:-"auto"} # auto | always | never
TETRIMMER_ARCHIVE_INCOMPLETE=${TETRIMMER_ARCHIVE_INCOMPLETE:-"TRUE"}
TETRIMMER_CLASSIFICATION=${TETRIMMER_CLASSIFICATION:-"classify_unknown"} # classify_unknown | classify_all | none
TETRIMMER_DEDUP=${TETRIMMER_DEDUP:-"TRUE"}
TETRIMMER_HMM=${TETRIMMER_HMM:-"TRUE"}
TETRIMMER_GENOME_ANNO=${TETRIMMER_GENOME_ANNO:-"FALSE"}
TETRIMMER_DEBUG=${TETRIMMER_DEBUG:-"FALSE"}
TETRIMMER_PFAM_DIR=${TETRIMMER_PFAM_DIR:-""}
TETRIMMER_CONS_THR=${TETRIMMER_CONS_THR:-"0.7"}
TETRIMMER_MINI_ORF=${TETRIMMER_MINI_ORF:-"200"}
TETRIMMER_MAX_MSA_LINES=${TETRIMMER_MAX_MSA_LINES:-"100"}
TETRIMMER_TOP_MSA_LINES=${TETRIMMER_TOP_MSA_LINES:-"100"}
TETRIMMER_MIN_SEQ_NUM=${TETRIMMER_MIN_SEQ_NUM:-"10"}
TETRIMMER_MIN_BLAST_LEN=${TETRIMMER_MIN_BLAST_LEN:-"150"}
TETRIMMER_MAX_CLUSTER_NUM=${TETRIMMER_MAX_CLUSTER_NUM:-"5"}
TETRIMMER_EXT_THR=${TETRIMMER_EXT_THR:-"0.7"}
TETRIMMER_EXT_CHECK_WIN=${TETRIMMER_EXT_CHECK_WIN:-"150"}
TETRIMMER_EXT_STEP=${TETRIMMER_EXT_STEP:-"1000"}
TETRIMMER_MAX_EXT=${TETRIMMER_MAX_EXT:-"7000"}
TETRIMMER_GAP_THR=${TETRIMMER_GAP_THR:-"0.4"}
TETRIMMER_GAP_NUL_THR=${TETRIMMER_GAP_NUL_THR:-"0.7"}
TETRIMMER_CROP_END_DIV_THR=${TETRIMMER_CROP_END_DIV_THR:-"0.7"}
TETRIMMER_CROP_END_DIV_WIN=${TETRIMMER_CROP_END_DIV_WIN:-"40"}
TETRIMMER_CROP_END_GAP_THR=${TETRIMMER_CROP_END_GAP_THR:-"0.1"}
TETRIMMER_CROP_END_GAP_WIN=${TETRIMMER_CROP_END_GAP_WIN:-"250"}
TETRIMMER_LTR_START_PATTERNS=${TETRIMMER_LTR_START_PATTERNS:-"TG"}
TETRIMMER_LTR_END_PATTERNS=${TETRIMMER_LTR_END_PATTERNS:-"CA"}
TETRIMMER_HELITRON_START_PATTERNS=${TETRIMMER_HELITRON_START_PATTERNS:-"ATC"}
TETRIMMER_HELITRON_END_PATTERNS=${TETRIMMER_HELITRON_END_PATTERNS:-"CTAAT,CTAGT,CTGAT,CTGGT"}
TETRIMMER_POLY_PATTERNS=${TETRIMMER_POLY_PATTERNS:-"A"}
TETRIMMER_POLY_LEN=${TETRIMMER_POLY_LEN:-"10"}

KNOWN_CDHIT_ID=${KNOWN_CDHIT_ID:-"0.98"}
KNOWN_CDHIT_WORD=${KNOWN_CDHIT_WORD:-"8"}
FINAL_CDHIT_ID=${FINAL_CDHIT_ID:-"0.95"}
FINAL_CDHIT_WORD=${FINAL_CDHIT_WORD:-"8"}

RM_ENGINE=${RM_ENGINE:-"rmblast"}
RM_SENSITIVITY=${RM_SENSITIVITY:-"-s"}              # "-s", "-q", "-qq", or empty
RM_XSMALL=${RM_XSMALL:-"yes"}
RM_GFF=${RM_GFF:-"yes"}
RM_ALIGN=${RM_ALIGN:-"yes"}
RM_PA=${RM_PA:-$(( REPEATS_THREADS / 4 ))}
[[ "$RM_PA" -lt 1 ]] && RM_PA=1

EARLGREY_REMOVE_SHORT=${EARLGREY_REMOVE_SHORT:-"yes"}
EARLGREY_SOFTMASK=${EARLGREY_SOFTMASK:-"yes"}
EARLGREY_HELIANO=${EARLGREY_HELIANO:-"yes"}
EARLGREY_REPEATMASKER_SPECIES=${EARLGREY_REPEATMASKER_SPECIES:-"none"}
EARLGREY_ENHANCED_PLOTS=${EARLGREY_ENHANCED_PLOTS:-"TRUE"}
EARLGREY_PLOT_DPI=${EARLGREY_PLOT_DPI:-600}
EARLGREY_PLOT_TOP_FAMILIES=${EARLGREY_PLOT_TOP_FAMILIES:-30}

# ===================================================================================================================
# ===================================================================================================================
# ---- DO NOT CHANGE ANYTHING BELOW THIS LINE -------------------------------------
# ===================================================================================================================
# ===================================================================================================================

# ===================================================================================================================
# 3) DIRECTORY LAYOUT (rarely changes) — create after config to keep things tidy
# ===================================================================================================================

declare -A dirs=(
  ["hifiasm"]="${working_dir}/01.hifiasm"
  ["fcs"]="${working_dir}/02.fcs_decontam"
  ["deepvariant"]="${working_dir}/03.deepvariant"
  ["inspector"]="${working_dir}/04.inspector"
  ["scaffolding"]="${working_dir}/05.scaffolding"
  ["final_asm"]="${working_dir}/06.final_asm"
  ["asm_stats"]="${working_dir}/07.asm_stats"
  ["quast_dir"]="${working_dir}/07.asm_stats/quast"
  ["busco_dir"]="${working_dir}/07.asm_stats/busco"
  ["qualimap_dir"]="${working_dir}/07.asm_stats/qualimap"
  ["btk_dir"]="${working_dir}/07.asm_stats/blobtoolkit"
  ["repeats"]="${working_dir}/08.repeats"
  ["summary"]="${working_dir}/09.summary"
)
for d in "${dirs[@]}"; do mkdir -p "$d"; done

# Short handles (exported)
export dir1=${dirs["hifiasm"]}
export dir2=${dirs["fcs"]}
export dir3=${dirs["deepvariant"]}
export dir4=${dirs["inspector"]}
export dir5=${dirs["scaffolding"]}
export dir6=${dirs["final_asm"]}
export dir7=${dirs["asm_stats"]}
export dir8=${dirs["summary"]}
export quast_dir=${dirs["quast_dir"]}
export busco_dir=${dirs["busco_dir"]}
export qualimap_dir=${dirs["qualimap_dir"]}
export btk_dir=${dirs["btk_dir"]}
export repeats_dir=${dirs["repeats"]}
export earlgrey_dir="${repeats_dir}/06.earlgrey_annotation_only"

# Early BlobToolKit location hint (no I/O yet)
early_btk_dir=${early_btk_dir:-"${btk_dir}"}

# Make sure our taxdump path is the canonical final one
export BLOB_TAXDUMP_DIR="${btk_dir}/taxdump"
mkdir -p "${BLOB_TAXDUMP_DIR}"

# If an alternate early path was provided and differs, link it for consistency
if [[ -n "${early_btk_dir}" && "${early_btk_dir}" != "${btk_dir}" ]]; then
  mkdir -p "${early_btk_dir}"
  ln -sfn "${BLOB_TAXDUMP_DIR}" "${early_btk_dir}/taxdump"
fi

# ===================================================================================================================
# 4) TMPDIR, LOGGING UTILS, PREFLIGHT — shared infra used by all steps
# ===================================================================================================================

# Temporary working dir for big tools. Clean up only directories created here;
# user-supplied tmpdir paths may be shared scratch locations.
if [[ -n "${tmpdir:-}" ]]; then
  mkdir -p "${tmpdir}"
  cleanup_tmpdir="FALSE"
else
  tmpdir="$(mktemp -d -p "${working_dir}" tmp.XXXXXXXXXX)"
  cleanup_tmpdir="TRUE"
fi
trap 'rv=$?; if [[ "${cleanup_tmpdir:-FALSE}" == "TRUE" && -n "${tmpdir:-}" ]]; then rm -rf "${tmpdir}"; fi; exit $rv' EXIT

gatk_spark() {
  singularity exec \
    -B "${working_dir}:${working_dir}" \
    -B "${tmpdir}:${tmpdir}" \
    --env SPARK_LOCAL_IP="${SPARK_LOCAL_IP:-127.0.0.1}" \
    --env SPARK_PUBLIC_DNS="${SPARK_PUBLIC_DNS:-localhost}" \
    "${sif_gatk4}" \
    gatk "$@"
}

# Shared log timestamp (override with LOG_TS to pin)
LOG_TS=${LOG_TS:-$(date +%Y%m%d_%H%M%S)}

# Utilities
# ---- Diagnostics (fail fast, clear) ----------------------------------------
error() {
  local msg="${1:-unknown error}"
  local code="${2:-1}"
  local src="${BASH_SOURCE[1]:-$(basename "$0")}"
  local line="${BASH_LINENO[0]:-0}"
  printf '%s [%s:%s] ERROR: %s\n' "$(date +'%F %T')" "${src}" "${line}" "${msg}" >&2
  return "${code}"
}

warn() { printf '%s [WARN] %s\n' "$(date +'%F %T')" "$*" >&2; }
info() { printf '%s [INFO] %s\n' "$(date +'%F %T')" "$*"; }
fatal() { printf '%s [FATAL] %s\n' "$(date +'%F %T')" "$*" >&2; exit 1; }
require() {
  local bin="${1:-}"
  [[ -n "${bin}" ]] || fatal "require: missing executable path"
  if [[ "${bin}" == */* ]]; then
    [[ -x "${bin}" ]] || fatal "Missing required executable: ${bin}"
  else
    command -v "${bin}" >/dev/null 2>&1 || fatal "Missing required executable: ${bin}"
  fi
}

is_true() {
  [[ "${1:-}" =~ ^(TRUE|True|true|1|YES|yes|Y|y)$ ]]
}

is_false() {
  [[ "${1:-}" =~ ^(FALSE|False|false|0|NO|no|N|n)$ ]]
}

activate_env() {
  local env="${1:-}"
  [[ -n "${env}" ]] || error "activate_env: missing conda environment name"
  "${ENV_RUNNER}" activate "${env}" || error "Failed to activate environment with ${ENV_RUNNER}: ${env}"
}

deactivate_env() {
  "${ENV_RUNNER}" deactivate || conda deactivate || true
}

env_run() {
  local env="${1:-}"
  shift || true
  [[ -n "${env}" ]] || error "env_run: missing environment name"
  "${ENV_RUNNER}" run -n "${env}" "$@"
}

done_file() {
  local f="${1:-}"
  [[ -s "${f}" ]]
}

# Print failing command, file, and line; keep set -Eeuo pipefail effective
err_trap() {
  local rc=$?
  printf "%s [ERROR] Command failed (rc=%s) at %s:%s: %s\n" \
    "$(date +"%F %T")" "$rc" "${BASH_SOURCE[1]:-${BASH_SOURCE[0]}}" "${BASH_LINENO[0]:-0}" "${BASH_COMMAND}" >&2
  exit "$rc"
}
trap err_trap ERR

# Preflight (core tools now; step-specific checked later)
for t in "$minimap2_bin" "$samtools_bin" "$bcftools_bin"; do require "$t"; done
[[ -f "$sif_gatk4"       ]] || fatal "Missing container: $sif_gatk4"
[[ -f "$sif_deepvariant" ]] || fatal "Missing container: $sif_deepvariant"

# Optional: FCS containers when early/late decontam enabled
if [[ "${run_fcs_early}" == "TRUE" || "${run_final_asm}" == "TRUE" ]]; then
  command -v singularity >/dev/null 2>&1 || fatal "'singularity' not on PATH"
  [[ -f "$sif_fcs_adaptor" ]] || fatal "Missing container (sif_fcs_adaptor): $sif_fcs_adaptor"
  [[ -f "$sif_fcs_gx"      ]] || fatal "Missing container (sif_fcs_gx): $sif_fcs_gx"
fi

if is_true "${FCS_TINI_SUBREAPER}"; then
  export TINI_SUBREAPER=1
  export SINGULARITYENV_TINI_SUBREAPER=1
  export APPTAINERENV_TINI_SUBREAPER=1
fi

# BLAST DB sanity (warn only)
if [[ ! -f "${BLAST_DB_ROOT}.nhr" && -z "$(compgen -G "${BLAST_DB_ROOT}.*.nhr" || true)" ]]; then
  echo "[WARN] Could not confirm BLAST DB shards at ${BLAST_DB_ROOT}*" >&2
fi

start_time=${start_time:-"$(date +%s)"}

# ===================================================================================================================
# 5) EARLY TAXONOMY BOOTSTRAP — download taxdump & resolve species/base now
# ===================================================================================================================
if [[ ! -s "${BLOB_TAXDUMP_DIR}/nodes.dmp" || ! -s "${BLOB_TAXDUMP_DIR}/names.dmp" ]]; then
  echo "[INFO] Downloading NCBI new_taxdump into ${BLOB_TAXDUMP_DIR}"
  command -v curl >/dev/null 2>&1 || { echo "[ERROR] curl not found on PATH"; exit 1; }
  mkdir -p "${btk_dir}" "${BLOB_TAXDUMP_DIR}"

  # Save the tarball inside the existing temporary directory
  taxdump_tar="${tmpdir}/new_taxdump.tar.gz"

  # Robust download
  set -x
  curl --fail --show-error --location \
       --retry 5 --retry-delay 3 \
       --connect-timeout 30 --max-time 900 \
       --continue-at - \
       "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz" \
       -o "${taxdump_tar}"
  set +x

  [[ -s "${taxdump_tar}" ]] || { echo "[ERROR] Downloaded file is empty: ${taxdump_tar}" >&2; exit 1; }
  tar -xzf "${taxdump_tar}" -C "${BLOB_TAXDUMP_DIR}"

  if [[ ! -s "${BLOB_TAXDUMP_DIR}/nodes.dmp" || ! -s "${BLOB_TAXDUMP_DIR}/names.dmp" ]]; then
    echo "[ERROR] Taxdump extraction incomplete; expected nodes.dmp and names.dmp in ${BLOB_TAXDUMP_DIR}" >&2
    ls -l "${BLOB_TAXDUMP_DIR}" || true
    exit 1
  fi

  rm -f "${taxdump_tar}" || true
  echo "[INFO] Taxdump ready: ${BLOB_TAXDUMP_DIR}/nodes.dmp, names.dmp"
else
  echo "[INFO] Using existing taxdump at ${BLOB_TAXDUMP_DIR}"
fi

# Extract scientific name for a TAXID from names.dmp (never trips 'set -e')
_species_from_taxdump() {
  local taxid="$1" names_file="$2"
  local sci=""
  if [[ ! -s "${names_file}" ]]; then
    echo "[WARN] names.dmp missing at ${names_file}" >&2
    printf ''
    return 0
  fi
  sci="$(
    awk -v id="${taxid}" '
      BEGIN { FS="|"; }
      {
        for (i=1; i<=NF; i++) gsub(/^[[:space:]]+|[[:space:]]+$/, "", $i)
      }
      $1 == id && $4 == "scientific name" { print $2; exit }
    ' "${names_file}" 2>/dev/null
  )" || true
  printf '%s' "${sci}"
}

# Resolve species/base (no hardcoded names; SPECIES_NAME overrides, SPECIES_HINT fallback, else taxid_XXXX)
species_name_raw="${SPECIES_NAME:-}"
if [[ -z "${species_name_raw}" ]]; then
  species_name_raw="$(_species_from_taxdump "${TAXID}" "${BLOB_TAXDUMP_DIR}/names.dmp")" || species_name_raw=""
fi
if [[ -z "${species_name_raw}" && -n "${SPECIES_HINT:-}" ]]; then
  species_name_raw="${SPECIES_HINT}"
fi
if [[ -z "${species_name_raw}" ]]; then
  species_name_raw="taxid_${TAXID}"
fi
# Filesystem-safe species identifier (spaces → underscores; strip odd chars)
species_name=$(printf '%s' "${species_name_raw}" | sed 's/ /_/g; s/[^A-Za-z0-9._-]/_/g')
export species_name SPECIES_NAME_CANON="${species_name}"

# Make species the assembly basename unless user provided an explicit base
base=${base:-"${species_name}"}
export base
echo "[INFO] Species resolved: ${species_name_raw} → ${species_name} (TAXID=${TAXID})"
echo "[INFO] Assembly base name: ${base}"
echo "[INFO] Configuration loaded. working_dir=${working_dir} base=${base} threads=${nthreads} jobs=${njobs}"

# ============================================================================================================================
# STEP 1: Hifiasm — primary assembly from PacBio HiFi reads (haplotype-resolved)
# Uses variables defined in the centralized configuration section.
# ============================================================================================================================
if [[ "${run_hifiasm}" == "TRUE" ]]; then
  echo "[INFO] Step 1: Hifiasm — starting"

  # Preflight
  [[ -s "${pacbio_hifi}" ]] || error "HiFi input not found: ${pacbio_hifi}"
  mkdir -p "${dir1}"

  # Activate env and verify binary
  activate_env "${ENV_HIFIASM}"
  command -v hifiasm >/dev/null 2>&1 || error "hifiasm binary not found in PATH (expected: hifiasm)"

  # Output base (Hifiasm writes multiple files with this prefix)
  asm_prefix="${dir1}/${base}"

  # Skip if final FASTA already present (idempotent resume)
  final_ctg_fa="${asm_prefix}.bp.p_ctg.fa"
  if [[ -s "${final_ctg_fa}" ]]; then
    echo "[INFO] Hifiasm output already exists: ${final_ctg_fa} — skipping assembly"
    deactivate_env
  else
    echo "[INFO] Running hifiasm on ${pacbio_hifi}"
    echo "       threads=$((nthreads * njobs)) tmpdir=${tmpdir} out=${asm_prefix}"

    # Core run (tuned flags preserved from prior version)
    # -u 1         : purge duplicated haplotigs (stringent); adjust if desired
    # --ctg-n 3    : aggressive contig break near potential errors (as in original)
    # Input        : single HiFi fq.gz file supported directly
    hifiasm \
      -o "${asm_prefix}" \
      -t "$((nthreads * njobs))" \
      -u 1 \
      --ctg-n 3 \
      "${pacbio_hifi}" \
      2>&1 | tee -a "${dir1}/hifiasm.log" || error "hifiasm failed"

    # Convert primary contig GFA -> FASTA (standard extraction for hifiasm)
    p_ctg_gfa="${asm_prefix}.bp.p_ctg.gfa"
    if [[ -s "${p_ctg_gfa}" ]]; then
      awk '/^S\t/{print ">"$2"\n"$3}' "${p_ctg_gfa}" > "${final_ctg_fa}" \
        || error "Failed to convert ${p_ctg_gfa} to FASTA"
      echo "[INFO] Generated FASTA: ${final_ctg_fa}"
    else
      error "Expected GFA not found: ${p_ctg_gfa}"
    fi

    # Optional: compress large intermediate GFAs to save space
    if [[ -s "${p_ctg_gfa}" ]]; then
      pigz -f --best "${p_ctg_gfa}" || echo "[WARN] pigz failed on ${p_ctg_gfa}"
    fi

    deactivate_env
  fi
fi

# ============================================================================================================================
# STEP 2: NCBI Foreign Contamination Screen (FCS) on primary raw assembly
# Uses centralized configuration (containers, DB paths, taxon settings) from the config block.
# Produces a cleaned assembly for downstream polishing, falling back to adaptor-cleaned/original input when no GX report is produced.
# Outputs:
#   - ${dir2}/adaptor/cleaned_sequences/<input_name> (preferred) or ${dir2}/adaptor/clean.fasta
#   - ${dir2}/gx/clean.fasta (GX-cleaned FASTA, or fallback link consumed by Step 3)
#   - ${dir2}/gx/contam.fasta (filtered sequences)
#   - ${dir2}/logs/* (tool logs)
#   - exports step2_cleaned_fasta (path to final cleaned FASTA)
# ============================================================================================================================
if [[ "${run_fcs_early}" == "TRUE" ]]; then
  echo "[INFO] Step 2: FCS (Adaptor + GX via Singularity) — starting"
  mkdir -p "${dir2}/adaptor" "${dir2}/gx" "${dir2}/logs"

  # --- lightweight log helpers if not already defined ---
  if ! type -t log >/dev/null 2>&1;  then log()  { echo "[INFO]  $*"; };  fi
  if ! type -t warn >/dev/null 2>&1; then warn() { echo "[WARN]  $*" >&2; }; fi
  if ! type -t fatal >/dev/null 2>&1; then fatal(){ echo "[FATAL] $*" >&2; exit 1; }; fi

  # --- Tool & input checks (containers already preflighted in config, but re-check lightly) ---
  command -v singularity >/dev/null 2>&1 || fatal "STEP 2: 'singularity' not on PATH"
  [[ -f "${sif_fcs_adaptor}" ]] || fatal "STEP 2: sif_fcs_adaptor missing: ${sif_fcs_adaptor}"
  [[ -f "${sif_fcs_gx}"      ]] || fatal "STEP 2: sif_fcs_gx missing: ${sif_fcs_gx}"

  input_fasta="${dir1}/${base}.bp.p_ctg.fa"
  [[ -s "${input_fasta}" ]] || fatal "STEP 2: Input FASTA not found: ${input_fasta} (did STEP 1 run?)"

  mkdir -p "${GX_DB}"
  [[ -w "${GX_DB}" ]] || fatal "STEP 2: GX_DB not writable: ${GX_DB}"

  # -------------------- 2A) FCS-adaptor: screen --------------------
  log "[STEP 2A] FCS-adaptor → ${dir2}/adaptor"
  "${ADAPTOR_RUNNER}" \
    --fasta-input "${input_fasta}" \
    --output-dir  "${dir2}/adaptor" \
    --container-engine singularity \
    --image "${sif_fcs_adaptor}" \
    --euk \
    --debug \
    2>&1 | tee -a "${dir2}/logs/fcs_adaptor.run.log"

  # Prefer tool's standard cleaned location; fallback to clean.fasta or original
  after_adaptor="${dir2}/adaptor/cleaned_sequences/$(basename "${input_fasta}")"
  if [[ ! -s "${after_adaptor}" && -s "${dir2}/adaptor/clean.fasta" ]]; then
    after_adaptor="${dir2}/adaptor/clean.fasta"
  fi
  if [[ ! -s "${after_adaptor}" ]]; then
    warn "[STEP 2A] Adaptor-cleaned FASTA not found; falling back to original assembly."
    after_adaptor="${input_fasta}"
  fi

  # -------------------- 2B) Ensure GX DB present (skip if complete) --------------------
  have_gxdb_step2() {
    local d="${GX_DB}"
    # Required files for a complete GX DB
    local -a req=(all.assemblies.tsv all.blast_div.tsv.gz all.gxi all.gxs all.manifest all.meta.jsonl all.README.txt all.seq_info.tsv.gz all.taxa.tsv)
    for f in "${req[@]}"; do [[ -s "${d}/$f" ]] || return 1; done
    # If an interrupted sync marker exists, treat as incomplete
    [[ -d "${d}/gxdb.in_progress" ]] && return 1
    return 0
  }

  mkdir -p "${GX_DB}"
  [[ -w "${GX_DB}" ]] || fatal "STEP 2: GX_DB not writable: ${GX_DB}"

  GX_BIND_MODE="rw"   # default writeable bind
  for gx_lockfile in "${GX_DB}.lockfile" "${GX_DB}/gxdb.lockfile"; do
    if [[ -f "${gx_lockfile}" ]] && ! head -n 1 "${gx_lockfile}" 2>/dev/null | grep -Eq '^[0-9]+,[^,[:space:]]+$'; then
      warn "[STEP 2B] Removing malformed GX DB lockfile: ${gx_lockfile}"
      rm -f "${gx_lockfile}" || fatal "STEP 2: failed to remove malformed lockfile: ${gx_lockfile}"
    fi
  done

  if have_gxdb_step2; then
    log "[STEP 2B] GX DB present → ${GX_DB} (skip download)"
    # Remove stale in-progress marker if any
    [[ -d "${GX_DB}/gxdb.in_progress" ]] && { warn "[STEP 2B] Removing stale gxdb.in_progress"; rm -rf "${GX_DB}/gxdb.in_progress"; }
    GX_BIND_MODE="ro"  # mount read-only so run_gx cannot re-sync
  else
    log "[STEP 2B] Sync GX DB → ${GX_DB}"
    # Clean any partial state before syncing
    rm -rf "${GX_DB}/gxdb.in_progress" 2>/dev/null || true
    singularity exec \
      --bind "${GX_ROOT}:${GX_ROOT}" \
      "${sif_fcs_gx}" \
      python3 /app/bin/sync_files get \
        --mft "${GX_MANIFEST}" \
        --dir "${GX_DB}" \
      2>&1 | tee -a "${dir2}/logs/fcs_gx.sync.log"
    # Re-assess; if still incomplete, bail
    have_gxdb_step2 || fatal "[STEP 2B] GX DB still incomplete after sync -> ${GX_DB}"
    GX_BIND_MODE="ro"
  fi

  # -------------------- 2C) GX screen --------------------
  log "[STEP 2C] GX screen → ${dir2}/gx"
  mkdir -p "${dir2}/gx"

  # Bind DB read-only when complete so run_gx cannot trigger a refresh
  singularity exec \
    --bind "${GX_DB}:/app/db/gxdb:${GX_BIND_MODE}" \
    --bind "$(dirname "${after_adaptor}"):/sample-volume" \
    --bind "${dir2}/gx:/output-volume" \
    "${sif_fcs_gx}" \
    python3 /app/bin/run_gx \
      --fasta "/sample-volume/$(basename "${after_adaptor}")" \
      --gx-db /app/db/gxdb \
      --out-dir /output-volume \
      ${TAXID:+--tax-id "${TAXID}"} \
      --generate-logfile True \
      --debug \
    2>&1 | tee -a "${dir2}/logs/fcs_gx.run.log"

  # -------------------- 2D) GX clean-genome (apply report) --------------------
  # Robustly pick the first GX action report if present (safe under set -e)
  gx_report="$(find "${dir2}/gx" -maxdepth 1 -type f -name '*fcs_gx_report.txt' -print -quit 2>/dev/null || true)"

  if [[ -z "${gx_report}" ]]; then
    log "[STEP 2D] No GX report found — using adaptor-cleaned (or original) as clean.fasta"
    if [[ -s "${after_adaptor}" ]]; then
      ln -sfn "${after_adaptor}" "${dir2}/gx/clean.fasta"
    else
      # Very defensive fallback: if adaptor stage produced nothing, pass through the original
      ln -sfn "${input_fasta}" "${dir2}/gx/clean.fasta"
    fi
    step2_cleaned_fasta="${dir2}/gx/clean.fasta"
  else
    log "[STEP 2D] GX clean-genome → ${dir2}/gx/clean.fasta"
    singularity exec \
      --bind "${dir2}/gx:/report-volume" \
      --bind "${dir2}/gx:/output-volume" \
      --bind "$(dirname "${after_adaptor}"):/input-volume" \
      "${sif_fcs_gx}" \
      /app/bin/gx clean-genome \
        --input "/input-volume/$(basename "${after_adaptor}")" \
        --action-report "/report-volume/$(basename "${gx_report}")" \
        --output "/output-volume/clean.fasta" \
        --contam-fasta-out "/output-volume/contam.fasta" \
        --min-seq-len 200 \
      2>&1 | tee -a "${dir2}/logs/fcs_gx.clean.log"
    step2_cleaned_fasta="${dir2}/gx/clean.fasta"
  fi

  export step2_cleaned_fasta
  log "[STEP 2] Done → ${step2_cleaned_fasta}"
fi

# ============================================================================================================================
# STEP 3: DeepVariant pre-polishing (3 rounds) on the NextPolish2 assembly
# Aligns with centralized configuration and Step 3 output.
# Produces: ${dir3}/${base}.round{1,2,3}.polished.fa  and logs/checkpoints for resuming.
# Final product of this step is typically the round3 polished FASTA.
# ============================================================================================================================
if [[ "${run_deepvariant}" == "TRUE" ]]; then
  echo "[INFO] Step 3: DeepVariant pre-polishing — starting"
  mkdir -p "${dir3}/logs" "${dir3}/checkpoints"

  # -------------------- Inputs --------------------
  initial_assembly="${dir2}/gx/clean.fasta"
  [[ -s "${initial_assembly}" ]] || error "Step 3 input missing: ${initial_assembly}"
  [[ -s "${pacbio_hifi}" ]] || error "PacBio HiFi input not found: ${pacbio_hifi}"

  # -------------------- Helper: ensure dict & fai for a given FASTA --------------------
  ensure_ref_index() {
    local ref="$1"
    local dict="${ref}.dict"
    if [[ ! -s "$dict" ]]; then
      echo "[INFO] Creating sequence dictionary for ${ref}"
      singularity exec -B "${working_dir}:${working_dir}" "${sif_gatk4}" \
        gatk CreateSequenceDictionary -R "$ref" -O "$dict" \
        || error "CreateSequenceDictionary failed for ${ref}"
    fi
    if [[ ! -s "${ref}.fai" ]]; then
      echo "[INFO] Indexing ${ref} with samtools faidx"
      "${samtools_bin}" faidx "$ref" || error "samtools faidx failed for ${ref}"
    fi
  }

  # -------------------- Rounds 1..3 --------------------
  for round in 1 2 3; do
    echo "----------------------------------------------"
    echo "[INFO] Step 3: Round ${round}"
    echo "----------------------------------------------"

    ck_prefix="${dir3}/checkpoints/round${round}"
    asm_in="${initial_assembly}"
    if (( round > 1 )); then
      prev=$(( round - 1 ))
      asm_in="${dir3}/${base}.round${prev}.polished.fa"
    fi
    ensure_ref_index "${asm_in}"

    # Use bam everywhere to avoid BAM QNAME limit (PacBio read IDs can be >255 bytes)
    sorted_bam="${dir3}/${base}.round${round}.sorted.bam"
    markdup_bam="${dir3}/${base}.round${round}.markdup.bam"
    dv_vcf="${dir3}/${base}.round${round}.dv.vcf.gz"
    dv_vcf_hom="${dir3}/${base}.round${round}.dv.hom.vcf.gz"
    asm_out="${dir3}/${base}.round${round}.polished.fa"

    # 3.1 Mapping → sorted bam (avoid BAM)
    if [[ ! -f "${ck_prefix}.mapping.done" ]]; then
      echo "[INFO] Round ${round}: Mapping PacBio HiFi reads to ${asm_in} (bam output)"
      # Write bam directly from sort; pass reference via -T for bam encoding
      "${minimap2_bin}" -H -Y --secondary=no -a -x map-hifi -K 16G -t "$((nthreads * njobs))" -R "@RG\\tID:${base}.HiFi\\tSM:${base}\\tPL:PACBIO\\tLB:SMRTbell" "${asm_in}" "${pacbio_hifi}" \
        | "${samtools_bin}" sort -@ "$((nthreads * njobs))" -m 4G -l 9 --reference "${asm_in}" -T "${dir3}/round${round}.tmp" -O BAM -o "${sorted_bam}" --write-index \
        || error "Round ${round}: Mapping HiFi reads and sorting to bam failed"
      touch "${ck_prefix}.mapping.done"
    else
      echo "[INFO] Round ${round}: using existing sorted bam"
    fi

    # 3.2 MarkDuplicatesSpark (BAM in/out, create CSI index)
    if [[ ! -f "${ck_prefix}.markdup.done" ]]; then
	      echo "[INFO] Marking duplicates with GATK4 (bam)"
	      spark_threads="$((nthreads * njobs))"
		      markdup_shard_tmp="${tmpdir}/gatk_markdup.step3.round${round}.shards"
		      markdup_tmp="${tmpdir}/gatk_markdup.step3.round${round}.tmp"
		      rm -rf "${markdup_shard_tmp}" "${markdup_tmp}"
		      mkdir -p "${markdup_tmp}"
		      rm -f "${markdup_bam}" "${markdup_bam}.bai" "${markdup_bam}.csi"
      gatk_spark \
        --java-options "-Xms128G -Xmx256G -Djava.io.tmpdir=${tmpdir}" \
        MarkDuplicatesSpark \
          --spark-master "local[${spark_threads}]" \
          --conf "spark.driver.host=127.0.0.1" \
          --conf "spark.driver.bindAddress=127.0.0.1" \
          --conf "spark.local.dir=${tmpdir}" \
          --input "${sorted_bam}" \
          --metrics-file "${sorted_bam}.metrics" \
          --output "${markdup_bam}" \
          --reference "${asm_in}" \
	          --remove-all-duplicates true \
	          --create-output-bam-index true \
          --create-output-bam-splitting-index false \
          --read-name-regex null \
          --read-validation-stringency SILENT \
	          --tmp-dir "${markdup_tmp}" \
	          --output-shard-tmp-dir "${markdup_shard_tmp}" \
	          2>&1 | tee -a "${dir3}/logs/round${round}.markduplicates.log" \
	          || error "Round ${round}: MarkDuplicates failed"
	      "${samtools_bin}" index -b -@ "${spark_threads}" "${markdup_bam}" \
	          >>"${dir3}/logs/round${round}.markduplicates.log" 2>&1 \
	          || error "Round ${round}: BAI indexing of MarkDuplicates output failed"
	      rm -rf "${markdup_shard_tmp}" "${markdup_tmp}"
      touch "${ck_prefix}.markdup.done"
    else
      echo "[INFO] Round ${round}: using existing markdup BAM"
    fi

    # 3.3 DeepVariant (reads=bam is OK; bam uses --ref for decoding)
    if [[ ! -f "${ck_prefix}.deepvariant.done" ]]; then
      echo "[INFO] Round ${round}: DeepVariant calling"
      dv_tmp=$(mktemp -d -p "${dir3}" dv.tmp.XXXXXXXX) || error "mktemp failed"
      singularity exec -B "${working_dir}:${working_dir}" "${sif_deepvariant}" /opt/deepvariant/bin/run_deepvariant \
        --model_type=PACBIO \
        --ref="${asm_in}" \
        --reads="${markdup_bam}" \
        --output_vcf="${dv_vcf}" \
        --num_shards="$((nthreads * njobs))" \
        --postprocess_cpus="$((nthreads * njobs))" \
        --intermediate_results_dir="${dv_tmp}" \
        --vcf_stats_report true \
        2>&1 | tee -a "${dir3}/logs/round${round}.deepvariant.log"
      rv=$?
      rm -rf "${dv_tmp}" || true
      [[ $rv -eq 0 ]] || error "Round ${round}: DeepVariant failed"
      touch "${ck_prefix}.deepvariant.done"
    else
      echo "[INFO] Round ${round}: using existing DeepVariant VCF"
    fi

    # 3.4 Filter → consensus
    if [[ ! -f "${ck_prefix}.consensus.done" ]]; then
      echo "[INFO] Round ${round}: filtering VCF and applying consensus"
      "${bcftools_bin}" view -l 0 -f 'PASS' -i 'GT="1/1"' "${dv_vcf}" \
        | "${bcftools_bin}" view --threads "$((nthreads * njobs))" -Oz9 -o "${dv_vcf_hom}" --write-index \
        || error "Round ${round}: bcftools filtering failed"

      "${bcftools_bin}" consensus -f "${asm_in}" -o "${asm_out}" "${dv_vcf_hom}" \
        || error "Round ${round}: consensus failed"
      touch "${ck_prefix}.consensus.done"
    else
      echo "[INFO] Round ${round}: consensus already complete"
    fi

    # 3.5 Cleanup
    if [[ -f "${ck_prefix}.consensus.done" ]]; then
      rm -f "${sorted_bam}"* "${markdup_bam}"* "${dv_vcf}"* "${dv_vcf_hom}"* || true
      rm -f "${asm_in}.fai" "${asm_in}.dict" || true
    fi

    echo "[INFO] Round ${round} complete: ${asm_out}"
  done
fi

# ============================================================================================================================
# STEP 4: Inspector — reference-free long-read assembly evaluation & correction (HiFi)
# Aligns with centralized configuration and Step 4 output.
# Consumes: ${dir3}/${base}.round3.polished.fa  (from DeepVariant pre-polish)
# Produces: ${dir4}/inspector_out/contig_corrected.fa
# ============================================================================================================================
if [[ "${run_inspector}" == "TRUE" ]]; then
  echo "[INFO] Step 4: Inspector — starting"
  mkdir -p "${dir4}/logs" "${dir4}/inspector_out" "${dir4}/checkpoints"

  # -------------------- Inputs --------------------
  asm_in="${dir3}/${base}.round3.polished.fa"
  [[ -s "${asm_in}" ]] || error "Step 4 input missing: ${asm_in}"
  [[ -s "${pacbio_hifi}" ]] || error "HiFi reads missing: ${pacbio_hifi}"

  # -------------------- Activate env --------------------
  activate_env "${ENV_INSPECTOR}"

  # -------------------- Paths & checkpoints --------------------
  outdir_ins="${dir4}/inspector_out"
  ck_eval="${dir4}/checkpoints/inspector.eval.done"
  ck_corr="${dir4}/checkpoints/inspector.correct.done"

  # -------------------- 5a) Inspector evaluation --------------------
  if [[ -s "${ck_eval}" ]]; then
    echo "[INFO] Inspector evaluation already done"
  else
    echo "[INFO] Running inspector.py on ${asm_in}"
    inspector.py \
      --contig "${asm_in}" \
      --read   "${pacbio_hifi}" \
      --outpath "${outdir_ins}" \
      --datatype hifi \
      --thread "$((nthreads * njobs))" \
      --min_contig_length 200 \
      --min_contig_length_assemblyerror 200 \
      >"${dir4}/logs/inspector.log" 2>"${dir4}/logs/inspector.err" \
      || error "inspector.py failed"
    touch "${ck_eval}"
  fi

  # -------------------- 5b) Inspector-based correction --------------------
  corrected_fa="${outdir_ins}/contig_corrected.fa"

  if [[ -s "${ck_corr}" && -s "${corrected_fa}" ]]; then
    echo "[INFO] Inspector correction already done: ${corrected_fa}"
  else
    echo "[INFO] Running inspector-correct.py"
    inspector-correct.py \
      --inspector "${outdir_ins}" \
      --datatype pacbio-hifi \
      --outpath  "${outdir_ins}" \
      --thread   "$((nthreads * njobs))" \
      >"${dir4}/logs/inspector-correct.log" 2>"${dir4}/logs/inspector-correct.err" \
      || error "inspector-correct.py failed"
    [[ -s "${corrected_fa}" ]] || error "Inspector-correct did not produce ${corrected_fa}"
    touch "${ck_corr}"
  fi

  deactivate_env
fi

# ============================================================================================================================
# STEP 5: Scaffolding — LongStitch grid, GenomeScope2, QUAST, optional yak, Merqury, robust selection
# Consumes: ${dir4}/inspector_out/contig_corrected.fa plus high-accuracy reads when available
# Produces: ${dir5}/${base}.scaffolds.fa (selected scaffold set) + scored selection artifacts
# ============================================================================================================================
if is_true "${run_scaffolding}"; then
  echo "[INFO] Step 5: Scaffolding — starting"
  mkdir -p \
    "${dir5}/logs" \
    "${dir5}/checkpoints" \
    "${dir5}/kmers" \
    "${dir5}/genomescope2" \
    "${dir5}/quast" \
    "${dir5}/yak" \
    "${dir5}/merqury" \
    "${dir5}/selection"

  scaff_ck_dir="${dir5}/checkpoints"
  scaff_logs_dir="${dir5}/logs"
  scaff_kmers_dir="${dir5}/kmers"
  scaff_gs_dir="${dir5}/genomescope2"
  scaff_quast_dir="${dir5}/quast"
  scaff_yak_dir="${dir5}/yak"
  scaff_merqury_dir="${dir5}/merqury"
  scaff_selection_dir="${dir5}/selection"

  scaffolding_checkpoint_file() {
    printf '%s/%s.ok' "$scaff_ck_dir" "$1"
  }

  scaffolding_merqury_metrics_ready() {
    local asm asm_id metrics_file found=0
    [[ -s "${scaff_selection_dir}/top_candidates.txt" ]] || return 1
    while IFS= read -r asm || [[ -n "$asm" ]]; do
      [[ -n "$asm" ]] || continue
      found=1
      asm_id="$(basename "${asm}")"
      asm_id="${asm_id%.fa}"; asm_id="${asm_id%.fasta}"; asm_id="${asm_id%.fna}"
      metrics_file="${scaff_merqury_dir}/${asm_id}/merqury_metrics.tsv"
      [[ -s "$metrics_file" ]] || return 1
      awk -F'\t' 'NR > 1 && NF >= 4 { ok = 1 } END { exit ok ? 0 : 1 }' "$metrics_file" || return 1
    done < "${scaff_selection_dir}/top_candidates.txt"
    (( found == 1 ))
  }

  scaffolding_genomescope2_ready() {
    [[ -s "${gs_summary:-}" && -s "${gs_model:-}" && -s "${gs_progress:-}" && -s "${gs_metrics:-}" ]] || return 1
    grep -Eq 'Genome|Heterozygosity|Model' "${gs_summary}" || return 1
  }

  scaffolding_longstitch_grid_ready() {
    local f found=0
    shopt -s nullglob
    for f in "${dir5}"/kw_k*_w*/*.scaffolds.fa; do
      [[ -s "$f" ]] || continue
      found=1
      break
    done
    shopt -u nullglob
    (( found == 1 ))
  }

  scaffolding_output_ready() {
    local step="${1:-}"
    case "$step" in
      05A.inputs) [[ -s "${input_manifest:-}" ]] ;;
      05B.kmers) [[ -s "${gs_hist:-}" && -e "${read_meryl_db:-}" ]] ;;
      05C.genomescope2) scaffolding_genomescope2_ready ;;
      05D.longstitch_grid) scaffolding_longstitch_grid_ready ;;
      05E.quast) [[ -s "${scaff_quast_dir}/report.tsv" ]] ;;
      05F.yak) [[ -s "${yak_metrics:-}" ]] ;;
      05G.preselect) [[ -s "${scaff_selection_dir}/top_candidates.txt" && -s "${scaff_selection_dir}/preselection_scores.tsv" ]] ;;
      05H.merqury) scaffolding_merqury_metrics_ready ;;
      05I.selection) [[ -s "${scaff_selection_dir}/best_assembly.path" && -s "${scaff_selection_dir}/longstitch_scores.tsv" ]] ;;
      05J.final_scaffolds) [[ -s "${dir5}/${base}.scaffolds.fa" && -s "${dir5}/best_assembly.path" ]] ;;
      *) return 1 ;;
    esac
  }

  scaffolding_checkpoint_ready() {
    local step="$1"
    local ck
    ck="$(scaffolding_checkpoint_file "$step")"
    if scaffolding_output_ready "$step"; then
      [[ -f "$ck" ]] || touch "$ck"
      return 0
    fi
    if [[ -f "$ck" ]]; then
      warn "Removing stale scaffolding checkpoint without valid outputs: $ck"
      rm -f "$ck"
    fi
    return 1
  }

  scaffolding_mark_checkpoint() {
    local step="$1"
    scaffolding_output_ready "$step" || error "Refusing to mark incomplete scaffolding checkpoint: $step"
    touch "$(scaffolding_checkpoint_file "$step")"
  }

  # -------------------- STEP 5A: Input resolution and preflight --------------------
  resolve_scaffolding_draft() {
    local cand
    for cand in \
      "${dir4}/inspector_out/contig_corrected.fa" \
      "${dir3}/${base}.round3.polished.fa" \
      "${dir2}/cleaned.fasta" \
      "${dir1}/${base}.bp.p_ctg.fa" \
      "${dir1}/${base}.asm.bp.p_ctg.fa"; do
      [[ -s "${cand}" ]] && { printf '%s\n' "${cand}"; return 0; }
    done
    return 1
  }

  resolve_longstitch_reads() {
    if [[ -n "${pacbio_clr:-}" && -s "${pacbio_clr}" ]]; then
      printf '%s\n' "${pacbio_clr}"
    elif [[ -n "${pacbio_hifi_fasta:-}" && -s "${pacbio_hifi_fasta}" ]]; then
      printf '%s\n' "${pacbio_hifi_fasta}"
    elif [[ -n "${pacbio_hifi:-}" && -s "${pacbio_hifi}" ]]; then
      printf '%s\n' "${pacbio_hifi}"
    else
      return 1
    fi
  }

  prepare_longstitch_reads() {
    local in="$1" out="$2"
    local tmp="${out}.tmp.$$"

    [[ -s "${in}" ]] || error "LongStitch read input missing: ${in}"
    mkdir -p "$(dirname "${out}")"

    if [[ -s "${out}" && "${out}" -nt "${in}" ]]; then
      printf '%s\n' "${out}"
      return 0
    fi

    rm -f "${tmp}"
    if [[ "${in}" == *.gz ]]; then
      require gzip
      gzip -dc "${in}" | awk '
        BEGIN { mode = ""; fqline = 0 }
        mode == "" {
          if ($0 ~ /^[[:space:]]*$/) next
          first = substr($0, 1, 1)
          if (first == ">") mode = "fasta"
          else if (first == "@") mode = "fastq"
          else {
            print "[ERROR] LongStitch reads are not FASTA/FASTQ" > "/dev/stderr"
            exit 2
          }
        }
        mode == "fasta" { print; next }
        mode == "fastq" {
          fqline++
          if (fqline % 4 == 1) { sub(/^@/, ">"); print; next }
          if (fqline % 4 == 2) { print; next }
        }
        END { if (mode == "") exit 2 }
      ' > "${tmp}" || { rm -f "${tmp}"; error "Failed to convert LongStitch reads to FASTA: ${in}"; }
    else
      awk '
        BEGIN { mode = ""; fqline = 0 }
        mode == "" {
          if ($0 ~ /^[[:space:]]*$/) next
          first = substr($0, 1, 1)
          if (first == ">") mode = "fasta"
          else if (first == "@") mode = "fastq"
          else {
            print "[ERROR] LongStitch reads are not FASTA/FASTQ" > "/dev/stderr"
            exit 2
          }
        }
        mode == "fasta" { print; next }
        mode == "fastq" {
          fqline++
          if (fqline % 4 == 1) { sub(/^@/, ">"); print; next }
          if (fqline % 4 == 2) { print; next }
        }
        END { if (mode == "") exit 2 }
      ' "${in}" > "${tmp}" || { rm -f "${tmp}"; error "Failed to convert LongStitch reads to FASTA: ${in}"; }
    fi

    [[ -s "${tmp}" ]] || { rm -f "${tmp}"; error "LongStitch FASTA conversion produced an empty file: ${out}"; }
    mv -f "${tmp}" "${out}"
    printf '%s\n' "${out}"
  }

  draft="$(resolve_scaffolding_draft)" || error "Step 5 input draft assembly not found"
  longstitch_raw_reads="$(resolve_longstitch_reads)" || error "No valid LongStitch reads found"
  longstitch_longmap="hifi"
  if [[ -n "${pacbio_clr:-}" && "${longstitch_raw_reads}" == "${pacbio_clr}" ]]; then
    longstitch_longmap="pacbio"
  fi
  longstitch_reads="${longstitch_raw_reads}"
  case "${longstitch_raw_reads}" in
    *.gz|*.fq|*.fastq)
      longstitch_reads="$(prepare_longstitch_reads "${longstitch_raw_reads}" "${dir5}/inputs/${base}.longstitch_reads.fa")"
      ;;
  esac

  read_inputs=()
  [[ -n "${illumina_r1:-}" && -s "${illumina_r1:-}" ]] && read_inputs+=("${illumina_r1}")
  [[ -n "${illumina_r2:-}" && -s "${illumina_r2:-}" ]] && read_inputs+=("${illumina_r2}")
  if (( ${#read_inputs[@]} == 0 )); then
    [[ -n "${pacbio_hifi:-}" && -s "${pacbio_hifi}" ]] && read_inputs+=("${pacbio_hifi}")
  fi
  (( ${#read_inputs[@]} > 0 )) || read_inputs+=("${longstitch_reads}")

  input_manifest="${scaff_selection_dir}/step5_input_manifest.tsv"
  {
    printf "key\tvalue\n"
    printf "draft\t%s\n" "${draft}"
    printf "longstitch_raw_reads\t%s\n" "${longstitch_raw_reads}"
    printf "longstitch_reads\t%s\n" "${longstitch_reads}"
    printf "longstitch_longmap\t%s\n" "${longstitch_longmap}"
    printf "longstitch_gap_fill\t%s\n" "${LONGSTITCH_GAP_FILL}"
    printf "kmer_reads\t%s\n" "${read_inputs[*]}"
    printf "longstitch_k_values\t%s\n" "${LONGSTITCH_K_VALUES}"
    printf "longstitch_w_values\t%s\n" "${LONGSTITCH_W_VALUES}"
  } > "${input_manifest}"
  scaffolding_mark_checkpoint 05A.inputs

  activate_env "${ENV_LONGSTITCH}"
  require longstitch
  require parallel
  deactivate_env

  activate_env "${ENV_GENOMESCOPE2}"
  if [[ -n "${GENOMESCOPE2_BIN}" ]]; then
    require "${GENOMESCOPE2_BIN}"
  else
    require genomescope2
  fi
  deactivate_env

  activate_env "${ENV_MERQURY}"
  require meryl
  if [[ -n "${MERQURY}" ]]; then
    require "${MERQURY}/merqury.sh"
  else
    require merqury.sh
  fi
  deactivate_env

  activate_env "${ENV_R}"
  require Rscript
  deactivate_env
  env_run "${ENV_R}" Rscript --vanilla -e 'if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Renv is missing required R package: ggplot2"); cat("Renv ggplot2 OK\n")' \
    > "${scaff_logs_dir}/renv.preflight.log" 2>&1 \
    || error "Renv preflight failed; install ggplot2 in ${ENV_R}"

  activate_env "${ENV_HIFIASM}"
  require yak
  deactivate_env

  # -------------------- STEP 5B: K-mer database preparation --------------------
  gs_meryl_db="${scaff_kmers_dir}/${base}.reads.k${GS_K}.meryl"
  gs_hist="${scaff_kmers_dir}/${base}.reads.k${GS_K}.hist"
  read_meryl_db="${scaff_kmers_dir}/${base}.reads.k${MERQURY_K}.meryl"

  if ! scaffolding_checkpoint_ready 05B.kmers; then
    activate_env "${ENV_MERQURY}"
    export OMP_NUM_THREADS="${nthreads}"
    if [[ ! -d "${gs_meryl_db}" ]]; then
      meryl k="${GS_K}" count output "${gs_meryl_db}" "${read_inputs[@]}" \
        > "${scaff_logs_dir}/meryl.genomescope2.k${GS_K}.log" 2>&1 \
        || error "meryl count failed for GenomeScope2 k=${GS_K}"
    fi
    if [[ ! -s "${gs_hist}" ]]; then
      meryl histogram "${gs_meryl_db}" > "${gs_hist}" \
        2> "${scaff_logs_dir}/meryl.genomescope2.hist.k${GS_K}.log" \
        || error "meryl histogram failed for GenomeScope2 k=${GS_K}"
    fi
    if [[ "${MERQURY_K}" == "${GS_K}" && -d "${gs_meryl_db}" && ! -e "${read_meryl_db}" && ! -L "${read_meryl_db}" ]]; then
      ln -s "${gs_meryl_db}" "${read_meryl_db}"
    fi
    if [[ ! -d "${read_meryl_db}" ]]; then
      meryl k="${MERQURY_K}" count output "${read_meryl_db}" "${read_inputs[@]}" \
        > "${scaff_logs_dir}/meryl.merqury.k${MERQURY_K}.log" 2>&1 \
        || error "meryl count failed for Merqury k=${MERQURY_K}"
    fi
    deactivate_env
    scaffolding_mark_checkpoint 05B.kmers
  fi

  # -------------------- STEP 5C: GenomeScope2 genome profiling --------------------
  gs_name_prefix="${GS_NAME_PREFIX:-${base}}"
  gs_output_prefix="${gs_name_prefix:+${gs_name_prefix}_}"
  gs_summary="${scaff_gs_dir}/${gs_output_prefix}summary.txt"
  gs_model="${scaff_gs_dir}/${gs_output_prefix}model.txt"
  gs_progress="${scaff_gs_dir}/${gs_output_prefix}progress.txt"
  gs_json="${scaff_gs_dir}/${gs_output_prefix}summary.json"
  gs_testing="${scaff_gs_dir}/${gs_output_prefix}testing.tsv"
  gs_metrics="${scaff_gs_dir}/genomescope2.metrics.tsv"
  if ! scaffolding_checkpoint_ready 05C.genomescope2; then
    activate_env "${ENV_GENOMESCOPE2}"
    genomescope_cmd="${GENOMESCOPE2_BIN:-genomescope2}"
    genomescope_args=(
      -i "${gs_hist}"
      -o "${scaff_gs_dir}"
      -p "${GS_PLOIDY}"
      -k "${GS_K}"
      -n "${gs_name_prefix}"
      --transform_exp "${GS_TRANSFORM_EXP}"
    )
    is_true "${GS_VERBOSE}" && genomescope_args+=(--verbose)
    is_true "${GS_NO_UNIQUE_SEQUENCE}" && genomescope_args+=(--no_unique_sequence)
    is_true "${GS_JSON_REPORT}" && genomescope_args+=(--json_report)
    is_true "${GS_FITTED_HIST}" && genomescope_args+=(--fitted_hist)
    [[ -n "${GS_MAX_KMERCOV}" ]] && genomescope_args+=(-m "${GS_MAX_KMERCOV}")
    [[ -n "${GS_LAMBDA}" ]] && genomescope_args+=(-l "${GS_LAMBDA}")
    [[ -n "${GS_TOPOLOGY}" ]] && genomescope_args+=(-t "${GS_TOPOLOGY}")
    [[ -n "${GS_INITIAL_REPETITIVENESS}" ]] && genomescope_args+=(--initial_repetitiveness "${GS_INITIAL_REPETITIVENESS}")
    [[ -n "${GS_INITIAL_HETEROZYGOSITIES}" ]] && genomescope_args+=(--initial_heterozygosities "${GS_INITIAL_HETEROZYGOSITIES}")
    is_true "${GS_TESTING}" && genomescope_args+=(--testing)
    [[ -n "${GS_TRUE_PARAMS}" ]] && genomescope_args+=(--true_params "${GS_TRUE_PARAMS}")
    is_true "${GS_TRACE_FLAG}" && genomescope_args+=(--trace_flag)
    [[ -n "${GS_NUM_ROUNDS}" ]] && genomescope_args+=(--num_rounds "${GS_NUM_ROUNDS}")
    [[ -n "${GS_START_SHIFT}" ]] && genomescope_args+=(--start_shift "${GS_START_SHIFT}")
    [[ -n "${GS_TYPICAL_ERROR}" ]] && genomescope_args+=(--typical_error "${GS_TYPICAL_ERROR}")
    "${genomescope_cmd}" "${genomescope_args[@]}" \
      > "${scaff_logs_dir}/genomescope2.log" 2>&1 \
      || error "GenomeScope2 failed"
    deactivate_env
    [[ -s "${gs_summary}" ]] || error "GenomeScope2 summary missing: ${gs_summary}"
    [[ -s "${gs_model}" ]] || error "GenomeScope2 model missing: ${gs_model}"
    [[ -s "${gs_progress}" ]] || error "GenomeScope2 progress missing: ${gs_progress}"
    {
      printf "metric\tvalue\n"
      printf "summary_file\t%s\n" "${gs_summary}"
      printf "model_file\t%s\n" "${gs_model}"
      printf "progress_file\t%s\n" "${gs_progress}"
      [[ -s "${gs_json}" ]] && printf "json_report\t%s\n" "${gs_json}"
      [[ -s "${gs_testing}" ]] && printf "testing_file\t%s\n" "${gs_testing}"
      awk '
        /^[[:space:]]*Genome Haploid Length/ {
          line = $0
          sub(/^[[:space:]]*Genome Haploid Length[[:space:]]+/, "", line)
          print "genome_haploid_length\t" line
        }
        /^[[:space:]]*Genome Repeat Length/ {
          line = $0
          sub(/^[[:space:]]*Genome Repeat Length[[:space:]]+/, "", line)
          print "genome_repeat_length\t" line
        }
        /^[[:space:]]*Genome Unique Length/ {
          line = $0
          sub(/^[[:space:]]*Genome Unique Length[[:space:]]+/, "", line)
          print "genome_unique_length\t" line
        }
        /^[[:space:]]*Heterozygosity/ {
          line = $0
          sub(/^[[:space:]]*Heterozygosity[[:space:]]+/, "", line)
          print "heterozygosity\t" line
        }
        /^[[:space:]]*Model Fit/ {
          line = $0
          sub(/^[[:space:]]*Model Fit[[:space:]]+/, "", line)
          print "model_fit\t" line
        }
        /^[[:space:]]*Read Error Rate/ {
          line = $0
          sub(/^[[:space:]]*Read Error Rate[[:space:]]+/, "", line)
          print "read_error_rate\t" line
        }
      ' "${gs_summary}"
    } > "${gs_metrics}"
    scaffolding_genomescope2_ready || error "GenomeScope2 outputs did not pass validation"
    scaffolding_mark_checkpoint 05C.genomescope2
  fi
  if [[ -z "${GENOME_SIZE_EST}" ]]; then
    gs_est="$(awk '
      /^[[:space:]]*Genome Haploid Length/ {
        line = $0
        gsub(/,/, "", line)
        if (match(line, /[0-9]+/)) {
          print substr(line, RSTART, RLENGTH)
          exit
        }
      }
    ' "${gs_summary}" 2>/dev/null || true)"
    [[ -n "${gs_est}" ]] && GENOME_SIZE_EST="${gs_est}"
  fi

  longstitch_genome_size="${GENOME_SIZE:-${GENOME_SIZE_EST:-}}"

  # -------------------- STEP 5D: LongStitch broad grid --------------------
  run_scaffolding_kw() {
    local k="${1:-}" w="${2:-}"
    local run_dir out_scaff reads_base draft_base ls_status
    local -a ls_args=()

    [[ -n "${k}" && -n "${w}" ]] || { echo "[ERROR] run_scaffolding_kw: missing k or w" >&2; return 2; }

    run_dir="${dir5}/kw_k${k}_w${w}"
    mkdir -p "${run_dir}"
    out_scaff="${run_dir}/${base}.k${k}.w${w}.scaffolds.fa"

    if [[ -s "${out_scaff}" && -s "${run_dir}/.longstitch.ok" ]]; then
      echo "[INFO] [k=${k}, w=${w}] LongStitch output exists; skipping"
      return 0
    fi

    find "${run_dir}" -maxdepth 1 \( -type f -o -type l \) -name "contig_corrected.fa.k${k}.w${w}*" -delete
    rm -f "${out_scaff}" "${run_dir}/.longstitch.ok"

    echo "[INFO] [k=${k}, w=${w}] Running LongStitch in ${run_dir}"
    pushd "${run_dir}" >/dev/null || { echo "[ERROR] cd failed for ${run_dir}" >&2; return 1; }

    longstitch make_links -C "${run_dir}" draft_path="${draft}" reads_path="${longstitch_reads}" \
      || { echo "[ERROR] make_links failed k=${k} w=${w}" >&2; popd >/dev/null; return 1; }

    draft_base="$(basename "${draft}")"
    draft_base="${draft_base%.gz}"; draft_base="${draft_base%.fasta}"; draft_base="${draft_base%.fa}"; draft_base="${draft_base%.fna}"
    reads_base="$(basename "${longstitch_reads}")"
    reads_base="${reads_base%.gz}"; reads_base="${reads_base%.fastq}"; reads_base="${reads_base%.fq}"
    reads_base="${reads_base%.fasta}"; reads_base="${reads_base%.fa}"; reads_base="${reads_base%.fna}"

    ls_args=(
      "${LONGSTITCH_MODE}"
      draft="${draft_base}"
      reads="${reads_base}"
      out_prefix="${base}.k${k}.w${w}"
      w="${w}"
      k_ntLink="${k}"
      longmap="${longstitch_longmap}"
      gap_fill=True
      rounds="${LONGSTITCH_ROUNDS}"
      t="${nthreads}"
    )
    [[ -n "${longstitch_genome_size}" ]] && ls_args+=(G="${longstitch_genome_size}" genome_size="${longstitch_genome_size}")

    longstitch "${ls_args[@]}" 2>&1 | tee -a "${run_dir}/longstitch.log"
    ls_status=${PIPESTATUS[0]}
    if (( ls_status != 0 )); then
      echo "[ERROR] LongStitch failed k=${k} w=${w} (rc=${ls_status})" >&2
      popd >/dev/null
      return "${ls_status}"
    fi

    [[ -s "${out_scaff}" ]] || { echo "[ERROR] LongStitch did not produce ${out_scaff}" >&2; popd >/dev/null; return 1; }
    touch "${run_dir}/.longstitch.ok"
    popd >/dev/null
  }
  export -f run_scaffolding_kw
  export dir5 base draft longstitch_reads longstitch_longmap longstitch_genome_size nthreads LONGSTITCH_MODE LONGSTITCH_ROUNDS

  read -r -a k_values <<< "${LONGSTITCH_K_VALUES}"
  read -r -a w_values <<< "${LONGSTITCH_W_VALUES}"
  (( ${#k_values[@]} > 0 && ${#w_values[@]} > 0 )) || error "LONGSTITCH_K_VALUES or LONGSTITCH_W_VALUES is empty"

  if ! scaffolding_checkpoint_ready 05D.longstitch_grid; then
    activate_env "${ENV_LONGSTITCH}"
    grid_log="${scaff_logs_dir}/longstitch.grid.log"
    parallel --jobs "${njobs}" --load 80% --memfree 128G --tmpdir "${tmpdir}" \
      --env run_scaffolding_kw --env dir5 --env base --env draft --env longstitch_reads --env longstitch_genome_size \
      --env nthreads --env longstitch_longmap --env LONGSTITCH_MODE --env LONGSTITCH_ROUNDS \
      run_scaffolding_kw "{1}" "{2}" ::: "${k_values[@]}" ::: "${w_values[@]}" \
      2>&1 | tee -a "${grid_log}"
    grid_rc=${PIPESTATUS[0]:-0}
    (( grid_rc == 0 )) || error "One or more LongStitch grid jobs failed (see ${grid_log})"
    deactivate_env
    scaffolding_mark_checkpoint 05D.longstitch_grid
  fi

  candidate_manifest="${scaff_selection_dir}/candidate_manifest.tsv"
  printf "assembly_id\tk\tw\tpath\n" > "${candidate_manifest}"
  shopt -s nullglob
  for f in "${dir5}"/kw_k*_w*/"${base}".k*.w*.scaffolds.fa; do
    [[ -s "${f}" ]] || continue
    run_leaf="$(basename "$(dirname "${f}")")"
    k_val="$(sed -E 's/.*kw_k([0-9]+)_w([0-9]+).*/\1/' <<< "${run_leaf}")"
    w_val="$(sed -E 's/.*kw_k([0-9]+)_w([0-9]+).*/\2/' <<< "${run_leaf}")"
    asm_id="$(basename "${f}")"
    asm_id="${asm_id%.fa}"; asm_id="${asm_id%.fasta}"; asm_id="${asm_id%.fna}"
    printf "%s\t%s\t%s\t%s\n" "${asm_id}" "${k_val}" "${w_val}" "${f}" >> "${candidate_manifest}"
  done
  shopt -u nullglob
  (( "$(wc -l < "${candidate_manifest}")" > 1 )) || error "No scaffold assemblies found after LongStitch grid"

  # -------------------- STEP 5E: QUAST structural evaluation --------------------
  if ! scaffolding_checkpoint_ready 05E.quast; then
    activate_env "${ENV_QUAST}"
    mapfile -t quast_inputs < <(awk -F'\t' 'NR > 1 {print $4}' "${candidate_manifest}")
    (( ${#quast_inputs[@]} > 0 )) || error "No scaffold assemblies found for QUAST evaluation"
    quast_args=(
      -o "${scaff_quast_dir}"
      --threads "$(( nthreads * njobs ))"
      --large
      --min-contig "${QUAST_MIN_CONTIG}"
      --report-all-metrics
      --k-mer-stats
      --k-mer-size "${QUAST_KMER_SIZE}"
      --contig-thresholds "1000,5000,10000,100000,500000,1000000,5000000,10000000"
      --x-for-Nx 90
      --no-snps
      --no-icarus
    )
    is_true "${QUAST_SPLIT_SCAFFOLDS}" && quast_args+=(--split-scaffolds)
    [[ -n "${GENOME_SIZE_EST}" ]] && quast_args+=(--est-ref-size "${GENOME_SIZE_EST}")
    if is_true "${QUAST_USE_READS}"; then
      quast_args+=(--pacbio "${longstitch_reads}")
    fi
    quast "${quast_args[@]}" "${quast_inputs[@]}" \
      2>&1 | tee -a "${scaff_logs_dir}/quast.longstitch.log"
    quast_rc=${PIPESTATUS[0]:-0}
    (( quast_rc == 0 )) || error "QUAST failed on Step 5 scaffold candidates"
    deactivate_env
    [[ -s "${scaff_quast_dir}/report.tsv" ]] || error "QUAST report.tsv not found in ${scaff_quast_dir}"
    scaffolding_mark_checkpoint 05E.quast
  fi

  # -------------------- STEP 5F: Optional yak fast screening --------------------
  yak_metrics="${scaff_yak_dir}/yak_metrics.tsv"
  if ! scaffolding_checkpoint_ready 05F.yak; then
    activate_env "${ENV_HIFIASM}"
    yak_db="${scaff_yak_dir}/${base}.reads.k${YAK_K}.yak"
    [[ -s "${yak_db}" ]] || yak count -k "${YAK_K}" -b "${YAK_BLOOM_BITS}" -t "$(( nthreads * njobs ))" \
      -o "${yak_db}" "${read_inputs[@]}" \
      > "${scaff_logs_dir}/yak.count.log" 2>&1 \
      || error "yak count failed"
    printf "assembly_id\tyak_qv\tyak_error_rate\tpath\n" > "${yak_metrics}"
    while IFS=$'\t' read -r asm_id _ _ asm_path; do
      [[ "${asm_id}" == "assembly_id" ]] && continue
      yak qv -t "${nthreads}" "${yak_db}" "${asm_path}" \
        > "${scaff_yak_dir}/${asm_id}.yak_qv.txt" 2>> "${scaff_logs_dir}/yak.qv.log" || true
      yak_line="$(tail -n1 "${scaff_yak_dir}/${asm_id}.yak_qv.txt" 2>/dev/null || true)"
      yak_qv="$(awk '{for (i=1;i<=NF;i++) if ($i ~ /^[0-9]+([.][0-9]+)?$/) {v=$i}} END{print v}' <<< "${yak_line}")"
      printf "%s\t%s\t\t%s\n" "${asm_id}" "${yak_qv:-NA}" "${asm_path}" >> "${yak_metrics}"
    done < "${candidate_manifest}"
    deactivate_env
    scaffolding_mark_checkpoint 05F.yak
  fi

  # -------------------- STEP 5G: Select top candidates before Merqury --------------------
  selector_r="${scaff_selection_dir}/select_step5_candidates.R"
  cat > "${selector_r}" <<'RSCRIPT'
#!/usr/bin/env Rscript

# ==============================================================================
# LongStitch candidate selector
# - Preselection mode: QUAST + optional yak
# - Final mode: QUAST + optional yak + Merqury
# - Outputs ranking tables, top candidates, best assembly path, and ggplot2 plots
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Missing required R package: ggplot2")
}

arg <- function(flag, default = NA_character_) {
  hit <- match(flag, args)
  if (is.na(hit) || hit == length(args)) return(default)
  args[[hit + 1]]
}

read_tsv <- function(path) {
  if (is.na(path) || !nzchar(path)) return(data.frame())
  if (!file.exists(path) || file.info(path)$size == 0) return(data.frame())
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

to_num <- function(x) {
  suppressWarnings(as.numeric(gsub("[,% ]", "", as.character(x))))
}

metric_row <- function(tab, labels) {
  if (nrow(tab) == 0) return(NULL)
  idx <- which(tolower(tab[[1]]) %in% tolower(labels))
  if (!length(idx)) return(NULL)
  tab[idx[1], -1, drop = FALSE]
}

score_col <- function(x, direction = 1) {
  x <- to_num(x)

  if (all(!is.finite(x))) {
    return(rep(0, length(x)))
  }

  med <- stats::median(x, na.rm = TRUE)
  madv <- stats::mad(x, constant = 1.4826, na.rm = TRUE)

  if (!is.finite(madv) || madv == 0) {
    madv <- stats::sd(x, na.rm = TRUE)
  }

  if (!is.finite(madv) || madv == 0) {
    madv <- 1
  }

  z <- direction * ((x - med) / madv)
  z[!is.finite(z)] <- 0
  z
}

normalize_id <- function(x) {
  x <- basename(as.character(x))
  x <- sub("\\.fasta$", "", x)
  x <- sub("\\.fa$", "", x)
  x <- sub("\\.fna$", "", x)
  x
}

long_metrics <- function(quast_path, manifest) {
  q <- read_tsv(quast_path)

  if (nrow(q) == 0 || ncol(q) < 2) {
    stop("QUAST report is empty or malformed: ", quast_path)
  }

  asm <- names(q)[-1]

  out <- data.frame(
    assembly_id = asm,
    assembly_id_norm = normalize_id(asm),
    stringsAsFactors = FALSE
  )

  get <- function(labels) {
    row <- metric_row(q, labels)
    if (is.null(row)) {
      rep(NA_real_, length(asm))
    } else {
      to_num(unlist(row, use.names = FALSE))
    }
  }

  out$n50 <- get(c("N50"))
  out$l50 <- get(c("L50"))
  out$n90 <- get(c("N90"))
  out$l90 <- get(c("L90"))
  out$ng50 <- get(c("NG50"))
  out$lg50 <- get(c("LG50"))
  out$largest_scaffold <- get(c("Largest contig", "Largest scaffold"))
  out$scaffold_count <- get(c("# contigs", "# scaffolds"))
  out$total_length <- get(c("Total length"))
  out$ns_per_100kbp <- get(c("# N's per 100 kbp", "# N's per 100kbp"))
  out$gc_percent <- get(c("GC (%)"))
  out$auN <- get(c("auN"))
  out$auNG <- get(c("auNG"))

  if (!"assembly_id_norm" %in% names(manifest)) {
    manifest$assembly_id_norm <- normalize_id(manifest$assembly_id)
  }

  merged <- merge(
    manifest,
    out,
    by = "assembly_id_norm",
    all.x = TRUE,
    suffixes = c("", ".quast")
  )

  if ("assembly_id.quast" %in% names(merged)) {
    merged$quast_assembly_id <- merged$assembly_id.quast
    merged$assembly_id.quast <- NULL
  }

  merged
}

safe_pdf <- function(outdir, filename, plot, width = 9, height = 6) {
  ggplot2::ggsave(
    filename = file.path(outdir, filename),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf
  )
}

safe_png <- function(outdir, filename, plot, width = 9, height = 6, dpi = 600) {
  ggplot2::ggsave(
    filename = file.path(outdir, filename),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi
  )
}

safe_svg <- function(outdir, filename, plot, width = 9, height = 6, dpi = 600) {
  ggplot2::ggsave(
    filename = file.path(outdir, filename),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = "svg",
    dpi = dpi
  )
}

save_selection_plot <- function(outdir, stem, plot, width = 9, height = 6) {
  safe_pdf(outdir, paste0(stem, ".pdf"), plot, width, height)
  safe_png(outdir, paste0(stem, ".png"), plot, width, height, plot_dpi)
  if (isTRUE(save_svg)) {
    safe_svg(outdir, paste0(stem, ".svg"), plot, width, height, plot_dpi)
  }
}

theme_selection <- function(base_size = 12) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
}

plot_selection_results <- function(metrics, top, outdir, final = FALSE, genome_size_est = NA_real_) {
  plot_df <- metrics[is.finite(metrics$score), , drop = FALSE]

  # 1. Ranking plot ------------------------------------------------------------
  if (nrow(top) > 0) {
    top$assembly_id_ranked <- stats::reorder(top$assembly_id, top$score)

    p_rank <- ggplot2::ggplot(
      top,
      ggplot2::aes(x = assembly_id_ranked, y = score)
    ) +
      ggplot2::geom_col(width = 0.75) +
      ggplot2::coord_flip() +
      ggplot2::labs(
        x = NULL,
        y = "Robust score",
        title = if (final) "Final LongStitch candidate ranking" else "LongStitch preselection ranking"
      ) +
      theme_selection()

    save_selection_plot(outdir, if (final) "ranking_ggplot" else "preselection_ranking_ggplot", p_rank, 10, 6)
  }

  # 2. k x w heatmap -----------------------------------------------------------
  if (nrow(plot_df) > 0 && all(c("k", "w") %in% names(plot_df))) {
    plot_df$k_num <- to_num(plot_df$k)
    plot_df$w_num <- to_num(plot_df$w)

    if (any(is.finite(plot_df$k_num)) && any(is.finite(plot_df$w_num))) {
      plot_df$k_factor <- factor(plot_df$k_num, levels = sort(unique(plot_df$k_num)))
      plot_df$w_factor <- factor(plot_df$w_num, levels = sort(unique(plot_df$w_num)))

      p_heat <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(x = w_factor, y = k_factor, fill = score)
      ) +
        ggplot2::geom_tile(color = "white", linewidth = 0.35) +
        ggplot2::geom_text(ggplot2::aes(label = round(score, 2)), size = 2.8) +
        ggplot2::labs(
          x = "LongStitch window size (w)",
          y = "LongStitch k-mer size (k)",
          fill = "Score",
          title = if (final) "Final LongStitch score heatmap" else "Preselection score heatmap"
        ) +
        theme_selection()

      save_selection_plot(outdir, if (final) "score_heatmap" else "preselection_score_heatmap", p_heat, 10, 7)
    }
  }

  # 3. N50 vs Ns ---------------------------------------------------------------
  if (nrow(plot_df) > 0 && all(c("n50", "ns_per_100kbp") %in% names(plot_df))) {
    p_n50_ns <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(
        x = n50,
        y = ns_per_100kbp,
        size = score
      )
    ) +
      ggplot2::geom_point(alpha = 0.75) +
      ggplot2::labs(
        x = "N50",
        y = "Ns per 100 kbp",
        size = "Score",
        title = "Contiguity versus gap burden"
      ) +
      theme_selection()

    save_selection_plot(outdir, if (final) "n50_vs_ns" else "preselection_n50_vs_ns", p_n50_ns, 8, 6)
  }

  # 4. Total assembly length vs N50 -------------------------------------------
  if (nrow(plot_df) > 0 && all(c("total_length", "n50") %in% names(plot_df))) {
    p_size <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(
        x = total_length,
        y = n50,
        size = score
      )
    ) +
      ggplot2::geom_point(alpha = 0.75) +
      ggplot2::labs(
        x = "Total assembly length",
        y = "N50",
        size = "Score",
        title = "Assembly size versus contiguity"
      ) +
      theme_selection()

    if (is.finite(genome_size_est) && genome_size_est > 0) {
      ymax <- max(plot_df$n50, na.rm = TRUE)

      if (is.finite(ymax)) {
        p_size <- p_size +
          ggplot2::geom_vline(xintercept = genome_size_est, linetype = "dashed") +
          ggplot2::annotate(
            "text",
            x = genome_size_est,
            y = ymax,
            label = "Genome-size estimate",
            angle = 90,
            vjust = -0.4,
            hjust = 1
          )
      }
    }

    save_selection_plot(outdir, if (final) "total_length_vs_n50" else "preselection_total_length_vs_n50", p_size, 8, 6)
  }

  # 5. Merqury QV vs completeness ---------------------------------------------
  if (
    final &&
    nrow(plot_df) > 0 &&
    all(c("merqury_qv", "merqury_completeness") %in% names(plot_df))
  ) {
    p_merqury <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(
        x = merqury_completeness,
        y = merqury_qv,
        size = score
      )
    ) +
      ggplot2::geom_point(alpha = 0.75) +
      ggplot2::labs(
        x = "Merqury k-mer completeness (%)",
        y = "Merqury QV",
        size = "Score",
        title = "Merqury completeness versus consensus accuracy"
      ) +
      theme_selection()

    save_selection_plot(outdir, "merqury_qv_vs_completeness", p_merqury, 8, 6)
  }

  # 6. yak QV vs score ---------------------------------------------------------
  if (nrow(plot_df) > 0 && "yak_qv" %in% names(plot_df)) {
    p_yak <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(
        x = yak_qv,
        y = score
      )
    ) +
      ggplot2::geom_point(alpha = 0.75, size = 3) +
      ggplot2::labs(
        x = "yak QV",
        y = "Robust score",
        title = "yak QV versus robust score"
      ) +
      theme_selection()

    save_selection_plot(outdir, if (final) "yak_qv_vs_score" else "preselection_yak_qv_vs_score", p_yak, 8, 6)
  }

  # 7. Score component plot ----------------------------------------------------
  component_cols <- intersect(
    c(
      "score_n50",
      "score_l50",
      "score_scaffold_count",
      "score_ns",
      "score_total_length",
      "score_largest_scaffold",
      "score_aun",
      "score_yak_qv",
      "score_merqury_qv",
      "score_merqury_completeness"
    ),
    names(top)
  )

  if (nrow(top) > 0 && length(component_cols) > 0) {
    comp_long <- data.frame()

    for (cc in component_cols) {
      comp_long <- rbind(
        comp_long,
        data.frame(
          assembly_id = top$assembly_id,
          component = cc,
          contribution = to_num(top[[cc]]),
          stringsAsFactors = FALSE
        )
      )
    }

    comp_long <- comp_long[is.finite(comp_long$contribution), , drop = FALSE]

    if (nrow(comp_long) > 0) {
      p_comp <- ggplot2::ggplot(
        comp_long,
        ggplot2::aes(
          x = assembly_id,
          y = contribution,
          fill = component
        )
      ) +
        ggplot2::geom_col() +
        ggplot2::coord_flip() +
        ggplot2::labs(
          x = NULL,
          y = "Weighted score contribution",
          fill = "Component",
          title = "Metric contributions to top candidate scores"
        ) +
        theme_selection()

      save_selection_plot(outdir, if (final) "score_contributions" else "preselection_score_contributions", p_comp, 11, 7)
    }
  }

  # 8. Metric distribution panel ----------------------------------------------
  metric_cols <- intersect(
    c(
      "n50",
      "l50",
      "n90",
      "l90",
      "ng50",
      "lg50",
      "largest_scaffold",
      "scaffold_count",
      "total_length",
      "ns_per_100kbp",
      "gc_percent",
      "auN",
      "auNG",
      "yak_qv",
      "merqury_qv",
      "merqury_completeness",
      "score"
    ),
    names(plot_df)
  )

  if (nrow(plot_df) > 0 && length(metric_cols) > 0) {
    metric_long <- data.frame()

    for (mc in metric_cols) {
      metric_long <- rbind(
        metric_long,
        data.frame(
          assembly_id = plot_df$assembly_id,
          metric = mc,
          value = to_num(plot_df[[mc]]),
          stringsAsFactors = FALSE
        )
      )
    }

    metric_long <- metric_long[is.finite(metric_long$value), , drop = FALSE]

    if (nrow(metric_long) > 0) {
      p_dist <- ggplot2::ggplot(
        metric_long,
        ggplot2::aes(x = value)
      ) +
        ggplot2::geom_histogram(bins = 25, alpha = 0.8) +
        ggplot2::facet_wrap(~ metric, scales = "free", ncol = 3) +
        ggplot2::labs(
          x = NULL,
          y = "Candidate count",
          title = "Distribution of candidate assembly metrics"
        ) +
        theme_selection(base_size = 11)

      save_selection_plot(outdir, if (final) "metric_distributions" else "preselection_metric_distributions", p_dist, 12, 8)
    }
  }
}

# ==============================================================================
# Parse command-line arguments
# ==============================================================================

manifest_path <- arg("--candidate_manifest")
quast_path <- arg("--quast")
yak_path <- arg("--yak", "")
merqury_dir <- arg("--merqury_dir", "")
outdir <- arg("--outdir")
top_n <- as.integer(arg("--top_n", "10"))
final <- identical(toupper(arg("--final", "FALSE")), "TRUE")
genome_size_est <- suppressWarnings(as.numeric(arg("--genome_size_est", "")))
save_svg <- identical(toupper(arg("--save_svg", "FALSE")), "TRUE")
plot_dpi <- suppressWarnings(as.numeric(arg("--plot_dpi", "600")))

if (is.na(manifest_path) || !nzchar(manifest_path)) stop("Missing --candidate_manifest")
if (is.na(quast_path) || !nzchar(quast_path)) stop("Missing --quast")
if (is.na(outdir) || !nzchar(outdir)) stop("Missing --outdir")
if (!is.finite(top_n) || top_n < 1) top_n <- 10
if (!is.finite(plot_dpi) || plot_dpi < 1) plot_dpi <- 600

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Read candidate manifest
# ==============================================================================

manifest <- read_tsv(manifest_path)

if (nrow(manifest) == 0) {
  stop("Candidate manifest is empty: ", manifest_path)
}

required_manifest <- c("assembly_id", "path")
missing_manifest <- setdiff(required_manifest, names(manifest))

if (length(missing_manifest) > 0) {
  stop(
    "Candidate manifest is missing required column(s): ",
    paste(missing_manifest, collapse = ", ")
  )
}

manifest$assembly_id_norm <- normalize_id(manifest$assembly_id)

# ==============================================================================
# Read QUAST metrics
# ==============================================================================

metrics <- long_metrics(quast_path, manifest)

# ==============================================================================
# Add yak metrics if available
# ==============================================================================

yak <- read_tsv(yak_path)

if (nrow(yak) > 0) {
  if (!"assembly_id_norm" %in% names(yak)) {
    if ("assembly_id" %in% names(yak)) {
      yak$assembly_id_norm <- normalize_id(yak$assembly_id)
    }
  }

  keep <- intersect(names(yak), c("assembly_id_norm", "yak_qv"))

  if (all(c("assembly_id_norm", "yak_qv") %in% keep)) {
    metrics <- merge(
      metrics,
      yak[, keep, drop = FALSE],
      by = "assembly_id_norm",
      all.x = TRUE
    )
  }
}

# ==============================================================================
# Add Merqury metrics if final mode
# ==============================================================================

if (final) {
  mq_files <- list.files(
    merqury_dir,
    pattern = "merqury_metrics.tsv$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(mq_files) > 0) {
    mq <- do.call(rbind, lapply(mq_files, read_tsv))

    if (nrow(mq) > 0) {
      if (!"assembly_id_norm" %in% names(mq)) {
        if ("assembly_id" %in% names(mq)) {
          mq$assembly_id_norm <- normalize_id(mq$assembly_id)
        }
      }

      keep <- intersect(
        names(mq),
        c("assembly_id_norm", "merqury_qv", "merqury_completeness")
      )

      if ("assembly_id_norm" %in% keep) {
        mq <- mq[, keep, drop = FALSE]
        metrics <- merge(metrics, mq, by = "assembly_id_norm", all.x = TRUE)
      }
    }
  }
}

# ==============================================================================
# Score candidates
# ==============================================================================

metrics$score <- 0

metrics$score_n50 <- 0.15 * score_col(metrics$n50, 1)
metrics$score_l50 <- 0.10 * score_col(metrics$l50, -1)
metrics$score_scaffold_count <- 0.10 * score_col(metrics$scaffold_count, -1)
metrics$score_ns <- 0.10 * score_col(metrics$ns_per_100kbp, -1)
metrics$score_largest_scaffold <- 0.05 * score_col(metrics$largest_scaffold, 1)

metrics$score <- metrics$score +
  metrics$score_n50 +
  metrics$score_l50 +
  metrics$score_scaffold_count +
  metrics$score_ns +
  metrics$score_largest_scaffold

# Optional auN/auNG contribution
if ("auNG" %in% names(metrics) && any(is.finite(to_num(metrics$auNG)))) {
  metrics$score_aun <- 0.05 * score_col(metrics$auNG, 1)
  metrics$score <- metrics$score + metrics$score_aun
} else if ("auN" %in% names(metrics) && any(is.finite(to_num(metrics$auN)))) {
  metrics$score_aun <- 0.05 * score_col(metrics$auN, 1)
  metrics$score <- metrics$score + metrics$score_aun
} else {
  metrics$score_aun <- 0
}

# Optional genome-size prior
metrics$score_total_length <- 0

if (any(is.finite(to_num(metrics$total_length)))) {
  if (is.finite(genome_size_est) && genome_size_est > 0) {
    metrics$total_length_penalty <- abs(to_num(metrics$total_length) - genome_size_est) / genome_size_est
    metrics$score_total_length <- 0.10 * score_col(metrics$total_length_penalty, -1)
    metrics$score <- metrics$score + metrics$score_total_length
  }
}

# Optional Merqury contributions
if ("merqury_qv" %in% names(metrics)) {
  metrics$score_merqury_qv <- 0.25 * score_col(metrics$merqury_qv, 1)
  metrics$score <- metrics$score + metrics$score_merqury_qv
} else {
  metrics$score_merqury_qv <- 0
}

if ("merqury_completeness" %in% names(metrics)) {
  metrics$score_merqury_completeness <- 0.20 * score_col(metrics$merqury_completeness, 1)
  metrics$score <- metrics$score + metrics$score_merqury_completeness
} else {
  metrics$score_merqury_completeness <- 0
}

# Optional yak contribution
if ("yak_qv" %in% names(metrics)) {
  metrics$score_yak_qv <- 0.05 * score_col(metrics$yak_qv, 1)
  metrics$score <- metrics$score + metrics$score_yak_qv
} else {
  metrics$score_yak_qv <- 0
}

# ==============================================================================
# Failure flags and hard penalties
# ==============================================================================

metrics$failed <- FALSE

metrics$failed <- metrics$failed |
  !file.exists(metrics$path) |
  is.na(metrics$n50) |
  is.na(metrics$scaffold_count)

# Mark severe gap outliers as soft-risk but not failed
metrics$gap_outlier <- FALSE

if (any(is.finite(to_num(metrics$ns_per_100kbp)))) {
  ns <- to_num(metrics$ns_per_100kbp)
  q95 <- stats::quantile(ns, 0.95, na.rm = TRUE)
  metrics$gap_outlier <- is.finite(ns) & ns > q95
}

metrics$score[metrics$gap_outlier] <- metrics$score[metrics$gap_outlier] - 1
metrics$score[metrics$failed] <- -Inf

# ==============================================================================
# Sort, write outputs, and select top candidates
# ==============================================================================

metrics <- metrics[
  order(
    -metrics$score,
    -to_num(metrics$n50),
    to_num(metrics$l50),
    to_num(metrics$scaffold_count)
  ),
  ,
  drop = FALSE
]

score_file <- if (final) "longstitch_scores.tsv" else "preselection_scores.tsv"

write.table(
  metrics,
  file.path(outdir, score_file),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

top <- head(metrics[is.finite(metrics$score), , drop = FALSE], top_n)

writeLines(top$path, file.path(outdir, "top_candidates.txt"))

if (nrow(top) > 0) {
  writeLines(top$path[[1]], file.path(outdir, "best_assembly.path"))
  writeLines(top$assembly_id[[1]], file.path(outdir, "best_assembly.txt"))
} else {
  warning("No valid top candidates were selected")
}

# ==============================================================================
# Plot outputs
# ==============================================================================

plot_selection_results(
  metrics = metrics,
  top = top,
  outdir = outdir,
  final = final,
  genome_size_est = genome_size_est
)

message("[INFO] Selection complete")
message("[INFO] Mode: ", if (final) "final" else "preselection")
message("[INFO] Candidates ranked: ", nrow(metrics))
message("[INFO] Valid top candidates: ", nrow(top))
if (nrow(top) > 0) {
  message("[INFO] Best assembly: ", top$assembly_id[[1]])
  message("[INFO] Best path: ", top$path[[1]])
}
RSCRIPT

  if ! scaffolding_checkpoint_ready 05G.preselect; then
    env_run "${ENV_R}" Rscript "${selector_r}" \
      --candidate_manifest "${candidate_manifest}" \
      --quast "${scaff_quast_dir}/report.tsv" \
      --yak "${yak_metrics}" \
      --merqury_dir "${scaff_merqury_dir}" \
      --outdir "${scaff_selection_dir}" \
      --top_n "${MERQURY_TOP_N}" \
      --genome_size_est "${GENOME_SIZE_EST}" \
      --save_svg "${STEP5_SAVE_SVG}" \
      --plot_dpi "${STEP5_PLOT_DPI}" \
      --final FALSE \
      > "${scaff_logs_dir}/preselect_longstitch_candidates.log" 2>&1 \
      || error "Step 5 preselection R script failed"
    [[ -s "${scaff_selection_dir}/top_candidates.txt" ]] || error "Step 5 preselection produced no top candidates"
    scaffolding_mark_checkpoint 05G.preselect
  fi

  # -------------------- STEP 5H: Merqury validation on top candidates --------------------
  run_merqury_candidate() {
      local asm="${1:-}" asm_id out_dir out_prefix merqury_bin qv_value comp_value candidate_qv candidate_qv_value candidate_comp candidate_comp_value
      [[ -s "${asm}" ]] || { echo "[ERROR] Merqury candidate missing: ${asm}" >&2; return 1; }
      asm_id="$(basename "${asm}")"
      asm_id="${asm_id%.fa}"; asm_id="${asm_id%.fasta}"; asm_id="${asm_id%.fna}"
      out_dir="${scaff_merqury_dir}/${asm_id}"
      out_prefix="${out_dir}/${asm_id}"
      mkdir -p "${out_dir}"
      if [[ -s "${out_dir}/merqury_metrics.tsv" ]]; then
        echo "[INFO] Merqury metrics exist for ${asm_id}; skipping"
        return 0
      fi
      export OMP_NUM_THREADS="${nthreads}"
      if [[ -n "${MERQURY:-}" && -x "${MERQURY}/merqury.sh" ]]; then
        merqury_bin="${MERQURY}/merqury.sh"
      else
        merqury_bin="merqury.sh"
      fi
      "${merqury_bin}" "${read_meryl_db}" "${asm}" "${out_prefix}" \
        > "${out_dir}/merqury.log" 2>&1 \
        || return 1
      qv_value=""
      while IFS= read -r candidate_qv; do
        candidate_qv_value="$(awk 'NF >= 4 && $4 ~ /^[-+]?[0-9.]+([eE][-+]?[0-9]+)?$/ {print $4; exit}' "${candidate_qv}" 2>/dev/null || true)"
        if [[ -n "${candidate_qv_value}" ]]; then
          qv_value="${candidate_qv_value}"
          break
        fi
      done < <(find "${out_dir}" -type f \( -name "*.qv" -o -name "*.qv.txt" -o -name "*qv*.tsv" \) -size +0c 2>/dev/null | sort)

      comp_value=""
      while IFS= read -r candidate_comp; do
        candidate_comp_value="$(awk 'NF >= 5 && $NF ~ /^[-+]?[0-9.]+([eE][-+]?[0-9]+)?$/ {v=$NF} END{print v}' "${candidate_comp}" 2>/dev/null || true)"
        if [[ -n "${candidate_comp_value}" ]]; then
          comp_value="${candidate_comp_value}"
          break
        fi
      done < <(find "${out_dir}" -type f -name "*.completeness.stats" -size +0c 2>/dev/null | sort)

      [[ -n "${qv_value}" ]] || echo "[WARN] Merqury QV file missing or unparsable for ${asm_id}; recording NA" >&2
      [[ -n "${comp_value}" ]] || echo "[WARN] Merqury completeness file missing or unparsable for ${asm_id}; recording NA" >&2
      {
        printf "assembly_id\tmerqury_qv\tmerqury_completeness\tpath\n"
        printf "%s\t%s\t%s\t%s\n" "${asm_id}" "${qv_value:-NA}" "${comp_value:-NA}" "${asm}"
      } > "${out_dir}/merqury_metrics.tsv"
    }
  export -f run_merqury_candidate
  export scaff_merqury_dir read_meryl_db MERQURY nthreads

  if ! scaffolding_checkpoint_ready 05H.merqury; then
    activate_env "${ENV_MERQURY}"
    parallel --jobs "${njobs}" --load 80% --memfree 64G --tmpdir "${tmpdir}" \
      --env run_merqury_candidate --env scaff_merqury_dir --env read_meryl_db --env MERQURY --env nthreads \
      run_merqury_candidate {} :::: "${scaff_selection_dir}/top_candidates.txt" \
      2>&1 | tee -a "${scaff_logs_dir}/merqury.parallel.log"
    merqury_rc=${PIPESTATUS[0]:-0}
    (( merqury_rc == 0 )) || error "Merqury failed for one or more Step 5 candidates"
    deactivate_env
    scaffolding_mark_checkpoint 05H.merqury
  fi

  # -------------------- STEP 5I: Robust R selection and plotting --------------------
  if ! scaffolding_checkpoint_ready 05I.selection; then
    env_run "${ENV_R}" Rscript "${selector_r}" \
      --candidate_manifest "${candidate_manifest}" \
      --quast "${scaff_quast_dir}/report.tsv" \
      --yak "${yak_metrics}" \
      --merqury_dir "${scaff_merqury_dir}" \
      --outdir "${scaff_selection_dir}" \
      --top_n "${MERQURY_TOP_N}" \
      --genome_size_est "${GENOME_SIZE_EST}" \
      --save_svg "${STEP5_SAVE_SVG}" \
      --plot_dpi "${STEP5_PLOT_DPI}" \
      --final TRUE \
      > "${scaff_logs_dir}/select_longstitch_best.log" 2>&1 \
      || error "Step 5 final R selector failed"
    [[ -s "${scaff_selection_dir}/best_assembly.path" ]] || error "Step 5 final selector did not write best_assembly.path"
    scaffolding_mark_checkpoint 05I.selection
  fi

  # -------------------- STEP 5J: Final scaffold FASTA creation --------------------
  best_path="$(head -n1 "${scaff_selection_dir}/best_assembly.path")"
  [[ -s "${best_path}" ]] || error "Selected Step 5 assembly missing: ${best_path}"
  final_scaff="${dir5}/${base}.scaffolds.fa"
  if [[ ! -s "${final_scaff}" || "${best_path}" -nt "${final_scaff}" ]]; then
    awk '/^>/{print ">scaff" ++i; next}{print}' "${best_path}" > "${final_scaff}" \
      || error "Failed to canonicalize selected scaffold FASTA"
  fi
  printf '%s\n' "${best_path}" > "${dir5}/best_assembly.path"
  scaffolding_mark_checkpoint 05J.final_scaffolds
  echo "[INFO] Step 5 complete. Selected scaffolds: ${final_scaff} (source: ${best_path})"
fi

# ==================================================================================================
# STEP 6: DeepVariant post-polishing on Step 5 scaffolds (+ FCS decontamination)
# Consumes : ${dir5}/${base}.scaffolds.fa  and ${pacbio_hifi}
# Produces : ${dir6}/${base}.round{1,2,3}.polished.fa
#            ${dir6}/final.polished.decontam.fa
#            ${dir6}/${species_name}.fasta (final published assembly with canonical headers)
# ==================================================================================================
if [[ "${run_final_asm:-FALSE}" == "TRUE" ]]; then
  echo "[INFO] Step 6: DeepVariant post-polishing + decontamination — starting"
  mkdir -p "${dir6}/logs" "${dir6}/checkpoints" "${dir6}/gx" "${dir6}/adaptor"

  # -------------------- Inputs --------------------
  initial_assembly="${dir5}/${base}.scaffolds.fa"
  [[ -s "${initial_assembly}" ]] || error "Step 6 input missing: ${initial_assembly}"
  [[ -s "${pacbio_hifi}" ]]      || error "PacBio HiFi input not found: ${pacbio_hifi}"

  # -------------------- Helper: ensure dict & fai --------------------
  ensure_ref_index() {
    local ref="$1"
    local dict="${ref}.dict"
    if [[ ! -s "${dict}" ]]; then
      echo "[INFO] Creating sequence dictionary for ${ref}"
      singularity exec -B "${working_dir}:${working_dir}" "${sif_gatk4}" \
        gatk CreateSequenceDictionary -R "${ref}" -O "${dict}" \
        >>"${dir6}/logs/gatk.dict.log" 2>&1 || error "CreateSequenceDictionary failed for ${ref}"
    fi
    if [[ ! -s "${ref}.fai" ]]; then
      echo "[INFO] Indexing ${ref} with samtools faidx"
      "${samtools_bin}" faidx "${ref}" >>"${dir6}/logs/samtools.faidx.log" 2>&1 \
        || error "samtools faidx failed for ${ref}"
    fi
  }

  # -------------------- Rounds 1..3 --------------------
  for round in 1 2 3; do
    echo "----------------------------------------------"
    echo "[INFO] Step 6: Round ${round}"
    echo "----------------------------------------------"

    ck_prefix="${dir6}/checkpoints/round${round}"
    asm_in="${initial_assembly}"
    if (( round > 1 )); then
      prev=$((round - 1))
      asm_in="${dir6}/${base}.round${prev}.polished.fa"
    fi
    ensure_ref_index "${asm_in}"

    sorted_bam="${dir6}/${base}.round${round}.sorted.bam"
    markdup_bam="${dir6}/${base}.round${round}.markdup.bam"
    dv_vcf="${dir6}/${base}.round${round}.dv.vcf.gz"
    dv_vcf_hom="${dir6}/${base}.round${round}.dv.hom.vcf.gz"
    asm_out="${dir6}/${base}.round${round}.polished.fa"

    # 6.1 Map HiFi → sorted BAM
    if [[ ! -f "${ck_prefix}.mapping.done" ]]; then
      echo "[INFO] Round ${round}: Mapping HiFi reads to ${asm_in} → sorted BAM"
      "${minimap2_bin}" -H -Y --secondary=no -a -x map-hifi -K 16G -t "$((nthreads * njobs))" -R "@RG\\tID:${base}.HiFi\\tSM:${base}\\tPL:PACBIO\\tLB:SMRTbell" "${asm_in}" "${pacbio_hifi}" \
      | "${samtools_bin}" sort -@ "$((nthreads * njobs))" -m 4G -l 9 -O BAM -o "${sorted_bam}" -T "${dir6}/round${round}.tmp" --write-index \
          >>"${dir6}/logs/round${round}.map_sort.log" 2>&1 \
        || error "Round ${round}: minimap2|samtools sort failed"
      touch "${ck_prefix}.mapping.done"
    else
      echo "[INFO] Round ${round}: using existing sorted BAM"
    fi

    # 6.2 MarkDuplicatesSpark (BAM in/out, create CSI index)
    if [[ ! -f "${ck_prefix}.markdup.done" ]]; then
	      echo "[INFO] Round ${round}: MarkDuplicatesSpark"
	      spark_threads="$((nthreads * njobs))"
		      markdup_shard_tmp="${tmpdir}/gatk_markdup.step6.round${round}.shards"
		      markdup_tmp="${tmpdir}/gatk_markdup.step6.round${round}.tmp"
		      rm -rf "${markdup_shard_tmp}" "${markdup_tmp}"
		      mkdir -p "${markdup_tmp}"
		      rm -f "${markdup_bam}" "${markdup_bam}.bai" "${markdup_bam}.csi"
      gatk_spark \
        --java-options "-Xms128G -Xmx256G -Djava.io.tmpdir=${tmpdir}" \
        MarkDuplicatesSpark \
          --spark-master "local[${spark_threads}]" \
          --conf "spark.driver.host=127.0.0.1" \
          --conf "spark.driver.bindAddress=127.0.0.1" \
          --conf "spark.local.dir=${tmpdir}" \
          --input "${sorted_bam}" \
          --metrics-file "${sorted_bam}.metrics" \
          --output "${markdup_bam}" \
          --reference "${asm_in}" \
	          --remove-all-duplicates true \
	          --create-output-bam-index true \
          --create-output-bam-splitting-index false \
          --read-name-regex null \
          --read-validation-stringency SILENT \
	          --tmp-dir "${markdup_tmp}" \
	          --output-shard-tmp-dir "${markdup_shard_tmp}" \
	        >>"${dir6}/logs/round${round}.markduplicates.log" 2>&1 \
	        || error "Round ${round}: MarkDuplicates failed"
	      "${samtools_bin}" index -b -@ "${spark_threads}" "${markdup_bam}" \
	        >>"${dir6}/logs/round${round}.markduplicates.log" 2>&1 \
	        || error "Round ${round}: BAI indexing of MarkDuplicates output failed"
	      rm -rf "${markdup_shard_tmp}" "${markdup_tmp}"
      touch "${ck_prefix}.markdup.done"
    else
      echo "[INFO] Round ${round}: using existing markdup BAM"
    fi

    # 6.3 DeepVariant
    if [[ ! -f "${ck_prefix}.deepvariant.done" ]]; then
      echo "[INFO] Round ${round}: DeepVariant calling"
      dv_tmp=$(mktemp -d -p "${dir6}" dv.tmp.XXXXXXXX) || error "mktemp failed"
      singularity exec -B "${working_dir}:${working_dir}" "${sif_deepvariant}" \
        /opt/deepvariant/bin/run_deepvariant \
          --model_type=PACBIO \
          --ref="${asm_in}" \
          --reads="${markdup_bam}" \
          --output_vcf="${dv_vcf}" \
          --num_shards="$((nthreads * njobs))" \
          --postprocess_cpus="$((nthreads * njobs))" \
          --intermediate_results_dir="${dv_tmp}" \
          --vcf_stats_report true \
        >>"${dir6}/logs/round${round}.deepvariant.log" 2>&1 || dv_fail=1
      rm -rf "${dv_tmp}" || true
      [[ -z "${dv_fail:-}" ]] || error "Round ${round}: DeepVariant failed"
      unset dv_fail
      touch "${ck_prefix}.deepvariant.done"
    else
      echo "[INFO] Round ${round}: using existing DeepVariant VCF"
    fi

    # 6.4 Filter VCF (PASS, hom) → consensus
    if [[ ! -f "${ck_prefix}.consensus.done" ]]; then
      echo "[INFO] Round ${round}: bcftools filter + consensus"
      "${bcftools_bin}" view -l 0 -f 'PASS' -i 'GT="1/1"' "${dv_vcf}" \
        | "${bcftools_bin}" view --threads "$((nthreads * njobs))" -Oz9 -o "${dv_vcf_hom}" --write-index \
        >>"${dir6}/logs/round${round}.bcftools.log" 2>&1 \
        || error "Round ${round}: bcftools filtering failed"

      "${bcftools_bin}" consensus -f "${asm_in}" -o "${asm_out}" "${dv_vcf_hom}" \
        >>"${dir6}/logs/round${round}.consensus.log" 2>&1 \
        || error "Round ${round}: consensus failed"
      touch "${ck_prefix}.consensus.done"
    else
      echo "[INFO] Round ${round}: consensus already complete"
    fi

    # 6.5 Cleanup intermediates
    if [[ -f "${ck_prefix}.consensus.done" ]]; then
      rm -f "${sorted_bam}"* "${markdup_bam}"* "${dv_vcf}"* "${dv_vcf_hom}"* || true
      rm -f "${asm_in}.fai" "${asm_in}.dict" || true
    fi

    echo "[INFO] Round ${round} complete: ${asm_out}"
  done

  ###############################################
  # Step 6 — Final Decontamination (mandatory)
  ###############################################
  echo "[INFO] Step 6: Final Decontamination — starting"

  # Input for decontam
  final_post="${dir6}/${base}.round3.polished.fa"
  [[ -s "${final_post}" ]] || error "Step 6 decontam: missing input ${final_post}"
  fcs_adaptor_input="${dir6}/${base}.round3.polished.fcs_input.fa"

  # Outputs
  decontam_out="${dir6}/final.polished.decontam.fa"
  canonical_out="${dir6}/final.polished.decontam.canonical.fa"
  species_fasta="${dir6}/${species_name}.fasta"

  # Checkpoints (kept as requested)
  ck_adaptor="${dir6}/checkpoints/step6.adaptor.done"
  ck_gx_screen="${dir6}/checkpoints/step6.gx_screen.done"
  ck_gx_clean="${dir6}/checkpoints/step6.gx_clean.done"
  ck_publish="${dir6}/checkpoints/step6.publish.done"

  # --- Tool & container checks ---
  command -v singularity >/dev/null 2>&1 || error "STEP 6: 'singularity' not on PATH"
  [[ -f "${sif_fcs_adaptor}" ]] || error "STEP 6: sif_fcs_adaptor missing: ${sif_fcs_adaptor}"
  [[ -f "${sif_fcs_gx}"      ]] || error "STEP 6: sif_fcs_gx missing: ${sif_fcs_gx}"
  [[ -n "${ADAPTOR_RUNNER:-}" ]] || error "STEP 6: ADAPTOR_RUNNER is not set"
  [[ -d "${GX_DB:-}" ]] || error "STEP 6: GX_DB directory missing"

  # ------------------ 6A) FCS-Adaptor ------------------
  adaptor_pref="${dir6}/adaptor"
  mkdir -p "${adaptor_pref}"

  if [[ ! -s "${fcs_adaptor_input}" || ! -s "${dir6}/logs/fcs_input.prepare.log" || "${final_post}" -nt "${fcs_adaptor_input}" ]]; then
    echo "[INFO] Step 6: preparing FASTA for FCS validation"
    awk -v min_len="${FCS_ADAPTOR_MIN_SEQ_LEN}" '
      BEGIN { header = ""; seq = ""; out = 0; short = 0; bad = 0 }
      function emit() {
        if (header == "") { return }
        gsub(/[[:space:]]/, "", seq)
        seq = toupper(seq)
        bad += gsub(/[^ACGTN]/, "N", seq)
        if (length(seq) < min_len) { short++; return }
        out++
        print header
        for (i = 1; i <= length(seq); i += 80) {
          print substr(seq, i, 80)
        }
      }
      /^>/ { emit(); header = $0; seq = ""; next }
      { seq = seq $0 }
      END {
        emit()
        if (out == 0) {
          print "[ERROR] No FASTA records >= " min_len " bp found" > "/dev/stderr"
          exit 2
        }
        if (bad > 0) {
          print "[WARN] Replaced " bad " non-ACGTN FASTA characters with N" > "/dev/stderr"
        }
        print "[INFO] FCS FASTA records written: " out > "/dev/stderr"
        print "[INFO] FCS FASTA records shorter than " min_len " bp skipped: " short > "/dev/stderr"
      }
    ' "${final_post}" > "${fcs_adaptor_input}" \
      2> "${dir6}/logs/fcs_input.prepare.log" \
      || error "Step 6: failed to prepare FASTA for FCS-adaptor"
  fi

  if [[ -f "${ck_adaptor}" && ( -s "${adaptor_pref}/cleaned_sequences/$(basename "${fcs_adaptor_input}")" || -s "${adaptor_pref}/cleaned_sequences/$(basename "${final_post}")" || -s "${adaptor_pref}/clean.fasta" ) ]]; then
    echo "[INFO] Step 6: FCS-Adaptor checkpoint present"
  else
    echo "[INFO] Step 6: FCS-Adaptor on ${fcs_adaptor_input}"
    "${ADAPTOR_RUNNER}" \
      --fasta-input "${fcs_adaptor_input}" \
      --output-dir  "${adaptor_pref}" \
      --container-engine singularity \
      --image "${sif_fcs_adaptor}" \
      --euk \
      --debug \
      2>&1 | tee -a "${dir6}/logs/fcs_adaptor.run.log"
    touch "${ck_adaptor}"
  fi

  # Choose adaptor-cleaned FASTA (robust fallback)
  adaptor_after="${adaptor_pref}/cleaned_sequences/$(basename "${fcs_adaptor_input}")"
  if [[ ! -s "${adaptor_after}" && -s "${adaptor_pref}/cleaned_sequences/$(basename "${final_post}")" ]]; then
    adaptor_after="${adaptor_pref}/cleaned_sequences/$(basename "${final_post}")"
  fi
  if [[ ! -s "${adaptor_after}" && -s "${adaptor_pref}/clean.fasta" ]]; then
    adaptor_after="${adaptor_pref}/clean.fasta"
  fi
  if [[ ! -s "${adaptor_after}" ]]; then
    echo "[WARN] [STEP 6A] Adaptor-cleaned FASTA not found; falling back to ${fcs_adaptor_input}"
    adaptor_after="${fcs_adaptor_input}"
  fi

  # ------------------ 6B) Ensure GX DB present (sync if incomplete) ------------------
  have_gxdb_step6() {
    local d="${GX_DB}"
    local -a req=(all.assemblies.tsv all.blast_div.tsv.gz all.gxi all.gxs all.manifest all.meta.jsonl all.README.txt all.seq_info.tsv.gz all.taxa.tsv)
    for f in "${req[@]}"; do [[ -s "${d}/$f" ]] || return 1; done
    [[ -d "${d}/gxdb.in_progress" ]] && return 1
    return 0
  }

  mkdir -p "${GX_DB}"
  [[ -w "${GX_DB}" ]] || error "STEP 6: GX_DB not writable: ${GX_DB}"

  GX_BIND_MODE="rw"
  for gx_lockfile in "${GX_DB}.lockfile" "${GX_DB}/gxdb.lockfile"; do
    if [[ -f "${gx_lockfile}" ]] && ! head -n 1 "${gx_lockfile}" 2>/dev/null | grep -Eq '^[0-9]+,[^,[:space:]]+$'; then
      echo "[WARN] [STEP 6B] Removing malformed GX DB lockfile: ${gx_lockfile}"
      rm -f "${gx_lockfile}" || error "STEP 6: failed to remove malformed lockfile: ${gx_lockfile}"
    fi
  done

  if have_gxdb_step6; then
    echo "[INFO] [STEP 6B] GX DB present → ${GX_DB} (skip download)"
    [[ -d "${GX_DB}/gxdb.in_progress" ]] && { echo "[WARN] [STEP 6B] Removing stale gxdb.in_progress"; rm -rf "${GX_DB}/gxdb.in_progress"; }
    GX_BIND_MODE="ro"
  else
    echo "[INFO] [STEP 6B] Sync GX DB → ${GX_DB}"
    rm -rf "${GX_DB}/gxdb.in_progress" 2>/dev/null || true
    singularity exec \
      --bind "${GX_ROOT}:${GX_ROOT}" \
      "${sif_fcs_gx}" \
      python3 /app/bin/sync_files get \
        --mft "${GX_MANIFEST}" \
        --dir "${GX_DB}" \
      2>&1 | tee -a "${dir6}/logs/fcs_gx.sync.log"
    have_gxdb_step6 || error "[STEP 6B] GX DB still incomplete after sync -> ${GX_DB}"
    GX_BIND_MODE="ro"
  fi

  # ------------------ 6C) GX screen ------------------
  if [[ -f "${ck_gx_screen}" && -n "$(find "${dir6}/gx" -maxdepth 1 -type f -name '*fcs_gx_report.txt' -print -quit 2>/dev/null || true)" ]]; then
    echo "[INFO] Step 6: FCS-GX screen checkpoint present"
  else
    echo "[INFO] Step 6: FCS-GX screen"
    singularity exec \
      --bind "${GX_DB}:/app/db/gxdb:${GX_BIND_MODE}" \
      --bind "$(dirname "${adaptor_after}"):/sample-volume" \
      --bind "${dir6}/gx:/output-volume" \
      "${sif_fcs_gx}" \
      python3 /app/bin/run_gx \
        --fasta "/sample-volume/$(basename "${adaptor_after}")" \
        --gx-db /app/db/gxdb \
        --out-dir /output-volume \
        ${TAXID:+--tax-id "${TAXID}"} \
        --generate-logfile True \
        --debug \
      2>&1 | tee -a "${dir6}/logs/fcs_gx.run.log"
    touch "${ck_gx_screen}"
  fi

  # ------------------ 6D) GX clean-genome ------------------
  if [[ -f "${ck_gx_clean}" && -s "${decontam_out}" ]]; then
    echo "[INFO] Step 6: FCS-GX clean-genome checkpoint present"
  else
    echo "[INFO] Step 6: FCS-GX clean-genome"
    gx_report="$(find "${dir6}/gx" -maxdepth 1 -type f -name '*fcs_gx_report.txt' -print -quit 2>/dev/null || true)"

    if [[ -z "${gx_report}" ]]; then
      echo "[INFO] [STEP 6D] No GX report found — using adaptor-cleaned (or post-polish) as clean.fasta"
      if [[ -s "${adaptor_after}" ]]; then
        ln -sfn "${adaptor_after}" "${dir6}/gx/clean.fasta"
      else
        ln -sfn "${final_post}" "${dir6}/gx/clean.fasta"
      fi
    else
      singularity exec \
        --bind "${dir6}/gx:/report-volume" \
        --bind "${dir6}/gx:/output-volume" \
        --bind "$(dirname "${adaptor_after}"):/input-volume" \
        "${sif_fcs_gx}" \
        /app/bin/gx clean-genome \
          --input "/input-volume/$(basename "${adaptor_after}")" \
          --action-report "/report-volume/$(basename "${gx_report}")" \
          --output "/output-volume/clean.fasta" \
          --contam-fasta-out "/output-volume/contam.fasta" \
          --min-seq-len 200 \
        2>&1 | tee -a "${dir6}/logs/fcs_gx.clean.log"
    fi

    # Promote to canonical decontam_out (materialize as a real file, not symlink)
    cp -f "${dir6}/gx/clean.fasta" "${decontam_out}" \
      || error "GX clean-genome did not produce ${decontam_out}"
    [[ -s "${decontam_out}" ]] || error "Decontaminated FASTA missing after GX clean-genome"
    touch "${ck_gx_clean}"
  fi

  # ------------------ Canonicalize headers + publish ------------------
  if [[ ! -f "${ck_publish}" ]]; then
    [[ -s "${decontam_out}" ]] || error "Step 6 publish: decontamination output missing: ${decontam_out}"

    awk '/^>/{print ">scaff" ++i; next}{print}' \
      "${decontam_out}" > "${canonical_out}" \
      || error "Canonicalization failed"

    cp -f "${canonical_out}" "${species_fasta}" \
      || error "Failed to publish ${species_fasta}"
    rm -f "${canonical_out}" || true
    echo "[INFO] Final canonical assembly published as: ${species_fasta}"
    touch "${ck_publish}"
  else
    echo "[INFO] Step 6: publish checkpoint present — ${species_fasta}"
  fi
fi

# ============================================================================================================================
# STEP 7: Assembly evaluation & reporting — QUAST + BUSCO + BlobToolKit + Qualimap + MultiQC
# Consumes: ${species_fasta} (final sanitized assembly). Falls back to earlier products if needed.
# Produces (in dir7):
#   - ${dir7}/quast/*
#   - ${dir7}/busco/<run>/* (+ plots)
#   - ${dir7}/blobtoolkit/*.blobdir and SVG plots under ${dir7}/blobtoolkit/plots
#   - ${dir7}/qualimap/* (bamqc on HiFi->assembly BAM)
#   - ${dir7}/multiqc_report.html
# ============================================================================================================================
if [[ "${run_stats:-FALSE}" == "TRUE" ]]; then
  echo "[INFO] Step 7: QUAST + BUSCO + BlobToolKit + Qualimap + MultiQC — starting"

  # Preferred final assembly; used when present
  : "${species_fasta:=${dir6}/${species_name}.fasta}"

  # Local report dirs (kept under dir7)
  quast_dir="${dir7}/quast"
  busco_dir="${dir7}/busco"
  btk_dir="${dir7}/blobtoolkit"
  qualimap_dir="${dir7}/qualimap"
  mkdir -p "${quast_dir}" "${busco_dir}" "${btk_dir}" "${qualimap_dir}" "${dir7}/checkpoints"
  ck_quast="${dir7}/checkpoints/step7.quast.done"
  ck_busco="${dir7}/checkpoints/step7.busco.done"
  ck_btk="${dir7}/checkpoints/step7.btk.done"
  ck7_btk_dataset="${dir7}/checkpoints/step7.btk.dataset.done"
  ck7_btk_plots="${dir7}/checkpoints/step7.btk.plots.done"
  ck_qualimap="${dir7}/checkpoints/step7.qualimap.done"
  ck_multiqc="${dir7}/checkpoints/step7.multiqc.done"

  # Shared BUSCO/BlobToolKit targets — compute once so substeps can run independently
  busco_run_name="${busco_run_name:-final.polished.busco}"
  busco_run_dir_base="${busco_dir}/${busco_run_name}"
  busco_run_dir="$(find "${busco_run_dir_base}" -maxdepth 1 -type d -name 'run_*' -print -quit 2>/dev/null || true)"
  [[ -n "${busco_run_dir}" ]] || busco_run_dir="${busco_run_dir_base}"
  busco_target="${species_fasta}"
  if [[ ! -s "${busco_target}" ]]; then
    echo "[WARN] species_fasta not found at ${species_fasta}; falling back to earlier outputs."
    for cand in \
      "${dir6}/final.polished.decontam.fa" \
      "${dir6}/final.polished.fa" \
      "${dir5}/${base}.scaffolds.fa" \
      "${dir4}/inspector_out/contig_corrected.fa" \
      "${dir3}/${base}.round3.polished.fa" \
      "${dir2}/gx/clean.fasta" \
      "${dir1}/${base}.bp.p_ctg.fa"; do
      if [[ -s "${cand}" ]]; then busco_target="${cand}"; break; fi
    done
  fi
  if [[ ! -s "${species_fasta}" && -s "${busco_target}" ]]; then
    species_fasta="${busco_target}"
  fi

  # BlobToolKit/Qualimap shared prefix (so later substeps don't depend on BTK block execution)
  bt_prefix="${btk_dir}/hifi"
  blobdir_path="${btk_dir}/blobdir"
  plots_out="${btk_dir}/plots"

  # --------------------------- 7a) QUAST ---------------------------
  echo "[INFO] Step 7a: QUAST"
  if [[ ! -f "${ck_quast}" ]]; then
    activate_env "${ENV_QUAST}"

    # Candidate inputs (existence-checked)
    declare -a quast_candidates=(
      "${dir1}/${base}.bp.p_ctg.fa"
      "${dir2}/gx/clean.fasta"
      "${dir3}/${base}.round3.polished.fa"
      "${dir4}/inspector_out/contig_corrected.fa"
      "${dir5}/${base}.scaffolds.fa"
      "${species_fasta}"                      # final sanitized assembly (preferred)
    )
    quast_inputs=()
    for f in "${quast_candidates[@]}"; do
      if [[ -s "$f" ]]; then
        quast_inputs+=("$f")
      else
        echo "[WARN] QUAST candidate missing (skipped): $f"
      fi
    done
    [[ ${#quast_inputs[@]} -gt 0 ]] || error "No inputs available for QUAST."

    quast \
      -o "${quast_dir}" \
      --threads "$(( nthreads * njobs ))" \
      --large \
      --min-contig "${QUAST_MIN_CONTIG}" \
      --report-all-metrics \
      --k-mer-stats \
      --k-mer-size "${QUAST_KMER_SIZE}" \
      --contig-thresholds "1000,5000,10000,100000,500000,1000000,5000000,10000000" \
      --x-for-Nx 90 \
      --no-snps \
      --no-icarus \
      "${quast_inputs[@]}" \
      || error "QUAST failed"
    deactivate_env
    touch "${ck_quast}"
  else
    echo "[INFO] Step 7a: QUAST checkpoint present — skipping"
  fi

  # --------------------------- 7b) BUSCO ---------------------------
  echo "[INFO] Step 7b: BUSCO (v6)"
  if [[ ! -f "${ck_busco}" ]]; then
    [[ -s "${busco_target}" ]] || error "No assembly found for BUSCO."

    # Keep all BUSCO outputs inside busco_dir (single top-level run folder)
    BUSCO_DOWNLOADS="${BUSCO_DOWNLOADS:-${busco_dir}/downloads}"
    mkdir -p "${BUSCO_DOWNLOADS}"

    export JAVA_TOOL_OPTIONS="-Xmx128g"
    env_run "${ENV_BUSCO}" busco \
      -i "${busco_target}" \
      -o "${busco_run_name}" \
      --out_path "${busco_dir}" \
      -l "${BUSCO_LINEAGE}" \
      -m genome \
      -c "$((nthreads * njobs))" \
      --download_path "${BUSCO_DOWNLOADS}" \
      --miniprot \
      --tar \
      --opt-out-run-stats \
      || error "BUSCO failed"

    # Resolve actual BUSCO run directory (BUSCO nests run_* under the run name)
    busco_run_dir_base="${busco_dir}/${busco_run_name}"
    busco_run_dir="$(find "${busco_run_dir_base}" -maxdepth 1 -type d -name 'run_*' -print -quit 2>/dev/null || true)"
    [[ -n "${busco_run_dir}" ]] || busco_run_dir="${busco_run_dir_base}"

    # BUSCO plots (if summaries exist)
    if compgen -G "${busco_run_dir}/short_summary*.json" > /dev/null; then
      env_run "${ENV_BUSCO}" busco --plot "${busco_run_dir}" --plot_percentages \
        || echo "[WARN] BUSCO plot generation failed (non-fatal)."
    else
      echo "[WARN] No short_summary JSON files in ${busco_run_dir}; skipping BUSCO plot."
    fi
    touch "${ck_busco}"
  else
    echo "[INFO] Step 7b: BUSCO checkpoint present — skipping"
  fi

  # BUSCO full_table.tsv for BlobToolKit
  busco_full_table="$(
    find "${busco_run_dir}" -maxdepth 2 -type f -name 'full_table.tsv' -print -quit 2>/dev/null || true
  )"
  if [[ -z "${busco_full_table}" ]]; then
    busco_full_table="$(
      find "${busco_dir}" -maxdepth 3 -type f -name 'full_table.tsv' -print -quit 2>/dev/null || true
    )"
  fi
  if [[ -n "${busco_full_table}" && -s "${busco_full_table}" ]]; then
    busco_full_table="$(readlink -f "${busco_full_table}")"
  else
    echo "[WARN] BUSCO full_table.tsv not found (expected under ${busco_run_dir}); BlobToolKit will use an empty table."
    : > "${btk_dir}/empty_busco.tsv"
    busco_full_table="${btk_dir}/empty_busco.tsv"
  fi

  # --------------------------- 7c) BlobToolKit ---------------------------
  echo "[INFO] Step 7c: BlobToolKit"
  if [[ ! -f "${ck_btk}" ]]; then
    mkdir -p "${btk_dir}" "${BLOB_TAXDUMP_DIR}"

    # Use the same assembly as BUSCO for BTK dataset
    bt_ref="${busco_target}"

    # Helper: ensure dict & fai for coverage mapping
    if [[ ! -s "${bt_ref}.dict" ]]; then
      echo "[INFO] Creating sequence dictionary for ${bt_ref}"
      singularity exec -B "${working_dir}:${working_dir}" "${sif_gatk4}" \
        gatk CreateSequenceDictionary -R "${bt_ref}" -O "${bt_ref}.dict" \
        || error "CreateSequenceDictionary failed for ${bt_ref}"
    fi
    if [[ ! -s "${bt_ref}.fai" ]]; then
      echo "[INFO] Indexing ${bt_ref} with samtools faidx"
      "${samtools_bin}" faidx "${bt_ref}" || error "samtools faidx failed for ${bt_ref}"
    fi

    # 1) Map HiFi reads for coverage
    echo "[INFO] Mapping PacBio HiFi reads to ${bt_ref} for coverage"
    "${minimap2_bin}" -ax map-hifi -H -Y --secondary=no -K 16G -t "$((nthreads * njobs))" -R "@RG\\tID:${base}.HiFi\\tSM:${base}\\tPL:PACBIO\\tLB:SMRTbell" "${bt_ref}" "${pacbio_hifi}" \
    | "${samtools_bin}" sort -@ "$((nthreads * njobs))" -m 8G -l 9 -O BAM -o "${bt_prefix}.bam" --write-index \
    || error "BlobToolKit mapping/sort failed"

	    echo "[INFO] Marking duplicates (BlobToolKit coverage)"
	    spark_threads="$((nthreads * njobs))"
		    markdup_shard_tmp="${tmpdir}/gatk_markdup.blobtoolkit.shards"
		    markdup_tmp="${tmpdir}/gatk_markdup.blobtoolkit.tmp"
		    rm -rf "${markdup_shard_tmp}" "${markdup_tmp}"
		    mkdir -p "${markdup_tmp}"
		    rm -f "${bt_prefix}.markdup.bam" "${bt_prefix}.markdup.bam.bai" "${bt_prefix}.markdup.bam.csi"
    gatk_spark \
      --java-options "-Xms128G -Xmx256G -Djava.io.tmpdir=${tmpdir}" \
      MarkDuplicatesSpark \
        --spark-master "local[${spark_threads}]" \
        --conf "spark.driver.host=127.0.0.1" \
        --conf "spark.driver.bindAddress=127.0.0.1" \
        --conf "spark.local.dir=${tmpdir}" \
        --input "${bt_prefix}.bam" \
        --metrics-file "${bt_prefix}.metrics" \
        --output "${bt_prefix}.markdup.bam" \
        --reference "${bt_ref}" \
	        --remove-all-duplicates true \
	        --create-output-bam-index true \
        --create-output-bam-splitting-index false \
        --read-name-regex null \
        --read-validation-stringency SILENT \
	        --tmp-dir "${markdup_tmp}" \
	        --output-shard-tmp-dir "${markdup_shard_tmp}" \
    2>&1 | tee -a "${btk_dir}/hifi.markduplicates.log" \
    || error "MarkDuplicates failed (BlobToolKit coverage)"

	    "${samtools_bin}" index -b -@ "${spark_threads}" "${bt_prefix}.markdup.bam" \
	      || error "Indexing markdup BAM failed"
	    rm -rf "${markdup_shard_tmp}" "${markdup_tmp}"

    # 2) BLAST hits
    hits_path="${btk_dir}/blast.hits"

    # Ensure we're in bash and be robust to empty/spacey vars
    if [[ -n ${blastn_bin-} ]] && [[ -x ${blastn_bin} ]]; then
      echo "[INFO] Running BLAST (megablast) for BlobToolKit hits"
      "${blastn_bin}" \
        -task megablast \
        -query "${bt_ref}" \
        -db "${BLAST_DB_ROOT}" \
        -outfmt '6 qseqid staxids bitscore std' \
        -num_threads "$((nthreads * njobs))" \
        -max_hsps 1 \
        -evalue 1e-25 \
        -out "${hits_path}" \
      || echo "[WARN] BLAST failed; continuing without hits"
    else
      echo "[WARN] blastn not found or not executable; proceeding without hits"
      : > "${hits_path}"
    fi

    # 3) Create BlobToolKit dataset (coverage + optional hits + optional BUSCO)
    info "Creating BlobToolKit dataset"
    mkdir -p "${blobdir_path}" "${plots_out}"

    # Sanity checks
    [[ -s "${bt_ref}" ]] || error "Missing FASTA for BlobToolKit: ${bt_ref}"
    [[ -s "${bt_prefix}.markdup.bam" ]] || error "Missing BAM for coverage: ${bt_prefix}.markdup.bam"
    [[ -s "${bt_prefix}.markdup.bam.csi" || -s "${bt_prefix}.markdup.bam.bai" ]] \
      || error "Missing BAM index for coverage: ${bt_prefix}.markdup.bam"

    # Optional BlobToolKit arguments
    btk_extra_args=()

    if [[ -s "${hits_path}" ]]; then
      btk_extra_args+=(--hits "${hits_path}")
    else
      warn "No BLAST hits file; taxonomy will rely on GC/coverage and taxid only"
    fi

    if [[ -s "${busco_full_table}" ]]; then
      btk_extra_args+=(--busco "${busco_full_table}")
    else
      warn "BUSCO full_table.tsv missing; BUSCO data will not be added"
    fi

    if [[ -s "${BLOB_TAXDUMP_DIR}/names.dmp" && -s "${BLOB_TAXDUMP_DIR}/nodes.dmp" ]]; then
      btk_extra_args+=(--taxdump "${BLOB_TAXDUMP_DIR}")
    else
      warn "Taxdump appears incomplete at ${BLOB_TAXDUMP_DIR}; taxonomy assignment may be incomplete"
    fi

    btk_threads="${nthreads:-1}"
    [[ "${btk_threads}" -ge 1 ]] || btk_threads=1

    if [[ -f "${ck7_btk_dataset}" && -s "${blobdir_path}/meta.json" ]]; then
      info "Step 7c.3: BlobToolKit dataset already complete"
    else
      rm -f "${ck7_btk_dataset}" "${ck_btk}"

      singularity exec \
        -B "${working_dir}:${working_dir}" \
        "${sif_blobtoolkit}" \
        blobtools add \
          --create \
          --fasta "${bt_ref}" \
          --cov "${bt_prefix}.markdup.bam" \
          "${btk_extra_args[@]}" \
          --taxrule bestsum \
          --taxid "${TAXID}" \
          --threads "${btk_threads}" \
          "${blobdir_path}" \
        || error "blobtools add --create failed"

      [[ -s "${blobdir_path}/meta.json" ]] \
        || error "BlobToolKit dataset creation failed: missing ${blobdir_path}/meta.json"

      touch "${ck7_btk_dataset}"
    fi

    if [[ -f "${ck7_btk_plots}" && -n "$(find "${plots_out}" -maxdepth 1 -type f -name '*.svg' -print -quit 2>/dev/null || true)" ]]; then
      info "Step 7c.4: BlobToolKit plots already complete"
    else
      info "Rendering BlobToolKit plots"
      rm -f "${ck7_btk_plots}"

      [[ -s "${blobdir_path}/meta.json" ]] \
        || error "Cannot render BlobToolKit plots: missing ${blobdir_path}/meta.json"

      run_blobtools_view() {
        local view="${1:?missing view}"
        local label="${2:-plot}"
        shift 2 || true

        local -a cmd=(
          singularity exec
          -B "${working_dir}:${working_dir}"
          "${sif_blobtoolkit}"
          xvfb-run -a
          blobtools view
          --plot
          --driver chromium
          --view "${view}"
          --format svg
          --out "${plots_out}"
        )

        (( $# > 0 )) && cmd+=("$@")
        cmd+=("${blobdir_path}")

        "${cmd[@]}" || warn "${label} failed"
      }

      run_blobtools_view snail "Snail plot"
      run_blobtools_view cumulative "Cumulative plot"
      run_blobtools_view blob "Blob plot" --param plotShape=circle

      if [[ -n "$(find "${plots_out}" -maxdepth 1 -type f -name '*.svg' -print -quit 2>/dev/null || true)" ]]; then
        touch "${ck7_btk_plots}"
      else
        warn "No BlobToolKit SVG plots were produced"
      fi
    fi

    [[ -s "${blobdir_path}/meta.json" ]] \
      || error "BlobToolKit dataset incomplete: missing ${blobdir_path}/meta.json"

	    touch "${ck_btk}"
	  else
	    echo "[INFO] Step 7c: BlobToolKit checkpoint present — skipping"
	  fi

  # --------------------------- 7d) Qualimap bamqc ---------------------------
  echo "[INFO] Step 7d: Qualimap bamqc (HiFi->assembly BAM)"
  if [[ ! -f "${ck_qualimap}" ]]; then
    activate_env "${ENV_QUALIMAP}"

  qualimap_logs="${qualimap_dir}/logs"
  mkdir -p "${qualimap_logs}"

  [[ -s "${bt_prefix}.markdup.bam" ]] || error "Qualimap input BAM missing: ${bt_prefix}.markdup.bam"
  [[ -s "${bt_prefix}.markdup.bam.csi" ]] || error "Qualimap BAM CSI index missing: ${bt_prefix}.markdup.bam.csi"

  qualimap_outdir="${qualimap_dir}/qualimap_${base}.HiFi"
  mkdir -p "${qualimap_outdir}"

  qualimap bamqc \
    -bam "${bt_prefix}.markdup.bam" \
    --paint-chromosome-limits \
    --skip-duplicated \
    --skip-dup-mode 2 \
    -nt "$((nthreads * njobs))" \
    -nr 1000000 \
    -outdir "${qualimap_outdir}" \
    -outformat PDF:HTML \
    --java-mem-size=256G \
  2>&1 | tee -a "${qualimap_logs}/qualimap.${LOG_TS}.log" \
  || error "Qualimap bamqc failed"
    deactivate_env
    touch "${ck_qualimap}"
  else
    echo "[INFO] Step 7d: Qualimap checkpoint present — skipping"
  fi

  # --------------------------- 7e) MultiQC (last) ---------------------------
  echo "[INFO] Step 7e: MultiQC"
  if [[ ! -f "${ck_multiqc}" ]]; then
    multiqc "${dir7}" -o "${dir7}" -p -ip -v || error "MultiQC failed"
    touch "${ck_multiqc}"
  else
    echo "[INFO] Step 7e: MultiQC checkpoint present — skipping"
  fi
  echo "[INFO] Step 7 complete. Reports at: ${dir7}"
fi

# ============================================================================================================================
# STEP 8: Repeats — RepeatModeler2, TEtrimmer, RepeatMasker, and Earl Grey annotation-only
# ============================================================================================================================
if is_true "${run_repeats:-FALSE}"; then
  echo "[INFO] Step 8: Repeat annotation — starting"

  repeats_target="${species_fasta:-${dir6}/${species_name}.fasta}"
  if [[ ! -s "${repeats_target}" ]]; then
    for cand in \
      "${dir6}/final.polished.decontam.fa" \
      "${dir6}/final.polished.fa" \
      "${dir5}/${base}.scaffolds.fa"; do
      if [[ -s "${cand}" ]]; then
        repeats_target="${cand}"
        break
      fi
    done
  fi
  [[ -s "${repeats_target}" ]] || error "Repeats input assembly not found (expected ${species_fasta:-${dir6}/${species_name}.fasta})"

  PROJECT_ROOT="${REPEATS_PROJECT_ROOT:-${repeats_dir}}"
  GENOME="${REPEATS_GENOME:-${repeats_target}}"
  SPECIES="${REPEATS_SPECIES:-${species_name}}"
  THREADS="${REPEATS_THREADS}"

  DFAM_FAMDB_DIR="${DFAM_FAMDB_DIR:-${PROJECT_ROOT}/Libraries/famdb}"
  REPBASE_EMBL="${REPBASE_EMBL:-${PROJECT_ROOT}/Libraries/RepBase-20240529.embl}"
  if [[ ! -s "$REPBASE_EMBL" && -s "$HOME/Documents/RepBase2024/RepBase-20240529.embl" ]]; then
    echo "[WARN] Project RepBase EMBL not found; using $HOME/Documents/RepBase2024/RepBase-20240529.embl"
    REPBASE_EMBL="$HOME/Documents/RepBase2024/RepBase-20240529.embl"
  fi

  if [[ "$ENV_RUNNER" == "mamba" ]]; then
    if [[ -z "${MAMBA_SH:-}" || ! -s "$MAMBA_SH" ]]; then
      echo "[ERROR] Mamba initialization script not found." >&2
      echo "[ERROR] Set MAMBA_SH globally or set ENV_RUNNER=conda." >&2
      exit 1
    fi
  fi

# ============================================================
# HELPER FUNCTIONS
# ============================================================

check_toggle() {
  local name="$1"
  local value="$2"

  case "${value^^}" in
    TRUE|FALSE) ;;
    *)
      echo "[ERROR] Invalid toggle for $name=$value. Use TRUE or FALSE." >&2
      exit 1
      ;;
  esac
}

validate_yes_no() {
  local name="$1"
  local value="$2"

  case "${value,,}" in
    yes|no) ;;
    *)
      echo "[ERROR] Invalid value for $name=$value. Use yes or no." >&2
      exit 1
      ;;
  esac
}

require_file() {
  local path="$1"
  local label="$2"

  [[ -s "$path" ]] || {
    echo "[ERROR] Required file missing for $label: $path" >&2
    exit 1
  }
}

require_dir() {
  local path="$1"
  local label="$2"

  [[ -d "$path" ]] || {
    echo "[ERROR] Required directory missing for $label: $path" >&2
    exit 1
  }
}

find_newest_file_by_name() {
  local dir="$1"
  shift

  [[ -d "$dir" ]] || return 0

  find "$dir" -type f "$@" -printf '%T@\t%p\n' 2>/dev/null \
    | sort -nr \
    | head -n 1 \
    | cut -f2- || true
}

file_has_fasta_records() {
  local path="$1"

  [[ -s "$path" ]] || return 1
  grep -q "^>" "$path"
}

dir_has_nonempty_file_or_link() {
  local dir="$1"
  local f

  [[ -d "$dir" ]] || return 1
  while IFS= read -r -d '' f; do
    [[ -s "$f" ]] && return 0
  done < <(find "$dir" \( -type f -o -type l \) -print0 2>/dev/null)

  return 1
}

auto_detect_outputs() {
  if ! is_true "$AUTO_RESUME"; then
    return 0
  fi

  echo "[INFO] Auto-resume enabled: scanning for existing outputs"

  if [[ -n "${KNOWN_LIBRARY_OVERRIDE:-}" ]]; then
    KNOWN_NR="$KNOWN_LIBRARY_OVERRIDE"
  elif [[ ! -s "${KNOWN_NR:-}" ]]; then
    detected_known="$(find_newest_file_by_name "$DIR_KNOWN" \( -name "known_repeats.nr.fa" -o -name "*known*nr*.fa" \))"
    [[ -n "$detected_known" ]] && KNOWN_NR="$detected_known"
  fi

  if [[ -n "${REPEATMODELER_LIBRARY_OVERRIDE:-}" ]]; then
    RM_LIBRARY="$REPEATMODELER_LIBRARY_OVERRIDE"
  elif [[ ! -s "${RM_LIBRARY:-}" ]]; then
    detected_rm="$(find_newest_file_by_name "$DIR_RM" \( -name "${GENOME_STEM}-families.fa" -o -name "*families.fa" \))"
    [[ -n "$detected_rm" ]] && RM_LIBRARY="$detected_rm"
  fi

  if [[ -n "${TETRIMMER_LIBRARY_OVERRIDE:-}" ]]; then
    TETRIMMER_LIBRARY="$TETRIMMER_LIBRARY_OVERRIDE"
  else
    if [[ -s "${DIR_TETRIMMER}/TEtrimmer_consensus_merged.fasta" ]]; then
      TETRIMMER_LIBRARY="${DIR_TETRIMMER}/TEtrimmer_consensus_merged.fasta"
    elif [[ -s "${DIR_TETRIMMER}/TEtrimmer_consensus.fasta" ]]; then
      echo "[WARN] Only TEtrimmer_consensus.fasta found; merged/deduplicated library missing."
      echo "[WARN] Using non-merged TEtrimmer library as fallback."
      TETRIMMER_LIBRARY="${DIR_TETRIMMER}/TEtrimmer_consensus.fasta"
    else
      detected_tetrimmer="$(find_newest_file_by_name "$DIR_TETRIMMER" \( -name "TEtrimmer_consensus_merged.fasta" -o -name "TEtrimmer_consensus.fasta" -o -name "*consensus*merged*.fa" -o -name "*consensus*merged*.fasta" \))"
      [[ -n "$detected_tetrimmer" ]] && TETRIMMER_LIBRARY="$detected_tetrimmer"
    fi
  fi

  if [[ -n "${FINAL_LIBRARY_OVERRIDE:-}" ]]; then
    FINAL_NR="$FINAL_LIBRARY_OVERRIDE"
  elif [[ ! -s "${FINAL_NR:-}" ]]; then
    detected_final="$(find_newest_file_by_name "$DIR_LIB" \( -name "${SPECIES}.TEtrimmer_Dfam_RepBase.nr.fa" -o -name "*TEtrimmer*Dfam*RepBase*nr*.fa" -o -name "*final*nr*.fa" \))"
    [[ -n "$detected_final" ]] && FINAL_NR="$detected_final"
  fi

  if [[ -s "${KNOWN_NR:-}" ]]; then
    echo "[INFO] Auto-detected known library: $KNOWN_NR"
  fi

  if [[ -s "${RM_LIBRARY:-}" ]]; then
    echo "[INFO] Auto-detected RepeatModeler library: $RM_LIBRARY"
  fi

  if [[ -s "${TETRIMMER_LIBRARY:-}" ]]; then
    echo "[INFO] Auto-detected TEtrimmer library: $TETRIMMER_LIBRARY"
  fi

  if [[ -s "${FINAL_NR:-}" ]]; then
    echo "[INFO] Auto-detected final TE library: $FINAL_NR"
  fi

  return 0
}

step_output_ready() {
  local step="$1"

  case "$step" in
    known)
      [[ -s "${KNOWN_NR:-}" ]] && file_has_fasta_records "$KNOWN_NR"
      ;;
    repeatmodeler)
      [[ -s "${RM_LIBRARY:-}" ]] && file_has_fasta_records "$RM_LIBRARY"
      ;;
    tetrimmer)
      tetrimmer_core_outputs_complete
      ;;
    final_library)
      [[ -s "${FINAL_NR:-}" ]] && file_has_fasta_records "$FINAL_NR"
      ;;
    repeatmasker)
      find "$DIR_REPEATMASKER" -maxdepth 1 -type f \( -name "*.out" -o -name "*.gff" -o -name "*.tbl" -o -name "*.masked" \) -size +0c 2>/dev/null | grep -q .
      ;;
    earlgrey)
      find "$DIR_EARLGREY" -type f \( -name "*.gff" -o -name "*.gff3" -o -name "*.summary" -o -name "*.tbl" -o -name "*.masked" -o -name "*.fa" -o -name "*.fasta" \) -size +0c 2>/dev/null | grep -q .
      ;;
    earlgrey_plots)
      if ! is_true "${EARLGREY_ENHANCED_PLOTS:-TRUE}"; then
        return 0
      fi
      [[ -d "${DIR_EARLGREY_ENHANCED_PLOTS:-}" ]] || return 1
      [[ -s "${DIR_EARLGREY_ENHANCED_PLOTS}/earlgrey_enhanced_plots_manifest.tsv" ]] || return 1
      find "$DIR_EARLGREY_ENHANCED_PLOTS" -maxdepth 1 -type f \( -name "*.png" -o -name "*.pdf" -o -name "*.svg" -o -name "*.vvg" \) -size +0c 2>/dev/null | grep -q .
      ;;
    collect_outputs)
      [[ -d "$DIR_FINAL/logs" ]] \
        && [[ -s "${DIR_FINAL}/logs/collect_outputs_manifest.tsv" ]] \
        && dir_has_nonempty_file_or_link "$DIR_FINAL"
      ;;
    non_repeats)
      [[ -s "${NON_REPEATS_BED:-}" && -s "${NON_REPEATS_GENOME_SIZES:-}" && -s "${NON_REPEATS_MERGED_BED:-}" ]]
      ;;
    *)
      return 1
      ;;
  esac
}

archive_unrecoverable_repeatmodeler_dirs() {
  [[ -d "$DIR_RM" ]] || return 0

  newest_rm_dir="$(find "$DIR_RM" -maxdepth 1 -type d -name 'RM_*' -printf '%T@\t%p\n' 2>/dev/null | sort -nr | head -n1 | cut -f2- || true)"
  [[ -n "$newest_rm_dir" && -d "$newest_rm_dir" ]] || return 0

  if grep -Rqs "did not get passed round-1\|Remove the -recoverDir option" "$newest_rm_dir" "${DIR_LOGS}/01_repeatmodeler.log" 2>/dev/null; then
    if is_true "$AUTO_ARCHIVE_UNRECOVERABLE_RM_DIRS"; then
      mkdir -p "$DIR_RM/failed_runs"
      base="$(basename "$newest_rm_dir")"
      stamp="$(date +%Y%m%d_%H%M%S)"
      dest="$DIR_RM/failed_runs/${base}.${stamp}"
      echo "[WARN] Archiving unrecoverable RepeatModeler directory: $newest_rm_dir -> $dest"
      mv "$newest_rm_dir" "$dest"
    else
      echo "[WARN] Unrecoverable RepeatModeler directory detected but AUTO_ARCHIVE_UNRECOVERABLE_RM_DIRS=FALSE: $newest_rm_dir"
    fi
  fi
}

tetrimmer_core_outputs_complete() {
  [[ -s "${DIR_TETRIMMER}/TEtrimmer_consensus_merged.fasta" || -s "${DIR_TETRIMMER}/TEtrimmer_consensus.fasta" ]] || return 1
  [[ -s "${DIR_TETRIMMER}/summary.txt" ]] || return 1
  [[ -s "${DIR_TETRIMMER}/Sequence_name_mapping.txt" ]] || return 1
  return 0
}

tetrimmer_dir_resumable() {
  # Conservative policy:
  # TEtrimmer --continue_analysis is strict and may abort if the directory
  # was not created by a compatible interrupted run. Therefore, automatic
  # mode does not infer resumability from logs or partial folders.
  # Use TETRIMMER_CONTINUE_MODE=always only when manually confirmed.
  return 1
}

prepare_tetrimmer_output_dir() {
  [[ -d "$DIR_TETRIMMER" ]] || return 0

  if tetrimmer_core_outputs_complete; then
    return 0
  fi

  # If no real TEtrimmer outputs exist but the directory contains only
  # logs/partial initialization from a failed run, archive it and start fresh.
  if is_true "$TETRIMMER_ARCHIVE_INCOMPLETE"; then
    mkdir -p "${PROJECT_ROOT}/failed_runs/tetrimmer"
    stamp="$(date +%Y%m%d_%H%M%S)"
    dest="${PROJECT_ROOT}/failed_runs/tetrimmer/$(basename "$DIR_TETRIMMER").${stamp}"
    echo "[WARN] TEtrimmer output directory is incomplete or stale."
    echo "[WARN] Archiving: $DIR_TETRIMMER -> $dest"
    mv "$DIR_TETRIMMER" "$dest"
    mkdir -p "$DIR_TETRIMMER"
  else
    echo "[ERROR] TEtrimmer output directory is incomplete: $DIR_TETRIMMER" >&2
    echo "[ERROR] Remove/archive it manually or set TETRIMMER_ARCHIVE_INCOMPLETE=TRUE." >&2
    exit 1
  fi

  return 0
}

audit_tetrimmer_outputs() {
  {
    echo "TEtrimmer output audit"
    echo "======================"
    echo "Output dir: $DIR_TETRIMMER"
    echo

    for f in \
      "$DIR_TETRIMMER/TEtrimmer_consensus.fasta" \
      "$DIR_TETRIMMER/TEtrimmer_consensus_merged.fasta" \
      "$DIR_TETRIMMER/summary.txt" \
      "$DIR_TETRIMMER/Sequence_name_mapping.txt" \
      "$DIR_TETRIMMER/error_file.txt"
    do
      if [[ -s "$f" ]]; then
        echo "[OK] $f"
      else
        echo "[MISSING] $f"
      fi
    done

    echo
    echo "Directories:"
    for d in \
      Classification \
      Single_fasta_files \
      TEtrimmer_for_proof_curation \
      Annotation_perfect \
      Annotation_good \
      Annotation_check_recommended \
      Annotation_check_required \
      Clustered_proof_curation \
      TE_low_copy \
      TE_skipped \
      HMM \
      Multiple_sequence_alignment
    do
      if [[ -d "$DIR_TETRIMMER/$d" ]]; then
        n=$(find "$DIR_TETRIMMER/$d" -type f 2>/dev/null | wc -l)
        echo "[OK] $d files=$n"
      else
        echo "[MISSING] $d"
      fi
    done
  } > "${DIR_LOGS}/02_tetrimmer.output_audit.txt"
}

known_library_needs_build() {
  if repeat_checkpoint_ready known; then
    return 1
  fi

  return 0
}

validate_toggles() {
  check_toggle AUTO_RESUME "$AUTO_RESUME"
  check_toggle AUTO_ARCHIVE_UNRECOVERABLE_RM_DIRS "$AUTO_ARCHIVE_UNRECOVERABLE_RM_DIRS"
  check_toggle TETRIMMER_ARCHIVE_INCOMPLETE "$TETRIMMER_ARCHIVE_INCOMPLETE"
  check_toggle REPEATMODELER_RECOVER "$REPEATMODELER_RECOVER"
  check_toggle REPEATMODELER_SANITIZE_FASTA_HEADERS "$REPEATMODELER_SANITIZE_FASTA_HEADERS"
}

validate_runner() {
  local name="$1"
  local value="$2"

  case "$value" in
    mamba|conda) ;;
    *)
      echo "[ERROR] Invalid runner for $name=$value. Use mamba or conda." >&2
      exit 1
      ;;
  esac
}

repeat_checkpoint_file() {
  printf '%s/%s.done' "$REPEATS_CHECKPOINT_DIR" "$1"
}

repeat_checkpoint_ready() {
  local step="$1"
  local ck
  ck="$(repeat_checkpoint_file "$step")"

  if ! is_true "$AUTO_RESUME"; then
    return 1
  fi

  if step_output_ready "$step"; then
    [[ -f "$ck" ]] || touch "$ck"
    return 0
  fi

  if [[ -f "$ck" ]]; then
    echo "[WARN] Removing stale repeat checkpoint without valid outputs: $ck" >&2
    rm -f "$ck"
  fi

  return 1
}

repeat_mark_checkpoint() {
  local step="$1"
  local ck
  ck="$(repeat_checkpoint_file "$step")"

  step_output_ready "$step" || {
    echo "[ERROR] Refusing to mark checkpoint for incomplete repeat step: $step" >&2
    exit 1
  }

  touch "$ck"
}

write_earlgrey_enhanced_plotter() {
  local script_path="$1"

  cat > "$script_path" <<'RSCRIPT'
#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(flag, default = NA_character_) {
  idx <- which(args == flag)
  if (length(idx) == 0 || idx[1] == length(args)) return(default)
  args[idx[1] + 1]
}

species <- arg_value("--species")
earlgrey_dir <- arg_value("--earlgrey_dir")
outdir <- arg_value("--outdir")
dpi <- suppressWarnings(as.integer(arg_value("--dpi", "600")))
top_n <- suppressWarnings(as.integer(arg_value("--top_families", "30")))

if (is.na(species) || !nzchar(species)) stop("--species is required")
if (is.na(earlgrey_dir) || !dir.exists(earlgrey_dir)) stop("--earlgrey_dir must exist")
if (is.na(outdir) || !nzchar(outdir)) stop("--outdir is required")
if (!is.finite(dpi) || dpi < 1) dpi <- 600L
if (!is.finite(top_n) || top_n < 1) top_n <- 30L

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

for (pkg in c("ggplot2")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Missing required R package: ", pkg)
  }
}

read_tsv <- function(path) {
  utils::read.delim(path, sep = "\t", header = TRUE, check.names = FALSE,
                    stringsAsFactors = FALSE, quote = "", comment.char = "")
}

num_col <- function(x) {
  suppressWarnings(as.numeric(gsub(",", "", as.character(x))))
}

find_first <- function(pattern, prefer = character()) {
  files <- list.files(earlgrey_dir, pattern = pattern, recursive = TRUE, full.names = TRUE)
  files <- files[file.exists(files) & file.info(files)$size > 0]
  if (length(files) == 0) return(NA_character_)
  if (length(prefer) > 0) {
    score <- rep(length(prefer) + 1L, length(files))
    for (i in seq_along(prefer)) {
      score[grepl(prefer[[i]], files, fixed = TRUE)] <- i
    }
    files <- files[order(score, files)]
  } else {
    files <- sort(files)
  }
  files[[1]]
}

class_base <- function(x) {
  x <- as.character(x)
  x <- sub("-nested$", "", x)
  x <- sub("/.*$", "", x)
  x <- sub("-.*$", "", x)
  x <- ifelse(grepl("^RC$|Rolling Circle", x), "Rolling Circle", x)
  x <- ifelse(grepl("Penelope|^PLE$", x), "Penelope", x)
  x <- ifelse(grepl("^Unknown|Unclassified|unknown", x), "Unclassified", x)
  x <- ifelse(grepl("Simple_repeat|Low_complexity|Satellite|RNA|Microsatellite|Other", x),
              "Other simple or low-complexity", x)
  x
}

extract_attr <- function(x, key) {
  pattern <- paste0("(^|;)", key, "=([^;]+)")
  hit <- regexec(pattern, x, ignore.case = TRUE)
  out <- regmatches(x, hit)
  vapply(out, function(z) if (length(z) >= 3) z[[3]] else NA_character_, character(1))
}

class_palette <- c(
  "DNA" = "#b2182b",
  "Rolling Circle" = "#ef8a62",
  "Penelope" = "#7b3294",
  "LINE" = "#2166ac",
  "SINE" = "#c51b7d",
  "LTR" = "#1b7837",
  "Unclassified" = "#636363",
  "Other simple or low-complexity" = "#a6dba0",
  "Non-Repeat" = "#000000",
  "Other" = "#bdbdbd"
)

theme_eg <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey92", colour = "grey70")
    )
}

save_plot <- function(plot, stem, width = 9, height = 6) {
  outputs <- file.path(outdir, paste0(stem, c(".png", ".pdf", ".svg", ".vvg")))
  ggplot2::ggsave(outputs[[1]], plot = plot, width = width, height = height,
                  units = "in", dpi = dpi)
  ggplot2::ggsave(outputs[[2]], plot = plot, width = width, height = height,
                  units = "in", device = grDevices::cairo_pdf)
  ggplot2::ggsave(outputs[[3]], plot = plot, width = width, height = height,
                  units = "in", device = "svg", dpi = dpi)
  ggplot2::ggsave(outputs[[4]], plot = plot, width = width, height = height,
                  units = "in", device = "svg", dpi = dpi)
  data.frame(plot = stem, file = outputs, stringsAsFactors = FALSE)
}

manifest <- data.frame(plot = character(), file = character(), stringsAsFactors = FALSE)

high_file <- find_first(paste0("(^|/)", species, "\\.highLevelCount\\.txt$"),
                        prefer = c("_summaryFiles"))
family_file <- find_first(paste0("(^|/)", species, "\\.familyLevelCount\\.txt$"),
                          prefer = c("_summaryFiles"))
div_file <- find_first(paste0("(^|/)", species, "(_divergence_summary_table|_summary_table)\\.tsv$"),
                       prefer = c("_summaryFiles", "_RepeatLandscape"))
gff_file <- find_first(paste0("(^|/)", species, "\\.filteredRepeats\\.withDivergence\\.gff$"),
                       prefer = c("_RepeatLandscape", "_summaryFiles", "looseMerge"))
if (is.na(gff_file)) {
  gff_file <- find_first(paste0("(^|/)", species, "\\.filteredRepeats\\.gff$"),
                         prefer = c("_summaryFiles", "looseMerge"))
}

if (!is.na(high_file)) {
  high <- read_tsv(high_file)
  required <- c("TE Classification", "Coverage (bp)", "Copy Number", "% Genome Coverage")
  if (all(required %in% names(high))) {
    high$coverage_bp <- num_col(high[["Coverage (bp)"]])
    high$genome_pct <- num_col(high[["% Genome Coverage"]])
    high$copy_number_num <- num_col(high[["Copy Number"]])
    high$is_nested <- grepl("-nested$", high[["TE Classification"]])
    high$class_base <- class_base(high[["TE Classification"]])
    high_plot <- high[!high[["TE Classification"]] %in% c("Total Interspersed Repeat", "Non-Repeat"), , drop = FALSE]

    p_comp <- ggplot2::ggplot(high_plot,
                              ggplot2::aes(x = stats::reorder(class_base, genome_pct),
                                           y = genome_pct, fill = class_base)) +
      ggplot2::geom_col(width = 0.75) +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values = class_palette, na.value = "#bdbdbd") +
      ggplot2::labs(title = paste(species, "EarlGrey repeat composition"),
                    x = NULL, y = "Genome coverage (%)") +
      theme_eg()
    manifest <- rbind(manifest, save_plot(p_comp, "earlgrey_repeat_composition", 9, 5.5))

    nested <- stats::aggregate(coverage_bp ~ class_base + is_nested, data = high_plot, sum, na.rm = TRUE)
    nested$nested_status <- ifelse(nested$is_nested, "Nested", "Non-nested")
    p_nested <- ggplot2::ggplot(nested,
                                ggplot2::aes(x = stats::reorder(class_base, coverage_bp),
                                             y = coverage_bp / 1e6, fill = nested_status)) +
      ggplot2::geom_col(width = 0.75) +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values = c("Non-nested" = "#4d4d4d", "Nested" = "#d95f02")) +
      ggplot2::labs(title = paste(species, "nested repeat burden"),
                    x = NULL, y = "Coverage (Mb)") +
      theme_eg()
    manifest <- rbind(manifest, save_plot(p_nested, "earlgrey_nested_repeat_burden", 9, 5.5))
  }
}

if (!is.na(family_file)) {
  fam <- read_tsv(family_file)
  if (all(c("TE Family", "Coverage (bp)", "Copy Number") %in% names(fam))) {
    fam$coverage_bp <- num_col(fam[["Coverage (bp)"]])
    fam$copy_number_num <- num_col(fam[["Copy Number"]])
    fam$family <- sub("#.*$", "", fam[["TE Family"]])
    fam$te_type <- ifelse(grepl("#", fam[["TE Family"]]), sub("^.*#", "", fam[["TE Family"]]), "Unknown")
    fam$class_base <- class_base(fam$te_type)
    fam <- fam[is.finite(fam$coverage_bp) & fam$coverage_bp > 0, , drop = FALSE]
    fam <- fam[order(-fam$coverage_bp), , drop = FALSE]
    top <- utils::head(fam, top_n)

    p_top <- ggplot2::ggplot(top,
                             ggplot2::aes(x = stats::reorder(family, coverage_bp),
                                          y = coverage_bp / 1e6, fill = class_base)) +
      ggplot2::geom_col(width = 0.75) +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values = class_palette, na.value = "#bdbdbd") +
      ggplot2::labs(title = paste(species, "top TE families by coverage"),
                    x = NULL, y = "Coverage (Mb)") +
      theme_eg(base_size = 10)
    manifest <- rbind(manifest, save_plot(p_top, "earlgrey_top_te_families", 10, 8))

    scatter <- fam[is.finite(fam$copy_number_num) & fam$copy_number_num > 0, , drop = FALSE]
    if (nrow(scatter) > 0) {
      p_scatter <- ggplot2::ggplot(scatter,
                                   ggplot2::aes(x = copy_number_num,
                                                y = coverage_bp / 1e6,
                                                colour = class_base)) +
        ggplot2::geom_point(alpha = 0.65, size = 1.8) +
        ggplot2::scale_x_log10() +
        ggplot2::scale_y_log10() +
        ggplot2::scale_colour_manual(values = class_palette, na.value = "#bdbdbd") +
        ggplot2::labs(title = paste(species, "TE family copy number versus coverage"),
                      x = "Copy number (log10)", y = "Coverage, Mb (log10)") +
        theme_eg()
      manifest <- rbind(manifest, save_plot(p_scatter, "earlgrey_family_copy_vs_coverage", 8, 6))
    }
  }
}

if (!is.na(div_file)) {
  div <- read_tsv(div_file)
  if (all(c("subclass", "superfamily", "family", "total_bp", "mean_div") %in% names(div))) {
    div$total_bp <- num_col(div$total_bp)
    div$mean_div <- num_col(div$mean_div)
    div <- div[is.finite(div$total_bp) & is.finite(div$mean_div) & div$total_bp > 0, , drop = FALSE]
    div$class_base <- class_base(div$subclass)
    if (nrow(div) > 0) {
      p_div <- ggplot2::ggplot(div,
                               ggplot2::aes(x = mean_div, y = total_bp / 1e6,
                                            colour = class_base)) +
        ggplot2::geom_point(alpha = 0.7, size = 1.8) +
        ggplot2::facet_wrap(stats::as.formula("~ class_base"), scales = "free_y") +
        ggplot2::scale_colour_manual(values = class_palette, na.value = "#bdbdbd") +
        ggplot2::labs(title = paste(species, "family-level Kimura divergence"),
                      x = "Mean Kimura 2-parameter distance", y = "Family coverage (Mb)") +
        theme_eg(base_size = 10)
      manifest <- rbind(manifest, save_plot(p_div, "earlgrey_family_divergence", 10, 7))
    }
  }
}

if (!is.na(gff_file)) {
  gff <- utils::read.delim(gff_file, sep = "\t", header = FALSE, quote = "",
                           comment.char = "", fill = TRUE, stringsAsFactors = FALSE,
                           col.names = paste0("V", seq_len(9)))
  gff <- gff[ncol(gff) >= 9 & grepl("kimura80=", gff$V9, ignore.case = TRUE), , drop = FALSE]
  if (nrow(gff) > 0) {
    gff$class_base <- class_base(gff$V3)
    gff$kimura80 <- suppressWarnings(as.numeric(extract_attr(gff$V9, "kimura80")))
    gff$width_bp <- num_col(gff$V5) - num_col(gff$V4) + 1
    gff <- gff[is.finite(gff$kimura80) & gff$kimura80 >= 0 & gff$kimura80 <= 0.6 &
                 is.finite(gff$width_bp) & gff$width_bp > 0, , drop = FALSE]
    if (nrow(gff) > 0) {
      p_kimura <- ggplot2::ggplot(gff,
                                  ggplot2::aes(x = kimura80, weight = width_bp,
                                               fill = class_base)) +
        ggplot2::geom_histogram(binwidth = 0.01, position = "identity", alpha = 0.65,
                                boundary = 0) +
        ggplot2::facet_wrap(stats::as.formula("~ class_base"), scales = "free_y") +
        ggplot2::scale_fill_manual(values = class_palette, na.value = "#bdbdbd") +
        ggplot2::labs(title = paste(species, "Kimura repeat landscape"),
                      x = "Kimura 2-parameter distance", y = "Repeat bp") +
        theme_eg(base_size = 10)
      manifest <- rbind(manifest, save_plot(p_kimura, "earlgrey_kimura_landscape", 10, 7))
    }
  }
}

if (nrow(manifest) == 0) {
  stop("No enhanced EarlGrey plots were created; required input files were missing or empty")
}

utils::write.table(manifest, file.path(outdir, "earlgrey_enhanced_plots_manifest.tsv"),
                   sep = "\t", row.names = FALSE, quote = FALSE)

input_manifest <- data.frame(
  key = c("species", "earlgrey_dir", "high_level_count", "family_level_count",
          "divergence_summary", "gff_with_kimura", "dpi", "top_families"),
  value = c(species, earlgrey_dir, high_file, family_file, div_file, gff_file,
            as.character(dpi), as.character(top_n)),
  stringsAsFactors = FALSE
)
utils::write.table(input_manifest, file.path(outdir, "earlgrey_enhanced_plots_inputs.tsv"),
                   sep = "\t", row.names = FALSE, quote = FALSE)
RSCRIPT

  chmod +x "$script_path"
}

# ============================================================
# VALIDATION
# ============================================================

validate_toggles

validate_runner ENV_RUNNER "$ENV_RUNNER"

validate_yes_no RM_XSMALL "$RM_XSMALL"
validate_yes_no RM_GFF "$RM_GFF"
validate_yes_no RM_ALIGN "$RM_ALIGN"

validate_yes_no EARLGREY_REMOVE_SHORT "$EARLGREY_REMOVE_SHORT"
validate_yes_no EARLGREY_SOFTMASK "$EARLGREY_SOFTMASK"
validate_yes_no EARLGREY_HELIANO "$EARLGREY_HELIANO"
check_toggle EARLGREY_ENHANCED_PLOTS "$EARLGREY_ENHANCED_PLOTS"

if ! [[ "$EARLGREY_PLOT_DPI" =~ ^[0-9]+$ ]] || [[ "$EARLGREY_PLOT_DPI" -lt 1 ]]; then
  echo "[ERROR] Invalid EARLGREY_PLOT_DPI=$EARLGREY_PLOT_DPI. Use a positive integer." >&2
  exit 1
fi

if ! [[ "$EARLGREY_PLOT_TOP_FAMILIES" =~ ^[0-9]+$ ]] || [[ "$EARLGREY_PLOT_TOP_FAMILIES" -lt 1 ]]; then
  echo "[ERROR] Invalid EARLGREY_PLOT_TOP_FAMILIES=$EARLGREY_PLOT_TOP_FAMILIES. Use a positive integer." >&2
  exit 1
fi

case "$TETRIMMER_CONTINUE_MODE" in
  auto|always|never) ;;
  *)
    echo "[ERROR] Invalid TETRIMMER_CONTINUE_MODE=$TETRIMMER_CONTINUE_MODE. Use auto, always, or never." >&2
    exit 1
    ;;
esac

# ============================================================
# PATH NORMALIZATION AND DIRECTORIES
# ============================================================

PROJECT_ROOT="$(readlink -f "$PROJECT_ROOT")"
GENOME="$(readlink -f "$GENOME")"
sif_repeatmasker="$(readlink -f "$sif_repeatmasker")"

DIR_GENOME="${PROJECT_ROOT}/00.genome"
DIR_KNOWN="${PROJECT_ROOT}/01.known_libraries"
DIR_RM="${PROJECT_ROOT}/02.repeatmodeler"
DIR_TETRIMMER="${PROJECT_ROOT}/03.tetrimmer"
DIR_LIB="${PROJECT_ROOT}/04.final_library"
DIR_REPEATMASKER="${PROJECT_ROOT}/05.repeatmasker"
DIR_EARLGREY="${PROJECT_ROOT}/06.earlgrey_annotation_only"
DIR_FINAL="${PROJECT_ROOT}/07.final_outputs"
DIR_EARLGREY_ENHANCED_PLOTS="${DIR_FINAL}/earlgrey/enhanced_plots"
DIR_NON_REPEATS="${PROJECT_ROOT}/08.non_repetitive_regions"
DIR_LOGS="${PROJECT_ROOT}/logs"
REPEATS_CHECKPOINT_DIR="${PROJECT_ROOT}/checkpoints"
NON_REPEATS_GENOME_SIZES="${DIR_NON_REPEATS}/genome.sizes"
NON_REPEATS_RAW_BED="${DIR_NON_REPEATS}/repeats.raw.bed"
NON_REPEATS_SORTED_BED="${DIR_NON_REPEATS}/repeats.sorted.bed"
NON_REPEATS_MERGED_BED="${DIR_NON_REPEATS}/repeats.merged.bed"
NON_REPEATS_BED="${DIR_NON_REPEATS}/non_repetitive_regions.bed"

mkdir -p \
  "$DIR_GENOME" \
  "$DIR_KNOWN" \
  "$DIR_RM" \
  "$DIR_TETRIMMER" \
  "$DIR_LIB" \
  "$DIR_REPEATMASKER" \
  "$DIR_EARLGREY" \
  "$DIR_FINAL" \
  "$DIR_EARLGREY_ENHANCED_PLOTS" \
  "$DIR_NON_REPEATS" \
  "$DIR_LOGS" \
  "$REPEATS_CHECKPOINT_DIR"

GENOME_BASENAME="$(basename "$GENOME")"
GENOME_STEM="${GENOME_BASENAME%%.*}"
GENOME_LOCAL="${DIR_GENOME}/${GENOME_BASENAME}"

if [[ "$REPEATMODELER_SANITIZE_FASTA_HEADERS" == "TRUE" ]]; then
  echo "[INFO] Creating sanitized local genome FASTA"

  awk '
    /^>/ {
      h=$0
      sub(/^>/, "", h)
      split(h, a, /[ \t]/)
      h=a[1]
      gsub(/[^A-Za-z0-9_.:|+-]/, "_", h)
      print ">" h
      next
    }
    {
      print
    }
  ' "$GENOME" > "$GENOME_LOCAL"
else
  cp -f "$GENOME" "$GENOME_LOCAL"
fi

# ============================================================
# OUTPUT PATHS
# ============================================================

DFAM_FASTA="${DIR_KNOWN}/Dfam.taxon_${DFAM_TAXON}.fa"
REPBASE_FASTA="${DIR_KNOWN}/RepBase.fa"
KNOWN_RAW="${DIR_KNOWN}/known_repeats.raw.fa"
KNOWN_NR="${DIR_KNOWN}/known_repeats.nr.fa"

RM_LIBRARY="${DIR_RM}/${GENOME_STEM}-families.fa"

FINAL_RAW="${DIR_LIB}/${SPECIES}.TEtrimmer_Dfam_RepBase.raw.fa"
FINAL_NR="${DIR_LIB}/${SPECIES}.TEtrimmer_Dfam_RepBase.nr.fa"

# Detect outputs from previous runs before deciding which steps to execute.
auto_detect_outputs || true


# ============================================================
# PREFLIGHT
# ============================================================

echo "[INFO] Running repeat annotation preflight checks"

  command -v singularity >/dev/null 2>&1 || {
    echo "[ERROR] singularity not found in PATH." >&2
    exit 1
  }

  command -v "$ENV_RUNNER" >/dev/null 2>&1 || {
    echo "[ERROR] $ENV_RUNNER not found in PATH." >&2
    exit 1
  }

  require_file "$GENOME_LOCAL" "genome"
  require_file "$sif_repeatmasker" "RepeatMasker/RepeatModeler SIF"
  require "$bedtools_bin"

  if known_library_needs_build; then
    require_dir "$DFAM_FAMDB_DIR" "Dfam FamDB directory"
    require_file "$REPBASE_EMBL" "RepBase EMBL file"
  else
    echo "[INFO] Step 0 raw-library inputs are not required because known library already exists or Step 0 is disabled"
  fi

  if [[ -n "$ARTEFACTS_FASTA" ]]; then
    require_file "$ARTEFACTS_FASTA" "artefacts FASTA"
  fi

  echo "[INFO] Auditing RepeatMasker/RepeatModeler container"

  # shellcheck disable=SC2016
  singularity exec \
    --no-home \
    -B "$PROJECT_ROOT":"$PROJECT_ROOT" \
    "$sif_repeatmasker" \
    bash -lc '
      set -Eeuo pipefail

      missing=""

      for tool in perl RepeatMasker RepeatModeler BuildDatabase cd-hit-est gunzip; do
        command -v "$tool" >/dev/null 2>&1 || missing="$missing $tool"
      done

      command -v makeblastdb >/dev/null 2>&1 || \
        [[ -x /opt/rmblast/bin/makeblastdb ]] || \
        missing="$missing makeblastdb"

      [[ -x /opt/RepeatMasker/util/buildRMLibFromEMBL.pl ]] || \
        missing="$missing buildRMLibFromEMBL.pl"

      if [[ -n "$missing" ]]; then
        echo "[ERROR] Missing tools inside container:$missing" >&2
        exit 1
      fi
    '

  if known_library_needs_build; then
    echo "[INFO] Auditing Dfam export tool inside container"

    singularity exec \
      --no-home \
      -B "$DFAM_FAMDB_DIR":"$DFAM_FAMDB_DIR" \
      "$sif_repeatmasker" \
      bash -lc '
        set -Eeuo pipefail

        command -v famdb.py >/dev/null 2>&1 || {
          echo "[ERROR] famdb.py not found inside container." >&2
          exit 1
        }
      '
  fi

  if ! { repeat_checkpoint_ready tetrimmer; }; then
    echo "[INFO] Auditing TEtrimmer environment"

    env_run "$ENV_TETRIMMER" bash -lc '
      set -Eeuo pipefail

      command -v TEtrimmer >/dev/null 2>&1 || {
        echo "[ERROR] TEtrimmer not found in environment." >&2
        exit 1
      }
    '
  fi

  if ! { repeat_checkpoint_ready earlgrey; }; then
    echo "[INFO] Auditing Earl Grey environment"

    env_run "$ENV_EARLGREY" bash -lc '
      set -Eeuo pipefail

      command -v earlGreyAnnotationOnly >/dev/null 2>&1 || {
        echo "[ERROR] earlGreyAnnotationOnly not found in environment." >&2
        exit 1
      }
    '
  fi

  if is_true "$EARLGREY_ENHANCED_PLOTS" && ! { repeat_checkpoint_ready earlgrey_plots; }; then
    echo "[INFO] Auditing R environment for Earl Grey enhanced plots"

    env_run "$ENV_R" Rscript --vanilla -e 'if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Renv is missing required R package: ggplot2"); cat("Renv ggplot2 OK\n")' \
      > "${DIR_LOGS}/05B_earlgrey_enhanced_plots.renv_preflight.log" \
      2>&1 \
      || error "Renv preflight failed for Earl Grey enhanced plots; install ggplot2 in ${ENV_R}"
  fi

# ============================================================
# STEP 0 — PREPARE DFAM + REPBASE KNOWN LIBRARY
# ============================================================

if repeat_checkpoint_ready known; then
    echo "[INFO] Step 0: existing known-repeat library found; skipping preparation: $KNOWN_NR"
  else
    echo "[INFO] Step 0: Preparing Dfam + RepBase known-repeat library"

  echo "[INFO] Exporting Dfam FamDB to FASTA"

  singularity exec \
    --no-home \
    -B "$PROJECT_ROOT":"$PROJECT_ROOT" \
    -B "$DFAM_FAMDB_DIR":"$DFAM_FAMDB_DIR" \
    "$sif_repeatmasker" \
    bash -lc "
      set -Eeuo pipefail

      famdb.py \
        -i '$DFAM_FAMDB_DIR' \
        families \
        -ad \
        $DFAM_FLAGS \
        --format fasta_name \
        --include-class-in-name \
        '$DFAM_TAXON' \
        > '$DFAM_FASTA'
    " \
    > "${DIR_LOGS}/00A_dfam_export.log" \
    2>&1

  require_file "$DFAM_FASTA" "Dfam FASTA export"

  echo "[INFO] Converting RepBase EMBL to FASTA using tolerant parser"

  embl_to_fa_tolerant_host() {
    awk -v LOG="/dev/stderr" -v OFS="" '
      BEGIN {
        IGNORECASE=1
        id=""
        de=""
        seq=""
        inseq=0
        out=0
      }

      function emit_record(    name, j) {
        if (seq != "") {
          gsub(/[ \t\r]/, "", seq)
          seq=toupper(seq)
          name=(id != "" ? id : "unknown")

          if (de != "") {
            print ">", name, " ", de
          } else {
            print ">", name
          }

          for (j=1; j<=length(seq); j+=80) {
            print substr(seq, j, 80)
          }
          out++
        } else if (id != "") {
          print "[WARN] Skipping RepBase record with no valid sequence: " id > LOG
        }

        id=""
        de=""
        seq=""
        inseq=0
      }

      /^\/\// {
        emit_record()
        next
      }

      /^ID[ \t]/ {
        line=$0
        sub(/^ID[ \t]+/, "", line)
        split(line, F, /[ \t;]+/)
        id=F[1]
        next
      }

      /^DE[ \t]/ {
        line=$0
        sub(/^DE[ \t]+/, "", line)
        if (de != "") de=de " "
        de=de line
        next
      }

      /^SQ[ \t]/ {
        inseq=1
        next
      }

      inseq == 1 {
        raw=$0
        clean=$0
        gsub(/[^A-Za-z]/, "", clean)

        if (clean == "") {
          next
        }

        if (clean ~ /^[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+$/) {
          seq = seq clean
        } else {
          if (id == "") id="unknown_record"
          print "[WARN] Ignoring malformed non-IUPAC RepBase SQ line in " id ": " raw > LOG
        }
      }

      END {
        emit_record()
        if (out < 1) {
          print "[ERROR] RepBase parser produced no FASTA records" > LOG
          exit 2
        }
      }
    '
  }

  mkdir -p "$(dirname "$REPBASE_FASTA")" "$DIR_LOGS"
  rm -f "$REPBASE_FASTA"
  REPBASE_FASTA_TMP="${REPBASE_FASTA}.tmp.$$"
  rm -f "$REPBASE_FASTA_TMP"

  if [[ "$REPBASE_EMBL" == *.gz ]]; then
    gunzip -c "$REPBASE_EMBL" | embl_to_fa_tolerant_host > "$REPBASE_FASTA_TMP" \
      2> "${DIR_LOGS}/00B_repbase_convert.log" \
      || error "RepBase FASTA conversion failed"
  else
    embl_to_fa_tolerant_host < "$REPBASE_EMBL" > "$REPBASE_FASTA_TMP" \
      2> "${DIR_LOGS}/00B_repbase_convert.log" \
      || error "RepBase FASTA conversion failed"
  fi
  require_file "$REPBASE_FASTA_TMP" "temporary RepBase FASTA conversion"
  mv -f "$REPBASE_FASTA_TMP" "$REPBASE_FASTA" \
    || error "Failed to promote RepBase FASTA conversion to $REPBASE_FASTA"

  require_file "$REPBASE_FASTA" "RepBase FASTA conversion"

  REPBASE_RECORDS=$(grep -c "^>" "$REPBASE_FASTA" || true)
  if [[ "$REPBASE_RECORDS" -lt 1 ]]; then
    echo "[ERROR] RepBase conversion produced no FASTA records: $REPBASE_FASTA" >&2
    echo "[ERROR] See log: ${DIR_LOGS}/00B_repbase_convert.log" >&2
    exit 1
  fi
  echo "[INFO] RepBase FASTA records: $REPBASE_RECORDS"

  echo "[INFO] Merging and deduplicating known-repeat library"

  cat "$DFAM_FASTA" "$REPBASE_FASTA" > "$KNOWN_RAW"

  if [[ -n "$ARTEFACTS_FASTA" && -s "$ARTEFACTS_FASTA" ]]; then
    cat "$ARTEFACTS_FASTA" >> "$KNOWN_RAW"
  fi

  singularity exec \
    --no-home \
    -B "$PROJECT_ROOT":"$PROJECT_ROOT" \
    "$sif_repeatmasker" \
    bash -lc "
      set -Eeuo pipefail

      cd-hit-est \
        -i '$KNOWN_RAW' \
        -o '$KNOWN_NR' \
        -c '$KNOWN_CDHIT_ID' \
        -n '$KNOWN_CDHIT_WORD' \
        -T '$THREADS' \
        -M 0 \
        -d 0

      if command -v makeblastdb >/dev/null 2>&1; then
        makeblastdb -dbtype nucl -in '$KNOWN_NR'
      else
        /opt/rmblast/bin/makeblastdb -dbtype nucl -in '$KNOWN_NR'
      fi
    " \
    > "${DIR_LOGS}/00C_known_library_cdhit.log" \
    2>&1

  require_file "$KNOWN_NR" "known-repeat library"
  repeat_mark_checkpoint known
fi

echo "[INFO] Known-repeat library: $KNOWN_NR"

# ============================================================
# STEP 1 — REPEATMODELER2
# ============================================================

if repeat_checkpoint_ready repeatmodeler; then
    echo "[INFO] Step 1: existing RepeatModeler library found; skipping RepeatModeler: $RM_LIBRARY"
  else
    echo "[INFO] Step 1: Running RepeatModeler2"

    archive_unrecoverable_repeatmodeler_dirs

  # Build RepeatModeler option string inside the host shell, then pass it to
  # the container. This allows recovery runs to use the exact same settings.
  RM_OPTS="-database '$GENOME_STEM' -threads '$THREADS'"

  if [[ "$USE_LTRSTRUCT" == "yes" ]]; then
    RM_OPTS="$RM_OPTS -LTRStruct"
  fi

  if [[ "$RM_ADDL_ROUNDS" -gt 0 ]]; then
    RM_OPTS="$RM_OPTS -numAddlRounds '$RM_ADDL_ROUNDS'"
  fi

  {
    echo "# RepeatModeler command logic"
    echo "BuildDatabase -name '$GENOME_STEM' '$GENOME_LOCAL'"
    echo "RepeatModeler $RM_OPTS"
    echo "REPEATMODELER_RECOVER=$REPEATMODELER_RECOVER"
    echo "REPEATMODELER_RECOVER_DIR=$REPEATMODELER_RECOVER_DIR"
  } > "${DIR_LOGS}/01_repeatmodeler.command.sh"

  singularity exec \
    --no-home \
    -B "$PROJECT_ROOT":"$PROJECT_ROOT" \
    "$sif_repeatmasker" \
    bash -lc "
      set -Eeuo pipefail
      cd '$DIR_RM'

      BuildDatabase \
        -name '$GENOME_STEM' \
        '$GENOME_LOCAL'

      newest_rm_dir=''

      if [[ -n '$REPEATMODELER_RECOVER_DIR' ]]; then
        newest_rm_dir='$REPEATMODELER_RECOVER_DIR'
      elif [[ '$REPEATMODELER_RECOVER' == 'TRUE' ]]; then
        newest_rm_dir=\$(find . -maxdepth 1 -type d -name 'RM_*' -printf '%T@\t%p\n' 2>/dev/null | sort -nr | head -n1 | cut -f2- || true)
      fi

      if [[ -n \"\$newest_rm_dir\" && -d \"\$newest_rm_dir\" ]]; then
        echo '[INFO] Attempting RepeatModeler recovery from:' \"\$newest_rm_dir\" >&2

        set +e
        RepeatModeler $RM_OPTS -recoverDir \"\$newest_rm_dir\"
        rc=\$?
        set -e

        if [[ \$rc -ne 0 ]]; then
          echo '[WARN] RepeatModeler recovery failed with exit code:' \$rc >&2
          echo '[WARN] Starting a fresh RepeatModeler run.' >&2
          RepeatModeler $RM_OPTS
        fi
      else
        RepeatModeler $RM_OPTS
      fi
    " \
    > "${DIR_LOGS}/01_repeatmodeler.log" \
    2>&1

  if [[ ! -s "$RM_LIBRARY" ]]; then
    RM_LIBRARY="$(find "$DIR_RM" -type f -name '*families.fa' | head -n 1 || true)"
  fi

  if [[ ! -s "$RM_LIBRARY" ]]; then
    echo "[ERROR] RepeatModeler consensus library not found." >&2
    echo "[ERROR] See log: ${DIR_LOGS}/01_repeatmodeler.log" >&2
    echo "[ERROR] Useful diagnostics:" >&2
    echo "        cd '$DIR_RM'" >&2
    echo "        find . -type f \\( -name '*.log' -o -name '*.err' -o -name 'stderr' -o -name 'stdout' \\) -exec grep -Hn 'eledef\\|RECON\\|failed\\|error\\|Exit code' {} \\; | tail -n 100" >&2
    exit 1
  fi
  repeat_mark_checkpoint repeatmodeler
fi

echo "[INFO] RepeatModeler library: $RM_LIBRARY"

# ============================================================
# STEP 2 — TETRIMMER
# ============================================================

if repeat_checkpoint_ready tetrimmer; then
    echo "[INFO] Step 2: existing TEtrimmer library found; skipping TEtrimmer: $TETRIMMER_LIBRARY"
  else
    echo "[INFO] Step 2: Running TEtrimmer"
  echo "[INFO] TEtrimmer mode:   $TETRIMMER_MODE"
  echo "[INFO] TEtrimmer preset: $TETRIMMER_PRESET"

  TETRIMMER_ARGS=(
    --input_file "$RM_LIBRARY"
    --genome_file "$GENOME_LOCAL"
    --output_dir "$DIR_TETRIMMER"
    --num_threads "$THREADS"
  )

  if [[ -n "$TETRIMMER_PRESET" ]]; then
    case "$TETRIMMER_PRESET" in
      conserved|divergent)
        TETRIMMER_ARGS+=(--preset "$TETRIMMER_PRESET")
        ;;
      *)
        echo "[ERROR] Invalid TETRIMMER_PRESET=$TETRIMMER_PRESET. Use conserved, divergent, or empty." >&2
        exit 1
        ;;
    esac
  fi

  prepare_tetrimmer_output_dir

  case "$TETRIMMER_CONTINUE_MODE" in
    auto)
      # Conservative auto mode: only complete outputs are skipped by AUTO_RESUME.
      # Incomplete directories are archived above and TEtrimmer starts fresh.
      echo "[INFO] TEtrimmer continue mode is auto; starting fresh unless manually set to always."
      ;;
    always)
      if is_true "$TETRIMMER_CONTINUE"; then
        echo "[INFO] TEtrimmer continue mode is always; using --continue_analysis"
        TETRIMMER_ARGS+=(--continue_analysis)
      fi
      ;;
    never)
      echo "[INFO] TEtrimmer continue mode is never; not using --continue_analysis"
      ;;
  esac

  if [[ "$TETRIMMER_MODE" == "curated" ]]; then
    EFFECTIVE_TETRIMMER_CLASSIFICATION="classify_all"
  else
    EFFECTIVE_TETRIMMER_CLASSIFICATION="$TETRIMMER_CLASSIFICATION"
  fi

  case "$EFFECTIVE_TETRIMMER_CLASSIFICATION" in
    classify_unknown)
      TETRIMMER_ARGS+=(--classify_unknown)
      ;;
    classify_all)
      TETRIMMER_ARGS+=(--classify_all)
      ;;
    none)
      ;;
    *)
      echo "[ERROR] Invalid TETRIMMER_CLASSIFICATION=$EFFECTIVE_TETRIMMER_CLASSIFICATION." >&2
      echo "[ERROR] Use classify_unknown, classify_all, or none." >&2
      exit 1
      ;;
  esac

  if is_true "$TETRIMMER_DEDUP"; then
    TETRIMMER_ARGS+=(--dedup)
  fi

  if is_true "$TETRIMMER_HMM"; then
    TETRIMMER_ARGS+=(--hmm)
  fi

  if is_true "$TETRIMMER_GENOME_ANNO"; then
    TETRIMMER_ARGS+=(--genome_anno)
  fi

  if [[ "$TETRIMMER_MODE" == "debug" ]] || is_true "$TETRIMMER_DEBUG"; then
    TETRIMMER_ARGS+=(--debug)
  fi

  if [[ -n "$TETRIMMER_PFAM_DIR" ]]; then
    TETRIMMER_ARGS+=(--pfam_dir "$TETRIMMER_PFAM_DIR")
  fi

  TETRIMMER_ARGS+=(
    --cons_thr "$TETRIMMER_CONS_THR"
    --mini_orf "$TETRIMMER_MINI_ORF"
    --max_msa_lines "$TETRIMMER_MAX_MSA_LINES"
    --top_msa_lines "$TETRIMMER_TOP_MSA_LINES"
    --min_seq_num "$TETRIMMER_MIN_SEQ_NUM"
    --min_blast_len "$TETRIMMER_MIN_BLAST_LEN"
    --max_cluster_num "$TETRIMMER_MAX_CLUSTER_NUM"
    --ext_thr "$TETRIMMER_EXT_THR"
    --ext_check_win "$TETRIMMER_EXT_CHECK_WIN"
    --ext_step "$TETRIMMER_EXT_STEP"
    --max_ext "$TETRIMMER_MAX_EXT"
    --gap_thr "$TETRIMMER_GAP_THR"
    --gap_nul_thr "$TETRIMMER_GAP_NUL_THR"
    --crop_end_div_thr "$TETRIMMER_CROP_END_DIV_THR"
    --crop_end_div_win "$TETRIMMER_CROP_END_DIV_WIN"
    --crop_end_gap_thr "$TETRIMMER_CROP_END_GAP_THR"
    --crop_end_gap_win "$TETRIMMER_CROP_END_GAP_WIN"
    --ltr_start_patterns "$TETRIMMER_LTR_START_PATTERNS"
    --ltr_end_patterns "$TETRIMMER_LTR_END_PATTERNS"
    --helitron_start_patterns "$TETRIMMER_HELITRON_START_PATTERNS"
    --helitron_end_patterns "$TETRIMMER_HELITRON_END_PATTERNS"
    --poly_patterns "$TETRIMMER_POLY_PATTERNS"
    --poly_len "$TETRIMMER_POLY_LEN"
  )

  case "$TETRIMMER_MODE" in
    fast|curated|debug)
      ;;
    *)
      echo "[ERROR] Invalid TETRIMMER_MODE=$TETRIMMER_MODE. Use fast, curated, or debug." >&2
      exit 1
      ;;
  esac

  {
    echo "# TEtrimmer command"
    printf '%q run -n %q TEtrimmer' "$ENV_RUNNER" "$ENV_TETRIMMER"
    printf ' %q' "${TETRIMMER_ARGS[@]}"
    echo
  } > "${DIR_LOGS}/02_tetrimmer.command.sh"

  env_run "$ENV_TETRIMMER" \
    TEtrimmer "${TETRIMMER_ARGS[@]}" \
    > "${DIR_LOGS}/02_tetrimmer.log" \
    2>&1

  audit_tetrimmer_outputs

  if [[ -s "${DIR_TETRIMMER}/TEtrimmer_consensus_merged.fasta" ]]; then
    TETRIMMER_LIBRARY="${DIR_TETRIMMER}/TEtrimmer_consensus_merged.fasta"
  elif [[ -s "${DIR_TETRIMMER}/TEtrimmer_consensus.fasta" ]]; then
    echo "[WARN] TEtrimmer_consensus_merged.fasta not found after TEtrimmer run."
    echo "[WARN] Falling back to TEtrimmer_consensus.fasta."
    TETRIMMER_LIBRARY="${DIR_TETRIMMER}/TEtrimmer_consensus.fasta"
  else
    TETRIMMER_LIBRARY=""
  fi

  require_file "$TETRIMMER_LIBRARY" "TEtrimmer refined library"
  repeat_mark_checkpoint tetrimmer
fi

echo "[INFO] TEtrimmer library: $TETRIMMER_LIBRARY"


# ============================================================
# STEP 3 — FINAL LIBRARY
# ============================================================

if repeat_checkpoint_ready final_library; then
    echo "[INFO] Step 3: existing final TE library found; skipping final-library construction: $FINAL_NR"
  else
    echo "[INFO] Step 3: Building final annotation library"

  cat "$TETRIMMER_LIBRARY" "$KNOWN_NR" > "$FINAL_RAW"

  singularity exec \
    --no-home \
    -B "$PROJECT_ROOT":"$PROJECT_ROOT" \
    "$sif_repeatmasker" \
    bash -lc "
      set -Eeuo pipefail

      cd-hit-est \
        -i '$FINAL_RAW' \
        -o '$FINAL_NR' \
        -c '$FINAL_CDHIT_ID' \
        -n '$FINAL_CDHIT_WORD' \
        -T '$THREADS' \
        -M 0 \
        -d 0

      if command -v makeblastdb >/dev/null 2>&1; then
        makeblastdb -dbtype nucl -in '$FINAL_NR'
      else
        /opt/rmblast/bin/makeblastdb -dbtype nucl -in '$FINAL_NR'
      fi
    " \
    > "${DIR_LOGS}/03_final_library.log" \
    2>&1

  require_file "$FINAL_NR" "final TE annotation library"
  repeat_mark_checkpoint final_library
fi

echo "[INFO] Final TE library: $FINAL_NR"


# ============================================================
# STEP 4 — REPEATMASKER
# ============================================================

if repeat_checkpoint_ready repeatmasker; then
    echo "[INFO] Step 4: existing RepeatMasker outputs found; skipping RepeatMasker"
  else
    echo "[INFO] Step 4: Running RepeatMasker"

  RM_ARGS=(
    -e "$RM_ENGINE"
    -pa "$RM_PA"
    -lib "$FINAL_NR"
    -dir "$DIR_REPEATMASKER"
  )

  if [[ "$RM_GFF" == "yes" ]]; then
    RM_ARGS+=(-gff)
  fi

  if [[ "$RM_XSMALL" == "yes" ]]; then
    RM_ARGS+=(-xsmall)
  fi

  if [[ "$RM_ALIGN" == "yes" ]]; then
    RM_ARGS+=(-a)
  fi

  if [[ -n "$RM_SENSITIVITY" ]]; then
    RM_ARGS+=("$RM_SENSITIVITY")
  fi

  {
    echo "# RepeatMasker command"
    printf 'singularity exec --no-home -B %q:%q %q RepeatMasker' \
      "$PROJECT_ROOT" "$PROJECT_ROOT" "$sif_repeatmasker"
    printf ' %q' "${RM_ARGS[@]}"
    printf ' %q\n' "$GENOME_LOCAL"
  } > "${DIR_LOGS}/04_repeatmasker.command.sh"

  singularity exec \
    --no-home \
    -B "$PROJECT_ROOT":"$PROJECT_ROOT" \
    "$sif_repeatmasker" \
    RepeatMasker "${RM_ARGS[@]}" "$GENOME_LOCAL" \
    > "${DIR_LOGS}/04_repeatmasker.log" \
    2>&1

  repeat_mark_checkpoint repeatmasker
fi

# ============================================================
# STEP 5 — EARL GREY ANNOTATION-ONLY
# ============================================================

if repeat_checkpoint_ready earlgrey; then
    echo "[INFO] Step 5: existing Earl Grey outputs found; skipping Earl Grey"
  else
    echo "[INFO] Step 5: Running Earl Grey annotation-only"

  EARLGREY_LIBRARY="$FINAL_NR"

  require_file "$EARLGREY_LIBRARY" "Earl Grey consensus library"

  EARLGREY_ARGS=(
    -g "$GENOME_LOCAL"
    -s "$SPECIES"
    -o "$DIR_EARLGREY"
    -l "$EARLGREY_LIBRARY"
    -t "$THREADS"
    -m "$EARLGREY_REMOVE_SHORT"
    -d "$EARLGREY_SOFTMASK"
    -e "$EARLGREY_HELIANO"
  )

  if [[ "$EARLGREY_REPEATMASKER_SPECIES" != "none" && -n "$EARLGREY_REPEATMASKER_SPECIES" ]]; then
    EARLGREY_ARGS+=(-r "$EARLGREY_REPEATMASKER_SPECIES")
  fi

  {
    echo "# Earl Grey annotation-only command"
    printf '%q run -n %q earlGreyAnnotationOnly' "$ENV_RUNNER" "$ENV_EARLGREY"
    printf ' %q' "${EARLGREY_ARGS[@]}"
    echo
  } > "${DIR_LOGS}/05_earlgrey_annotation_only.command.sh"

  env_run "$ENV_EARLGREY" \
    earlGreyAnnotationOnly "${EARLGREY_ARGS[@]}" \
    > "${DIR_LOGS}/05_earlgrey_annotation_only.log" \
    2>&1

  repeat_mark_checkpoint earlgrey
fi

# ============================================================
# STEP 6 — EARL GREY ENHANCED GGPLOT2 PLOTS
# ============================================================

if is_true "$EARLGREY_ENHANCED_PLOTS"; then
  if repeat_checkpoint_ready earlgrey_plots; then
    echo "[INFO] Step 6: existing Earl Grey enhanced plots found; skipping"
  else
    echo "[INFO] Step 6: Creating Earl Grey enhanced ggplot2 plots"

    EARLGREY_PLOTTER="${DIR_LOGS}/05B_earlgrey_enhanced_plots.R"
    write_earlgrey_enhanced_plotter "$EARLGREY_PLOTTER"

    {
      echo "# Earl Grey enhanced ggplot2 plotting command"
      printf '%q run -n %q Rscript --vanilla %q' "$ENV_RUNNER" "$ENV_R" "$EARLGREY_PLOTTER"
      printf ' --species %q --earlgrey_dir %q --outdir %q --dpi %q --top_families %q\n' \
        "$SPECIES" "$DIR_EARLGREY" "$DIR_EARLGREY_ENHANCED_PLOTS" "$EARLGREY_PLOT_DPI" "$EARLGREY_PLOT_TOP_FAMILIES"
    } > "${DIR_LOGS}/05B_earlgrey_enhanced_plots.command.sh"

    env_run "$ENV_R" Rscript --vanilla "$EARLGREY_PLOTTER" \
      --species "$SPECIES" \
      --earlgrey_dir "$DIR_EARLGREY" \
      --outdir "$DIR_EARLGREY_ENHANCED_PLOTS" \
      --dpi "$EARLGREY_PLOT_DPI" \
      --top_families "$EARLGREY_PLOT_TOP_FAMILIES" \
      > "${DIR_LOGS}/05B_earlgrey_enhanced_plots.log" \
      2>&1

    repeat_mark_checkpoint earlgrey_plots
  fi
else
  echo "[INFO] Step 6: Earl Grey enhanced plots disabled"
fi

# ============================================================
# STEP 7 — COLLECT OUTPUTS
# ============================================================

if repeat_checkpoint_ready collect_outputs; then
    echo "[INFO] Step 7: existing collected repeat outputs found; skipping collection"
  else
  echo "[INFO] Step 7: Collecting final outputs"

  mkdir -p \
    "$DIR_FINAL/known_libraries" \
    "$DIR_FINAL/repeatmodeler" \
    "$DIR_FINAL/tetrimmer" \
    "$DIR_FINAL/final_library" \
    "$DIR_FINAL/repeatmasker" \
    "$DIR_FINAL/earlgrey" \
    "$DIR_FINAL/logs"

  if [[ -s "${DFAM_FASTA:-}" ]]; then cp -f "$DFAM_FASTA" "$DIR_FINAL/known_libraries/"; fi
  if [[ -s "${REPBASE_FASTA:-}" ]]; then cp -f "$REPBASE_FASTA" "$DIR_FINAL/known_libraries/"; fi
  if [[ -s "${KNOWN_RAW:-}" ]]; then cp -f "$KNOWN_RAW" "$DIR_FINAL/known_libraries/"; fi
  if [[ -s "${KNOWN_NR:-}" ]]; then cp -f "$KNOWN_NR"* "$DIR_FINAL/known_libraries/"; fi

  if [[ -s "${RM_LIBRARY:-}" ]]; then cp -f "$RM_LIBRARY" "$DIR_FINAL/repeatmodeler/"; fi
  # TEtrimmer main libraries and reports
  if [[ -s "${TETRIMMER_LIBRARY:-}" ]]; then cp -f "$TETRIMMER_LIBRARY" "$DIR_FINAL/tetrimmer/"; fi

  for f in \
    "$DIR_TETRIMMER/TEtrimmer_consensus.fasta" \
    "$DIR_TETRIMMER/TEtrimmer_consensus_merged.fasta" \
    "$DIR_TETRIMMER/summary.txt" \
    "$DIR_TETRIMMER/Sequence_name_mapping.txt" \
    "$DIR_TETRIMMER/error_file.txt" \
    "$DIR_LOGS/02_tetrimmer.output_audit.txt"
  do
    if [[ -s "$f" ]]; then cp -f "$f" "$DIR_FINAL/tetrimmer/"; fi
  done

  # TEtrimmer important output directories
  for d in \
    Classification \
    Single_fasta_files \
    TEtrimmer_for_proof_curation \
    Annotation_perfect \
    Annotation_good \
    Annotation_check_recommended \
    Annotation_check_required \
    Clustered_proof_curation \
    TE_low_copy \
    TE_skipped \
    HMM \
    Multiple_sequence_alignment
  do
    if [[ -d "$DIR_TETRIMMER/$d" ]]; then
      mkdir -p "$DIR_FINAL/tetrimmer/$d"
      find "$DIR_TETRIMMER/$d" \
        -maxdepth 2 \
        -type f \
        \( \
          -name "*.fa" \
          -o -name "*.fasta" \
          -o -name "*.pdf" \
          -o -name "*.txt" \
          -o -name "*.tsv" \
          -o -name "*.csv" \
          -o -name "*.hmm" \
          -o -name "*.log" \
        \) \
        -exec cp -f {} "$DIR_FINAL/tetrimmer/$d/" \; 2>/dev/null || true
    fi
  done

  if [[ -s "${FINAL_RAW:-}" ]]; then cp -f "$FINAL_RAW" "$DIR_FINAL/final_library/"; fi
  if [[ -s "${FINAL_NR:-}" ]]; then cp -f "$FINAL_NR"* "$DIR_FINAL/final_library/"; fi

  find "$DIR_REPEATMASKER" \
    -maxdepth 1 \
    -type f \
    \( \
      -name "*.out" \
      -o -name "*.gff" \
      -o -name "*.tbl" \
      -o -name "*.masked" \
      -o -name "*.align" \
      -o -name "*.cat" \
    \) \
    -exec cp -f {} "$DIR_FINAL/repeatmasker/" \; 2>/dev/null || true

  find "$DIR_EARLGREY" \
    -type f \
    \( \
      -name "*.gff" \
      -o -name "*.gff3" \
      -o -name "*.tbl" \
      -o -name "*.summary" \
      -o -name "*.html" \
      -o -name "*.png" \
      -o -name "*.pdf" \
      -o -name "*.txt" \
      -o -name "*.tsv" \
      -o -name "*.csv" \
      -o -name "*.masked" \
      -o -name "*.fa" \
      -o -name "*.fasta" \
    \) \
    -exec cp -f {} "$DIR_FINAL/earlgrey/" \; 2>/dev/null || true

	  cp -f "$DIR_LOGS"/*.log "$DIR_FINAL/logs/" 2>/dev/null || true
	  cp -f "$DIR_LOGS"/*.command.sh "$DIR_FINAL/logs/" 2>/dev/null || true

	  collect_manifest="${DIR_FINAL}/logs/collect_outputs_manifest.tsv"
	  {
	    printf "relative_path\tsize_bytes\n"
	    while IFS= read -r -d '' f; do
	      if [[ -s "$f" ]]; then
	        rel="${f#"${DIR_FINAL}"/}"
	        size="$(wc -c < "$f" 2>/dev/null || printf "0")"
	        printf "%s\t%s\n" "$rel" "$size"
	      fi
	    done < <(find "$DIR_FINAL" \( -type f -o -type l \) -print0 2>/dev/null | sort -z)
	  } > "$collect_manifest"

	  awk 'NR > 1 { found = 1 } END { exit found ? 0 : 1 }' "$collect_manifest" \
	    || error "Collect outputs did not copy any non-empty files into ${DIR_FINAL}"

	  repeat_mark_checkpoint collect_outputs
	fi

	# ============================================================
	# STEP 8 — NON-REPETITIVE REGIONS BED
	# ============================================================

	if repeat_checkpoint_ready non_repeats; then
	    echo "[INFO] Step 8: existing non-repetitive regions BED found; skipping"
	  else
	  echo "[INFO] Step 8: Creating non-repetitive regions BED"

	  require_file "$GENOME_LOCAL" "genome for non-repetitive BED"
	  require "$bedtools_bin"

	  REPEATS_BED="${NON_REPEATS_REPEATS_BED:-}"
	  if [[ -z "$REPEATS_BED" ]]; then
	    expected_earlgrey_bed="${DIR_EARLGREY}/${SPECIES}_EarlGrey/${SPECIES}_summaryFiles/${SPECIES}.filteredRepeats.bed"
	    if [[ -s "$expected_earlgrey_bed" ]]; then
	      REPEATS_BED="$expected_earlgrey_bed"
	    else
	      REPEATS_BED="$(find "$DIR_EARLGREY" -type f -name "*.filteredRepeats.bed" -size +0c -print -quit 2>/dev/null || true)"
	    fi
	  fi

	  require_file "$REPEATS_BED" "Earl Grey filtered repeats BED"

	  if [[ ! -s "${GENOME_LOCAL}.fai" ]]; then
	    "$samtools_bin" faidx "$GENOME_LOCAL" \
	      > "${DIR_LOGS}/07_non_repeats.samtools_faidx.log" \
	      2>&1 \
	      || error "samtools faidx failed for non-repetitive BED input: $GENOME_LOCAL"
	  fi

	  cut -f1,2 "${GENOME_LOCAL}.fai" > "$NON_REPEATS_GENOME_SIZES"

	  awk 'BEGIN{OFS="\t"} NF >= 3 && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ && $2 < $3 {print $1,$2,$3}' "$REPEATS_BED" \
	    > "$NON_REPEATS_RAW_BED"

	  "$bedtools_bin" sort \
	    -i "$NON_REPEATS_RAW_BED" \
	    -g "$NON_REPEATS_GENOME_SIZES" \
	    > "$NON_REPEATS_SORTED_BED"

	  "$bedtools_bin" merge \
	    -i "$NON_REPEATS_SORTED_BED" \
	    > "$NON_REPEATS_MERGED_BED"

	  "$bedtools_bin" complement \
	    -i "$NON_REPEATS_MERGED_BED" \
	    -g "$NON_REPEATS_GENOME_SIZES" \
	    > "$NON_REPEATS_BED"

	  {
	    echo "# Non-repetitive region BED command"
	    printf 'cut -f1,2 %q > %q\n' "${GENOME_LOCAL}.fai" "$NON_REPEATS_GENOME_SIZES"
	    # shellcheck disable=SC2016
	    printf 'awk %q %q > %q\n' 'BEGIN{OFS="\t"} NF >= 3 && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ && $2 < $3 {print $1,$2,$3}' "$REPEATS_BED" "$NON_REPEATS_RAW_BED"
	    printf '%q sort -i %q -g %q > %q\n' "$bedtools_bin" "$NON_REPEATS_RAW_BED" "$NON_REPEATS_GENOME_SIZES" "$NON_REPEATS_SORTED_BED"
	    printf '%q merge -i %q > %q\n' "$bedtools_bin" "$NON_REPEATS_SORTED_BED" "$NON_REPEATS_MERGED_BED"
	    printf '%q complement -i %q -g %q > %q\n' "$bedtools_bin" "$NON_REPEATS_MERGED_BED" "$NON_REPEATS_GENOME_SIZES" "$NON_REPEATS_BED"
	  } > "${DIR_LOGS}/07_non_repeats.command.sh"

	  mkdir -p "$DIR_FINAL/non_repetitive_regions"
	  cp -f "$NON_REPEATS_GENOME_SIZES" "$NON_REPEATS_RAW_BED" "$NON_REPEATS_SORTED_BED" "$NON_REPEATS_MERGED_BED" "$NON_REPEATS_BED" \
	    "$DIR_FINAL/non_repetitive_regions/"
	  cp -f "${DIR_LOGS}/07_non_repeats.command.sh" "$DIR_FINAL/logs/" 2>/dev/null || true

	  repeat_mark_checkpoint non_repeats
	fi

	  echo "[INFO] Step 8 complete. Repeat outputs: ${DIR_FINAL}"
	fi

# ============================================================================================================================
# STEP 9: Final Summary & Run Report — execution stats, resources, key metrics, artifacts
# Produces:
#   - ${dir7}/run_summary.txt
#   - ${dir7}/run_summary.tsv
#   - ${dir7}/run_summary.json
#   - ${dir7}/logs/summary_tail.log
# ============================================================================================================================

echo "[INFO] Step 9: Final Summary — starting"
mkdir -p "${dir7}" "${dir7}/logs"

# --------------------------- helpers ---------------------------
now_ts() { date +%s; }
human_ts() { date -d "@$1" '+%F %T'; }
duration_hms() {
  local s="$1"; local d=$(( s/86400 )); local h=$(( (s%86400)/3600 )); local m=$(( (s%3600)/60 )); local sec=$(( s%60 ));
  printf "%d days, %02d:%02d:%02d" "$d" "$h" "$m" "$sec"
}
filesize() {
  local f="$1"
  [[ -s "$f" ]] && du -h --apparent-size "$f" 2>/dev/null | awk 'NR==1{print $1}' || echo "-"
}
contig_stats() {
  # Works for .fa or .fa.gz if zcat exists
  local fa="$1"
  [[ -s "$fa" ]] || { echo -e "-\t-\t-"; return; }
  local reader="cat"
  [[ "$fa" =~ \.gz$ ]] && reader="zcat"
  local n bp
  n=$($reader "$fa" | awk '/^>/{c++} END{print c+0}')
  bp=$($reader "$fa" | awk 'BEGIN{bp=0} /^>/{next} {bp+=length($0)} END{print bp}')
  echo -e "$n\t$bp\t$(filesize "$fa")"
}
# QUAST metric extractor: returns "-" if not found/empty; strips commas.
quast_metric() {
  local report="$1" asm_base="$2" key="$3"
  [[ -s "$report" ]] || { echo "-"; return; }
  # header fields → names array (skip column 1)
  IFS=$'\t' read -r -a _hdr < <(head -n1 "$report")
  local names=()
  for ((i=1;i<${#_hdr[@]};i++)); do names+=("${_hdr[$i]}"); done
  # find assembly column (try exact, then stem without .fa/.fasta)
  local idx=-1 stem="${asm_base%.fa}"; stem="${stem%.fasta}"
  for i in "${!names[@]}"; do [[ "${names[$i]}" == "$asm_base" ]] && idx=$i && break; done
  if (( idx < 0 )); then
    for i in "${!names[@]}"; do [[ "${names[$i]}" == "$stem" ]] && idx=$i && break; done
  fi
  (( idx < 0 )) && idx=$(( ${#names[@]} - 1 ))  # fallback: last column
  local v
  v=$(awk -v k="$key" -v col=$((idx+2)) -F'\t' '($1==k){print $col; exit}' "$report" | sed 's/,//g' || true)
  [[ -n "$v" ]] && echo "$v" || echo "-"
}
# BUSCO short summary parser (portable)
busco_summary() {
  local run_dir="$1" ss
  ss=$(compgen -G "$run_dir/short_summary*.txt" | head -n1 || true)
  [[ -n "$ss" && -s "$ss" ]] || { echo -e "-\t-\t-\t-\t-\t-"; return; }
  awk '
    { gsub(/,/, " "); for(i=1;i<=NF;i++){
        if ($i ~ /^C:/) C=$i; if ($i ~ /^S:/) S=$i; if ($i ~ /^D:/) D=$i;
        if ($i ~ /^F:/) F=$i; if ($i ~ /^M:/) M=$i; if ($i ~ /^n:/) n=$i;
      }}
    END{
      sub(/.*:/,"",C); sub(/.*:/,"",S); sub(/.*:/,"",D);
      sub(/.*:/,"",F); sub(/.*:/,"",M); sub(/.*:/,"",n);
      gsub(/[%\]]/,"",C); gsub(/[%\]]/,"",S); gsub(/[%\]]/,"",D);
      gsub(/[%\]]/,"",F); gsub(/[%\]]/,"",M);
      printf "%s\t%s\t%s\t%s\t%s\t%s\n", C,S,D,F,M,n;
    }' "$ss"
}

# --------------------------- timing & system ---------------------------
[[ -n "${start_time:-}" ]] || start_time=$(now_ts)
end_time=$(now_ts)
runtime=$(( end_time - start_time ))
start_h=$(human_ts "$start_time")
end_h=$(human_ts "$end_time")
runtime_h=$(duration_hms "$runtime")

host=$(hostname)
kernel=$(uname -srm)
cpu_model=$(awk -F: '/model name/{print $2; exit}' /proc/cpuinfo | sed 's/^ //')
cpu_threads=$(nproc 2>/dev/null || echo "-")
mem_total_gb=$(awk '/MemTotal/ {printf "%.2f", $2/1024/1024}' /proc/meminfo 2>/dev/null || echo "0")

# Disk usage (human)
fs_usage=$(df -h "$working_dir" | awk 'NR==1 || NR==2')
run_du=$(du -sh "$working_dir" 2>/dev/null | awk '{print $1}')

# Per-step disk footprints
step_du=$(
  for d in "$dir1" "$dir2" "$dir3" "$dir4" "$dir5" "$dir6" "$dir7" "${repeats_dir}" "$dir8"; do
    [[ -d "$d" ]] && du -sh "$d"
  done | sort -h
)

# --------------------------- key artifacts / paths ---------------------------
# These dirs were created in Step 7
quast_dir="${dir7}/quast"
busco_dir="${dir7}/busco"
btk_dir="${dir7}/blobtoolkit"

: "${species_fasta:=${dir7}/${species_name}.fasta}"      # final sanitized assembly
quast_report="${quast_dir}/report.tsv"
busco_run_name="${busco_run_name:-final.polished.busco}"
busco_run_dir_base="${busco_dir}/${busco_run_name}"
busco_run_dir="$(find "${busco_run_dir_base}" -maxdepth 1 -type d -name 'run_*' -print -quit 2>/dev/null || true)"
[[ -n "${busco_run_dir}" ]] || busco_run_dir="${busco_run_dir_base}"
blobdir_path="${btk_dir}/blobdir"

# Final target (default to species-named FASTA if present)
final_target="${species_fasta}"
[[ -s "$final_target" ]] || final_target="${dir6}/final.polished.decontam.canonical.fa"
[[ -s "$final_target" ]] || final_target="${dir6}/final.polished.decontam.fa"
[[ -s "$final_target" ]] || final_target="${dir6}/final.polished.fa"

final_base="$(basename "$final_target")"
final_base="${final_base%.fa}"

# QUAST metrics against the column matching the final target
N50=$(quast_metric "$quast_report" "$final_base" "N50")
L50=$(quast_metric "$quast_report" "$final_base" "L50")
CONTIGS=$(quast_metric "$quast_report" "$final_base" "# contigs")
GC=$(quast_metric "$quast_report" "$final_base" "GC (%)")
Ns100=$(quast_metric "$quast_report" "$final_base" "# N's per 100 kbp")

# BUSCO headline stats
read -r BUSCO_C BUSCO_S BUSCO_D BUSCO_F BUSCO_M BUSCO_n < <(busco_summary "$busco_run_dir")

# Final assembly contig counts + size (human-readable)
read -r final_ctg_n final_ctg_bp final_ctg_size < <(contig_stats "$final_target")

summary_txt="${dir7}/run_summary.txt"
summary_tsv="${dir7}/run_summary.tsv"
summary_json="${dir7}/run_summary.json"

# Step 5 scaffolding selection artifacts.
scaffolding_scores="${dir5}/selection/longstitch_scores.tsv"
scaffolding_preselection="${dir5}/selection/preselection_scores.tsv"
scaffolding_best_path_file="${dir5}/best_assembly.path"
scaffolding_best_source="-"
[[ -s "${scaffolding_best_path_file}" ]] && scaffolding_best_source="$(head -n1 "${scaffolding_best_path_file}")"
tsv_first_value() {
  local file="${1:-}" key="${2:-}"
  [[ -s "${file}" && -n "${key}" ]] || { printf '%s\n' ""; return 0; }
  awk -F'\t' -v key="${key}" '
    NR == 1 { for (i = 1; i <= NF; i++) h[$i] = i; next }
    NR == 2 { if (key in h) print $h[key]; exit }
  ' "${file}"
}
scaffolding_best_id="$(tsv_first_value "${scaffolding_scores}" "assembly_id")"
scaffolding_best_k="$(tsv_first_value "${scaffolding_scores}" "k")"
scaffolding_best_w="$(tsv_first_value "${scaffolding_scores}" "w")"
scaffolding_score="$(tsv_first_value "${scaffolding_scores}" "score")"
scaffolding_n50="$(tsv_first_value "${scaffolding_scores}" "n50")"
scaffolding_l50="$(tsv_first_value "${scaffolding_scores}" "l50")"
scaffolding_contigs="$(tsv_first_value "${scaffolding_scores}" "scaffold_count")"
scaffolding_ns100="$(tsv_first_value "${scaffolding_scores}" "ns_per_100kbp")"
scaffolding_merqury_qv="$(tsv_first_value "${scaffolding_scores}" "merqury_qv")"
scaffolding_merqury_comp="$(tsv_first_value "${scaffolding_scores}" "merqury_completeness")"
scaffolding_yak_qv="$(tsv_first_value "${scaffolding_scores}" "yak_qv")"

# Step 8 repeat-annotation artifacts. These paths are stable even when Step 8 is skipped.
repeat_known_library="${repeats_dir}/01.known_libraries/known_repeats.nr.fa"
repeat_rm_library="$(find "${repeats_dir}/02.repeatmodeler" -type f -name '*families.fa' -print -quit 2>/dev/null || true)"
repeat_tetrimmer_library="${repeats_dir}/03.tetrimmer/TEtrimmer_consensus_merged.fasta"
[[ -s "$repeat_tetrimmer_library" ]] || repeat_tetrimmer_library="${repeats_dir}/03.tetrimmer/TEtrimmer_consensus.fasta"
repeat_final_library="$(find "${repeats_dir}/04.final_library" -type f \( -name '*TEtrimmer*Dfam*RepBase*nr*.fa' -o -name '*final*nr*.fa' \) -print -quit 2>/dev/null || true)"
repeat_repeatmasker_dir="${repeats_dir}/05.repeatmasker"
repeat_earlgrey_dir="${earlgrey_dir}"
repeat_final_outputs_dir="${repeats_dir}/07.final_outputs"
repeat_logs_dir="${repeats_dir}/logs"
repeat_checkpoints_dir="${repeats_dir}/checkpoints"
repeat_summary_status="not_run"
if is_true "${run_repeats:-FALSE}"; then
  repeat_summary_status="enabled"
  [[ -d "${repeat_final_outputs_dir}" ]] && repeat_summary_status="outputs_collected"
fi

# --------------------------- write human-readable (TOPLINE) ---------------------------
{
  echo "==================== PIPELINE RUN SUMMARY (Topline) ===================="
  echo "Final Assembly     : $(basename "$final_target")"
  echo "Location           : $final_target"
  echo "Size               : $final_ctg_size"
  echo "Contigs            : ${final_ctg_n:-'-'}"
  echo "Total bp           : ${final_ctg_bp:-'-'}"
  echo "Scaffolding        : best=${scaffolding_best_id:-'-'} | k=${scaffolding_best_k:-'-'} | w=${scaffolding_best_w:-'-'} | score=${scaffolding_score:-'-'}"
  echo "Scaffolding metrics: N50=${scaffolding_n50:-'-'} | L50=${scaffolding_l50:-'-'} | contigs=${scaffolding_contigs:-'-'} | Ns/100kb=${scaffolding_ns100:-'-'} | Merqury_QV=${scaffolding_merqury_qv:-'-'} | Merqury_comp=${scaffolding_merqury_comp:-'-'}"
  echo "QUAST              : N50=${N50:-'-'} | L50=${L50:-'-'} | #contigs=${CONTIGS:-'-'} | GC%=${GC:-'-'} | Ns/100kb=${Ns100:-'-'}"
  echo "BUSCO (${BUSCO_LINEAGE:-'-'}): C=${BUSCO_C:-'-'}% [S=${BUSCO_S:-'-'}% | D=${BUSCO_D:-'-'}%] | F=${BUSCO_F:-'-'}% | M=${BUSCO_M:-'-'}% | n=${BUSCO_n:-'-'}"
  echo "BlobToolKit        : ${blobdir_path} (plots: ${btk_dir}/plots)"
  echo "Repeats            : ${repeat_summary_status} | root=${repeats_dir}"
  echo "Repeat libraries   : known=${repeat_known_library:-'-'} | de_novo=${repeat_rm_library:-'-'} | refined=${repeat_tetrimmer_library:-'-'} | final=${repeat_final_library:-'-'}"
  echo "Repeat annotation  : RepeatMasker=${repeat_repeatmasker_dir} | EarlGrey=${repeat_earlgrey_dir}"
  echo "MultiQC            : ${dir7}/multiqc_report.html"
  echo "======================================================================="
  echo
  echo "========================= EXECUTION RESOURCES =========================="
  echo "Host               : $host"
  echo "OS/Kernel          : $kernel"
  echo "CPU                : $cpu_model"
  echo "CPU Threads (avail): $cpu_threads"
  echo "RAM (total)        : ${mem_total_gb} GB"
  echo "Start Time         : $start_h"
  echo "End Time           : $end_h"
  echo "Duration           : $runtime_h"
  echo
  echo "Working dir        : $working_dir"
  echo "Working dir size   : ${run_du:--}"
  echo "Filesystem usage   :"
  echo "$fs_usage"
  echo
  echo "Per-step disk usage:"
  echo "$step_du"
  echo
  echo "============================ ARTIFACT PATHS ============================"
  echo "QUAST report       : ${quast_report}"
  echo "Scaffolding scores : ${scaffolding_scores}"
  echo "Scaffolding source : ${scaffolding_best_source}"
  echo "BUSCO run dir      : ${busco_run_dir}"
  echo "BlobToolKit blobdir: ${blobdir_path}"
  echo "Repeats root       : ${repeats_dir}"
  echo "Repeats final out  : ${repeat_final_outputs_dir}"
  echo "Repeats checkpoints: ${repeat_checkpoints_dir}"
  echo "Repeats logs       : ${repeat_logs_dir}"
  echo "Logs (tail)        : ${dir7}/logs/summary_tail.log"
  echo "======================================================================="
} > "$summary_txt"

# --------------------------- write TSV ---------------------------
{
  printf "key\tvalue\n"
  printf "host\t%s\n" "$host"
  printf "kernel\t%s\n" "$kernel"
  printf "cpu_model\t%s\n" "$cpu_model"
  printf "cpu_threads\t%s\n" "$cpu_threads"
  printf "ram_total_gb\t%s\n" "$mem_total_gb"
  printf "start_time\t%s\n" "$start_h"
  printf "end_time\t%s\n" "$end_h"
  printf "duration\t%s\n" "$runtime_h"
  printf "working_dir\t%s\n" "$working_dir"
  printf "working_dir_size\t%s\n" "${run_du:-}"
  printf "final_asm\t%s\n" "$final_target"
  printf "final_contigs\t%s\n" "$final_ctg_n"
  printf "final_total_bp\t%s\n" "$final_ctg_bp"
  printf "scaffolding_best_id\t%s\n" "$scaffolding_best_id"
  printf "scaffolding_best_k\t%s\n" "$scaffolding_best_k"
  printf "scaffolding_best_w\t%s\n" "$scaffolding_best_w"
  printf "scaffolding_score\t%s\n" "$scaffolding_score"
  printf "scaffolding_N50\t%s\n" "$scaffolding_n50"
  printf "scaffolding_L50\t%s\n" "$scaffolding_l50"
  printf "scaffolding_contigs\t%s\n" "$scaffolding_contigs"
  printf "scaffolding_Ns_per_100kb\t%s\n" "$scaffolding_ns100"
  printf "scaffolding_merqury_qv\t%s\n" "$scaffolding_merqury_qv"
  printf "scaffolding_merqury_completeness\t%s\n" "$scaffolding_merqury_comp"
  printf "scaffolding_yak_qv\t%s\n" "$scaffolding_yak_qv"
  printf "scaffolding_scores\t%s\n" "$scaffolding_scores"
  printf "scaffolding_preselection_scores\t%s\n" "$scaffolding_preselection"
  printf "scaffolding_best_source\t%s\n" "$scaffolding_best_source"
  printf "quast_N50\t%s\n" "$N50"
  printf "quast_L50\t%s\n" "$L50"
  printf "quast_contigs\t%s\n" "$CONTIGS"
  printf "quast_GC_percent\t%s\n" "$GC"
  printf "quast_Ns_per_100kb\t%s\n" "$Ns100"
  printf "busco_lineage\t%s\n" "${BUSCO_LINEAGE:-}"
  printf "busco_C\t%s\n" "$BUSCO_C"
  printf "busco_S\t%s\n" "$BUSCO_S"
  printf "busco_D\t%s\n" "$BUSCO_D"
  printf "busco_F\t%s\n" "$BUSCO_F"
  printf "busco_M\t%s\n" "$BUSCO_M"
  printf "busco_n\t%s\n" "$BUSCO_n"
  printf "blobdir\t%s\n" "$blobdir_path"
  printf "blobplots\t%s\n" "${btk_dir}/plots"
  printf "repeats_status\t%s\n" "$repeat_summary_status"
  printf "repeats_dir\t%s\n" "$repeats_dir"
  printf "repeats_known_library\t%s\n" "$repeat_known_library"
  printf "repeats_repeatmodeler_library\t%s\n" "$repeat_rm_library"
  printf "repeats_tetrimmer_library\t%s\n" "$repeat_tetrimmer_library"
  printf "repeats_final_library\t%s\n" "$repeat_final_library"
  printf "repeats_repeatmasker_dir\t%s\n" "$repeat_repeatmasker_dir"
  printf "repeats_earlgrey_dir\t%s\n" "$repeat_earlgrey_dir"
  printf "repeats_final_outputs_dir\t%s\n" "$repeat_final_outputs_dir"
  printf "repeats_checkpoints_dir\t%s\n" "$repeat_checkpoints_dir"
  printf "repeats_logs_dir\t%s\n" "$repeat_logs_dir"
  printf "multiqc\t%s\n" "${dir7}/multiqc_report.html"
} > "$summary_tsv"

# --------------------------- write JSON (no jq) ---------------------------
{
  printf '{\n'
  printf '  "host": "%s",\n' "${host//\"/\\\"}"
  printf '  "kernel": "%s",\n' "${kernel//\"/\\\"}"
  printf '  "cpu_model": "%s",\n' "${cpu_model//\"/\\\"}"
  printf '  "cpu_threads": "%s",\n' "${cpu_threads}"
  printf '  "ram_total_gb": %.2f,\n' "${mem_total_gb:-0}"
  printf '  "start_time": "%s",\n' "${start_h//\"/\\\"}"
  printf '  "end_time": "%s",\n' "${end_h//\"/\\\"}"
  printf '  "duration_hms": "%s",\n' "${runtime_h//\"/\\\"}"
  printf '  "working_dir": "%s",\n' "${working_dir//\"/\\\"}"
  printf '  "working_dir_size": "%s",\n' "${run_du:-}"
  printf '  "final_asm": "%s",\n' "${final_target//\"/\\\"}"
  printf '  "final_contigs": %s,\n' "${final_ctg_n:-0}"
  printf '  "final_total_bp": %s,\n' "${final_ctg_bp:-0}"
  printf '  "scaffolding": {"best_id": "%s", "k": "%s", "w": "%s", "score": "%s", "N50": "%s", "L50": "%s", "contigs": "%s", "Ns_per_100kb": "%s", "merqury_qv": "%s", "merqury_completeness": "%s", "yak_qv": "%s", "scores": "%s", "preselection_scores": "%s", "best_source": "%s"},\n' \
         "${scaffolding_best_id//\"/\\\"}" \
         "${scaffolding_best_k//\"/\\\"}" \
         "${scaffolding_best_w//\"/\\\"}" \
         "${scaffolding_score//\"/\\\"}" \
         "${scaffolding_n50//\"/\\\"}" \
         "${scaffolding_l50//\"/\\\"}" \
         "${scaffolding_contigs//\"/\\\"}" \
         "${scaffolding_ns100//\"/\\\"}" \
         "${scaffolding_merqury_qv//\"/\\\"}" \
         "${scaffolding_merqury_comp//\"/\\\"}" \
         "${scaffolding_yak_qv//\"/\\\"}" \
         "${scaffolding_scores//\"/\\\"}" \
         "${scaffolding_preselection//\"/\\\"}" \
         "${scaffolding_best_source//\"/\\\"}"
  printf '  "quast": {"N50": "%s", "L50": "%s", "contigs": "%s", "GC_percent": "%s", "Ns_per_100kb": "%s"},\n' \
         "${N50:-}" "${L50:-}" "${CONTIGS:-}" "${GC:-}" "${Ns100:-}"
  printf '  "busco": {"lineage": "%s", "C": "%s", "S": "%s", "D": "%s", "F": "%s", "M": "%s", "n": "%s"},\n' \
         "${BUSCO_LINEAGE:-}" "${BUSCO_C:-}" "${BUSCO_S:-}" "${BUSCO_D:-}" "${BUSCO_F:-}" "${BUSCO_M:-}" "${BUSCO_n:-}"
  printf '  "blobtoolkit": {"blobdir": "%s", "plots_dir": "%s"},\n' \
         "${blobdir_path//\"/\\\"}" "${btk_dir}/plots"
  printf '  "repeats": {"status": "%s", "root": "%s", "known_library": "%s", "repeatmodeler_library": "%s", "tetrimmer_library": "%s", "final_library": "%s", "repeatmasker_dir": "%s", "earlgrey_dir": "%s", "final_outputs_dir": "%s", "checkpoints_dir": "%s", "logs_dir": "%s"},\n' \
         "${repeat_summary_status//\"/\\\"}" \
         "${repeats_dir//\"/\\\"}" \
         "${repeat_known_library//\"/\\\"}" \
         "${repeat_rm_library//\"/\\\"}" \
         "${repeat_tetrimmer_library//\"/\\\"}" \
         "${repeat_final_library//\"/\\\"}" \
         "${repeat_repeatmasker_dir//\"/\\\"}" \
         "${repeat_earlgrey_dir//\"/\\\"}" \
         "${repeat_final_outputs_dir//\"/\\\"}" \
         "${repeat_checkpoints_dir//\"/\\\"}" \
         "${repeat_logs_dir//\"/\\\"}"
  printf '  "multiqc": "%s"\n' "${dir7}/multiqc_report.html"
  printf '}\n'
} > "$summary_json"

# --------------------------- log tail aggregation ---------------------------
{
  # Collect recent tails from all logs under the pipeline tree
  find "$dir1" "$dir2" "$dir3" "$dir4" "$dir5" "$dir6" "$dir7" "${repeats_dir}" "$dir8" -type f -name "*.log" -print0 2>/dev/null \
    | sort -z \
    | while IFS= read -r -d '' f; do
        printf '\n### %s\n' "$(realpath --relative-to="$working_dir" "$f" 2>/dev/null || printf '%s' "$f")"
        tail -n 20 "$f" || true
      done
} > "${dir7}/logs/summary_tail.log" || true

# --------------------------- stderr topline ---------------------------
{
  printf "RUN SUMMARY :: Final=%s | Size=%s | N50=%s | L50=%s | BUSCO_C=%s%% | Repeats=%s | Host=%s | Threads=%s | RAM=%sGB | Start=%s | End=%s | Duration=%s\n" \
    "$(basename "$final_target")" "$final_ctg_size" "${N50:-'-'}" "${L50:-'-'}" "${BUSCO_C:-'-'}" \
    "$repeat_summary_status" "$host" "$cpu_threads" "$mem_total_gb" "$start_h" "$end_h" "$runtime_h"
} >&2

echo "[INFO] Pipeline complete. Summary: ${summary_txt} | ${summary_tsv} | ${summary_json}"
