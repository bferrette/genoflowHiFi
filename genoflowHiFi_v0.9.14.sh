#!/usr/bin/env bash
# shellcheck disable=SC2016 # literal inner-shell programs expand inside their target runtime

# Re-exec with bash if not already running under bash
if [ -z "${BASH_VERSION:-}" ]; then
  exec /usr/bin/env bash "$0" "$@"
fi

# Release metadata. Author identity follows the repository's published citation.
SCRIPT_NAME="genoflowHiFi"
SCRIPT_VERSION="0.9.14"
SCRIPT_RELEASE_DATE="2026-11-23"
SCRIPT_AUTHOR="B. L. S. Ferrette"
SCRIPT_CREATED="Not recorded in source history"
GENOFLOW_SCRIPT_DIR=${GENOFLOW_SCRIPT_DIR:-"$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"}

# Usage: bash genoflowHiFi_v0.9.14.sh 2>&1 | tee -a genoflowHiFi_v0.9.14.$(date +%Y%m%d_%H%M%S).log

# =============================== GLOBAL SHELL MODE ===============================
set -Eeuo pipefail
shopt -s extglob
umask 0022
export LC_ALL=C
ulimit -c 0   # production default (disable core dumps); set to 'unlimited' when debugging

# ==============================================================================================================
# HiFi-only workflow
# ==============================================================================================================
cat <<'EOF'
##############################################################################################################
#                                         🧬 genoflowHiFi v0.9.14 🧬                                         #
#                               PacBio HiFi reference-genome assembly                                        #
#                                      Release: 2026-11-23                                                   #
#                                                                                                            #
#  🎯 Purpose                                                                                                #
#    End-to-end reference-genome assembly from PacBio HiFi reads. The workflow performs haplotype-aware      #
#    primary assembly, optional decontamination, iterative consensus polishing, reference-free evaluation    #
#    and correction, HiFi-read scaffolding and selection, final refinement, validation, repeat analysis,     #
#    and consolidated reporting.                                                                             #
#                                                                                                            #
#  🧭 Main Workflow                                                                                          #
#                                                                                                            #
#    1) 🧩 Hifiasm                                                                                           #
#       Primary haplotype-aware contig assembly from PacBio HiFi reads.                                      #
#                                                                                                            #
#    2) 🧼 FCS-Decontam                                                                                      #
#       Optional NCBI FCS-Adaptor and FCS-GX screening, adaptor treatment, and contaminant removal.          #
#                                                                                                            #
#    3) 🛠️ DeepVariant Pre-Polish                                                                            #
#       Iterative HiFi-read consensus polishing with minimap2, SAMtools, DeepVariant, and BCFtools.          #
#                                                                                                            #
#    4) 🔍 Inspector                                                                                         #
#       Enhanced reference-free HiFi evaluation, correction, and independent replacement validation.         #
#                                                                                                            #
#    5) 🧱 LongStitch Scaffolding and Selection                                                              #
#       Optional HiFi-only Tigmint correction, Inspector-direct ntLink/ARKS scaffolding and mandatory        #
#       gap-filling, followed by genoflowSTATS metrics, yak, Merqury, and validated final candidate          #
#       selection.                                                                                           #
#                                                                                                            #
#    6) 💎 Final Refinement                                                                                  #
#       Post-scaffold DeepVariant polishing, final FCS cleanup, ≥1-kb scaffold filtering, and canonical      #
#       publication.                                                                                         #
#                                                                                                            #
#    7) 📊 Final Evaluation                                                                                  #
#       Native FASTA statistics, BUSCO, HiFi coverage/QC, DIAMOND, BLASTN, BlobToolKit, GenomeScope2,        #
#       genoflowSTATS longitudinal comparison, and MultiQC reporting.                                        #
#                                                                                                            #
#    8) 🧬 Repeat Analysis                                                                                   #
#       Dfam 4.0 and optional RepBase combine with RepeatModeler 2.0.9. Terrier 0.4 classifies unknown       #
#       families before TEtrimmer 1.7.4 and Pfam curation. Earl Grey 7.3.1 and RepeatMasker 4.2.4 perform    #
#       final annotation, enhanced plotting, soft masking, and non-repetitive interval derivation.           #
#                                                                                                            #
#    9) 📦 Final Summary                                                                                     #
#       Consolidated metrics, provenance, resource use, logs, checkpoints, and artifact manifests.           #
#                                                                                                            #
#  ✅ Design Principles                                                                                      #
#    ✅ Project inputs, resources, and step toggles are centralized in Section 1.                            #
#    ✅ Software paths, conda environments, containers, and tool policies are centralized in Section 2.      #
#    ✅ nthreads defines the per-process CPU budget; njobs controls supported parallel candidate work.       #
#    ✅ DeepVariant uses coordinate-sorted HiFi alignments without GATK duplicate marking.                   #
#    ✅ Checkpoints validate expected outputs and configuration fingerprints before reuse.                   #
#    ✅ Downstream stages consume only validated canonical assemblies.                                       #
#    ✅ Assembly taxonomy and repeat-library taxonomy are configured independently.                          #
#                                                                                                            #
#  🧰 Key Tool Groups                                                                                        #
#    🧩 Assembly / polishing     : Hifiasm, DeepVariant, Inspector                                           #
#    🧱 Scaffolding / selection  : LongStitch, GenomeScope2, genoflowSTATS, yak, Merqury                     #
#    📊 Assembly evaluation      : genoflowSTATS, BUSCO, DIAMOND, BLASTN, BlobToolKit, MultiQC               #
#    🧬 Repeat analysis          : RepeatModeler2, Terrier, TEtrimmer, Earl Grey, RepeatMasker               #
#    ⚙️  Infrastructure          : minimap2, SAMtools, BCFtools, GNU Parallel, Apptainer/Singularity         #
#    📚 Repeat libraries         : Dfam 4.0, optional RepBase, Pfam, and custom consensus                    #
#                                                                                                            #
#  📄 Primary Deliverable                                                                                    #
#    06.final_asm/<Scientific_Name>.fasta                                                                    #
#                                                                                                            #
#  📥 Required Inputs                                                                                        #
#    🧬 PacBio HiFi reads : required for assembly, polishing, validation, and repeat annotation              #
#    🧱 PacBio CLR reads  : optional; preferred for ntLink/ARKS scaffolding and gap filling when configured  #
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
pacbio_hifi=${pacbio_hifi:-"/data/bferrette/snook/rawdata/pacbio/snook_hifi.fq.gz"}
# Optional: pre-converted HiFi FASTA for scaffolding (skip if not set)
pacbio_hifi_fasta=${pacbio_hifi_fasta:-"/data/bferrette/snook/rawdata/pacbio/snook_hifi.fasta"}
# PacBio CLR reads for LongStitch scaffolding (skip if not set)
pacbio_clr=${pacbio_clr:-""}
# Hi-C raw data (FASTQ/FASTA or gz) for Hifiasm scaffolding (skip if not set)
hic_r1=${hic_r1:-""}  # optional; used only when HIFIASM_USE_HIC=TRUE
hic_r2=${hic_r2:-""}  # optional; used only when HIFIASM_USE_HIC=TRUE

# ---- Species / Taxonomy (edit per project) --------------------------------------
TAXID=${TAXID:-28805}            # NCBI Taxonomy ID for species
SPECIES_NAME=${SPECIES_NAME:-""} # explicit override, e.g. "Makaira nigricans"
SPECIES_HINT=${SPECIES_HINT:-""} # fallback label if taxdump unavailable

# ---- Resources (edit per project) -----------------------------------------------
nthreads=${nthreads:-32}            # threads per heavy task
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
# Completed-stage cleanup is limited to regenerable scratch. Canonical FASTAs,
# checkpoints, logs, VCFs, BAMs, repeat libraries, and provenance are retained.
CLEANUP_COMPLETED_STEPS=${CLEANUP_COMPLETED_STEPS:-"TRUE"}

# ===================================================================================================================
# 2) SOFTWARE / ENVIRONMENTS / CONTAINERS
# ===================================================================================================================

# ---- Conda bootstrap -------------------------------------------------------------
CONDA_SH=${CONDA_SH:-""}
if [[ -z "$CONDA_SH" || ! -s "$CONDA_SH" ]]; then
  candidate="$HOME/miniforge3/etc/profile.d/conda.sh"
  [[ -s "$candidate" ]] && CONDA_SH="$candidate"
fi

if [[ -s "$CONDA_SH" ]]; then
  # shellcheck source=/dev/null
  source "$CONDA_SH"
else
  echo "[ERROR] Cannot source conda from $HOME miniforge3 paths" >&2
  exit 1
fi

# Root of the Miniforge installation. Derive it from the selected conda.sh
# so `set -u` never depends on an undeclared MINIFORGE_PREFIX.
if [[ -z "${MINIFORGE_PREFIX:-}" ]]; then
  MINIFORGE_PREFIX="$(
    cd "$(dirname "${CONDA_SH}")/../.." && pwd -P
  )" || {
    echo "[ERROR] Could not resolve Miniforge prefix from ${CONDA_SH}" >&2
    exit 1
  }
fi

[[ -d "${MINIFORGE_PREFIX}" ]] || {
  echo "[ERROR] Miniforge prefix does not exist: ${MINIFORGE_PREFIX}" >&2
  exit 1
}

MAMBA_SH=${MAMBA_SH:-""}
if [[ -z "$MAMBA_SH" || ! -s "$MAMBA_SH" ]]; then
  candidate="$HOME/miniforge3/etc/profile.d/mamba.sh"
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
ENV_GENOFLOWSTATS=${ENV_GENOFLOWSTATS:-"GENOFLOWSTATSenv"}
ENV_LONGSTITCH=${ENV_LONGSTITCH:-"LONGSTITCHenv"}
ENV_MERQURY=${ENV_MERQURY:-"MERQURYenv"}
ENV_GENOMESCOPE2=${ENV_GENOMESCOPE2:-"GENOMESCOPE2env"}
GENOMESCOPE2_ENV_PREFIX=${GENOMESCOPE2_ENV_PREFIX:-"${MINIFORGE_PREFIX}/envs/${ENV_GENOMESCOPE2}"}
ENV_R=${ENV_R:-"Renv"}
ENV_BUSCO=${ENV_BUSCO:-"BUSCOenv"}
ENV_DIAMOND=${ENV_DIAMOND:-"DIAMONDenv"}
ENV_TERRIER=${ENV_TERRIER:-"TERRIERenv"}
ENV_TETRIMMER=${ENV_TETRIMMER:-"TETRIMMERenv"}
TETRIMMER_SCRIPT=${TETRIMMER_SCRIPT:-""}  # resolved lazily by resolve_tetrimmer_script()
ENV_EARLGREY=${ENV_EARLGREY:-"EARLGREYenv"}

# ---- Tool binaries (compiled installs) ------------------------------------------
bcftools_bin=${bcftools_bin:-"$HOME/software/bcftools/bcftools"}
bedtools_bin=${bedtools_bin:-"$HOME/software/bedtools2/bin/bedtools"}
blastn_bin=${blastn_bin:-"$HOME/software/ncbi-blast-2.17.0+/bin/blastn"}
minimap2_bin=${minimap2_bin:-"$HOME/software/minimap2/minimap2"}
python_bin=${python_bin:-"$HOME/miniforge3/bin/python3"}
samtools_bin=${samtools_bin:-"$HOME/software/samtools/samtools"}
seqkit_bin=${seqkit_bin:-"$HOME/software/seqkit-2.13.0/seqkit"}

# ---- Containers (Singularity/Apptainer) -----------------------------------------
sif_blobtoolkit=${sif_blobtoolkit:-"$HOME/containers/blobtoolkit.sif"}
sif_deepvariant=${sif_deepvariant:-"$HOME/containers/deepvariant.sif"}
sif_fcs_adaptor=${sif_fcs_adaptor:-${ADAPTOR_SIF:-"$HOME/containers/fcs-adaptor.sif"}}
sif_fcs_gx=${sif_fcs_gx:-${GX_SIF:-"$HOME/containers/fcs-gx.sif"}}
sif_repeatmasker=${sif_repeatmasker:-${REPEATMASKER_SIF:-"$HOME/containers/repeatmasker_v4.sif"}}

# ---- NCBI FCS (central knobs) ----------------------------------------------------
ADAPTOR_RUNNER=${ADAPTOR_RUNNER:-"$HOME/software/fcs/dist/run_fcsadaptor.sh"}
GX_ROOT=${GX_ROOT:-"/storage/bferrette/ncbi/refseq/fcs-gx"}
GX_DB=${GX_DB:-"$GX_ROOT/gxdb"}
GX_MANIFEST=${GX_MANIFEST:-"https://ncbi-fcs-gx.s3.amazonaws.com/gxdb/latest/all.manifest"}
FCS_TINI_SUBREAPER=${FCS_TINI_SUBREAPER:-"TRUE"}
FCS_ADAPTOR_MIN_SEQ_LEN=${FCS_ADAPTOR_MIN_SEQ_LEN:-11}
FCS_ADAPTOR_MAX_SEQ_LEN=${FCS_ADAPTOR_MAX_SEQ_LEN:-2147483647}

# ---- Steps 1-4 enhanced assembly controls ---------------------------------------
# Hifiasm: unset expert parameters deliberately defer to upstream defaults.
HIFIASM_USE_HIC=${HIFIASM_USE_HIC:-FALSE}
HIFIASM_PRIMARY_ALT=${HIFIASM_PRIMARY_ALT:-"FALSE"}
HIFIASM_COMPRESS_GFA=${HIFIASM_COMPRESS_GFA:-"TRUE"}
# Store Hifiasm resume caches as .bin.zst (zstd -9, multithreaded).
# If zstd is unavailable, use pigz --best -p and .bin.gz.
HIFIASM_COMPRESS_BIN=${HIFIASM_COMPRESS_BIN:-"TRUE"}
HIFIASM_HG_SIZE=${HIFIASM_HG_SIZE:-""}
HIFIASM_HOM_COV=${HIFIASM_HOM_COV:-""}
HIFIASM_PURGE_LEVEL=${HIFIASM_PURGE_LEVEL:-""}
HIFIASM_PURGE_SIMILARITY=${HIFIASM_PURGE_SIMILARITY:-""}
HIFIASM_POST_JOIN=${HIFIASM_POST_JOIN:-"1"}
HIFIASM_CTG_N=${HIFIASM_CTG_N:-""}
HIFIASM_ADAPTER_TRIM=${HIFIASM_ADAPTER_TRIM:-""}
HIFIASM_BLOOM_BITS=${HIFIASM_BLOOM_BITS:-""}
HIFIASM_CORRECTION_ROUNDS=${HIFIASM_CORRECTION_ROUNDS:-""}
HIFIASM_MAX_KOCC=${HIFIASM_MAX_KOCC:-""}
HIFIASM_HIGH_KMER_FACTOR=${HIFIASM_HIGH_KMER_FACTOR:-""}
HIFIASM_MAX_OVERLAPS=${HIFIASM_MAX_OVERLAPS:-""}
PACBIO_HIFI_READS=${PACBIO_HIFI_READS:-"${pacbio_hifi}"}

# FCS: assembly-organism taxonomy is deliberately independent of repeat taxonomy.
FCS_RUNNER=${FCS_RUNNER:-"$(dirname "${ADAPTOR_RUNNER}")/fcs.py"}
FCS_GX_TAXID=${FCS_GX_TAXID:-"${TAXID}"}
FCS_GX_DB=${FCS_GX_DB:-"${GX_DB}"}
FCS_GX_MANIFEST=${FCS_GX_MANIFEST:-"${GX_MANIFEST}"}
# Database lifecycle policy:
#   check          validate the installed DB, compare its manifest with NCBI latest, report updates, continue locally
#   update         synchronize a missing/outdated DB to the configured manifest, then validate it
#   require_latest validate locally and abort when the configured remote manifest is newer/different
#   offline        perform local integrity validation only; never access the network
#   download       bootstrap a missing DB and otherwise synchronize it (equivalent to update once present)
#   validate       backward-compatible alias for check
FCS_GX_DB_POLICY=${FCS_GX_DB_POLICY:-"check"}
# Full FCS-GX database integrity scan. FALSE skips only `fcs.py db check`;
# required-file checks, manifest freshness checks, synchronization, and screening remain active.
FCS_GX_DB_CHECK=${FCS_GX_DB_CHECK:-"FALSE"}
FCS_GX_CORES=${FCS_GX_CORES:-"${nthreads}"}
FCS_GX_ALIGN_EXCLUDE_TAXA=${FCS_GX_ALIGN_EXCLUDE_TAXA:-""}
FCS_GX_EXTRA_CONTAM_DIVS=${FCS_GX_EXTRA_CONTAM_DIVS:-""}
FCS_ADAPTOR_TAXONOMY=${FCS_ADAPTOR_TAXONOMY:-"euk"}  # euk | prok
FCS_ADAPTOR_INTERNAL_ACTION=${FCS_ADAPTOR_INTERNAL_ACTION:-"split"}  # split | mask
FCS_GX_INTERNAL_ACTION=${FCS_GX_INTERNAL_ACTION:-"fix"}  # fix | split
FCS_DEBUG=${FCS_DEBUG:-"FALSE"}

# Legacy lowercase round names are accepted only as input aliases here.
# DeepVariant uses nthreads as its per-process CPU budget.
DEEPVARIANT_PRE_ROUNDS=${DEEPVARIANT_PRE_ROUNDS:-${deepvariant_pre_rounds:-3}}
DEEPVARIANT_SHARDS=${DEEPVARIANT_SHARDS:-"${nthreads}"}
DEEPVARIANT_POSTPROCESS_CPUS=${DEEPVARIANT_POSTPROCESS_CPUS:-"${nthreads}"}
DEEPVARIANT_VCF_STATS=${DEEPVARIANT_VCF_STATS:-"true"}
DEEPVARIANT_STOP_ON_ZERO=${DEEPVARIANT_STOP_ON_ZERO:-"TRUE"}
DEEPVARIANT_WRITE_CHAIN=${DEEPVARIANT_WRITE_CHAIN:-"TRUE"}
DEEPVARIANT_CLEAN_INTERMEDIATES=${DEEPVARIANT_CLEAN_INTERMEDIATES:-"TRUE"}

# DeepVariant/HiFi BAM compatibility.
#
# DeepVariant r1.10 reads QUAL into Read.aligned_quality and AlleleCounter
# assumes that quality values exist for every aligned query base. Records with
# QUAL=* can therefore abort make_examples with:
#   Check failed: offset + len <= read.aligned_quality().size()
#
# minimap2 --secondary=no suppresses secondary mappings, but supplementary
# alignments are separate SAM records. DeepVariant's read-requirement filter is
# applied before QUAL decoding, so keep_supplementary_alignments=false is valid
# defense in depth. The BAM itself is nevertheless made stricter: exclude
# secondary/supplementary records and require a real QUAL array before sorting.
readonly DEEPVARIANT_EXCLUDE_FLAGS="0x900"       # 0x100 secondary + 0x800 supplementary
readonly DEEPVARIANT_VIEW_EXPR="length(seq) > 0 && length(qual) == length(seq)"   # require non-empty SEQ and one QUAL value per base
readonly DEEPVARIANT_MAKE_EXAMPLES_EXTRA_ARGS="keep_supplementary_alignments=false"

# minimap2 -K controls the query mini-batch size. MINIMAP2_INDEX_BATCH is kept
# only as a backward-compatible input alias from older pipeline versions.
MINIMAP2_INDEX_BATCH=${MINIMAP2_INDEX_BATCH:-"16G"}
MINIMAP2_BATCH_SIZE=${MINIMAP2_BATCH_SIZE:-"${MINIMAP2_INDEX_BATCH}"}
SAMTOOLS_SORT_MEM=${SAMTOOLS_SORT_MEM:-"8G"}
SAMTOOLS_BAM_COMPRESSION=${SAMTOOLS_BAM_COMPRESSION:-9}

# Shared post-scaffolding refinement (DeepVariant + final FCS).
DEEPVARIANT_POST_ROUNDS=${DEEPVARIANT_POST_ROUNDS:-${deepvariant_post_rounds:-3}}
DEEPVARIANT_FINAL_POSTPROCESS_CPUS=${DEEPVARIANT_FINAL_POSTPROCESS_CPUS:-${DEEPVARIANT_POSTPROCESS_CPUS}}
POST_POLISH_HASH_READS=${POST_POLISH_HASH_READS:-"FALSE"}
POST_POLISH_KEEP_VCF=${POST_POLISH_KEEP_VCF:-"TRUE"}
RUN_FINAL_FCS=${RUN_FINAL_FCS:-"TRUE"}
FCS_FINAL_ADAPTOR_INTERNAL_ACTION=${FCS_FINAL_ADAPTOR_INTERNAL_ACTION:-"split"}  # mask | split
FCS_FINAL_GX_INTERNAL_ACTION=${FCS_FINAL_GX_INTERNAL_ACTION:-"fix"}            # fix | split
FINAL_FASTA_ID_PREFIX=${FINAL_FASTA_ID_PREFIX:-"scaff"}

# Inspector: enhanced Python 3 replacements are expected inside ENV_INSPECTOR.
# Empty length thresholds intentionally defer to Inspector's documented defaults (10 kb / 1 Mb).
INSPECTOR_THREADS=${INSPECTOR_THREADS:-"${nthreads}"}
INSPECTOR_MIN_CONTIG_LENGTH=${INSPECTOR_MIN_CONTIG_LENGTH:-""}
INSPECTOR_MIN_STRUCTURAL_CONTIG=${INSPECTOR_MIN_STRUCTURAL_CONTIG:-""}
INSPECTOR_MINIMAP2_PRESET=${INSPECTOR_MINIMAP2_PRESET:-"auto"}
INSPECTOR_MINIMAP2_HPC=${INSPECTOR_MINIMAP2_HPC:-"auto"}
INSPECTOR_MINIMAP2_INDEX_BATCH=${INSPECTOR_MINIMAP2_INDEX_BATCH:-"16G"}
INSPECTOR_SAMTOOLS_SORT_MEMORY=${INSPECTOR_SAMTOOLS_SORT_MEMORY:-"8G"}
INSPECTOR_BAM_COMPRESSION=${INSPECTOR_BAM_COMPRESSION:-9}
INSPECTOR_EVALUATE_STRUCTURAL=${INSPECTOR_EVALUATE_STRUCTURAL:-"TRUE"}
INSPECTOR_EVALUATE_BASE=${INSPECTOR_EVALUATE_BASE:-"TRUE"}
INSPECTOR_CORRECT=${INSPECTOR_CORRECT:-"TRUE"}
INSPECTOR_CORRECT_STRUCTURAL=${INSPECTOR_CORRECT_STRUCTURAL:-"TRUE"}
INSPECTOR_CORRECT_BASE=${INSPECTOR_CORRECT_BASE:-"TRUE"}
INSPECTOR_CORRECTION_JOBS=${INSPECTOR_CORRECTION_JOBS:-1}
INSPECTOR_FLYE_TIMEOUT=${INSPECTOR_FLYE_TIMEOUT:-1200}
INSPECTOR_KEEP_INTERMEDIATES=${INSPECTOR_KEEP_INTERMEDIATES:-"FALSE"}
INSPECTOR_VALIDATE_CORRECTION=${INSPECTOR_VALIDATE_CORRECTION:-"TRUE"}

# ---- External DBs / References ---------------------------------------------------
BLAST_DB_ROOT=${BLAST_DB_ROOT:-"/data/bferrette/ncbi/refseq/blastdb/core_nt/core_nt"}
BLASTX_DB_ROOT=${BLASTX_DB_ROOT:-"/data/bferrette/databases/diamond/refseq_protein.dmnd"}
BUSCO_LINEAGE=${BUSCO_LINEAGE:-"actinopterygii_odb12.2"}

# ---- Shared final assembly publication / statistics -----------------------------
# Final FASTA cleanup is intentionally conservative: retain only scaffolds >=1 kb,
# characterize the pre/post-filter length distribution with SeqKit + Python, then
# publish simple sequential scaffold IDs (scaff1, scaff2, ...).
FINAL_FASTA_MIN_LENGTH=${FINAL_FASTA_MIN_LENGTH:-1000}
FINAL_FASTA_QC_DPI=${FINAL_FASTA_QC_DPI:-600}
FINAL_FASTA_QC_MAX_CURVE_POINTS=${FINAL_FASTA_QC_MAX_CURVE_POINTS:-100000}
FINAL_FASTA_QC_SAVE_SVG=${FINAL_FASTA_QC_SAVE_SVG:-"TRUE"}

# ---- Shared final assembly statistics / validation ------------------------------
BUSCO_DATA_POLICY=${BUSCO_DATA_POLICY:-"check"}        # check | update | offline
# Online preparation is bounded; scoring uses the validated local lineage offline.
BUSCO_NETWORK_TIMEOUT=${BUSCO_NETWORK_TIMEOUT:-60}   # seconds per socket operation
BUSCO_DOWNLOAD_TIMEOUT=${BUSCO_DOWNLOAD_TIMEOUT:-3600} # seconds per download attempt
BUSCO_DOWNLOAD_ATTEMPTS=${BUSCO_DOWNLOAD_ATTEMPTS:-3}
BUSCO_DOWNLOAD_RETRY_DELAY=${BUSCO_DOWNLOAD_RETRY_DELAY:-10}
BUSCO_DOWNLOAD_BASE_URL=${BUSCO_DOWNLOAD_BASE_URL:-"https://busco-data.ezlab.org/v5/data/"}
BUSCO_DOWNLOAD_FALLBACK_URL=${BUSCO_DOWNLOAD_FALLBACK_URL:-"https://busco-data2.ezlab.org/v5/data/"}
RUN_DIAMOND_STATS=${RUN_DIAMOND_STATS:-"TRUE"}
RUN_BLAST_STATS=${RUN_BLAST_STATS:-"TRUE"}             # BLASTN fallback on DIAMOND no-hit chunks
DIAMOND_SENSITIVITY=${DIAMOND_SENSITIVITY:-"sensitive"} # fast | sensitive | more-sensitive | very-sensitive | ultra-sensitive
DIAMOND_EVALUE=${DIAMOND_EVALUE:-"1e-10"}
DIAMOND_MAX_TARGET_SEQS=${DIAMOND_MAX_TARGET_SEQS:-10}
DIAMOND_BLOCK_SIZE=${DIAMOND_BLOCK_SIZE:-2.0}
DIAMOND_INDEX_CHUNKS=${DIAMOND_INDEX_CHUNKS:-4}
BTK_TAXONOMY_CHUNK_BP=${BTK_TAXONOMY_CHUNK_BP:-100000}
BLAST_TASK=${BLAST_TASK:-"megablast"}                  # megablast | dc-megablast
BLAST_MAX_TARGET_SEQS=${BLAST_MAX_TARGET_SEQS:-10}
BLAST_MAX_HSPS=${BLAST_MAX_HSPS:-1}
BLAST_EVALUE=${BLAST_EVALUE:-"1e-25"}
BTK_ENABLE_TAXONOMY=${BTK_ENABLE_TAXONOMY:-"TRUE"}
BTK_HIT_COUNT=${BTK_HIT_COUNT:-10}
STATS_SAMTOOLS_SORT_MEM=${STATS_SAMTOOLS_SORT_MEM:-"8G"}
STATS_BAM_COMPRESSION=${STATS_BAM_COMPRESSION:-9}
STATS_GENOMESCOPE_K=${STATS_GENOMESCOPE_K:-31}
STATS_GENOMESCOPE_PLOIDY=${STATS_GENOMESCOPE_PLOIDY:-2}
STATS_GENOMESCOPE_TRANSFORM_EXP=${STATS_GENOMESCOPE_TRANSFORM_EXP:-1}

# ---- Step 5 enhanced scaffolding and assembly selection --------------------------
GENOME_SIZE=${GENOME_SIZE:-""}              # optional; avoid hard-coded false precision
GENOME_SIZE_EST=${GENOME_SIZE_EST:-""}

LONGSTITCH_K_VALUES=${LONGSTITCH_K_VALUES:-"21 24 27 32 35 40 45 51"}
LONGSTITCH_W_VALUES=${LONGSTITCH_W_VALUES:-"50 100 150 250 350 500 750 1000"}
# Tigmint correction is optional. By default, Step 5 skips Tigmint and feeds the
# Inspector-corrected FASTA directly to the ntLink grid. When enabled, Tigmint is
# run once with HiFi reads because its correction is invariant across the k x w grid.
LONGSTITCH_USE_TIGMINT=${LONGSTITCH_USE_TIGMINT:-"FALSE"}
LONGSTITCH_TIGMINT_CUT=${LONGSTITCH_TIGMINT_CUT:-250}
LONGSTITCH_TIGMINT_SPAN=${LONGSTITCH_TIGMINT_SPAN:-"auto"}
LONGSTITCH_TIGMINT_DIST=${LONGSTITCH_TIGMINT_DIST:-"auto"}
NTLINK_Z=${NTLINK_Z:-1000}
NTLINK_N=${NTLINK_N:-2}
NTLINK_A=${NTLINK_A:-1}
NTLINK_OVERLAP=${NTLINK_OVERLAP:-"True"}
NTLINK_SENSITIVE=${NTLINK_SENSITIVE:-"False"}

# Post-grid refinement. This HiFi-only/no-Hi-C pipeline always performs a final
# ntLink gap-fill pass. Optional refinement before gap filling is controlled by
# SCAFFOLD_REFINEMENT. Allowed modes all include mandatory gap filling:
#   gapfill | rounds-gapfill | arks-gapfill | rounds-arks-gapfill
SCAFFOLD_REFINEMENT=${SCAFFOLD_REFINEMENT:-"arks-gapfill"}
SCAFFOLD_REFINEMENT_TOP_N=${SCAFFOLD_REFINEMENT_TOP_N:-5}
NTLINK_REFINEMENT_ROUNDS=${NTLINK_REFINEMENT_ROUNDS:-3}
ARKS_M=${ARKS_M:-"8-10000"}
ARKS_J=${ARKS_J:-0.05}
ARKS_K=${ARKS_K:-20}
ARKS_C=${ARKS_C:-4}
ARKS_L=${ARKS_L:-4}
ARKS_A=${ARKS_A:-0.3}
ARKS_CUT=${ARKS_CUT:-250}
ARKS_R=${ARKS_R:-0.05}
ARKS_E=${ARKS_E:-30000}
ARKS_D=${ARKS_D:-true}


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
GENOMESCOPE2_BIN=${GENOMESCOPE2_BIN:-"genomescope.R"}

MERQURY_K=${MERQURY_K:-31}
MERQURY_TOP_N=${MERQURY_TOP_N:-10}
MERQURY=${MERQURY:-""}                       # optional path to Merqury repo if merqury.sh is not on PATH

YAK_K=${YAK_K:-31}
YAK_BLOOM_BITS=${YAK_BLOOM_BITS:-37}
STEP5_SAVE_SVG=${STEP5_SAVE_SVG:-"TRUE"}
STEP5_PLOT_DPI=${STEP5_PLOT_DPI:-600}

# ---- genoflowSTATS assembly assessment ------------------------------------
GENOFLOWSTATS_GAP_MIN=${GENOFLOWSTATS_GAP_MIN:-10}
GENOFLOWSTATS_PLOT_DPI=${GENOFLOWSTATS_PLOT_DPI:-600}
GENOFLOWSTATS_PALETTE=${GENOFLOWSTATS_PALETTE:-"colorblind"}
GENOFLOWSTATS_COMPARE_ALL_STAGES=${GENOFLOWSTATS_COMPARE_ALL_STAGES:-"TRUE"}
GENOFLOWSTATS_COMPARE_READ_SUPPORT=${GENOFLOWSTATS_COMPARE_READ_SUPPORT:-"TRUE"}
GENOFLOWSTATS_KEEP_STAGE_BAMS=${GENOFLOWSTATS_KEEP_STAGE_BAMS:-"FALSE"}
GENOFLOWSTATS_COV_MIN=${GENOFLOWSTATS_COV_MIN:-1}
GENOFLOWSTATS_COV_MAX=${GENOFLOWSTATS_COV_MAX:-1000}
GENOFLOWSTATS_COV_STEP=${GENOFLOWSTATS_COV_STEP:-1}
GENOFLOWSTATS_COMPARE_THREADS=${GENOFLOWSTATS_COMPARE_THREADS:-"${nthreads}"}

# ---- Repeat annotation knobs (edit per project) ---------------------------------
# Empty REPEATS_* path/name overrides are resolved at Step 8 after final assembly paths exist.
REPEATS_PROJECT_ROOT=${REPEATS_PROJECT_ROOT:-""}   # default: ${working_dir}/08.repeats
REPEATS_GENOME=${REPEATS_GENOME:-""}               # default: final species FASTA from Step 6
REPEAT_TAXON=${REPEAT_TAXON:-"7777"}  # authoritative repeat-annotation taxon; may be a broader NCBI Taxonomy ID such as 7777 or 1489904
REPEAT_TAXON_USE_GENOME_SPECIES=${REPEAT_TAXON_USE_GENOME_SPECIES:-"TRUE"}   # fallback only; explicit REPEAT_TAXON takes precedence
REPEATS_SPECIES=${REPEATS_SPECIES:-""}  # repeat-annotation target; may differ from the assembled organism
REPEATS_THREADS=${REPEATS_THREADS:-${nthreads}} # CPU budget for one repeat-annotation process

AUTO_RESUME=${AUTO_RESUME:-"TRUE"}
AUTO_ARCHIVE_UNRECOVERABLE_RM_DIRS=${AUTO_ARCHIVE_UNRECOVERABLE_RM_DIRS:-"TRUE"}

KNOWN_LIBRARY_OVERRIDE=${KNOWN_LIBRARY_OVERRIDE:-""}
REPEATMODELER_LIBRARY_OVERRIDE=${REPEATMODELER_LIBRARY_OVERRIDE:-""}
TERRIER_LIBRARY_OVERRIDE=${TERRIER_LIBRARY_OVERRIDE:-""}
TETRIMMER_LIBRARY_OVERRIDE=${TETRIMMER_LIBRARY_OVERRIDE:-""}
FINAL_LIBRARY_OVERRIDE=${FINAL_LIBRARY_OVERRIDE:-""}

REPEATS_LIB_DIR=${REPEATS_LIB_DIR:-"$HOME/miniforge3/envs/EARLGREYenv/share/RepeatMasker/Libraries"}              # default: ${REPEATS_PROJECT_ROOT}/Libraries
# Host directory containing external Dfam FamDB HDF5 partitions.
# It is mounted inside the container at /opt/RepeatMasker/Libraries/famdb.
DFAM_FAMDB_DIR=${DFAM_FAMDB_DIR:-""}               # default: ${REPEATS_LIB_DIR}/famdb
FAMDB_CONTAINER_DIR=${FAMDB_CONTAINER_DIR:-"/opt/FamDB"}
FAMDB_CONTAINER_DATA_DIR=${FAMDB_CONTAINER_DATA_DIR:-"/opt/RepeatMasker/Libraries/famdb"}
REPEATMODELER_CONTAINER_DIR=${REPEATMODELER_CONTAINER_DIR:-"/opt/RepeatModeler"}
DFAM_TAXON=${DFAM_TAXON:-""}  # Dfam/FamDB query taxon; may differ from genome TAXID/SPECIES_NAME
DFAM_FLAGS=${DFAM_FLAGS:-"--curated"}               # curated Dfam families by default; may be overridden explicitly
DFAM_AUTODECOMPRESS_GZ=${DFAM_AUTODECOMPRESS_GZ:-"TRUE"}
DFAM_VALIDATE_MD5=${DFAM_VALIDATE_MD5:-"TRUE"}
DFAM_RELEASE_TAG=${DFAM_RELEASE_TAG:-"dfam40"}
DFAM_REQUIRED_PARTITIONS=${DFAM_REQUIRED_PARTITIONS:-"dfam40.0 dfam40.curated.consensus.0 dfam40.curated.hmm.0 dfam40.curated.hmm.1 dfam40.uncurated.consensus.0 dfam40.uncurated.consensus.1 dfam40.uncurated.hmm.0 dfam40.uncurated.hmm.83 dfam40.uncurated.hmm.84 dfam40.uncurated.hmm.85 dfam40.uncurated.hmm.91 dfam40.uncurated.hmm.92 dfam40.uncurated.hmm.93 dfam40.uncurated.hmm.94 dfam40.uncurated.hmm.95 dfam40.uncurated.hmm.96"}
REPEAT_USE_REPBASE=${REPEAT_USE_REPBASE:-"TRUE"}     # optional complementary known-repeat source
REPBASE_DIR=${REPBASE_DIR:-""}                     # default: ${REPEATS_LIB_DIR}
REPBASE_EMBL=${REPBASE_EMBL:-""}                   # optional explicit RepBase EMBL/archive path
REPBASE_EMBL_PATTERN=${REPBASE_EMBL_PATTERN:-"RepBase*.embl.tar.gz"}
ARTEFACTS_FASTA=${ARTEFACTS_FASTA:-""}

USE_LTRSTRUCT=${USE_LTRSTRUCT:-"yes"}
REPEATMODELER_RECOVER=${REPEATMODELER_RECOVER:-"TRUE"}
REPEATMODELER_RECOVER_DIR=${REPEATMODELER_RECOVER_DIR:-""}
REPEATMODELER_SANITIZE_FASTA_HEADERS=${REPEATMODELER_SANITIZE_FASTA_HEADERS:-"TRUE"}
RM_QUICK=${RM_QUICK:-"FALSE"}
RM_SRAND=${RM_SRAND:-""}
RM_GENOME_SAMPLE_SIZE_MAX=${RM_GENOME_SAMPLE_SIZE_MAX:-""}
RM_ADDL_ROUNDS=${RM_ADDL_ROUNDS:-0}

# Terrier classifies RepeatModeler Unknown families before TEtrimmer curation.
TERRIER_ENABLE=${TERRIER_ENABLE:-"TRUE"}
TERRIER_REQUIRED_VERSION=${TERRIER_REQUIRED_VERSION:-"0.4"}
TERRIER_THRESHOLD=${TERRIER_THRESHOLD:-"0.90"}
TERRIER_MIN_LENGTH=${TERRIER_MIN_LENGTH:-128}
TERRIER_MAX_LENGTH=${TERRIER_MAX_LENGTH:-5000}
TERRIER_BATCH_SIZE=${TERRIER_BATCH_SIZE:-1}
TERRIER_UNKNOWN_ONLY=${TERRIER_UNKNOWN_ONLY:-"TRUE"}
TERRIER_IMAGE_FORMAT=${TERRIER_IMAGE_FORMAT:-"png"}
TERRIER_IMAGE_THRESHOLD=${TERRIER_IMAGE_THRESHOLD:-"0.005"}
TERRIER_THREADS=${TERRIER_THREADS:-${REPEATS_THREADS}}

TETRIMMER_REQUIRED_VERSION=${TETRIMMER_REQUIRED_VERSION:-"1.7.4"}
TETRIMMER_MODE=${TETRIMMER_MODE:-"curated"}              # fast | curated | debug
TETRIMMER_PRESET=${TETRIMMER_PRESET:-"divergent"}        # conserved | divergent | empty
TETRIMMER_THREAD_CAP=${TETRIMMER_THREAD_CAP:-16}         # safety cap; hard-limited to 16 to avoid RAM overflow
TETRIMMER_CONTINUE=${TETRIMMER_CONTINUE:-"TRUE"}
TETRIMMER_CONTINUE_MODE=${TETRIMMER_CONTINUE_MODE:-"auto"} # auto | always | never
TETRIMMER_ARCHIVE_INCOMPLETE=${TETRIMMER_ARCHIVE_INCOMPLETE:-"TRUE"}
TETRIMMER_CLASSIFICATION=${TETRIMMER_CLASSIFICATION:-"none"} # none | classify_all; legacy classify_unknown also disables classification
TETRIMMER_DEDUP=${TETRIMMER_DEDUP:-"FALSE"}
TETRIMMER_HMM=${TETRIMMER_HMM:-"FALSE"}
TETRIMMER_GENOME_ANNO=${TETRIMMER_GENOME_ANNO:-"FALSE"}
TETRIMMER_DEBUG=${TETRIMMER_DEBUG:-"FALSE"}
TETRIMMER_COMPRESS_OUTPUT=${TETRIMMER_COMPRESS_OUTPUT:-"FALSE"}
TETRIMMER_EXPORT_COVERAGE=${TETRIMMER_EXPORT_COVERAGE:-"FALSE"}
TETRIMMER_MAX_THREAD_TIME=${TETRIMMER_MAX_THREAD_TIME:-3600}
TETRIMMER_SKIP_REST_SEQ=${TETRIMMER_SKIP_REST_SEQ:-"FALSE"} # rescue only; requires continuation
TETRIMMER_PFAM_DIR=${TETRIMMER_PFAM_DIR:-""}        # resolved later to ${PROJECT_ROOT}/resources/pfam; kept outside TEtrimmer output
# Keep FALSE to let --preset supply its internally consistent parameter set.
# Expert values below are passed only when TETRIMMER_EXPERT_OVERRIDES=TRUE.
TETRIMMER_EXPERT_OVERRIDES=${TETRIMMER_EXPERT_OVERRIDES:-"FALSE"}
TETRIMMER_CONS_THR=${TETRIMMER_CONS_THR:-"0.7"}
TETRIMMER_MINI_ORF=${TETRIMMER_MINI_ORF:-"200"}
TETRIMMER_MAX_MSA_LINES=${TETRIMMER_MAX_MSA_LINES:-"100"}
TETRIMMER_TOP_MSA_LINES=${TETRIMMER_TOP_MSA_LINES:-"100"}
TETRIMMER_MIN_SEQ_NUM=${TETRIMMER_MIN_SEQ_NUM:-"10"}
TETRIMMER_MIN_BLAST_LEN=${TETRIMMER_MIN_BLAST_LEN:-"150"}
TETRIMMER_MAX_CLUSTER_NUM=${TETRIMMER_MAX_CLUSTER_NUM:-"5"}
TETRIMMER_EXT_THR=${TETRIMMER_EXT_THR:-"0.6"}
TETRIMMER_EXT_CHECK_WIN=${TETRIMMER_EXT_CHECK_WIN:-"150"}
TETRIMMER_EXT_STEP=${TETRIMMER_EXT_STEP:-"1500"}
TETRIMMER_MAX_EXT=${TETRIMMER_MAX_EXT:-"8000"}
TETRIMMER_GAP_THR=${TETRIMMER_GAP_THR:-"0.4"}
TETRIMMER_GAP_NUL_THR=${TETRIMMER_GAP_NUL_THR:-"0.65"}
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

FINAL_CDHIT_ID=${FINAL_CDHIT_ID:-"0.95"}
FINAL_CDHIT_WORD=${FINAL_CDHIT_WORD:-"9"}
FINAL_CDHIT_AS=${FINAL_CDHIT_AS:-"0.90"}
FINAL_CDHIT_AL=${FINAL_CDHIT_AL:-"0.80"}
FINAL_CDHIT_MEMORY_MB=${FINAL_CDHIT_MEMORY_MB:-16000}

EARLGREY_REQUIRED_VERSION=${EARLGREY_REQUIRED_VERSION:-"7.3.1"}
EARLGREY_REQUIRED_REPEATMASKER=${EARLGREY_REQUIRED_REPEATMASKER:-"4.2.4"}
EARLGREY_REQUIRED_REPEATMODELER=${EARLGREY_REQUIRED_REPEATMODELER:-"2.0.9"}
EARLGREY_REQUIRED_DFAM=${EARLGREY_REQUIRED_DFAM:-"4.0"}
EARLGREY_REMOVE_SHORT=${EARLGREY_REMOVE_SHORT:-"no"}
EARLGREY_SOFTMASK=${EARLGREY_SOFTMASK:-"yes"}
# Earl Grey documents HELIANO integration as experimental; enable explicitly when required.
EARLGREY_HELIANO=${EARLGREY_HELIANO:-"yes"}
# FINAL_NR already contains the selected Dfam/RepBase consensus records.
EARLGREY_REPEATMASKER_SPECIES=${EARLGREY_REPEATMASKER_SPECIES:-"none"}
EARLGREY_ENHANCED_PLOTS=${EARLGREY_ENHANCED_PLOTS:-"TRUE"}
EARLGREY_ARCHIVE_INCOMPLETE=${EARLGREY_ARCHIVE_INCOMPLETE:-"TRUE"}
EARLGREY_MIGRATE_LEGACY_LAYOUT=${EARLGREY_MIGRATE_LEGACY_LAYOUT:-"TRUE"}
EARLGREY_RUNTIME_PATCH_ANNOTATION_CLI=${EARLGREY_RUNTIME_PATCH_ANNOTATION_CLI:-"TRUE"}
EARLGREY_VALIDATE_RESUME=${EARLGREY_VALIDATE_RESUME:-"yes"}  # Earl Grey -v yes: checksum-validated resume
# How to handle incomplete Earl Grey outputs: resume | archive | error
# resume keeps expensive RepeatMasker/HELIANO outputs and removes only downstream stale merge/defragmentation artifacts.
EARLGREY_INCOMPLETE_POLICY=${EARLGREY_INCOMPLETE_POLICY:-"resume"}
EARLGREY_PLOT_DPI=${EARLGREY_PLOT_DPI:-600}
EARLGREY_PLOT_TOP_FAMILIES=${EARLGREY_PLOT_TOP_FAMILIES:-30}
EARLGREY_PLOT_BINWIDTH=${EARLGREY_PLOT_BINWIDTH:-"0.01"}
EARLGREY_PLOT_MAX_KIMURA=${EARLGREY_PLOT_MAX_KIMURA:-"0.5"}
EARLGREY_COMPOSITION_TOLERANCE=${EARLGREY_COMPOSITION_TOLERANCE:-"0.005"}
NON_REPEATS_REPEATS_BED=${NON_REPEATS_REPEATS_BED:-""}

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
  ["busco_dir"]="${working_dir}/07.asm_stats/busco"
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
export busco_dir=${dirs["busco_dir"]}
export btk_dir=${dirs["btk_dir"]}
export repeats_dir=${dirs["repeats"]}
export earlgrey_dir="${repeats_dir}/05.earlgrey_annotation_only"

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
declare -a owned_lock_dirs=()
cleanup_pipeline_tmpdir() {
  local exit_status=$?
  local lock_dir
  trap - EXIT
  for lock_dir in "${owned_lock_dirs[@]}"; do
    rmdir -- "${lock_dir}" 2>/dev/null || true
  done
  if is_true "${cleanup_tmpdir:-FALSE}" && [[ -n "${tmpdir:-}" ]]; then
    rm -rf -- "${tmpdir:?}"
  fi
  exit "${exit_status}"
}
trap cleanup_pipeline_tmpdir EXIT

# ---- Embedded genoflowHiFi assessment engine -------------------------------------
# Keep the pipeline self-contained: the Python assessment code is stored in this
# Bash script and materialized only inside the pipeline-owned temporary directory.
# It is executed with Python from GENOFLOWSTATSenv and is removed automatically when the
# pipeline-owned TMPDIR is cleaned up. No external selector script is required.
ASSESS_SCRIPT="${tmpdir}/genoflow_assess.embedded.py"
cat > "${ASSESS_SCRIPT}" <<'__GENOFLOW_ASSESS_PY__'
#!/usr/bin/env python3
"""genoflowHiFi assembly assessment utilities.

Reference-free assembly and read-support metrics intended to replace the uses of
legacy external assembly/QC reporting in genoflowHiFi. Heavy alignment/sorting remains delegated
to minimap2/SAMtools; this module performs deterministic FASTA metrics, robust
coverage summaries, BAM alignment QC, candidate ranking and SVG plotting.

Python >=3.14; core dependencies: numpy, scipy, pysam, polars, pyarrow,
statsmodels, scikit-learn, matplotlib, rich.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import re
import shutil
import statistics
import subprocess
import sys
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable, Iterator, Sequence

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

try:
    import pysam
except Exception:
    pysam = None

try:
    from rich.console import Console
    console = Console(stderr=True)
except Exception:
    class _Console:
        def log(self, *args: Any, **kwargs: Any) -> None:
            print(*args, file=sys.stderr)
    console = _Console()

VERSION = "0.1.0"
HEAT = LinearSegmentedColormap.from_list(
    "genoflow_heat", ["#2b0a3d", "#7a1f59", "#d94841", "#f98e52", "#ffd166"]
)
HEAT_R = LinearSegmentedColormap.from_list(
    "genoflow_heat_r", ["#ffd166", "#f98e52", "#d94841", "#7a1f59", "#2b0a3d"]
)


def open_text(path: str | Path):
    p = str(path)
    return gzip.open(p, "rt", encoding="utf-8") if p.endswith(".gz") else open(p, "rt", encoding="utf-8")


def fasta_records(path: str | Path) -> Iterator[tuple[str, str]]:
    name: str | None = None
    chunks: list[str] = []
    with open_text(path) as fh:
        for raw in fh:
            if raw.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = raw[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(raw.strip())
        if name is not None:
            yield name, "".join(chunks)


def nxx(lengths: Sequence[int], denominator: float, pct: float) -> tuple[int, int]:
    if denominator <= 0 or not lengths:
        return 0, 0
    target = denominator * (pct / 100.0)
    c = 0
    for i, length in enumerate(sorted(lengths, reverse=True), start=1):
        c += length
        if c >= target:
            return int(length), i
    return 0, 0


def exact_quantile_from_hist(hist: np.ndarray, q: float) -> float:
    total = int(hist.sum())
    if total <= 0:
        return math.nan
    # numpy default quantile uses interpolation. For integer depth QC we use the
    # inverse empirical CDF, which is exact and directly interpretable.
    rank = max(1, int(math.ceil(q * total)))
    cs = np.cumsum(hist, dtype=np.int64)
    return float(np.searchsorted(cs, rank, side="left"))


def exact_mad_from_hist(hist: np.ndarray, median: int) -> float:
    total = int(hist.sum())
    if total <= 0:
        return math.nan
    maxdev = max(median, len(hist) - 1 - median)
    dev = np.zeros(maxdev + 1, dtype=np.int64)
    nz = np.flatnonzero(hist)
    for d in nz:
        dev[abs(int(d) - median)] += int(hist[d])
    return exact_quantile_from_hist(dev, 0.5)


def write_kv_tsv(path: Path, values: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["metric", "value"])
        for k, v in values.items():
            w.writerow([k, "NA" if v is None or (isinstance(v, float) and not math.isfinite(v)) else v])


def write_json(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=2, sort_keys=True, allow_nan=False)
        fh.write("\n")


def style_ax(ax: plt.Axes) -> None:
    ax.set_facecolor("#fffaf0")
    ax.figure.set_facecolor("white")
    ax.grid(True, color="#efe6d5", linewidth=0.8, alpha=0.9)
    ax.set_axisbelow(True)
    for s in ax.spines.values():
        s.set_visible(False)


def save_svg(fig: plt.Figure, path: Path, dpi: int = 600) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, format="svg", dpi=dpi, bbox_inches="tight")
    plt.close(fig)


# -----------------------------------------------------------------------------
# FASTA metrics
# -----------------------------------------------------------------------------

def assembly_metrics(path: str | Path, genome_size: int | None, gap_min: int) -> tuple[dict[str, Any], list[dict[str, Any]], list[tuple[int, int, int]]]:
    lengths: list[int] = []
    contig_lengths: list[int] = []
    scaffold_rows: list[dict[str, Any]] = []
    gap_rows: list[tuple[int, int, int]] = []
    total_gc = total_atgc = total_n = total_ambig = 0

    gap_re = re.compile(rf"N{{{gap_min},}}", re.I)
    for name, seq in fasta_records(path):
        s = seq.upper()
        L = len(s)
        lengths.append(L)
        n = s.count("N")
        gc = s.count("G") + s.count("C")
        atgc = s.count("A") + s.count("T") + gc
        ambig = L - atgc - n
        total_gc += gc; total_atgc += atgc; total_n += n; total_ambig += max(0, ambig)

        gaps = list(gap_re.finditer(s))
        for m in gaps:
            gap_rows.append((len(lengths) - 1, m.start(), m.end()))
        pieces = [x for x in gap_re.split(s) if x]
        contig_lengths.extend(len(x) for x in pieces)
        scaffold_rows.append({
            "scaffold": name, "length": L, "gc_percent": (100.0 * gc / atgc if atgc else math.nan),
            "n_bases": n, "gap_count": len(gaps), "gap_bp": sum(m.end()-m.start() for m in gaps),
        })

    total = int(sum(lengths))
    if not lengths:
        raise ValueError(f"No FASTA records found: {path}")
    largest = int(max(lengths))
    n50, l50 = nxx(lengths, total, 50); n90, l90 = nxx(lengths, total, 90)
    g = int(genome_size) if genome_size and genome_size > 0 else None
    ng50, lg50 = nxx(lengths, g, 50) if g else (0, 0)
    ng90, lg90 = nxx(lengths, g, 90) if g else (0, 0)
    aun = float(sum(x*x for x in lengths) / total) if total else 0.0
    aung = float(sum(x*x for x in lengths) / g) if g else math.nan
    cn50, cl50 = nxx(contig_lengths, sum(contig_lengths), 50) if contig_lengths else (0, 0)
    cn90, cl90 = nxx(contig_lengths, sum(contig_lengths), 90) if contig_lengths else (0, 0)
    gap_bp = int(sum(e-s for _, s, e in gap_rows))
    gap_lengths = np.array([e-s for _, s, e in gap_rows], dtype=np.int64)
    sorted_lengths = sorted(lengths, reverse=True)

    def topfrac(n: int) -> float:
        return float(sum(sorted_lengths[:n]) / total) if total else math.nan

    vals: dict[str, Any] = {
        "assembly": Path(path).name,
        "scaffold_count": len(lengths), "total_length": total, "canonical_bases": total_atgc,
        "n_bases": total_n, "ambiguous_non_N_bases": total_ambig,
        "gc_percent": 100.0 * total_gc / total_atgc if total_atgc else math.nan,
        "largest_scaffold": largest, "mean_scaffold_length": float(np.mean(lengths)),
        "median_scaffold_length": float(np.median(lengths)),
        "n50": n50, "l50": l50, "n90": n90, "l90": l90,
        "ng50": ng50 if g else None, "lg50": lg50 if g else None,
        "ng90": ng90 if g else None, "lg90": lg90 if g else None,
        "auN": aun, "auNG": aung if g else None,
        "genome_size_estimate": g,
        "assembly_to_genome_ratio": total / g if g else None,
        "assembly_vs_genome_percent": 100.0 * (total-g) / g if g else None,
        "contig_count_gap_split": len(contig_lengths), "contig_n50": cn50, "contig_l50": cl50,
        "contig_n90": cn90, "contig_l90": cl90,
        "gap_min_run": gap_min, "gap_count": len(gap_rows), "gap_bp": gap_bp,
        "gap_fraction": gap_bp / total if total else math.nan,
        "gaps_per_100mb": len(gap_rows) * 1e8 / total if total else math.nan,
        "ns_per_100kbp": total_n * 1e5 / total if total else math.nan,
        "median_gap": float(np.median(gap_lengths)) if gap_lengths.size else 0.0,
        "mean_gap": float(np.mean(gap_lengths)) if gap_lengths.size else 0.0,
        "p95_gap": float(np.percentile(gap_lengths, 95)) if gap_lengths.size else 0.0,
        "largest_gap": int(gap_lengths.max()) if gap_lengths.size else 0,
        "top1_fraction": topfrac(1), "top5_fraction": topfrac(5), "top10_fraction": topfrac(10),
        "top20_fraction": topfrac(20), "top50_fraction": topfrac(50),
    }
    return vals, scaffold_rows, gap_rows


def plot_fasta_metrics(vals: dict[str, Any], scaffold_rows: list[dict[str, Any]], outdir: Path, prefix: str, dpi: int) -> None:
    lens = np.array([r["length"] for r in scaffold_rows], dtype=np.int64)
    s = np.sort(lens)[::-1]
    cum = np.cumsum(s) / s.sum()
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(np.arange(1, len(s)+1), cum, color="#d94841", linewidth=2)
    ax.fill_between(np.arange(1, len(s)+1), cum, color="#f98e52", alpha=0.25)
    ax.set_xscale("log"); ax.set_ylim(0, 1.01)
    ax.set_xlabel("Scaffold rank (log scale)"); ax.set_ylabel("Cumulative assembly fraction")
    ax.set_title(f"{prefix}: cumulative scaffold length")
    style_ax(ax); save_svg(fig, outdir / f"{prefix}.scaffold_cumulative.svg", dpi)

    x = np.arange(1, 101)
    nx = np.array([nxx(lens.tolist(), int(lens.sum()), p)[0] for p in x], dtype=float)
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(x, nx, color="#7a1f59", linewidth=2, label="Nx")
    if vals.get("genome_size_estimate"):
        ng = np.array([nxx(lens.tolist(), int(vals["genome_size_estimate"]), p)[0] for p in x], dtype=float)
        ax.plot(x, ng, color="#f98e52", linewidth=2, label="NGx")
    ax.set_xlabel("x (%)"); ax.set_ylabel("Length (bp)"); ax.set_title(f"{prefix}: Nx/NGx profile")
    ax.legend(frameon=False); style_ax(ax); save_svg(fig, outdir / f"{prefix}.nx_ngx.svg", dpi)


# -----------------------------------------------------------------------------
# Streaming depth histogram
# -----------------------------------------------------------------------------

def iter_depths(path: str | Path, mode: str) -> Iterator[int]:
    with open_text(path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                vals = [int(x) for x in parts[2:]]
            except ValueError:
                continue
            if not vals:
                continue
            if mode == "single":
                yield vals[0]
            elif mode == "sum":
                yield sum(vals)
            else:
                raise ValueError("depth mode 'each' is not supported by this scalar histogram command")


def build_depth_hist(path: str | Path, mode: str, include_zero: bool) -> tuple[np.ndarray, int]:
    counts: list[int] = [0] * 128
    n = 0
    for d in iter_depths(path, mode):
        if d < 0 or (d == 0 and not include_zero):
            continue
        if d >= len(counts):
            new = len(counts)
            while new <= d:
                new *= 2
            counts.extend([0] * (new - len(counts)))
        counts[d] += 1; n += 1
    arr = np.asarray(counts, dtype=np.int64)
    nz = np.flatnonzero(arr)
    if nz.size:
        arr = arr[: int(nz[-1]) + 1]
    return arr, n


def depth_metrics_from_hist(hist: np.ndarray, scale_mad: bool, mad_factor: float, low_pct: float, high_pct: float) -> dict[str, Any]:
    n = int(hist.sum())
    if n == 0:
        raise ValueError("No valid depth values")
    med = int(exact_quantile_from_hist(hist, 0.5))
    mad_raw = float(exact_mad_from_hist(hist, med))
    mad = mad_raw * 1.4826 if scale_mad else mad_raw
    mode = int(np.argmax(hist))
    q = {p: exact_quantile_from_hist(hist, p/100.0) for p in (1, 5, 25, 50, 75, 95, 99)}
    mean = float(np.dot(np.arange(hist.size, dtype=np.float64), hist.astype(np.float64)) / n)
    total_zero = int(hist[0]) if hist.size else 0
    vals: dict[str, Any] = {
        "n_sites": n, "zero_depth_sites": total_zero, "zero_fraction": total_zero/n,
        "mean": mean, "median": med, "mode": mode, "mad_raw": mad_raw, "mad_used": mad,
        "iqr": q[75]-q[25], "robust_cv": (mad/med if med > 0 else None),
        "p01": q[1], "p05": q[5], "p25": q[25], "p50": q[50], "p75": q[75], "p95": q[95], "p99": q[99],
        "min_cutoff_mad": max(0.0, med-mad_factor*mad), "max_cutoff_mad": med+mad_factor*mad,
        "min_cutoff_percentile": exact_quantile_from_hist(hist, low_pct/100.0),
        "max_cutoff_percentile": exact_quantile_from_hist(hist, high_pct/100.0),
    }
    cs_rev = np.cumsum(hist[::-1], dtype=np.int64)[::-1]
    for t in (1, 5, 10, 20, 30, 50, 100):
        vals[f"fraction_ge_{t}x"] = float(cs_rev[t]/n) if t < cs_rev.size else 0.0
    for label, mult in (("025", .25), ("050", .5), ("150", 1.5), ("200", 2.0)):
        t = int(math.ceil(med * mult))
        vals[f"threshold_{label}_median"] = t
        vals[f"fraction_ge_{label}_median"] = float(cs_rev[t]/n) if t < cs_rev.size else 0.0
    vals["fraction_lt_half_median"] = float(hist[: max(0, int(math.floor(.5*med)))+1].sum()/n) if med else None
    vals["fraction_gt_twice_median"] = float(hist[int(math.floor(2*med))+1:].sum()/n) if med else None
    return vals


def plot_depth_hist(hist: np.ndarray, vals: dict[str, Any], outdir: Path, prefix: str, dpi: int) -> None:
    x = np.arange(hist.size)
    n = hist.sum(); density = hist / n
    xmax = max(5, int(vals["p99"] * 1.15), int(vals["max_cutoff_mad"] * 1.1))
    xmax = min(xmax, hist.size-1)
    mask = x <= xmax
    fig, ax = plt.subplots(figsize=(10, 6))
    colors = HEAT(np.linspace(0.20, 0.95, int(mask.sum())))
    ax.bar(x[mask], density[mask], color=colors, width=1.0, linewidth=0)
    ax.axvline(vals["median"], color="#2b0a3d", linewidth=1.8, label=f"Median={vals['median']}")
    ax.axvspan(vals["min_cutoff_mad"], vals["max_cutoff_mad"], color="#ffd166", alpha=.22, label="MAD interval")
    ax.set_xlabel("Depth"); ax.set_ylabel("Genome fraction"); ax.set_title(f"{prefix}: HiFi depth distribution")
    ax.legend(frameon=False); style_ax(ax); save_svg(fig, outdir / f"{prefix}.depth_hist.svg", dpi)

    cdf = np.cumsum(hist, dtype=np.int64) / n
    surv = np.cumsum(hist[::-1], dtype=np.int64)[::-1] / n
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(x, cdf, color="#7a1f59", linewidth=2, label="CDF")
    ax.plot(x, surv, color="#f98e52", linewidth=2, label="Coverage fraction ≥ depth")
    ax.set_xlim(0, xmax); ax.set_ylim(0, 1.01); ax.set_xlabel("Depth"); ax.set_ylabel("Fraction")
    ax.set_title(f"{prefix}: cumulative HiFi depth"); ax.legend(frameon=False)
    style_ax(ax); save_svg(fig, outdir / f"{prefix}.depth_cumulative.svg", dpi)


# -----------------------------------------------------------------------------
# BAM QC: one streaming pass, no depth calculation
# -----------------------------------------------------------------------------

def bam_metrics(path: str | Path, include_secondary: bool = False) -> tuple[dict[str, Any], np.ndarray, np.ndarray]:
    if pysam is None:
        raise RuntimeError("pysam is required for BAM metrics")
    totals = Counter(); mapq_hist = np.zeros(256, dtype=np.int64); indel_hist = np.zeros(1002, dtype=np.int64)
    aligned_bases = mismatch_nm = ins_bases = del_bases = softclip = hardclip = 0
    with pysam.AlignmentFile(str(path), "rb") as bam:
        for r in bam.fetch(until_eof=True):
            totals["records"] += 1
            if r.is_unmapped:
                totals["unmapped"] += 1; continue
            totals["mapped"] += 1
            if r.is_secondary:
                totals["secondary"] += 1
                if not include_secondary:
                    continue
            if r.is_supplementary:
                totals["supplementary"] += 1
            if r.is_duplicate: totals["duplicate_flag"] += 1
            if r.is_qcfail: totals["qcfail"] += 1
            mq = min(255, max(0, int(r.mapping_quality))); mapq_hist[mq] += 1
            if r.cigartuples:
                for op, ln in r.cigartuples:
                    if op in (0, 7, 8): aligned_bases += ln
                    elif op == 1:
                        ins_bases += ln; indel_hist[min(1001, ln)] += 1
                    elif op == 2:
                        del_bases += ln; indel_hist[min(1001, ln)] += 1
                    elif op == 4: softclip += ln
                    elif op == 5: hardclip += ln
            try:
                mismatch_nm += int(r.get_tag("NM"))
            except KeyError:
                totals["missing_nm"] += 1
    primary_denom = max(1, totals["mapped"])
    considered = int(mapq_hist.sum())
    vals = {
        "alignment_records": totals["records"], "mapped_records": totals["mapped"], "unmapped_records": totals["unmapped"],
        "mapped_fraction": totals["mapped"]/max(1, totals["records"]),
        "secondary_records": totals["secondary"], "secondary_fraction_mapped": totals["secondary"]/primary_denom,
        "supplementary_records": totals["supplementary"], "supplementary_fraction_mapped": totals["supplementary"]/primary_denom,
        "duplicate_flag_records": totals["duplicate_flag"], "qcfail_records": totals["qcfail"],
        "mapq0_fraction": mapq_hist[0]/max(1, considered),
        "mapq_lt10_fraction": mapq_hist[:10].sum()/max(1, considered),
        "mapq_lt20_fraction": mapq_hist[:20].sum()/max(1, considered),
        "mapq_ge30_fraction": mapq_hist[30:].sum()/max(1, considered),
        "mapq_ge50_fraction": mapq_hist[50:].sum()/max(1, considered),
        "mapq_median": weighted_quantile_hist(mapq_hist, .5) if considered else None,
        "aligned_bases": aligned_bases, "nm_total": mismatch_nm, "inserted_bases": ins_bases,
        "deleted_bases": del_bases, "softclipped_bases": softclip, "hardclipped_bases": hardclip,
        "missing_nm_records": totals["missing_nm"],
    }
    # NM includes mismatches plus indels for common aligners; report conservative NM/aligned concordance separately.
    e = mismatch_nm / aligned_bases if aligned_bases else math.nan
    vals["nm_rate_per_aligned_base"] = e
    vals["alignment_concordance_qv_nm"] = -10*math.log10(e) if e and e > 0 else (math.inf if e == 0 else None)
    return vals, mapq_hist, indel_hist


def weighted_quantile_hist(hist: np.ndarray, q: float) -> float:
    return exact_quantile_from_hist(hist, q)


def plot_bam_metrics(mapq: np.ndarray, indel: np.ndarray, outdir: Path, prefix: str, dpi: int) -> None:
    if mapq.sum():
        x = np.arange(mapq.size); y = mapq/mapq.sum(); nz = y > 0
        fig, ax = plt.subplots(figsize=(10, 6)); colors = HEAT(np.linspace(.15,.95,int(nz.sum())))
        ax.bar(x[nz], y[nz], color=colors, width=.9); ax.set_xlabel("MAPQ"); ax.set_ylabel("Alignment fraction")
        ax.set_title(f"{prefix}: mapping-quality distribution"); style_ax(ax); save_svg(fig, outdir/f"{prefix}.mapq.svg", dpi)
    if indel.sum():
        x=np.arange(indel.size); mask=(indel>0)&(x<=100); y=indel/max(1,indel.sum())
        fig,ax=plt.subplots(figsize=(10,6)); colors=HEAT(np.linspace(.15,.95,int(mask.sum())))
        ax.bar(x[mask],y[mask],color=colors,width=.9); ax.set_xlabel("Indel length (bp; ≤100)"); ax.set_ylabel("Event fraction")
        ax.set_title(f"{prefix}: indel-size spectrum"); style_ax(ax); save_svg(fig,outdir/f"{prefix}.indel_spectrum.svg",dpi)


# -----------------------------------------------------------------------------
# Candidate ranking
# -----------------------------------------------------------------------------

def read_table(path: str | Path) -> list[dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def fnum(v: Any) -> float:
    try: return float(str(v).replace(",", ""))
    except Exception: return math.nan


def robust_z(v: np.ndarray, direction: float = 1.0) -> np.ndarray:
    out=np.zeros(v.size,float); finite=np.isfinite(v)
    if not finite.any(): return out
    x=v[finite]; med=np.median(x); mad=np.median(np.abs(x-med))*1.4826
    if not np.isfinite(mad) or mad==0:
        mad=np.std(x)
    if not np.isfinite(mad) or mad==0: mad=1.0
    out[finite]=direction*(x-med)/mad
    return out


def load_optional_simple(path: str | Path | None, id_col: str, wanted: Sequence[str]) -> dict[str, dict[str, float]]:
    if not path or not Path(path).is_file(): return {}
    rows=read_table(path); out={}
    for r in rows:
        k=Path(r.get(id_col,"")).stem
        out[k]={w:fnum(r.get(w)) for w in wanted}
    return out


def rank_candidates(metrics_path: str | Path, manifest_path: str | Path, outdir: Path, top_n: int,
                    genome_size: float | None, yak_path: str | None, merqury_dir: str | None, final: bool, dpi: int) -> None:
    metrics=read_table(metrics_path); manifest=read_table(manifest_path)
    mm={Path(r["assembly_id"]).stem:r for r in metrics}
    yak=load_optional_simple(yak_path,"assembly_id",["yak_qv"])
    mq: dict[str,dict[str,float]]={}
    if merqury_dir and Path(merqury_dir).is_dir():
        for p in Path(merqury_dir).rglob("merqury_metrics.tsv"):
            for r in read_table(p):
                mq[Path(r.get("assembly_id","")).stem]={"merqury_qv":fnum(r.get("merqury_qv")),"merqury_completeness":fnum(r.get("merqury_completeness"))}
    rows=[]
    for m in manifest:
        key=Path(m["assembly_id"]).stem; r=dict(m); r.update(mm.get(key,{})); r.update(yak.get(key,{})); r.update(mq.get(key,{})); rows.append(r)
    if not rows: raise ValueError("No candidates to rank")

    def col(name:str)->np.ndarray: return np.array([fnum(r.get(name)) for r in rows],float)
    # Evidence groups reduce double-counting of correlated metrics.
    contig_base = col("auNG") if np.isfinite(col("auNG")).any() else col("auN")
    score_contig=0.20*robust_z(contig_base,1)
    gap_metric=col("gap_fraction") if np.isfinite(col("gap_fraction")).any() else col("ns_per_100kbp")
    score_gap=0.15*robust_z(gap_metric,-1)
    if genome_size and genome_size>0:
        size_pen=np.abs(col("total_length")-genome_size)/genome_size
        score_size=0.15*robust_z(size_pen,-1)
    else: score_size=np.zeros(len(rows))
    score_kmer=np.zeros(len(rows))
    if final and np.isfinite(col("merqury_qv")).any(): score_kmer += .25*robust_z(col("merqury_qv"),1)
    if final and np.isfinite(col("merqury_completeness")).any(): score_kmer += .20*robust_z(col("merqury_completeness"),1)
    if np.isfinite(col("yak_qv")).any(): score_kmer += .05*robust_z(col("yak_qv"),1)
    score=score_contig+score_gap+score_size+score_kmer
    for i,r in enumerate(rows):
        r.update({"score_contiguity":score_contig[i],"score_gap":score_gap[i],"score_size":score_size[i],"score_kmer":score_kmer[i],"score":score[i]})
        if not Path(r.get("path","")).is_file() or not math.isfinite(fnum(r.get("n50"))): r["score"]=-math.inf
    rows.sort(key=lambda r:(-fnum(r.get("score")),-fnum(r.get("auNG")),-fnum(r.get("n50"))))
    outdir.mkdir(parents=True,exist_ok=True)
    cols=list(dict.fromkeys(k for r in rows for k in r.keys()))
    name="longstitch_scores.tsv" if final else "preselection_scores.tsv"
    with (outdir/name).open("w",newline="",encoding="utf-8") as fh:
        w=csv.DictWriter(fh,fieldnames=cols,delimiter="\t",lineterminator="\n",extrasaction="ignore"); w.writeheader(); w.writerows(rows)
    valid=[r for r in rows if math.isfinite(fnum(r.get("score")))][:top_n]
    (outdir/"top_candidates.txt").write_text("".join(f"{r['path']}\n" for r in valid),encoding="utf-8")
    if valid:
        (outdir/"best_assembly.path").write_text(valid[0]["path"]+"\n",encoding="utf-8")
        (outdir/"best_assembly.txt").write_text(valid[0]["assembly_id"]+"\n",encoding="utf-8")
    plot_candidate_ranking(rows,valid,outdir,final,dpi)


def plot_candidate_ranking(rows:list[dict[str,Any]], top:list[dict[str,Any]], outdir:Path, final:bool, dpi:int)->None:
    if top:
        labels=[r["assembly_id"] for r in top][::-1]; vals=np.array([fnum(r["score"]) for r in top][::-1])
        fig,ax=plt.subplots(figsize=(11,max(5,.38*len(labels)+2))); colors=HEAT(np.linspace(.25,.95,len(labels)))
        ax.barh(labels,vals,color=colors); ax.set_xlabel("Robust evidence-group score"); ax.set_title("Final scaffold candidate ranking" if final else "Scaffold preselection ranking")
        style_ax(ax); save_svg(fig,outdir/("ranking.svg" if final else "preselection_ranking.svg"),dpi)
    # k/w heatmap where metadata exist
    plot=[r for r in rows if math.isfinite(fnum(r.get("score"))) and math.isfinite(fnum(r.get("k"))) and math.isfinite(fnum(r.get("w")))]
    if plot:
        ks=sorted(set(fnum(r["k"]) for r in plot)); ws=sorted(set(fnum(r["w"]) for r in plot)); mat=np.full((len(ks),len(ws)),np.nan)
        for r in plot: mat[ks.index(fnum(r["k"])),ws.index(fnum(r["w"]))]=fnum(r["score"])
        fig,ax=plt.subplots(figsize=(10,7)); im=ax.imshow(mat,aspect="auto",origin="lower",cmap=HEAT)
        ax.set_xticks(range(len(ws)),[str(int(x)) for x in ws]); ax.set_yticks(range(len(ks)),[str(int(x)) for x in ks])
        ax.set_xlabel("LongStitch w"); ax.set_ylabel("LongStitch k"); ax.set_title("Candidate score heatmap")
        fig.colorbar(im,ax=ax,label="Score"); style_ax(ax); save_svg(fig,outdir/("score_heatmap.svg" if final else "preselection_score_heatmap.svg"),dpi)


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def cmd_fasta(a: argparse.Namespace) -> None:
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    if a.manifest:
        man=read_table(a.manifest); rows=[]
        for r in man:
            v,_,_=assembly_metrics(r["path"],a.genome_size,a.gap_min); v["assembly_id"]=r["assembly_id"]; v["path"]=r["path"]
            for k in ("k","w","round","mode"):
                if k in r: v[k]=r[k]
            rows.append(v)
        cols=list(dict.fromkeys(k for r in rows for k in r.keys()))
        with (out/a.table_name).open("w",newline="",encoding="utf-8") as fh:
            w=csv.DictWriter(fh,fieldnames=cols,delimiter="\t",lineterminator="\n"); w.writeheader(); w.writerows(rows)
        return
    vals,scaffolds,_=assembly_metrics(a.fasta,a.genome_size,a.gap_min); prefix=a.prefix or Path(a.fasta).stem
    write_kv_tsv(out/f"{prefix}.assembly_metrics.tsv",vals); write_json(out/f"{prefix}.assembly_metrics.json",vals)
    with (out/f"{prefix}.scaffolds.tsv").open("w",newline="",encoding="utf-8") as fh:
        w=csv.DictWriter(fh,fieldnames=list(scaffolds[0]),delimiter="\t",lineterminator="\n"); w.writeheader(); w.writerows(scaffolds)
    if not a.no_plots: plot_fasta_metrics(vals,scaffolds,out,prefix,a.dpi)


def cmd_depth(a: argparse.Namespace) -> None:
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True); prefix=a.prefix or Path(a.input).name.removesuffix(".gz")
    hist,n=build_depth_hist(a.input,a.depth_mode,a.include_zero)
    vals=depth_metrics_from_hist(hist,a.scale_mad,a.mad_factor,a.low_percentile,a.high_percentile)
    write_kv_tsv(out/f"{prefix}.depth_metrics.tsv",vals); write_json(out/f"{prefix}.depth_metrics.json",vals)
    with (out/f"{prefix}.depth_histogram.tsv").open("w",encoding="utf-8") as fh:
        fh.write("depth\tsites\n");
        for d,c in enumerate(hist):
            if c: fh.write(f"{d}\t{int(c)}\n")
    if not a.no_plots: plot_depth_hist(hist,vals,out,prefix,a.dpi)


def cmd_bam(a: argparse.Namespace) -> None:
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True); prefix=a.prefix or Path(a.bam).stem
    vals,mapq,indel=bam_metrics(a.bam,a.include_secondary); write_kv_tsv(out/f"{prefix}.bam_metrics.tsv",vals); write_json(out/f"{prefix}.bam_metrics.json",vals)
    if not a.no_plots: plot_bam_metrics(mapq,indel,out,prefix,a.dpi)


def build_parser()->argparse.ArgumentParser:
    p=argparse.ArgumentParser(prog="genoflow_assess",description="genoflowHiFi assembly QC metrics")
    p.add_argument("--version",action="version",version=f"%(prog)s {VERSION}")
    sub=p.add_subparsers(dest="cmd",required=True)
    f=sub.add_parser("fasta"); f.add_argument("--fasta"); f.add_argument("--manifest"); f.add_argument("--outdir",required=True); f.add_argument("--prefix"); f.add_argument("--table-name",default="candidate_metrics.tsv"); f.add_argument("--genome-size",type=int); f.add_argument("--gap-min",type=int,default=10); f.add_argument("--dpi",type=int,default=600); f.add_argument("--no-plots",action="store_true"); f.set_defaults(func=cmd_fasta)
    d=sub.add_parser("depth"); d.add_argument("--input",required=True); d.add_argument("--outdir",required=True); d.add_argument("--prefix"); d.add_argument("--depth-mode",choices=["single","sum"],default="single"); d.add_argument("--include-zero",action="store_true"); d.add_argument("--scale-mad",action="store_true"); d.add_argument("--mad-factor",type=float,default=3.0); d.add_argument("--low-percentile",type=float,default=1.0); d.add_argument("--high-percentile",type=float,default=99.0); d.add_argument("--dpi",type=int,default=600); d.add_argument("--no-plots",action="store_true"); d.set_defaults(func=cmd_depth)
    b=sub.add_parser("bam"); b.add_argument("--bam",required=True); b.add_argument("--outdir",required=True); b.add_argument("--prefix"); b.add_argument("--include-secondary",action="store_true"); b.add_argument("--dpi",type=int,default=600); b.add_argument("--no-plots",action="store_true"); b.set_defaults(func=cmd_bam)
    r=sub.add_parser("candidates"); r.add_argument("--metrics",required=True); r.add_argument("--manifest",required=True); r.add_argument("--outdir",required=True); r.add_argument("--top-n",type=int,default=10); r.add_argument("--genome-size",type=float); r.add_argument("--yak"); r.add_argument("--merqury-dir"); r.add_argument("--final",action="store_true"); r.add_argument("--dpi",type=int,default=600)
    def _rank(a): rank_candidates(a.metrics,a.manifest,Path(a.outdir),a.top_n,a.genome_size,a.yak,a.merqury_dir,a.final,a.dpi)
    r.set_defaults(func=_rank)
    return p


def main(argv:Sequence[str]|None=None)->int:
    a=build_parser().parse_args(argv)
    if getattr(a,"cmd",None)=="fasta" and not a.fasta and not a.manifest:
        raise SystemExit("fasta command requires --fasta or --manifest")
    try:
        a.func(a); return 0
    except Exception as e:
        console.log(f"[ERROR] {e}"); return 1

if __name__=="__main__": raise SystemExit(main())
__GENOFLOW_ASSESS_PY__
chmod 600 "${ASSESS_SCRIPT}"


release_lock_dir() {
  local target="${1:?missing lock directory}"
  local item
  local -a retained=()
  rmdir -- "${target}" 2>/dev/null || true
  for item in "${owned_lock_dirs[@]}"; do
    [[ "${item}" == "${target}" ]] || retained+=("${item}")
  done
  owned_lock_dirs=("${retained[@]}")
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
# Verify the installed implementation, not just its version: 0.2.0 gained
# reference-FASTA and comparison support after its initial release.
genoflowstats_implementation_id() {
  local python_bin="${1:?missing genoflowSTATS Python}"
  "${python_bin}" - <<'PY_GENOFLOWSTATS_ID'
import hashlib, inspect
from pathlib import Path
import genoflowstats
from genoflowstats import cli

parameters = inspect.signature(cli.samtools_stats_metrics).parameters
if not {'reference_fasta', 'known_fasta_length'} <= set(parameters):
    raise SystemExit('genoflowSTATS lacks automatic reference FASTA support')
if not all(callable(getattr(cli, name, None)) for name in ('fasta_total_length', 'stats_reference_command')):
    raise SystemExit('genoflowSTATS lacks automatic reference FASTA support')
parser = cli.build_parser()
subcommands = next(a.choices for a in parser._actions if getattr(a, 'dest', None) == 'cmd')
for name in ('fasta', 'stats', 'compare'):
    if '--palette' not in subcommands[name]._option_string_actions:
        raise SystemExit('genoflowSTATS lacks palette support')
if '--reference-fasta' not in subcommands['stats']._option_string_actions:
    raise SystemExit('genoflowSTATS lacks --reference-fasta')
hash_value = hashlib.sha256()
for module in (genoflowstats, cli):
    path = Path(module.__file__)
    hash_value.update(path.name.encode('utf-8') + b'\0')
    hash_value.update(path.read_bytes())
print(hash_value.hexdigest())
PY_GENOFLOWSTATS_ID
}

require() {
  local bin="${1:-}"
  [[ -n "${bin}" ]] || fatal "require: missing executable path"
  if [[ "${bin}" == */* ]]; then
    [[ -x "${bin}" ]] || fatal "Missing required executable: ${bin}"
  else
    command -v "${bin}" >/dev/null 2>&1 || fatal "Missing required executable: ${bin}"
  fi
}

# Resolve the GenomeScope2 wrapper directly from the environment prefix.
# This avoids `mamba run`/Conda activation and guarantees use of the wrapper
# that pins GENOMESCOPE2env/bin/Rscript.
resolve_genomescope2_bin() {
  local configured="${GENOMESCOPE2_BIN:-genomescope.R}"
  local env_prefix="${GENOMESCOPE2_ENV_PREFIX}"
  local candidate=""

  if [[ "${configured}" == */* ]]; then
    [[ -x "${configured}" ]] || return 1
    printf '%s\n' "${configured}"
    return 0
  fi

  candidate="${env_prefix}/bin/${configured}"
  if [[ -x "${candidate}" ]]; then
    printf '%s\n' "${candidate}"
    return 0
  fi

  # Canonical wrapper installed by create_genomescope2_env().
  candidate="${env_prefix}/bin/genomescope.R"
  if [[ -x "${candidate}" ]]; then
    printf '%s\n' "${candidate}"
    return 0
  fi

  return 1
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

# ---- BUSCO dataset preparation --------------------------------------------------
# Run the installed entry point without modifying BUSCO or its environment.
# urlretrieve() inherits the socket timeout, including TLS and response reads.
busco_bounded_run() {
  local deadline="${1:?missing deadline (0 disables the overall deadline)}"
  shift
  local -a runner_args=(run -n "${ENV_BUSCO}")
  [[ "${ENV_RUNNER##*/}" != conda ]] || runner_args+=(--no-capture-output)
  local -a launch=("${ENV_RUNNER}" "${runner_args[@]}" python -u - "${BUSCO_NETWORK_TIMEOUT}" "$@")
  if [[ "${deadline}" != 0 ]]; then
    launch=(timeout --signal=TERM --kill-after=30s "${deadline}s" "${launch[@]}")
  fi
  "${launch[@]}" <<'PY_BUSCO_BOUNDED'
import runpy, shutil, socket, sys
# Offline scoring keeps Python's default behavior for local IPC sockets.
if "--offline" not in sys.argv[2:]:
    socket.setdefaulttimeout(float(sys.argv[1]))
entry = shutil.which("busco")
if not entry:
    raise SystemExit("[ERROR] BUSCO entry point not found in its environment")
sys.argv = [entry, *sys.argv[2:]]
try:
    runpy.run_path(entry, run_name="__main__")
except TimeoutError as exc:
    raise SystemExit(f"[ERROR] BUSCO network operation timed out: {exc}")
PY_BUSCO_BOUNDED
}

# A directory/dataset.cfg alone is insufficient after interrupted extraction.
# Check every scoring marker and the Miniprot reference, not just one HMM.
busco_lineage_ready() {
  env_run "${ENV_BUSCO}" python - "${1:?missing lineage directory}" <<'PY_BUSCO_LINEAGE'
import gzip, pathlib, sys, zlib
p = pathlib.Path(sys.argv[1])
try:
    cfg = dict(line.strip().split("=", 1) for line in (p / "dataset.cfg").read_text().splitlines()
               if line.strip() and not line.lstrip().startswith("#"))
    n = int(cfg["number_of_BUSCOs"])
    if n < 1 or not cfg.get("creation_date") or not cfg.get("name") or not cfg.get("domain"):
        raise ValueError("incomplete dataset.cfg")
    markers = set()
    for filename in ("scores_cutoff", "lengths_cutoff"):
        rows = [line.split() for line in (p / filename).read_text().splitlines() if line.strip()]
        ids = {row[0] for row in rows}
        if len(rows) != n or len(ids) != n or (markers and ids != markers):
            raise ValueError(f"incomplete or inconsistent {filename}: expected {n} markers")
        markers = ids
    for marker in markers:
        if pathlib.Path(marker).name != marker or (p / "hmms" / f"{marker}.hmm").stat().st_size == 0:
            raise ValueError(f"missing/empty HMM: {marker}")
    ref = p / "refseq_db.faa"
    if ref.is_file() and ref.stat().st_size:
        pass
    else:
        with gzip.open(p / "refseq_db.faa.gz", "rb") as handle:
            size = 0
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                size += len(chunk)
            if not size:
                raise ValueError("empty Miniprot reference")
except (OSError, EOFError, ValueError, KeyError, zlib.error) as exc:
    print(f"[ERROR] Incomplete BUSCO lineage {p}: {exc}", file=sys.stderr)
    raise SystemExit(1)
PY_BUSCO_LINEAGE
}

# Sets BUSCO_PREPARED_LINEAGE in the caller. Keep explicit offline/local-path
# contracts; online named lineages must finish a verified BUSCO download first.
prepare_busco_lineage() {
  local downloads="${1:?missing download directory}"
  local download_log="${2:?missing download log}"
  local setting attempt rc url lock_dir archive
  for setting in BUSCO_NETWORK_TIMEOUT BUSCO_DOWNLOAD_TIMEOUT BUSCO_DOWNLOAD_ATTEMPTS; do
    [[ "${!setting}" =~ ^[1-9][0-9]*$ ]] || { error "${setting} must be a positive integer"; return 1; }
  done
  [[ "${BUSCO_DOWNLOAD_RETRY_DELAY}" =~ ^(0|[1-9][0-9]*)$ ]] \
    || { error "BUSCO_DOWNLOAD_RETRY_DELAY must be a non-negative integer"; return 1; }
  mkdir -p "${downloads}" || return 1
  downloads="$(cd -- "${downloads}" && pwd -P)" || return 1
  if [[ "${BUSCO_LINEAGE}" == */* ]]; then
    BUSCO_PREPARED_LINEAGE="${BUSCO_LINEAGE}"
    busco_lineage_ready "${BUSCO_PREPARED_LINEAGE}" \
      || { error "Explicit BUSCO lineage path is incomplete: ${BUSCO_LINEAGE}"; return 1; }
    BUSCO_PREPARED_LINEAGE="$(cd -- "${BUSCO_PREPARED_LINEAGE}" && pwd -P)" || return 1
    return 0
  fi
  [[ "${BUSCO_LINEAGE}" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*_odb[0-9]+(\.[0-9]+)?$ ]] \
    || { error "BUSCO_LINEAGE must be one lineage name or a local dataset path: ${BUSCO_LINEAGE}"; return 1; }
  BUSCO_PREPARED_LINEAGE="${downloads}/lineages/${BUSCO_LINEAGE}"
  if [[ "${BUSCO_DATA_POLICY}" == offline ]]; then
    busco_lineage_ready "${BUSCO_PREPARED_LINEAGE}" \
      || { error "BUSCO offline cache is incomplete; populate ${downloads} with ${BUSCO_LINEAGE} first"; return 1; }
    return 0
  fi
  command -v timeout >/dev/null 2>&1 || { error "GNU timeout is required for BUSCO dataset downloads"; return 1; }
  lock_dir="${downloads}/.genoflowHiFi.busco.lock"
  mkdir "${lock_dir}" 2>/dev/null || { error "BUSCO download cache is locked: ${lock_dir}"; return 1; }
  owned_lock_dirs+=("${lock_dir}")
  local -a download_args=(--download_path "${downloads}")
  [[ "${BUSCO_DATA_POLICY}" != update ]] || download_args+=(--update-data)
  for ((attempt=1; attempt<=BUSCO_DOWNLOAD_ATTEMPTS; attempt++)); do
    # Preserve incomplete caches; BUSCO otherwise mistakes a partial extraction
    # with a current dataset.cfg for a complete dataset and skips the download.
    if [[ -e "${BUSCO_PREPARED_LINEAGE}" ]] && ! busco_lineage_ready "${BUSCO_PREPARED_LINEAGE}" >> "${download_log}" 2>&1; then
      archive="${BUSCO_PREPARED_LINEAGE}.incomplete.${LOG_TS}.$$.${attempt}"
      mv -- "${BUSCO_PREPARED_LINEAGE}" "${archive}" || return 1
      warn "Preserved incomplete BUSCO lineage as ${archive}"
    fi
    url="${BUSCO_DOWNLOAD_BASE_URL}"
    if (( attempt % 2 == 0 )); then url="${BUSCO_DOWNLOAD_FALLBACK_URL}"; fi
    info "BUSCO dataset ${BUSCO_LINEAGE}: download attempt ${attempt}/${BUSCO_DOWNLOAD_ATTEMPTS}, server=${url}, deadline=${BUSCO_DOWNLOAD_TIMEOUT}s" \
      | tee -a "${download_log}"
    # --download consumes all following arguments and must remain last.
    if busco_bounded_run "${BUSCO_DOWNLOAD_TIMEOUT}" "${download_args[@]}" \
        --download_base_url "${url}" --download "${BUSCO_LINEAGE}" 2>&1 | tee -a "${download_log}"; then
      if busco_lineage_ready "${BUSCO_PREPARED_LINEAGE}" 2>&1 | tee -a "${download_log}"; then
        release_lock_dir "${lock_dir}"
        info "BUSCO lineage ready: ${BUSCO_PREPARED_LINEAGE}"
        return 0
      fi
      rc=1
    else
      rc=$?
    fi
    warn "BUSCO dataset preparation failed (rc=${rc}, attempt ${attempt}/${BUSCO_DOWNLOAD_ATTEMPTS}); see ${download_log}"
    if (( attempt < BUSCO_DOWNLOAD_ATTEMPTS )); then sleep "${BUSCO_DOWNLOAD_RETRY_DELAY}"; fi
  done
  release_lock_dir "${lock_dir}"
  error "BUSCO dataset ${BUSCO_LINEAGE} unavailable after ${BUSCO_DOWNLOAD_ATTEMPTS} bounded attempts; see ${download_log}"
  return 1
}

# Print failing command, file, and line; keep set -Eeuo pipefail effective
err_trap() {
  local rc=$?
  printf "%s [ERROR] Command failed (rc=%s) at %s:%s: %s\n" \
    "$(date +"%F %T")" "$rc" "${BASH_SOURCE[1]:-${BASH_SOURCE[0]}}" "${BASH_LINENO[0]:-0}" "${BASH_COMMAND}" >&2
  exit "$rc"
}
trap err_trap ERR

# Global preflight is minimal; each stage validates its optional dependencies.
command -v "$ENV_RUNNER" >/dev/null 2>&1 || fatal "$ENV_RUNNER not found in PATH"

# Optional: FCS containers when early/late decontam enabled
if is_true "${run_fcs_early}" || is_true "${run_final_asm}"; then
  command -v singularity >/dev/null 2>&1 || fatal "'singularity' not on PATH"
  [[ -f "$sif_fcs_adaptor" ]] || fatal "Missing container (sif_fcs_adaptor): $sif_fcs_adaptor"
  [[ -f "$sif_fcs_gx"      ]] || fatal "Missing container (sif_fcs_gx): $sif_fcs_gx"
fi

if is_true "${FCS_TINI_SUBREAPER:-FALSE}"; then
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

  # Robust download; preserve the caller's xtrace state.
  curl --fail --show-error --location \
       --retry 5 --retry-delay 3 \
       --connect-timeout 30 --max-time 900 \
       --continue-at - \
       "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz" \
       -o "${taxdump_tar}"

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

# Repeat annotation uses one authoritative taxonomic assignment.
# REPEAT_TAXON may intentionally represent a broader clade than the assembled
# organism (for example, NCBI Taxonomy ID 7777 for Chondrichthyes or 1489904
# for Carangaria) so non-model species can use a richer curated repeat set.
#
# REPEATS_SPECIES and DFAM_TAXON are synchronized deliberately.  This avoids
# mixing an exact assembly species with a broader Dfam/Earl Grey taxonomy.
#
# Resolution order:
#   1. explicit REPEAT_TAXON
#   2. existing REPEATS_SPECIES
#   3. existing DFAM_TAXON
#   4. TAXID
#   5. genome species only when REPEAT_TAXON_USE_GENOME_SPECIES=TRUE
if [[ -z "${REPEAT_TAXON:-}" ]]; then
  if [[ -n "${REPEATS_SPECIES:-}" ]]; then
    REPEAT_TAXON="${REPEATS_SPECIES}"
  elif [[ -n "${DFAM_TAXON:-}" ]]; then
    REPEAT_TAXON="${DFAM_TAXON}"
  elif [[ -n "${TAXID:-}" ]]; then
    REPEAT_TAXON="${TAXID}"
  elif is_true "${REPEAT_TAXON_USE_GENOME_SPECIES:-TRUE}"; then
    REPEAT_TAXON="${species_name_raw}"
  else
    fatal "REPEAT_TAXON is required when no repeat taxonomy can be inferred"
  fi
fi

[[ -n "${REPEAT_TAXON:-}" ]] || fatal "REPEAT_TAXON resolved to an empty value"
[[ "$REPEAT_TAXON" != *$'\n'* && "$REPEAT_TAXON" != *$'\r'* && "$REPEAT_TAXON" != *$'\t'* ]] || fatal \
  "REPEAT_TAXON cannot contain tab/newline characters"

REPEATS_SPECIES="${REPEAT_TAXON}"
DFAM_TAXON="${REPEAT_TAXON}"

export REPEAT_TAXON REPEATS_SPECIES DFAM_TAXON
echo "[INFO] Species resolved: ${species_name_raw} → ${species_name} (TAXID=${TAXID})"
echo "[INFO] Assembly base name: ${base}"
echo "[INFO] Configuration loaded. working_dir=${working_dir} base=${base} threads=${nthreads} jobs=${njobs}"

# Hifiasm can reuse its .bin caches only in uncompressed form. Restore them
# before a real rebuild; a valid Step 1 checkpoint never needs decompression.
hifiasm_restore_bin_caches() {
  local archive bin tmp
  for archive in "${dir1}/${base}.asm."*.bin.zst "${dir1}/${base}.asm."*.bin.gz; do
    [[ -s "$archive" ]] || continue
    if [[ "$archive" == *.zst ]]; then
      bin="${archive%.zst}"
    else
      bin="${archive%.gz}"
    fi
    [[ -e "$bin" ]] && continue
    tmp="${bin}.restore.$$"
    rm -f -- "$tmp"
    if [[ "$archive" == *.zst ]]; then
      require zstd
      zstd -q -t "$archive" || fatal "Invalid Hifiasm cache archive: $archive"
      zstd -q -d -c "$archive" > "$tmp" ||
        { rm -f -- "$tmp"; fatal "Could not restore Hifiasm cache: $archive"; }
    else
      require pigz
      pigz -t "$archive" || fatal "Invalid Hifiasm cache archive: $archive"
      pigz -dc "$archive" > "$tmp" ||
        { rm -f -- "$tmp"; fatal "Could not restore Hifiasm cache: $archive"; }
    fi
    [[ -s "$tmp" ]] || { rm -f -- "$tmp"; fatal "Restored Hifiasm cache is empty: $archive"; }
    mv -f -- "$tmp" "$bin"
    info "Restored Hifiasm cache: $bin"
  done
}

hifiasm_compress_bin_caches() {
  local THREADS="${nthreads}" bin archive tmp codec
  [[ "$THREADS" =~ ^[0-9]+$ ]] && (( THREADS > 0 )) ||
    fatal "Hifiasm cache compression requires a positive thread count"
  if command -v zstd >/dev/null 2>&1; then
    codec="zstd"
  elif command -v pigz >/dev/null 2>&1; then
    codec="pigz"
    warn "zstd unavailable; compressing Hifiasm caches with pigz --best"
  else
    fatal "Hifiasm cache compression requires zstd or pigz"
  fi

  for bin in "${dir1}"/*.bin; do
    [[ -e "$bin" ]] || continue
    [[ -s "$bin" ]] || fatal "Hifiasm cache is empty and cannot be compressed: $bin"
    if [[ "$codec" == "zstd" ]]; then
      archive="${bin}.zst"
    else
      archive="${bin}.gz"
    fi
    tmp="${archive}.tmp.$$"
    rm -f -- "$tmp"
    info "Compressing Hifiasm cache: $bin -> $archive"
    if [[ "$codec" == "zstd" ]]; then
      zstd -q -9 -T"${THREADS}" -c "$bin" > "$tmp" ||
        { rm -f -- "$tmp"; fatal "zstd compression failed: $bin"; }
      zstd -q -t "$tmp" ||
        { rm -f -- "$tmp"; fatal "zstd archive verification failed: $tmp"; }
    else
      pigz --best -p "$THREADS" -c "$bin" > "$tmp" ||
        { rm -f -- "$tmp"; fatal "pigz compression failed: $bin"; }
      pigz -t "$tmp" ||
        { rm -f -- "$tmp"; fatal "pigz archive verification failed: $tmp"; }
    fi
    [[ -s "$tmp" ]] || { rm -f -- "$tmp"; fatal "Compressed Hifiasm cache is empty: $bin"; }
    mv -f -- "$tmp" "$archive"
    rm -f -- "$bin"
    if [[ "$codec" == "zstd" ]]; then rm -f -- "${bin}.gz"; else rm -f -- "${bin}.zst"; fi
  done
}

# ============================================================================================================================
# STEP 1: Hifiasm — primary haplotype-aware assembly from PacBio HiFi reads
# Shared implementation for genoflowHiFi and genoflowHiFi_HiC. Hi-C phasing is controlled only by HIFIASM_USE_HIC.
# Produces canonical pipeline output: ${dir1}/${base}.hifiasm.primary.fa
# ============================================================================================================================
step1_primary_fasta="${dir1}/${base}.hifiasm.primary.fa"
if is_true "${run_hifiasm}"; then
  info "Step 1: Hifiasm — starting"
  mkdir -p "${dir1}/logs" "${dir1}/checkpoints" "${dir1}/haplotypes"

  require "${ENV_RUNNER}"
  env_run "${ENV_HIFIASM}" hifiasm --version > "${dir1}/logs/hifiasm.version.txt" 2>&1 \
    || fatal "Step 1: hifiasm is unavailable in ${ENV_HIFIASM}"
  require "${samtools_bin}"
  command -v pigz >/dev/null 2>&1 || warn "Step 1: pigz not found; GFA compression will use gzip"

  # Support one or more whitespace-delimited HiFi inputs while preserving the legacy pacbio_hifi variable.
  hifi_input_spec="${PACBIO_HIFI_READS:-${pacbio_hifi}}"
  read -r -a hifi_reads <<< "${hifi_input_spec}"
  (( ${#hifi_reads[@]} > 0 )) || fatal "Step 1: no PacBio HiFi input files configured"
  for read_file in "${hifi_reads[@]}"; do
    [[ -s "${read_file}" ]] || fatal "Step 1: HiFi input missing or empty: ${read_file}"
  done

  if is_true "${HIFIASM_USE_HIC}"; then
    [[ -s "${hic_r1:-}" ]] || fatal "Step 1: HIFIASM_USE_HIC=TRUE but hic_r1 is missing"
    [[ -s "${hic_r2:-}" ]] || fatal "Step 1: HIFIASM_USE_HIC=TRUE but hic_r2 is missing"
  fi

  asm_prefix="${dir1}/${base}.asm"
  if is_true "${HIFIASM_USE_HIC}"; then
    primary_gfa="${asm_prefix}.hic.p_ctg.gfa"
    hap1_gfa="${asm_prefix}.hic.hap1.p_ctg.gfa"
    hap2_gfa="${asm_prefix}.hic.hap2.p_ctg.gfa"
  else
    primary_gfa="${asm_prefix}.bp.p_ctg.gfa"
    hap1_gfa="${asm_prefix}.bp.hap1.p_ctg.gfa"
    hap2_gfa="${asm_prefix}.bp.hap2.p_ctg.gfa"
  fi

  ck_step1="${dir1}/checkpoints/step1.done"
  if [[ -f "${ck_step1}" && -s "${step1_primary_fasta}" ]]; then
    info "Step 1 checkpoint validated — reusing ${step1_primary_fasta}"
  else
    hifiasm_cmd=(hifiasm -o "${asm_prefix}" -t "${nthreads}")

    [[ -n "${HIFIASM_HG_SIZE}" ]] && hifiasm_cmd+=(--hg-size "${HIFIASM_HG_SIZE}")
    [[ -n "${HIFIASM_HOM_COV}" ]] && hifiasm_cmd+=(--hom-cov "${HIFIASM_HOM_COV}")
    [[ -n "${HIFIASM_PURGE_LEVEL}" ]] && hifiasm_cmd+=(-l "${HIFIASM_PURGE_LEVEL}")
    [[ -n "${HIFIASM_PURGE_SIMILARITY}" ]] && hifiasm_cmd+=(-s "${HIFIASM_PURGE_SIMILARITY}")
    [[ -n "${HIFIASM_POST_JOIN}" ]] && hifiasm_cmd+=(-u "${HIFIASM_POST_JOIN}")
    [[ -n "${HIFIASM_CTG_N}" ]] && hifiasm_cmd+=(--ctg-n "${HIFIASM_CTG_N}")
    [[ -n "${HIFIASM_ADAPTER_TRIM}" ]] && hifiasm_cmd+=(-z "${HIFIASM_ADAPTER_TRIM}")
    [[ -n "${HIFIASM_BLOOM_BITS}" ]] && hifiasm_cmd+=(-f "${HIFIASM_BLOOM_BITS}")
    [[ -n "${HIFIASM_CORRECTION_ROUNDS}" ]] && hifiasm_cmd+=(-r "${HIFIASM_CORRECTION_ROUNDS}")
    [[ -n "${HIFIASM_MAX_KOCC}" ]] && hifiasm_cmd+=(--max-kocc "${HIFIASM_MAX_KOCC}")
    [[ -n "${HIFIASM_HIGH_KMER_FACTOR}" ]] && hifiasm_cmd+=(-D "${HIFIASM_HIGH_KMER_FACTOR}")
    [[ -n "${HIFIASM_MAX_OVERLAPS}" ]] && hifiasm_cmd+=(-N "${HIFIASM_MAX_OVERLAPS}")
    is_true "${HIFIASM_PRIMARY_ALT}" && hifiasm_cmd+=(--primary)

    if is_true "${HIFIASM_USE_HIC}"; then
      hifiasm_cmd+=(--h1 "${hic_r1}" --h2 "${hic_r2}")
    fi
    hifiasm_cmd+=("${hifi_reads[@]}")

    {
      echo "# Hifiasm command"
      printf '%q ' "${hifiasm_cmd[@]}"
      echo
    } > "${dir1}/logs/hifiasm.command.sh"

    hifiasm_restore_bin_caches
    env_run "${ENV_HIFIASM}" "${hifiasm_cmd[@]}" \
      2>&1 | tee "${dir1}/logs/hifiasm.${LOG_TS}.log" \
      || fatal "Step 1: Hifiasm failed"

    [[ -s "${primary_gfa}" ]] || fatal "Step 1: expected primary Hifiasm GFA not found: ${primary_gfa}"

    awk '/^S\t/ {print ">"$2"\n"$3}' "${primary_gfa}" > "${step1_primary_fasta}.tmp"
    [[ -s "${step1_primary_fasta}.tmp" ]] || fatal "Step 1: failed to convert primary GFA to FASTA"
    mv -f "${step1_primary_fasta}.tmp" "${step1_primary_fasta}"
    "${samtools_bin}" faidx "${step1_primary_fasta}" || fatal "Step 1: primary FASTA failed samtools faidx validation"

    for hap in 1 2; do
      if [[ "${hap}" == 1 ]]; then gfa="${hap1_gfa}"; else gfa="${hap2_gfa}"; fi
      [[ -s "${gfa}" ]] || continue
      hap_fa="${dir1}/haplotypes/${base}.hifiasm.hap${hap}.fa"
      awk '/^S\t/ {print ">"$2"\n"$3}' "${gfa}" > "${hap_fa}.tmp"
      if [[ -s "${hap_fa}.tmp" ]]; then
        mv -f "${hap_fa}.tmp" "${hap_fa}"
      else
        rm -f "${hap_fa}.tmp"
      fi
    done

    # Hifiasm *.bin files are deliberately retained because they are reusable caches.
    if is_true "${HIFIASM_COMPRESS_GFA}"; then
      shopt -s nullglob
      for gfa in "${asm_prefix}".*.gfa; do
        if command -v pigz >/dev/null 2>&1; then pigz -f -p "${nthreads}" "${gfa}"; else gzip -f "${gfa}"; fi
      done
      shopt -u nullglob
    fi

    {
      echo -e "key\tvalue"
      echo -e "mode\t$([[ ${HIFIASM_USE_HIC} =~ ^(TRUE|True|true|1|YES|yes|Y|y)$ ]] && echo hifi_hic || echo hifi)"
      echo -e "primary_fasta\t${step1_primary_fasta}"
      echo -e "hifiasm_version\t$(tr '\n' ' ' < "${dir1}/logs/hifiasm.version.txt")"
      echo -e "threads\t${nthreads}"
    } > "${dir1}/hifiasm.provenance.tsv"

    touch "${ck_step1}"
  fi
else
  [[ -s "${step1_primary_fasta}" ]] || {
    # Backward-compatible discovery of v0.3.0 outputs when Step 1 is disabled.
    for candidate in "${dir1}/${base}.hic.p_ctg.fa" "${dir1}/${base}.bp.p_ctg.fa" "${dir1}/${base}.p_ctg.fa"; do
      [[ -s "${candidate}" ]] && { step1_primary_fasta="${candidate}"; break; }
    done
  }
fi

if is_true "${HIFIASM_COMPRESS_BIN}" &&
   [[ -f "${dir1}/checkpoints/step1.done" && -s "${dir1}/${base}.hifiasm.primary.fa" ]]; then
  hifiasm_compress_bin_caches
fi

# ============================================================================================================================
# STEP 2: NCBI FCS-Adaptor + FCS-GX — complete adaptor and foreign-contamination cleaning
# Uses the documented FCS action-report workflow; internal FCS-Adaptor ACTION_TRIMs are explicitly resolved.
# Produces canonical pipeline output: ${dir2}/clean.fasta
# ============================================================================================================================
step2_cleaned_fasta="${dir2}/clean.fasta"
if is_true "${run_fcs_early}"; then
  info "Step 2: NCBI FCS-Adaptor + FCS-GX — starting"
  mkdir -p "${dir2}/adaptor" "${dir2}/gx" "${dir2}/logs" "${dir2}/checkpoints"

  [[ -s "${step1_primary_fasta}" ]] || fatal "Step 2: Step 1 primary FASTA missing: ${step1_primary_fasta}"
  require "${python_bin}"
  require "${ADAPTOR_RUNNER}"
  [[ -s "${FCS_RUNNER}" ]] || fatal "Step 2: FCS runner missing or empty: ${FCS_RUNNER}"
  command -v singularity >/dev/null 2>&1 || fatal "Step 2: singularity is required for the configured FCS images"
  [[ -s "${sif_fcs_adaptor}" ]] || fatal "Step 2: FCS-Adaptor image missing: ${sif_fcs_adaptor}"
  [[ -s "${sif_fcs_gx}" ]] || fatal "Step 2: FCS-GX image missing: ${sif_fcs_gx}"
  [[ "${FCS_GX_TAXID}" =~ ^[0-9]+$ ]] || fatal "Step 2: FCS_GX_TAXID must be numeric: ${FCS_GX_TAXID}"
  (( FCS_GX_CORES >= 1 )) || fatal "Step 2: FCS_GX_CORES must be >=1"
  (( FCS_GX_CORES <= 48 )) || warn "Step 2: FCS-GX documentation describes utilization up to 48 cores; configured FCS_GX_CORES=${FCS_GX_CORES}"

  # FCS organism identity is intentionally independent from repeat-annotation taxonomy.
  info "Step 2: assembly TAXID=${TAXID}; FCS-GX TAXID=${FCS_GX_TAXID}; repeat taxonomy is configured separately"

  fcs_input="${dir2}/${base}.fcs.input.fa"
  "${python_bin}" - "${step1_primary_fasta}" "${fcs_input}.tmp" "${FCS_ADAPTOR_MIN_SEQ_LEN}" <<'PY_FCS_VALIDATE' \
    || fatal "Step 2: FCS FASTA validation failed"
import re
import sys
from pathlib import Path

src = Path(sys.argv[1])
dst = Path(sys.argv[2])
min_len = int(sys.argv[3])
max_len = 2_000_000_000
valid_id = re.compile(r"^[A-Za-z0-9._:*-]+$")
valid_seq = re.compile(r"^[ACGTRYSWKMBDHVN.-]+$")
seen = set()
records = 0

with src.open("r", encoding="utf-8") as inp, dst.open("w", encoding="utf-8") as out:
    name = None
    seq_parts = []

    def emit(record_name, parts):
        global records
        if record_name is None:
            return
        seq = "".join(parts).upper()
        if len(seq) < min_len:
            raise SystemExit(f"FCS-invalid FASTA record {record_name}: length {len(seq)} < {min_len}")
        if len(seq) > max_len:
            raise SystemExit(f"FCS-invalid FASTA record {record_name}: length {len(seq)} > {max_len}")
        if not valid_seq.fullmatch(seq):
            raise SystemExit(f"FCS-invalid FASTA sequence alphabet in {record_name}")
        n_fraction = seq.count("N") / len(seq)
        if n_fraction > 0.50:
            raise SystemExit(f"FCS-invalid FASTA record {record_name}: N fraction {n_fraction:.4f} > 0.50")
        out.write(f">{record_name}\n")
        for i in range(0, len(seq), 80):
            out.write(seq[i:i+80] + "\n")
        records += 1

    for raw in inp:
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            emit(name, seq_parts)
            header = line[1:].strip()
            record_id = header.split()[0] if header else ""
            if not record_id or len(record_id) >= 50 or not valid_id.fullmatch(record_id):
                raise SystemExit(f"FCS-invalid FASTA identifier: {record_id!r}")
            if record_id in seen:
                raise SystemExit(f"Duplicate FASTA identifier: {record_id}")
            seen.add(record_id)
            name = record_id
            seq_parts = []
        else:
            if name is None:
                raise SystemExit("FASTA sequence encountered before the first header")
            seq_parts.append(line)
    emit(name, seq_parts)

if records == 0:
    raise SystemExit("No FASTA records found")
PY_FCS_VALIDATE
  mv -f "${fcs_input}.tmp" "${fcs_input}"
  [[ -s "${fcs_input}" ]] || fatal "Step 2: validated FCS input is empty"

  adaptor_report="${dir2}/adaptor/fcs_adaptor_report.txt"
  adaptor_clean="${dir2}/adaptor/clean.fasta"
  adaptor_contam="${dir2}/adaptor/contam.fasta"
  ck_adaptor="${dir2}/checkpoints/adaptor.done"
  if [[ ! -f "${ck_adaptor}" || ! -s "${adaptor_report}" || ! -s "${adaptor_clean}" ]]; then
    adaptor_tax_flag="--euk"
    [[ "${FCS_ADAPTOR_TAXONOMY}" == "prok" ]] && adaptor_tax_flag="--prok"
    [[ "${FCS_ADAPTOR_TAXONOMY}" =~ ^(euk|prok)$ ]] || fatal "Step 2: FCS_ADAPTOR_TAXONOMY must be euk or prok"

    adaptor_cmd=("${ADAPTOR_RUNNER}" --fasta-input "${fcs_input}" --output-dir "${dir2}/adaptor" "${adaptor_tax_flag}" \
      --container-engine singularity --image "${sif_fcs_adaptor}")
    is_true "${FCS_DEBUG}" && adaptor_cmd+=(--debug)
    "${adaptor_cmd[@]}" \
      2>&1 | tee "${dir2}/logs/fcs_adaptor.${LOG_TS}.log" \
      || fatal "Step 2: FCS-Adaptor failed"
    [[ -s "${adaptor_report}" ]] || fatal "Step 2: FCS-Adaptor action report missing: ${adaptor_report}"

    adaptor_actions="${adaptor_report}"
    if [[ "${FCS_ADAPTOR_INTERNAL_ACTION}" == "mask" ]]; then
      adaptor_actions="${dir2}/adaptor/fcs_adaptor_report.mask-internal.txt"
      awk -F'\t' 'BEGIN{OFS="\t"} /^#/ {print; next} {
        if ($3=="ACTION_TRIM") { split($4,r,/\.\./); if ((r[1]+0)>1 && (r[2]+0)<($2+0)) $3="FIX" }
        print
      }' "${adaptor_report}" > "${adaptor_actions}"
    elif [[ "${FCS_ADAPTOR_INTERNAL_ACTION}" != "split" ]]; then
      fatal "Step 2: FCS_ADAPTOR_INTERNAL_ACTION must be split or mask"
    fi

    FCS_DEFAULT_IMAGE="${sif_fcs_gx}" "${python_bin}" "${FCS_RUNNER}" clean genome \
      --action-report "${adaptor_actions}" \
      --output "${adaptor_clean}" \
      --contam-fasta-out "${adaptor_contam}" \
      < "${fcs_input}" \
      2>&1 | tee "${dir2}/logs/fcs_adaptor.clean.${LOG_TS}.log" \
      || fatal "Step 2: applying FCS-Adaptor actions failed"
    [[ -s "${adaptor_clean}" ]] || fatal "Step 2: FCS-Adaptor clean FASTA missing"
    touch "${ck_adaptor}"
  fi

  # --------------------------------------------------------------------------------------------------------------
  # FCS-GX database lifecycle: integrity and freshness are deliberately separate checks.
  # NCBI distributes GX through a manifest-driven synchronizer.  We first validate the installed release against
  # its own local all.manifest, then (unless offline) retrieve only the configured remote manifest and compare it.
  # Full database synchronization is performed only when policy=update/download.
  # --------------------------------------------------------------------------------------------------------------
  mkdir -p "${FCS_GX_DB}"
  gx_local_manifest="${FCS_GX_DB}/all.manifest"
  gx_remote_manifest="${dir2}/gx_database_remote.manifest"
  gx_status_file="${dir2}/gx_database_status.tsv"
  gx_update_lock="${FCS_GX_DB}.genoflow-update.lock"
  gx_policy="${FCS_GX_DB_POLICY,,}"
  [[ "${gx_policy}" == "validate" ]] && gx_policy="check"  # backward compatibility with v0.4.0-draft1
  gx_integrity_mode="FULL"
  is_true "${FCS_GX_DB_CHECK}" || gx_integrity_mode="SKIPPED"

  case "${gx_policy}" in
    check|update|require_latest|offline|download) ;;
    *) fatal "Step 2: FCS_GX_DB_POLICY must be check, update, require_latest, offline, download, or validate (alias for check)" ;;
  esac

  fcs_gx_required_files=(
    all.assemblies.tsv all.blast_div.tsv.gz all.gxi all.gxs all.manifest
    all.meta.jsonl all.README.txt all.seq_info.tsv.gz all.taxa.tsv
  )

  fcs_gx_db_complete() {
    local db_file
    for db_file in "${fcs_gx_required_files[@]}"; do
      [[ -s "${FCS_GX_DB}/${db_file}" ]] || return 1
    done
    return 0
  }

  fcs_gx_validate_local() {
    fcs_gx_db_complete || return 1
    if is_true "${FCS_GX_DB_CHECK}"; then
      "${python_bin}" "${FCS_RUNNER}" db check \
        --mft "${gx_local_manifest}" \
        --dir "${FCS_GX_DB}" \
        > "${dir2}/logs/fcs_gx.db_check.${LOG_TS}.log" 2>&1
    else
      info "Step 2: skipping full FCS-GX database integrity scan (FCS_GX_DB_CHECK=FALSE)"
      return 0
    fi
  }

  fcs_gx_fetch_remote_manifest() {
    rm -f "${gx_remote_manifest}"
    if [[ "${FCS_GX_MANIFEST}" =~ ^https?:// ]]; then
      command -v curl >/dev/null 2>&1 || fatal "Step 2: curl is required to check the remote FCS-GX manifest"
      curl --fail --location --silent --show-error \
        --retry 3 --retry-delay 5 \
        --output "${gx_remote_manifest}.tmp" \
        "${FCS_GX_MANIFEST}" \
        || { rm -f "${gx_remote_manifest}.tmp"; return 1; }
      mv -f "${gx_remote_manifest}.tmp" "${gx_remote_manifest}"
    elif [[ -s "${FCS_GX_MANIFEST}" ]]; then
      cp -f "${FCS_GX_MANIFEST}" "${gx_remote_manifest}"
    else
      return 1
    fi
    [[ -s "${gx_remote_manifest}" ]]
  }

  gx_local_state="MISSING_OR_INCOMPLETE"
  gx_freshness="NOT_CHECKED"
  gx_action="NONE"
  gx_local_sha="NA"
  gx_remote_sha="NA"

  if fcs_gx_db_complete; then
    if fcs_gx_validate_local; then
      gx_local_state="VALID"
    else
      gx_local_state="INVALID"
    fi
  fi

  # Fetch only the small current manifest when network access is allowed.  A transient freshness-check failure does
  # not invalidate an otherwise sound local DB under policy=check, but update/require_latest/download require it.
  if [[ "${gx_policy}" != "offline" ]]; then
    if fcs_gx_fetch_remote_manifest; then
      gx_remote_sha="$(sha256sum "${gx_remote_manifest}" | awk '{print $1}')"
      if [[ -s "${gx_local_manifest}" ]]; then
        gx_local_sha="$(sha256sum "${gx_local_manifest}" | awk '{print $1}')"
        if cmp -s "${gx_local_manifest}" "${gx_remote_manifest}"; then
          gx_freshness="UP_TO_DATE"
        else
          gx_freshness="UPDATE_AVAILABLE"
        fi
      else
        gx_freshness="LOCAL_MANIFEST_MISSING"
      fi
    else
      gx_freshness="REMOTE_CHECK_FAILED"
      if [[ "${gx_policy}" == "check" ]]; then
        warn "Step 2: could not retrieve the configured FCS-GX manifest; continuing only if the local database validates"
      else
        fatal "Step 2: unable to retrieve FCS-GX manifest required by policy=${gx_policy}: ${FCS_GX_MANIFEST}"
      fi
    fi
  else
    gx_freshness="OFFLINE_NOT_CHECKED"
  fi

  case "${gx_policy}" in
    offline)
      [[ "${gx_local_state}" == "VALID" ]] \
        || fatal "Step 2: offline policy requires a complete FCS-GX database that passes its local manifest integrity check: ${FCS_GX_DB}"
      gx_action="USE_VALIDATED_LOCAL"
      ;;

    check)
      [[ "${gx_local_state}" == "VALID" ]] \
        || fatal "Step 2: FCS-GX database is missing/incomplete/corrupt. Use FCS_GX_DB_POLICY=download or update to synchronize it: ${FCS_GX_DB}"
      if [[ "${gx_freshness}" == "UPDATE_AVAILABLE" ]]; then
        warn "Step 2: a newer/different FCS-GX database manifest is available; policy=check keeps the validated local release"
        gx_action="USE_VALIDATED_LOCAL_UPDATE_AVAILABLE"
      elif [[ "${gx_freshness}" == "UP_TO_DATE" ]]; then
        if is_true "${FCS_GX_DB_CHECK}"; then
          info "Step 2: FCS-GX database is valid and up to date"
        else
          info "Step 2: GX required files are present and the manifest is current; full integrity scan skipped"
        fi
        gx_action="USE_UP_TO_DATE_LOCAL"
      else
        gx_action="USE_VALIDATED_LOCAL_FRESHNESS_UNKNOWN"
      fi
      ;;

    require_latest)
      [[ "${gx_local_state}" == "VALID" ]] \
        || fatal "Step 2: require_latest policy requires a valid installed FCS-GX database: ${FCS_GX_DB}"
      [[ "${gx_freshness}" == "UP_TO_DATE" ]] \
        || fatal "Step 2: installed FCS-GX database is not identical to the configured current manifest (status=${gx_freshness}); use FCS_GX_DB_POLICY=update"
      gx_action="USE_UP_TO_DATE_LOCAL"
      ;;

    update|download)
      # download and update both use NCBI's manifest-driven synchronizer.  Existing valid files are reused, so this is
      # incremental rather than an unconditional ~470-GiB redownload.  A genoflow-specific lock prevents concurrent
      # pipeline runs from updating the same persistent database simultaneously.
      if [[ "${gx_policy}" == "download" && "${gx_local_state}" == "VALID" && "${gx_freshness}" == "UP_TO_DATE" ]]; then
        info "Step 2: FCS-GX database already exists, validates, and is current; no download required"
        gx_action="USE_UP_TO_DATE_LOCAL"
      elif [[ "${gx_policy}" == "update" && "${gx_local_state}" == "VALID" && "${gx_freshness}" == "UP_TO_DATE" ]]; then
        info "Step 2: FCS-GX database already validates against the current manifest; no update required"
        gx_action="USE_UP_TO_DATE_LOCAL"
      else
        if ! mkdir "${gx_update_lock}" 2>/dev/null; then
          fatal "Step 2: another process may be updating the FCS-GX database (lock exists: ${gx_update_lock})"
        fi
        owned_lock_dirs+=("${gx_update_lock}")
        if [[ -s "${gx_local_manifest}" ]]; then
          cp -f "${gx_local_manifest}" "${dir2}/gx_database_manifest.previous"
        fi

        info "Step 2: synchronizing FCS-GX database using policy=${gx_policy}"
        if ! "${python_bin}" "${FCS_RUNNER}" db get \
          --mft "${FCS_GX_MANIFEST}" \
          --dir "${FCS_GX_DB}" \
          2>&1 | tee "${dir2}/logs/fcs_gx.db_get.${LOG_TS}.log"; then
          release_lock_dir "${gx_update_lock}"
          fatal "Step 2: FCS-GX database synchronization failed"
        fi

        release_lock_dir "${gx_update_lock}"

        fcs_gx_validate_local \
          || fatal "Step 2: FCS-GX database failed integrity validation after synchronization"
        gx_local_state="VALID"
        gx_local_sha="$(sha256sum "${gx_local_manifest}" | awk '{print $1}')"

        # The remote manifest was fetched immediately before synchronization.  Confirm that the synchronized local
        # database release now corresponds to it; this also protects against incomplete/interrupted updates.
        if [[ -s "${gx_remote_manifest}" ]] && cmp -s "${gx_local_manifest}" "${gx_remote_manifest}"; then
          gx_freshness="UP_TO_DATE"
        elif [[ -s "${gx_remote_manifest}" ]]; then
          fatal "Step 2: FCS-GX synchronization completed but the installed all.manifest does not match the requested manifest"
        fi
        gx_action="SYNCHRONIZED"
      fi
      ;;
  esac

  # Revalidate the final database state before screening regardless of policy.
  fcs_gx_validate_local \
    || fatal "Step 2: final FCS-GX database integrity check failed: ${FCS_GX_DB}"
  gx_local_sha="$(sha256sum "${gx_local_manifest}" | awk '{print $1}')"

  cp -f "${gx_local_manifest}" "${dir2}/gx_database_manifest.used"
  sha256sum "${gx_local_manifest}" > "${dir2}/gx_database_manifest.sha256"
  {
    printf 'key\tvalue\n'
    printf 'database_path\t%s\n' "${FCS_GX_DB}"
    printf 'policy\t%s\n' "${gx_policy}"
    printf 'db_check_enabled\t%s\n' "${FCS_GX_DB_CHECK}"
    printf 'integrity_check\t%s\n' "${gx_integrity_mode}"
    printf 'local_integrity\t%s\n' "${gx_local_state}"
    printf 'freshness\t%s\n' "${gx_freshness}"
    printf 'action\t%s\n' "${gx_action}"
    printf 'local_manifest_sha256\t%s\n' "${gx_local_sha}"
    printf 'remote_manifest_sha256\t%s\n' "${gx_remote_sha}"
    printf 'manifest_source\t%s\n' "${FCS_GX_MANIFEST}"
    printf 'checked_at\t%s\n' "$(date -Is)"
  } > "${gx_status_file}"

  info "Step 2: FCS-GX database status: integrity=${gx_local_state}; freshness=${gx_freshness}; action=${gx_action}"

  gx_env="${dir2}/gx/fcs_gx.env"
  {
    echo "GX_NUM_CORES=${FCS_GX_CORES}"
    [[ -n "${FCS_GX_ALIGN_EXCLUDE_TAXA}" ]] && echo "GX_ALIGN_EXCLUDE_TAXA=${FCS_GX_ALIGN_EXCLUDE_TAXA}"
    [[ -n "${FCS_GX_EXTRA_CONTAM_DIVS}" ]] && echo "GX_EXTRA_CONTAM_DIVS=${FCS_GX_EXTRA_CONTAM_DIVS}"
  } > "${gx_env}"

  gx_ck="${dir2}/checkpoints/gx_screen.done"
  if [[ ! -f "${gx_ck}" ]]; then
    FCS_DEFAULT_IMAGE="${sif_fcs_gx}" "${python_bin}" "${FCS_RUNNER}" --env-file "${gx_env}" screen genome \
      --fasta "${adaptor_clean}" \
      --out-dir "${dir2}/gx" \
      --gx-db "${FCS_GX_DB}" \
      --tax-id "${FCS_GX_TAXID}" \
      2>&1 | tee "${dir2}/logs/fcs_gx.screen.${LOG_TS}.log" \
      || fatal "Step 2: FCS-GX screen failed"
    touch "${gx_ck}"
  fi

  mapfile -t gx_reports < <(find "${dir2}/gx" -maxdepth 1 -type f -name '*.fcs_gx_report.txt' -size +0c | sort)
  (( ${#gx_reports[@]} == 1 )) || fatal "Step 2: expected exactly one FCS-GX action report, found ${#gx_reports[@]}"
  gx_report="${gx_reports[0]}"
  gx_actions="${gx_report}"
  if [[ "${FCS_GX_INTERNAL_ACTION}" == "split" ]]; then
    gx_actions="${dir2}/gx/fcs_gx_report.split-internal.txt"
    sed 's/\tFIX\t/\tSPLIT\t/g' "${gx_report}" > "${gx_actions}"
  elif [[ "${FCS_GX_INTERNAL_ACTION}" != "fix" ]]; then
    fatal "Step 2: FCS_GX_INTERNAL_ACTION must be fix or split"
  fi

  gx_contam="${dir2}/gx/contam.fasta"
  FCS_DEFAULT_IMAGE="${sif_fcs_gx}" "${python_bin}" "${FCS_RUNNER}" clean genome \
    --action-report "${gx_actions}" \
    --output "${step2_cleaned_fasta}.tmp" \
    --contam-fasta-out "${gx_contam}" \
    < "${adaptor_clean}" \
    2>&1 | tee "${dir2}/logs/fcs_gx.clean.${LOG_TS}.log" \
    || fatal "Step 2: applying FCS-GX actions failed"
  [[ -s "${step2_cleaned_fasta}.tmp" ]] || fatal "Step 2: cleaned assembly is empty"
  mv -f "${step2_cleaned_fasta}.tmp" "${step2_cleaned_fasta}"
  "${samtools_bin}" faidx "${step2_cleaned_fasta}" || fatal "Step 2: cleaned FASTA failed validation"
  touch "${dir2}/checkpoints/step2.done"
  info "Step 2: complete → ${step2_cleaned_fasta}"
else
  step2_cleaned_fasta="${step1_primary_fasta}"
  info "Step 2: disabled — carrying Step 1 assembly forward unchanged"
fi

# The validated cleaned FASTA is the Step 2 output. The FCS input is rebuilt
# from Step 1 whenever FCS must run again.
if is_true "${CLEANUP_COMPLETED_STEPS}" &&
   [[ -f "${dir2}/checkpoints/step2.done" && -s "${dir2}/clean.fasta" ]]; then
  rm -f -- "${dir2}/${base}.fcs.input.fa"
fi

# ============================================================================================================================
# STEP 3: Three-round maximum PacBio HiFi polishing — minimap2 + SAMtools + DeepVariant + BCFtools
# DeepVariant corrections are restricted to PASS homozygous-alt genotypes (GT=AA); rounds can stop early at convergence.
# Produces canonical pipeline output: ${dir3}/${base}.prepolish.final.fa
# ============================================================================================================================
step3_polished_fasta="${dir3}/${base}.prepolish.final.fa"
if is_true "${run_deepvariant}"; then
  info "Step 3: DeepVariant HiFi polishing — starting"
  mkdir -p "${dir3}/logs" "${dir3}/checkpoints" "${dir3}/metrics"
  [[ -s "${step2_cleaned_fasta}" ]] || fatal "Step 3: input assembly missing: ${step2_cleaned_fasta}"
  require "${minimap2_bin}"
  require "${samtools_bin}"
  require "${bcftools_bin}"
  command -v singularity >/dev/null 2>&1 || fatal "Step 3: singularity is required for DeepVariant"
  [[ -s "${sif_deepvariant}" ]] || fatal "Step 3: DeepVariant image missing: ${sif_deepvariant}"

  "${minimap2_bin}" --version > "${dir3}/logs/minimap2.version.txt" 2>&1 || fatal "Step 3: minimap2 version check failed"
  "${samtools_bin}" --version > "${dir3}/logs/samtools.version.txt" 2>&1 || fatal "Step 3: samtools version check failed"
  "${bcftools_bin}" --version > "${dir3}/logs/bcftools.version.txt" 2>&1 || fatal "Step 3: bcftools version check failed"

  hifi_input_spec="${PACBIO_HIFI_READS:-${pacbio_hifi}}"
  read -r -a hifi_reads <<< "${hifi_input_spec}"
  for read_file in "${hifi_reads[@]}"; do [[ -s "${read_file}" ]] || fatal "Step 3: HiFi input missing: ${read_file}"; done


  asm_in="${step2_cleaned_fasta}"
  completed_round=0
  for ((round=1; round<=DEEPVARIANT_PRE_ROUNDS; round++)); do
    prefix="${dir3}/${base}.round${round}"
    sorted_bam="${prefix}.sorted.bam"
    dv_vcf="${prefix}.dv.vcf.gz"
    correction_vcf="${prefix}.corrections.vcf.gz"
    asm_out="${prefix}.polished.fa"
    chain_file="${prefix}.consensus.chain"
    metrics_file="${dir3}/metrics/round${round}.corrections.tsv"
    round_ck="${dir3}/checkpoints/round${round}.done"

    [[ -s "${asm_in}.fai" ]] || "${samtools_bin}" faidx "${asm_in}" || fatal "Round ${round}: faidx failed"

    if [[ -f "${round_ck}" && -s "${asm_out}" ]]; then
      info "Round ${round}: validated checkpoint found"
      completed_round="${round}"
      asm_in="${asm_out}"
      continue
    fi

    mm2_args=("${minimap2_bin}" -a -x map-hifi --secondary=no --cs=long -Y -K "${MINIMAP2_BATCH_SIZE}" -t "${nthreads}" \
      -R "@RG\\tID:${base}.HiFi\\tSM:${base}\\tPL:PACBIO\\tLB:SMRTbell")
    mm2_args+=("${asm_in}" "${hifi_reads[@]}")

    info "Round ${round}: mapping → DeepVariant-compatible BAM → coordinate sort + BAI"
    "${mm2_args[@]}" \
      | "${samtools_bin}" view -u -F "${DEEPVARIANT_EXCLUDE_FLAGS}" -e "${DEEPVARIANT_VIEW_EXPR}" - \
      | "${samtools_bin}" sort --write-index -@ "${nthreads}" -m "${SAMTOOLS_SORT_MEM}" -l "${SAMTOOLS_BAM_COMPRESSION}" \
          -T "${dir3}/round${round}.coord.tmp" -O BAM \
          -o "${sorted_bam}##idx##${sorted_bam}.bai" - \
      || fatal "Round ${round}: minimap2/SAMtools DeepVariant BAM preparation failed"
    bam_for_dv="${sorted_bam}"
    "${samtools_bin}" quickcheck -v "${bam_for_dv}" || fatal "Round ${round}: BAM quickcheck failed"
    [[ -s "${bam_for_dv}.bai" ]] || fatal "Round ${round}: DeepVariant BAI index missing"
    if [[ "$("${samtools_bin}" view -c -f "${DEEPVARIANT_EXCLUDE_FLAGS}" "${bam_for_dv}")" -ne 0 ]]; then
      fatal "Round ${round}: DeepVariant BAM still contains secondary/supplementary alignments"
    fi
    if [[ "$("${samtools_bin}" view -c -e 'length(seq) == 0 || length(qual) != length(seq)' "${bam_for_dv}")" -ne 0 ]]; then
      fatal "Round ${round}: DeepVariant BAM still contains records without base qualities"
    fi

    dv_tmp=$(mktemp -d -p "${dir3}" "dv.round${round}.XXXXXXXX") || fatal "Round ${round}: mktemp failed"
    singularity exec --cleanenv -B "${working_dir}:${working_dir}" -B "${tmpdir}:${tmpdir}" "${sif_deepvariant}" \
      /opt/deepvariant/bin/run_deepvariant \
        --model_type=PACBIO \
        --ref="${asm_in}" \
        --reads="${bam_for_dv}" \
        --output_vcf="${dv_vcf}" \
        --num_shards="${DEEPVARIANT_SHARDS}" \
        --postprocess_cpus="${DEEPVARIANT_POSTPROCESS_CPUS}" \
        --intermediate_results_dir="${dv_tmp}" \
        --make_examples_extra_args="${DEEPVARIANT_MAKE_EXAMPLES_EXTRA_ARGS}" \
        --vcf_stats_report="${DEEPVARIANT_VCF_STATS}" \
      2>&1 | tee "${dir3}/logs/round${round}.deepvariant.${LOG_TS}.log" \
      || { rm -rf "${dv_tmp}"; fatal "Round ${round}: DeepVariant failed"; }
    rm -rf "${dv_tmp}"
    [[ -s "${dv_vcf}" ]] || fatal "Round ${round}: DeepVariant VCF missing"
    "${bcftools_bin}" view -h "${dv_vcf}" >/dev/null || fatal "Round ${round}: DeepVariant VCF validation failed"

    mapfile -t dv_samples < <("${bcftools_bin}" query -l "${dv_vcf}")
    (( ${#dv_samples[@]} == 1 )) || fatal "Round ${round}: expected exactly one DeepVariant sample, found ${#dv_samples[@]}"
    dv_sample="${dv_samples[0]}"

    "${bcftools_bin}" view --threads "${nthreads}" -f PASS -i 'GT="AA"' -Oz -o "${correction_vcf}" --write-index "${dv_vcf}" \
      || fatal "Round ${round}: filtering PASS homozygous-alt corrections failed"
    corrections=$("${bcftools_bin}" view -H "${correction_vcf}" | wc -l | awk '{print $1}')
    snps=$("${bcftools_bin}" view -H -v snps "${correction_vcf}" | wc -l | awk '{print $1}')
    indels=$("${bcftools_bin}" view -H -v indels "${correction_vcf}" | wc -l | awk '{print $1}')
    {
      echo -e "metric\tvalue"
      echo -e "corrections_total\t${corrections}"
      echo -e "snps\t${snps}"
      echo -e "indels\t${indels}"
      echo -e "sample\t${dv_sample}"
    } > "${metrics_file}"

    if (( corrections == 0 )); then
      cp -f "${asm_in}" "${asm_out}"
      info "Round ${round}: zero applicable corrections; polishing converged"
    else
      consensus_args=("${bcftools_bin}" consensus --fasta-ref "${asm_in}" --samples "${dv_sample}" --haplotype A --output "${asm_out}.tmp")
      is_true "${DEEPVARIANT_WRITE_CHAIN}" && consensus_args+=(--chain "${chain_file}")
      consensus_args+=("${correction_vcf}")
      "${consensus_args[@]}" || fatal "Round ${round}: bcftools consensus failed"
      mv -f "${asm_out}.tmp" "${asm_out}"
    fi
    [[ -s "${asm_out}" ]] || fatal "Round ${round}: polished FASTA missing"
    "${samtools_bin}" faidx "${asm_out}" || fatal "Round ${round}: polished FASTA validation failed"

    touch "${round_ck}"
    completed_round="${round}"
    asm_in="${asm_out}"

    if is_true "${DEEPVARIANT_CLEAN_INTERMEDIATES}"; then
      rm -f "${sorted_bam}" "${sorted_bam}.bai" "${sorted_bam}.csi"
    fi
    if is_true "${DEEPVARIANT_STOP_ON_ZERO}" && (( corrections == 0 )); then break; fi
  done

  (( completed_round > 0 )) || fatal "Step 3: no polishing round completed"
  cp -f "${asm_in}" "${step3_polished_fasta}.tmp"
  mv -f "${step3_polished_fasta}.tmp" "${step3_polished_fasta}"
  "${samtools_bin}" faidx "${step3_polished_fasta}" || fatal "Step 3: final polished FASTA failed validation"
  printf 'completed_round\t%s\nfinal_fasta\t%s\n' "${completed_round}" "${step3_polished_fasta}" > "${dir3}/polishing.summary.tsv"
  touch "${dir3}/checkpoints/step3.done"
  info "Step 3: complete after ${completed_round} round(s) → ${step3_polished_fasta}"
else
  step3_polished_fasta="${step2_cleaned_fasta}"
  info "Step 3: disabled — carrying Step 2 assembly forward unchanged"
fi

# ============================================================================================================================
# STEP 4: Inspector — enhanced Python 3 reference-free HiFi evaluation, correction, and post-correction validation
# genoflowHiFi expects the enhanced inspector.py / inspector-correct.py replacements inside ENV_INSPECTOR.
# Short contigs excluded by Inspector's evaluation threshold are carried forward unchanged.
# Produces canonical pipeline output: ${dir4}/${base}.inspector.corrected.fa
# ============================================================================================================================
step4_output_fasta="${dir4}/${base}.inspector.corrected.fa"
if is_true "${run_inspector}"; then
  info "Step 4: Inspector enhanced evaluation/correction — starting"
  eval_dir="${dir4}/evaluation"
  corr_dir="${dir4}/correction"
  valid_dir="${dir4}/validation"
  mkdir -p "${eval_dir}" "${corr_dir}" "${valid_dir}" "${dir4}/logs" "${dir4}/checkpoints"
  [[ -s "${step3_polished_fasta}" ]] || fatal "Step 4: input assembly missing: ${step3_polished_fasta}"

  # Verify that the environment contains the enhanced wrappers, not merely upstream Inspector.
  env_run "${ENV_INSPECTOR}" inspector.py --version > "${dir4}/logs/inspector.version.txt" 2>&1 \
    || fatal "Step 4: inspector.py unavailable in ${ENV_INSPECTOR}"
  inspector_help="$(env_run "${ENV_INSPECTOR}" inspector.py --help 2>&1)" \
    || fatal "Step 4: inspector.py help failed in ${ENV_INSPECTOR}"
  grep -Fq -- '--strict_commands' <<< "${inspector_help}" \
    || fatal "Step 4: inspector.py does not expose --strict_commands; install/replace with the enhanced Python 3 version"
  inspector_correct_help="$(env_run "${ENV_INSPECTOR}" inspector-correct.py --help 2>&1)" \
    || fatal "Step 4: inspector-correct.py help failed in ${ENV_INSPECTOR}"
  grep -Fq -- '--strict' <<< "${inspector_correct_help}" \
    || fatal "Step 4: inspector-correct.py does not expose --strict; install/replace with the enhanced Python 3 version"

  hifi_input_spec="${PACBIO_HIFI_READS:-${pacbio_hifi}}"
  read -r -a hifi_reads <<< "${hifi_input_spec}"
  for read_file in "${hifi_reads[@]}"; do [[ -s "${read_file}" ]] || fatal "Step 4: HiFi input missing: ${read_file}"; done

  eval_ck="${dir4}/checkpoints/evaluation.done"
  if [[ ! -f "${eval_ck}" || ! -s "${eval_dir}/summary_statistics" || ! -s "${eval_dir}/valid_contig.fa" ]]; then
    inspector_eval=(inspector.py --contig "${step3_polished_fasta}" --read "${hifi_reads[@]}" --outpath "${eval_dir}" \
      --datatype hifi --thread "${INSPECTOR_THREADS}" --strict_commands --strict_workers \
      --minimap2_read_preset "${INSPECTOR_MINIMAP2_PRESET}" --minimap2_hpc "${INSPECTOR_MINIMAP2_HPC}" \
      --minimap2_index_batch "${INSPECTOR_MINIMAP2_INDEX_BATCH}" --samtools_sort_level "${INSPECTOR_BAM_COMPRESSION}")
    [[ -n "${INSPECTOR_SAMTOOLS_SORT_MEMORY}" ]] && inspector_eval+=(--samtools_sort_memory "${INSPECTOR_SAMTOOLS_SORT_MEMORY}")
    [[ -n "${INSPECTOR_MIN_CONTIG_LENGTH}" ]] && inspector_eval+=(--min_contig_length "${INSPECTOR_MIN_CONTIG_LENGTH}")
    [[ -n "${INSPECTOR_MIN_STRUCTURAL_CONTIG}" ]] && inspector_eval+=(--min_contig_length_assemblyerror "${INSPECTOR_MIN_STRUCTURAL_CONTIG}")
    is_false "${INSPECTOR_EVALUATE_STRUCTURAL}" && inspector_eval+=(--skip_structural_error)
    is_false "${INSPECTOR_EVALUATE_BASE}" && inspector_eval+=(--skip_base_error)

    env_run "${ENV_INSPECTOR}" "${inspector_eval[@]}" \
      2>&1 | tee "${dir4}/logs/inspector.evaluation.${LOG_TS}.log" \
      || fatal "Step 4: Inspector evaluation failed"
    [[ -s "${eval_dir}/summary_statistics" ]] || fatal "Step 4: Inspector summary_statistics missing"
    [[ -s "${eval_dir}/valid_contig.fa" ]] || fatal "Step 4: Inspector valid_contig.fa missing"
    [[ -s "${eval_dir}/read_to_contig.bam" ]] || fatal "Step 4: Inspector read_to_contig.bam missing"
    touch "${eval_ck}"
  fi

  if is_true "${INSPECTOR_CORRECT}"; then
    corr_ck="${dir4}/checkpoints/correction.done"
    if [[ ! -f "${corr_ck}" || ! -s "${corr_dir}/contig_corrected.fa" ]]; then
      inspector_corr=(inspector-correct.py --inspector "${eval_dir}" --outpath "${corr_dir}" --datatype pacbio-hifi \
        --thread "${INSPECTOR_THREADS}" --jobs "${INSPECTOR_CORRECTION_JOBS}" --flyetimeout "${INSPECTOR_FLYE_TIMEOUT}" --strict)
      is_false "${INSPECTOR_CORRECT_STRUCTURAL}" && inspector_corr+=(--skip_structural)
      is_false "${INSPECTOR_CORRECT_BASE}" && inspector_corr+=(--skip_baseerror)
      is_true "${INSPECTOR_KEEP_INTERMEDIATES}" && inspector_corr+=(--keep_intermediates)

      env_run "${ENV_INSPECTOR}" "${inspector_corr[@]}" \
        2>&1 | tee "${dir4}/logs/inspector.correction.${LOG_TS}.log" \
        || fatal "Step 4: Inspector correction failed"
      [[ -s "${corr_dir}/contig_corrected.fa" ]] || fatal "Step 4: corrected Inspector FASTA missing"
      touch "${corr_ck}"
    fi
    corrected_valid="${corr_dir}/contig_corrected.fa"
  else
    corrected_valid="${eval_dir}/valid_contig.fa"
  fi

  # Reconstruct the full assembly: corrected Inspector-valid contigs + all original contigs excluded by Inspector's
  # general evaluation threshold. This keeps Inspector's documented default threshold without silently deleting sequence.
  "${python_bin}" - "${step3_polished_fasta}" "${eval_dir}/valid_contig.fa" "${corrected_valid}" "${step4_output_fasta}.tmp" <<'PY'
import sys
from pathlib import Path

def fasta(path):
    out=[]; name=None; seq=[]
    with open(path) as fh:
        for raw in fh:
            line=raw.rstrip("\n\r")
            if not line: continue
            if line.startswith(">"):
                if name is not None: out.append((name,"".join(seq)))
                name=line[1:].split()[0]; seq=[]
            else: seq.append(line.strip())
    if name is not None: out.append((name,"".join(seq)))
    return out

original, valid, corrected, output = map(Path, sys.argv[1:])
orig = fasta(original); valid_ids={n for n,_ in fasta(valid)}; corr=fasta(corrected)
if not corr: raise SystemExit("corrected/valid Inspector FASTA is empty")
with open(output,"w") as out:
    for n,s in corr: out.write(f">{n}\n{s}\n")
    for n,s in orig:
        if n not in valid_ids: out.write(f">{n}\n{s}\n")
PY
  [[ -s "${step4_output_fasta}.tmp" ]] || fatal "Step 4: failed to reconstruct full Inspector output assembly"
  mv -f "${step4_output_fasta}.tmp" "${step4_output_fasta}"
  "${samtools_bin}" faidx "${step4_output_fasta}" || fatal "Step 4: reconstructed Inspector FASTA failed validation"

  if is_true "${INSPECTOR_VALIDATE_CORRECTION}"; then
    validation_ck="${dir4}/checkpoints/validation.done"
    if [[ ! -f "${validation_ck}" || ! -s "${valid_dir}/summary_statistics" ]]; then
      inspector_validate=(inspector.py --contig "${step4_output_fasta}" --read "${hifi_reads[@]}" --outpath "${valid_dir}" \
        --datatype hifi --thread "${INSPECTOR_THREADS}" --strict_commands --strict_workers \
        --minimap2_read_preset "${INSPECTOR_MINIMAP2_PRESET}" --minimap2_hpc "${INSPECTOR_MINIMAP2_HPC}" \
        --minimap2_index_batch "${INSPECTOR_MINIMAP2_INDEX_BATCH}" --samtools_sort_level "${INSPECTOR_BAM_COMPRESSION}")
      [[ -n "${INSPECTOR_SAMTOOLS_SORT_MEMORY}" ]] && inspector_validate+=(--samtools_sort_memory "${INSPECTOR_SAMTOOLS_SORT_MEMORY}")
      [[ -n "${INSPECTOR_MIN_CONTIG_LENGTH}" ]] && inspector_validate+=(--min_contig_length "${INSPECTOR_MIN_CONTIG_LENGTH}")
      [[ -n "${INSPECTOR_MIN_STRUCTURAL_CONTIG}" ]] && inspector_validate+=(--min_contig_length_assemblyerror "${INSPECTOR_MIN_STRUCTURAL_CONTIG}")
      is_false "${INSPECTOR_EVALUATE_STRUCTURAL}" && inspector_validate+=(--skip_structural_error)
      is_false "${INSPECTOR_EVALUATE_BASE}" && inspector_validate+=(--skip_base_error)

      env_run "${ENV_INSPECTOR}" "${inspector_validate[@]}" \
        2>&1 | tee "${dir4}/logs/inspector.validation.${LOG_TS}.log" \
        || fatal "Step 4: post-correction Inspector validation failed"
      [[ -s "${valid_dir}/summary_statistics" ]] || fatal "Step 4: validation summary_statistics missing"
      touch "${validation_ck}"
    fi
  fi

  {
    echo "# Inspector baseline"
    cat "${eval_dir}/summary_statistics"
    if [[ -s "${valid_dir}/summary_statistics" ]]; then
      echo
      echo "# Inspector post-correction validation"
      cat "${valid_dir}/summary_statistics"
    fi
  } > "${dir4}/inspector.before_after.txt"

  touch "${dir4}/checkpoints/step4.done"
  info "Step 4: complete → ${step4_output_fasta}"
else
  step4_output_fasta="${step3_polished_fasta}"
  info "Step 4: disabled — carrying Step 3 assembly forward unchanged"
fi


# ============================================================================================================================
# STEP 5: Scaffolding — optional Tigmint correction, ntLink grid, CLR/HiFi refinement, native metrics/yak/Merqury selection
# Consumes: Inspector-corrected draft directly by default; optional HiFi-only Tigmint correction; optional CLR reads for ntLink/ARKS/gap filling
# Produces: ${dir5}/${base}.scaffolds.fa (selected refined, gap-filled scaffold set) + selection artifacts
# ============================================================================================================================
if is_true "${run_scaffolding}"; then
  echo "[INFO] Step 5: Scaffolding — starting"
  mkdir -p \
    "${dir5}/logs" \
    "${dir5}/checkpoints" \
    "${dir5}/kmers" \
    "${dir5}/genomescope2" \
    "${dir5}/tigmint" \
    "${dir5}/metrics" \
    "${dir5}/yak" \
    "${dir5}/merqury" \
    "${dir5}/selection" \
    "${dir5}/refinement" \
    "${dir5}/refinement/metrics" \
    "${dir5}/refinement/yak" \
    "${dir5}/refinement/merqury"

  scaff_ck_dir="${dir5}/checkpoints"
  scaff_logs_dir="${dir5}/logs"
  scaff_kmers_dir="${dir5}/kmers"
  scaff_gs_dir="${dir5}/genomescope2"
  scaff_tigmint_dir="${dir5}/tigmint"
  scaff_metrics_dir="${dir5}/metrics"
  scaff_yak_dir="${dir5}/yak"
  scaff_merqury_dir="${dir5}/merqury"
  scaff_selection_dir="${dir5}/selection"
  scaff_refine_dir="${dir5}/refinement"
  scaff_refined_metrics_dir="${scaff_refine_dir}/metrics"
  scaff_refined_yak_dir="${scaff_refine_dir}/yak"
  scaff_refined_merqury_dir="${scaff_refine_dir}/merqury"
  refined_manifest="${scaff_selection_dir}/refined_candidate_manifest.tsv"

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

  longstitch_cell_output() {
    local k="$1" w="$2"
    printf '%s/kw_k%s_w%s/%s.k%s.w%s.scaffolds.fa' "${dir5}" "$k" "$w" "${base}" "$k" "$w"
  }

  longstitch_cell_checkpoint() {
    local k="$1" w="$2"
    printf '%s/kw_k%s_w%s/.ntlink.ok' "${dir5}" "$k" "$w"
  }

  tigmint_output="${scaff_tigmint_dir}/${base}.tigmint.fa"
  tigmint_checkpoint="${scaff_tigmint_dir}/.tigmint.ok"
  tigmint_signature="${scaff_tigmint_dir}/tigmint.signature.tsv"
  grid_signature="${scaff_selection_dir}/ntlink_grid.signature.tsv"

  scaffolding_tigmint_ready() {
    [[ -s "${tigmint_output}" && -f "${tigmint_checkpoint}" ]]
  }

  scaffolding_longstitch_cell_ready() {
    local k="$1" w="$2" out ck
    out="$(longstitch_cell_output "$k" "$w")"
    ck="$(longstitch_cell_checkpoint "$k" "$w")"
    [[ -s "$out" ]] || return 1
    [[ -f "$ck" ]] || touch "$ck"
    return 0
  }

  scaffolding_longstitch_grid_ready() {
    local k w found=0
    read -r -a _grid_k_values <<< "${LONGSTITCH_K_VALUES}"
    read -r -a _grid_w_values <<< "${LONGSTITCH_W_VALUES}"
    (( ${#_grid_k_values[@]} > 0 && ${#_grid_w_values[@]} > 0 )) || return 1
    for k in "${_grid_k_values[@]}"; do
      for w in "${_grid_w_values[@]}"; do
        found=1
        scaffolding_longstitch_cell_ready "$k" "$w" || return 1
      done
    done
    (( found == 1 ))
  }

  # Native metrics must contain exactly one row per canonical LongStitch grid candidate.
  scaffolding_candidate_metrics_ready() {
    local report="${scaff_metrics_dir}/candidate_metrics.tsv"
    local manifest="${scaff_selection_dir}/candidate_manifest.tsv"
    local expected actual
    [[ -s "${report}" && -s "${manifest}" ]] || return 1
    expected="$(awk -F'\t' 'NR>1 && NF>=4 {n++} END{print n+0}' "${manifest}")"
    actual="$(awk -F'\t' 'NR>1 {n++} END{print n+0}' "${report}")"
    [[ "${expected}" =~ ^[0-9]+$ && "${actual}" =~ ^[0-9]+$ ]] || return 1
    (( expected > 0 && actual == expected ))
  }

  scaffolding_refinement_ready() {
    local manifest="${refined_manifest}" path found=0
    [[ -s "${manifest}" ]] || return 1
    while IFS=$'\t' read -r asm_id k w path; do
      [[ "${asm_id}" == "assembly_id" ]] && continue
      [[ -n "${asm_id}" && -s "${path}" ]] || return 1
      found=1
    done < "${manifest}"
    (( found == 1 ))
  }

  scaffolding_refined_metrics_ready() {
    local report="${scaff_refined_metrics_dir}/candidate_metrics.tsv" expected actual
    [[ -s "${report}" && -s "${refined_manifest}" ]] || return 1
    expected="$(awk -F'\t' 'NR>1 && NF>=4 {n++} END{print n+0}' "${refined_manifest}")"
    actual="$(awk -F'\t' 'NR>1 {n++} END{print n+0}' "${report}")"
    [[ "${expected}" =~ ^[0-9]+$ && "${actual}" =~ ^[0-9]+$ ]] || return 1
    (( expected > 0 && actual == expected ))
  }

  scaffolding_refined_merqury_ready() {
    local asm_id metrics_file found=0
    [[ -s "${refined_manifest}" ]] || return 1
    while IFS=$'\t' read -r asm_id _ _ _; do
      [[ "${asm_id}" == "assembly_id" ]] && continue
      [[ -n "${asm_id}" ]] || continue
      found=1
      metrics_file="${scaff_refined_merqury_dir}/${asm_id}/merqury_metrics.tsv"
      [[ -s "${metrics_file}" ]] || return 1
      awk -F'\t' 'NR > 1 && NF >= 4 { ok = 1 } END { exit ok ? 0 : 1 }' "${metrics_file}" || return 1
    done < "${refined_manifest}"
    (( found == 1 ))
  }

  scaffolding_output_ready() {
    local step="${1:-}"
    case "$step" in
      05A.inputs) [[ -s "${input_manifest:-}" ]] ;;
      05B.kmers) [[ -s "${gs_hist:-}" && -e "${read_meryl_db:-}" ]] ;;
      05C.genomescope2) scaffolding_genomescope2_ready ;;
      05D.longstitch_grid) scaffolding_longstitch_grid_ready ;;
      05E.metrics) scaffolding_candidate_metrics_ready ;;
      05F.yak) [[ -s "${yak_metrics:-}" ]] ;;
      05G.preselect) [[ -s "${scaff_selection_dir}/preselected_candidates.txt" && -s "${scaff_selection_dir}/preselection_scores.tsv" ]] ;;
      05H.refinement) scaffolding_refinement_ready ;;
      05I.refined_metrics) scaffolding_refined_metrics_ready ;;
      05J.refined_yak) [[ -s "${scaff_refined_yak_dir}/yak_metrics.tsv" ]] ;;
      05K.refined_merqury) scaffolding_refined_merqury_ready ;;
      05L.selection) [[ -s "${scaff_selection_dir}/best_assembly.path" && -s "${scaff_selection_dir}/longstitch_scores.tsv" ]] ;;
      05M.final_scaffolds) [[ -s "${dir5}/${base}.scaffolds.fa" && -s "${dir5}/best_assembly.path" ]] ;;
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
      "${dir4}/${base}.inspector.corrected.fa" \
      "${dir3}/${base}.round3.polished.fa" \
      "${dir2}/cleaned.fasta" \
      "${dir1}/${base}.bp.p_ctg.fa" \
      "${dir1}/${base}.asm.bp.p_ctg.fa"; do
      [[ -s "${cand}" ]] && { printf '%s\n' "${cand}"; return 0; }
    done
    return 1
  }

  resolve_tigmint_reads() {
    # Tigmint is a correction stage: always use HiFi reads, never CLR.
    if [[ -n "${pacbio_hifi_fasta:-}" && -s "${pacbio_hifi_fasta}" ]]; then
      printf '%s\n' "${pacbio_hifi_fasta}"
    elif [[ -n "${pacbio_hifi:-}" && -s "${pacbio_hifi}" ]]; then
      printf '%s\n' "${pacbio_hifi}"
    else
      return 1
    fi
  }

  resolve_scaffolding_reads() {
    # ntLink/ARKS/gap filling use CLR when available for long-range support.
    # Fall back to HiFi when CLR reads are absent.
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

  prepare_long_reads_fasta() {
    local in="$1" out="$2" label="${3:-long-read}"
    local tmp="${out}.tmp.$$"

    [[ -s "${in}" ]] || error "${label} input missing: ${in}"
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
          else { print "[ERROR] reads are not FASTA/FASTQ" > "/dev/stderr"; exit 2 }
        }
        mode == "fasta" { print; next }
        mode == "fastq" {
          fqline++
          if (fqline % 4 == 1) { sub(/^@/, ">"); print; next }
          if (fqline % 4 == 2) { print; next }
        }
        END { if (mode == "") exit 2 }
      ' > "${tmp}" || { rm -f "${tmp}"; error "Failed to convert ${label} reads to FASTA: ${in}"; }
    else
      awk '
        BEGIN { mode = ""; fqline = 0 }
        mode == "" {
          if ($0 ~ /^[[:space:]]*$/) next
          first = substr($0, 1, 1)
          if (first == ">") mode = "fasta"
          else if (first == "@") mode = "fastq"
          else { print "[ERROR] reads are not FASTA/FASTQ" > "/dev/stderr"; exit 2 }
        }
        mode == "fasta" { print; next }
        mode == "fastq" {
          fqline++
          if (fqline % 4 == 1) { sub(/^@/, ">"); print; next }
          if (fqline % 4 == 2) { print; next }
        }
        END { if (mode == "") exit 2 }
      ' "${in}" > "${tmp}" || { rm -f "${tmp}"; error "Failed to convert ${label} reads to FASTA: ${in}"; }
    fi

    [[ -s "${tmp}" ]] || { rm -f "${tmp}"; error "${label} FASTA conversion produced an empty file: ${out}"; }
    mv -f "${tmp}" "${out}"
    printf '%s\n' "${out}"
  }

  draft="$(resolve_scaffolding_draft)" || error "Step 5 input draft assembly not found"
  if is_false "${LONGSTITCH_USE_TIGMINT}"; then
    inspector_draft="${dir4}/${base}.inspector.corrected.fa"
    [[ -s "${inspector_draft}" ]] \
      || error "Step 5 Tigmint bypass requires the Inspector-corrected FASTA: ${inspector_draft}"
    draft="${inspector_draft}"
  fi
  scaffold_raw_reads="$(resolve_scaffolding_reads)" || error "No valid CLR/HiFi reads found for ntLink/ARKS scaffolding"

  tigmint_raw_reads=""
  tigmint_reads=""
  if is_true "${LONGSTITCH_USE_TIGMINT}"; then
    tigmint_raw_reads="$(resolve_tigmint_reads)" || error "Step 5 Tigmint is enabled but no valid PacBio HiFi reads were found"
    tigmint_reads="${tigmint_raw_reads}"
    case "${tigmint_raw_reads}" in
      *.gz|*.fq|*.fastq)
        tigmint_reads="$(prepare_long_reads_fasta "${tigmint_raw_reads}" "${dir5}/inputs/${base}.tigmint_hifi.fa" "Tigmint HiFi")"
        ;;
    esac
  fi

  scaffold_reads="${scaffold_raw_reads}"
  case "${scaffold_raw_reads}" in
    *.gz|*.fq|*.fastq)
      scaffold_reads="$(prepare_long_reads_fasta "${scaffold_raw_reads}" "${dir5}/inputs/${base}.scaffold_reads.fa" "scaffolding")"
      ;;
  esac

  # Backward-compatible alias used by a few downstream optional native-metric paths.
  longstitch_reads="${scaffold_reads}"
  longstitch_longmap="hifi"
  scaffold_read_type="hifi"
  if [[ -n "${pacbio_clr:-}" && "${scaffold_raw_reads}" == "${pacbio_clr}" ]]; then
    scaffold_read_type="clr"
  fi

  case "${SCAFFOLD_REFINEMENT}" in
    gapfill|rounds-gapfill|arks-gapfill|rounds-arks-gapfill) ;;
    *) error "Invalid SCAFFOLD_REFINEMENT=${SCAFFOLD_REFINEMENT}. Allowed: gapfill, rounds-gapfill, arks-gapfill, rounds-arks-gapfill" ;;
  esac
  [[ "${SCAFFOLD_REFINEMENT_TOP_N}" =~ ^[0-9]+$ ]] && (( SCAFFOLD_REFINEMENT_TOP_N > 0 )) \
    || error "SCAFFOLD_REFINEMENT_TOP_N must be a positive integer"
  [[ "${NTLINK_REFINEMENT_ROUNDS}" =~ ^[0-9]+$ ]] && (( NTLINK_REFINEMENT_ROUNDS >= 2 )) \
    || { [[ "${SCAFFOLD_REFINEMENT}" != *rounds* ]] || error "NTLINK_REFINEMENT_ROUNDS must be >=2 when rounds refinement is enabled"; }
  is_true "${NTLINK_OVERLAP}" \
    || error "NTLINK_OVERLAP must be TRUE because final ntLink gap filling requires overlap=True"

  if is_true "${LONGSTITCH_USE_TIGMINT}"; then
    echo "[INFO] Step 5 pre-scaffolding correction: Tigmint ENABLED with HiFi reads (${tigmint_raw_reads})"
  else
    echo "[INFO] Step 5 pre-scaffolding correction: Tigmint DISABLED; using Inspector-corrected FASTA directly (${draft})"
  fi
  echo "[INFO] Step 5 scaffold/refine reads: ${scaffold_read_type} (${scaffold_raw_reads})"
  echo "[INFO] Step 5 refinement mode: ${SCAFFOLD_REFINEMENT}; final gap filling is mandatory"

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
    printf "use_tigmint\t%s\n" "${LONGSTITCH_USE_TIGMINT}"
    printf "tigmint_raw_reads\t%s\n" "${tigmint_raw_reads}"
    printf "tigmint_reads\t%s\n" "${tigmint_reads}"
    printf "tigmint_longmap\t%s\n" "${longstitch_longmap}"
    printf "scaffold_raw_reads\t%s\n" "${scaffold_raw_reads}"
    printf "scaffold_reads\t%s\n" "${scaffold_reads}"
    printf "scaffold_read_type\t%s\n" "${scaffold_read_type}"
    printf "refinement_mode\t%s\n" "${SCAFFOLD_REFINEMENT}"
    printf "refinement_top_n\t%s\n" "${SCAFFOLD_REFINEMENT_TOP_N}"
    printf "refinement_rounds\t%s\n" "${NTLINK_REFINEMENT_ROUNDS}"
    printf "arks_m\t%s\n" "${ARKS_M}"
    printf "arks_j\t%s\n" "${ARKS_J}"
    printf "arks_k\t%s\n" "${ARKS_K}"
    printf "arks_c\t%s\n" "${ARKS_C}"
    printf "arks_l\t%s\n" "${ARKS_L}"
    printf "arks_a\t%s\n" "${ARKS_A}"
    printf "arks_cut\t%s\n" "${ARKS_CUT}"
    printf "tigmint_cut\t%s\n" "${LONGSTITCH_TIGMINT_CUT}"
    printf "tigmint_span\t%s\n" "${LONGSTITCH_TIGMINT_SPAN}"
    printf "tigmint_dist\t%s\n" "${LONGSTITCH_TIGMINT_DIST}"
    printf "ntlink_z\t%s\n" "${NTLINK_Z}"
    printf "ntlink_n\t%s\n" "${NTLINK_N}"
    printf "ntlink_a\t%s\n" "${NTLINK_A}"
    printf "ntlink_overlap\t%s\n" "${NTLINK_OVERLAP}"
    printf "ntlink_sensitive\t%s\n" "${NTLINK_SENSITIVE}"
    printf "kmer_reads\t%s\n" "${read_inputs[*]}"
    printf "longstitch_k_values\t%s\n" "${LONGSTITCH_K_VALUES}"
    printf "longstitch_w_values\t%s\n" "${LONGSTITCH_W_VALUES}"
  } > "${input_manifest}"
  scaffolding_mark_checkpoint 05A.inputs

  activate_env "${ENV_LONGSTITCH}"
  if is_true "${LONGSTITCH_USE_TIGMINT}"; then
    require longstitch
  fi
  require ntLink
  require parallel
  if [[ "${SCAFFOLD_REFINEMENT}" == *rounds* ]]; then
    require ntLink_rounds
  fi
  if [[ "${SCAFFOLD_REFINEMENT}" == *arks* ]]; then
    require arcs-make
  fi
  deactivate_env

  [[ -d "${GENOMESCOPE2_ENV_PREFIX}" ]] \
    || error "GenomeScope2 environment not found: ${GENOMESCOPE2_ENV_PREFIX}"

  genomescope_cmd="$(resolve_genomescope2_bin)" \
    || error "GenomeScope2 CLI not found or not executable under ${GENOMESCOPE2_ENV_PREFIX}/bin"

  [[ -x "${genomescope_cmd}" ]] \
    || error "GenomeScope2 CLI is not executable: ${genomescope_cmd}"

  "${genomescope_cmd}" --help \
    > "${scaff_logs_dir}/genomescope2.preflight.log" 2>&1 \
    || error "GenomeScope2 preflight failed: ${genomescope_cmd}"

  activate_env "${ENV_MERQURY}"
  require meryl
  if [[ -n "${MERQURY}" ]]; then
    require "${MERQURY}/merqury.sh"
  else
    require merqury.sh
  fi
  deactivate_env

  # Renv is invoked by absolute path rather than activated in the pipeline shell.
  # Conda-built R environments can install activate.d hooks that reference build-time
  # variables such as NATIVE_CONDA_TRIPLET. Under this script's `set -u`, sourcing
  # those hooks can abort the pipeline even though R itself is functional.
  assess_python="${MINIFORGE_PREFIX}/envs/${ENV_GENOFLOWSTATS}/bin/python"
  genoflowstats_bin="${MINIFORGE_PREFIX}/envs/${ENV_GENOFLOWSTATS}/bin/genoflowSTATS"
  [[ -x "${assess_python}" ]] || error "Python not found in ${ENV_GENOFLOWSTATS}: ${assess_python}"
  [[ -x "${genoflowstats_bin}" ]] || error "genoflowSTATS not found in ${ENV_GENOFLOWSTATS}: ${genoflowstats_bin}"
  [[ -s "${ASSESS_SCRIPT}" ]] || error "Embedded scaffold-selector helper not found: ${ASSESS_SCRIPT}"
  "${assess_python}" -c 'import sys,numpy,scipy,pysam,polars,pyarrow,statsmodels,sklearn,matplotlib; assert sys.version_info >= (3,14)' \
    >/dev/null || error "${ENV_GENOFLOWSTATS} integrity/import preflight failed"
  "${genoflowstats_bin}" --version >/dev/null 2>&1 || error "genoflowSTATS CLI preflight failed"
  genoflowstats_implementation_id "${assess_python}" >/dev/null \
    || error "Install the current genoflowSTATS 0.2.0 implementation with automatic FASTA length and palettes"

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
  gs_json=""
  gs_testing="${scaff_gs_dir}/${gs_output_prefix}testing.tsv"
  gs_metrics="${scaff_gs_dir}/genomescope2.metrics.tsv"
  if ! scaffolding_checkpoint_ready 05C.genomescope2; then
    genomescope_cmd="$(resolve_genomescope2_bin)" || error "GenomeScope2 CLI could not be resolved"
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
    [[ -s "${gs_summary}" ]] || error "GenomeScope2 summary missing: ${gs_summary}"
    [[ -s "${gs_model}" ]] || error "GenomeScope2 model missing: ${gs_model}"
    [[ -s "${gs_progress}" ]] || error "GenomeScope2 progress missing: ${gs_progress}"

    if is_true "${GS_JSON_REPORT}"; then
      gs_json="$(find "${scaff_gs_dir}" -maxdepth 1 -type f -name '*report.json' -print -quit 2>/dev/null || true)"
      [[ -s "${gs_json}" ]] || error "GenomeScope2 JSON report missing in ${scaff_gs_dir}"
    fi
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
  longstitch_genome_size="${longstitch_genome_size//,/}"
  [[ "${longstitch_genome_size}" =~ ^[0-9]+$ ]] \
    || error "Step 5 requires a positive integer haploid genome-size estimate for NG metrics/scaffolding assessment; got: ${longstitch_genome_size:-<empty>}"
  (( longstitch_genome_size > 0 )) \
    || error "Step 5 requires a positive haploid genome-size estimate for NG metrics/scaffolding assessment; got: ${longstitch_genome_size}"

  read -r -a k_values <<< "${LONGSTITCH_K_VALUES}"
  read -r -a w_values <<< "${LONGSTITCH_W_VALUES}"
  (( ${#k_values[@]} > 0 && ${#w_values[@]} > 0 )) \
    || error "LONGSTITCH_K_VALUES or LONGSTITCH_W_VALUES is empty"

  # -------------------- STEP 5D1: optional Tigmint correction --------------------
  # The default path bypasses Tigmint and uses the Inspector-corrected assembly
  # directly. When enabled, Tigmint uses HiFi reads only and runs once before the
  # ntLink grid because it is independent of minimizer k/w values.
  scaffold_input_assembly="${draft}"

  if is_true "${LONGSTITCH_USE_TIGMINT}"; then
    tigmint_signature_tmp="${scaff_selection_dir}/tigmint.signature.tmp.$$"
    {
      printf "key\tvalue\n"
      printf "enabled\tTRUE\n"
      printf "draft\t%s\n" "$(readlink -f "${draft}")"
      printf "draft_stat\t%s\n" "$(stat -c '%s:%Y' "${draft}")"
      printf "reads\t%s\n" "$(readlink -f "${tigmint_reads}")"
      printf "reads_stat\t%s\n" "$(stat -c '%s:%Y' "${tigmint_reads}")"
      printf "genome_size\t%s\n" "${longstitch_genome_size}"
      printf "longmap\t%s\n" "${longstitch_longmap}"
      printf "cut\t%s\n" "${LONGSTITCH_TIGMINT_CUT}"
      printf "span\t%s\n" "${LONGSTITCH_TIGMINT_SPAN}"
      printf "dist\t%s\n" "${LONGSTITCH_TIGMINT_DIST}"
    } > "${tigmint_signature_tmp}"

    tigmint_signature_changed=0
    if [[ ! -s "${tigmint_signature}" ]] || ! cmp -s "${tigmint_signature_tmp}" "${tigmint_signature}"; then
      tigmint_signature_changed=1
      echo "[INFO] Tigmint inputs/parameters changed; invalidating cached Tigmint correction"
      rm -f "${tigmint_checkpoint}" "${tigmint_output}"
    fi

    if (( tigmint_signature_changed == 1 )) || ! scaffolding_tigmint_ready; then
      activate_env "${ENV_LONGSTITCH}"
      require longstitch

      find "${scaff_tigmint_dir}" -mindepth 1 -maxdepth 1 -exec rm -rf {} +

      draft_link_name="${base}.draft.fa"
      reads_link_name="${base}.reads.fa"
      ln -s "$(readlink -f "${draft}")" "${scaff_tigmint_dir}/${draft_link_name}"
      ln -s "$(readlink -f "${tigmint_reads}")" "${scaff_tigmint_dir}/${reads_link_name}"

      echo "[INFO] Running optional Tigmint correction once before the ntLink k x w grid"
      pushd "${scaff_tigmint_dir}" >/dev/null || error "Cannot enter Tigmint directory: ${scaff_tigmint_dir}"
      set +e
      longstitch tigmint \
        draft="${base}.draft" \
        reads="${base}.reads" \
        cut="${LONGSTITCH_TIGMINT_CUT}" \
        t="${nthreads}" \
        G="${longstitch_genome_size}" \
        span="${LONGSTITCH_TIGMINT_SPAN}" \
        dist="${LONGSTITCH_TIGMINT_DIST}" \
        longmap="${longstitch_longmap}" \
        2>&1 | tee -a "${scaff_logs_dir}/tigmint.log"
      tigmint_rc=${PIPESTATUS[0]:-1}
      set -e
      (( tigmint_rc == 0 )) || { popd >/dev/null; error "Tigmint failed (rc=${tigmint_rc})"; }

      tigmint_native="${base}.draft.cut${LONGSTITCH_TIGMINT_CUT}.tigmint.fa"
      [[ -s "${tigmint_native}" ]] \
        || { popd >/dev/null; error "Tigmint did not produce expected corrected assembly: ${scaff_tigmint_dir}/${tigmint_native}"; }
      ln -s "${tigmint_native}" "$(basename "${tigmint_output}")"
      popd >/dev/null

      [[ -s "${tigmint_output}" ]] || error "Canonical Tigmint output is missing: ${tigmint_output}"
      mv -f "${tigmint_signature_tmp}" "${tigmint_signature}"
      touch "${tigmint_checkpoint}"
      deactivate_env
    else
      rm -f "${tigmint_signature_tmp}"
      echo "[INFO] Tigmint corrected assembly exists; reusing ${tigmint_output}"
    fi
    scaffold_input_assembly="${tigmint_output}"
  else
    # Tigmint artifacts are intentionally ignored in bypass mode. The input is
    # the complete Inspector-reconstructed corrected assembly from Step 4.
    echo "[INFO] Bypassing Tigmint; ntLink grid input is ${scaffold_input_assembly}"
  fi

  [[ -s "${scaffold_input_assembly}" ]] || error "Step 5 scaffold input assembly is missing: ${scaffold_input_assembly}"

  # Signature for the actual ntLink grid input. This deliberately includes the
  # Tigmint switch and resolved source assembly so resume cannot mix candidates
  # generated from Inspector-direct and Tigmint-corrected inputs.
  grid_signature_tmp="${scaff_selection_dir}/ntlink_grid.signature.tmp.$$"
  {
    printf "key\tvalue\n"
    printf "use_tigmint\t%s\n" "${LONGSTITCH_USE_TIGMINT}"
    printf "assembly\t%s\n" "$(readlink -f "${scaffold_input_assembly}")"
    printf "assembly_stat\t%s\n" "$(stat -c '%s:%Y' "${scaffold_input_assembly}")"
    printf "scaffold_reads\t%s\n" "$(readlink -f "${scaffold_reads}")"
    printf "scaffold_reads_stat\t%s\n" "$(stat -c '%s:%Y' "${scaffold_reads}")"
    printf "k_values\t%s\n" "${LONGSTITCH_K_VALUES}"
    printf "w_values\t%s\n" "${LONGSTITCH_W_VALUES}"
    printf "ntlink_z\t%s\n" "${NTLINK_Z}"
    printf "ntlink_n\t%s\n" "${NTLINK_N}"
    printf "ntlink_a\t%s\n" "${NTLINK_A}"
    printf "ntlink_overlap\t%s\n" "${NTLINK_OVERLAP}"
    printf "ntlink_sensitive\t%s\n" "${NTLINK_SENSITIVE}"
  } > "${grid_signature_tmp}"

  grid_signature_changed=0
  if [[ ! -s "${grid_signature}" ]] || ! cmp -s "${grid_signature_tmp}" "${grid_signature}"; then
    grid_signature_changed=1
    echo "[INFO] ntLink grid input/parameters changed; invalidating cached grid candidates and downstream selection"
    rm -f "$(scaffolding_checkpoint_file 05D.longstitch_grid)"
    for k in "${k_values[@]}"; do
      for w in "${w_values[@]}"; do
        rm -f "$(longstitch_cell_output "${k}" "${w}")" "$(longstitch_cell_checkpoint "${k}" "${w}")"
      done
    done
  fi

  # -------------------- STEP 5D2: ntLink broad k x w grid --------------------
  # The primary search intentionally uses one scaffolding round and no gap filling.
  # Iterative rounds/gap filling are refinement operations and should not multiply
  # the cost of the full grid search.
  run_scaffolding_kw() {
    local k="${1:-}" w="${2:-}"
    local run_dir out_scaff cell_ck ntlink_native ntlink_status

    [[ -n "${k}" && -n "${w}" ]] \
      || { echo "[ERROR] run_scaffolding_kw: missing k or w" >&2; return 2; }

    run_dir="${dir5}/kw_k${k}_w${w}"
    mkdir -p "${run_dir}"
    out_scaff="$(longstitch_cell_output "${k}" "${w}")"
    cell_ck="$(longstitch_cell_checkpoint "${k}" "${w}")"

    if [[ -s "${out_scaff}" && -f "${cell_ck}" ]]; then
      echo "[INFO] [k=${k}, w=${w}] canonical ntLink output exists; skipping"
      return 0
    fi
    if [[ -s "${out_scaff}" ]]; then
      echo "[INFO] [k=${k}, w=${w}] canonical ntLink FASTA exists; marking checkpoint"
      touch "${cell_ck}"
      return 0
    fi

    rm -f "${out_scaff}" "${cell_ck}"
    find "${run_dir}" -maxdepth 1 \( -type f -o -type l \) \
      -name 'target.fa.k*.ntLink*' -delete
    rm -f "${run_dir}/target.fa" "${run_dir}/reads.fa"
    ln -s "$(readlink -f "${scaffold_input_assembly}")" "${run_dir}/target.fa"
    ln -s "$(readlink -f "${scaffold_reads}")" "${run_dir}/reads.fa"

    echo "[INFO] [k=${k}, w=${w}] Running ntLink scaffold in ${run_dir}"
    pushd "${run_dir}" >/dev/null || { echo "[ERROR] cd failed for ${run_dir}" >&2; return 1; }

    set +e
    ntLink scaffold \
      target=target.fa \
      reads=reads.fa \
      t="${nthreads}" \
      k="${k}" \
      w="${w}" \
      n="${NTLINK_N}" \
      a="${NTLINK_A}" \
      z="${NTLINK_Z}" \
      overlap="${NTLINK_OVERLAP}" \
      sensitive="${NTLINK_SENSITIVE}" \
      > "ntlink.log" 2>&1
    ntlink_status=$?
    set -e
    if (( ntlink_status != 0 )); then
      echo "[ERROR] ntLink failed k=${k} w=${w} (rc=${ntlink_status}); see ${run_dir}/ntlink.log" >&2
      popd >/dev/null
      return "${ntlink_status}"
    fi

    ntlink_native="target.fa.k${k}.w${w}.z${NTLINK_Z}.ntLink.scaffolds.fa"
    if [[ ! -s "${ntlink_native}" ]]; then
      echo "[ERROR] ntLink did not produce expected scaffold FASTA: ${run_dir}/${ntlink_native}" >&2
      popd >/dev/null
      return 1
    fi

    ln -s "${ntlink_native}" "$(basename "${out_scaff}")"
    [[ -s "$(basename "${out_scaff}")" ]] \
      || { echo "[ERROR] canonical ntLink output link is invalid: ${out_scaff}" >&2; popd >/dev/null; return 1; }
    touch "$(basename "${cell_ck}")"
    popd >/dev/null
  }
  export -f run_scaffolding_kw
  export -f longstitch_cell_output longstitch_cell_checkpoint
  export dir5 base scaffold_reads scaffold_input_assembly nthreads
  export NTLINK_Z NTLINK_N NTLINK_A NTLINK_OVERLAP NTLINK_SENSITIVE

  grid_recomputed=0
  if ! scaffolding_checkpoint_ready 05D.longstitch_grid; then
    activate_env "${ENV_LONGSTITCH}"
    require ntLink
    grid_log="${scaff_logs_dir}/ntlink.grid.log"
    parallel --jobs "${njobs}" --load 80% --memfree 128G --tmpdir "${tmpdir}" \
      --env run_scaffolding_kw --env dir5 --env base --env scaffold_reads --env scaffold_input_assembly \
      --env nthreads --env NTLINK_Z --env NTLINK_N --env NTLINK_A \
      --env NTLINK_OVERLAP --env NTLINK_SENSITIVE \
      --env longstitch_cell_output --env longstitch_cell_checkpoint \
      run_scaffolding_kw "{1}" "{2}" ::: "${k_values[@]}" ::: "${w_values[@]}" \
      2>&1 | tee -a "${grid_log}"
    grid_rc=${PIPESTATUS[0]:-0}
    (( grid_rc == 0 )) || error "One or more ntLink grid jobs failed (see ${grid_log})"
    deactivate_env
    scaffolding_mark_checkpoint 05D.longstitch_grid
    mv -f "${grid_signature_tmp}" "${grid_signature}"
    grid_recomputed=1
  else
    rm -f "${grid_signature_tmp}"
  fi

  candidate_manifest="${scaff_selection_dir}/candidate_manifest.tsv"
  candidate_manifest_tmp="${candidate_manifest}.tmp.$$"
  printf "assembly_id\tk\tw\tpath\n" > "${candidate_manifest_tmp}"

  # Build the candidate set deterministically from the configured LongStitch
  # grid.  Never glob for *.scaffolds.fa here: each LongStitch run directory
  # contains many intermediate ntLink/gap-fill scaffold FASTAs that are not
  # independent grid candidates and must not enter native metrics or selection.
  for k_val in "${k_values[@]}"; do
    for w_val in "${w_values[@]}"; do
      f="$(longstitch_cell_output "${k_val}" "${w_val}")"
      [[ -s "${f}" ]] || error "Canonical LongStitch final scaffold missing for k=${k_val}, w=${w_val}: ${f}"

      asm_id="${base}.k${k_val}.w${w_val}.scaffolds"
      printf "%s\t%s\t%s\t%s\n" \
        "${asm_id}" "${k_val}" "${w_val}" "${f}" >> "${candidate_manifest_tmp}"
    done
  done

  (( "$(wc -l < "${candidate_manifest_tmp}")" > 1 )) \
    || { rm -f "${candidate_manifest_tmp}"; error "No canonical scaffold assemblies found after LongStitch grid"; }

  # Preserve manifest mtime when the grid candidate set is unchanged so that
  # downstream checkpoint logic remains stable across resumed runs.
  candidate_manifest_changed=0
  if [[ -s "${candidate_manifest}" ]] && cmp -s "${candidate_manifest_tmp}" "${candidate_manifest}"; then
    rm -f "${candidate_manifest_tmp}"
  else
    mv -f "${candidate_manifest_tmp}" "${candidate_manifest}"
    candidate_manifest_changed=1
  fi

  # Recomputed ntLink FASTAs invalidate all downstream metrics even when the
  # canonical candidate names are unchanged. A changed manifest also invalidates
  # them. This prevents native metrics/yak/Merqury checkpoints from referring to stale
  # assembly contents after a Tigmint or ntLink parameter/input change.
  if (( candidate_manifest_changed == 1 || grid_recomputed == 1 )); then
    echo "[INFO] Step 5 scaffold candidates changed or were recomputed; invalidating Step 5E-5M evaluation/refinement/selection checkpoints"
    rm -f \
      "$(scaffolding_checkpoint_file 05E.metrics)" \
      "$(scaffolding_checkpoint_file 05F.yak)" \
      "$(scaffolding_checkpoint_file 05G.preselect)" \
      "$(scaffolding_checkpoint_file 05H.refinement)" \
      "$(scaffolding_checkpoint_file 05I.refined_metrics)" \
      "$(scaffolding_checkpoint_file 05J.refined_yak)" \
      "$(scaffolding_checkpoint_file 05K.refined_merqury)" \
      "$(scaffolding_checkpoint_file 05L.selection)" \
      "$(scaffolding_checkpoint_file 05M.final_scaffolds)" \
      "${scaff_yak_dir}/yak_metrics.tsv" \
      "${scaff_selection_dir}/top_candidates.txt" \
      "${scaff_selection_dir}/preselected_candidates.txt" \
      "${scaff_selection_dir}/preselection_scores.tsv" \
      "${scaff_selection_dir}/refined_candidate_manifest.tsv" \
      "${scaff_selection_dir}/best_assembly.path" \
      "${scaff_selection_dir}/best_assembly.txt" \
      "${scaff_selection_dir}/longstitch_scores.tsv" \
      "${dir5}/best_assembly.path"
  fi

  expected_candidates=$((${#k_values[@]} * ${#w_values[@]}))
  observed_candidates="$(awk -F'\t' 'NR > 1 && NF >= 4 {n++} END {print n+0}' "${candidate_manifest}")"
  (( observed_candidates == expected_candidates )) \
    || error "LongStitch candidate manifest mismatch: expected ${expected_candidates}, found ${observed_candidates}"
  echo "[INFO] Step 5 native metric candidate set: ${observed_candidates} canonical LongStitch final scaffold assemblies"


  genoflowstats_candidate_metrics() {
    local source_manifest="${1:?missing source manifest}"
    local out_dir="${2:?missing output directory}"
    local table_name="${3:?missing table name}"
    local log_file="${4:?missing log file}"
    local gf_manifest="${out_dir}/genoflowSTATS.manifest.tsv"

    mkdir -p "${out_dir}"
    printf 'name\tfasta\n' > "${gf_manifest}.tmp"
    awk -F'\t' 'NR>1 && NF>=4 {print $1 "\t" $4}' "${source_manifest}" >> "${gf_manifest}.tmp"
    mv -f "${gf_manifest}.tmp" "${gf_manifest}"

    local -a args=(compare --manifest "${gf_manifest}" --outdir "${out_dir}" \
      --threads "${GENOFLOWSTATS_COMPARE_THREADS}" --gap-min "${GENOFLOWSTATS_GAP_MIN}" \
      --dpi "${GENOFLOWSTATS_PLOT_DPI}" --palette "${GENOFLOWSTATS_PALETTE}")
    [[ -n "${GENOME_SIZE_EST:-}" ]] && args+=(--genome-size "${GENOME_SIZE_EST}")
    "${genoflowstats_bin}" "${args[@]}" > "${log_file}" 2>&1 || return 1

    "${assess_python}" - "${source_manifest}" "${out_dir}/comparison.long.tsv" "${out_dir}/${table_name}" <<'PY_GF_CANDIDATE_TABLE'
import csv, math, sys
manifest, long_tsv, out = sys.argv[1:4]
meta={}
with open(manifest, newline='', encoding='utf-8') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        meta[r['assembly_id']] = r
vals={k:{} for k in meta}
with open(long_tsv, newline='', encoding='utf-8') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        a=r['assembly']
        if a not in vals or not r['metric'].startswith('assembly.'):
            continue
        key=r['metric'].split('.',1)[1]
        vals[a][key]=r['value']
cols=['assembly_id','k','w','path']
metric_cols=sorted({k for d in vals.values() for k in d})
with open(out,'w',newline='',encoding='utf-8') as fh:
    w=csv.writer(fh,delimiter='\t',lineterminator='\n')
    w.writerow(cols+metric_cols)
    for a,m in meta.items():
        w.writerow([a,m.get('k',''),m.get('w',''),m.get('path','')]+[vals[a].get(k,'NA') for k in metric_cols])
PY_GF_CANDIDATE_TABLE
  }

  # -------------------- STEP 5E: genoflowSTATS candidate metrics --------------------
  if ! scaffolding_checkpoint_ready 05E.metrics; then
    rm -rf "${scaff_metrics_dir:?}"/*
    genoflowstats_candidate_metrics "${candidate_manifest}" "${scaff_metrics_dir}" "candidate_metrics.tsv" \
      "${scaff_logs_dir}/genoflowSTATS.longstitch.log" \
      || error "Step 5 genoflowSTATS candidate assessment failed"
    scaffolding_candidate_metrics_ready || error "Step 5 candidate metrics are incomplete"
    scaffolding_mark_checkpoint 05E.metrics
  fi

  # -------------------- STEP 5F: Optional yak fast screening --------------------
  # Build the read k-mer database once, then evaluate independent scaffold
  # candidates concurrently.  nthreads is the CPU budget for each yak qv
  # process; njobs is the number of candidate-level yak qv processes allowed
  # to run simultaneously (maximum CPU demand ~= nthreads * njobs).
  yak_metrics="${scaff_yak_dir}/yak_metrics.tsv"
  if ! scaffolding_checkpoint_ready 05F.yak; then
    activate_env "${ENV_HIFIASM}"
    yak_db="${scaff_yak_dir}/${base}.reads.k${YAK_K}.yak"
    [[ -s "${yak_db}" ]] || yak count -k "${YAK_K}" -b "${YAK_BLOOM_BITS}" -t "${nthreads}" \
      -o "${yak_db}" "${read_inputs[@]}" \
      > "${scaff_logs_dir}/yak.count.log" 2>&1 \
      || error "yak count failed"

    yak_parts_dir="${scaff_yak_dir}/parallel_parts"
    rm -rf "${yak_parts_dir}"
    mkdir -p "${yak_parts_dir}"
    : > "${scaff_logs_dir}/yak.qv.log"

    run_yak_qv_candidate() {
      local asm_id="${1:?assembly ID required}"
      local asm_path="${2:?assembly path required}"
      local qv_out="${scaff_yak_dir}/${asm_id}.yak_qv.txt"
      local err_out="${yak_parts_dir}/${asm_id}.stderr.log"
      local row_out="${yak_parts_dir}/${asm_id}.tsv"
      local yak_line yak_qv

      if yak qv -t "${nthreads}" "${yak_db}" "${asm_path}" > "${qv_out}" 2> "${err_out}"; then
        yak_line="$(tail -n1 "${qv_out}" 2>/dev/null || true)"
        yak_qv="$(awk '{for (i=1;i<=NF;i++) if ($i ~ /^[0-9]+([.][0-9]+)?$/) {v=$i}} END{print v}' <<< "${yak_line}")"
      else
        yak_qv=""
      fi

      printf "%s\t%s\t\t%s\n" \
        "${asm_id}" "${yak_qv:-NA}" "${asm_path}" > "${row_out}"
    }

    export -f run_yak_qv_candidate
    export yak_db scaff_yak_dir yak_parts_dir nthreads

    yak_jobs_file="${yak_parts_dir}/jobs.tsv"
    awk -F '\t' 'NR > 1 {print $1 "\t" $4}' "${candidate_manifest}" > "${yak_jobs_file}"
    [[ -s "${yak_jobs_file}" ]] || error "No canonical LongStitch candidates available for yak evaluation"

    if ! parallel --jobs "${njobs}" --load 80% --memfree 128G --tmpdir "${tmpdir}" \
      --env run_yak_qv_candidate --env yak_db --env scaff_yak_dir --env yak_parts_dir --env nthreads \
      --colsep '\t' \
      run_yak_qv_candidate '{1}' '{2}' :::: "${yak_jobs_file}"; then
      error "GNU Parallel failed while dispatching yak QV jobs"
    fi

    printf "assembly_id\tyak_qv\tyak_error_rate\tpath\n" > "${yak_metrics}"
    while IFS=$'\t' read -r asm_id _ _ asm_path; do
      [[ "${asm_id}" == "assembly_id" ]] && continue
      row_file="${yak_parts_dir}/${asm_id}.tsv"
      [[ -s "${row_file}" ]] \
        || printf "%s\tNA\t\t%s\n" "${asm_id}" "${asm_path}" > "${row_file}"
      cat "${row_file}" >> "${yak_metrics}"
      if [[ -s "${yak_parts_dir}/${asm_id}.stderr.log" ]]; then
        {
          printf '[%s]\n' "${asm_id}"
          cat "${yak_parts_dir}/${asm_id}.stderr.log"
        } >> "${scaff_logs_dir}/yak.qv.log"
      fi
    done < "${candidate_manifest}"

    deactivate_env
    scaffolding_mark_checkpoint 05F.yak
  fi

  # -------------------- STEP 5G: Select top candidates before Merqury --------------------
  selector_py="${ASSESS_SCRIPT}"

  if ! scaffolding_checkpoint_ready 05G.preselect; then
    "${assess_python}" "${selector_py}" candidates \
      --metrics "${scaff_metrics_dir}/candidate_metrics.tsv" \
      --manifest "${candidate_manifest}" --yak "${yak_metrics}" \
      --merqury-dir "${scaff_merqury_dir}" --outdir "${scaff_selection_dir}" \
      --top-n "${SCAFFOLD_REFINEMENT_TOP_N}" --genome-size "${GENOME_SIZE_EST}" \
      --dpi "${GENOFLOWSTATS_PLOT_DPI}" \
      > "${scaff_logs_dir}/preselect_longstitch_candidates.log" 2>&1 \
      || error "Step 5 Python preselection failed"
    [[ -s "${scaff_selection_dir}/top_candidates.txt" ]] || error "Step 5 preselection produced no top candidates"
    head -n "${SCAFFOLD_REFINEMENT_TOP_N}" "${scaff_selection_dir}/top_candidates.txt" \
      > "${scaff_selection_dir}/preselected_candidates.txt"
    [[ -s "${scaff_selection_dir}/preselected_candidates.txt" ]] \
      || error "Step 5 preselection produced no candidates for refinement"
    scaffolding_mark_checkpoint 05G.preselect
  fi

  # -------------------- STEP 5H: Refine preselected candidates --------------------
  # Read-role policy:
  #   - Optional Tigmint correction, when enabled, always uses HiFi reads.
  #   - When Tigmint is disabled, the Inspector-corrected FASTA is scaffolded directly.
  #   - ntLink/ARKS/refinement use CLR reads when available, otherwise HiFi.
  #   - final ntLink gap filling is mandatory in this non-Hi-C pipeline.
  refinement_inputs="${scaff_selection_dir}/preselected_candidates.txt"
  refinement_signature="${scaff_refine_dir}/refinement.signature.tsv"
  refinement_signature_tmp="${scaff_refine_dir}/refinement.signature.tmp.$$"

  {
    printf "key\tvalue\n"
    printf "mode\t%s\n" "${SCAFFOLD_REFINEMENT}"
    printf "top_n\t%s\n" "${SCAFFOLD_REFINEMENT_TOP_N}"
    printf "scaffold_reads\t%s\n" "$(readlink -f "${scaffold_reads}")"
    printf "scaffold_reads_stat\t%s\n" "$(stat -c '%s:%Y' "${scaffold_reads}")"
    printf "rounds\t%s\n" "${NTLINK_REFINEMENT_ROUNDS}"
    printf "ntlink_z\t%s\n" "${NTLINK_Z}"
    printf "ntlink_n\t%s\n" "${NTLINK_N}"
    printf "ntlink_a\t%s\n" "${NTLINK_A}"
    printf "ntlink_overlap\t%s\n" "${NTLINK_OVERLAP}"
    printf "ntlink_sensitive\t%s\n" "${NTLINK_SENSITIVE}"
    printf "arks_m\t%s\n" "${ARKS_M}"
    printf "arks_j\t%s\n" "${ARKS_J}"
    printf "arks_k\t%s\n" "${ARKS_K}"
    printf "arks_c\t%s\n" "${ARKS_C}"
    printf "arks_l\t%s\n" "${ARKS_L}"
    printf "arks_a\t%s\n" "${ARKS_A}"
    printf "arks_cut\t%s\n" "${ARKS_CUT}"
    while IFS= read -r asm || [[ -n "${asm}" ]]; do
      [[ -n "${asm}" ]] || continue
      printf "candidate\t%s\t%s\n" "$(readlink -f "${asm}")" "$(stat -c '%s:%Y' "${asm}")"
    done < "${refinement_inputs}"
  } > "${refinement_signature_tmp}"

  refinement_signature_changed=0
  if [[ ! -s "${refinement_signature}" ]] || ! cmp -s "${refinement_signature_tmp}" "${refinement_signature}"; then
    refinement_signature_changed=1
    echo "[INFO] Step 5 refinement inputs/parameters changed; invalidating Step 5H-5M"
    rm -f \
      "$(scaffolding_checkpoint_file 05H.refinement)" \
      "$(scaffolding_checkpoint_file 05I.refined_metrics)" \
      "$(scaffolding_checkpoint_file 05J.refined_yak)" \
      "$(scaffolding_checkpoint_file 05K.refined_merqury)" \
      "$(scaffolding_checkpoint_file 05L.selection)" \
      "$(scaffolding_checkpoint_file 05M.final_scaffolds)" \
      "${refined_manifest}" \
      "${scaff_selection_dir}/best_assembly.path" \
      "${scaff_selection_dir}/best_assembly.txt" \
      "${scaff_selection_dir}/longstitch_scores.tsv" \
      "${dir5}/best_assembly.path"
    find "${scaff_refine_dir}" -mindepth 1 -maxdepth 1 -type d \
      ! -name metrics ! -name yak ! -name merqury -exec rm -rf {} +
    rm -rf "${scaff_refined_metrics_dir:?}"/* "${scaff_refined_yak_dir:?}"/* "${scaff_refined_merqury_dir:?}"/*
  fi

  run_refine_candidate() {
    local asm="${1:?candidate assembly required}"
    local asm_file asm_id k w out_dir current rounds_native arks_native gap_native refined
    local status

    [[ -s "${asm}" ]] || { echo "[ERROR] refinement candidate missing: ${asm}" >&2; return 1; }
    asm_file="$(basename "${asm}")"
    asm_id="${asm_file%.fa}"
    if [[ ! "${asm_id}" =~ \.k([0-9]+)\.w([0-9]+)\.scaffolds$ ]]; then
      echo "[ERROR] cannot parse ntLink k/w from candidate name: ${asm_file}" >&2
      return 1
    fi
    k="${BASH_REMATCH[1]}"
    w="${BASH_REMATCH[2]}"
    out_dir="${scaff_refine_dir}/${asm_id}"
    refined="${out_dir}/${asm_id}.refined.fa"
    mkdir -p "${out_dir}"

    if [[ -s "${refined}" && -f "${out_dir}/.refinement.ok" ]]; then
      echo "[INFO] [${asm_id}] refined candidate exists; skipping"
      return 0
    fi

    rm -rf "${out_dir:?}"/*
    ln -s "$(readlink -f "${asm}")" "${out_dir}/input.fa"
    ln -s "$(readlink -f "${scaffold_reads}")" "${out_dir}/reads.fa"
    current="input.fa"

    pushd "${out_dir}" >/dev/null || return 1

    if [[ "${SCAFFOLD_REFINEMENT}" == *rounds* ]]; then
      echo "[INFO] [${asm_id}] ntLink iterative refinement: ${NTLINK_REFINEMENT_ROUNDS} rounds"
      set +e
      ntLink_rounds run_rounds \
        target="${current}" \
        reads=reads.fa \
        rounds="${NTLINK_REFINEMENT_ROUNDS}" \
        t="${nthreads}" \
        k="${k}" \
        w="${w}" \
        n="${NTLINK_N}" \
        a="${NTLINK_A}" \
        z="${NTLINK_Z}" \
        overlap="${NTLINK_OVERLAP}" \
        sensitive="${NTLINK_SENSITIVE}" \
        > ntlink_rounds.log 2>&1
      status=$?
      set -e
      if (( status != 0 )); then
        echo "[ERROR] [${asm_id}] ntLink_rounds failed (rc=${status})" >&2
        popd >/dev/null
        return "${status}"
      fi
      rounds_native="${current}.k${k}.w${w}.z${NTLINK_Z}.ntLink.${NTLINK_REFINEMENT_ROUNDS}rounds.fa"
      [[ -s "${rounds_native}" ]] || { echo "[ERROR] [${asm_id}] missing rounds output: ${rounds_native}" >&2; popd >/dev/null; return 1; }
      current="${rounds_native}"
    fi

    if [[ "${SCAFFOLD_REFINEMENT}" == *arks* ]]; then
      rm -f arks_input.fa arks_result.scaffolds.fa
      ln -s "$(readlink -f "${current}")" arks_input.fa
      echo "[INFO] [${asm_id}] ARKS-long refinement using ${scaffold_read_type} reads"
      set +e
      arcs-make arks-long \
        draft=arks_input \
        reads=reads \
        m="${ARKS_M}" \
        cut="${ARKS_CUT}" \
        j="${ARKS_J}" \
        k="${ARKS_K}" \
        c="${ARKS_C}" \
        l="${ARKS_L}" \
        a="${ARKS_A}" \
        r="${ARKS_R}" \
        e="${ARKS_E}" \
        z="${NTLINK_Z}" \
        D="${ARKS_D}" \
        t="${nthreads}" \
        prefix=arks_result \
        > arks.log 2>&1
      status=$?
      set -e
      if (( status != 0 )); then
        echo "[ERROR] [${asm_id}] ARKS-long failed (rc=${status})" >&2
        popd >/dev/null
        return "${status}"
      fi
      arks_native="arks_result.scaffolds.fa"
      [[ -s "${arks_native}" ]] || { echo "[ERROR] [${asm_id}] missing ARKS output: ${arks_native}" >&2; popd >/dev/null; return 1; }
      current="${arks_native}"
    fi

    # Mandatory final gap-closing pass for the no-Hi-C workflow. ntLink gap_fill
    # fills scaffold gaps with raw long-read sequence, so later polishing is required.
    rm -f gap_input.fa
    ln -s "$(readlink -f "${current}")" gap_input.fa
    echo "[INFO] [${asm_id}] mandatory ntLink gap filling using ${scaffold_read_type} reads"
    set +e
    ntLink scaffold gap_fill \
      target=gap_input.fa \
      reads=reads.fa \
      t="${nthreads}" \
      k="${k}" \
      w="${w}" \
      n="${NTLINK_N}" \
      a="${NTLINK_A}" \
      z="${NTLINK_Z}" \
      overlap=True \
      sensitive="${NTLINK_SENSITIVE}" \
      > ntlink_gapfill.log 2>&1
    status=$?
    set -e
    if (( status != 0 )); then
      echo "[ERROR] [${asm_id}] ntLink gap filling failed (rc=${status})" >&2
      popd >/dev/null
      return "${status}"
    fi
    gap_native="gap_input.fa.k${k}.w${w}.z${NTLINK_Z}.ntLink.scaffolds.gap_fill.fa"
    [[ -s "${gap_native}" ]] || { echo "[ERROR] [${asm_id}] missing gap-filled output: ${gap_native}" >&2; popd >/dev/null; return 1; }

    ln -s "${gap_native}" "$(basename "${refined}")"
    [[ -s "$(basename "${refined}")" ]] || { popd >/dev/null; return 1; }
    touch .refinement.ok
    popd >/dev/null
  }
  export -f run_refine_candidate
  export scaff_refine_dir scaffold_reads scaffold_read_type nthreads SCAFFOLD_REFINEMENT
  export NTLINK_REFINEMENT_ROUNDS NTLINK_Z NTLINK_N NTLINK_A NTLINK_OVERLAP NTLINK_SENSITIVE
  export ARKS_M ARKS_J ARKS_K ARKS_C ARKS_L ARKS_A ARKS_CUT ARKS_R ARKS_E ARKS_D

  if ! scaffolding_checkpoint_ready 05H.refinement; then
    activate_env "${ENV_LONGSTITCH}"
    require ntLink
    [[ "${SCAFFOLD_REFINEMENT}" == *rounds* ]] && require ntLink_rounds
    [[ "${SCAFFOLD_REFINEMENT}" == *arks* ]] && require arcs-make

    parallel --jobs "${njobs}" --load 80% --memfree 128G --tmpdir "${tmpdir}" \
      --env run_refine_candidate --env scaff_refine_dir --env scaffold_reads --env scaffold_read_type \
      --env nthreads --env SCAFFOLD_REFINEMENT --env NTLINK_REFINEMENT_ROUNDS \
      --env NTLINK_Z --env NTLINK_N --env NTLINK_A --env NTLINK_OVERLAP --env NTLINK_SENSITIVE \
      --env ARKS_M --env ARKS_J --env ARKS_K --env ARKS_C --env ARKS_L --env ARKS_A \
      --env ARKS_CUT --env ARKS_R --env ARKS_E --env ARKS_D \
      run_refine_candidate {} :::: "${refinement_inputs}" \
      2>&1 | tee -a "${scaff_logs_dir}/refinement.parallel.log"
    refinement_rc=${PIPESTATUS[0]:-0}
    (( refinement_rc == 0 )) || error "One or more Step 5 refinement jobs failed"
    deactivate_env

    printf "assembly_id\tk\tw\tpath\n" > "${refined_manifest}"
    while IFS= read -r asm || [[ -n "${asm}" ]]; do
      [[ -n "${asm}" ]] || continue
      asm_file="$(basename "${asm}")"
      asm_id="${asm_file%.fa}"
      [[ "${asm_id}" =~ \.k([0-9]+)\.w([0-9]+)\.scaffolds$ ]] \
        || error "Cannot parse k/w for refined candidate: ${asm_file}"
      k="${BASH_REMATCH[1]}"; w="${BASH_REMATCH[2]}"
      refined_path="${scaff_refine_dir}/${asm_id}/${asm_id}.refined.fa"
      [[ -s "${refined_path}" ]] || error "Refined candidate missing: ${refined_path}"
      printf "%s.refined\t%s\t%s\t%s\n" "${asm_id}" "${k}" "${w}" "${refined_path}" >> "${refined_manifest}"
    done < "${refinement_inputs}"

    mv -f "${refinement_signature_tmp}" "${refinement_signature}"
    scaffolding_mark_checkpoint 05H.refinement
  else
    rm -f "${refinement_signature_tmp}"
  fi

  # -------------------- STEP 5I: genoflowSTATS refined-candidate metrics --------------------
  if ! scaffolding_checkpoint_ready 05I.refined_metrics; then
    rm -rf "${scaff_refined_metrics_dir:?}"/*
    genoflowstats_candidate_metrics "${refined_manifest}" "${scaff_refined_metrics_dir}" "candidate_metrics.tsv" \
      "${scaff_logs_dir}/genoflowSTATS.refined.log" \
      || error "Step 5 genoflowSTATS refined-candidate assessment failed"
    scaffolding_refined_metrics_ready || error "Refined candidate metrics are incomplete"
    scaffolding_mark_checkpoint 05I.refined_metrics
  fi

  # -------------------- STEP 5J: yak reevaluation of refined candidates --------------------
  refined_yak_metrics="${scaff_refined_yak_dir}/yak_metrics.tsv"
  if ! scaffolding_checkpoint_ready 05J.refined_yak; then
    activate_env "${ENV_HIFIASM}"
    yak_db="${scaff_yak_dir}/${base}.reads.k${YAK_K}.yak"
    [[ -s "${yak_db}" ]] || yak count -k "${YAK_K}" -b "${YAK_BLOOM_BITS}" -t "${nthreads}" \
      -o "${yak_db}" "${read_inputs[@]}" \
      > "${scaff_logs_dir}/yak.count.refined.log" 2>&1 \
      || error "yak count failed before refined-candidate evaluation"

    refined_yak_parts="${scaff_refined_yak_dir}/parallel_parts"
    rm -rf "${refined_yak_parts}"
    mkdir -p "${refined_yak_parts}"

    run_refined_yak_qv() {
      local asm_id="${1:?assembly ID required}" asm_path="${2:?assembly path required}"
      local qv_out="${scaff_refined_yak_dir}/${asm_id}.yak_qv.txt"
      local row_out="${refined_yak_parts}/${asm_id}.tsv" line qv
      if yak qv -t "${nthreads}" "${yak_db}" "${asm_path}" > "${qv_out}" 2> "${refined_yak_parts}/${asm_id}.stderr.log"; then
        line="$(tail -n1 "${qv_out}" 2>/dev/null || true)"
        qv="$(awk '{for (i=1;i<=NF;i++) if ($i ~ /^[0-9]+([.][0-9]+)?$/) {v=$i}} END{print v}' <<< "${line}")"
      else
        qv=""
      fi
      printf "%s\t%s\t\t%s\n" "${asm_id}" "${qv:-NA}" "${asm_path}" > "${row_out}"
    }
    export -f run_refined_yak_qv
    export scaff_refined_yak_dir refined_yak_parts yak_db nthreads

    refined_yak_jobs="${refined_yak_parts}/jobs.tsv"
    awk -F'\t' 'NR > 1 {print $1 "\t" $4}' "${refined_manifest}" > "${refined_yak_jobs}"
    parallel --jobs "${njobs}" --load 80% --memfree 128G --tmpdir "${tmpdir}" \
      --env run_refined_yak_qv --env scaff_refined_yak_dir --env refined_yak_parts --env yak_db --env nthreads \
      --colsep '\t' run_refined_yak_qv '{1}' '{2}' :::: "${refined_yak_jobs}" \
      || error "GNU Parallel failed during refined yak evaluation"

    printf "assembly_id\tyak_qv\tyak_error_rate\tpath\n" > "${refined_yak_metrics}"
    while IFS=$'\t' read -r asm_id _ _ asm_path; do
      [[ "${asm_id}" == "assembly_id" ]] && continue
      row="${refined_yak_parts}/${asm_id}.tsv"
      [[ -s "${row}" ]] || printf "%s\tNA\t\t%s\n" "${asm_id}" "${asm_path}" > "${row}"
      cat "${row}" >> "${refined_yak_metrics}"
    done < "${refined_manifest}"
    deactivate_env
    scaffolding_mark_checkpoint 05J.refined_yak
  fi

  # -------------------- STEP 5K: Merqury validation of refined candidates --------------------
  run_refined_merqury_candidate() {
    local asm_id="${1:?assembly ID required}" asm="${2:?assembly path required}"
    local out_dir="${scaff_refined_merqury_dir}/${asm_id}" out_prefix merqury_bin
    local qv_value="" comp_value="" f v
    [[ -s "${asm}" ]] || { echo "[ERROR] Merqury candidate missing: ${asm}" >&2; return 1; }
    mkdir -p "${out_dir}"
    if [[ -s "${out_dir}/merqury_metrics.tsv" ]]; then return 0; fi
    out_prefix="${out_dir}/${asm_id}"
    if [[ -n "${MERQURY:-}" && -x "${MERQURY}/merqury.sh" ]]; then merqury_bin="${MERQURY}/merqury.sh"; else merqury_bin=merqury.sh; fi
    "${merqury_bin}" "${read_meryl_db}" "${asm}" "${out_prefix}" > "${out_dir}/merqury.log" 2>&1 || return 1
    while IFS= read -r f; do
      v="$(awk 'NF >= 4 && $4 ~ /^[-+]?[0-9.]+([eE][-+]?[0-9]+)?$/ {print $4; exit}' "${f}" 2>/dev/null || true)"
      [[ -n "${v}" ]] && { qv_value="${v}"; break; }
    done < <(find "${out_dir}" -type f \( -name '*.qv' -o -name '*.qv.txt' -o -name '*qv*.tsv' \) -size +0c 2>/dev/null | sort)
    while IFS= read -r f; do
      v="$(awk 'NF >= 5 && $NF ~ /^[-+]?[0-9.]+([eE][-+]?[0-9]+)?$/ {x=$NF} END{print x}' "${f}" 2>/dev/null || true)"
      [[ -n "${v}" ]] && { comp_value="${v}"; break; }
    done < <(find "${out_dir}" -type f -name '*.completeness.stats' -size +0c 2>/dev/null | sort)
    {
      printf "assembly_id\tmerqury_qv\tmerqury_completeness\tpath\n"
      printf "%s\t%s\t%s\t%s\n" "${asm_id}" "${qv_value:-NA}" "${comp_value:-NA}" "${asm}"
    } > "${out_dir}/merqury_metrics.tsv"
  }
  export -f run_refined_merqury_candidate
  export scaff_refined_merqury_dir read_meryl_db MERQURY nthreads

  if ! scaffolding_checkpoint_ready 05K.refined_merqury; then
    activate_env "${ENV_MERQURY}"
    refined_merqury_jobs="${scaff_refine_dir}/merqury_jobs.tsv"
    awk -F'\t' 'NR > 1 {print $1 "\t" $4}' "${refined_manifest}" > "${refined_merqury_jobs}"
    parallel --jobs "${njobs}" --load 80% --memfree 64G --tmpdir "${tmpdir}" \
      --env run_refined_merqury_candidate --env scaff_refined_merqury_dir --env read_meryl_db --env MERQURY --env nthreads \
      --colsep '\t' run_refined_merqury_candidate '{1}' '{2}' :::: "${refined_merqury_jobs}" \
      2>&1 | tee -a "${scaff_logs_dir}/merqury.refined.parallel.log"
    merqury_rc=${PIPESTATUS[0]:-0}
    (( merqury_rc == 0 )) || error "Merqury failed for one or more refined Step 5 candidates"
    deactivate_env
    scaffolding_mark_checkpoint 05K.refined_merqury
  fi

  # -------------------- STEP 5L: Final selection among refined candidates --------------------
  if ! scaffolding_checkpoint_ready 05L.selection; then
    "${assess_python}" "${selector_py}" candidates \
      --metrics "${scaff_refined_metrics_dir}/candidate_metrics.tsv" \
      --manifest "${refined_manifest}" --yak "${refined_yak_metrics}" \
      --merqury-dir "${scaff_refined_merqury_dir}" --outdir "${scaff_selection_dir}" \
      --top-n "${SCAFFOLD_REFINEMENT_TOP_N}" --final --genome-size "${GENOME_SIZE_EST}" \
      --dpi "${GENOFLOWSTATS_PLOT_DPI}" \
      > "${scaff_logs_dir}/select_refined_best.log" 2>&1 \
      || error "Step 5 final refined-candidate selector failed"
    [[ -s "${scaff_selection_dir}/best_assembly.path" ]] || error "Final selector did not write best_assembly.path"
    scaffolding_mark_checkpoint 05L.selection
  fi

  # -------------------- STEP 5M: Final scaffold FASTA creation --------------------
  best_path="$(head -n1 "${scaff_selection_dir}/best_assembly.path")"
  [[ -s "${best_path}" ]] || error "Selected refined Step 5 assembly missing: ${best_path}"
  final_scaff="${dir5}/${base}.scaffolds.fa"
  if [[ ! -s "${final_scaff}" || "${best_path}" -nt "${final_scaff}" ]]; then
    awk '/^>/{print ">scaff" ++i; next}{print}' "${best_path}" > "${final_scaff}" \
      || error "Failed to canonicalize selected refined scaffold FASTA"
  fi
  printf '%s\n' "${best_path}" > "${dir5}/best_assembly.path"
  scaffolding_mark_checkpoint 05M.final_scaffolds
  echo "[INFO] Step 5 complete. Selected refined/gap-filled scaffolds: ${final_scaff} (source: ${best_path})"

fi

# Only remove ntLink/ARKS inputs and mapping scratch after their canonical
# FASTAs and cell checkpoints exist. These files are regenerated on a rebuild.
cleanup_completed_scaffolding_scratch() {
  local run_dir cell k w canonical scratch refined_dir refined_name
  [[ -f "${dir5}/checkpoints/05D.longstitch_grid.ok" ]] || return 0

  for run_dir in "${dir5}"/kw_k*_w*; do
    [[ -d "$run_dir" ]] || continue
    cell="${run_dir##*/}"
    [[ "$cell" =~ ^kw_k([0-9]+)_w([0-9]+)$ ]] || continue
    k="${BASH_REMATCH[1]}"
    w="${BASH_REMATCH[2]}"
    canonical="${run_dir}/${base}.k${k}.w${w}.scaffolds.fa"
    [[ -f "${run_dir}/.ntlink.ok" && -s "$canonical" ]] || continue
    rm -f -- "${run_dir}/reads.fa" "${run_dir}/target.fa" \
      "${run_dir}/target.fa.k${k}.w${w}.tsv"
    for scratch in "${run_dir}"/target.fa.k"${k}".w"${w}".z*.verbose_mapping.tsv \
                   "${run_dir}"/target.fa.k"${k}".w"${w}".z*.stitch.abyss-scaffold.fa \
                   "${run_dir}"/target.fa.k"${k}".w"${w}".z*.trimmed_scafs.fa; do
      [[ -e "$scratch" ]] || continue
      rm -f -- "$scratch"
    done
  done

  [[ -f "${dir5}/checkpoints/05H.refinement.ok" ]] || return 0
  for refined_dir in "${dir5}/refinement"/${base}.k*.w*.scaffolds; do
    [[ -d "$refined_dir" ]] || continue
    refined_name="${refined_dir##*/}"
    [[ -f "${refined_dir}/.refinement.ok" &&
       -s "${refined_dir}/${refined_name}.refined.fa" ]] || continue
    rm -f -- "${refined_dir}/reads.fa" "${refined_dir}/input.fa" \
      "${refined_dir}/arks_input.fa" "${refined_dir}/gap_input.fa"
    for scratch in "${refined_dir}"/gap_input.fa.k*.w*.tsv \
                   "${refined_dir}"/gap_input.fa.k*.w*.z*.verbose_mapping.tsv; do
      [[ -e "$scratch" ]] || continue
      rm -f -- "$scratch"
    done
  done
}

if is_true "${CLEANUP_COMPLETED_STEPS}"; then
  cleanup_completed_scaffolding_scratch
fi

# ==================================================================================================
# SHARED FINAL-FASTA LENGTH QC / FILTERING
# Runs SeqKit statistics and per-sequence length extraction, generates memory-efficient
# Python/Matplotlib length-distribution plots, filters scaffolds shorter than the configured
# minimum (default 1,000 bp), and returns the filtered FASTA path through the caller.
# ==================================================================================================
run_final_fasta_length_qc() {
  local stage="$1" source_fasta="$2" qc_dir="$3" filtered_fasta="$4"
  local min_len="${FINAL_FASTA_MIN_LENGTH}"
  local pre_stats="${qc_dir}/prefilter.seqkit_stats.tsv"
  local pre_lengths="${qc_dir}/prefilter.lengths.tsv"
  local pre_summary="${qc_dir}/prefilter.length_summary.tsv"
  local retention="${qc_dir}/filter_retention.tsv"
  local post_stats="${qc_dir}/postfilter.seqkit_stats.tsv"
  local post_lengths="${qc_dir}/postfilter.lengths.tsv"
  local py_log="${qc_dir}/python_length_qc.log"
  local py_resource="${qc_dir}/python_length_qc.resources.log"

  [[ -s "${source_fasta}" ]] || fatal "${stage}: final FASTA length-QC input missing: ${source_fasta}"
  command -v "${seqkit_bin}" >/dev/null 2>&1 || fatal "${stage}: SeqKit not found: ${seqkit_bin}"
  [[ "${min_len}" =~ ^[0-9]+$ ]] && (( min_len > 0 )) \
    || fatal "${stage}: FINAL_FASTA_MIN_LENGTH must be a positive integer"
  mkdir -p "${qc_dir}"

  info "${stage}: SeqKit pre-filter assembly statistics"
  "${seqkit_bin}" stats -a -T -j "${nthreads}" "${source_fasta}" > "${pre_stats}" \
    || fatal "${stage}: SeqKit pre-filter stats failed"
  "${seqkit_bin}" fx2tab -n -l -j "${nthreads}" "${source_fasta}" > "${pre_lengths}" \
    || fatal "${stage}: SeqKit per-scaffold length extraction failed"
  [[ -s "${pre_lengths}" ]] || fatal "${stage}: SeqKit produced no scaffold-length table"

  info "${stage}: plotting scaffold-length distribution and fixed ${min_len}-bp retention cutoff"
  local -a pycmd=("${python_bin}")
  if [[ -x /usr/bin/time ]]; then
    pycmd=(/usr/bin/time -v -o "${py_resource}" "${python_bin}")
  fi
  "${pycmd[@]}" - \
      "${pre_lengths}" "${qc_dir}" "${pre_summary}" "${retention}" \
      "${min_len}" "${FINAL_FASTA_QC_DPI}" "${FINAL_FASTA_QC_MAX_CURVE_POINTS}" \
      "${FINAL_FASTA_QC_SAVE_SVG}" \
      >"${py_log}" 2>&1 <<'PY_FINAL_FASTA_QC'
import sys
from pathlib import Path

import numpy as np

length_path = Path(sys.argv[1])
outdir = Path(sys.argv[2])
summary_path = Path(sys.argv[3])
retention_path = Path(sys.argv[4])
min_len = int(sys.argv[5])
dpi = max(600, int(sys.argv[6]))
max_curve_points = max(1000, int(sys.argv[7]))
save_svg = sys.argv[8].strip().lower() in {"1", "true", "yes", "y", "on"}

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception as exc:
    raise SystemExit(f"Matplotlib is required for final FASTA QC plots: {exc}")

vals = []
with length_path.open("rt", encoding="utf-8", errors="replace") as handle:
    for ln, raw in enumerate(handle, 1):
        if not raw.strip():
            continue
        fields = raw.rstrip("\n\r").split("\t")
        if len(fields) < 2:
            raise SystemExit(f"Malformed SeqKit fx2tab row {ln}: expected name and length")
        try:
            v = int(fields[-1])
        except ValueError as exc:
            raise SystemExit(f"Invalid scaffold length at row {ln}: {fields[-1]!r}") from exc
        if v < 0:
            raise SystemExit(f"Negative scaffold length at row {ln}")
        vals.append(v)

if not vals:
    raise SystemExit("No scaffold lengths were parsed")

lengths = np.asarray(vals, dtype=np.uint64)
del vals
lengths.sort()
n = int(lengths.size)
total_bp = int(lengths.sum(dtype=np.uint64))
if total_bp <= 0:
    raise SystemExit("Assembly contains zero sequence bases")

prefix = np.empty(n + 1, dtype=np.uint64)
prefix[0] = 0
np.cumsum(lengths, dtype=np.uint64, out=prefix[1:])

def nx_lx(frac):
    target = total_bp * frac
    desc_cum = total_bp - prefix[:-1]
    # prefix is ascending; search from largest end without materializing another full vector.
    lo, hi = 0, n - 1
    while lo < hi:
        mid = (lo + hi) // 2
        bp_from_mid = total_bp - int(prefix[mid])
        if bp_from_mid >= target:
            lo = mid + 1
        else:
            hi = mid
    idx = max(0, lo - 1 if total_bp - int(prefix[lo]) < target else lo)
    nx = int(lengths[idx])
    lx = n - idx
    return nx, lx

# Robust exact N50/L50 using a reverse cumulative loop over the sorted vector.
def exact_n50_l50():
    target = total_bp / 2.0
    acc = 0
    for j in range(n - 1, -1, -1):
        acc += int(lengths[j])
        if acc >= target:
            return int(lengths[j]), n - j
    return int(lengths[0]), n

n50, l50 = exact_n50_l50()
idx = int(np.searchsorted(lengths, min_len, side="left"))
kept_n = n - idx
kept_bp = total_bp - int(prefix[idx])
removed_n = idx
removed_bp = total_bp - kept_bp

median = float(np.median(lengths))
mean = float(total_bp / n)
min_bp = int(lengths[0])
max_bp = int(lengths[-1])
q1, q3 = np.percentile(lengths, [25, 75])

with summary_path.open("wt", encoding="utf-8") as out:
    out.write("metric\tvalue\n")
    for k, v in [
        ("scaffolds", n), ("total_bp", total_bp), ("min_bp", min_bp),
        ("mean_bp", f"{mean:.2f}"), ("median_bp", f"{median:.2f}"),
        ("q1_bp", f"{q1:.2f}"), ("q3_bp", f"{q3:.2f}"),
        ("max_bp", max_bp), ("N50_bp", n50), ("L50_scaffolds", l50),
        ("filter_min_length_bp", min_len),
    ]:
        out.write(f"{k}\t{v}\n")

with retention_path.open("wt", encoding="utf-8") as out:
    out.write("minimum_length_bp\tscaffolds_before\tscaffolds_retained\tscaffolds_removed\treads_retained_pct\tbases_before\tbases_retained\tbases_removed\tbases_retained_pct\n")
    out.write(
        f"{min_len}\t{n}\t{kept_n}\t{removed_n}\t{100.0*kept_n/n:.6f}\t"
        f"{total_bp}\t{kept_bp}\t{removed_bp}\t{100.0*kept_bp/total_bp:.6f}\n"
    )

# Histogram edges are computed once; plotting aggregated counts avoids passing the full
# sequence vector repeatedly to Matplotlib.
HEAT_DARK = "#7a1f59"
HEAT_RED = "#d94841"
HEAT_ORANGE = "#f98e52"
HEAT_GOLD = "#ffd166"

def finish_axes(ax):
    ax.set_facecolor("#fffaf0")
    ax.grid(True, color="#efe6d5", linewidth=0.8, alpha=0.9)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_visible(False)

def save_both(fig, stem):
    # SVG is the canonical publication output; PNG is retained for compatibility.
    fig.savefig(outdir / f"{stem}.svg", format="svg", dpi=dpi, bbox_inches="tight")
    fig.savefig(outdir / f"{stem}.png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)

if min_bp == max_bp:
    linear_edges = np.array([max(0, min_bp - 1), max_bp + 1], dtype=float)
else:
    linear_edges = np.linspace(min_bp, max_bp, min(201, max(31, int(np.sqrt(n)) + 1)))
counts, edges = np.histogram(lengths, bins=linear_edges)
fig, ax = plt.subplots(figsize=(10, 6))
ax.stairs(counts, edges, fill=True, color=HEAT_RED, alpha=0.78)
ax.axvline(min_len, linestyle="-.", linewidth=1.6, label=f"Minimum retained = {min_len:,} bp")
ax.axvline(median, linestyle="--", linewidth=1.4, label=f"Median = {median:,.0f} bp")
ax.axvline(n50, linestyle=":", linewidth=1.6, label=f"N50 = {n50:,} bp")
ax.set_xlabel("Scaffold length (bp)"); ax.set_ylabel("Number of scaffolds")
ax.set_title("Final assembly scaffold-length distribution")
ax.legend(frameon=False); finish_axes(ax); fig.tight_layout(); save_both(fig, "length_histogram")

positive = lengths[lengths > 0]
if positive.size:
    lo = max(1, int(positive[0])); hi = max(lo + 1, int(positive[-1]))
    log_edges = np.logspace(np.log10(lo), np.log10(hi), 100)
    if np.unique(log_edges).size >= 2:
        counts_log, edges_log = np.histogram(positive, bins=log_edges)
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.stairs(counts_log, edges_log, fill=True, color=HEAT_ORANGE, alpha=0.80)
        ax.set_xscale("log")
        ax.axvline(min_len, linestyle="-.", linewidth=1.6, label=f"Minimum retained = {min_len:,} bp")
        ax.axvline(median, linestyle="--", linewidth=1.4, label=f"Median = {median:,.0f} bp")
        ax.axvline(n50, linestyle=":", linewidth=1.6, label=f"N50 = {n50:,} bp")
        ax.set_xlabel("Scaffold length (bp, log scale)"); ax.set_ylabel("Number of scaffolds")
        ax.set_title("Final assembly scaffold-length distribution")
        ax.legend(frameon=False); finish_axes(ax); fig.tight_layout(); save_both(fig, "length_histogram_logx")

    # Yield-weighted histogram: sum scaffold bases per length bin.
    weights = positive.astype(np.float64, copy=False)
    yield_counts, yield_edges = np.histogram(positive, bins=log_edges, weights=weights)
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.stairs(yield_counts, yield_edges, fill=True, color=HEAT_GOLD, alpha=0.82)
    ax.set_xscale("log")
    ax.axvline(min_len, linestyle="-.", linewidth=1.6, label=f"Minimum retained = {min_len:,} bp")
    ax.set_xlabel("Scaffold length (bp, log scale)"); ax.set_ylabel("Assembly span (bp)")
    ax.set_title("Yield-weighted scaffold-length distribution")
    ax.legend(frameon=False); finish_axes(ax); fig.tight_layout(); save_both(fig, "yield_weighted_length_histogram")

# Downsample ECDF/survival plotting only; exact statistics still use all scaffolds.
step = max(1, int(np.ceil(n / max_curve_points)))
plot_idx = np.arange(0, n, step, dtype=np.int64)
plot_x = lengths[plot_idx]
ecdf = (plot_idx + 1).astype(np.float64) / n
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(plot_x, ecdf, linewidth=1.8, color=HEAT_DARK)
ax.axvline(min_len, linestyle="-.", linewidth=1.4)
ax.set_xscale("log" if min_bp > 0 else "linear")
ax.set_xlabel("Scaffold length (bp)"); ax.set_ylabel("Cumulative fraction of scaffolds")
ax.set_title("Scaffold-length ECDF"); finish_axes(ax); fig.tight_layout(); save_both(fig, "length_ecdf")

survival = 1.0 - ecdf + (1.0 / n)
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(plot_x, survival, linewidth=1.8, color=HEAT_RED)
ax.axvline(min_len, linestyle="-.", linewidth=1.4)
ax.set_xscale("log" if min_bp > 0 else "linear")
ax.set_xlabel("Scaffold length (bp)"); ax.set_ylabel("Fraction of scaffolds ≥ length")
ax.set_title("Scaffold-length survival curve"); finish_axes(ax); fig.tight_layout(); save_both(fig, "length_survival")

print(f"Scaffolds before filtering: {n:,}")
print(f"Assembly span before filtering: {total_bp:,} bp")
print(f"N50/L50 before filtering: {n50:,} bp / {l50:,}")
print(f"Cutoff: {min_len:,} bp")
print(f"Scaffolds retained: {kept_n:,} ({100.0*kept_n/n:.2f}%)")
print(f"Bases retained: {kept_bp:,} ({100.0*kept_bp/total_bp:.4f}%)")
PY_FINAL_FASTA_QC
  local py_rc=$?
  (( py_rc == 0 )) || { cat "${py_log}" >&2 || true; fatal "${stage}: Python final-FASTA length QC failed"; }

  info "${stage}: removing scaffolds shorter than ${min_len} bp with SeqKit"
  "${seqkit_bin}" seq -m "${min_len}" -w 80 "${source_fasta}" > "${filtered_fasta}.tmp" \
    || fatal "${stage}: SeqKit minimum-length filtering failed"
  [[ -s "${filtered_fasta}.tmp" ]] || fatal "${stage}: minimum-length filtering removed every scaffold"
  mv -f "${filtered_fasta}.tmp" "${filtered_fasta}"

  "${seqkit_bin}" stats -a -T -j "${nthreads}" "${filtered_fasta}" > "${post_stats}" \
    || fatal "${stage}: SeqKit post-filter stats failed"
  "${seqkit_bin}" fx2tab -n -l -j "${nthreads}" "${filtered_fasta}" > "${post_lengths}" \
    || fatal "${stage}: SeqKit post-filter length extraction failed"

  local residual_short
  residual_short=$(awk -F'\t' -v m="${min_len}" '$NF < m {n++} END{print n+0}' "${post_lengths}")
  (( residual_short == 0 )) || fatal "${stage}: ${residual_short} scaffolds shorter than ${min_len} bp remain after SeqKit filtering"

  info "${stage}: final FASTA length filter retained $(awk 'NR==2{print $3}' "${retention}" 2>/dev/null || echo '?') scaffolds"
}

# ==================================================================================================
# SHARED POST-SCAFFOLD REFINEMENT — DeepVariant HiFi polishing + final NCBI FCS
# This implementation is intentionally identical in genoflowHiFi and genoflowHiFi_HiC.
# Arguments to run_post_scaffold_refinement:
#   1) display step number, 2) output directory, 3) canonical scaffold FASTA input.
# Produces:
#   <outdir>/<base>.postpolish.final.fa
#   <outdir>/final.polished.decontam.fa (when RUN_FINAL_FCS=TRUE)
#   <outdir>/<species_name>.fasta      (canonical published FASTA)
# ==================================================================================================
run_post_scaffold_refinement() {
  local step_no="$1" outdir="$2" initial_assembly="$3"
  local stage="Step ${step_no}"
  local final_polished="${outdir}/${base}.postpolish.final.fa"
  local published="${outdir}/${species_name}.fasta"
  local fcs_final="${outdir}/final.polished.decontam.fa"

  info "${stage}: shared post-scaffold DeepVariant + FCS refinement — starting"
  mkdir -p "${outdir}/logs" "${outdir}/checkpoints" "${outdir}/metrics" \
           "${outdir}/adaptor" "${outdir}/gx" "${outdir}/provenance"

  [[ -s "${initial_assembly}" ]] || fatal "${stage}: canonical scaffold input missing: ${initial_assembly}"
  require "${minimap2_bin}"
  require "${samtools_bin}"
  require "${bcftools_bin}"
  command -v singularity >/dev/null 2>&1 || fatal "${stage}: singularity is required for DeepVariant/FCS"
  [[ -s "${sif_deepvariant}" ]] || fatal "${stage}: DeepVariant image missing: ${sif_deepvariant}"

  "${minimap2_bin}" --version > "${outdir}/provenance/minimap2.version.txt" 2>&1 || fatal "${stage}: minimap2 version check failed"
  "${samtools_bin}" --version > "${outdir}/provenance/samtools.version.txt" 2>&1 || fatal "${stage}: samtools version check failed"
  "${bcftools_bin}" --version > "${outdir}/provenance/bcftools.version.txt" 2>&1 || fatal "${stage}: bcftools version check failed"
  singularity exec --cleanenv "${sif_deepvariant}" /opt/deepvariant/bin/run_deepvariant --version \
    > "${outdir}/provenance/deepvariant.version.txt" 2>&1 \
    || warn "${stage}: DeepVariant version string could not be captured; execution will still be tested"

  local hifi_input_spec="${PACBIO_HIFI_READS:-${pacbio_hifi}}"
  local -a hifi_reads=()
  read -r -a hifi_reads <<< "${hifi_input_spec}"
  (( ${#hifi_reads[@]} > 0 )) || fatal "${stage}: no PacBio HiFi reads configured"
  local rf
  for rf in "${hifi_reads[@]}"; do [[ -s "${rf}" ]] || fatal "${stage}: HiFi input missing: ${rf}"; done

  # Assembly SHA256 is exact; huge read files use an inexpensive path/size/mtime signature unless explicitly requested.
  post_file_signature() {
    local f="$1"
    if is_true "${POST_POLISH_HASH_READS}"; then
      sha256sum "${f}" | awk '{print $1}'
    else
      local rp sz mt
      rp=$(readlink -f "${f}" 2>/dev/null || printf '%s' "${f}")
      sz=$(stat -c '%s' "${f}")
      mt=$(stat -c '%Y' "${f}")
      printf '%s|%s|%s' "${rp}" "${sz}" "${mt}" | sha256sum | awk '{print $1}'
    fi
  }


  post_fasta_stats() {
    local fa="$1"
    "${samtools_bin}" faidx "${fa}" >/dev/null 2>&1 || return 1
    awk 'BEGIN{n=0;bp=0} {n++;bp+=$2} END{printf "%d\t%d\n",n,bp}' "${fa}.fai"
  }

  local initial_sha
  initial_sha=$(sha256sum "${initial_assembly}" | awk '{print $1}')
  {
    printf 'key\tvalue\n'
    printf 'input_fasta\t%s\n' "${initial_assembly}"
    printf 'input_sha256\t%s\n' "${initial_sha}"
    printf 'post_rounds\t%s\n' "${DEEPVARIANT_POST_ROUNDS}"
    printf 'deepvariant_exclude_flags\t%s\n' "${DEEPVARIANT_EXCLUDE_FLAGS}"
    printf 'deepvariant_view_expr\t%s\n' "${DEEPVARIANT_VIEW_EXPR}"
    printf 'deepvariant_make_examples_extra_args\t%s\n' "${DEEPVARIANT_MAKE_EXAMPLES_EXTRA_ARGS}"
    printf 'minimap2_index_batch\t%s\n' "${MINIMAP2_INDEX_BATCH}"
    printf 'samtools_sort_mem\t%s\n' "${SAMTOOLS_SORT_MEM}"
    printf 'bam_compression\t%s\n' "${SAMTOOLS_BAM_COMPRESSION}"
    printf 'checked_at\t%s\n' "$(date -Is)"
  } > "${outdir}/provenance/post_refinement.input.tsv"

  local asm_in="${initial_assembly}" completed_round=0 round
  for ((round=1; round<=DEEPVARIANT_POST_ROUNDS; round++)); do
    local prefix="${outdir}/${base}.round${round}"
    local sorted_bam="${prefix}.sorted.bam"
    local dv_vcf="${prefix}.dv.vcf.gz"
    local correction_vcf="${prefix}.corrections.vcf.gz"
    local asm_out="${prefix}.polished.fa"
    local chain_file="${prefix}.consensus.chain"
    local metrics_file="${outdir}/metrics/round${round}.corrections.tsv"
    local round_ck="${outdir}/checkpoints/round${round}.done"
    local fp_file="${outdir}/checkpoints/round${round}.fingerprint.sha256"
    local fp_tmp="${fp_file}.tmp"
    local ref_sha read_sig config_fp

    [[ -s "${asm_in}.fai" ]] || "${samtools_bin}" faidx "${asm_in}" || fatal "${stage} round ${round}: faidx failed"
    ref_sha=$(sha256sum "${asm_in}" | awk '{print $1}')
    read_sig=""
    for rf in "${hifi_reads[@]}"; do read_sig+="$(post_file_signature "${rf}")|"; done
    config_fp=$(printf '%s\n' \
      "ref=${ref_sha}" "reads=${read_sig}" "model=PACBIO" "shards=${DEEPVARIANT_SHARDS}" \
      "postcpus=${DEEPVARIANT_FINAL_POSTPROCESS_CPUS}" \
      "exclude_flags=${DEEPVARIANT_EXCLUDE_FLAGS}" \
      "view_expr=${DEEPVARIANT_VIEW_EXPR}" \
      "make_examples_extra=${DEEPVARIANT_MAKE_EXAMPLES_EXTRA_ARGS}" \
      "mm2I=${MINIMAP2_INDEX_BATCH}" "mm2K=${MINIMAP2_BATCH_SIZE}" \
      "sortmem=${SAMTOOLS_SORT_MEM}" "bamlevel=${SAMTOOLS_BAM_COMPRESSION}" \
      "writechain=${DEEPVARIANT_WRITE_CHAIN}" | sha256sum | awk '{print $1}')
    printf '%s\n' "${config_fp}" > "${fp_tmp}"

    if [[ -f "${round_ck}" && -s "${asm_out}" && -s "${fp_file}" ]] && cmp -s "${fp_tmp}" "${fp_file}"; then
      rm -f "${fp_tmp}"
      info "${stage} round ${round}: fingerprinted checkpoint validated"
      completed_round="${round}"
      asm_in="${asm_out}"
      continue
    fi
    rm -f "${round_ck}" "${fp_file}"

    local -a mm2_args=("${minimap2_bin}" -a -x map-hifi --secondary=no -Y -K "${MINIMAP2_BATCH_SIZE}" -t "${nthreads}" \
      -R "@RG\\tID:${base}.HiFi\\tSM:${base}\\tPL:PACBIO\\tLB:SMRTbell")
    mm2_args+=("${asm_in}" "${hifi_reads[@]}")

    local bam_for_dv
    info "${stage} round ${round}: HiFi mapping → DeepVariant-compatible BAM → coordinate sort + BAI"
    "${mm2_args[@]}" \
      | "${samtools_bin}" view -u -F "${DEEPVARIANT_EXCLUDE_FLAGS}" -e "${DEEPVARIANT_VIEW_EXPR}" - \
      | "${samtools_bin}" sort --write-index -@ "${nthreads}" -m "${SAMTOOLS_SORT_MEM}" -l "${SAMTOOLS_BAM_COMPRESSION}" \
          -T "${outdir}/round${round}.coord.tmp" -O BAM \
          -o "${sorted_bam}##idx##${sorted_bam}.bai" - \
      || fatal "${stage} round ${round}: minimap2/SAMtools DeepVariant BAM preparation failed"
    bam_for_dv="${sorted_bam}"
    "${samtools_bin}" quickcheck -v "${bam_for_dv}" || fatal "${stage} round ${round}: BAM quickcheck failed"
    [[ -s "${bam_for_dv}.bai" ]] || fatal "${stage} round ${round}: DeepVariant BAI index missing"
    if [[ "$("${samtools_bin}" view -c -f "${DEEPVARIANT_EXCLUDE_FLAGS}" "${bam_for_dv}")" -ne 0 ]]; then
      fatal "${stage} round ${round}: DeepVariant BAM still contains secondary/supplementary alignments"
    fi
    if [[ "$("${samtools_bin}" view -c -e 'length(seq) == 0 || length(qual) != length(seq)' "${bam_for_dv}")" -ne 0 ]]; then
      fatal "${stage} round ${round}: DeepVariant BAM still contains records without base qualities"
    fi

    local dv_tmp="${outdir}/deepvariant_tmp.round${round}"
    local dv_logdir="${outdir}/logs/deepvariant.round${round}"
    rm -rf "${dv_tmp}" "${dv_logdir}"; mkdir -p "${dv_tmp}" "${dv_logdir}"
    local -a dv_binds=(-B "${working_dir}:${working_dir}" -B "${tmpdir}:${tmpdir}" -B "${outdir}:${outdir}")
    local d
    for rf in "${hifi_reads[@]}" "${asm_in}" "${bam_for_dv}"; do
      d=$(dirname "$(readlink -f "${rf}")")
      dv_binds+=(-B "${d}:${d}")
    done
    singularity exec --cleanenv "${dv_binds[@]}" "${sif_deepvariant}" \
      /opt/deepvariant/bin/run_deepvariant \
        --model_type=PACBIO --ref="${asm_in}" --reads="${bam_for_dv}" \
        --output_vcf="${dv_vcf}" --num_shards="${DEEPVARIANT_SHARDS}" \
        --postprocess_cpus="${DEEPVARIANT_FINAL_POSTPROCESS_CPUS}" \
        --intermediate_results_dir="${dv_tmp}" --logging_dir="${dv_logdir}" \
        --make_examples_extra_args="${DEEPVARIANT_MAKE_EXAMPLES_EXTRA_ARGS}" \
        --vcf_stats_report="${DEEPVARIANT_VCF_STATS}" \
      2>&1 | tee "${outdir}/logs/round${round}.deepvariant.${LOG_TS}.log" \
      || { rm -rf "${dv_tmp}"; fatal "${stage} round ${round}: DeepVariant failed"; }
    rm -rf "${dv_tmp}"
    [[ -s "${dv_vcf}" ]] || fatal "${stage} round ${round}: DeepVariant VCF missing"
    "${bcftools_bin}" view -h "${dv_vcf}" >/dev/null || fatal "${stage} round ${round}: DeepVariant VCF validation failed"

    local -a dv_samples=()
    mapfile -t dv_samples < <("${bcftools_bin}" query -l "${dv_vcf}")
    (( ${#dv_samples[@]} == 1 )) || fatal "${stage} round ${round}: expected exactly one DeepVariant sample, found ${#dv_samples[@]}"
    local dv_sample="${dv_samples[0]}"

    "${bcftools_bin}" view --threads "${nthreads}" -f PASS -i 'GT="AA"' -Oz -o "${correction_vcf}" --write-index "${dv_vcf}" \
      || fatal "${stage} round ${round}: PASS homozygous-alt filtering failed"
    local total_calls pass_calls corrections snps indels
    total_calls=$("${bcftools_bin}" view -H "${dv_vcf}" | wc -l | awk '{print $1}')
    pass_calls=$("${bcftools_bin}" view -H -f PASS "${dv_vcf}" | wc -l | awk '{print $1}')
    corrections=$("${bcftools_bin}" view -H "${correction_vcf}" | wc -l | awk '{print $1}')
    snps=$("${bcftools_bin}" view -H -v snps "${correction_vcf}" | wc -l | awk '{print $1}')
    indels=$("${bcftools_bin}" view -H -v indels "${correction_vcf}" | wc -l | awk '{print $1}')
    {
      printf 'metric\tvalue\n'
      printf 'total_calls\t%s\n' "${total_calls}"
      printf 'pass_calls\t%s\n' "${pass_calls}"
      printf 'corrections_total\t%s\n' "${corrections}"
      printf 'snps\t%s\n' "${snps}"
      printf 'indels\t%s\n' "${indels}"
      printf 'sample\t%s\n' "${dv_sample}"
    } > "${metrics_file}"

    if (( corrections == 0 )); then
      cp -f "${asm_in}" "${asm_out}"
      info "${stage} round ${round}: zero applicable corrections; polishing converged"
    else
      local -a consensus_args=("${bcftools_bin}" consensus --fasta-ref "${asm_in}" --samples "${dv_sample}" \
        --haplotype A --output "${asm_out}.tmp")
      is_true "${DEEPVARIANT_WRITE_CHAIN}" && consensus_args+=(--chain "${chain_file}")
      consensus_args+=("${correction_vcf}")
      "${consensus_args[@]}" \
        2>&1 | tee "${outdir}/logs/round${round}.consensus.${LOG_TS}.log" \
        || fatal "${stage} round ${round}: bcftools consensus failed"
      mv -f "${asm_out}.tmp" "${asm_out}"
    fi
    [[ -s "${asm_out}" ]] || fatal "${stage} round ${round}: polished FASTA missing"
    "${samtools_bin}" faidx "${asm_out}" || fatal "${stage} round ${round}: polished FASTA validation failed"

    # Consensus polishing may alter sequence lengths but must never lose/reorder FASTA records.
    local in_n out_n in_bp out_bp
    read -r in_n in_bp < <(post_fasta_stats "${asm_in}")
    read -r out_n out_bp < <(post_fasta_stats "${asm_out}")
    [[ "${in_n}" == "${out_n}" ]] || fatal "${stage} round ${round}: FASTA record count changed (${in_n} -> ${out_n})"
    printf 'input_records\t%s\noutput_records\t%s\ninput_bp\t%s\noutput_bp\t%s\n' \
      "${in_n}" "${out_n}" "${in_bp}" "${out_bp}" >> "${metrics_file}"

    mv -f "${fp_tmp}" "${fp_file}"
    touch "${round_ck}"
    completed_round="${round}"
    asm_in="${asm_out}"

    if is_true "${DEEPVARIANT_CLEAN_INTERMEDIATES}"; then
      rm -f "${sorted_bam}" "${sorted_bam}.bai" "${sorted_bam}.csi"
      if is_false "${POST_POLISH_KEEP_VCF}"; then
        rm -f "${dv_vcf}" "${dv_vcf}.tbi" "${dv_vcf}.csi" \
              "${correction_vcf}" "${correction_vcf}.tbi" "${correction_vcf}.csi"
      fi
    fi
    is_true "${DEEPVARIANT_STOP_ON_ZERO}" && (( corrections == 0 )) && break
  done

  (( completed_round > 0 )) || fatal "${stage}: no post-scaffold polishing round completed"
  cp -f "${asm_in}" "${final_polished}.tmp" && mv -f "${final_polished}.tmp" "${final_polished}"
  "${samtools_bin}" faidx "${final_polished}" || fatal "${stage}: final polished FASTA failed validation"
  printf 'completed_round\t%s\nfinal_polished_fasta\t%s\n' "${completed_round}" "${final_polished}" \
    > "${outdir}/post_polishing.summary.tsv"

  local final_source="${final_polished}"

  # Dedicated final-FCS checkpoint. This prevents rerunning FCS-Adaptor/FCS-GX
  # when only the downstream SeqKit final-FASTA QC/publication needs to be rerun.
  local final_fcs_ck="${outdir}/checkpoints/final_fcs.done"
  local final_fcs_fp_file="${outdir}/checkpoints/final_fcs.fingerprint.sha256"
  local final_polished_sha final_fcs_fp
  final_polished_sha="$(sha256sum "${final_polished}" | awk '{print $1}')"
  final_fcs_fp="$(printf '%s\n' \
    "final_polished_sha256=${final_polished_sha}" \
    "run_final_fcs=${RUN_FINAL_FCS}" \
    "fcs_adaptor_min_seq_len=${FCS_ADAPTOR_MIN_SEQ_LEN}" \
    "fcs_adaptor_max_seq_len=${FCS_ADAPTOR_MAX_SEQ_LEN}" \
    "fcs_gx_taxid=${FCS_GX_TAXID}" \
    "fcs_gx_internal_action=${FCS_FINAL_GX_INTERNAL_ACTION}" \
    "fcs_gx_align_exclude_taxa=${FCS_GX_ALIGN_EXCLUDE_TAXA}" \
    "fcs_gx_extra_contam_divs=${FCS_GX_EXTRA_CONTAM_DIVS}" \
    "fcs_gx_db=${FCS_GX_DB}" \
    "fcs_gx_db_policy=${FCS_GX_DB_POLICY}" \
    "fcs_gx_manifest=${FCS_GX_MANIFEST}" \
    "fcs_adaptor_image=${sif_fcs_adaptor}" \
    "fcs_gx_image=${sif_fcs_gx}" \
    | sha256sum | awk '{print $1}')"

  if is_true "${RUN_FINAL_FCS}"; then
    local final_fcs_checkpoint_valid="FALSE"
    if [[ -s "${final_fcs_fp_file}" && -e "${final_fcs_ck}" && -s "${fcs_final}" ]]; then
      local saved_final_fcs_fp
      saved_final_fcs_fp="$(awk 'NF{print $1; exit}' "${final_fcs_fp_file}")"
      if [[ "${saved_final_fcs_fp}" == "${final_fcs_fp}" ]]; then
        if "${samtools_bin}" faidx "${fcs_final}" >/dev/null 2>&1; then
          final_fcs_checkpoint_valid="TRUE"
        else
          warn "${stage}: final FCS checkpoint exists but cleaned FASTA validation failed; FCS will rerun"
        fi
      else
        info "${stage}: final FCS checkpoint fingerprint changed; FCS will rerun"
      fi
    fi

    if is_true "${final_fcs_checkpoint_valid}"; then
      info "${stage}: final FCS checkpoint valid — reusing ${fcs_final}"
      final_source="${fcs_final}"
    else
      # Migration path for successful v0.7.1 runs that predate this checkpoint.
      # Only adopt the old FCS output when all major completion evidence exists,
      # the cleaned FASTA is newer than its polished input, and faidx validation passes.
      local legacy_gx_report=""
      legacy_gx_report="$(find "${outdir}/gx" -maxdepth 1 -type f -name '*.fcs_gx_report.txt' -size +0c 2>/dev/null | sort | head -n 1 || true)"
      if [[ ! -e "${final_fcs_ck}" && ! -s "${final_fcs_fp_file}" \
            && -s "${fcs_final}" && "${fcs_final}" -nt "${final_polished}" \
            && -s "${outdir}/adaptor/clean.fasta" \
            && -s "${outdir}/gx_database_status.tsv" \
            && -n "${legacy_gx_report}" ]] \
         && "${samtools_bin}" faidx "${fcs_final}" >/dev/null 2>&1; then
        info "${stage}: adopting existing validated v0.7.1 final FCS output and creating checkpoint"
        printf '%s\n' "${final_fcs_fp}" > "${final_fcs_fp_file}.tmp"
        mv -f "${final_fcs_fp_file}.tmp" "${final_fcs_fp_file}"
        touch "${final_fcs_ck}"
        final_source="${fcs_final}"
      else
        info "${stage}: final FCS-Adaptor + FCS-GX screening/cleaning — starting"
    [[ -s "${sif_fcs_adaptor}" ]] || fatal "${stage}: FCS-Adaptor image missing: ${sif_fcs_adaptor}"
    [[ -s "${sif_fcs_gx}" ]] || fatal "${stage}: FCS-GX image missing: ${sif_fcs_gx}"
    [[ -s "${FCS_RUNNER}" ]] || fatal "${stage}: fcs.py runner missing: ${FCS_RUNNER}"
    [[ -x "${ADAPTOR_RUNNER}" || -s "${ADAPTOR_RUNNER}" ]] || fatal "${stage}: FCS-Adaptor runner missing: ${ADAPTOR_RUNNER}"

    local fcs_input="${outdir}/final.polished.fcs_input.fa"
    "${python_bin}" - "${final_polished}" "${fcs_input}.tmp" "${FCS_ADAPTOR_MIN_SEQ_LEN}" "${FCS_ADAPTOR_MAX_SEQ_LEN}" <<'PY_FCS_FINAL'
import re, sys
src, dst, min_len, max_len = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
valid_id = re.compile(r'^[A-Za-z0-9_.:-]+$')
valid_seq = re.compile(r'^[ACGTN]+$', re.I)
seen=set(); records=0
with open(src) as inp, open(dst,'w') as out:
    name=None; seq=[]
    def emit(n, parts):
        if n is None: return 0
        s=''.join(parts).upper()
        if len(s) < min_len:
            sys.stderr.write(f'FCS-skip short sequence {n}: {len(s)} bp < minimum {min_len} bp\n')
            return 0
        if len(s) > max_len:
            raise SystemExit(f'FCS-invalid length for {n}: {len(s)} bp > maximum {max_len} bp')
        if not valid_seq.fullmatch(s):
            raise SystemExit(f'FCS-invalid sequence alphabet in {n}')
        if s.count('N')/len(s) > 0.50:
            raise SystemExit(f'FCS-invalid N fraction in {n}')
        out.write(f'>{n}\n')
        for i in range(0,len(s),80):
            out.write(s[i:i+80]+'\n')
        return 1
    for raw in inp:
        line=raw.strip()
        if not line: continue
        if line.startswith('>'):
            records += emit(name, seq)
            hdr=line[1:].strip(); rid=hdr.split()[0] if hdr else ''
            if not rid or len(rid) >= 50 or not valid_id.fullmatch(rid): raise SystemExit(f'FCS-invalid FASTA identifier: {rid!r}')
            if rid in seen: raise SystemExit(f'Duplicate FASTA identifier: {rid}')
            seen.add(rid); name=rid; seq=[]
        else:
            if name is None: raise SystemExit('Sequence before first FASTA header')
            seq.append(line)
    records += emit(name, seq)
if records == 0: raise SystemExit('No FASTA records found')
PY_FCS_FINAL
    mv -f "${fcs_input}.tmp" "${fcs_input}"
    "${samtools_bin}" faidx "${fcs_input}" || fatal "${stage}: prepared FCS FASTA failed validation"

    local adaptor_report="${outdir}/adaptor/fcs_adaptor_report.txt"
    local adaptor_clean="${outdir}/adaptor/clean.fasta"
    local adaptor_contam="${outdir}/adaptor/contam.fasta"
    local adaptor_tax_flag="--euk"
    [[ "${FCS_ADAPTOR_TAXONOMY}" == "prok" ]] && adaptor_tax_flag="--prok"
    [[ "${FCS_ADAPTOR_TAXONOMY}" =~ ^(euk|prok)$ ]] || fatal "${stage}: FCS_ADAPTOR_TAXONOMY must be euk or prok"
    local -a adaptor_cmd=("${ADAPTOR_RUNNER}" --fasta-input "${fcs_input}" --output-dir "${outdir}/adaptor" \
      "${adaptor_tax_flag}" --container-engine singularity --image "${sif_fcs_adaptor}")
    is_true "${FCS_DEBUG}" && adaptor_cmd+=(--debug)
    "${adaptor_cmd[@]}" 2>&1 | tee "${outdir}/logs/fcs_adaptor.${LOG_TS}.log" \
      || fatal "${stage}: FCS-Adaptor failed"
    [[ -s "${adaptor_report}" ]] || fatal "${stage}: FCS-Adaptor report missing"

    local adaptor_actions="${adaptor_report}"
    if [[ "${FCS_FINAL_ADAPTOR_INTERNAL_ACTION}" == "mask" ]]; then
      adaptor_actions="${outdir}/adaptor/fcs_adaptor_report.mask-internal.txt"
      awk -F'\t' 'BEGIN{OFS="\t"} /^#/ {print; next} {
        if ($3=="ACTION_TRIM") { split($4,r,/\.\./); if ((r[1]+0)>1 && (r[2]+0)<($2+0)) $3="FIX" }
        print
      }' "${adaptor_report}" > "${adaptor_actions}"
    elif [[ "${FCS_FINAL_ADAPTOR_INTERNAL_ACTION}" != "split" ]]; then
      fatal "${stage}: FCS_FINAL_ADAPTOR_INTERNAL_ACTION must be mask or split"
    fi
    FCS_DEFAULT_IMAGE="${sif_fcs_gx}" "${python_bin}" "${FCS_RUNNER}" clean genome \
      --action-report "${adaptor_actions}" --output "${adaptor_clean}" --contam-fasta-out "${adaptor_contam}" \
      < "${fcs_input}" 2>&1 | tee "${outdir}/logs/fcs_adaptor.clean.${LOG_TS}.log" \
      || fatal "${stage}: applying FCS-Adaptor actions failed"
    [[ -s "${adaptor_clean}" ]] || fatal "${stage}: adaptor-cleaned FASTA missing"

    # GX DB lifecycle check: local integrity and remote freshness remain separate, exactly as in early FCS.
    local gx_policy="${FCS_GX_DB_POLICY,,}"
    [[ "${gx_policy}" == "validate" ]] && gx_policy="check"
    case "${gx_policy}" in check|update|require_latest|offline|download) ;; *) fatal "${stage}: invalid FCS_GX_DB_POLICY=${FCS_GX_DB_POLICY}" ;; esac
    local local_mft="${FCS_GX_DB}/all.manifest"
    local remote_mft="${outdir}/gx_database_remote.manifest"
    local db_lock="${FCS_GX_DB}.genoflow-update.lock"
    local local_state="INVALID" freshness="NOT_CHECKED" action="NONE" local_sha="NA" remote_sha="NA"
    local integrity_mode="FULL"
    is_true "${FCS_GX_DB_CHECK}" || integrity_mode="SKIPPED"

    local -a gx_required_files=(
      all.assemblies.tsv all.blast_div.tsv.gz all.gxi all.gxs all.manifest
      all.meta.jsonl all.README.txt all.seq_info.tsv.gz all.taxa.tsv
    )
    local gx_file
    local gx_complete="TRUE"
    for gx_file in "${gx_required_files[@]}"; do
      if [[ ! -s "${FCS_GX_DB}/${gx_file}" ]]; then
        gx_complete="FALSE"
        break
      fi
    done

    if [[ "${gx_complete}" == "TRUE" && -s "${local_mft}" ]]; then
      if is_true "${FCS_GX_DB_CHECK}"; then
        if "${python_bin}" "${FCS_RUNNER}" db check --mft "${local_mft}" --dir "${FCS_GX_DB}" \
          > "${outdir}/logs/fcs_gx.db_check.local.${LOG_TS}.log" 2>&1; then
          local_state="VALID"
        fi
      else
        info "${stage}: skipping full FCS-GX database integrity scan (FCS_GX_DB_CHECK=FALSE)"
        local_state="VALID"
      fi
      [[ "${local_state}" == "VALID" ]] && local_sha=$(sha256sum "${local_mft}" | awk '{print $1}')
    fi
    if [[ "${gx_policy}" != "offline" ]]; then
      if [[ "${FCS_GX_MANIFEST}" =~ ^https?:// ]]; then
        if curl --fail --location --silent --show-error --retry 3 --retry-delay 5 \
          -o "${remote_mft}.tmp" "${FCS_GX_MANIFEST}"; then
          mv -f "${remote_mft}.tmp" "${remote_mft}"
        else
          rm -f "${remote_mft}.tmp"
        fi
      elif [[ -s "${FCS_GX_MANIFEST}" ]]; then
        cp -f "${FCS_GX_MANIFEST}" "${remote_mft}"
      fi
      if [[ -s "${remote_mft}" ]]; then
        remote_sha=$(sha256sum "${remote_mft}" | awk '{print $1}')
        if [[ -s "${local_mft}" ]] && cmp -s "${local_mft}" "${remote_mft}"; then freshness="UP_TO_DATE"; else freshness="UPDATE_AVAILABLE"; fi
      else
        freshness="REMOTE_CHECK_FAILED"
        [[ "${gx_policy}" == "check" ]] || fatal "${stage}: current GX manifest could not be retrieved for policy=${gx_policy}"
      fi
    else
      freshness="OFFLINE_NOT_CHECKED"
    fi

    case "${gx_policy}" in
      offline) [[ "${local_state}" == "VALID" ]] || fatal "${stage}: offline FCS policy requires a valid local GX DB"; action="USE_VALIDATED_LOCAL" ;;
      check)
        [[ "${local_state}" == "VALID" ]] || fatal "${stage}: GX DB invalid/missing; use FCS_GX_DB_POLICY=update or download"
        [[ "${freshness}" == "UPDATE_AVAILABLE" ]] && warn "${stage}: newer/different GX manifest available; policy=check retains validated local DB"
        action="USE_VALIDATED_LOCAL"
        ;;
      require_latest)
        [[ "${local_state}" == "VALID" && "${freshness}" == "UP_TO_DATE" ]] || fatal "${stage}: require_latest requires a valid current GX DB"
        action="USE_UP_TO_DATE_LOCAL"
        ;;
      update|download)
        if [[ "${local_state}" != "VALID" || "${freshness}" != "UP_TO_DATE" ]]; then
          mkdir "${db_lock}" 2>/dev/null || fatal "${stage}: GX DB update lock exists: ${db_lock}"
          owned_lock_dirs+=("${db_lock}")
          if ! "${python_bin}" "${FCS_RUNNER}" db get --mft "${FCS_GX_MANIFEST}" --dir "${FCS_GX_DB}" \
            2>&1 | tee "${outdir}/logs/fcs_gx.db_get.${LOG_TS}.log"; then
            release_lock_dir "${db_lock}"
            fatal "${stage}: GX DB synchronization failed"
          fi
          release_lock_dir "${db_lock}"
          [[ -s "${local_mft}" ]] || fatal "${stage}: GX DB synchronization did not install all.manifest"
          for gx_file in "${gx_required_files[@]}"; do
            [[ -s "${FCS_GX_DB}/${gx_file}" ]] || fatal "${stage}: GX DB synchronization left required file missing/empty: ${gx_file}"
          done
          if is_true "${FCS_GX_DB_CHECK}"; then
            "${python_bin}" "${FCS_RUNNER}" db check --mft "${local_mft}" --dir "${FCS_GX_DB}" \
              > "${outdir}/logs/fcs_gx.db_check.after_update.${LOG_TS}.log" 2>&1 \
              || fatal "${stage}: GX DB failed validation after synchronization"
          else
            info "${stage}: post-update full GX integrity scan skipped (FCS_GX_DB_CHECK=FALSE)"
          fi
          local_state="VALID"; local_sha=$(sha256sum "${local_mft}" | awk '{print $1}')
          if [[ -s "${remote_mft}" ]] && ! cmp -s "${local_mft}" "${remote_mft}"; then
            fatal "${stage}: synchronized GX local manifest does not match requested current manifest"
          fi
          freshness="UP_TO_DATE"; action="SYNCHRONIZED"
        else action="USE_UP_TO_DATE_LOCAL"; fi
        ;;
    esac
    [[ "${local_state}" == "VALID" ]] || fatal "${stage}: final GX DB integrity state is not valid"
    cp -f "${local_mft}" "${outdir}/gx_database_manifest.used"
    {
      printf 'key\tvalue\n'; printf 'database_path\t%s\n' "${FCS_GX_DB}"; printf 'policy\t%s\n' "${gx_policy}"
      printf 'db_check_enabled\t%s\n' "${FCS_GX_DB_CHECK}"; printf 'integrity_check\t%s\n' "${integrity_mode}"
      printf 'local_integrity\t%s\n' "${local_state}"; printf 'freshness\t%s\n' "${freshness}"; printf 'action\t%s\n' "${action}"
      printf 'local_manifest_sha256\t%s\n' "${local_sha}"; printf 'remote_manifest_sha256\t%s\n' "${remote_sha}"
      printf 'manifest_source\t%s\n' "${FCS_GX_MANIFEST}"; printf 'checked_at\t%s\n' "$(date -Is)"
    } > "${outdir}/gx_database_status.tsv"

    local gx_env="${outdir}/gx/fcs_gx.env"
    {
      echo "GX_NUM_CORES=${FCS_GX_CORES}"
      [[ -n "${FCS_GX_ALIGN_EXCLUDE_TAXA}" ]] && echo "GX_ALIGN_EXCLUDE_TAXA=${FCS_GX_ALIGN_EXCLUDE_TAXA}"
      [[ -n "${FCS_GX_EXTRA_CONTAM_DIVS}" ]] && echo "GX_EXTRA_CONTAM_DIVS=${FCS_GX_EXTRA_CONTAM_DIVS}"
    } > "${gx_env}"
    FCS_DEFAULT_IMAGE="${sif_fcs_gx}" "${python_bin}" "${FCS_RUNNER}" --env-file "${gx_env}" screen genome \
      --fasta "${adaptor_clean}" --out-dir "${outdir}/gx" --gx-db "${FCS_GX_DB}" --tax-id "${FCS_GX_TAXID}" \
      2>&1 | tee "${outdir}/logs/fcs_gx.screen.${LOG_TS}.log" || fatal "${stage}: FCS-GX screen failed"
    local -a gx_reports=()
    mapfile -t gx_reports < <(find "${outdir}/gx" -maxdepth 1 -type f -name '*.fcs_gx_report.txt' -size +0c | sort)
    (( ${#gx_reports[@]} == 1 )) || fatal "${stage}: expected exactly one GX action report; found ${#gx_reports[@]}"
    local gx_report="${gx_reports[0]}" gx_actions="${gx_reports[0]}"
    if [[ "${FCS_FINAL_GX_INTERNAL_ACTION}" == "split" ]]; then
      gx_actions="${outdir}/gx/fcs_gx_report.split-internal.txt"
      sed 's/\tFIX\t/\tSPLIT\t/g' "${gx_report}" > "${gx_actions}"
    elif [[ "${FCS_FINAL_GX_INTERNAL_ACTION}" != "fix" ]]; then
      fatal "${stage}: FCS_FINAL_GX_INTERNAL_ACTION must be fix or split"
    fi
    FCS_DEFAULT_IMAGE="${sif_fcs_gx}" "${python_bin}" "${FCS_RUNNER}" clean genome \
      --action-report "${gx_actions}" --output "${fcs_final}.tmp" --contam-fasta-out "${outdir}/gx/contam.fasta" \
      < "${adaptor_clean}" 2>&1 | tee "${outdir}/logs/fcs_gx.clean.${LOG_TS}.log" \
      || fatal "${stage}: applying FCS-GX actions failed"
    [[ -s "${fcs_final}.tmp" ]] || fatal "${stage}: final FCS-cleaned FASTA is empty"
    mv -f "${fcs_final}.tmp" "${fcs_final}"
    "${samtools_bin}" faidx "${fcs_final}" || fatal "${stage}: FCS-cleaned FASTA failed validation"
    final_source="${fcs_final}"

        [[ -s "${fcs_final}" ]] || fatal "${stage}: final FCS-cleaned FASTA missing after FCS completion"
        "${samtools_bin}" faidx "${fcs_final}" >/dev/null 2>&1 \
          || fatal "${stage}: final FCS-cleaned FASTA failed checkpoint validation"
        printf '%s\n' "${final_fcs_fp}" > "${final_fcs_fp_file}.tmp"
        mv -f "${final_fcs_fp_file}.tmp" "${final_fcs_fp_file}"
        touch "${final_fcs_ck}"
        final_source="${fcs_final}"
        info "${stage}: final FCS checkpoint written"
      fi
    fi
  else
    info "${stage}: RUN_FINAL_FCS=FALSE — publishing polished assembly without final FCS"
  fi

  # Final assembly cleanup: characterize scaffold lengths, remove <1-kb sequences (configurable),
  # then publish simple sequential identifiers with awk. Filtering occurs after FCS so the published
  # assembly, its length statistics, and downstream annotations all use the same retained sequence set.
  local final_qc_dir="${outdir}/final_fasta_qc"
  local length_filtered="${outdir}/${species_name}.min${FINAL_FASTA_MIN_LENGTH}bp.fasta"
  run_final_fasta_length_qc "${stage}" "${final_source}" "${final_qc_dir}" "${length_filtered}"

  local map_tsv="${outdir}/final.header_map.tsv"
  awk -v prefix="${FINAL_FASTA_ID_PREFIX}" -v map="${map_tsv}.tmp" '
    BEGIN {
      OFS="\t"
      print "new_id", "original_header" > map
    }
    /^>/ {
      original=substr($0,2)
      ++i
      id=sprintf("%s%d", prefix, i)
      print id, original >> map
      print ">" id
      next
    }
    { print }
    END {
      if (i == 0) exit 42
    }
  ' "${length_filtered}" > "${published}.tmp" \
    || fatal "${stage}: failed to rename final FASTA headers sequentially"
  [[ -s "${published}.tmp" ]] || fatal "${stage}: header-renamed final FASTA is empty"
  mv -f "${published}.tmp" "${published}"
  mv -f "${map_tsv}.tmp" "${map_tsv}"
  "${samtools_bin}" faidx "${published}" || fatal "${stage}: published FASTA failed samtools faidx validation"
  sha256sum "${published}" > "${published}.sha256"
  printf 'final_source\t%s\nlength_filtered_fasta\t%s\nminimum_scaffold_length_bp\t%s\nfinal_fasta_qc_dir\t%s\npublished_fasta\t%s\nheader_map\t%s\n' \
    "${final_source}" "${length_filtered}" "${FINAL_FASTA_MIN_LENGTH}" "${final_qc_dir}" "${published}" "${map_tsv}" \
    >> "${outdir}/post_polishing.summary.tsv"
  touch "${outdir}/checkpoints/post_scaffold_refinement.done"
  export final_assembly="${published}"
  info "${stage}: complete → ${published}"
}

# ==================================================================================================
# STEP 6: Shared post-scaffold DeepVariant polishing + final FCS
# Consumes canonical header-sanitized Step-5 scaffolds.
# ==================================================================================================
if is_true "${run_final_asm:-FALSE}"; then
  post_scaffold_input="${dir5}/${base}.scaffolds.fa"
  [[ -s "${post_scaffold_input}" ]] || fatal "Step 6: canonical Step-5 scaffolds missing: ${post_scaffold_input}"
  run_post_scaffold_refinement "6" "${dir6}" "${post_scaffold_input}"
fi

# The prepared FCS input is reproducible from the retained post-polish FASTA.
if is_true "${CLEANUP_COMPLETED_STEPS}" &&
   [[ -f "${dir6}/checkpoints/final_fcs.done" &&
      -s "${dir6}/final.polished.decontam.fa" &&
      -s "${dir6}/${species_name}.fasta" ]]; then
  rm -f -- "${dir6}/final.polished.fcs_input.fa" \
    "${dir6}/final.polished.fcs_input.fa.fai"
fi

# ==================================================================================================
# SHARED FINAL-ASSEMBLY VALIDATION / STATISTICS
# Native Python metrics + BUSCO + HiFi coverage/QC + DIAMOND + BLASTN + BlobToolKit + GenomeScope2 + MultiQC
# The final FCS-cleaned assembly is mandatory; no fallback to earlier assembly generations is allowed.
# ==================================================================================================
run_final_assembly_stats() {
  local stage="${1:?missing stage number}"
  local outdir="${2:?missing statistics output directory}"
  local final_fasta="${3:?missing final assembly FASTA}"

  [[ -s "${final_fasta}" ]] || fatal "Step ${stage}: final FCS-cleaned assembly missing: ${final_fasta}"
  command -v sha256sum >/dev/null 2>&1 || fatal "Step ${stage}: sha256sum is required"

  local assess_dir="${outdir}/assessment"
  local busco_dir="${outdir}/busco"
  local coverage_dir="${outdir}/coverage"
  local blast_dir="${outdir}/blast"
  local diamond_dir="${outdir}/diamond"
  local btk_dir="${outdir}/blobtoolkit"
  local gs_dir="${outdir}/genomescope2"
  local logs_dir="${outdir}/logs"
  local ck_dir="${outdir}/checkpoints"
  local provenance_dir="${outdir}/provenance"
  mkdir -p "${assess_dir}" "${busco_dir}" "${coverage_dir}" "${blast_dir}" "${diamond_dir}" "${btk_dir}" \
           "${gs_dir}" "${logs_dir}" "${ck_dir}" "${provenance_dir}"

  local final_sha
  final_sha="$(sha256sum "${final_fasta}" | awk '{print $1}')"
  printf '%s  %s\n' "${final_sha}" "${final_fasta}" > "${provenance_dir}/final_assembly.sha256"

  stats_read_signature() {
    local f="${1:?missing file}"
    local rp sz mt
    rp="$(readlink -f "${f}")"
    sz="$(stat -c '%s' "${f}")"
    mt="$(stat -c '%Y' "${f}")"
    printf '%s|%s|%s' "${rp}" "${sz}" "${mt}" | sha256sum | awk '{print $1}'
  }

  stats_fp() {
    printf '%s\n' "$@" | sha256sum | awk '{print $1}'
  }

  stats_ck_matches() {
    local tag="${1:?missing checkpoint tag}"
    local fp="${2:?missing fingerprint}"
    [[ -s "${ck_dir}/${tag}.fingerprint.sha256" ]] || return 1
    [[ "$(cat "${ck_dir}/${tag}.fingerprint.sha256")" == "${fp}" ]] || return 1
    [[ -f "${ck_dir}/${tag}.done" ]]
  }

  stats_ck_write() {
    local tag="${1:?missing checkpoint tag}"
    local fp="${2:?missing fingerprint}"
    printf '%s\n' "${fp}" > "${ck_dir}/${tag}.fingerprint.sha256"
    touch "${ck_dir}/${tag}.done"
  }

  stats_ck_drop() {
    local tag="${1:?missing checkpoint tag}"
    rm -f "${ck_dir}/${tag}.done" "${ck_dir}/${tag}.fingerprint.sha256"
  }

  stats_nonempty_files() {
    local f
    for f in "$@"; do
      [[ -s "${f}" ]] || return 1
    done
    return 0
  }

  stats_bam_ready() {
    local bam="${1:?missing BAM}"
    [[ -s "${bam}" && -s "${bam}.bai" && -s "${bam}.csi" ]] || return 1
    "${samtools_bin}" quickcheck -v "${bam}" >/dev/null 2>&1
  }

  stats_hits7_valid() {
    local hits="${1:?missing taxonomy hit table}"
    [[ -s "${hits}" ]] || return 1
    awk -F'\t' '
      NF != 7 { bad=1; exit }
      $1 == "" || $2 == "" ||
      $3 !~ /^[0-9]+$/ || $4 !~ /^[0-9]+$/ ||
      $5 == "" || $6 == "" ||
      ($7 != "" && $7 !~ /^[0-9]+(;[0-9]+)*$/) {
        bad=1
        exit
      }
      END { exit bad ? 1 : 0 }
    ' "${hits}"
  }

  stats_hits_concat_matches() {
    local diamond_hits="${1:?missing DIAMOND normalized table}"
    local blast_hits="${2:?missing BLASTN normalized table}"
    local combined_hits="${3:?missing combined taxonomy table}"
    [[ -s "${diamond_hits}" && -s "${blast_hits}" && -s "${combined_hits}" ]] || return 1
    cmp -s <(cat "${diamond_hits}" "${blast_hits}") "${combined_hits}"
  }

  stats_recover_diamond_hits() {
    local combined_hits="${1:?missing combined taxonomy table}"
    local blast_hits="${2:?missing BLASTN normalized table}"
    local diamond_hits="${3:?missing DIAMOND normalized table}"

    [[ -s "${combined_hits}" && -s "${blast_hits}" ]] || return 1
    stats_hits7_valid "${combined_hits}" || return 1
    stats_hits7_valid "${blast_hits}" || return 1

    local n_total n_blast n_diamond
    n_total="$(wc -l < "${combined_hits}")"
    n_blast="$(wc -l < "${blast_hits}")"
    [[ "${n_total}" =~ ^[0-9]+$ && "${n_blast}" =~ ^[0-9]+$ ]] || return 1
    (( n_total > n_blast )) || return 1
    n_diamond=$(( n_total - n_blast ))

    # v0.7.2 created combined_hits as:
    #   cat "${diamond_hits}" "${blast_hits}" > "${combined_hits}"
    # Require an exact byte-for-byte tail match before reconstructing DIAMOND.
    tail -n "${n_blast}" "${combined_hits}" | cmp -s - "${blast_hits}" || return 1

    local recovered="${diamond_hits}.recovered.$$"
    head -n "${n_diamond}" "${combined_hits}" > "${recovered}"
    stats_hits7_valid "${recovered}" || {
      rm -f "${recovered}"
      return 1
    }

    # Strongest consistency test: recovered DIAMOND + BLASTN must exactly
    # reproduce the preserved combined table.
    cmp -s <(cat "${recovered}" "${blast_hits}") "${combined_hits}" || {
      rm -f "${recovered}"
      return 1
    }

    mv -f "${recovered}" "${diamond_hits}"
    return 0
  }

  # ---------------------- genoflowSTATS final FASTA metrics ----------------------
  info "Step ${stage}a: genoflowSTATS final-assembly statistics"
  local genoflowstats_env_prefix="${MINIFORGE_PREFIX}/envs/${ENV_GENOFLOWSTATS}"
  local genoflowstats_bin="${genoflowstats_env_prefix}/bin/genoflowSTATS"
  local genoflowstats_python="${genoflowstats_env_prefix}/bin/python"
  [[ -x "${genoflowstats_bin}" ]] || fatal "Step ${stage}a: genoflowSTATS executable missing: ${genoflowstats_bin}"
  [[ -x "${genoflowstats_python}" ]] || fatal "Step ${stage}a: Python missing from ${ENV_GENOFLOWSTATS}: ${genoflowstats_python}"
  "${genoflowstats_python}" -c 'import sys,numpy,scipy,pysam,polars,pyarrow,statsmodels,sklearn,matplotlib; assert sys.version_info >= (3,14)' \
    || fatal "Step ${stage}a: ${ENV_GENOFLOWSTATS} integrity/import check failed"

  local assess_genome_size=""
  if [[ "${GENOME_SIZE_EST:-}" =~ ^[0-9]+$ ]] && (( GENOME_SIZE_EST > 0 )); then
    assess_genome_size="${GENOME_SIZE_EST}"
  elif [[ "${GENOME_SIZE:-}" =~ ^[0-9]+$ ]] && (( GENOME_SIZE > 0 )); then
    assess_genome_size="${GENOME_SIZE}"
  fi

  local genoflowstats_version genoflowstats_impl genoflowstats_fasta_fp assembly_metrics_tsv
  genoflowstats_version="$("${genoflowstats_bin}" --version 2>&1 | head -n 1 || true)"
  [[ -n "${genoflowstats_version}" ]] || fatal "Step ${stage}a: unable to obtain genoflowSTATS version"
  genoflowstats_impl="$(genoflowstats_implementation_id "${genoflowstats_python}")" \
    || fatal "Step ${stage}a: install the current genoflowSTATS 0.2.0 implementation with automatic FASTA length and palettes"
  genoflowstats_fasta_fp="$(stats_fp "fasta=${final_sha}" "genoflowstats=${genoflowstats_version}" "implementation=${genoflowstats_impl}" "genome_size=${assess_genome_size}" "gap_min=${GENOFLOWSTATS_GAP_MIN}" "dpi=${GENOFLOWSTATS_PLOT_DPI}" "palette=${GENOFLOWSTATS_PALETTE}")"
  assembly_metrics_tsv="${assess_dir}/${base}.assembly_metrics.tsv"

  if ! stats_ck_matches genoflowstats_fasta "${genoflowstats_fasta_fp}" || [[ ! -s "${assembly_metrics_tsv}" ]]; then
    rm -f "${assess_dir}/${base}."{assembly_metrics.tsv,assembly_metrics.json,scaffolds.tsv,scaffold_cumulative.svg,nx_ngx.svg}
    local -a gf_fasta_args=(fasta --fasta "${final_fasta}" --outdir "${assess_dir}" --prefix "${base}" \
      --gap-min "${GENOFLOWSTATS_GAP_MIN}" --dpi "${GENOFLOWSTATS_PLOT_DPI}" --palette "${GENOFLOWSTATS_PALETTE}")
    [[ -n "${assess_genome_size}" ]] && gf_fasta_args+=(--genome-size "${assess_genome_size}")
    "${genoflowstats_bin}" "${gf_fasta_args[@]}" > "${logs_dir}/genoflowSTATS.fasta.log" 2>&1 \
      || fatal "Step ${stage}a: genoflowSTATS FASTA assessment failed"
    [[ -s "${assembly_metrics_tsv}" ]] || fatal "Step ${stage}a: genoflowSTATS assembly metrics missing"
    stats_ck_write genoflowstats_fasta "${genoflowstats_fasta_fp}"
  fi
  cp -f "${assembly_metrics_tsv}" "${outdir}/assembly.basic_stats.tsv"

  # -------------------------------- BUSCO --------------------------------  # -------------------------------- BUSCO --------------------------------
  info "Step ${stage}b: BUSCO genome completeness"
  case "${BUSCO_DATA_POLICY}" in
    check|update|offline) ;;
    *) fatal "BUSCO_DATA_POLICY must be check, update, or offline (got ${BUSCO_DATA_POLICY})" ;;
  esac
  local busco_run_name="final.assembly.busco"
  local busco_downloads="${BUSCO_DOWNLOADS:-${busco_dir}/downloads}"
  mkdir -p "${busco_downloads}"
  local busco_version
  busco_version="$(env_run "${ENV_BUSCO}" busco --version 2>&1 | head -n 1 || true)"
  local busco_fp
  busco_fp="$(stats_fp "fasta=${final_sha}" "busco=${busco_version}" "lineage=${BUSCO_LINEAGE}" \
    "policy=${BUSCO_DATA_POLICY}" "mode=genome" "mapper=miniprot")"

  local busco_run_dir busco_full_table busco_summary_json
  busco_run_dir="$(find "${busco_dir}/${busco_run_name}" -maxdepth 1 -type d -name 'run_*' -print -quit 2>/dev/null || true)"
  [[ -n "${busco_run_dir}" ]] || busco_run_dir="${busco_dir}/${busco_run_name}"
  busco_full_table="$(find "${busco_run_dir}" -maxdepth 2 -type f -name 'full_table.tsv' -print -quit 2>/dev/null || true)"
  busco_summary_json="$(find "${busco_run_dir}" -maxdepth 2 -type f -name 'short_summary*.json' -print -quit 2>/dev/null || true)"

  local busco_ready="FALSE"
  if stats_ck_matches busco "${busco_fp}" \
     && [[ -n "${busco_full_table}" && -n "${busco_summary_json}" ]] \
     && stats_nonempty_files "${busco_full_table}" "${busco_summary_json}"; then
    busco_ready="TRUE"
    info "Step ${stage}b: BUSCO checkpoint valid — reusing existing BUSCO results"
  fi

  if ! is_true "${busco_ready}"; then
    local BUSCO_PREPARED_LINEAGE=""
    prepare_busco_lineage "${busco_downloads}" "${logs_dir}/busco.download.log" \
      || fatal "Step ${stage}b: BUSCO dataset preparation failed; analysis was not started"
    stats_ck_drop busco
    rm -rf "${busco_dir:?}/${busco_run_name:?}"
    local -a busco_args=(
      -i "${final_fasta}"
      -o "${busco_run_name}"
      --out_path "${busco_dir}"
      -l "${BUSCO_PREPARED_LINEAGE}"
      -m genome
      -c "${nthreads}"
      --download_path "${busco_downloads}"
      --miniprot
      --opt-out-run-stats
      --offline
    )
    busco_bounded_run 0 "${busco_args[@]}" 2>&1 | tee "${logs_dir}/busco.log" \
      || fatal "Step ${stage}b: BUSCO failed; see ${logs_dir}/busco.log"

    busco_run_dir="$(find "${busco_dir}/${busco_run_name}" -maxdepth 1 -type d -name 'run_*' -print -quit 2>/dev/null || true)"
    [[ -n "${busco_run_dir}" ]] || busco_run_dir="${busco_dir}/${busco_run_name}"
    busco_full_table="$(find "${busco_run_dir}" -maxdepth 2 -type f -name 'full_table.tsv' -print -quit 2>/dev/null || true)"
    busco_summary_json="$(find "${busco_run_dir}" -maxdepth 2 -type f -name 'short_summary*.json' -print -quit 2>/dev/null || true)"
    [[ -s "${busco_full_table}" ]] || fatal "Step ${stage}b: BUSCO full_table.tsv missing"
    [[ -s "${busco_summary_json}" ]] || fatal "Step ${stage}b: BUSCO short-summary JSON missing"
    stats_ck_write busco "${busco_fp}"
  fi

  [[ -s "${busco_full_table}" ]] || fatal "Step ${stage}b: BUSCO full_table.tsv missing"
  [[ -s "${busco_summary_json}" ]] || fatal "Step ${stage}b: BUSCO short-summary JSON missing"
  if [[ ! -s "${busco_dir}/busco.metrics.tsv" || "${busco_summary_json}" -nt "${busco_dir}/busco.metrics.tsv" ]]; then
    "${python_bin}" - "${busco_summary_json}" "${busco_dir}/busco.metrics.tsv" <<'PY_BUSCO'
import json,sys
j=json.load(open(sys.argv[1])); out=sys.argv[2]
r=j.get('results',{})
keys=['Complete percentage','Single copy percentage','Multi copy percentage','Fragmented percentage','Missing percentage','n_markers']
with open(out,'w') as f:
    f.write('metric\tvalue\n')
    f.write(f"dataset\t{j.get('lineage_dataset', j.get('lineage', 'NA'))}\n")
    f.write(f"dataset_creation_date\t{j.get('dataset_creation_date','NA')}\n")
    for k in keys: f.write(f"{k}\t{r.get(k,'NA')}\n")
PY_BUSCO
  fi

  # ----------------------- HiFi mapping for coverage -----------------------
  info "Step ${stage}c: HiFi-to-final-assembly mapping"
  [[ -s "${pacbio_hifi}" ]] || fatal "Step ${stage}c: PacBio HiFi reads missing: ${pacbio_hifi}"
  local hifi_sig coverage_bam minimap2_version samtools_version coverage_fp
  hifi_sig="$(stats_read_signature "${pacbio_hifi}")"
  coverage_bam="${coverage_dir}/${base}.final.hifi.bam"
  minimap2_version="$(${minimap2_bin} --version 2>&1 | head -n 1 || true)"
  samtools_version="$(${samtools_bin} --version 2>&1 | head -n 1 || true)"
  coverage_fp="$(stats_fp "fasta=${final_sha}" "reads=${hifi_sig}" "minimap2=${minimap2_version}" \
    "samtools=${samtools_version}" "preset=map-hifi" "secondary=no" "cs=long" \
    "sortmem=${STATS_SAMTOOLS_SORT_MEM}" "compression=${STATS_BAM_COMPRESSION}")"

  local coverage_ready="FALSE"
  local coverage_rebuilt="FALSE"
  if stats_ck_matches coverage "${coverage_fp}" && stats_bam_ready "${coverage_bam}"; then
    coverage_ready="TRUE"
    info "Step ${stage}c: HiFi coverage checkpoint valid — reusing existing BAM and indexes"
  fi

  if ! is_true "${coverage_ready}"; then
    stats_ck_drop coverage
    rm -f "${coverage_bam}" "${coverage_bam}.bai" "${coverage_bam}.csi"
    "${minimap2_bin}" -a -x map-hifi --secondary=no -Y -K "${MINIMAP2_BATCH_SIZE}" -t "${nthreads}" \
      -R "@RG\\tID:${base}.HiFi\\tSM:${base}\\tPL:PACBIO\\tLB:SMRTbell" \
      "${final_fasta}" "${pacbio_hifi}" \
      2> "${logs_dir}/coverage.minimap2.log" \
    | "${samtools_bin}" sort --write-index -@ "${nthreads}" -m "${STATS_SAMTOOLS_SORT_MEM}" \
        -l "${STATS_BAM_COMPRESSION}" -O BAM \
        -o "${coverage_bam}##idx##${coverage_bam}.bai" - \
      2> "${logs_dir}/coverage.samtools-sort.log" \
    || fatal "Step ${stage}c: HiFi mapping/sort failed"
    "${samtools_bin}" quickcheck -v "${coverage_bam}" || fatal "Step ${stage}c: BAM validation failed"
    [[ -s "${coverage_bam}.bai" ]] || fatal "Step ${stage}c: BAI index was not created by samtools sort"
    "${samtools_bin}" index -c -@ "${nthreads}" "${coverage_bam}" || fatal "Step ${stage}c: CSI creation failed"
    [[ -s "${coverage_bam}.bai" && -s "${coverage_bam}.csi" ]] || fatal "Step ${stage}c: BAM indexes missing"
    stats_ck_write coverage "${coverage_fp}"
    coverage_rebuilt="TRUE"
  fi

  # genoflowSTATS consumes the compact COV/MAPQ/ID/GCD/RFS records from
  # samtools stats, so these summaries are required even when the BAM itself is reused.
  if is_true "${coverage_rebuilt}" || [[ ! -s "${coverage_dir}/flagstat.txt" ]]; then
    "${samtools_bin}" flagstat -@ "${nthreads}" "${coverage_bam}" > "${coverage_dir}/flagstat.txt"
  fi
  if is_true "${coverage_rebuilt}" || [[ ! -s "${coverage_dir}/samtools.stats.txt" ]] \
     || ! grep -q '^COV[[:space:]]' "${coverage_dir}/samtools.stats.txt"; then
    "${samtools_bin}" stats -@ "${nthreads}" \
      --coverage "${GENOFLOWSTATS_COV_MIN},${GENOFLOWSTATS_COV_MAX},${GENOFLOWSTATS_COV_STEP}" \
      --ref-seq "${final_fasta}" "${coverage_bam}" > "${coverage_dir}/samtools.stats.txt" \
      || fatal "Step ${stage}c: samtools stats failed"
  fi

  # -------------------- DIAMOND + BLASTN taxonomy evidence --------------------
  # Expensive taxonomy work is checkpointed independently:
  #   taxonomy_chunks -> diamond -> blastn -> taxonomy_merge
  # A legacy v0.7.2 taxonomy checkpoint is automatically adopted when all
  # existing outputs are present and valid.
  info "Step ${stage}d: DIAMOND BLASTX + BLASTN fallback taxonomy evidence"

  local taxonomy_chunks="${diamond_dir}/${base}.taxonomy_chunks.fa"
  local diamond_raw="${diamond_dir}/${base}.diamond.blastx.raw.tsv"
  local diamond_hits="${diamond_dir}/${base}.diamond.blastx.blobtools.tsv"
  local nohit_chunks="${diamond_dir}/${base}.diamond.nohit_chunks.fa"
  local blast_raw="${blast_dir}/core_nt.nohit.raw.tsv"
  local blast_hits="${blast_dir}/core_nt.nohit.blobtools.tsv"
  local combined_hits="${btk_dir}/taxonomy.combined.hits.tsv"
  local diamond_info="${diamond_dir}/database.info.txt"
  local blast_info="${blast_dir}/database.info.txt"

  case "${DIAMOND_SENSITIVITY}" in
    fast|sensitive|more-sensitive|very-sensitive|ultra-sensitive) ;;
    *) fatal "DIAMOND_SENSITIVITY must be fast, sensitive, more-sensitive, very-sensitive, or ultra-sensitive (got ${DIAMOND_SENSITIVITY})" ;;
  esac
  [[ "${BTK_TAXONOMY_CHUNK_BP}" =~ ^[1-9][0-9]*$ ]] \
    || fatal "BTK_TAXONOMY_CHUNK_BP must be a positive integer"

  local diamond_version="disabled"
  local diamond_db_sig="disabled"
  if is_true "${RUN_DIAMOND_STATS}"; then
    [[ -s "${BLASTX_DB_ROOT}" ]] || fatal "Step ${stage}d: DIAMOND DB missing: ${BLASTX_DB_ROOT}"
    diamond_version="$(env_run "${ENV_DIAMOND}" diamond version 2>&1 | head -n 1 || true)"
    if [[ ! -s "${diamond_info}" ]]; then
      env_run "${ENV_DIAMOND}" diamond dbinfo --db "${BLASTX_DB_ROOT}" \
        > "${diamond_info}" 2>&1 || fatal "Step ${stage}d: DIAMOND database unavailable: ${BLASTX_DB_ROOT}"
    fi
    diamond_db_sig="$(stats_read_signature "${BLASTX_DB_ROOT}")"
  else
    printf 'DIAMOND disabled by RUN_DIAMOND_STATS=%s\n' "${RUN_DIAMOND_STATS}" > "${diamond_info}"
  fi

  local blast_version="disabled"
  local blast_db_sig="disabled"
  if is_true "${RUN_BLAST_STATS}"; then
    [[ -x "${blastn_bin}" ]] || fatal "Step ${stage}d: blastn not executable: ${blastn_bin}"
    local blastdbcmd_bin
    blastdbcmd_bin="$(dirname "${blastn_bin}")/blastdbcmd"
    [[ -x "${blastdbcmd_bin}" ]] || fatal "Step ${stage}d: blastdbcmd not found beside blastn"
    if [[ ! -s "${blast_info}" ]]; then
      "${blastdbcmd_bin}" -db "${BLAST_DB_ROOT}" -info > "${blast_info}" 2>&1 \
        || fatal "Step ${stage}d: BLAST database unavailable: ${BLAST_DB_ROOT}"
    fi
    blast_version="$(${blastn_bin} -version 2>&1 | head -n 1 || true)"
    blast_db_sig="$(sha256sum "${blast_info}" | awk '{print $1}')"
  else
    printf 'BLASTN disabled by RUN_BLAST_STATS=%s\n' "${RUN_BLAST_STATS}" > "${blast_info}"
  fi

  local chunks_fp diamond_fp blast_fp taxonomy_merge_fp legacy_taxonomy_fp
  chunks_fp="$(stats_fp "fasta=${final_sha}" "chunk_bp=${BTK_TAXONOMY_CHUNK_BP}")"
  diamond_fp="$(stats_fp \
    "chunks=${chunks_fp}" \
    "diamond_enabled=${RUN_DIAMOND_STATS}" \
    "diamond=${diamond_version}" \
    "diamond_db=${diamond_db_sig}" \
    "diamond_sensitivity=${DIAMOND_SENSITIVITY}" \
    "diamond_evalue=${DIAMOND_EVALUE}" \
    "diamond_targets=${DIAMOND_MAX_TARGET_SEQS}" \
    "diamond_block=${DIAMOND_BLOCK_SIZE}" \
    "diamond_index_chunks=${DIAMOND_INDEX_CHUNKS}")"
  blast_fp="$(stats_fp \
    "diamond=${diamond_fp}" \
    "blast_enabled=${RUN_BLAST_STATS}" \
    "blast=${blast_version}" \
    "blast_db=${blast_db_sig}" \
    "blast_task=${BLAST_TASK}" \
    "blast_evalue=${BLAST_EVALUE}" \
    "blast_targets=${BLAST_MAX_TARGET_SEQS}" \
    "blast_hsps=${BLAST_MAX_HSPS}")"
  taxonomy_merge_fp="$(stats_fp "diamond=${diamond_fp}" "blast=${blast_fp}" "schema=blobtools7-v1")"

  # Exact legacy fingerprint used by v0.7.2. This enables safe migration of a
  # completed taxonomy run into the new fine-grained checkpoints without rerunning searches.
  legacy_taxonomy_fp="$(stats_fp \
    "fasta=${final_sha}" \
    "chunk_bp=${BTK_TAXONOMY_CHUNK_BP}" \
    "diamond_enabled=${RUN_DIAMOND_STATS}" \
    "diamond=${diamond_version}" \
    "diamond_db=${diamond_db_sig}" \
    "diamond_sensitivity=${DIAMOND_SENSITIVITY}" \
    "diamond_evalue=${DIAMOND_EVALUE}" \
    "diamond_targets=${DIAMOND_MAX_TARGET_SEQS}" \
    "diamond_block=${DIAMOND_BLOCK_SIZE}" \
    "diamond_index_chunks=${DIAMOND_INDEX_CHUNKS}" \
    "blast_enabled=${RUN_BLAST_STATS}" \
    "blast=${blast_version}" \
    "blast_db=${blast_db_sig}" \
    "blast_task=${BLAST_TASK}" \
    "blast_evalue=${BLAST_EVALUE}" \
    "blast_targets=${BLAST_MAX_TARGET_SEQS}" \
    "blast_hsps=${BLAST_MAX_HSPS}")"

  # Migration/recovery path for completed v0.7.2 taxonomy searches.
  #
  # The normalized tables and the preserved combined table are the durable
  # downstream products. The raw DIAMOND/BLASTN tables are caches and must not
  # force multi-day recomputation when the normalized evidence is complete.
  #
  # This also recovers a DIAMOND normalized table that was deleted by an
  # interrupted v0.7.3 rerun, provided that:
  #   1. the combined table is valid,
  #   2. the BLASTN table is valid,
  #   3. BLASTN matches the exact tail of combined_hits, and
  #   4. reconstructed DIAMOND + BLASTN exactly reproduces combined_hits.
  local legacy_taxonomy_adopted="FALSE"

  if [[ -s "${combined_hits}" && -s "${blast_hits}" ]] \
     && stats_hits7_valid "${combined_hits}" \
     && stats_hits7_valid "${blast_hits}"; then

    if [[ ! -s "${diamond_hits}" ]] || ! stats_hits7_valid "${diamond_hits}"; then
      if stats_recover_diamond_hits "${combined_hits}" "${blast_hits}" "${diamond_hits}"; then
        info "Step ${stage}d: recovered normalized DIAMOND evidence from preserved combined taxonomy table"
      fi
    fi

    if [[ -s "${diamond_hits}" ]] \
       && stats_hits7_valid "${diamond_hits}" \
       && stats_hits_concat_matches "${diamond_hits}" "${blast_hits}" "${combined_hits}"; then
      info "Step ${stage}d: adopting validated existing DIAMOND/BLASTN taxonomy evidence; expensive searches will not rerun"

      [[ ! -s "${taxonomy_chunks}" ]] || stats_ck_write taxonomy_chunks "${chunks_fp}"
      stats_ck_write diamond "${diamond_fp}"
      stats_ck_write blastn "${blast_fp}"
      stats_ck_write taxonomy_merge "${taxonomy_merge_fp}"
      stats_ck_write taxonomy "${legacy_taxonomy_fp}"
      legacy_taxonomy_adopted="TRUE"
    fi
  fi

  # ---- Taxonomy chunks ----------------------------------------------------
  if stats_ck_matches taxonomy_chunks "${chunks_fp}" && [[ -s "${taxonomy_chunks}" ]]; then
    info "Step ${stage}d1: taxonomy-chunk checkpoint valid — reusing existing chunks"
  else
    stats_ck_drop taxonomy_chunks
    rm -f "${taxonomy_chunks}"
    "${python_bin}" - "${final_fasta}" "${taxonomy_chunks}" "${BTK_TAXONOMY_CHUNK_BP}" <<'PY_BTK_CHUNKS'
import sys
from pathlib import Path
fa = Path(sys.argv[1]); out = Path(sys.argv[2]); chunk = int(sys.argv[3])
def records(path):
    name=None; seq=[]
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name=line[1:].strip().split()[0]; seq=[]
            else:
                seq.append(line.strip())
        if name is not None:
            yield name, "".join(seq)
with out.open("w") as w:
    for name, seq in records(fa):
        for start0 in range(0, len(seq), chunk):
            end0=min(start0+chunk, len(seq))
            frag=seq[start0:end0]
            if not frag:
                continue
            w.write(f">{name}__BTK__{start0+1}__{end0}\n")
            for i in range(0, len(frag), 80):
                w.write(frag[i:i+80] + "\n")
PY_BTK_CHUNKS
    [[ -s "${taxonomy_chunks}" ]] || fatal "Step ${stage}d1: taxonomy chunk FASTA is empty"
    stats_ck_write taxonomy_chunks "${chunks_fp}"
  fi

  # ---- DIAMOND ------------------------------------------------------------
  local diamond_ready="FALSE"
  if stats_ck_matches diamond "${diamond_fp}"; then
    if ! is_true "${RUN_DIAMOND_STATS}"; then
      diamond_ready="TRUE"
      info "Step ${stage}d2: DIAMOND disabled — checkpoint valid"
    elif [[ -s "${diamond_hits}" ]] && stats_hits7_valid "${diamond_hits}"; then
      diamond_ready="TRUE"
      info "Step ${stage}d2: DIAMOND checkpoint valid — reusing normalized BLASTX results"
    fi
  fi

  if ! is_true "${diamond_ready}"; then
    stats_ck_drop diamond
    # Preserve any existing normalized evidence until recomputation succeeds.
    # Only the raw temporary output is replaced up front.
    rm -f "${diamond_raw}"
    if is_true "${RUN_DIAMOND_STATS}"; then
      local diamond_sensitivity_flag="--${DIAMOND_SENSITIVITY}"
      [[ "${DIAMOND_SENSITIVITY}" == "fast" ]] && diamond_sensitivity_flag="--fast"
      env_run "${ENV_DIAMOND}" diamond blastx \
        --query "${taxonomy_chunks}" \
        --db "${BLASTX_DB_ROOT}" \
        --out "${diamond_raw}" \
        --threads "${nthreads}" \
        "${diamond_sensitivity_flag}" \
        --evalue "${DIAMOND_EVALUE}" \
        --max-target-seqs "${DIAMOND_MAX_TARGET_SEQS}" \
        --block-size "${DIAMOND_BLOCK_SIZE}" \
        --index-chunks "${DIAMOND_INDEX_CHUNKS}" \
        --outfmt 6 qseqid sseqid qstart qend evalue bitscore staxids \
        > "${logs_dir}/diamond.stdout.log" 2> "${logs_dir}/diamond.stderr.log" \
        || fatal "Step ${stage}d2: DIAMOND BLASTX failed"

      if [[ -s "${diamond_raw}" ]] && ! awk -F'\t' '$7 != "" && $7 != "0" {ok=1; exit} END{exit !ok}' "${diamond_raw}"; then
        fatal "Step ${stage}d2: DIAMOND hits lack usable staxids; rebuild ${BLASTX_DB_ROOT} with taxonomy metadata"
      fi

      local diamond_hits_tmp="${diamond_hits}.tmp.$$"
      rm -f "${diamond_hits_tmp}"
      "${python_bin}" - "${diamond_raw}" "${diamond_hits_tmp}" <<'PY_BTK_REMAP'
import sys
inp,out=sys.argv[1:]
with open(inp) as r, open(out,"w") as w:
    for line in r:
        f=line.rstrip("\n").split("\t")
        if len(f) < 7: continue
        qid,sid,qs,qe,ev,bs,tax=f[:7]
        if "__BTK__" not in qid: continue
        parent=qid.rsplit("__BTK__",1)[0]
        parts=qid.rsplit("__BTK__",1)[1].split("__")
        if len(parts) != 2: continue
        chunk_start=int(parts[0])
        gqs=chunk_start+int(qs)-1; gqe=chunk_start+int(qe)-1
        w.write("\t".join([parent,sid,str(gqs),str(gqe),ev,bs,tax])+"\n")
PY_BTK_REMAP
      [[ ! -s "${diamond_hits_tmp}" ]] || stats_hits7_valid "${diamond_hits_tmp}" \
        || fatal "Step ${stage}d2: normalized DIAMOND table failed 7-column validation"
      mv -f "${diamond_hits_tmp}" "${diamond_hits}"
    else
      : > "${diamond_raw}"
      : > "${diamond_hits}"
    fi
    stats_ck_write diamond "${diamond_fp}"
  fi

  # ---- BLASTN fallback ----------------------------------------------------
  local blast_ready="FALSE"
  if stats_ck_matches blastn "${blast_fp}"; then
    if ! is_true "${RUN_BLAST_STATS}"; then
      blast_ready="TRUE"
      info "Step ${stage}d3: BLASTN disabled — checkpoint valid"
    elif [[ -s "${blast_hits}" ]] && stats_hits7_valid "${blast_hits}"; then
      blast_ready="TRUE"
      info "Step ${stage}d3: BLASTN checkpoint valid — reusing normalized fallback results"
    fi
  fi

  if ! is_true "${blast_ready}"; then
    stats_ck_drop blastn
    # Keep any existing normalized BLASTN table until a replacement succeeds.
    rm -f "${nohit_chunks}" "${blast_raw}"

    "${python_bin}" - "${taxonomy_chunks}" "${diamond_raw}" "${nohit_chunks}" <<'PY_BTK_NOHIT'
import sys
fa,hits,out=sys.argv[1:]
seen=set()
try:
    with open(hits) as h:
        for line in h:
            if line.strip():
                seen.add(line.split("\t",1)[0])
except FileNotFoundError:
    pass
with open(fa) as r, open(out,"w") as w:
    keep=False
    for line in r:
        if line.startswith(">"):
            qid=line[1:].strip().split()[0]
            keep=qid not in seen
        if keep:
            w.write(line)
PY_BTK_NOHIT

    if is_true "${RUN_BLAST_STATS}" && [[ -s "${nohit_chunks}" ]]; then
      "${blastn_bin}" \
        -task "${BLAST_TASK}" \
        -query "${nohit_chunks}" \
        -db "${BLAST_DB_ROOT}" \
        -outfmt '6 qseqid sseqid qstart qend evalue bitscore staxids' \
        -num_threads "${nthreads}" \
        -max_target_seqs "${BLAST_MAX_TARGET_SEQS}" \
        -max_hsps "${BLAST_MAX_HSPS}" \
        -evalue "${BLAST_EVALUE}" \
        -out "${blast_raw}" \
        > "${logs_dir}/blast.stdout.log" 2> "${logs_dir}/blast.stderr.log" \
        || fatal "Step ${stage}d3: BLASTN fallback failed"

      local blast_hits_tmp="${blast_hits}.tmp.$$"
      rm -f "${blast_hits_tmp}"
      "${python_bin}" - "${blast_raw}" "${blast_hits_tmp}" <<'PY_BTK_REMAP_BLAST'
import sys
inp,out=sys.argv[1:]
with open(inp) as r, open(out,"w") as w:
    for line in r:
        f=line.rstrip("\n").split("\t")
        if len(f) < 7: continue
        qid,sid,qs,qe,ev,bs,tax=f[:7]
        if "__BTK__" not in qid: continue
        parent=qid.rsplit("__BTK__",1)[0]
        parts=qid.rsplit("__BTK__",1)[1].split("__")
        if len(parts) != 2: continue
        chunk_start=int(parts[0])
        gqs=chunk_start+int(qs)-1; gqe=chunk_start+int(qe)-1
        w.write("\t".join([parent,sid,str(gqs),str(gqe),ev,bs,tax])+"\n")
PY_BTK_REMAP_BLAST
      [[ ! -s "${blast_hits_tmp}" ]] || stats_hits7_valid "${blast_hits_tmp}" \
        || fatal "Step ${stage}d3: normalized BLASTN table failed 7-column validation"
      mv -f "${blast_hits_tmp}" "${blast_hits}"
    else
      : > "${blast_raw}"
      : > "${blast_hits}"
    fi
    stats_ck_write blastn "${blast_fp}"
  fi

  # ---- Merge normalized taxonomy evidence --------------------------------
  if stats_ck_matches taxonomy_merge "${taxonomy_merge_fp}" \
     && stats_hits7_valid "${combined_hits}"; then
    info "Step ${stage}d4: taxonomy-merge checkpoint valid — reusing combined hit table"
  elif [[ -s "${combined_hits}" && -s "${diamond_hits}" && -s "${blast_hits}" ]] \
       && stats_hits7_valid "${combined_hits}" \
       && stats_hits_concat_matches "${diamond_hits}" "${blast_hits}" "${combined_hits}"; then
    info "Step ${stage}d4: adopting validated preserved combined taxonomy table"
    stats_ck_write taxonomy_merge "${taxonomy_merge_fp}"
  else
    stats_ck_drop taxonomy_merge
    local combined_tmp="${combined_hits}.tmp.$$"
    cat "${diamond_hits}" "${blast_hits}" > "${combined_tmp}"
    if is_true "${BTK_ENABLE_TAXONOMY}"; then
      stats_hits7_valid "${combined_tmp}" \
        || fatal "Step ${stage}d4: combined taxonomy table is empty or malformed"
    fi
    mv -f "${combined_tmp}" "${combined_hits}"
    stats_ck_write taxonomy_merge "${taxonomy_merge_fp}"
  fi

  # Retain the legacy umbrella checkpoint for backward compatibility.
  stats_ck_write taxonomy "${legacy_taxonomy_fp}"

  if is_true "${BTK_ENABLE_TAXONOMY}"; then
    [[ -s "${combined_hits}" ]] \
      || fatal "Step ${stage}d: BlobToolKit taxonomy is enabled but DIAMOND/BLASTN produced no taxonomic hits"
  fi

  # ----------------------------- BlobToolKit -----------------------------
  #
  # BlobToolKit 4.5.4 notes:
  #   * --hits-cols uses 1-based column numbering.
  #   * an empty terminal staxids field is unsafe because BlobToolKit strips the
  #     trailing tab before splitting the row; convert empty staxids to 0 in a
  #     derived BlobToolKit-only table while preserving the canonical hit table.
  #   * static plots are generated directly by blobtools/blobtk; no browser,
  #     Xvfb, or --driver option is required.
  #   * valid --view values used here are snail, cumulative, and blob.
  #   * blob shape is controlled with --param shape=... (not plotShape=...).
  #   * BlobToolKit 4.5.4 may name shaped blob plots inconsistently, so each
  #     shape is rendered in an isolated temporary directory and then renamed
  #     deterministically.
  info "Step ${stage}e: BlobToolKit BUSCO + coverage + taxonomy dataset and SVG plots"

  local blobdir_path="${btk_dir}/blobdir"
  local plots_out="${btk_dir}/plots"
  local btk_hits="${btk_dir}/taxonomy.combined.hits.blobtools.tsv"
  local blob_stem
  blob_stem="$(basename "${blobdir_path}")"

  # Prepare a BlobToolKit-safe taxonomy table. The canonical combined_hits file
  # remains unchanged. Empty staxids are represented as taxid 0 solely for
  # BlobToolKit parsing.
  if is_true "${BTK_ENABLE_TAXONOMY}"; then
    [[ -s "${combined_hits}" ]] \
      || fatal "Step ${stage}e: taxonomy is enabled but ${combined_hits} is missing or empty"

    local btk_hits_tmp="${btk_hits}.tmp.$$"
    awk -F'\t' -v OFS='\t' '
      NF != 7 {
        printf "Malformed taxonomy row %d: expected 7 columns, found %d\n", NR, NF > "/dev/stderr"
        bad=1
        next
      }
      {
        if ($7 == "") $7="0"
        print
      }
      END { exit bad ? 1 : 0 }
    ' "${combined_hits}" > "${btk_hits_tmp}" \
      || { rm -f "${btk_hits_tmp}"; fatal "Step ${stage}e: failed to prepare BlobToolKit-safe taxonomy hits"; }

    awk -F'\t' '
      NF != 7 || $1 == "" || $2 == "" ||
      $3 !~ /^[0-9]+$/ || $4 !~ /^[0-9]+$/ ||
      $5 == "" || $6 == "" || $7 == "" ||
      $7 !~ /^[0-9]+(;[0-9]+)*$/ {
        bad=1
        exit
      }
      END { exit bad ? 1 : 0 }
    ' "${btk_hits_tmp}" \
      || { rm -f "${btk_hits_tmp}"; fatal "Step ${stage}e: BlobToolKit-safe taxonomy table failed validation"; }

    mv -f "${btk_hits_tmp}" "${btk_hits}"
  else
    rm -f "${btk_hits}"
  fi

  local -a btk_extra_args=(--busco "${busco_full_table}")
  if is_true "${BTK_ENABLE_TAXONOMY}"; then
    [[ -s "${BLOB_TAXDUMP_DIR}/names.dmp" && -s "${BLOB_TAXDUMP_DIR}/nodes.dmp" ]] \
      || fatal "Step ${stage}e: BlobToolKit taxonomy requires names.dmp and nodes.dmp in ${BLOB_TAXDUMP_DIR}"
    btk_extra_args+=(
      --hits "${btk_hits}"
      --hits-cols "1=qseqid,2=sseqid,3=qstart,4=qend,5=evalue,6=bitscore,7=staxids"
      --hit-count "${BTK_HIT_COUNT}"
      --taxdump "${BLOB_TAXDUMP_DIR}"
    )
  fi

  local btk_version btk_dataset_fp btk_plots_fp btk_fp
  btk_version="$(
    singularity exec \
      -B "${working_dir}:${working_dir}" \
      "${sif_blobtoolkit}" \
      blobtools --version 2>&1 | head -n 1 || true
  )"

  local btk_hits_sha="disabled"
  if is_true "${BTK_ENABLE_TAXONOMY}"; then
    btk_hits_sha="$(sha256sum "${btk_hits}" | awk '{print $1}')"
  fi

  btk_dataset_fp="$(stats_fp \
    "fasta=${final_sha}" \
    "bam=${coverage_fp}" \
    "busco=${busco_fp}" \
    "taxonomy_merge=${taxonomy_merge_fp}" \
    "btk_hits=${btk_hits_sha}" \
    "btk=${btk_version}" \
    "taxonomy_enabled=${BTK_ENABLE_TAXONOMY}" \
    "hits_cols=1,2,3,4,5,6,7" \
    "empty_staxids=0" \
    "hitcount=${BTK_HIT_COUNT}" \
    "taxid=${TAXID}")"

  btk_plots_fp="$(stats_fp \
    "dataset=${btk_dataset_fp}" \
    "format=svg" \
    "views=snail,cumulative,blob" \
    "blob_shapes=circle,hex,square" \
    "renderer=blobtk")"

  btk_fp="$(stats_fp "dataset=${btk_dataset_fp}" "plots=${btk_plots_fp}")"

  # Dataset creation is independent of rendering. A plotting failure must never
  # invalidate BUSCO, read mapping, DIAMOND, BLASTN, or a valid BlobDir.
  if stats_ck_matches blobtoolkit_dataset "${btk_dataset_fp}" \
     && [[ -s "${blobdir_path}/meta.json" ]]; then
    info "Step ${stage}e1: BlobToolKit dataset checkpoint valid — reusing BlobDir"
  else
    info "Step ${stage}e1: creating BlobToolKit dataset"
    stats_ck_drop blobtoolkit_dataset
    stats_ck_drop blobtoolkit_plots
    stats_ck_drop blobtoolkit
    rm -rf "${blobdir_path}"
    mkdir -p "${blobdir_path}"

    local -a btk_bind_args=(-B "${working_dir}:${working_dir}")
    if is_true "${BTK_ENABLE_TAXONOMY}"; then
      btk_bind_args+=(-B "${BLOB_TAXDUMP_DIR}:${BLOB_TAXDUMP_DIR}:ro")
    fi

    singularity exec \
      "${btk_bind_args[@]}" \
      "${sif_blobtoolkit}" \
      blobtools add --create \
        --fasta "${final_fasta}" \
        --cov "${coverage_bam}" \
        "${btk_extra_args[@]}" \
        --taxrule bestsum \
        --taxid "${TAXID}" \
        --threads "${nthreads}" \
        "${blobdir_path}" \
      > "${logs_dir}/blobtools.add.log" 2>&1 \
      || fatal "Step ${stage}e1: BlobToolKit dataset creation failed"

    [[ -s "${blobdir_path}/meta.json" ]] \
      || fatal "Step ${stage}e1: BlobToolKit meta.json missing"

    stats_ck_write blobtoolkit_dataset "${btk_dataset_fp}"
  fi

  stats_btk_plots_ready() {
    local dir="${1:?missing BlobToolKit plot directory}"
    local stem="${2:?missing BlobToolKit plot stem}"
    local f
    for f in \
      "${dir}/${stem}.snail.svg" \
      "${dir}/${stem}.cumulative.svg" \
      "${dir}/${stem}.blob.circle.svg" \
      "${dir}/${stem}.blob.hex.svg" \
      "${dir}/${stem}.blob.square.svg"; do
      [[ -s "${f}" ]] || return 1
    done
    return 0
  }

  run_stats_blob_view() {
    local view="${1:?missing BlobToolKit view}"
    local outdir="${2:?missing BlobToolKit plot output directory}"
    shift 2

    mkdir -p "${outdir}"
    singularity exec \
      -B "${working_dir}:${working_dir}" \
      "${sif_blobtoolkit}" \
      blobtools view \
        --plot \
        --view "${view}" \
        --format svg \
        --out "${outdir}" \
        "$@" \
        "${blobdir_path}" \
      >> "${logs_dir}/blobtools.plots.log" 2>&1
  }

  run_stats_blob_shape() {
    local shape="${1:?missing BlobToolKit blob shape}"
    case "${shape}" in
      circle|hex|square) ;;
      *) fatal "Step ${stage}e2: unsupported BlobToolKit blob shape: ${shape}" ;;
    esac

    local shape_tmp="${plots_out}/.tmp.${shape}.$$"
    local rendered target
    rm -rf "${shape_tmp}"
    mkdir -p "${shape_tmp}"

    if ! run_stats_blob_view blob "${shape_tmp}" --param "shape=${shape}"; then
      rm -rf "${shape_tmp}"
      return 1
    fi

    rendered="$(
      find "${shape_tmp}" -maxdepth 1 -type f -name '*.svg' -size +0c -print -quit 2>/dev/null || true
    )"
    [[ -n "${rendered}" ]] || {
      rm -rf "${shape_tmp}"
      return 1
    }

    target="${plots_out}/${blob_stem}.blob.${shape}.svg"
    mv -f "${rendered}" "${target}"
    rm -rf "${shape_tmp}"
    [[ -s "${target}" ]]
  }

  if stats_ck_matches blobtoolkit_plots "${btk_plots_fp}" \
     && stats_btk_plots_ready "${plots_out}" "${blob_stem}"; then
    info "Step ${stage}e2: BlobToolKit SVG plot checkpoint valid — reusing existing plots"
  else
    info "Step ${stage}e2: rendering BlobToolKit SVG plots"
    stats_ck_drop blobtoolkit_plots
    stats_ck_drop blobtoolkit
    rm -rf "${plots_out}"
    mkdir -p "${plots_out}"
    : > "${logs_dir}/blobtools.plots.log"

    local plot_fail=0

    run_stats_blob_view snail "${plots_out}" \
      || { warn "BlobToolKit snail SVG plot failed"; plot_fail=1; }

    run_stats_blob_view cumulative "${plots_out}" \
      || { warn "BlobToolKit cumulative SVG plot failed"; plot_fail=1; }

    run_stats_blob_shape circle \
      || { warn "BlobToolKit circle blob SVG plot failed"; plot_fail=1; }

    run_stats_blob_shape hex \
      || { warn "BlobToolKit hex blob SVG plot failed"; plot_fail=1; }

    run_stats_blob_shape square \
      || { warn "BlobToolKit square blob SVG plot failed"; plot_fail=1; }

    # Remove obsolete files from older implementations if they somehow survived.
    rm -f "${plots_out}/${blob_stem}.grid.svg"

    if (( plot_fail == 0 )) \
       && stats_btk_plots_ready "${plots_out}" "${blob_stem}"; then
      stats_ck_write blobtoolkit_plots "${btk_plots_fp}"
      info "Step ${stage}e2: BlobToolKit SVG plots completed and checkpointed"
    else
      warn "Step ${stage}e2: one or more BlobToolKit SVG plots failed; plot checkpoint not written"
    fi
  fi

  # Umbrella checkpoint is written only when both dataset and plots are valid.
  if stats_ck_matches blobtoolkit_dataset "${btk_dataset_fp}" \
     && stats_ck_matches blobtoolkit_plots "${btk_plots_fp}"; then
    stats_ck_write blobtoolkit "${btk_fp}"
  else
    stats_ck_drop blobtoolkit
  fi

  # Preserve a compact machine-readable taxonomy evidence manifest.
  {
    printf 'key\tvalue\n'
    printf 'taxonomy_chunk_bp\t%s\n' "${BTK_TAXONOMY_CHUNK_BP}"
    printf 'diamond_enabled\t%s\n' "${RUN_DIAMOND_STATS}"
    printf 'diamond_database\t%s\n' "${BLASTX_DB_ROOT}"
    printf 'diamond_hits\t%s\n' "${diamond_hits}"
    printf 'blastn_fallback_enabled\t%s\n' "${RUN_BLAST_STATS}"
    printf 'blastn_database\t%s\n' "${BLAST_DB_ROOT}"
    printf 'blastn_hits\t%s\n' "${blast_hits}"
    printf 'combined_hits\t%s\n' "${combined_hits}"
    printf 'blobtoolkit_hits\t%s\n' "${btk_hits}"
    printf 'blobtoolkit_version\t%s\n' "${btk_version}"
    printf 'blobdir\t%s\n' "${blobdir_path}"
    printf 'plot_format\tsvg\n'
    printf 'blob_shapes\tcircle,hex,square\n'
  } > "${btk_dir}/taxonomy_evidence.manifest.tsv"

  # ----------------------------- GenomeScope2 -----------------------------
  # Intentionally placed after BlobToolKit in the reporting workflow. It profiles the HiFi reads,
  # then compares the independent k-mer genome-size estimate with the final assembly length.
  info "Step ${stage}f: GenomeScope2 HiFi k-mer profiling"
  local gs_meryl_db="${gs_dir}/${base}.reads.k${STATS_GENOMESCOPE_K}.meryl"
  local gs_hist="${gs_dir}/${base}.reads.k${STATS_GENOMESCOPE_K}.hist"
  local gs_name="${base}.final"
  local gs_version meryl_version gs_fp gs_cli
  meryl_version="$(env_run "${ENV_MERQURY}" meryl --version 2>&1 | head -n 1 || true)"
  gs_cli="$(resolve_genomescope2_bin || true)"
  gs_version="$([[ -x "${gs_cli}" ]] && "${gs_cli}" --version 2>&1 | head -n 1 || true)"
  gs_fp="$(stats_fp "reads=${hifi_sig}" "meryl=${meryl_version}" "genomescope=${gs_version}" \
    "k=${STATS_GENOMESCOPE_K}" "ploidy=${STATS_GENOMESCOPE_PLOIDY}" "transform=${STATS_GENOMESCOPE_TRANSFORM_EXP}")"
  if ! stats_ck_matches genomescope2 "${gs_fp}"; then
    rm -rf "${gs_meryl_db}"; rm -f "${gs_hist}"
    activate_env "${ENV_MERQURY}"
    export OMP_NUM_THREADS="${nthreads}"
    meryl k="${STATS_GENOMESCOPE_K}" count output "${gs_meryl_db}" "${pacbio_hifi}" \
      > "${logs_dir}/meryl.genomescope.log" 2>&1 || fatal "Step ${stage}f: meryl count failed"
    meryl histogram "${gs_meryl_db}" > "${gs_hist}" 2> "${logs_dir}/meryl.histogram.log" \
      || fatal "Step ${stage}f: meryl histogram failed"
    deactivate_env
    gs_cli="$(resolve_genomescope2_bin)" || fatal "Step ${stage}f: GenomeScope2 CLI could not be resolved"
    "${gs_cli}" -i "${gs_hist}" -o "${gs_dir}" \
      -k "${STATS_GENOMESCOPE_K}" -p "${STATS_GENOMESCOPE_PLOIDY}" -n "${gs_name}" \
      --transform_exp "${STATS_GENOMESCOPE_TRANSFORM_EXP}" --json_report --fitted_hist --verbose \
      > "${logs_dir}/genomescope2.log" 2>&1 || fatal "Step ${stage}f: GenomeScope2 failed"
    stats_ck_write genomescope2 "${gs_fp}"
  fi
  local gs_json
  gs_json="$(find "${gs_dir}" -maxdepth 1 -type f -name '*report.json' -print -quit 2>/dev/null || true)"
  [[ -s "${gs_json}" ]] || fatal "Step ${stage}f: GenomeScope2 JSON report missing"
  "${python_bin}" - "${gs_json}" "${outdir}/assembly.basic_stats.tsv" "${gs_dir}/genomescope2.metrics.tsv" <<'PY_GS'
import json,sys
j=json.load(open(sys.argv[1])); basic=sys.argv[2]; out=sys.argv[3]
b={}
for line in open(basic):
    if line.startswith('metric'): continue
    k,v=line.rstrip().split('\t',1); b[k]=v
L=float(b.get('total_length', b.get('total_bp', 0)) or 0)
g=j.get('genome_haploid_length',{})
avg=float(g.get('avg',0) or 0)
vals={
 'kcov':j.get('kcov','NA'),
 'genome_haploid_length_avg':g.get('avg','NA'),
 'genome_haploid_length_min':g.get('min','NA'),
 'genome_haploid_length_max':g.get('max','NA'),
 'heterozygous_avg':j.get('heterozygous',{}).get('avg','NA'),
 'homozygous_avg':j.get('homozygous',{}).get('avg','NA'),
 'genome_repeat_length_avg':j.get('genome_repeat_length',{}).get('avg','NA'),
 'genome_unique_length_avg':j.get('genome_unique_length',{}).get('avg','NA'),
 'model_fit_all':j.get('model_fit',{}).get('all','NA'),
 'model_fit_full':j.get('model_fit',{}).get('full','NA'),
 'model_fit_unique':j.get('model_fit',{}).get('unique','NA'),
 'read_error_rate':j.get('read_error_rate','NA'),
 'final_assembly_length':int(L),
 'assembly_to_genomescope_ratio':(L/avg if avg else 'NA'),
 'assembly_vs_genomescope_percent':(((L-avg)/avg)*100 if avg else 'NA'),
}
with open(out,'w') as f:
    f.write('metric\tvalue\n')
    for k,v in vals.items(): f.write(f'{k}\t{v}\n')
PY_GS

  # ---------------- genoflowSTATS alignment QC + longitudinal comparison ----------------
  info "Step ${stage}g: genoflowSTATS alignment QC and cross-stage comparison"

  local final_stats_file="${coverage_dir}/samtools.stats.txt"
  local final_stats_metrics="${assess_dir}/${base}.final.hifi.samtools_metrics.tsv"
  local native_qc_fp
  native_qc_fp="$(stats_fp "coverage=${coverage_fp}" "genoflowstats=${genoflowstats_version}" "implementation=${genoflowstats_impl}" \
    "cov=${GENOFLOWSTATS_COV_MIN},${GENOFLOWSTATS_COV_MAX},${GENOFLOWSTATS_COV_STEP}" "dpi=${GENOFLOWSTATS_PLOT_DPI}" "palette=${GENOFLOWSTATS_PALETTE}")"

  if ! stats_ck_matches genoflowstats_final_stats "${native_qc_fp}" || [[ ! -s "${final_stats_metrics}" ]]; then
    "${genoflowstats_bin}" stats --stats "${final_stats_file}" --reference-fasta "${final_fasta}" --outdir "${assess_dir}" \
      --prefix "${base}.final.hifi" --dpi "${GENOFLOWSTATS_PLOT_DPI}" --palette "${GENOFLOWSTATS_PALETTE}" \
      > "${logs_dir}/genoflowSTATS.samtools.log" 2>&1 \
      || fatal "Step ${stage}g: genoflowSTATS samtools-stats assessment failed"
    [[ -s "${final_stats_metrics}" ]] || fatal "Step ${stage}g: genoflowSTATS final read-support metrics missing"
    stats_ck_write genoflowstats_final_stats "${native_qc_fp}"
  fi

  local comparison_dir="${outdir}/assembly_comparison"
  local comparison_reads_dir="${comparison_dir}/read_support"
  local comparison_manifest="${comparison_dir}/assemblies.tsv"
  mkdir -p "${comparison_dir}" "${comparison_reads_dir}"

  local comparison_fp="disabled"
  if is_true "${GENOFLOWSTATS_COMPARE_ALL_STAGES}"; then
    declare -a compare_labels=()
    declare -a compare_fastas=()
    add_compare_stage() {
      local label="${1:?missing stage label}" fasta="${2:-}"
      [[ -n "${fasta}" && -s "${fasta}" ]] || return 0
      compare_labels+=("${label}")
      compare_fastas+=("${fasta}")
    }

    add_compare_stage "01_hifiasm_primary" "${step1_primary_fasta:-${dir1}/${base}.hifiasm.primary.fa}"
    add_compare_stage "02_fcs_cleaned" "${step2_cleaned_fasta:-${dir2}/clean.fasta}"
    add_compare_stage "03_prepolished" "${step3_polished_fasta:-${dir3}/${base}.prepolish.final.fa}"
    add_compare_stage "04_inspector_corrected" "${step4_output_fasta:-${dir4}/${base}.inspector.corrected.fa}"
    add_compare_stage "05_scaffolded" "${dir5}/${base}.scaffolds.fa"
    add_compare_stage "06_final_refined" "${final_fasta}"

    (( ${#compare_labels[@]} > 0 )) || fatal "Step ${stage}g: no stage assemblies available for genoflowSTATS comparison"

    printf 'name	fasta	samtools_stats
' > "${comparison_manifest}.tmp"
    declare -A stats_by_sha=()
    declare -A fasta_by_sha=()
    local hifi_signature
    hifi_signature="$(stats_read_signature "${pacbio_hifi}")"

    local idx label asm asm_sha stage_stats stage_dir stage_bam stage_fp stage_fp_file
    for ((idx=0; idx<${#compare_labels[@]}; idx++)); do
      label="${compare_labels[$idx]}"
      asm="${compare_fastas[$idx]}"
      asm_sha="$(sha256sum "${asm}" | awk '{print $1}')"
      stage_stats=""

      if is_true "${GENOFLOWSTATS_COMPARE_READ_SUPPORT}"; then
        if [[ -n "${stats_by_sha[$asm_sha]:-}" && -s "${stats_by_sha[$asm_sha]}" ]]; then
          stage_stats="${stats_by_sha[$asm_sha]}"
        elif [[ "${asm_sha}" == "${final_sha}" && -s "${final_stats_file}" ]]; then
          stage_stats="${final_stats_file}"
          stats_by_sha["${asm_sha}"]="${stage_stats}"
        else
          stage_dir="${comparison_reads_dir}/${label}"
          mkdir -p "${stage_dir}"
          stage_stats="${stage_dir}/${label}.samtools.stats.txt"
          stage_bam="${stage_dir}/${label}.hifi.bam"
          stage_fp_file="${stage_dir}/read_support.fingerprint.sha256"
          stage_fp="$(stats_fp "assembly=${asm_sha}" "reads=${hifi_signature}" "minimap2=$(${minimap2_bin} --version 2>/dev/null | head -n1)" \
            "samtools=$(${samtools_bin} --version 2>/dev/null | head -n1)" \
            "cov=${GENOFLOWSTATS_COV_MIN},${GENOFLOWSTATS_COV_MAX},${GENOFLOWSTATS_COV_STEP}")"

          if [[ ! -s "${stage_stats}" || ! -s "${stage_fp_file}" || "$(cat "${stage_fp_file}" 2>/dev/null || true)" != "${stage_fp}" ]]; then
            rm -f "${stage_bam}" "${stage_stats}"
            info "Step ${stage}g: mapping HiFi reads to ${label} for comparative read-support metrics"
            local -a stage_mm2_args=(-a -x map-hifi --secondary=no -Y --cs=long -t "${nthreads}" \
              -R "@RG\tID:${base}.${label}.HiFi\tSM:${base}\tPL:PACBIO\tLB:SMRTbell")
            [[ -n "${MINIMAP2_BATCH_SIZE:-}" ]] && stage_mm2_args+=(-K "${MINIMAP2_BATCH_SIZE}")
            "${minimap2_bin}" "${stage_mm2_args[@]}" "${asm}" "${pacbio_hifi}" 2> "${stage_dir}/minimap2.log" \
            | "${samtools_bin}" sort -@ "${nthreads}" -m "${STATS_SAMTOOLS_SORT_MEM}" -l 1 -O BAM \
                -o "${stage_bam}" - 2> "${stage_dir}/samtools-sort.log" \
              || fatal "Step ${stage}g: comparative HiFi mapping failed for ${label}"
            "${samtools_bin}" quickcheck -v "${stage_bam}" \
              || fatal "Step ${stage}g: comparative BAM validation failed for ${label}"
            "${samtools_bin}" stats -@ "${nthreads}" \
              --coverage "${GENOFLOWSTATS_COV_MIN},${GENOFLOWSTATS_COV_MAX},${GENOFLOWSTATS_COV_STEP}" \
              --ref-seq "${asm}" "${stage_bam}" > "${stage_stats}" \
              || fatal "Step ${stage}g: samtools stats failed for ${label}"
            printf '%s\n' "${stage_fp}" > "${stage_fp_file}"
            if ! is_true "${GENOFLOWSTATS_KEEP_STAGE_BAMS}"; then
              rm -f "${stage_bam}" "${stage_bam}.bai" "${stage_bam}.csi"
            fi
          fi
          stats_by_sha["${asm_sha}"]="${stage_stats}"
        fi
      fi

      printf '%s	%s	%s
' "${label}" "${asm}" "${stage_stats}" >> "${comparison_manifest}.tmp"
      fasta_by_sha["${asm_sha}"]="${asm}"
    done
    mv -f "${comparison_manifest}.tmp" "${comparison_manifest}"

    comparison_fp="$(sha256sum "${comparison_manifest}" | awk '{print $1}')"
    comparison_fp="$(stats_fp "manifest=${comparison_fp}" "genoflowstats=${genoflowstats_version}" "implementation=${genoflowstats_impl}" \
      "genome_size=${assess_genome_size}" "gap_min=${GENOFLOWSTATS_GAP_MIN}" "dpi=${GENOFLOWSTATS_PLOT_DPI}" "palette=${GENOFLOWSTATS_PALETTE}")"

    if ! stats_ck_matches genoflowstats_compare "${comparison_fp}" \
       || [[ ! -s "${comparison_dir}/comparison.tsv" || ! -s "${comparison_dir}/comparison.long.tsv" ]]; then
      rm -f "${comparison_dir}/comparison.tsv" "${comparison_dir}/comparison.long.tsv" \
            "${comparison_dir}/comparison.json" "${comparison_dir}/comparison.nx_ngx.svg" "${comparison_dir}/comparison.heatmap.svg"
      local -a compare_args=(compare --manifest "${comparison_manifest}" --outdir "${comparison_dir}" \
        --threads "${GENOFLOWSTATS_COMPARE_THREADS}" --gap-min "${GENOFLOWSTATS_GAP_MIN}" --dpi "${GENOFLOWSTATS_PLOT_DPI}" --palette "${GENOFLOWSTATS_PALETTE}")
      [[ -n "${assess_genome_size}" ]] && compare_args+=(--genome-size "${assess_genome_size}")
      "${genoflowstats_bin}" "${compare_args[@]}" > "${logs_dir}/genoflowSTATS.compare.log" 2>&1 \
        || fatal "Step ${stage}g: genoflowSTATS cross-stage comparison failed"
      stats_ck_write genoflowstats_compare "${comparison_fp}"
    fi

    # Explicit per-step deltas make the direction and magnitude of change visible,
    # while the original wide and long genoflowSTATS tables remain authoritative.
    "${genoflowstats_python}" - "${comparison_dir}/comparison.long.tsv" "${comparison_dir}/stage_improvements.tsv" <<'PY_STAGE_DELTAS'
import csv, math, sys
src, dst = sys.argv[1], sys.argv[2]
rows=list(csv.DictReader(open(src, newline='', encoding='utf-8'), delimiter='	'))
order=[]
by={}
for r in rows:
    a=r['assembly']
    if a not in by:
        by[a]={}
        order.append(a)
    try:
        by[a][r['metric']]=(float(r['value']),r.get('unit',''),r.get('direction',''))
    except Exception:
        pass
with open(dst,'w',newline='',encoding='utf-8') as fh:
    w=csv.writer(fh,delimiter='\t',lineterminator='\n')
    w.writerow(['stage','previous_stage','metric','value','previous_value','delta','percent_change','unit','direction'])
    for i,a in enumerate(order):
        prev=order[i-1] if i else ''
        metrics=sorted(by[a])
        for m in metrics:
            v,unit,direction=by[a][m]
            pv=by.get(prev,{}).get(m,(math.nan,'',''))[0] if prev else math.nan
            delta=v-pv if math.isfinite(pv) else math.nan
            pct=(delta/abs(pv)*100.0) if math.isfinite(pv) and pv!=0 else math.nan
            w.writerow([a,prev,m,v,'' if not math.isfinite(pv) else pv,
                        '' if not math.isfinite(delta) else delta,
                        '' if not math.isfinite(pct) else pct,unit,direction])
PY_STAGE_DELTAS
    [[ -s "${comparison_dir}/stage_improvements.tsv" ]] \
      || fatal "Step ${stage}g: stage-improvement table was not generated"
  fi

  # -------------------------------- MultiQC --------------------------------  # -------------------------------- MultiQC --------------------------------
  info "Step ${stage}h: MultiQC consolidated report"
  local multiqc_fp
  multiqc_fp="$(stats_fp "fasta=${final_sha}" "assessment=${genoflowstats_fasta_fp}" "busco=${busco_fp}" \
    "coverage=${coverage_fp}" "native_qc=${native_qc_fp}" "comparison=${comparison_fp}" "blobtoolkit=${btk_fp}" "genomescope=${gs_fp}")"
  if ! stats_ck_matches multiqc "${multiqc_fp}" || [[ ! -s "${outdir}/multiqc_report.html" ]]; then
    rm -f "${outdir}/multiqc_report.html"
    multiqc "${outdir}" -o "${outdir}" -p -ip -v > "${logs_dir}/multiqc.log" 2>&1 \
      || fatal "Step ${stage}h: MultiQC failed"
    [[ -s "${outdir}/multiqc_report.html" ]] || fatal "Step ${stage}h: MultiQC report missing"
    stats_ck_write multiqc "${multiqc_fp}"
  fi

  {
    printf 'key\tvalue\n'
    printf 'final_fasta\t%s\n' "${final_fasta}"
    printf 'final_fasta_sha256\t%s\n' "${final_sha}"
    printf 'assembly_metrics\t%s\n' "${assembly_metrics_tsv}"
    printf 'busco_metrics\t%s\n' "${busco_dir}/busco.metrics.tsv"
    printf 'coverage_bam\t%s\n' "${coverage_bam}"
    printf 'diamond_hits\t%s\n' "${diamond_hits}"
    printf 'blastn_fallback_hits\t%s\n' "${blast_hits}"
    printf 'combined_taxonomy_hits\t%s\n' "${combined_hits}"
    printf 'blobtoolkit_dataset\t%s\n' "${blobdir_path}"
    printf 'genomescope_metrics\t%s\n' "${gs_dir}/genomescope2.metrics.tsv"
    printf 'samtools_metrics\t%s\n' "${final_stats_metrics}"
    printf 'assembly_comparison\t%s\n' "${comparison_dir}/comparison.tsv"
    printf 'stage_improvements\t%s\n' "${comparison_dir}/stage_improvements.tsv"
    printf 'multiqc_report\t%s\n' "${outdir}/multiqc_report.html"
  } > "${outdir}/final_statistics.manifest.tsv"

  info "Step ${stage}: final assembly validation/statistics complete"
}

# ==================================================================================================
# STEP 7: Final assembly validation and reporting
# Consumes only the canonical published Step-6 assembly.
# ==================================================================================================
if is_true "${run_stats:-FALSE}"; then
  run_final_assembly_stats "7" "${dir7}" "${dir6}/${species_name}.fasta"
fi

# ============================================================================================================================
# STEP 8: Repeats — RepeatModeler2, TEtrimmer, final-library construction,
# Earl Grey annotation, enhanced plotting, and non-repeat intervals
# ============================================================================================================================
if is_true "${run_repeats:-FALSE}"; then
  echo "[INFO] Step 8: Repeat annotation — starting"

  repeats_target="${species_fasta:-${dir6}/${species_name}.fasta}"
  [[ -s "${repeats_target}" ]] || fatal "Final assembly required for repeat annotation but not found: ${repeats_target}"

  PROJECT_ROOT="${REPEATS_PROJECT_ROOT:-${repeats_dir}}"
  GENOME="${REPEATS_GENOME:-${repeats_target}}"

  # Re-resolve the authoritative repeat taxonomy at the repeat stage so
  # Step 8/9 can be resumed safely in a fresh shell.  The same taxon is used
  # for Earl Grey and Dfam/FamDB by design.
  if [[ -z "${REPEAT_TAXON:-}" ]]; then
    if [[ -n "${REPEATS_SPECIES:-}" ]]; then
      REPEAT_TAXON="${REPEATS_SPECIES}"
    elif [[ -n "${DFAM_TAXON:-}" ]]; then
      REPEAT_TAXON="${DFAM_TAXON}"
    elif [[ -n "${TAXID:-}" ]]; then
      REPEAT_TAXON="${TAXID}"
    elif is_true "${REPEAT_TAXON_USE_GENOME_SPECIES:-TRUE}"; then
      REPEAT_TAXON="${species_name_raw}"
    else
      fatal "REPEAT_TAXON is required when no repeat taxonomy can be inferred"
    fi
  fi

  REPEATS_SPECIES="${REPEAT_TAXON}"
  DFAM_TAXON="${REPEAT_TAXON}"

  SPECIES="${REPEATS_SPECIES}"
  THREADS="${REPEATS_THREADS}"

  echo "[INFO] Repeat taxonomy: REPEAT_TAXON='${REPEAT_TAXON}' REPEATS_SPECIES='${REPEATS_SPECIES}' DFAM_TAXON='${DFAM_TAXON}'"

  echo "[INFO] Repeat taxonomy: REPEATS_SPECIES='${REPEATS_SPECIES}' DFAM_TAXON='${DFAM_TAXON}' (reuse_genome_species=${REPEAT_TAXON_USE_GENOME_SPECIES})"

  REPEATS_LIB_DIR="${REPEATS_LIB_DIR:-${PROJECT_ROOT}/Libraries}"
  DFAM_FAMDB_DIR="${DFAM_FAMDB_DIR:-${REPEATS_LIB_DIR}/famdb}"
  REPBASE_DIR="${REPBASE_DIR:-${REPEATS_LIB_DIR}}"
  # REPBASE_EMBL may be set explicitly. If empty, the newest file matching
  # REPBASE_EMBL_PATTERN, RepBase*.embl.gz, or RepBase*.embl in REPBASE_DIR is used.

  if [[ "$ENV_RUNNER" == "mamba" ]]; then
    if [[ -z "${MAMBA_SH:-}" || ! -s "$MAMBA_SH" ]]; then
      echo "[ERROR] Mamba initialization script not found." >&2
      echo "[ERROR] Set MAMBA_SH globally or set ENV_RUNNER=conda." >&2
      exit 1
    fi
  fi

env_run_clean_r() {
  # Run a command in a conda environment while removing R variables inherited
  # from the login shell or other conda environments. This prevents Earl Grey's
  # internal Rscript calls from using the wrong R library/configuration.
  local env="${1:-}"
  shift || true
  [[ -n "${env}" ]] || error "env_run_clean_r: missing environment name"
  OPENBLAS_NUM_THREADS=1 \
  OMP_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 \
  NUMEXPR_NUM_THREADS=1 \
  VECLIB_MAXIMUM_THREADS=1 \
  BLIS_NUM_THREADS=1 \
  env \
    -u R_ARCH \
    -u R_ARCHS \
    -u R_HOME \
    -u R_LIBS \
    -u R_LIBS_USER \
    -u R_LIBS_SITE \
    -u R_PROFILE_USER \
    -u R_ENVIRON_USER \
    "${ENV_RUNNER}" run -n "${env}" "$@"
}

EARLGREY_ANNOTATION_CLI=""
EARLGREY_DIVERGENCE_CALC=""
EARLGREY_ENV_PREFIX=""
EARLGREY_SHARE_DIR=""
EARLGREY_SCRIPT_DIR=""
EARLGREY_DETECTED_VERSION=""

version_ge() {
  # Return success when version $1 is greater than or equal to version $2.
  # GNU sort -V is already a pipeline requirement on supported Linux systems.
  local have="${1:?}"
  local need="${2:?}"
  [[ "$(printf '%s\n%s\n' "$need" "$have" | sort -V | head -n 1)" == "$need" ]]
}

earlgrey_resolve_installation() {
  # Resolve the active Earl Grey package layout from the selected environment.
  # Never construct share/earlgrey-<version>-<build> paths from a hard-coded
  # version string: conda build directories change between releases.
  [[ -n "$EARLGREY_ENV_PREFIX" &&
     -d "$EARLGREY_SHARE_DIR" &&
     -d "$EARLGREY_SCRIPT_DIR" ]] && return 0

  local layout=""
  local key=""
  local value=""

  # shellcheck disable=SC2016
  layout="$(env_run_clean_r "$ENV_EARLGREY" bash -lc '
    set -Eeuo pipefail

    prefix="${CONDA_PREFIX:?CONDA_PREFIX is not set}"
    share_root="${prefix}/share"
    [[ -d "$share_root" ]] || {
      echo "Earl Grey share directory does not exist: $share_root" >&2
      exit 1
    }

    shopt -s nullglob
    candidates=("${share_root}"/earlgrey-*)
    valid=()

    for candidate in "${candidates[@]}"; do
      [[ -d "${candidate}/scripts" ]] || continue
      [[ -s "${candidate}/scripts/divergenceCalc/divergence_calc.py" ]] || continue
      [[ -f "${candidate}/scripts/rcMergeRepeatsLoose" ]] || continue
      valid+=("$candidate")
    done

    [[ "${#valid[@]}" -ge 1 ]] || {
      echo "No complete Earl Grey installation was found under ${share_root}/earlgrey-*" >&2
      exit 1
    }

    # A conda environment should normally contain one package directory. If
    # stale directories remain after an update, select the highest versioned
    # complete installation rather than an obsolete hard-coded path.
    share_dir="$(
      printf "%s\n" "${valid[@]}" |
        sort -V |
        tail -n 1
    )"

    script_dir="${share_dir}/scripts"
    package_tag="$(basename "$share_dir")"
    package_tag="${package_tag#earlgrey-}"
    version="${package_tag%%-*}"

    printf "prefix=%s\n" "$prefix"
    printf "share=%s\n" "$share_dir"
    printf "scripts=%s\n" "$script_dir"
    printf "version=%s\n" "$version"
  ')" || fatal "Could not resolve the active Earl Grey installation in environment: $ENV_EARLGREY"

  while IFS='=' read -r key value; do
    case "$key" in
      prefix)  EARLGREY_ENV_PREFIX="$value" ;;
      share)   EARLGREY_SHARE_DIR="$value" ;;
      scripts) EARLGREY_SCRIPT_DIR="$value" ;;
      version) EARLGREY_DETECTED_VERSION="$value" ;;
    esac
  done <<< "$layout"

  [[ -d "$EARLGREY_ENV_PREFIX" ]] ||
    fatal "Resolved Earl Grey environment prefix is invalid: $EARLGREY_ENV_PREFIX"
  [[ -d "$EARLGREY_SHARE_DIR" ]] ||
    fatal "Resolved Earl Grey package directory is invalid: $EARLGREY_SHARE_DIR"
  [[ -d "$EARLGREY_SCRIPT_DIR" ]] ||
    fatal "Resolved Earl Grey script directory is invalid: $EARLGREY_SCRIPT_DIR"
  [[ -n "$EARLGREY_DETECTED_VERSION" ]] ||
    fatal "Could not determine the installed Earl Grey version from: $EARLGREY_SHARE_DIR"

  require_file \
    "${EARLGREY_SCRIPT_DIR}/divergenceCalc/divergence_calc.py" \
    "Earl Grey divergence calculator"
  [[ -f "${EARLGREY_SCRIPT_DIR}/rcMergeRepeatsLoose" ]] ||
    fatal "Earl Grey helper missing: ${EARLGREY_SCRIPT_DIR}/rcMergeRepeatsLoose"

  if ! version_ge "$EARLGREY_DETECTED_VERSION" "$EARLGREY_REQUIRED_VERSION"; then
    fatal \
      "Earl Grey ${EARLGREY_DETECTED_VERSION} is older than the required minimum " \
      "${EARLGREY_REQUIRED_VERSION}"
  fi
}

earlgrey_prepare_divergence_calc() {
  local source_py=""
  local patched_py="${DIR_LOGS}/divergence_calc.runtime.py"
  local patched_tmp="${patched_py}.tmp.$$"
  local syntax_log="${DIR_LOGS}/04_earlgrey_divergence_calc.runtime.syntax.log"
  local patch_audit="${DIR_LOGS}/04_earlgrey_divergence_calc.runtime_patch.txt"
  local split_bug_count=0
  local lookup_bug_count=0

  earlgrey_resolve_installation
  source_py="${EARLGREY_SCRIPT_DIR}/divergenceCalc/divergence_calc.py"

  [[ -s "$source_py" ]] ||
    fatal "Resolved Earl Grey divergence calculator is missing or empty: $source_py"

  split_bug_count="$(
    grep -Fc \
      'file_name = os.path.join(temp_dir, "split_library", repeat_name + ".fasta")' \
      "$source_py" || true
  )"
  lookup_bug_count="$(
    grep -Fc \
      'subject_path = os.path.join(temp_dir, "split_library", repeat_family + ".fasta")' \
      "$source_py" || true
  )"

  if [[ "$split_bug_count" -eq 0 && "$lookup_bug_count" -eq 0 ]]; then
    EARLGREY_DIVERGENCE_CALC="$source_py"
    {
      echo "mode=original-valid"
      echo "installed_version=$EARLGREY_DETECTED_VERSION"
      echo "source_py=$source_py"
      echo "resolved_py=$EARLGREY_DIVERGENCE_CALC"
      echo "repair_safe_repeat_filenames=0"
      sha256sum "$source_py"
    } > "$patch_audit"
    return 0
  fi

  [[ "$split_bug_count" -eq 1 && "$lookup_bug_count" -eq 1 ]] || {
    fatal \
      "Refusing to patch Earl Grey divergence calculator with unexpected " \
      "split/lookup signature counts (${split_bug_count}/${lookup_bug_count}): $source_py"
  }

  is_true "$EARLGREY_RUNTIME_PATCH_ANNOTATION_CLI" || {
    fatal \
      "Earl Grey divergence filename repair is required but " \
      "EARLGREY_RUNTIME_PATCH_ANNOTATION_CLI is disabled"
  }

  # Compatibility is determined by the exact source signatures above and by
  # the dynamically resolved installation path, not by an Earl Grey version
  # string. This keeps the repair safe across compatible future releases.
  [[ "$source_py" == "${EARLGREY_SCRIPT_DIR}/divergenceCalc/divergence_calc.py" ]] ||
    fatal "Resolved divergence calculator escaped active Earl Grey script directory: $source_py"

  [[ "$(grep -Ec '^import os[[:space:]]*$' "$source_py")" -eq 1 ]] || {
    fatal "Refusing to patch divergence_calc.py without exactly one import os statement: $source_py"
  }
  [[ "$(grep -Ec '^def splitter\(in_seq, temp_dir\):[[:space:]]*$' "$source_py")" -eq 1 ]] || {
    fatal "Refusing to patch divergence_calc.py without exactly one splitter function: $source_py"
  }

  if ! awk '
    BEGIN { imported = 0; helper = 0; split_fixed = 0; lookup_fixed = 0 }
    /^import os[[:space:]]*$/ {
      print
      print "import hashlib"
      imported++
      next
    }
    /^def splitter\(in_seq, temp_dir\):[[:space:]]*$/ {
      print "def repeat_library_filename(repeat_name):"
      print "    return hashlib.sha256(repeat_name.encode(\"utf-8\")).hexdigest() + \".fasta\""
      print ""
      print
      helper++
      next
    }
    /file_name = os\.path\.join\(temp_dir, "split_library", repeat_name \+ "\.fasta"\)/ {
      match($0, /^[[:space:]]*/)
      print substr($0, RSTART, RLENGTH) \
        "file_name = os.path.join(temp_dir, \"split_library\", repeat_library_filename(repeat_name))"
      split_fixed++
      next
    }
    /subject_path = os\.path\.join\(temp_dir, "split_library", repeat_family \+ "\.fasta"\)/ {
      match($0, /^[[:space:]]*/)
      print substr($0, RSTART, RLENGTH) \
        "subject_path = os.path.join(temp_dir, \"split_library\", repeat_library_filename(repeat_family))"
      lookup_fixed++
      next
    }
    { print }
    END {
      if (imported != 1 || helper != 1 || split_fixed != 1 || lookup_fixed != 1) exit 42
    }
  ' "$source_py" > "$patched_tmp"; then
    rm -f "$patched_tmp"
    fatal "Failed to create runtime-patched Earl Grey divergence calculator"
  fi

  if ! env_run_clean_r "$ENV_EARLGREY" python -c \
    'import pathlib, sys; p = pathlib.Path(sys.argv[1]); compile(p.read_text(), str(p), "exec")' \
    "$patched_tmp" > "$syntax_log" 2>&1; then
    rm -f "$patched_tmp"
    fatal \
      "Runtime-patched Earl Grey divergence calculator failed Python syntax validation: " \
      "$syntax_log"
  fi

  chmod 0644 "$patched_tmp"
  mv -f "$patched_tmp" "$patched_py"
  EARLGREY_DIVERGENCE_CALC="$patched_py"

  {
    echo "mode=runtime-safe-repeat-filename-repair"
    echo "installed_version=$EARLGREY_DETECTED_VERSION"
    echo "source_py=$source_py"
    echo "resolved_py=$EARLGREY_DIVERGENCE_CALC"
    echo "repair_safe_repeat_filenames=1"
    echo "filename_mapping=sha256-lowercase-repeat-family"
    echo "syntax_log=$syntax_log"
    sha256sum "$source_py" "$EARLGREY_DIVERGENCE_CALC"
  } > "$patch_audit"

  echo \
    "[WARN] Earl Grey ${EARLGREY_DETECTED_VERSION} divergence calculator contains " \
    "the known unsafe repeat-filename implementation." >&2
  echo "[WARN] Using validated runtime copy: $EARLGREY_DIVERGENCE_CALC" >&2
}

earlgrey_prepare_annotation_cli() {
  local source_cli=""
  local patched_cli="${DIR_LOGS}/earlGreyAnnotationOnly.runtime.sh"
  local patched_tmp="${patched_cli}.tmp.$$"
  local original_syntax_log="${DIR_LOGS}/04_earlgrey_cli.original.bash_n.log"
  local patched_syntax_log="${DIR_LOGS}/04_earlgrey_cli.runtime.bash_n.log"
  local patch_audit="${DIR_LOGS}/04_earlgrey_cli.runtime_patch.txt"
  local marker_count=0
  local buggy_library_condition_count=0
  local divergence_call_count=0
  local routed_divergence_call_count=0
  local legacy_divcalc_tmp_count=0
  local process_divcalc_tmp_count=0
  local heliano_glob_bug_count=0
  local heliano_glob_fixed_count=0
  local script_dir_count=0
  local routed_script_dir_count=0
  local script_dir_ref_count=0
  local needs_brace_repair=0
  local needs_library_repair=0
  local needs_divergence_routing=0
  local needs_divcalc_tmp_repair=0
  local needs_heliano_resume_repair=0
  local needs_script_dir_routing=0
  local patch_mode=""

  earlgrey_resolve_installation

  source_cli="$(env_run_clean_r "$ENV_EARLGREY" bash -lc 'command -v earlGreyAnnotationOnly')" || {
    fatal "Could not resolve earlGreyAnnotationOnly in environment: $ENV_EARLGREY"
  }
  source_cli="${source_cli%%$'\n'*}"
  [[ -s "$source_cli" ]] || {
    fatal "Resolved Earl Grey annotation-only CLI is missing or empty: $source_cli"
  }

  if ! env_run_clean_r "$ENV_EARLGREY" bash -n "$source_cli" > "$original_syntax_log" 2>&1; then
    needs_brace_repair=1
  fi

  buggy_library_condition_count="$(
    grep -Ec \
      '^[[:space:]]*if \[ -z "\$RepSpec" \] && \[ -z "\$startCust" \]; then[[:space:]]*$' \
      "$source_cli" || true
  )"
  divergence_call_count="$(
    grep -Fc 'python ${SCRIPT_DIR}/divergenceCalc/divergence_calc.py' "$source_cli" || true
  )"
  routed_divergence_call_count="$(
    grep -Fc \
      'python "${EARLGREY_DIVERGENCE_CALC:-${SCRIPT_DIR}/divergenceCalc/divergence_calc.py}"' \
      "$source_cli" || true
  )"
  legacy_divcalc_tmp_count="$(
    grep -Fc 'divcalc_tmp="/tmp/egdiv_${species}"' "$source_cli" || true
  )"
  process_divcalc_tmp_count="$(
    grep -Fc 'divcalc_tmp="${TMPDIR:-/tmp}/egdiv_${species}_$$"' "$source_cli" || true
  )"
  heliano_glob_bug_count="$(
    grep -Fc \
      'if [ ! -s ${OUTDIR}/${species}_heliano/*/RC.representative.gff ]; then' \
      "$source_cli" || true
  )"
  heliano_glob_fixed_count="$(
    grep -Fc \
      'if ! compgen -G "${OUTDIR}/${species}_heliano/*/RC.representative.gff" > /dev/null; then' \
      "$source_cli" || true
  )"

  # Route Earl Grey helper calls through the dynamically resolved package
  # script directory. This specifically prevents a conda update from leaving
  # runtime calls pinned to an obsolete path such as share/earlgrey-7.3.0-1.
  script_dir_count="$(
    grep -Ec '^[[:space:]]*SCRIPT_DIR=' "$source_cli" || true
  )"
  routed_script_dir_count="$(
    grep -Fc 'SCRIPT_DIR="${EARLGREY_SCRIPT_DIR:-${SCRIPT_DIR}}"' "$source_cli" || true
  )"
  script_dir_ref_count="$(
    grep -Fc '${SCRIPT_DIR}' "$source_cli" || true
  )"

  if [[ "$script_dir_count" -eq 1 && "$routed_script_dir_count" -eq 0 ]]; then
    needs_script_dir_routing=1
  elif [[ "$routed_script_dir_count" -gt 1 || "$script_dir_count" -gt 2 ]]; then
    fatal \
      "Refusing to patch Earl Grey CLI with unexpected SCRIPT_DIR assignment counts " \
      "(${script_dir_count}/${routed_script_dir_count}): $source_cli"
  elif [[ "$script_dir_count" -eq 0 && "$script_dir_ref_count" -gt 0 ]]; then
    fatal \
      "Earl Grey CLI references SCRIPT_DIR but no SCRIPT_DIR assignment could be found: " \
      "$source_cli"
  fi

  if [[ "$buggy_library_condition_count" -eq 1 ]]; then
    needs_library_repair=1
  elif [[ "$buggy_library_condition_count" -gt 1 ]]; then
    fatal \
      "Refusing to patch Earl Grey CLI with an unexpected calcDivRL " \
      "library-selection condition: $source_cli"
  fi

  if [[ "$divergence_call_count" -eq 2 && "$routed_divergence_call_count" -eq 0 ]]; then
    needs_divergence_routing=1
  elif [[ "$divergence_call_count" -ne 0 || "$routed_divergence_call_count" -ne 2 ]]; then
    fatal \
      "Refusing to patch Earl Grey CLI with unexpected divergence command counts " \
      "(${divergence_call_count}/${routed_divergence_call_count}): $source_cli"
  fi

  if [[ "$legacy_divcalc_tmp_count" -eq 1 && "$process_divcalc_tmp_count" -eq 0 ]]; then
    needs_divcalc_tmp_repair=1
  elif [[ "$legacy_divcalc_tmp_count" -ne 0 || "$process_divcalc_tmp_count" -ne 1 ]]; then
    fatal \
      "Refusing to patch Earl Grey CLI with unexpected divergence temp-directory counts " \
      "(${legacy_divcalc_tmp_count}/${process_divcalc_tmp_count}): $source_cli"
  fi

  if [[ "$heliano_glob_bug_count" -eq 1 && "$heliano_glob_fixed_count" -eq 0 ]]; then
    needs_heliano_resume_repair=1
  elif [[ "$heliano_glob_bug_count" -ne 0 || "$heliano_glob_fixed_count" -ne 1 ]]; then
    fatal \
      "Refusing to patch Earl Grey CLI with unexpected HELIANO resume-condition counts " \
      "(${heliano_glob_bug_count}/${heliano_glob_fixed_count}): $source_cli"
  fi

  if [[ "$needs_brace_repair" -eq 0 &&
        "$needs_library_repair" -eq 0 &&
        "$needs_divergence_routing" -eq 0 &&
        "$needs_divcalc_tmp_repair" -eq 0 &&
        "$needs_heliano_resume_repair" -eq 0 &&
        "$needs_script_dir_routing" -eq 0 ]]; then
    EARLGREY_ANNOTATION_CLI="$source_cli"
    {
      echo "mode=original-valid"
      echo "installed_version=$EARLGREY_DETECTED_VERSION"
      echo "share_dir=$EARLGREY_SHARE_DIR"
      echo "script_dir=$EARLGREY_SCRIPT_DIR"
      echo "source_cli=$source_cli"
      echo "resolved_cli=$EARLGREY_ANNOTATION_CLI"
      echo "repair_missing_brace=0"
      echo "repair_divergence_library_selection=0"
      echo "repair_divergence_runtime_routing=0"
      echo "repair_process_unique_divergence_tmp=0"
      echo "repair_heliano_resume_glob=0"
      echo "repair_dynamic_script_dir=0"
      sha256sum "$source_cli"
    } > "$patch_audit"
    return 0
  fi

  is_true "$EARLGREY_RUNTIME_PATCH_ANNOTATION_CLI" || {
    echo "[ERROR] Earl Grey annotation-only CLI requires runtime repair: $source_cli" >&2
    echo "[ERROR] Runtime repair is disabled by EARLGREY_RUNTIME_PATCH_ANNOTATION_CLI." >&2
    echo "[ERROR] Check: $original_syntax_log" >&2
    exit 1
  }

  # No exact Earl Grey version string is used as a patch gate. Compatibility is
  # established through the active environment path and the exact source
  # signatures/counts checked above. Future releases that remove these defects
  # simply bypass the relevant repairs.
  if [[ "$needs_brace_repair" -eq 1 ]]; then
    [[ "$(grep -Ec '^[[:space:]]*Checks\(\)[[:space:]]*$' "$source_cli")" -eq 1 ]] || {
      fatal "Refusing to patch Earl Grey CLI without exactly one Checks() function: $source_cli"
    }
    [[ "$(grep -Ec '^[[:space:]]*dfamCheck\(\)[[:space:]]*$' "$source_cli")" -eq 1 ]] || {
      fatal "Refusing to patch Earl Grey CLI without exactly one dfamCheck() function: $source_cli"
    }
    marker_count="$(
      grep -Ec '^[[:space:]]*# Subprocess Dfam Checks[[:space:]]*$' "$source_cli"
    )"
    [[ "$marker_count" -eq 1 ]] || {
      fatal \
        "Refusing to patch Earl Grey CLI with unexpected Dfam marker count " \
        "($marker_count): $source_cli"
    }
  fi

  # Keep awk function calls on a single source line. Some awk builds used by
  # the host reject a newline immediately after sub()/gsub() opening parenthesis.
  if ! awk \
    -v add_brace="$needs_brace_repair" \
    -v fix_library="$needs_library_repair" \
    -v route_divergence="$needs_divergence_routing" \
    -v fix_divcalc_tmp="$needs_divcalc_tmp_repair" \
    -v fix_heliano_resume="$needs_heliano_resume_repair" \
    -v route_script_dir="$needs_script_dir_routing" '
    BEGIN {
      inserted = 0
      fixed = 0
      routed = 0
      tmp_fixed = 0
      heliano_fixed = 0
      script_dir_routed = 0
    }
    add_brace == 1 &&
      /^[[:space:]]*# Subprocess Dfam Checks[[:space:]]*$/ &&
      inserted == 0 {
      print "}"
      inserted = 1
    }
    route_script_dir == 1 &&
      /^[[:space:]]*SCRIPT_DIR=/ &&
      script_dir_routed == 0 {
      print
      match($0, /^[[:space:]]*/)
      print substr($0, RSTART, RLENGTH) \
        "SCRIPT_DIR=\"${EARLGREY_SCRIPT_DIR:-${SCRIPT_DIR}}\""
      script_dir_routed = 1
      next
    }
    fix_library == 1 &&
      /^[[:space:]]*if \[ -z "\$RepSpec" \] && \[ -z "\$startCust" \]; then[[:space:]]*$/ {
      match($0, /^[[:space:]]*/)
      print substr($0, RSTART, RLENGTH) "if [ -z \"$RepSpec\" ]; then"
      fixed = 1
      next
    }
    route_divergence == 1 &&
      index($0, "python ${SCRIPT_DIR}/divergenceCalc/divergence_calc.py") {
      sub(/python \$\{SCRIPT_DIR\}\/divergenceCalc\/divergence_calc\.py/, "python \"${EARLGREY_DIVERGENCE_CALC:-${SCRIPT_DIR}/divergenceCalc/divergence_calc.py}\"")
      routed++
    }
    fix_divcalc_tmp == 1 &&
      /^[[:space:]]*divcalc_tmp="\/tmp\/egdiv_\$\{species\}"[[:space:]]*$/ {
      match($0, /^[[:space:]]*/)
      print substr($0, RSTART, RLENGTH) \
        "divcalc_tmp=\"${TMPDIR:-/tmp}/egdiv_${species}_$$\""
      tmp_fixed++
      next
    }
    fix_heliano_resume == 1 &&
      /^[[:space:]]*if \[ ! -s \$\{OUTDIR\}\/\$\{species\}_heliano\/\*\/RC\.representative\.gff \]; then[[:space:]]*$/ {
      match($0, /^[[:space:]]*/)
      print substr($0, RSTART, RLENGTH) \
        "if ! compgen -G \"${OUTDIR}/${species}_heliano/*/RC.representative.gff\" > /dev/null; then"
      heliano_fixed++
      next
    }
    { print }
    END {
      if (add_brace == 1 && inserted != 1) exit 42
      if (fix_library == 1 && fixed != 1) exit 43
      if (route_divergence == 1 && routed != 2) exit 44
      if (fix_divcalc_tmp == 1 && tmp_fixed != 1) exit 45
      if (fix_heliano_resume == 1 && heliano_fixed != 1) exit 46
      if (route_script_dir == 1 && script_dir_routed != 1) exit 47
    }
  ' "$source_cli" > "$patched_tmp"; then
    rm -f "$patched_tmp"
    fatal "Failed to create runtime-patched Earl Grey annotation-only CLI"
  fi

  if ! bash -n "$patched_tmp" > "$patched_syntax_log" 2>&1; then
    rm -f "$patched_tmp"
    fatal \
      "Runtime-patched Earl Grey annotation-only CLI still fails bash syntax validation: " \
      "$patched_syntax_log"
  fi

  chmod 0755 "$patched_tmp"
  mv -f "$patched_tmp" "$patched_cli"
  EARLGREY_ANNOTATION_CLI="$patched_cli"
  patch_mode="runtime-compatible-repair"

  {
    echo "mode=$patch_mode"
    echo "installed_version=$EARLGREY_DETECTED_VERSION"
    echo "share_dir=$EARLGREY_SHARE_DIR"
    echo "script_dir=$EARLGREY_SCRIPT_DIR"
    echo "source_cli=$source_cli"
    echo "resolved_cli=$EARLGREY_ANNOTATION_CLI"
    echo "repair_missing_brace=$needs_brace_repair"
    echo "repair_divergence_library_selection=$needs_library_repair"
    echo "repair_divergence_runtime_routing=$needs_divergence_routing"
    echo "repair_process_unique_divergence_tmp=$needs_divcalc_tmp_repair"
    echo "repair_heliano_resume_glob=$needs_heliano_resume_repair"
    echo "repair_dynamic_script_dir=$needs_script_dir_routing"
    echo "original_bash_n_log=$original_syntax_log"
    echo "patched_bash_n_log=$patched_syntax_log"
    sha256sum "$source_cli" "$EARLGREY_ANNOTATION_CLI"
  } > "$patch_audit"

  echo \
    "[WARN] Earl Grey ${EARLGREY_DETECTED_VERSION} annotation-only CLI required runtime repair." >&2
  echo \
    "[WARN] Repairs: missing_brace=${needs_brace_repair}, " \
    "divergence_library_selection=${needs_library_repair}, " \
    "divergence_runtime_routing=${needs_divergence_routing}, " \
    "process_unique_divergence_tmp=${needs_divcalc_tmp_repair}, " \
    "heliano_resume_glob=${needs_heliano_resume_repair}, " \
    "dynamic_script_dir=${needs_script_dir_routing}." >&2
  echo "[WARN] Using validated runtime copy: $EARLGREY_ANNOTATION_CLI" >&2
}

earlgrey_write_environment_provenance() {
  local packages="${DIR_LOGS}/earlgrey.environment.packages.txt"
  local explicit="${DIR_LOGS}/earlgrey.environment.explicit.txt"
  local yaml="${DIR_LOGS}/earlgrey.environment.yml"

  if [[ -s "$packages" && -s "$explicit" && -s "$yaml" ]]; then
    return 0
  fi

  "$ENV_RUNNER" list -n "$ENV_EARLGREY" > "${packages}.tmp.$$"
  "$ENV_RUNNER" list -n "$ENV_EARLGREY" --explicit > "${explicit}.tmp.$$"
  "$ENV_RUNNER" env export -n "$ENV_EARLGREY" --no-builds > "${yaml}.tmp.$$"
  mv -f "${packages}.tmp.$$" "$packages"
  mv -f "${explicit}.tmp.$$" "$explicit"
  mv -f "${yaml}.tmp.$$" "$yaml"
}

earlgrey_preflight() {
  local audit="${DIR_LOGS}/04_earlgrey.preflight.txt"

  earlgrey_resolve_installation
  prepare_earlgrey_famdb_runtime
  earlgrey_prepare_divergence_calc
  earlgrey_prepare_annotation_cli
  earlgrey_write_environment_provenance

  # shellcheck disable=SC2016
  env_run_clean_r "$ENV_EARLGREY" env \
    EARLGREY_HELIANO_REQUIRED="$EARLGREY_HELIANO" \
    EARLGREY_ANNOTATION_CLI="$EARLGREY_ANNOTATION_CLI" \
    EARLGREY_SCRIPT_DIR="$EARLGREY_SCRIPT_DIR" \
    EARLGREY_DIVERGENCE_CALC="$EARLGREY_DIVERGENCE_CALC" \
    REQUIRED_DFAM_VERSION="$EARLGREY_REQUIRED_DFAM" \
    FAMDB_DIR="$EARLGREY_FAMDB_RUNTIME_DIR" \
    FAMDB_DATA_DIR="$DFAM_FAMDB_DIR" \
    FAMDB_REAL_PY="$EARLGREY_FAMDB_REAL_PY" \
    EARLGREY_ENV_PREFIX="$EARLGREY_ENV_PREFIX" \
    EARLGREY_SHARE_DIR="$EARLGREY_SHARE_DIR" \
    EARLGREY_DETECTED_VERSION="$EARLGREY_DETECTED_VERSION" \
    bash -lc '
    set -Eeuo pipefail

    echo "Executables"
    echo "==========="
    printf "%-28s%s\n" "earlGreyAnnotationOnly" "${EARLGREY_ANNOTATION_CLI:?}"
    for tool in RepeatMasker RepeatModeler famdb.py bedtools samtools blastn makeblastdb python Rscript; do
      printf "%-28s" "$tool"
      command -v "$tool"
    done

    case "${EARLGREY_HELIANO_REQUIRED:-no}" in
      yes|YES|Yes|true|TRUE|True|1)
        printf "%-28s" "heliano"
        command -v heliano
        ;;
    esac

    echo
    echo "CLI contract"
    echo "============"
    earlgrey_bin="${EARLGREY_ANNOTATION_CLI:?}"
    bash -n "$earlgrey_bin"
    help_text="$("$earlgrey_bin" -h 2>&1)"
    for required_option in "-g" "-s" "-o" "-l" "-t" "-r" "-m" "-d" "-e"; do
      grep -Fq -- "$required_option" <<< "$help_text" || {
        echo "[ERROR] Earl Grey annotation-only CLI is missing option: $required_option" >&2
        exit 1
      }
    done

    echo
    echo "Divergence runtime"
    echo "=================="
    echo "Divergence calculator: ${EARLGREY_DIVERGENCE_CALC:?}"
    [[ -s "$EARLGREY_DIVERGENCE_CALC" ]] || {
      echo "[ERROR] Earl Grey divergence calculator is missing or empty: $EARLGREY_DIVERGENCE_CALC" >&2
      exit 1
    }
    python "$EARLGREY_DIVERGENCE_CALC" -h >/dev/null

    echo
    echo "FamDB binding"
    echo "============="
    echo "FamDB real executable: $FAMDB_REAL_PY"
    echo "FamDB runtime wrapper: $FAMDB_DIR/famdb.py"
    echo "FamDB data directory: $FAMDB_DATA_DIR"
    [[ -x "$FAMDB_DIR/famdb.py" ]] || {
      echo "[ERROR] Earl Grey FamDB runtime wrapper is missing or not executable: $FAMDB_DIR/famdb.py" >&2
      exit 1
    }
    famdb_output="$("$FAMDB_DIR/famdb.py" info 2>&1)"
    printf "%s\n" "$famdb_output"
    grep -Eq "^Version[[:space:]]*:[[:space:]]*${REQUIRED_DFAM_VERSION//./\\.}([[:space:]]|$)" <<< "$famdb_output" || {
      echo "[ERROR] Required Dfam version ${REQUIRED_DFAM_VERSION} was not confirmed in the Earl Grey environment." >&2
      exit 1
    }

    echo
    echo "Earl Grey installation"
    echo "======================"
    echo "Environment prefix: ${EARLGREY_ENV_PREFIX:-unknown}"
    echo "Package directory: ${EARLGREY_SHARE_DIR:-unknown}"
    echo "Script directory: ${EARLGREY_SCRIPT_DIR:-unknown}"
    echo "Detected version: ${EARLGREY_DETECTED_VERSION:-unknown}"

    echo
    echo "Versions"
    echo "========"
    sed -n '1,40p' <<< "$help_text"
    RepeatMasker -v
    RepeatModeler -version
  ' > "$audit" 2>&1

  if ! version_ge "$EARLGREY_DETECTED_VERSION" "$EARLGREY_REQUIRED_VERSION"; then
    echo "[ERROR] Earl Grey ${EARLGREY_DETECTED_VERSION} is older than required minimum ${EARLGREY_REQUIRED_VERSION}." >&2
    echo "[ERROR] Check: $audit" >&2
    exit 1
  fi
  grep -Eq "RepeatMasker version[[:space:]]+${EARLGREY_REQUIRED_REPEATMASKER//./\\.}|RepeatMasker.*${EARLGREY_REQUIRED_REPEATMASKER//./\\.}" "$audit" || {
    echo "[ERROR] RepeatMasker ${EARLGREY_REQUIRED_REPEATMASKER} was not confirmed in ${ENV_EARLGREY}." >&2
    echo "[ERROR] Check: $audit" >&2
    exit 1
  }
  grep -Eq "RepeatModeler version[[:space:]]+${EARLGREY_REQUIRED_REPEATMODELER//./\\.}|RepeatModeler.*${EARLGREY_REQUIRED_REPEATMODELER//./\\.}" "$audit" || {
    echo "[ERROR] RepeatModeler ${EARLGREY_REQUIRED_REPEATMODELER} was not confirmed in ${ENV_EARLGREY}." >&2
    echo "[ERROR] Check: $audit" >&2
    exit 1
  }
}

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

find_has_match() {
  local match
  match="$(find "$@" -print -quit 2>/dev/null)" || return 1
  [[ -n "$match" ]]
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

resolve_repbase_embl_input() {
  # v1.0: RepBase is expected under Libraries/ as RepBase*.embl.tar.gz.
  # Plain *.embl.gz, *.embl, *.ref, and extracted RepBase*.embl directories
  # are also supported for backward compatibility and newer RepBase layouts.
  local detected=""

  if [[ -n "${REPBASE_EMBL:-}" ]]; then
    [[ -f "$REPBASE_EMBL" || -d "$REPBASE_EMBL" ]] || fatal "REPBASE_EMBL was set but is missing: $REPBASE_EMBL"
    [[ -d "$REPBASE_EMBL" || -s "$REPBASE_EMBL" ]] || fatal "REPBASE_EMBL was set but is empty: $REPBASE_EMBL"
    REPBASE_EMBL="$(readlink -f "$REPBASE_EMBL")"
  else
    require_dir "$REPBASE_DIR" "RepBase directory"
    detected="$(
      find "$REPBASE_DIR" -maxdepth 2 \( -type f -o -type d \) \
        \( -name "$REPBASE_EMBL_PATTERN" -o -name "RepBase*.embl.gz" -o -name "RepBase*.embl" \) \
        -printf '%T@\t%p\n' 2>/dev/null \
        | sort -nr \
        | head -n 1 \
        | cut -f2- || true
    )"
    [[ -n "$detected" ]] || fatal "No RepBase EMBL input found in $REPBASE_DIR. Expected ${REPBASE_EMBL_PATTERN}, RepBase*.embl.gz, RepBase*.embl, or an extracted RepBase*.embl directory"
    REPBASE_EMBL="$(readlink -f "$detected")"
  fi

  if [[ ! -d "$REPBASE_EMBL" ]]; then
    case "$REPBASE_EMBL" in
      *.embl|*.embl.gz|*.embl.tar.gz|*.tgz|*.ref|*.ref.gz) ;;
      *)
        fatal "Unsupported RepBase input extension: $REPBASE_EMBL. Use *.embl, *.embl.gz, *.embl.tar.gz, *.tgz, *.ref, *.ref.gz, or an extracted RepBase*.embl directory"
        ;;
    esac
  fi

  echo "[INFO] RepBase input: $REPBASE_EMBL"
}

prepare_dfam_famdb_dir() {
  # Each required Dfam partition may be supplied as an uncompressed *.h5
  # working file or as a *.h5.gz archive that can be validated and decompressed.
  require_dir "$DFAM_FAMDB_DIR" "Dfam FamDB directory"
  mkdir -p "$DIR_LOGS"

  local -a required_partitions=()
  read -r -a required_partitions <<< "$DFAM_REQUIRED_PARTITIONS"
  [[ "${#required_partitions[@]}" -gt 0 ]] || fatal "DFAM_REQUIRED_PARTITIONS is empty"

  local -A seen_partitions=()
  local -a available_archives=()
  local -a available_md5_files=()
  local -a working_partitions=()
  local -a partition_sources=()
  local partition archive md5_file working tmp source

  for partition in "${required_partitions[@]}"; do
    [[ "$partition" == "${DFAM_RELEASE_TAG}."* ]] || {
      fatal "Dfam partition '$partition' does not match DFAM_RELEASE_TAG=$DFAM_RELEASE_TAG"
    }
    [[ -z "${seen_partitions[$partition]:-}" ]] || {
      fatal "Duplicate Dfam partition in DFAM_REQUIRED_PARTITIONS: $partition"
    }
    seen_partitions["$partition"]=1

    archive="$DFAM_FAMDB_DIR/$partition.h5.gz"
    md5_file="$archive.md5"
    working="$DFAM_FAMDB_DIR/$partition.h5"
    source=""

    [[ -s "$archive" ]] && available_archives+=( "$archive" )
    [[ -s "$md5_file" ]] && available_md5_files+=( "$md5_file" )

    if [[ -s "$working" ]]; then
      source="uncompressed"
      if [[ -s "$archive" ]]; then
        source="compressed+uncompressed"
      fi
    else
      [[ -s "$archive" ]] || {
        fatal "Required Dfam partition not found in compressed or uncompressed form: $archive or $working"
      }

      if is_true "$DFAM_VALIDATE_MD5"; then
        [[ -s "$md5_file" ]] || {
          fatal "Checksum required for compressed Dfam archive but not found or empty: $md5_file"
        }
        echo "[INFO] Validating Dfam archive MD5: $md5_file"
        (
          cd "$DFAM_FAMDB_DIR"
          md5sum -c "$(basename "$md5_file")"
        ) || fatal "Dfam archive MD5 validation failed: $md5_file"
      fi

      is_true "$DFAM_AUTODECOMPRESS_GZ" || {
        fatal "Required uncompressed FamDB partition is missing and DFAM_AUTODECOMPRESS_GZ is disabled: $working"
      }
      tmp="${working}.tmp.$$"
      echo "[INFO] Decompressing FamDB partition: $archive -> $working"
      if ! gunzip -c "$archive" > "$tmp"; then
        rm -f "$tmp"
        fatal "Failed to decompress FamDB partition: $archive"
      fi
      [[ -s "$tmp" ]] || {
        rm -f "$tmp"
        fatal "Dfam archive produced an empty working partition: $archive"
      }
      mv -f "$tmp" "$working"
      source="compressed->uncompressed"
    fi

    [[ -s "$working" ]] || {
      fatal "Required FamDB working partition not found or empty: $working"
    }
    working_partitions+=( "$working" )
    partition_sources+=( "$partition=$source" )
  done

  {
    echo "Dfam FamDB audit"
    echo "================"
    echo "DFAM_FAMDB_DIR=$DFAM_FAMDB_DIR"
    echo "DFAM_RELEASE_TAG=$DFAM_RELEASE_TAG"
    echo "DFAM_TAXON=$DFAM_TAXON"
    echo "DFAM_FLAGS=$DFAM_FLAGS"
    echo "DFAM_AUTODECOMPRESS_GZ=$DFAM_AUTODECOMPRESS_GZ"
    echo "DFAM_VALIDATE_MD5=$DFAM_VALIDATE_MD5"
    echo "DFAM_REQUIRED_PARTITIONS=$DFAM_REQUIRED_PARTITIONS"
    echo "Partition source forms:"
    printf '  %s\n' "${partition_sources[@]}"
    echo "Detected compressed archives:"
    if [[ "${#available_archives[@]}" -gt 0 ]]; then
      printf '  %s\n' "${available_archives[@]}"
    else
      echo "  none"
    fi
    echo "Detected archive checksum files:"
    if [[ "${#available_md5_files[@]}" -gt 0 ]]; then
      printf '  %s\n' "${available_md5_files[@]}"
    else
      echo "  none"
    fi
    echo "FamDB HDF5 working partitions:"
    printf '  %s\n' "${working_partitions[@]}"
  } > "${DIR_LOGS}/00A_dfam_famdb_audit.txt"
}

resolve_tetrimmer_script() {
  # Resolve TEtrimmer lazily at the point of use.  This deliberately does not
  # rely on the global TETRIMMER_SCRIPT assignment because repeat-only resumes
  # and strict `set -u` execution may reach this stage with that variable unset.
  local env_prefix="${MINIFORGE_PREFIX}/envs/${ENV_TETRIMMER}"
  local candidate="${TETRIMMER_SCRIPT:-}"

  if [[ -z "$candidate" ]]; then
    candidate="${env_prefix}/share/TEtrimmer/tetrimmer/TEtrimmer.py"
  fi

  if [[ ! -s "$candidate" ]]; then
    candidate="$(
      find "${env_prefix}/share" \
        -maxdepth 5 \
        -type f \
        -path '*/TEtrimmer*/tetrimmer/TEtrimmer.py' \
        -print \
        -quit 2>/dev/null || true
    )"
  fi

  [[ -n "$candidate" && -s "$candidate" ]] || fatal \
    "Could not resolve TEtrimmer.py in environment ${ENV_TETRIMMER}. Expected under ${env_prefix}/share/TEtrimmer/tetrimmer/TEtrimmer.py"

  TETRIMMER_SCRIPT="$(readlink -f "$candidate")"
  export TETRIMMER_SCRIPT
}

prepare_earlgrey_famdb_runtime() {
  local real_famdb_py=""
  local wrapper_tmp="${EARLGREY_FAMDB_WRAPPER}.tmp.$$"

  # shellcheck disable=SC2016
  real_famdb_py="$(env_run_clean_r "$ENV_EARLGREY" bash -lc '
    set -Eeuo pipefail
    famdb_cmd="$(command -v famdb.py 2>/dev/null || true)"
    if [[ -n "$famdb_cmd" && -f "$famdb_cmd" ]]; then
      readlink -f "$famdb_cmd"
      exit 0
    fi

    # RepeatMasker commonly installs FamDB here.
    famdb_cmd="${CONDA_PREFIX:-}/share/RepeatMasker/famdb.py"
    if [[ -f "$famdb_cmd" ]]; then
      readlink -f "$famdb_cmd"
      exit 0
    fi

    # Standalone FamDB conda layouts.
    shopt -s nullglob
    for famdb_cmd in "${CONDA_PREFIX:-}"/share/famdb-*/famdb.py; do
      if [[ -f "$famdb_cmd" ]]; then
        readlink -f "$famdb_cmd"
        exit 0
      fi
    done

    # Final compatibility fallback for future/alternate conda layouts.
    famdb_cmd="$(
      find "${CONDA_PREFIX:-}/share"         -maxdepth 4         -type f         -name famdb.py         -print         -quit 2>/dev/null || true
    )"
    if [[ -n "$famdb_cmd" && -f "$famdb_cmd" ]]; then
      readlink -f "$famdb_cmd"
      exit 0
    fi

    echo "famdb.py not found in ${CONDA_PREFIX:-unset}; checked PATH, share/RepeatMasker/famdb.py, share/famdb-*/famdb.py, and share/**/famdb.py" >&2
    exit 1
  ')" || fatal "Could not resolve famdb.py in conda environment: $ENV_EARLGREY"

  [[ -f "$real_famdb_py" ]] || fatal "Resolved Earl Grey FamDB executable does not exist: $real_famdb_py"
  require_dir "$DFAM_FAMDB_DIR" "Dfam FamDB directory for Earl Grey"
  find_has_match "$DFAM_FAMDB_DIR" -maxdepth 1 -type f -name '*.h5' -size +0c \
    || fatal "No uncompressed Dfam *.h5 partitions found for Earl Grey: $DFAM_FAMDB_DIR"

  mkdir -p "$EARLGREY_FAMDB_RUNTIME_DIR"
  {
    echo '#!/usr/bin/env bash'
    echo 'set -Eeuo pipefail'
    printf 'exec %q -i %q "$@"\n' "$real_famdb_py" "$DFAM_FAMDB_DIR"
  } > "$wrapper_tmp"
  chmod 0755 "$wrapper_tmp"
  mv -f "$wrapper_tmp" "$EARLGREY_FAMDB_WRAPPER"

  EARLGREY_FAMDB_REAL_PY="$real_famdb_py"
}

convert_repbase_to_fasta() {
  local input="$1"
  local output="$2"
  local tmp="${output}.tmp.$$"
  local member=""
  local n_members=0
  local converted_members=0
  local member_tmp=""

  mkdir -p "$(dirname "$output")" "$DIR_LOGS"
  rm -f "$output" "$tmp"

  if [[ -d "$input" ]]; then
    echo "[INFO] Extracting RepBase EMBL records from directory: $input"
    find "$input" -type f \( -name "*.embl" -o -name "*.embl.gz" -o -name "*.ref" -o -name "*.ref.gz" \) \
      | sort > "${DIR_LOGS}/00B_repbase_directory_members.txt"
    n_members="$(wc -l < "${DIR_LOGS}/00B_repbase_directory_members.txt" | tr -d ' ')"
    [[ "$n_members" -gt 0 ]] || fatal "No .embl, .embl.gz, .ref, or .ref.gz file found inside RepBase directory: $input"

    while IFS= read -r member; do
      [[ -n "$member" ]] || continue
      echo "[INFO] Converting RepBase directory member: $member" >&2
      member_tmp="${tmp}.${converted_members}.member"
      rm -f "$member_tmp"
      if [[ "$member" == *.gz ]]; then
        if gunzip -c "$member" | embl_to_fa_tolerant_host > "$member_tmp"; then
          cat "$member_tmp" >> "$tmp"
          converted_members=$(( converted_members + 1 ))
        else
          echo "[WARN] Skipping directory member that did not parse as gzipped RepBase EMBL: $member" >&2
        fi
      else
        if embl_to_fa_tolerant_host < "$member" > "$member_tmp"; then
          cat "$member_tmp" >> "$tmp"
          converted_members=$(( converted_members + 1 ))
        else
          echo "[WARN] Skipping directory member that did not parse as RepBase EMBL: $member" >&2
        fi
      fi
      rm -f "$member_tmp"
    done < "${DIR_LOGS}/00B_repbase_directory_members.txt" \
      2> "${DIR_LOGS}/00B_repbase_convert.log" \
      || error "RepBase FASTA conversion failed"

    [[ "$converted_members" -gt 0 ]] || fatal "No RepBase EMBL records could be parsed from directory: $input"
    require_file "$tmp" "temporary RepBase FASTA conversion"
    mv -f "$tmp" "$output" || error "Failed to promote RepBase FASTA conversion to $output"
    require_file "$output" "RepBase FASTA conversion"
    return 0
  fi

  case "$input" in
    *.embl.tar.gz|*.tgz)
      echo "[INFO] Extracting RepBase EMBL records from tar archive: $input"
      tar -tzf "$input" > "${DIR_LOGS}/00B_repbase_archive_all_members.txt" \
        || fatal "Cannot list RepBase tar archive: $input"
      grep -E '\.(embl|ref)(\.gz)?$' "${DIR_LOGS}/00B_repbase_archive_all_members.txt" > "${DIR_LOGS}/00B_repbase_archive_members.txt" || true
      n_members="$(wc -l < "${DIR_LOGS}/00B_repbase_archive_members.txt" | tr -d ' ')"
      if [[ "$n_members" -lt 1 ]]; then
        echo "[WARN] No member ending in .embl, .embl.gz, .ref, or .ref.gz found; scanning non-directory archive members by content"
        sed '/\/$/d' "${DIR_LOGS}/00B_repbase_archive_all_members.txt" > "${DIR_LOGS}/00B_repbase_archive_members.txt"
        n_members="$(wc -l < "${DIR_LOGS}/00B_repbase_archive_members.txt" | tr -d ' ')"
        [[ "$n_members" -gt 0 ]] || fatal "No non-directory member found inside RepBase archive: $input"
      fi

      while IFS= read -r member; do
        [[ -n "$member" ]] || continue
        echo "[INFO] Converting RepBase archive member: $member" >&2
        member_tmp="${tmp}.${converted_members}.member"
        rm -f "$member_tmp"
        if [[ "$member" == *.gz ]]; then
          if tar -xOzf "$input" "$member" | gunzip -c | embl_to_fa_tolerant_host > "$member_tmp"; then
            cat "$member_tmp" >> "$tmp"
            converted_members=$(( converted_members + 1 ))
          else
            echo "[WARN] Skipping archive member that did not parse as gzipped RepBase EMBL: $member" >&2
          fi
        else
          if tar -xOzf "$input" "$member" | embl_to_fa_tolerant_host > "$member_tmp"; then
            cat "$member_tmp" >> "$tmp"
            converted_members=$(( converted_members + 1 ))
          else
            echo "[WARN] Skipping archive member that did not parse as RepBase EMBL: $member" >&2
          fi
        fi
        rm -f "$member_tmp"
      done < "${DIR_LOGS}/00B_repbase_archive_members.txt" \
        2> "${DIR_LOGS}/00B_repbase_convert.log" \
        || error "RepBase FASTA conversion failed"

      [[ "$converted_members" -gt 0 ]] || fatal "No RepBase EMBL records could be parsed from archive: $input"
      ;;
    *.embl.gz)
      gunzip -c "$input" | embl_to_fa_tolerant_host > "$tmp" \
        2> "${DIR_LOGS}/00B_repbase_convert.log" \
        || error "RepBase FASTA conversion failed"
      ;;
    *.embl|*.ref)
      embl_to_fa_tolerant_host < "$input" > "$tmp" \
        2> "${DIR_LOGS}/00B_repbase_convert.log" \
        || error "RepBase FASTA conversion failed"
      ;;
    *.ref.gz)
      gunzip -c "$input" | embl_to_fa_tolerant_host > "$tmp" \
        2> "${DIR_LOGS}/00B_repbase_convert.log" \
        || error "RepBase FASTA conversion failed"
      ;;
    *)
      fatal "Unsupported RepBase input extension: $input"
      ;;
  esac

  require_file "$tmp" "temporary RepBase FASTA conversion"
  mv -f "$tmp" "$output" || error "Failed to promote RepBase FASTA conversion to $output"
  require_file "$output" "RepBase FASTA conversion"
}

file_has_fasta_records() {
  local path="$1"

  [[ -s "$path" ]] || return 1
  grep -q "^>" "$path"
}

dir_has_nonempty_file_or_link() {
  local dir="$1"
  local match

  [[ -d "$dir" ]] || return 1
  # Let find stop itself after the first nonempty file. Returning early from a
  # process-substitution reader closes its pipe and makes find fail with SIGPIPE
  # (rc=141) under the global ERR trap, even though the directory is valid.
  match="$(find "$dir" \( -type f -o -type l \) -exec test -s {} \; -print -quit 2>/dev/null)" || return 1
  [[ -n "$match" ]]
}

auto_detect_outputs() {
  if ! is_true "$AUTO_RESUME"; then
    return 0
  fi

  echo "[INFO] Auto-resume enabled: scanning for existing outputs"

  if [[ -n "${KNOWN_LIBRARY_OVERRIDE:-}" ]]; then
    KNOWN_NR="$KNOWN_LIBRARY_OVERRIDE"
  elif [[ ! -s "${KNOWN_NR:-}" ]]; then
    detected_known="$(find_newest_file_by_name "$DIR_KNOWN" \( -name "known_repeats.Dfam4_RepBase.nr.fa" -o -name "known_repeats.nr.fa" -o -name "*known*nr*.fa" \))"
    [[ -n "$detected_known" ]] && KNOWN_NR="$detected_known"
  fi

  if [[ -n "${REPEATMODELER_LIBRARY_OVERRIDE:-}" ]]; then
    RM_LIBRARY="$REPEATMODELER_LIBRARY_OVERRIDE"
  fi

  if [[ -n "${TERRIER_LIBRARY_OVERRIDE:-}" ]]; then
    TERRIER_LIBRARY="$TERRIER_LIBRARY_OVERRIDE"
  elif [[ -s "${DIR_TERRIER:-}/repeatmodeler.terrier.fa" ]]; then
    TERRIER_LIBRARY="${DIR_TERRIER}/repeatmodeler.terrier.fa"
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
    detected_final="$(find_newest_file_by_name "$DIR_LIB" \( -name "${SPECIES}.TEtrimmer_Dfam4_RepBase.nr.fa" -o -name "${SPECIES}.TEtrimmer_Dfam_RepBase.nr.fa" -o -name "*TEtrimmer*Dfam*RepBase*nr*.fa" -o -name "*final*nr*.fa" \))"
    [[ -n "$detected_final" ]] && FINAL_NR="$detected_final"
  fi

  if [[ -s "${KNOWN_NR:-}" ]]; then
    echo "[INFO] Auto-detected known library: $KNOWN_NR"
  fi

  if [[ -s "${RM_LIBRARY:-}" ]]; then
    echo "[INFO] Auto-detected RepeatModeler library: $RM_LIBRARY"
  fi

  if [[ -s "${TERRIER_LIBRARY:-}" ]]; then
    echo "[INFO] Auto-detected Terrier library: $TERRIER_LIBRARY"
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
    terrier)
      if ! is_true "${TERRIER_ENABLE:-TRUE}"; then return 0; fi
      [[ -s "${TERRIER_LIBRARY:-}" ]] && file_has_fasta_records "$TERRIER_LIBRARY" \
        && [[ -s "${DIR_TERRIER}/terrier_predictions.csv" ]] \
        && [[ -s "${DIR_TERRIER}/terrier_tip_probabilities.csv" ]]
      ;;
    tetrimmer)
      tetrimmer_core_outputs_complete
      ;;
    final_library)
      [[ -s "${FINAL_NR:-}" ]] && file_has_fasta_records "$FINAL_NR"
      ;;
    earlgrey)
      earlgrey_core_outputs_complete
      ;;
    earlgrey_plots)
      if ! is_true "${EARLGREY_ENHANCED_PLOTS:-TRUE}"; then
        return 0
      fi
      [[ -d "${DIR_EARLGREY_ENHANCED_PLOTS:-}" ]] || return 1
      local required
      for required in \
        enhanced_plots_manifest.tsv \
        enhanced_plots_inputs.tsv \
        enhanced_plots_schema.tsv \
        palette_audit_subclass.tsv \
        palette_audit_superfamily.tsv
      do
        [[ -s "${DIR_EARLGREY_ENHANCED_PLOTS}/${required}" ]] || return 1
      done
      local ext
      for ext in png pdf svg; do
        find_has_match "$DIR_EARLGREY_ENHANCED_PLOTS" -maxdepth 1 -type f -name "*.${ext}" -size +0c || return 1
      done
      while IFS= read -r plot_file; do
        [[ -s "$plot_file" ]] || return 1
      done < <(awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) h[$i] = i; next} {print $h["file"]}' "${DIR_EARLGREY_ENHANCED_PLOTS}/enhanced_plots_manifest.tsv")
      ;;
    collect_outputs)
      [[ -d "$DIR_FINAL/logs" ]] \
        && [[ -s "${DIR_FINAL}/logs/collect_outputs_manifest.tsv" ]] \
        && [[ "$(head -n 1 "${DIR_FINAL}/logs/collect_outputs_manifest.tsv")" == $'category\tsource\tdestination\tsize_bytes\tlink_or_copy' ]] \
        && [[ -L "${DIR_FINAL}/earlgrey/repeatmasker_internal" ]] \
        && [[ -L "${DIR_FINAL}/earlgrey/summaryFiles" ]] \
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
  local base
  [[ -d "$DIR_RM" ]] || return 0

  newest_rm_dir="$(find "$DIR_RM" -maxdepth 1 -type d -name 'RM_*' -printf '%T@\t%p\n' 2>/dev/null | sort -nr | head -n1 | cut -f2- || true)"
  [[ -n "$newest_rm_dir" && -d "$newest_rm_dir" ]] || return 0

  if [[ -s "$newest_rm_dir/consensi.fa.classified" ]]; then
    echo "[INFO] Complete classified RepeatModeler directory: $newest_rm_dir"
    return 0
  fi

  if [[ -s "$newest_rm_dir/consensi.fa" && -s "$newest_rm_dir/families.stk" ]]; then
    echo "[INFO] Preserving discovery-complete directory for RepeatClassifier: $newest_rm_dir"
    return 0
  fi

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

tetrimmer_famdb_errors_present() {
  local log_file="${DIR_LOGS}/02_tetrimmer.log"

  [[ -s "$log_file" ]] || return 1
  LC_ALL=C grep -Fq \
    -e 'Could not determine FamDB version' \
    -e 'FamDB data directory not found' \
    "$log_file"
}

resolve_tetrimmer_result_dir() {
  # TEtrimmer 1.7.4 creates TEtrimmer_output_<timestamp>/ whenever --output_dir
  # is not empty. Prefer a complete direct layout, otherwise the newest complete
  # nested layout, then the newest resumable layout.
  TETRIMMER_RESULT_DIR=""

  if [[ -d "$DIR_TETRIMMER" ]]; then
    if [[ -s "${DIR_TETRIMMER}/TEtrimmer_consensus_merged.fasta" || -s "${DIR_TETRIMMER}/TEtrimmer_consensus.fasta" ]]; then
      TETRIMMER_RESULT_DIR="$DIR_TETRIMMER"
    else
      local candidate=""
      while IFS=$'\t' read -r _mtime candidate; do
        [[ -d "$candidate" ]] || continue
        if [[ -s "${candidate}/TEtrimmer_consensus_merged.fasta" || -s "${candidate}/TEtrimmer_consensus.fasta" ]]; then
          TETRIMMER_RESULT_DIR="$candidate"
          break
        fi
      done < <(
        find "$DIR_TETRIMMER" -mindepth 1 -maxdepth 1 -type d \
          -name 'TEtrimmer_output_*' -printf '%T@\t%p\n' 2>/dev/null \
          | sort -nr
      )

      if [[ -z "$TETRIMMER_RESULT_DIR" ]]; then
        while IFS=$'\t' read -r _mtime candidate; do
          [[ -d "$candidate" ]] || continue
          if [[ -s "${candidate}/summary.txt" &&
                -s "${candidate}/Sequence_name_mapping.txt" &&
                -d "${candidate}/Multiple_sequence_alignment" ]]; then
            TETRIMMER_RESULT_DIR="$candidate"
            break
          fi
        done < <(
          find "$DIR_TETRIMMER" -mindepth 1 -maxdepth 1 -type d \
            -name 'TEtrimmer_output_*' -printf '%T@\t%p\n' 2>/dev/null \
            | sort -nr
        )
      fi
    fi
  fi

  export TETRIMMER_RESULT_DIR
  [[ -n "$TETRIMMER_RESULT_DIR" ]]
}

tetrimmer_core_outputs_complete() {
  local result_dir=""
  resolve_tetrimmer_result_dir || return 1
  result_dir="$TETRIMMER_RESULT_DIR"

  [[ -s "${result_dir}/TEtrimmer_consensus_merged.fasta" || -s "${result_dir}/TEtrimmer_consensus.fasta" ]] || return 1
  [[ -s "${result_dir}/summary.txt" ]] || return 1
  [[ -s "${result_dir}/Sequence_name_mapping.txt" ]] || return 1
  tetrimmer_famdb_errors_present && return 1
  return 0
}

tetrimmer_dir_resumable() {
  local result_dir=""
  local singles=""
  local summary=""

  if resolve_tetrimmer_result_dir; then
    result_dir="$TETRIMMER_RESULT_DIR"
  else
    result_dir="$DIR_TETRIMMER"
  fi

  singles="${result_dir}/Single_fasta_files"
  summary="${result_dir}/summary.txt"

  [[ -d "$singles" ]] || return 1
  [[ -s "$summary" ]] || return 1
  head -n 1 "$summary" | grep -q '^input_name,output_name,.*status$' || return 1
  if ! find_has_match "$singles" -maxdepth 1 -type f -name "*.fasta" -size +0c; then
    awk -F',' 'NR > 1 && $NF == "processed" {found = 1; exit} END {exit found ? 0 : 1}' "$summary" || return 1
  fi
  [[ -s "${result_dir}/Sequence_name_mapping.txt" ]] || return 1
  [[ -d "${result_dir}/Multiple_sequence_alignment" ]] || return 1
  return 0
}

prepare_tetrimmer_output_dir() {
  [[ -d "$DIR_TETRIMMER" ]] || {
    mkdir -p "$DIR_TETRIMMER"
    return 0
  }

  if tetrimmer_core_outputs_complete; then
    return 0
  fi

  if [[ "$TETRIMMER_CONTINUE_MODE" != "never" ]] && is_true "$TETRIMMER_CONTINUE" && tetrimmer_dir_resumable; then
    echo "[INFO] Compatible interrupted TEtrimmer v${TETRIMMER_REQUIRED_VERSION} output detected: $TETRIMMER_RESULT_DIR"
    echo "[INFO] Preserving it for continuation."
    return 0
  fi

  if [[ "$TETRIMMER_CONTINUE_MODE" == "always" ]]; then
    echo "[ERROR] TETRIMMER_CONTINUE_MODE=always but no resumable TEtrimmer result directory was found under: $DIR_TETRIMMER" >&2
    exit 1
  fi

  if [[ -z "$(find "$DIR_TETRIMMER" -mindepth 1 -print -quit 2>/dev/null)" ]]; then
    return 0
  fi

  if is_true "$TETRIMMER_ARCHIVE_INCOMPLETE"; then
    mkdir -p "${PROJECT_ROOT}/failed_runs/tetrimmer"
    local stamp dest
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
  local result_dir="$DIR_TETRIMMER"
  resolve_tetrimmer_result_dir && result_dir="$TETRIMMER_RESULT_DIR"

  {
    echo "TEtrimmer output audit"
    echo "======================"
    echo "Requested output dir: $DIR_TETRIMMER"
    echo "Resolved result dir:  $result_dir"
    echo

    for f in \
      "$result_dir/TEtrimmer_consensus.fasta" \
      "$result_dir/TEtrimmer_consensus_merged.fasta" \
      "$result_dir/summary.txt" \
      "$result_dir/Sequence_name_mapping.txt" \
      "$result_dir/Multiple_sequence_alignment/error_file.txt"
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
      Classification_and_deduplication \
      Single_fasta_files \
      TEtrimmer_for_proof_curation \
      TEtrimmer_for_proof_curation/Annotations_perfect \
      TEtrimmer_for_proof_curation/Annotations_good \
      TEtrimmer_for_proof_curation/Annotations_check_recommended \
      TEtrimmer_for_proof_curation/Annotations_check_required \
      TEtrimmer_for_proof_curation/Clustered_proof_curation \
      TEtrimmer_for_proof_curation/TE_low_copy \
      TEtrimmer_for_proof_curation/TE_skipped \
      HMM_files \
      Multiple_sequence_alignment
    do
      if [[ -d "$result_dir/$d" ]]; then
        n=$(find "$result_dir/$d" -type f 2>/dev/null | wc -l)
        echo "[OK] $d files=$n"
      else
        echo "[MISSING] $d"
      fi
    done
  } > "${DIR_LOGS}/02_tetrimmer.output_audit.txt"
}

migrate_legacy_tetrimmer_pfam() {
  local legacy="${DIR_TETRIMMER}/pfam"

  # Only migrate when the user did not explicitly override the new resource path.
  if [[ "$TETRIMMER_PFAM_DIR" == "${PROJECT_ROOT}/resources/pfam" &&
        -d "$legacy" &&
        "$legacy" != "$TETRIMMER_PFAM_DIR" ]]; then
    if [[ ! -e "$TETRIMMER_PFAM_DIR" ]]; then
      mkdir -p "$(dirname "$TETRIMMER_PFAM_DIR")"
      echo "[INFO] Migrating legacy TEtrimmer Pfam directory outside the output tree:"
      echo "       $legacy -> $TETRIMMER_PFAM_DIR"
      mv "$legacy" "$TETRIMMER_PFAM_DIR"
    elif [[ -d "$TETRIMMER_PFAM_DIR" ]]; then
      echo "[INFO] External TEtrimmer Pfam resource already exists: $TETRIMMER_PFAM_DIR"
    fi
  fi
}

earlgrey_famdb_errors_present() {
  local log_file="${DIR_LOGS}/04_earlgrey_annotation_only.log"

  [[ -s "$log_file" ]] || return 1
  LC_ALL=C grep -Fq \
    -e 'Could not determine FamDB version' \
    -e 'FamDB data directory not found' \
    -e 'Could not execute famdb.py using' \
    "$log_file"
}

earlgrey_divergence_filename_errors_present() {
  local log_file="${DIR_LOGS}/04_earlgrey_annotation_only.log"

  [[ -s "$log_file" ]] || return 1
  LC_ALL=C grep -Fq 'FileNotFoundError:' "$log_file" \
    && LC_ALL=C grep -Fq '/split_library/' "$log_file"
}


earlgrey_prepare_library_ids() {
  # RepeatCraft encodes groups as chromosome|family|number and splits on every
  # pipe. Source prefixes such as DFAM| and REPBASE| therefore crash its merge.
  # Escape only the FASTA identifier; keep the class and description intact.
  local source_fa="$1"
  local output_fa="${DIR_LOGS}/${SPECIES}.earlgrey_safe_ids.fa"
  local id_map="${DIR_LOGS}/${SPECIES}.earlgrey_id_map.tsv"
  local tmp_fa="${output_fa}.tmp.$$"
  local tmp_map="${id_map}.tmp.$$"
  local cache_manifest="${DIR_LOGS}/${SPECIES}.earlgrey_safe_ids.source.tsv"
  local source_sha

  source_sha="$(sha256sum -- "$source_fa" | cut -d ' ' -f 1)" ||
    fatal "Could not fingerprint EarlGrey source library: $source_fa"
  if [[ -s "$output_fa" && -s "$id_map" && -s "$cache_manifest" ]] &&
     cmp -s "$cache_manifest" <(printf 'format\tsafe-id-v1\nsource_sha256\t%s\nsafe_bytes\t%s\nmap_bytes\t%s\n' \
       "$source_sha" "$(stat -c '%s' "$output_fa")" "$(stat -c '%s' "$id_map")"); then
    EARLGREY_LIBRARY="$output_fa"
    echo "[INFO] Reusing EarlGrey-safe ID library for unchanged source: $output_fa"
    echo "[INFO] EarlGrey ID mapping: $id_map"
    return 0
  fi

  if ! python3 - "$source_fa" "$tmp_fa" "$tmp_map" <<'PYIDS'
import re
import sys

source, output, mapping = sys.argv[1:]
seen = set()
count = escaped = disambiguated = 0
with open(source, encoding="utf-8") as inp, \
     open(output, "w", encoding="utf-8") as out, \
     open(mapping, "w", encoding="utf-8") as names:
    names.write("record_index\toriginal_id\tearlgrey_id\n")
    for line in inp:
        if not line.startswith(">"):
            out.write(line)
            continue
        header = line[1:].rstrip("\r\n")
        match = re.fullmatch(r"(\S+)(.*)", header)
        if match is None:
            raise ValueError(f"Invalid FASTA header at record {count + 1}")
        token, description = match.groups()
        identifier, separator, repeat_class = token.partition("#")
        count += 1
        safe_base = identifier.replace("|", "_PIPE_")
        safe_id = safe_base
        suffix = 0
        while safe_id in seen:
            suffix += 1
            safe_id = f"{safe_base}_EG{count}_{suffix}"
        seen.add(safe_id)
        escaped += "|" in identifier
        disambiguated += suffix > 0
        names.write(f"{count}\t{identifier}\t{safe_id}\n")
        out.write(">" + safe_id + (separator + repeat_class if separator else "") + description + "\n")
if not count:
    raise ValueError("EarlGrey library contains no FASTA records")
print(f"[INFO] EarlGrey library IDs: {count} records, {escaped} escaped pipes, {disambiguated} duplicate/colliding names disambiguated")
PYIDS
  then
    rm -f "$tmp_fa" "$tmp_map"
    fatal "Could not create EarlGrey-safe identifiers from: $source_fa"
  fi

  mv -f "$tmp_fa" "$output_fa"
  mv -f "$tmp_map" "$id_map"
  printf 'format\tsafe-id-v1\nsource_sha256\t%s\nsafe_bytes\t%s\nmap_bytes\t%s\n' \
    "$source_sha" "$(stat -c '%s' "$output_fa")" "$(stat -c '%s' "$id_map")" > "${cache_manifest}.tmp.$$"
  mv -f "${cache_manifest}.tmp.$$" "$cache_manifest"
  EARLGREY_LIBRARY="$output_fa"
  echo "[INFO] EarlGrey ID mapping: $id_map"
}

earlgrey_repeatcraft_merge_failed() {
  local log_file="${DIR_LOGS}/04_earlgrey_annotation_only.log"
  [[ -s "$log_file" ]] || return 1
  LC_ALL=C grep -Fq 'extraTrueMergeTEm.py' "$log_file" &&
    LC_ALL=C grep -Fq 'ValueError: too many values to unpack (expected 3)' "$log_file"
}

# Empty HELIANO files are valid only when the original log establishes that
# all reported subgenomes had zero candidates and both sorting phases finished.
earlgrey_heliano_confirmed_zero() {
  local hel_dir="${1:?missing HELIANO directory}"
  local log_file="${DIR_LOGS}/04_earlgrey_annotation_only.log"
  local gff bed
  [[ -s "$log_file" ]] || return 1
  gff="$(find "$hel_dir" -mindepth 2 -maxdepth 2 -type f -name RC.representative.gff -empty -print -quit 2>/dev/null)" || return 1
  [[ -n "$gff" ]] || return 1
  bed="${gff%.gff}.bed"
  [[ -f "$bed" && ! -s "$bed" ]] || return 1
  LC_ALL=C awk '
    /<<< Running HELIANO to Detect Helitrons >>>/ { active=1; next }
    active && /<<< (Defragmenting Repeats|Repeats already defragmented, skipping)/ { finished=1; exit }
    active && /^Find [0-9]+ rep-hel blocks in subgenome[0-9]+\./ {
      if ($2 != 0) bad=1; count++
    }
    active && /Sort representative files finished/ { representative=1 }
    active && /Sort alternative files finished/ { alternative=1 }
    active && /(Traceback|(^|[[:space:]])([Ee][Rr][Rr][Oo][Rr]|[Kk]illed|[Ff]atal)([[:space:]:]|$))/ { bad=1 }
    END { exit !(active && finished && count > 0 && representative && alternative && !bad) }
  ' "$log_file"
}

earlgrey_core_outputs_complete() {
  # Earl Grey creates many large intermediate files before the final report is
  # complete. Use the observed completed Earl Grey annotation-only layout as
  # the success contract: summaryFiles + RepeatLandscape + looseMerge finals.
  local eg_root="${DIR_EARLGREY}/${SPECIES}_EarlGrey"
  local merged="${eg_root}/${SPECIES}_mergedRepeats"
  local loose="${merged}/looseMerge"
  local summary="${eg_root}/${SPECIES}_summaryFiles"
  local landscape="${eg_root}/${SPECIES}_RepeatLandscape"
  local rm_dir="${eg_root}/${SPECIES}_RepeatMasker_Against_Custom_Library"

  earlgrey_famdb_errors_present && return 1
  earlgrey_divergence_filename_errors_present && return 1
  [[ -d "$eg_root" ]] || return 1
  [[ -d "$summary" ]] || return 1
  [[ -d "$landscape" ]] || return 1
  [[ -d "$rm_dir" ]] || return 1
  find_has_match "$rm_dir" -maxdepth 1 -type f -name "*.out" -size +0c || return 1
  find_has_match "$rm_dir" -maxdepth 1 -type f -name "*.tbl" -size +0c || return 1
  find_has_match "$rm_dir" -maxdepth 1 -type f -name "*.align" -size +0c || return 1

  # Final annotation products copied by Earl Grey into summaryFiles.
  [[ -s "${summary}/${SPECIES}.filteredRepeats.bed" ]] || return 1
  [[ -s "${summary}/${SPECIES}.filteredRepeats.gff" ]] || return 1

  # Final abundance summaries.
  [[ -s "${summary}/${SPECIES}.highLevelCount.txt" ]] || return 1
  [[ -s "${summary}/${SPECIES}.familyLevelCount.txt" ]] || return 1

  # Divergence table and Kimura-bearing GFF.
  [[ -s "${summary}/${SPECIES}_divergence_summary_table.tsv" ]] || return 1
  [[ -s "${landscape}/${SPECIES}.filteredRepeats.withDivergence.gff" ]] || return 1

  # Optional softmask and HELIANO products are required only when requested.
  if is_true "$EARLGREY_SOFTMASK"; then
    [[ -s "${summary}/${SPECIES}.softmasked.fasta" ]] || return 1
  fi
  if is_true "$EARLGREY_HELIANO"; then
    if ! find_has_match "${eg_root}/${SPECIES}_heliano" -mindepth 2 -maxdepth 2 -type f -name "RC.representative.gff" -size +0c; then
      earlgrey_heliano_confirmed_zero "${eg_root}/${SPECIES}_heliano" || return 1
    fi
  fi

  # Earl Grey native plots.
  [[ -s "${summary}/${SPECIES}_classification_landscape.pdf" || -s "${landscape}/${SPECIES}_classification_landscape.pdf" ]] || return 1
  [[ -s "${summary}/${SPECIES}_split_class_landscape.pdf" || -s "${landscape}/${SPECIES}_split_class_landscape.pdf" ]] || return 1
  [[ -s "${summary}/${SPECIES}_superfamily_div_plot.pdf" || -s "${landscape}/${SPECIES}_superfamily_div_plot.pdf" ]] || return 1
  [[ -s "${summary}/${SPECIES}.summaryPie.pdf" ]] || return 1

  # looseMerge finals are useful to confirm that the defragmentation/filtering
  # stage really completed. These are not enough alone, but they should exist
  # in a complete annotation-only run.
  [[ -s "${loose}/${SPECIES}.rmerge.gff" ]] || return 1
  [[ -s "${loose}/${SPECIES}.rmerge.gff.filtered" ]] || return 1
  [[ -s "${loose}/${SPECIES}.filteredRepeats.bed" ]] || return 1
  [[ -s "${loose}/${SPECIES}.filteredRepeats.gff" ]] || return 1

  return 0
}

audit_earlgrey_outputs() {
  local eg_root="${DIR_EARLGREY}/${SPECIES}_EarlGrey"
  local summary="${eg_root}/${SPECIES}_summaryFiles"
  local landscape="${eg_root}/${SPECIES}_RepeatLandscape"
  local loose="${eg_root}/${SPECIES}_mergedRepeats/looseMerge"
  local rm_dir="${eg_root}/${SPECIES}_RepeatMasker_Against_Custom_Library"
  local audit="${DIR_LOGS}/04_earlgrey_annotation_only.output_audit.tsv"
  local f
  local stage

  mkdir -p "${DIR_LOGS}"

  {
    printf 'stage\tlabel\tpath\texists\tsize_bytes\n'

    while IFS= read -r -d '' f; do
      printf 'repeatmasker_internal\t%s\t%s\tTRUE\t%s\n' "$(basename "$f")" "$f" "$(stat -c '%s' "$f")"
    done < <(find "$rm_dir" -maxdepth 1 -type f \( -name "*.out" -o -name "*.tbl" -o -name "*.align" -o -name "*.masked" \) -size +0c -print0 2>/dev/null)

    for f in \
      "${summary}/${SPECIES}.filteredRepeats.bed" \
      "${summary}/${SPECIES}.filteredRepeats.gff" \
      "${summary}/${SPECIES}.highLevelCount.txt" \
      "${summary}/${SPECIES}.familyLevelCount.txt" \
      "${summary}/${SPECIES}_divergence_summary_table.tsv" \
      "${summary}/${SPECIES}.softmasked.fasta" \
      "${summary}/${SPECIES}_classification_landscape.pdf" \
      "${summary}/${SPECIES}_split_class_landscape.pdf" \
      "${summary}/${SPECIES}_superfamily_div_plot.pdf" \
      "${summary}/${SPECIES}.summaryPie.pdf" \
      "${landscape}/${SPECIES}.filteredRepeats.withDivergence.gff" \
      "${landscape}/${SPECIES}_summary_table.tsv" \
      "${loose}/${SPECIES}.rmerge.gff.filtered" \
      "${loose}/${SPECIES}.filteredRepeats.bed" \
      "${loose}/${SPECIES}.filteredRepeats.gff"
    do
      if [[ -e "$f" ]]; then
        case "$f" in
          "$summary"/*) stage="earlgrey_summary" ;;
          "$landscape"/*) stage="earlgrey_landscape" ;;
          *) stage="earlgrey_merge" ;;
        esac
        printf '%s\t%s\t%s\tTRUE\t%s\n' "$stage" "$(basename "$f")" "$f" "$(stat -c '%s' "$f")"
      else
        printf 'earlgrey_expected\t%s\t%s\tFALSE\t0\n' "$(basename "$f")" "$f"
      fi
    done
  } > "$audit"
}
earlgrey_prepare_resume_from_failed_merge() {
  # Preserve expensive upstream results; archive stale merge products that
  # prevent Earl Grey from resuming safely.
  local eg_root="${DIR_EARLGREY}/${SPECIES}_EarlGrey"
  local merged="${eg_root}/${SPECIES}_mergedRepeats"
  local summary="${eg_root}/${SPECIES}_summaryFiles"
  local landscape="${eg_root}/${SPECIES}_RepeatLandscape"
  local stamp
  local rescue_dir

  [[ -d "$eg_root" ]] || return 0

  stamp="$(date +%Y%m%d_%H%M%S)"
  rescue_dir="${PROJECT_ROOT}/failed_runs/earlgrey_resume_cleanup/${stamp}"
  mkdir -p "$rescue_dir"

  echo "[WARN] Earl Grey output directory is incomplete; attempting resume-safe cleanup."
  echo "[WARN] Keeping completed upstream outputs such as RepeatMasker and HELIANO."
  echo "[WARN] Moving stale downstream files to: $rescue_dir"

  mkdir -p "$merged" "$summary"

  if [[ -e "${merged}/looseMerge" ]]; then
    mv "${merged}/looseMerge" "${rescue_dir}/looseMerge.${stamp}"
  fi

  for f in \
    "${merged}/${SPECIES}.filteredRepeats.bed" \
    "${merged}/${SPECIES}.filteredRepeats.bed.2" \
    "${merged}/${SPECIES}.rmerge.gff.filtered" \
    "${merged}/${SPECIES}.rmerge.gff.filtered.2" \
    "${summary}/${SPECIES}.filteredRepeats.bed" \
    "${summary}/${SPECIES}.filteredRepeats.gff" \
    "${summary}/${SPECIES}.filteredRepeats.withDivergence.gff" \
    "${summary}/${SPECIES}.highLevelCount.txt" \
    "${summary}/${SPECIES}.familyLevelCount.txt" \
    "${landscape}/${SPECIES}.filteredRepeats.withDivergence.gff"
  do
    if [[ -e "$f" ]]; then
      mkdir -p "${rescue_dir}/$(dirname "${f#"${eg_root}/"}")"
      mv "$f" "${rescue_dir}/${f#"${eg_root}/"}.${stamp}"
    fi
  done

  rm -f \
    "$(repeat_checkpoint_file earlgrey)" \
    "$(repeat_checkpoint_file earlgrey_plots)" \
    "$(repeat_checkpoint_file collect_outputs)" \
    "$(repeat_checkpoint_file non_repeats)"

  return 0
}

earlgrey_archive_incompatible_masker_cache() {
  # Earl Grey skips RepeatMasker when any *.tbl exists. Require its matching
  # nonempty .out and .align files before trusting that resume marker. Also
  # reject annotations made with the original pipe-bearing library IDs.
  local rm_dir="${DIR_EARLGREY}/${SPECIES}_EarlGrey/${SPECIES}_RepeatMasker_Against_Custom_Library"
  local rm_out="" rm_tbl="" rm_align="" reason=""
  local archive_root="${PROJECT_ROOT}/failed_runs/earlgrey_incompatible_repeatmasker"
  local archive_dir=""

  [[ -d "$rm_dir" ]] || return 0
  rm_tbl="$(find "$rm_dir" -maxdepth 1 -type f -name '*.tbl' -print -quit)" || return 1
  if [[ -n "$rm_tbl" ]]; then
    rm_out="${rm_tbl%.tbl}.out"
    rm_align="${rm_tbl%.tbl}.align"
  else
    rm_out="$(find "$rm_dir" -maxdepth 1 -type f -name '*.out' -print -quit)" || return 1
  fi

  if [[ -n "$rm_tbl" && ( ! -s "$rm_tbl" || ! -s "$rm_out" || ! -s "$rm_align" ) ]]; then
    reason="RepeatMasker resume table exists without complete .out/.align outputs"
  elif [[ -s "$rm_out" ]] &&
       awk 'NR > 3 && NF >= 11 && index($10, "|") { found = 1; exit }
            END { exit !found }' "$rm_out"; then
    reason="RepeatMasker cache contains source IDs with pipes"
  else
    return 0
  fi

  mkdir -p "$archive_root"
  archive_dir="$(mktemp -d "${archive_root}/${SPECIES}.XXXXXXXX")" ||
    fatal "Could not create archive for incompatible Earl Grey RepeatMasker cache"
  echo "[WARN] $reason; archiving before reannotation."
  echo "[WARN] Archived RepeatMasker directory: $rm_dir -> $archive_dir/$(basename "$rm_dir")"
  mv "$rm_dir" "$archive_dir/" ||
    fatal "Could not archive incompatible Earl Grey RepeatMasker directory: $rm_dir"
  mkdir -p "$rm_dir"
}

prepare_earlgrey_output_dir() {
  [[ -d "$DIR_EARLGREY" ]] || return 0

  if earlgrey_core_outputs_complete; then
    return 0
  fi

  if dir_has_nonempty_file_or_link "$DIR_EARLGREY"; then
    case "$EARLGREY_INCOMPLETE_POLICY" in
      resume)
        earlgrey_prepare_resume_from_failed_merge
        ;;
      archive)
        if is_true "$EARLGREY_ARCHIVE_INCOMPLETE"; then
          mkdir -p "${PROJECT_ROOT}/failed_runs/earlgrey"
          stamp="$(date +%Y%m%d_%H%M%S)"
          dest="${PROJECT_ROOT}/failed_runs/earlgrey/$(basename "$DIR_EARLGREY").${stamp}"
          echo "[WARN] Earl Grey output directory is incomplete or stale."
          echo "[WARN] Archiving: $DIR_EARLGREY -> $dest"
          mv "$DIR_EARLGREY" "$dest"
          mkdir -p "$DIR_EARLGREY"
          rm -f "$(repeat_checkpoint_file earlgrey)" "$(repeat_checkpoint_file earlgrey_plots)" "$(repeat_checkpoint_file collect_outputs)" "$(repeat_checkpoint_file non_repeats)"
        else
          echo "[ERROR] Earl Grey output directory is incomplete: $DIR_EARLGREY" >&2
          echo "[ERROR] Set EARLGREY_INCOMPLETE_POLICY=resume, archive, or error." >&2
          exit 1
        fi
        ;;
      error)
        echo "[ERROR] Earl Grey output directory is incomplete: $DIR_EARLGREY" >&2
        echo "[ERROR] Set EARLGREY_INCOMPLETE_POLICY=resume to keep upstream outputs and retry downstream stages." >&2
        exit 1
        ;;
    esac
  fi

  return 0
}

known_library_needs_build() {
  if repeat_checkpoint_ready known; then
    return 1
  fi

  return 0
}

migrate_earlgrey_v730_layout() {
  is_true "$EARLGREY_MIGRATE_LEGACY_LAYOUT" || return 0

  local stamp
  stamp="$(date +%Y%m%d_%H%M%S)"

  migrate_repeat_dir() {
    local source="$1"
    local destination="$2"
    if [[ -d "$source" && ! -e "$destination" ]]; then
      echo "[INFO] Migrating legacy repeat output: $source -> $destination"
      mv "$source" "$destination"
    elif [[ -d "$source" && -e "$destination" ]]; then
      echo "[WARN] Legacy and current repeat outputs both exist; preserving both: $source | $destination" >&2
    fi
  }

  migrate_repeat_dir "${PROJECT_ROOT}/06.earlgrey_annotation_only" "${PROJECT_ROOT}/05.earlgrey_annotation_only"
  migrate_repeat_dir "${PROJECT_ROOT}/07.final_outputs" "${PROJECT_ROOT}/06.final_outputs"
  migrate_repeat_dir "${PROJECT_ROOT}/08.non_repetitive_regions" "${PROJECT_ROOT}/07.non_repetitive_regions"

  if [[ -d "${PROJECT_ROOT}/05.repeatmasker" ]]; then
    mkdir -p "${PROJECT_ROOT}/retired_outputs"
    echo "[INFO] Archiving retired standalone RepeatMasker output."
    mv "${PROJECT_ROOT}/05.repeatmasker" \
      "${PROJECT_ROOT}/retired_outputs/standalone_repeatmasker.${stamp}"
  fi

  rm -f "${PROJECT_ROOT}/checkpoints/repeatmasker.done"
}

validate_toggles() {
  check_toggle AUTO_RESUME "$AUTO_RESUME"
  check_toggle AUTO_ARCHIVE_UNRECOVERABLE_RM_DIRS "$AUTO_ARCHIVE_UNRECOVERABLE_RM_DIRS"
  check_toggle TETRIMMER_ARCHIVE_INCOMPLETE "$TETRIMMER_ARCHIVE_INCOMPLETE"
  check_toggle REPEATMODELER_RECOVER "$REPEATMODELER_RECOVER"
  check_toggle REPEATMODELER_SANITIZE_FASTA_HEADERS "$REPEATMODELER_SANITIZE_FASTA_HEADERS"
  check_toggle EARLGREY_MIGRATE_LEGACY_LAYOUT "$EARLGREY_MIGRATE_LEGACY_LAYOUT"
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

# ============================================================
# earlgrey_enhanced_plotter.R
#
# Enhanced publication plotter for Earl Grey v7.3.x final outputs.
# Earl Grey native PDFs remain authoritative for its standard plots.
# This palette-fixed version follows Earl Grey autoLand.R/autoPie.R
# class collapsing and London Underground colour assignments where possible.
# v5 consolidates all enhanced Earl Grey plotting in one script:
#   - annotation-level K2P landscapes from the divergence GFF;
#   - family-level K2P summaries from <species>_divergence_summary_table.tsv;
#   - high-level genome/repeat composition from highLevelCount;
#   - family abundance summaries from familyLevelCount;
#   - robust log-table fallback, explicit schema audits, and publication exports.
# Family-mean K2P landscapes are explicitly labelled as family-level summaries and
# are never presented as substitutes for annotation-level Earl Grey landscapes.
# It keeps SINE/tRNA/7SL/5S/MIR in SINE while matching standalone RNA
# labels to Earl Grey's Other class:
#   <species>_summaryFiles/<species>.highLevelCount.txt
#   <species>_summaryFiles/<species>.familyLevelCount.txt
#   <species>_summaryFiles/<species>_divergence_summary_table.tsv
#   <species>_RepeatLandscape/<species>.filteredRepeats.withDivergence.gff
#
# Usage:
#   Rscript --vanilla earlgrey_enhanced_plotter.R \
#     --species 7777 \
#     --earlgrey_dir /path/to/7777_EarlGrey \
#     --outdir /path/to/enhanced_plots \
#     --dpi 600 \
#     --top_families 30 \
#     --top_superfamilies 12 \
#     --y_scale mb \
#     --axis_flip TRUE \
#     --divergence_summary /path/to/7777_divergence_summary_table.tsv
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

script_version <- "5.0.0"

usage_text <- function() {
  cat(
"earlgrey_enhanced_plotter.R v", script_version, "

Usage:
  Rscript --vanilla earlgrey_enhanced_plotter.R \\
    --species <name> \\
    --earlgrey_dir <EarlGrey output directory> \\
    --outdir <output directory> [options]

Required:
  --species <str>                 Species/sample identifier used by Earl Grey
  --earlgrey_dir <path>           Earl Grey <species>_EarlGrey directory
  --outdir <path>                 Output directory for enhanced plots/tables

Optional input overrides:
  --high_level <path>             highLevelCount.txt
  --family_level <path>           familyLevelCount.txt
  --divergence_summary <path>     <species>_divergence_summary_table.tsv
  --gff <path>                    filteredRepeats.withDivergence.gff
  --summary_table <path>          high-level TSV or markdown summary table
  --earlgrey_log <path>           Earl Grey log containing the final markdown table

Plot controls:
  --dpi <int>                     PNG resolution [600]
  --top_families <int>            Number of families in ranked plots [30]
  --top_superfamilies <int>       Number of superfamilies retained per class [12]
  --binwidth <float>              K2P bin width [0.01]
  --max_kimura <float>            Maximum K2P shown in landscapes [0.5]
  --y_scale <bp|mb|percent>       Landscape y-axis representation [mb]
  --axis_flip <TRUE|FALSE>        Reverse K2P axis to Earl Grey convention [TRUE]
  --composition_tolerance <float> Numeric composition tolerance [0.005]
  --debug_palette <TRUE|FALSE>    Use conspicuous fallback colour for unmapped classes [FALSE]

Other:
  --help                          Show this help
  --version                       Print script version

Main output groups:
  - annotation-level K2P landscapes from the divergence GFF
  - family-mean K2P landscapes and divergence plots from divergence_summary_table
  - genome/repeat composition from highLevelCount
  - family abundance and copy-number plots from familyLevelCount
  - normalized TSV summaries, schema audit, input manifest, plot manifest
", sep = "")
}

if ("--help" %in% args || "-h" %in% args) {
  usage_text()
  quit(save = "no", status = 0)
}
if ("--version" %in% args) {
  cat(script_version, "\n")
  quit(save = "no", status = 0)
}

arg_value <- function(flag, default = NA_character_) {
  idx <- which(args == flag)
  if (length(idx) == 0 || idx[1] == length(args)) return(default)
  args[idx[1] + 1]
}

is_true <- function(x) {
  tolower(as.character(x)) %in% c("true", "t", "yes", "y", "1")
}

species <- arg_value("--species")
earlgrey_dir <- arg_value("--earlgrey_dir")
outdir <- arg_value("--outdir")
dpi <- suppressWarnings(as.integer(arg_value("--dpi", "600")))
top_n <- suppressWarnings(as.integer(arg_value("--top_families", "30")))
binwidth <- suppressWarnings(as.numeric(arg_value("--binwidth", "0.01")))
max_kimura <- suppressWarnings(as.numeric(arg_value("--max_kimura", "0.5")))
debug_palette <- is_true(arg_value("--debug_palette", "FALSE"))
composition_tolerance <- suppressWarnings(as.numeric(arg_value("--composition_tolerance", "0.005")))
top_superfamilies <- suppressWarnings(as.integer(arg_value("--top_superfamilies", "12")))
y_scale <- tolower(arg_value("--y_scale", "mb"))
axis_flip <- is_true(arg_value("--axis_flip", "TRUE"))
summary_table_arg <- arg_value("--summary_table", NA_character_)
earlgrey_log_arg <- arg_value("--earlgrey_log", NA_character_)
high_level_arg <- arg_value("--high_level", NA_character_)
family_level_arg <- arg_value("--family_level", NA_character_)
divergence_summary_arg <- arg_value("--divergence_summary", NA_character_)
gff_arg <- arg_value("--gff", NA_character_)

if (is.na(species) || !nzchar(species)) stop("--species is required")
if (is.na(earlgrey_dir) || !dir.exists(earlgrey_dir)) stop("--earlgrey_dir must exist")
if (is.na(outdir) || !nzchar(outdir)) stop("--outdir is required")
if (!is.finite(dpi) || dpi < 1) dpi <- 600L
if (!is.finite(top_n) || top_n < 1) top_n <- 30L
if (!is.finite(binwidth) || binwidth <= 0) binwidth <- 0.01
if (!is.finite(max_kimura) || max_kimura <= 0) max_kimura <- 0.5
if (!is.finite(top_superfamilies) || top_superfamilies < 1) top_superfamilies <- 12L
if (!y_scale %in% c("bp", "mb", "percent")) stop("--y_scale must be one of: bp, mb, percent")
if (!is.finite(composition_tolerance) || composition_tolerance < 0 || composition_tolerance >= 1) {
  stop("--composition_tolerance must be numeric and in [0,1)")
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

required_pkgs <- c("ggplot2")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Missing required R package(s): ", paste(missing_pkgs, collapse = ", "))

has_data_table <- requireNamespace("data.table", quietly = TRUE)
has_scales <- requireNamespace("scales", quietly = TRUE)
has_grid <- requireNamespace("grid", quietly = TRUE)
has_patchwork <- requireNamespace("patchwork", quietly = TRUE)

read_tsv <- function(path) {
  if (is.na(path) || !file.exists(path) || file.info(path)$size == 0) {
    return(data.frame())
  }
  if (has_data_table) {
    return(as.data.frame(data.table::fread(
      path, sep = "\t", header = TRUE, quote = "", data.table = FALSE,
      showProgress = FALSE
    )))
  }
  utils::read.delim(
    path, sep = "\t", header = TRUE, check.names = FALSE,
    stringsAsFactors = FALSE, quote = "", comment.char = ""
  )
}

schema_audit <- data.frame(
  input = character(), path = character(), canonical = character(),
  source = character(), observed_columns = character(), status = character(),
  stringsAsFactors = FALSE
)

resolve_schema <- function(x, input, path, aliases) {
  observed <- paste(names(x), collapse = ";")
  if (nrow(x) == 0) stop("Empty required Earl Grey input: ", path)
  for (canonical in names(aliases)) {
    candidates <- aliases[[canonical]]
    source <- candidates[candidates %in% names(x)][1]
    if (is.na(source)) {
      stop(
        "Missing required field '", canonical, "' in ", input,
        ". Observed columns: ", observed
      )
    }
    if (source != canonical) names(x)[names(x) == source] <- canonical
    schema_audit <<- rbind(
      schema_audit,
      data.frame(
        input = input, path = path, canonical = canonical, source = source,
        observed_columns = observed, status = "mapped",
        stringsAsFactors = FALSE
      )
    )
  }
  x
}

read_gff_minimal <- function(path) {
  if (is.na(path) || !file.exists(path) || file.info(path)$size == 0) {
    return(data.frame())
  }
  if (has_data_table) {
    x <- data.table::fread(
      path, sep = "\t", header = FALSE, quote = "", comment.char = "",
      data.table = FALSE, select = c(3, 4, 5, 9), showProgress = FALSE
    )
    names(x) <- c("repeat_class", "start", "end", "metadata")
    return(x)
  }
  x <- utils::read.delim(
    path, sep = "\t", header = FALSE, quote = "", comment.char = "",
    fill = TRUE, stringsAsFactors = FALSE, col.names = paste0("V", seq_len(9))
  )
  data.frame(
    repeat_class = x$V3, start = x$V4, end = x$V5, metadata = x$V9,
    stringsAsFactors = FALSE
  )
}

num_col <- function(x) {
  suppressWarnings(as.numeric(gsub(",", "", as.character(x))))
}

label_bp <- function(x) {
  if (has_scales) return(scales::label_number(scale_cut = scales::cut_short_scale())(x))
  out <- ifelse(abs(x) >= 1e9, paste0(format(round(x / 1e9, 1), trim = TRUE), "G"),
         ifelse(abs(x) >= 1e6, paste0(format(round(x / 1e6, 1), trim = TRUE), "M"),
         ifelse(abs(x) >= 1e3, paste0(format(round(x / 1e3, 1), trim = TRUE), "K"),
                format(round(x, 0), trim = TRUE))))
  out
}

find_first <- function(pattern, prefer = character()) {
  files <- list.files(earlgrey_dir, pattern = pattern, recursive = TRUE, full.names = TRUE)
  files <- files[file.exists(files) & file.info(files)$size > 0]
  if (length(files) == 0) return(NA_character_)
  if (length(prefer) > 0) {
    score <- rep(length(prefer) + 1L, length(files))
    for (i in seq_along(prefer)) score[grepl(prefer[[i]], files, fixed = TRUE)] <- i
    files <- files[order(score, files)]
  } else {
    files <- sort(files)
  }
  files[[1]]
}

extract_attr <- function(x, key) {
  pattern <- paste0("(^|;)", key, "=([^;]+)")
  hit <- regexec(pattern, x, ignore.case = TRUE)
  out <- regmatches(x, hit)
  vapply(out, function(z) if (length(z) >= 3) z[[3]] else NA_character_, character(1))
}

clean_token <- function(x) {
  x <- as.character(x)
  x <- gsub("%20", " ", x, fixed = TRUE)
  x <- gsub("%2F", "/", x, fixed = TRUE)
  x <- gsub("%23", "#", x, fixed = TRUE)
  x <- sub("-nested$", "", x)
  x <- sub("_merged$", "", x)
  trimws(x)
}

parse_nested <- function(x) {
  grepl("-nested$", as.character(x), ignore.case = TRUE)
}

strip_nested_suffix <- function(x) {
  sub("-nested$", "", as.character(x), ignore.case = TRUE)
}

subclass_abbrev <- function(x) {
  # Earl Grey-style class collapsing with explicit SINE protection.
  # Important: do not use a broad "rna" regex before/after SINE matching,
  # because tRNA contains the substring "rna" and should remain SINE.
  z <- tolower(clean_token(x))
  z <- gsub("%2f", "/", z, fixed = TRUE)
  z <- gsub("%23", "#", z, fixed = TRUE)
  sep <- "(^|[/#;:_\\. -])"
  end <- "($|[/#;:_\\. -])"

  out <- rep("Unknown", length(z))

  # Earl Grey puts Unspecified, Retroposon, and standalone RNA in Other.
  other_pat <- paste0(
    "simple[_ -]?repeat|satellite|low[_ -]?complexity|microsatellite|",
    sep, "rna", end, "|^other|^unspecified$|^retroposon$|^t-?rna$|^rrna$|^snrna$",
    "|unknown_other|artefact|segmenta"
  )
  out[grepl(other_pat, z, perl = TRUE)] <- "Other"

  # Specific TE classes. These deliberately run after Other so SINE/tRNA can
  # override an accidental broad label, matching the biological class rather
  # than the generic RNA bucket.
  out[grepl("(^|[/#;:_\\. -])(dna|tir|cmc|hat|tc1|tcmar|pif|zisupton|maverick|polinton|piggybac|academ|crypton)($|[/#;:_\\. -])|dna[ _-]?transposon", z, perl = TRUE)] <- "DNA"
  out[grepl("(^|[/#;:_\\. -])(line|cr1|l1|l2|r2|rex|rte|jockey|i-jockey|ingi|tx1)($|[/#;:_\\. -])", z, perl = TRUE)] <- "LINE"
  out[grepl("(^|[/#;:_\\. -])(ltr|gypsy|copia|dirs|erv|erv1|ervl|pao|bel)($|[/#;:_\\. -])|ltr[ _-]?retrotransposon", z, perl = TRUE)] <- "LTR"

  # SINE/tRNA stays SINE via its SINE token; standalone tRNA stays Other.
  out[grepl("(^|[/#;:_\\. -])(sine|7sl|5s|mir|alu|b2|b4|id)($|[/#;:_\\. -])", z, perl = TRUE)] <- "SINE"

  out[grepl("penelope|(^|[/#;:_\\. -])ple($|[/#;:_\\. -])", z, perl = TRUE)] <- "PLE"
  out[grepl("rolling|helitron|(^|[/#;:_\\. -])rc($|[/#;:_\\. -])", z, perl = TRUE)] <- "RC"

  out
}

subclass_label <- function(x) {
  y <- subclass_abbrev(x)
  map <- c(
    DNA = "DNA Transposon", LINE = "LINE", LTR = "LTR Retrotransposon",
    SINE = "SINE", PLE = "Penelope", RC = "Rolling Circle",
    Other = "Other", Unknown = "Unknown"
  )
  unname(map[y])
}

class_base <- function(x) subclass_label(x)

superfamily_label <- function(te_type, repeat_class = NA_character_, metadata = NA_character_) {
  # This function intentionally avoids using raw family names as superfamily
  # labels. Earl Grey's native plotting scripts collapse repeat labels into a
  # small set of biologically interpretable classes, and only use explicit
  # class/superfamily tokens when they are present. If NAME lacks a "#" or the
  # token after "#" lacks "/", we fall back to the broad Earl Grey class.
  te_type <- clean_token(te_type)
  repeat_class <- clean_token(repeat_class)
  name_attr <- extract_attr(metadata, "NAME")

  n <- max(length(te_type), length(repeat_class), length(name_attr))
  te_type <- rep(te_type, length.out = n)
  repeat_class <- rep(repeat_class, length.out = n)
  name_attr <- rep(name_attr, length.out = n)

  explicit_type <- rep(NA_character_, n)
  has_hash <- !is.na(name_attr) & nzchar(name_attr) & grepl("#", name_attr, fixed = TRUE)
  explicit_type[has_hash] <- sub("^.*#", "", name_attr[has_hash])
  explicit_type <- clean_token(explicit_type)
  explicit_type <- sub(";.*$", "", explicit_type)

  # Source used for broad Earl Grey class assignment.
  class_source <- ifelse(!is.na(explicit_type) & nzchar(explicit_type), explicit_type,
                         ifelse(!is.na(te_type) & nzchar(te_type), te_type, repeat_class))
  broad <- subclass_abbrev(class_source)

  # Only derive superfamily from explicit class/superfamily strings, e.g.
  # DNA/TcMar, LINE/CR1, LTR/Gypsy. Do not treat family names such as
  # ACADEM-1_SP or tandem-repeat labels as superfamilies.
  sf <- rep(NA_character_, n)
  has_slash <- !is.na(explicit_type) & nzchar(explicit_type) & grepl("/", explicit_type, fixed = TRUE)
  sf[has_slash] <- sub("^.*/", "", explicit_type[has_slash])
  sf <- clean_token(sf)

  # Normalize specific superfamily labels.
  sf[grepl("cmc", sf, ignore.case = TRUE)] <- "CMC"
  sf[grepl("hat", sf, ignore.case = TRUE)] <- "hAT"
  sf[grepl("pif|harbinger", sf, ignore.case = TRUE)] <- "PIF"
  sf[grepl("tcmar|tc1|mariner", sf, ignore.case = TRUE)] <- "TcMar"
  sf[grepl("zisupton", sf, ignore.case = TRUE)] <- "Zisupton"
  sf[grepl("maverick|polinton", sf, ignore.case = TRUE)] <- "Maverick"
  sf[grepl("crypton", sf, ignore.case = TRUE)] <- "Crypton"
  sf[grepl("piggy", sf, ignore.case = TRUE)] <- "piggyBac"
  sf[grepl("academ", sf, ignore.case = TRUE)] <- "Academ"

  sf[grepl("^cr1$", sf, ignore.case = TRUE)] <- "CR1"
  sf[grepl("^l1$", sf, ignore.case = TRUE)] <- "L1"
  sf[grepl("^l2$", sf, ignore.case = TRUE)] <- "L2"
  sf[grepl("^r2$", sf, ignore.case = TRUE)] <- "R2"
  sf[grepl("^rex$", sf, ignore.case = TRUE)] <- "Rex"
  sf[grepl("^rte$", sf, ignore.case = TRUE)] <- "RTE"
  sf[grepl("jockey", sf, ignore.case = TRUE)] <- "Jockey"
  sf[grepl("^i$|i-jockey", sf, ignore.case = TRUE)] <- "I-Jockey"
  sf[grepl("tx1", sf, ignore.case = TRUE)] <- "Tx1"
  sf[grepl("ingi", sf, ignore.case = TRUE)] <- "Ingi"

  sf[grepl("dirs", sf, ignore.case = TRUE)] <- "DIRS"
  sf[grepl("erv1", sf, ignore.case = TRUE)] <- "ERV1"
  sf[grepl("ervl", sf, ignore.case = TRUE)] <- "ERVL"
  sf[grepl("^erv$", sf, ignore.case = TRUE)] <- "ERV"
  sf[grepl("gypsy", sf, ignore.case = TRUE)] <- "Gypsy"
  sf[grepl("copia", sf, ignore.case = TRUE)] <- "Copia"
  sf[grepl("pao|bel", sf, ignore.case = TRUE)] <- "Pao"

  sf[grepl("trna|t-rna", sf, ignore.case = TRUE)] <- "tRNA"
  sf[grepl("7sl", sf, ignore.case = TRUE)] <- "7SL"
  sf[grepl("5s", sf, ignore.case = TRUE)] <- "5S"
  sf[grepl("mir", sf, ignore.case = TRUE)] <- "MIR"
  sf[grepl("^sine$", sf, ignore.case = TRUE)] <- "SINE"
  sf[grepl("^alu$", sf, ignore.case = TRUE)] <- "Alu"
  sf[grepl("^b2$", sf, ignore.case = TRUE)] <- "B2"
  sf[grepl("^b4$", sf, ignore.case = TRUE)] <- "B4"

  sf[grepl("^rc$|rolling|helitron", sf, ignore.case = TRUE)] <- "Helitron"
  sf[grepl("^ple$|penelope", sf, ignore.case = TRUE)] <- "Penelope"

  # If no explicit superfamily was found, use the broad Earl Grey class.
  broad_fallback <- c(
    DNA = "DNA", LINE = "LINE", LTR = "LTR", SINE = "SINE",
    PLE = "Penelope", RC = "Helitron", Other = "Other", Unknown = "Unclassified"
  )
  use_fallback <- is.na(sf) | !nzchar(sf) | sf == "NA"
  sf[use_fallback] <- unname(broad_fallback[broad[use_fallback]])
  sf[is.na(sf) | !nzchar(sf)] <- "Unclassified"

  sf
}

subclass_levels <- c("DNA", "LINE", "LTR", "SINE", "PLE", "RC", "Other", "Unknown")
subclass_label_levels <- c(
  "DNA Transposon", "LINE", "LTR Retrotransposon", "SINE",
  "Penelope", "Rolling Circle", "Other", "Unknown"
)

# Earl Grey native class colours. These reproduce autoLand.R/autoPie.R:
# DNA=#E32017, Rolling Circle=#EE7C0E, Penelope=#7156A5,
# LINE=#0098D4, SINE=#9B0056, LTR=#00782A,
# Other=#F3A9BB, Unclassified=#A0A5A9, Non-Repeat=#000000.
# Extra aliases are deliberately included so ggplot never falls back to grey
# just because labels differ between Earl Grey tables, GFF fields, and enhanced plots.
class_palette <- c(
  "DNA" = "#E32017", "DNA Transposon" = "#E32017", "TIR" = "#E32017",
  "Rolling Circle" = "#EE7C0E", "RC" = "#EE7C0E", "Helitron" = "#EE7C0E",
  "Penelope" = "#7156A5", "PLE" = "#7156A5",
  "LINE" = "#0098D4",
  "SINE" = "#9B0056",
  "LTR" = "#00782A", "LTR Retrotransposon" = "#00782A",
  "Other (Simple Repeat, Microsatellite, RNA)" = "#F3A9BB",
  "Other" = "#F3A9BB", "Other simple or low-complexity" = "#F3A9BB",
  "Unclassified" = "#A0A5A9", "Unknown" = "#A0A5A9",
  "Non-Repeat" = "#000000"
)

# Superfamily colours. The named superfamilies use related palettes within
# each TE class, while broad fallback labels reuse Earl Grey native class colours.
# Additional aliases are used for robust parsing of RepeatMasker/Earl Grey labels.
superfamily_palette <- c(
  "CMC" = "#FB8072", "hAT" = "#B3DE69", "PIF" = "#66C2A5", "TcMar" = "#00B0F6",
  "Tc1" = "#00B0F6", "Tc1-Mariner" = "#00B0F6", "Zisupton" = "#E76BF3",
  "Maverick" = "#FDB462", "Polinton" = "#FDB462", "Crypton" = "#FCCDE5",
  "piggyBac" = "#BC80BD", "Academ" = "#CCEBC5", "DNA" = "#E32017", "DNA Transposon" = "#E32017",
  "CR1" = "#084594", "L1" = "#2171B5", "L2" = "#4292C6", "LINE" = "#6BAED6",
  "R2" = "#9ECAE1", "Rex" = "#C6DBEF", "RTE" = "#EFF3FF", "Jockey" = "#08519C",
  "I" = "#3182BD", "I-Jockey" = "#3182BD", "Tx1" = "#BDD7E7", "Ingi" = "#6BAED6",
  "DIRS" = "#005A32", "ERV" = "#238B45", "ERV1" = "#41AB5D", "ERVL" = "#74C476",
  "Gypsy" = "#A1D99B", "Copia" = "#31A354", "LTR" = "#C7E9C0", "LTR Retrotransposon" = "#00782A",
  "Pao" = "#EDF8E9", "BEL" = "#EDF8E9",
  "SINE" = "#9B0056", "tRNA" = "#9B0056", "7SL" = "#C51B7D", "5S" = "#DE77AE", "MIR" = "#F1B6DA", "Alu" = "#8E0152", "B2" = "#C51B7D", "B4" = "#DE77AE",
  "Penelope" = "#ECEAF7", "PLE" = "#7156A5",
  "Helitron" = "#F6DFC8", "RC" = "#EE7C0E", "Rolling Circle" = "#EE7C0E",
  "Other (Simple Repeat, Microsatellite, RNA)" = "#F3A9BB", "Other" = "#F3A9BB",
  "Unclassified" = "#A0A5A9", "Unknown" = "#A0A5A9"
)

na_debug_color <- if (debug_palette) "#ff00ff" else "#bdbdbd"

make_extra_palette <- function(labels, existing = character()) {
  labels <- labels[!is.na(labels) & nzchar(labels)]
  labels <- setdiff(sort(unique(labels)), existing)
  if (length(labels) == 0) return(setNames(character(), character()))
  cols <- if ("hcl.colors" %in% getNamespaceExports("grDevices")) {
    grDevices::hcl.colors(length(labels), palette = "Dark 3")
  } else {
    grDevices::rainbow(length(labels), s = 0.65, v = 0.85)
  }
  names(cols) <- labels
  cols
}

complete_palette <- function(values, palette, label = "palette") {
  values <- sort(unique(as.character(values)))
  values <- values[!is.na(values) & nzchar(values)]
  missing <- setdiff(values, names(palette))
  if (length(missing) > 0) {
    warning("Adding automatic colours for ", label, ": ", paste(missing, collapse = ", "))
    palette <- c(palette, make_extra_palette(missing, names(palette)))
  }
  palette
}

audit_palette <- function(values, palette, label, outdir) {
  values <- sort(unique(as.character(values)))
  values <- values[!is.na(values) & nzchar(values)]
  missing <- setdiff(values, names(palette))
  audit <- data.frame(
    label = label,
    value = values,
    in_palette = values %in% names(palette),
    stringsAsFactors = FALSE
  )
  utils::write.table(
    audit,
    file.path(outdir, paste0("palette_audit_", label, ".tsv")),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  if (length(missing) > 0) {
    warning("Missing palette colors for ", label, ": ", paste(missing, collapse = ", "))
  }
  invisible(missing)
}

enhanced_theme <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(linewidth = 0.25, colour = "grey90"),
      panel.border = ggplot2::element_rect(linewidth = 0.5, colour = "grey45"),
      axis.title = ggplot2::element_text(face = "bold", colour = "grey15"),
      axis.text = ggplot2::element_text(colour = "grey20"),
      plot.title = ggplot2::element_text(hjust = 0, face = "bold", size = ggplot2::rel(1.22), colour = "grey10"),
      plot.subtitle = ggplot2::element_text(hjust = 0, colour = "grey30", margin = ggplot2::margin(b = 8)),
      plot.caption = ggplot2::element_text(hjust = 0, colour = "grey40", size = ggplot2::rel(0.82), margin = ggplot2::margin(t = 8)),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.text = ggplot2::element_text(size = ggplot2::rel(0.9)),
      legend.key.size = if (has_grid) grid::unit(0.42, "cm") else ggplot2::rel(0.8),
      strip.background = ggplot2::element_rect(fill = "grey96", colour = "grey75", linewidth = 0.4),
      strip.text = ggplot2::element_text(face = "bold", colour = "grey15"),
      plot.margin = ggplot2::margin(10, 14, 10, 10)
    )
}

species_label <- gsub("_", " ", species)

format_mb <- function(x, accuracy = 0.1) {
  if (has_scales) return(scales::label_number(scale = 1e-6, accuracy = accuracy)(x))
  format(round(x / 1e6, ifelse(accuracy < 1, 1, 0)), trim = TRUE)
}

format_pct <- function(x, accuracy = 0.1) {
  if (has_scales) return(scales::label_number(accuracy = accuracy, suffix = "%")(x))
  paste0(format(round(x, 1), trim = TRUE), "%")
}

kimura_scale <- function(reverse = axis_flip, axis_name = "Kimura 2-parameter divergence (K2P)") {
  breaks <- seq(0, max_kimura, by = 0.05)
  minor <- seq(0, max_kimura, by = binwidth)
  if (reverse) {
    ggplot2::scale_x_reverse(
      limits = c(max_kimura, 0), breaks = breaks, minor_breaks = minor,
      expand = ggplot2::expansion(mult = c(0, 0)), name = axis_name
    )
  } else {
    ggplot2::scale_x_continuous(
      limits = c(0, max_kimura), breaks = breaks, minor_breaks = minor,
      expand = ggplot2::expansion(mult = c(0, 0)), name = axis_name
    )
  }
}

landscape_y_scale <- function(genome_size = NA_real_) {
  if (y_scale == "percent" && is.finite(genome_size) && genome_size > 0) {
    return(ggplot2::scale_y_continuous(
      name = "Genome coverage per K2P bin (%)",
      labels = function(x) format_pct(100 * x / genome_size),
      expand = ggplot2::expansion(mult = c(0, 0.05))
    ))
  }
  if (y_scale == "bp") {
    return(ggplot2::scale_y_continuous(
      name = "Repeat sequence (bp)", labels = label_bp,
      expand = ggplot2::expansion(mult = c(0, 0.05))
    ))
  }
  ggplot2::scale_y_continuous(
    name = "Repeat sequence (Mb)",
    labels = function(x) format_mb(x),
    expand = ggplot2::expansion(mult = c(0, 0.05))
  )
}

save_plot <- function(plot, stem, width = 9, height = 6) {
  outputs <- file.path(outdir, paste0(stem, c(".png", ".pdf", ".svg")))
  ggplot2::ggsave(outputs[[1]], plot = plot, width = width, height = height,
                  units = "in", dpi = dpi, limitsize = FALSE, bg = "white")
  pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else "pdf"
  ggplot2::ggsave(outputs[[2]], plot = plot, width = width, height = height,
                  units = "in", device = pdf_device, limitsize = FALSE, bg = "white")
  svg_device <- if (requireNamespace("svglite", quietly = TRUE)) svglite::svglite else "svg"
  ggplot2::ggsave(outputs[[3]], plot = plot, width = width, height = height,
                  units = "in", device = svg_device, limitsize = FALSE, bg = "white")
  data.frame(plot = stem, format = c("png", "pdf", "svg"), file = outputs,
             stringsAsFactors = FALSE)
}

read_markdown_table <- function(path) {
  if (is.na(path) || !file.exists(path) || file.info(path)$size == 0) return(data.frame())
  lines <- readLines(path, warn = FALSE)
  header_idx <- grep("\\|\\s*TE Classification\\s*\\|.*Coverage \\(bp\\).*% Genome Coverage", lines)
  if (!length(header_idx)) return(data.frame())
  i <- tail(header_idx, 1)
  block <- lines[i:length(lines)]
  block <- block[grepl("^\\s*\\|", block)]
  if (length(block) < 3) return(data.frame())
  # stop after the contiguous markdown table
  source_lines <- lines[i:length(lines)]
  keep <- c()
  for (ln in source_lines) {
    if (!grepl("^\\s*\\|", ln)) {
      if (length(keep)) break else next
    }
    keep <- c(keep, ln)
  }
  if (length(keep) < 3) return(data.frame())
  split_line <- function(x) {
    x <- sub("^\\s*\\|", "", x)
    x <- sub("\\|\\s*$", "", x)
    trimws(strsplit(x, "\\|", fixed = FALSE)[[1]])
  }
  header <- split_line(keep[[1]])
  data_lines <- keep[-c(1, 2)]
  rows <- lapply(data_lines, split_line)
  rows <- rows[vapply(rows, length, integer(1)) == length(header)]
  if (!length(rows)) return(data.frame())
  out <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE, check.names = FALSE)
  names(out) <- header
  out
}

read_summary_table <- function() {
  # Explicit table takes precedence, then native highLevelCount, then Earl Grey log.
  if (!is.na(summary_table_arg) && file.exists(summary_table_arg)) {
    x <- read_tsv(summary_table_arg)
    if (nrow(x)) return(list(data = x, path = summary_table_arg, source = "explicit_table"))
    x <- read_markdown_table(summary_table_arg)
    if (nrow(x)) return(list(data = x, path = summary_table_arg, source = "explicit_markdown"))
  }
  if (!is.na(high_file) && file.exists(high_file)) {
    x <- read_tsv(high_file)
    if (nrow(x)) return(list(data = x, path = high_file, source = "highLevelCount"))
  }
  log_file <- earlgrey_log_arg
  if (is.na(log_file) || !file.exists(log_file)) {
    logs <- list.files(earlgrey_dir, pattern = "\\.log$", recursive = TRUE, full.names = TRUE)
    logs <- logs[file.exists(logs) & file.info(logs)$size > 0]
    if (length(logs)) log_file <- logs[order(file.info(logs)$mtime, decreasing = TRUE)][1]
  }
  x <- read_markdown_table(log_file)
  if (nrow(x)) return(list(data = x, path = log_file, source = "earlgrey_log"))
  list(data = data.frame(), path = NA_character_, source = "missing")
}

bin_landscape <- function(gff) {
  if (!nrow(gff)) return(data.frame())
  gff$kimura_bin <- round(gff$kimura80 / binwidth) * binwidth
  gff$kimura_bin <- pmin(max_kimura, pmax(0, gff$kimura_bin))
  agg <- stats::aggregate(
    width_bp ~ subclass + subclass_label + superfamily + kimura_bin,
    data = gff, FUN = sum, na.rm = TRUE
  )
  names(agg)[names(agg) == "width_bp"] <- "bp"
  agg
}

weighted_mean_safe <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  stats::weighted.mean(x[ok], w[ok])
}

weighted_median_safe <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  x <- x[ok]
  w <- w[ok]
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  x[which(cumsum(w) >= sum(w) / 2)[1]]
}

bin_numeric <- function(x, width = binwidth, max_value = max_kimura) {
  out <- round(x / width) * width
  pmin(max_value, pmax(0, out))
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  stats::sd(x)
}

prepare_gff_landscape <- function(gff_file) {
  if (is.na(gff_file) || !file.exists(gff_file) || file.info(gff_file)$size == 0) {
    return(data.frame())
  }
  gff <- read_gff_minimal(gff_file)
  if (nrow(gff) == 0) return(data.frame())
  name_attr <- extract_attr(gff$metadata, "NAME")
  te_type <- ifelse(!is.na(name_attr) & grepl("#", name_attr), sub("^.*#", "", name_attr), gff$repeat_class)

  gff$start <- num_col(gff$start)
  gff$end <- num_col(gff$end)
  gff$width_bp <- gff$end - gff$start + 1
  gff <- gff[
    is.finite(gff$start) & is.finite(gff$end) &
      gff$start >= 1 & gff$end >= gff$start,
    , drop = FALSE
  ]
  kimura80 <- extract_attr(gff$metadata, "kimura80")
  kimura80_alt <- extract_attr(gff$metadata, "Kimura80")
  divergence <- extract_attr(gff$metadata, "divergence")
  kimura_value <- ifelse(!is.na(kimura80), kimura80, ifelse(!is.na(kimura80_alt), kimura80_alt, divergence))
  gff$kimura80 <- suppressWarnings(as.numeric(kimura_value))
  percent_scale <- is.finite(gff$kimura80) & gff$kimura80 > 1
  gff$kimura80[percent_scale] <- gff$kimura80[percent_scale] / 100
  gff$subclass <- subclass_abbrev(te_type)
  gff$subclass_label <- subclass_label(te_type)
  gff$superfamily <- superfamily_label(te_type, gff$repeat_class, gff$metadata)

  gff <- gff[
    is.finite(gff$kimura80) & gff$kimura80 >= 0 & gff$kimura80 <= max_kimura &
      is.finite(gff$width_bp) & gff$width_bp > 0,
    , drop = FALSE
  ]
  if (nrow(gff) == 0) return(data.frame())

  gff$subclass <- factor(gff$subclass, levels = subclass_levels)
  gff$subclass_label <- factor(gff$subclass_label, levels = subclass_label_levels)
  gff
}

manifest <- data.frame(plot = character(), format = character(), file = character(), stringsAsFactors = FALSE)

high_file <- if (!is.na(high_level_arg) && file.exists(high_level_arg)) high_level_arg else
  find_first(paste0("(^|/)", species, "\\.highLevelCount\\.txt$"), prefer = c("_summaryFiles"))
family_file <- if (!is.na(family_level_arg) && file.exists(family_level_arg)) family_level_arg else
  find_first(paste0("(^|/)", species, "\\.familyLevelCount\\.txt$"), prefer = c("_summaryFiles"))
div_file <- if (!is.na(divergence_summary_arg) && file.exists(divergence_summary_arg)) divergence_summary_arg else
  find_first(
    paste0("(^|/)", species, "(_divergence_summary_table|_summary_table)\\.(tsv|txt)$"),
    prefer = c("_summaryFiles", "_RepeatLandscape")
  )
gff_file <- if (!is.na(gff_arg) && file.exists(gff_arg)) gff_arg else
  find_first(
    paste0("(^|/)", species, "\\.filteredRepeats\\.withDivergence\\.gff$"),
    prefer = c("_RepeatLandscape", "_summaryFiles", "looseMerge")
  )
if (is.na(gff_file)) {
  gff_file <- find_first(
    paste0("(^|/)", species, "\\.filteredRepeats\\.gff$"),
    prefer = c("_summaryFiles", "looseMerge")
  )
}

landscape_gff <- prepare_gff_landscape(gff_file)
schema_audit <- rbind(
  schema_audit,
  data.frame(
    input = "gff", path = gff_file,
    canonical = c("type", "start", "end", "attributes"),
    source = c("column_3", "column_4", "column_5", "column_9"),
    observed_columns = "GFF3 columns 3,4,5,9",
    status = if (nrow(landscape_gff) > 0) "mapped" else "unreadable",
    stringsAsFactors = FALSE
  )
)
if (nrow(landscape_gff) > 0) {
  class_palette <- complete_palette(c(landscape_gff$subclass_label, landscape_gff$subclass), class_palette, "TE subclass")
  superfamily_palette <- complete_palette(landscape_gff$superfamily, superfamily_palette, "TE superfamily")
  audit_palette(landscape_gff$subclass_label, class_palette, "subclass", outdir)
  audit_palette(landscape_gff$superfamily, superfamily_palette, "superfamily", outdir)
}

# ============================================================
# Enhanced summary/composition plots
# ============================================================
summary_info <- read_summary_table()
high <- data.frame()
genome_size <- NA_real_
summary_source <- summary_info$source
summary_source_path <- summary_info$path

if (nrow(summary_info$data) > 0) {
  high <- resolve_schema(
    summary_info$data,
    "high_level",
    summary_info$path,
    list(
      "TE Classification" = c("TE Classification", "TE_Classification", "Classification"),
      "Coverage (bp)" = c("Coverage (bp)", "Coverage_bp", "Coverage"),
      "Copy Number" = c("Copy Number", "Copy_Number", "Copies"),
      "% Genome Coverage" = c("% Genome Coverage", "Genome Coverage (%)", "Percent Genome Coverage")
    )
  )

  # Optional columns present in Earl Grey's printed kable summary.
  genome_col <- c("Genome Size", "Genome_Size", "Genome size")
  genome_col <- genome_col[genome_col %in% names(high)][1]
  family_count_col <- c("TE Family Count", "TE_Family_Count", "Family Count")
  family_count_col <- family_count_col[family_count_col %in% names(high)][1]

  high$coverage_bp <- num_col(high[["Coverage (bp)"]])
  high$genome_pct <- num_col(gsub("%", "", high[["% Genome Coverage"]], fixed = TRUE))
  high$copy_number_num <- num_col(high[["Copy Number"]])
  high$family_count_num <- if (!is.na(family_count_col)) num_col(high[[family_count_col]]) else NA_real_
  high$is_nested <- parse_nested(high[["TE Classification"]])
  high$class_token <- strip_nested_suffix(high[["TE Classification"]])
  high$class_base <- class_base(high$class_token)
  high$class_base[high[["TE Classification"]] == "Non-Repeat"] <- "Non-Repeat"
  high$is_total <- grepl("Total Interspersed Repeat", high[["TE Classification"]], ignore.case = TRUE)
  high$is_nonrepeat <- high[["TE Classification"]] == "Non-Repeat"

  if (!is.na(genome_col)) {
    genome_candidates <- num_col(high[[genome_col]])
    genome_candidates <- genome_candidates[is.finite(genome_candidates) & genome_candidates > 0]
    if (length(genome_candidates)) genome_size <- stats::median(genome_candidates)
  }
  if (!is.finite(genome_size)) {
    stable <- is.finite(high$coverage_bp) & high$coverage_bp > 0 & is.finite(high$genome_pct) & high$genome_pct > 0.1
    candidates <- high$coverage_bp[stable] * 100 / high$genome_pct[stable]
    candidates <- candidates[is.finite(candidates) & candidates > 0]
    if (length(candidates)) genome_size <- stats::median(candidates)
  }

  # Preserve the parsed table as a pipeline artifact.
  parsed_summary <- high
  utils::write.table(parsed_summary, file.path(outdir, paste0(species, "_earlgrey_summary_table.tsv")),
                     sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")

  repeat_rows <- high[!high$is_total & !high$is_nonrepeat & is.finite(high$coverage_bp) & high$coverage_bp >= 0, , drop = FALSE]
  unique_rows <- repeat_rows[!repeat_rows$is_nested, , drop = FALSE]
  nested_rows <- repeat_rows[repeat_rows$is_nested, , drop = FALSE]
  nonrepeat_row <- high[high$is_nonrepeat & is.finite(high$coverage_bp), , drop = FALSE]

  total_repeat_bp <- if (nrow(unique_rows)) sum(unique_rows$coverage_bp, na.rm = TRUE) else NA_real_
  nested_repeat_bp <- if (nrow(nested_rows)) sum(nested_rows$coverage_bp, na.rm = TRUE) else 0
  nonrepeat_bp <- if (nrow(nonrepeat_row)) sum(nonrepeat_row$coverage_bp, na.rm = TRUE) else if (is.finite(genome_size)) genome_size - total_repeat_bp else NA_real_
  repeat_pct <- if (is.finite(genome_size) && genome_size > 0) 100 * total_repeat_bp / genome_size else sum(unique_rows$genome_pct, na.rm = TRUE)
  family_total <- if (any(is.finite(unique_rows$family_count_num))) sum(unique_rows$family_count_num, na.rm = TRUE) else NA_real_
  copy_total <- if (any(is.finite(unique_rows$copy_number_num))) sum(unique_rows$copy_number_num, na.rm = TRUE) else NA_real_

  summary_subtitle <- paste0(
    if (is.finite(genome_size)) sprintf("Genome %.2f Gb", genome_size / 1e9) else "Genome size unavailable",
    "  ·  ", if (is.finite(repeat_pct)) sprintf("%.1f%% repetitive", repeat_pct) else "repeat fraction unavailable",
    "  ·  ", if (is.finite(family_total) && family_total > 0) paste(format(round(family_total), big.mark = ","), "families") else "family count unavailable"
  )

  # 1. Genome composition: one proportional bar, explicitly repeat vs non-repeat.
  genome_comp <- data.frame(
    component = c("Annotated repeats", "Non-repeat"),
    coverage_bp = c(total_repeat_bp, nonrepeat_bp),
    stringsAsFactors = FALSE
  )
  genome_comp <- genome_comp[is.finite(genome_comp$coverage_bp) & genome_comp$coverage_bp >= 0, , drop = FALSE]
  if (nrow(genome_comp) == 2) {
    genome_comp$percent <- 100 * genome_comp$coverage_bp / sum(genome_comp$coverage_bp)
    genome_palette <- c("Annotated repeats" = "#666666", "Non-repeat" = "#E6E6E6")
    p <- ggplot2::ggplot(genome_comp, ggplot2::aes(x = "Genome", y = percent, fill = component)) +
      ggplot2::geom_col(width = 0.62, colour = "white", linewidth = 0.5) +
      ggplot2::geom_text(
        ggplot2::aes(label = ifelse(percent >= 4, sprintf("%.1f%%", percent), "")),
        position = ggplot2::position_stack(vjust = 0.5), fontface = "bold", size = 4
      ) +
      ggplot2::scale_fill_manual(values = genome_palette) +
      ggplot2::scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = c(0, 0)) +
      ggplot2::labs(
        title = paste("Genome composition of", species_label), subtitle = summary_subtitle,
        x = NULL, y = "Genome coverage (%)", fill = NULL,
        caption = "Repeat coverage follows Earl Grey's non-nested class totals; nested intervals are summarized separately."
      ) + enhanced_theme() +
      ggplot2::theme(legend.position = "bottom", axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
    manifest <- rbind(manifest, save_plot(p, "enhanced_genome_repeat_fraction", 8, 5.4))
  }

  # 2. Class composition as genome coverage (%), preserving Earl Grey class colours.
  # Aggregate each metric independently so optional all-NA columns (especially
  # TE Family Count in older Earl Grey outputs) do not remove valid classes.
  class_groups <- split(unique_rows, unique_rows$class_base)
  class_summary <- do.call(rbind, lapply(class_groups, function(d) {
    data.frame(
      class_base = d$class_base[1],
      coverage_bp = sum(d$coverage_bp, na.rm = TRUE),
      genome_pct = sum(d$genome_pct, na.rm = TRUE),
      copy_number_num = if (any(is.finite(d$copy_number_num))) sum(d$copy_number_num, na.rm = TRUE) else NA_real_,
      family_count_num = if (any(is.finite(d$family_count_num))) sum(d$family_count_num, na.rm = TRUE) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
  class_summary <- class_summary[class_summary$class_base != "Non-Repeat", , drop = FALSE]
  class_summary$class_base <- factor(class_summary$class_base, levels = subclass_label_levels)
  class_summary <- class_summary[order(class_summary$genome_pct), , drop = FALSE]
  if (nrow(class_summary)) {
    p <- ggplot2::ggplot(class_summary, ggplot2::aes(x = class_base, y = genome_pct, fill = class_base)) +
      ggplot2::geom_col(width = 0.72) +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f%%", genome_pct)), hjust = -0.12, size = 3.3) +
      ggplot2::coord_flip(clip = "off") +
      ggplot2::scale_fill_manual(values = class_palette, drop = FALSE, na.value = na_debug_color) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.13))) +
      ggplot2::labs(
        title = paste("Repeat composition of", species_label), subtitle = summary_subtitle,
        x = NULL, y = "Genome coverage (%)", fill = "Repeat subclass",
        caption = "Values exclude Earl Grey '-nested' rows to avoid double-counting genomic coverage."
      ) + enhanced_theme() + ggplot2::theme(legend.position = "none")
    manifest <- rbind(manifest, save_plot(p, "enhanced_repeat_class_composition", 9.5, 5.8))

    # 3. One 100% stacked proportional bar showing the repeatome itself.
    prop <- class_summary
    prop$repeat_pct <- 100 * prop$coverage_bp / sum(prop$coverage_bp, na.rm = TRUE)
    p <- ggplot2::ggplot(prop, ggplot2::aes(x = "Repeatome", y = repeat_pct, fill = class_base)) +
      ggplot2::geom_col(width = 0.62, colour = "white", linewidth = 0.3) +
      ggplot2::scale_fill_manual(values = class_palette, drop = FALSE, na.value = na_debug_color) +
      ggplot2::scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = c(0, 0)) +
      ggplot2::labs(
        title = paste("Proportional repeatome of", species_label), subtitle = "Relative contribution of each repeat subclass to annotated repeat sequence",
        x = NULL, y = "Share of annotated repeats (%)", fill = "Repeat subclass",
        caption = sprintf("Annotated repeat sequence: %.1f Mb", total_repeat_bp / 1e6)
      ) + enhanced_theme() +
      ggplot2::theme(legend.position = "bottom", axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank()) +
      ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2, byrow = TRUE))
    manifest <- rbind(manifest, save_plot(p, "enhanced_repeatome_proportional_bar", 9.5, 5.5))
  }

  # 4. Nested versus non-nested class burden.
  if (nrow(repeat_rows)) {
    nn <- repeat_rows
    nn$state <- ifelse(nn$is_nested, "Nested", "Primary annotation")
    nn_sum <- stats::aggregate(coverage_bp ~ class_base + state, data = nn, FUN = sum, na.rm = TRUE)
    nn_sum$class_base <- factor(nn_sum$class_base, levels = subclass_label_levels)
    p <- ggplot2::ggplot(nn_sum, ggplot2::aes(x = class_base, y = coverage_bp / 1e6, fill = state)) +
      ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.76), width = 0.68) +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values = c("Primary annotation" = "grey35", "Nested" = "grey75")) +
      ggplot2::labs(
        title = paste("Nested repeat burden of", species_label),
        subtitle = sprintf("Nested repeat intervals: %.1f Mb (reported separately from unique genome composition)", nested_repeat_bp / 1e6),
        x = NULL, y = "Annotated sequence (Mb)", fill = NULL,
        caption = "Nested values may overlap primary repeat annotations and should not be added to genome-wide coverage totals."
      ) + enhanced_theme() + ggplot2::theme(legend.position = "bottom")
    manifest <- rbind(manifest, save_plot(p, "enhanced_nested_vs_primary_repeats", 9.5, 6))
  }

  # 5. Copy number and TE family count summary from the printed Earl Grey table.
  metrics_long <- rbind(
    data.frame(class_base = class_summary$class_base, metric = "Copy number", value = class_summary$copy_number_num),
    data.frame(class_base = class_summary$class_base, metric = "TE family count", value = class_summary$family_count_num)
  )
  metrics_long <- metrics_long[is.finite(metrics_long$value) & metrics_long$value > 0, , drop = FALSE]
  if (nrow(metrics_long)) {
    p <- ggplot2::ggplot(metrics_long, ggplot2::aes(x = class_base, y = value, fill = class_base)) +
      ggplot2::geom_col(width = 0.72) +
      ggplot2::coord_flip() +
      ggplot2::facet_wrap(~ metric, scales = "free_x", ncol = 2) +
      ggplot2::scale_fill_manual(values = class_palette, drop = FALSE, na.value = na_debug_color) +
      ggplot2::scale_y_continuous(labels = label_bp, expand = ggplot2::expansion(mult = c(0, 0.06))) +
      ggplot2::labs(
        title = paste("Repeat abundance and library diversity of", species_label),
        subtitle = "Copy number quantifies annotated occurrences; family count summarizes represented TE families",
        x = NULL, y = NULL, fill = "Repeat subclass"
      ) + enhanced_theme() + ggplot2::theme(legend.position = "none")
    manifest <- rbind(manifest, save_plot(p, "enhanced_repeat_copy_and_family_counts", 10, 6))
  }

  utils::write.table(class_summary, file.path(outdir, paste0(species, "_repeat_class_summary.tsv")),
                     sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
}

# ============================================================
# Enhanced family-level plots
# ============================================================
if (!is.na(family_file)) {
  fam <- resolve_schema(
    read_tsv(family_file), "family_level", family_file,
    list(
      "TE Family" = c("TE Family", "TE_Family", "Family"),
      "Coverage (bp)" = c("Coverage (bp)", "Coverage_bp", "Coverage"),
      "Copy Number" = c("Copy Number", "Copy_Number", "Copies")
    )
  )
  if (nrow(fam)) {
    fam$coverage_bp <- num_col(fam[["Coverage (bp)"]])
    fam$copy_number_num <- num_col(fam[["Copy Number"]])
    fam$family <- sub("#.*$", "", fam[["TE Family"]])
    fam$te_type <- ifelse(grepl("#", fam[["TE Family"]]), sub("^.*#", "", fam[["TE Family"]]), "Unknown")
    fam$class_base <- class_base(fam$te_type)
    fam <- fam[is.finite(fam$coverage_bp) & fam$coverage_bp > 0, , drop = FALSE]
    fam <- fam[order(-fam$coverage_bp), , drop = FALSE]
    top <- utils::head(fam, top_n)

    if (nrow(top)) {
      top$family <- factor(top$family, levels = rev(top$family))
      p <- ggplot2::ggplot(top, ggplot2::aes(x = family, y = coverage_bp / 1e6, fill = class_base)) +
        ggplot2::geom_col(width = 0.72) + ggplot2::coord_flip() +
        ggplot2::scale_fill_manual(values = class_palette, na.value = na_debug_color, drop = FALSE) +
        ggplot2::labs(
          title = paste("Top", top_n, "repeat families of", species_label),
          subtitle = "Families ranked by total Earl Grey genomic coverage",
          x = NULL, y = "Coverage (Mb)", fill = "Repeat subclass"
        ) + enhanced_theme(base_size = 10)
      manifest <- rbind(manifest, save_plot(p, "enhanced_top_te_families", 10.5, 8.5))
    }

    top_copy <- fam[is.finite(fam$copy_number_num) & fam$copy_number_num > 0, , drop = FALSE]
    top_copy <- top_copy[order(-top_copy$copy_number_num), , drop = FALSE]
    top_copy <- utils::head(top_copy, top_n)
    if (nrow(top_copy)) {
      top_copy$family <- factor(top_copy$family, levels = rev(top_copy$family))
      p <- ggplot2::ggplot(top_copy, ggplot2::aes(x = family, y = copy_number_num, fill = class_base)) +
        ggplot2::geom_col(width = 0.72) + ggplot2::coord_flip() +
        ggplot2::scale_fill_manual(values = class_palette, na.value = na_debug_color, drop = FALSE) +
        ggplot2::scale_y_continuous(labels = label_bp, expand = ggplot2::expansion(mult = c(0, 0.05))) +
        ggplot2::labs(
          title = paste("Top", top_n, "repeat families by copy number in", species_label),
          subtitle = "Families ranked by Earl Grey annotated copy count",
          x = NULL, y = "Copy number", fill = "Repeat subclass"
        ) + enhanced_theme(base_size = 10)
      manifest <- rbind(manifest, save_plot(p, "enhanced_top_te_families_by_copy_number", 10.5, 8.5))
    }

    scatter <- fam[is.finite(fam$copy_number_num) & fam$copy_number_num > 0, , drop = FALSE]
    if (nrow(scatter)) {
      p <- ggplot2::ggplot(scatter, ggplot2::aes(x = copy_number_num, y = coverage_bp / 1e6, colour = class_base)) +
        ggplot2::geom_point(alpha = 0.62, size = 2) +
        ggplot2::scale_x_log10(labels = label_bp) + ggplot2::scale_y_log10() +
        ggplot2::scale_colour_manual(values = class_palette, na.value = na_debug_color, drop = FALSE) +
        ggplot2::labs(
          title = paste("Repeat-family abundance of", species_label),
          subtitle = "Coverage integrates copy number and the lengths of annotated family fragments",
          x = "Copy number (log scale)", y = "Coverage (Mb, log scale)", colour = "Repeat subclass"
        ) + enhanced_theme()
      manifest <- rbind(manifest, save_plot(p, "enhanced_family_copy_vs_coverage", 8.8, 6.4))
    }
  }
}

# ============================================================
# Enhanced Earl Grey divergence landscapes
# ============================================================
if (nrow(landscape_gff) > 0) {
  landscape_bins <- bin_landscape(landscape_gff)
  class_bins <- stats::aggregate(bp ~ subclass + subclass_label + kimura_bin, data = landscape_bins, FUN = sum, na.rm = TRUE)
  class_bins$subclass_label <- factor(class_bins$subclass_label, levels = subclass_label_levels)

  total_landscape_bp <- sum(landscape_gff$width_bp, na.rm = TRUE)
  n_fragments <- nrow(landscape_gff)
  n_families_gff <- length(unique(extract_attr(landscape_gff$metadata, "NAME")))
  landscape_subtitle <- paste0(
    format(n_fragments, big.mark = ","), " divergence-scored fragments  ·  ",
    sprintf("%.1f Mb", total_landscape_bp / 1e6), " represented",
    if (n_families_gff > 0) paste0("  ·  ", format(n_families_gff, big.mark = ","), " families") else ""
  )

  # Peak and weighted K2P summaries by subclass.
  peak_rows <- do.call(rbind, lapply(split(class_bins, class_bins$subclass_label), function(d) {
    d[which.max(d$bp), c("subclass_label", "kimura_bin", "bp"), drop = FALSE]
  }))
  names(peak_rows)[names(peak_rows) == "kimura_bin"] <- "peak_k2p"
  names(peak_rows)[names(peak_rows) == "bp"] <- "peak_bin_bp"
  subclass_stats <- do.call(rbind, lapply(split(landscape_gff, landscape_gff$subclass_label), function(d) {
    data.frame(
      subclass = as.character(d$subclass[1]), subclass_label = as.character(d$subclass_label[1]),
      fragments = nrow(d), total_bp = sum(d$width_bp, na.rm = TRUE),
      mean_k2p = mean(d$kimura80, na.rm = TRUE), median_k2p = stats::median(d$kimura80, na.rm = TRUE),
      weighted_mean_k2p = weighted_mean_safe(d$kimura80, d$width_bp), stringsAsFactors = FALSE
    )
  }))
  subclass_stats <- merge(subclass_stats, peak_rows, by = "subclass_label", all.x = TRUE)
  utils::write.table(subclass_stats, file.path(outdir, paste0(species, "_divergence_subclass_summary.tsv")),
                     sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")

  # Family-level K2P summary reconstructed directly from the divergence GFF.
  family_name <- extract_attr(landscape_gff$metadata, "NAME")
  family_name[is.na(family_name) | !nzchar(family_name)] <- "Unidentified"
  gff_family <- landscape_gff
  gff_family$family <- sub("#.*$", "", family_name)
  family_stats <- do.call(rbind, lapply(split(gff_family, gff_family$family), function(d) {
    data.frame(
      family = d$family[1], subclass = as.character(d$subclass[1]),
      subclass_label = as.character(d$subclass_label[1]), superfamily = as.character(d$superfamily[1]),
      fragments = nrow(d), total_bp = sum(d$width_bp, na.rm = TRUE),
      mean_width = mean(d$width_bp, na.rm = TRUE), median_width = stats::median(d$width_bp, na.rm = TRUE),
      min_width = min(d$width_bp, na.rm = TRUE), max_width = max(d$width_bp, na.rm = TRUE),
      sd_width = stats::sd(d$width_bp, na.rm = TRUE), mean_k2p = mean(d$kimura80, na.rm = TRUE),
      median_k2p = stats::median(d$kimura80, na.rm = TRUE), weighted_mean_k2p = weighted_mean_safe(d$kimura80, d$width_bp),
      min_k2p = min(d$kimura80, na.rm = TRUE), max_k2p = max(d$kimura80, na.rm = TRUE),
      sd_k2p = stats::sd(d$kimura80, na.rm = TRUE), stringsAsFactors = FALSE
    )
  }))
  family_stats <- family_stats[order(family_stats$subclass_label, family_stats$superfamily, family_stats$family), , drop = FALSE]
  utils::write.table(family_stats, file.path(outdir, paste0(species, "_enhanced_family_divergence_summary.tsv")),
                     sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")

  # A. Main stacked subclass landscape.
  p <- ggplot2::ggplot(class_bins, ggplot2::aes(x = kimura_bin, y = bp, fill = subclass_label)) +
    ggplot2::geom_col(width = binwidth * 0.94, linewidth = 0) +
    kimura_scale() + landscape_y_scale(genome_size) +
    ggplot2::scale_fill_manual(values = class_palette, na.value = na_debug_color, drop = FALSE) +
    ggplot2::labs(
      title = paste("Repeat landscape of", species_label), subtitle = landscape_subtitle,
      fill = "Repeat subclass",
      caption = paste0("K2P values are explicitly binned at width ", format(binwidth), ". Lower divergence indicates less-diverged repeat copies.")
    ) + enhanced_theme() + ggplot2::theme(legend.position = "bottom") +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2, byrow = TRUE))
  manifest <- rbind(manifest, save_plot(p, "enhanced_classification_landscape", 12, 7.4))

  # B. Split-class landscape: free-y emphasizes temporal structure.
  p <- ggplot2::ggplot(class_bins, ggplot2::aes(x = kimura_bin, y = bp, fill = subclass_label)) +
    ggplot2::geom_col(width = binwidth * 0.94, linewidth = 0) +
    ggplot2::facet_wrap(~ subclass_label, ncol = 1, scales = "free_y") +
    kimura_scale() + landscape_y_scale(genome_size) +
    ggplot2::scale_fill_manual(values = class_palette, na.value = na_debug_color, drop = FALSE) +
    ggplot2::labs(
      title = paste("Subclass repeat landscapes of", species_label), subtitle = "Free y-scales emphasize divergence structure within each repeat subclass",
      fill = "Repeat subclass", caption = "Panel heights are independently scaled; compare peak positions, not absolute panel heights."
    ) + enhanced_theme(base_size = 9.5) + ggplot2::theme(legend.position = "none")
  manifest <- rbind(manifest, save_plot(p, "enhanced_split_class_landscape", 11.5, 10.5))

  # C. Fixed-y split landscape enables direct abundance comparison.
  p <- ggplot2::ggplot(class_bins, ggplot2::aes(x = kimura_bin, y = bp, fill = subclass_label)) +
    ggplot2::geom_col(width = binwidth * 0.94, linewidth = 0) +
    ggplot2::facet_wrap(~ subclass_label, ncol = 2, scales = "fixed") +
    kimura_scale() + landscape_y_scale(genome_size) +
    ggplot2::scale_fill_manual(values = class_palette, na.value = na_debug_color, drop = FALSE) +
    ggplot2::labs(
      title = paste("Absolute subclass landscapes of", species_label), subtitle = "Common y-scale permits direct comparison of repeat abundance among subclasses",
      fill = "Repeat subclass"
    ) + enhanced_theme(base_size = 9.5) + ggplot2::theme(legend.position = "none")
  manifest <- rbind(manifest, save_plot(p, "enhanced_split_class_landscape_fixed", 11.5, 8.4))

  # D. Superfamily landscape with top-N superfamilies retained per subclass.
  super_bins <- stats::aggregate(bp ~ subclass + subclass_label + superfamily + kimura_bin, data = landscape_bins, FUN = sum, na.rm = TRUE)
  sf_totals <- stats::aggregate(bp ~ subclass_label + superfamily, data = super_bins, FUN = sum, na.rm = TRUE)
  keep_keys <- do.call(rbind, lapply(split(sf_totals, sf_totals$subclass_label), function(d) {
    d <- d[order(-d$bp), , drop = FALSE]
    utils::head(d, top_superfamilies)[, c("subclass_label", "superfamily"), drop = FALSE]
  }))
  key <- paste(super_bins$subclass_label, super_bins$superfamily, sep = "\r")
  keep_key <- paste(keep_keys$subclass_label, keep_keys$superfamily, sep = "\r")
  super_bins$superfamily_plot <- as.character(super_bins$superfamily)
  drop_idx <- !key %in% keep_key
  super_bins$superfamily_plot[drop_idx] <- paste("Other", as.character(super_bins$subclass_label[drop_idx]))
  super_plot <- stats::aggregate(bp ~ subclass + subclass_label + superfamily_plot + kimura_bin, data = super_bins, FUN = sum, na.rm = TRUE)
  sf_palette <- complete_palette(super_plot$superfamily_plot, superfamily_palette, "TE superfamily")
  # Make aggregated Other-* entries subtle versions of their parent class.
  for (lab in unique(super_plot$subclass_label)) {
    other_lab <- paste("Other", lab)
    if (other_lab %in% super_plot$superfamily_plot) sf_palette[other_lab] <- "#D9D9D9"
  }

  p <- ggplot2::ggplot(super_plot, ggplot2::aes(x = kimura_bin, y = bp, fill = superfamily_plot)) +
    ggplot2::geom_col(width = binwidth * 0.94, linewidth = 0) +
    ggplot2::facet_grid(subclass_label ~ ., scales = "free_y", switch = "y") +
    kimura_scale() + landscape_y_scale(genome_size) +
    ggplot2::scale_fill_manual(values = sf_palette, na.value = na_debug_color, drop = FALSE) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 3, byrow = TRUE)) +
    ggplot2::labs(
      title = paste("Superfamily repeat landscape of", species_label),
      subtitle = paste("Top", top_superfamilies, "superfamilies per subclass; remaining groups combined as Other"),
      fill = "Superfamily", caption = "Free y-scales emphasize within-subclass activity profiles."
    ) + enhanced_theme(base_size = 8.8) +
    ggplot2::theme(strip.placement = "outside", legend.position = "bottom")
  manifest <- rbind(manifest, save_plot(p, "enhanced_superfamily_divergence_landscape", 12, 10))

  # E. Peak-K2P dot plot summarizing the dominant divergence mode of each class.
  peak_rows$subclass_label <- factor(peak_rows$subclass_label, levels = subclass_label_levels)
  p <- ggplot2::ggplot(peak_rows, ggplot2::aes(x = peak_k2p, y = subclass_label, colour = subclass_label, size = peak_bin_bp)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = peak_k2p, yend = subclass_label), colour = "grey82", linewidth = 0.8) +
    ggplot2::geom_point(alpha = 0.95) +
    ggplot2::scale_colour_manual(values = class_palette, drop = FALSE, na.value = na_debug_color) +
    ggplot2::scale_size_continuous(range = c(3, 9), labels = function(x) paste0(format_mb(x), " Mb")) +
    ggplot2::scale_x_continuous(limits = c(0, max_kimura), breaks = seq(0, max_kimura, 0.05), expand = ggplot2::expansion(mult = c(0.01, 0.03))) +
    ggplot2::labs(
      title = paste("Peak repeat divergence of", species_label),
      subtitle = "Peak K2P identifies the divergence bin containing the largest annotated sequence contribution",
      x = "Peak Kimura 2-parameter divergence (K2P)", y = NULL, colour = "Repeat subclass", size = "Peak-bin sequence"
    ) + enhanced_theme() + ggplot2::theme(legend.position = "bottom")
  manifest <- rbind(manifest, save_plot(p, "enhanced_peak_k2p_summary", 9, 5.8))
}

# ============================================================
# Family-level divergence summaries from Earl Grey TSV
# ============================================================
if (!is.na(div_file)) {
  div_raw <- read_tsv(div_file)
  if (nrow(div_raw)) {
    aliases <- list(
      "subclass" = c("subclass", "Subclass", "TE subclass"),
      "superfamily" = c("superfamily", "Superfamily", "TE superfamily"),
      "family" = c("family", "Family", "TE family"),
      "total_bp" = c("total_bp", "total bp", "Coverage (bp)"),
      "mean_width" = c("mean_width", "mean width", "Mean width"),
      "mean_div" = c("weighted_mean_k2p", "mean_div", "mean divergence", "Mean divergence", "mean_k2p"),
      "min_div" = c("min_div", "min divergence", "Min divergence", "min_k2p"),
      "max_div" = c("max_div", "max divergence", "Max divergence", "max_k2p"),
      "sd_div" = c("sd_div", "sd divergence", "SD divergence", "sd_k2p")
    )

    # Only the first five fields are essential; optional dispersion fields are
    # preserved when available without making older Earl Grey tables fail.
    essential_aliases <- aliases[c("subclass", "superfamily", "family", "total_bp", "mean_div")]
    div <- resolve_schema(div_raw, "divergence", div_file, essential_aliases)

    optional_map <- aliases[c("mean_width", "min_div", "max_div", "sd_div")]
    for (canonical in names(optional_map)) {
      if (!canonical %in% names(div)) {
        candidates <- optional_map[[canonical]]
        source <- candidates[candidates %in% names(div)][1]
        if (!is.na(source)) {
          names(div)[names(div) == source] <- canonical
        } else {
          div[[canonical]] <- NA_real_
        }
      }
    }

    numeric_fields <- c("total_bp", "mean_width", "mean_div", "min_div", "max_div", "sd_div")
    for (nm in numeric_fields) div[[nm]] <- num_col(div[[nm]])

    div$subclass_abbrev <- subclass_abbrev(div$subclass)
    div$class_base <- factor(subclass_label(div$subclass), levels = subclass_label_levels)
    div$superfamily <- clean_token(div$superfamily)
    div$family <- clean_token(div$family)

    div_valid <- div[
      is.finite(div$total_bp) & div$total_bp > 0 &
        is.finite(div$mean_div) & div$mean_div >= 0,
      , drop = FALSE
    ]
    div_valid$coverage_mb <- div_valid$total_bp / 1e6

    # Persist the cleaned family-level table used by all v5 divergence figures.
    utils::write.table(
      div_valid,
      file.path(outdir, paste0(species, "_family_divergence_clean.tsv")),
      sep = "\t", row.names = FALSE, quote = FALSE, na = "NA"
    )

    if (nrow(div_valid)) {
      # F1. Family-level abundance vs divergence.
      point_size <- ifelse(
        is.finite(div_valid$total_bp),
        pmax(1.2, pmin(6.2, 1.2 + 1.25 * log10(pmax(div_valid$total_bp, 1)) - 5)),
        1.8
      )
      p <- ggplot2::ggplot(
        div_valid,
        ggplot2::aes(x = mean_div, y = coverage_mb, colour = class_base)
      ) +
        ggplot2::geom_point(
          ggplot2::aes(size = total_bp),
          alpha = 0.72, stroke = 0.25
        ) +
        ggplot2::facet_wrap(~ class_base, scales = "free_y", ncol = 3, drop = TRUE) +
        ggplot2::scale_colour_manual(values = class_palette, na.value = na_debug_color, drop = FALSE) +
        ggplot2::scale_size_continuous(
          range = c(1.4, 6.5),
          labels = function(x) paste0(format_mb(x), " Mb")
        ) +
        ggplot2::scale_x_continuous(
          limits = c(0, max(max_kimura, min(0.8, max(div_valid$mean_div, na.rm = TRUE)))),
          expand = ggplot2::expansion(mult = c(0.01, 0.02))
        ) +
        ggplot2::labs(
          title = paste("Family-level Kimura divergence of", species_label),
          subtitle = "Each point is a repeat family; coverage and point size represent total annotated sequence",
          x = "Mean Kimura 2-parameter divergence (K2P)",
          y = "Family coverage (Mb)",
          colour = "Repeat subclass",
          size = "Family coverage",
          caption = "Rows lacking a finite mean divergence are excluded rather than treated as K2P = 0."
        ) +
        enhanced_theme(base_size = 9.8) +
        ggplot2::theme(legend.position = "right")
      manifest <- rbind(
        manifest,
        save_plot(p, "enhanced_family_divergence", 12.2, 8.4)
      )

      # F2. Coverage-weighted family-mean K2P distribution by subclass.
      div_binned <- div_valid[
        div_valid$mean_div <= max_kimura,
        , drop = FALSE
      ]
      if (nrow(div_binned)) {
        div_binned$kimura_bin <- bin_numeric(div_binned$mean_div)
        family_bins <- stats::aggregate(
          total_bp ~ class_base + kimura_bin,
          data = div_binned,
          FUN = sum,
          na.rm = TRUE
        )
        names(family_bins)[names(family_bins) == "total_bp"] <- "bp"

        utils::write.table(
          family_bins,
          file.path(outdir, paste0(species, "_family_mean_k2p_binned_summary.tsv")),
          sep = "\t", row.names = FALSE, quote = FALSE, na = "NA"
        )

        p <- ggplot2::ggplot(
          family_bins,
          ggplot2::aes(x = kimura_bin, y = bp, fill = class_base)
        ) +
          ggplot2::geom_col(width = binwidth * 0.94, linewidth = 0) +
          kimura_scale() +
          landscape_y_scale(genome_size) +
          ggplot2::scale_fill_manual(values = class_palette, na.value = na_debug_color, drop = FALSE) +
          ggplot2::labs(
            title = paste("Family-mean repeat landscape of", species_label),
            subtitle = "Family coverage summarized by each family's mean K2P; complementary to the annotation-level Earl Grey landscape",
            fill = "Repeat subclass",
            caption = paste0(
              "Family means are binned at K2P width ", format(binwidth),
              ". This is a family-level summary and must not be interpreted as a per-copy K2P landscape."
            )
          ) +
          enhanced_theme() +
          ggplot2::theme(legend.position = "bottom") +
          ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2, byrow = TRUE))
        manifest <- rbind(
          manifest,
          save_plot(p, "enhanced_family_mean_classification_landscape", 12, 7.3)
        )

        # F3. Faceted family-mean Kimura landscape, matching Earl Grey visual grammar.
        p <- ggplot2::ggplot(
          family_bins,
          ggplot2::aes(x = kimura_bin, y = bp, fill = class_base)
        ) +
          ggplot2::geom_col(width = binwidth * 0.94, linewidth = 0, alpha = 0.86) +
          ggplot2::facet_wrap(~ class_base, ncol = 3, scales = "free_y", drop = TRUE) +
          kimura_scale() +
          landscape_y_scale(genome_size) +
          ggplot2::scale_fill_manual(values = class_palette, na.value = na_debug_color, drop = FALSE) +
          ggplot2::labs(
            title = paste("Family-mean Kimura repeat landscape of", species_label),
            subtitle = "Summed family coverage within mean-K2P bins; free y-scales emphasize subclass-specific divergence structure",
            fill = "Repeat subclass",
            caption = "These family-level summaries complement, but do not replace, annotation-level K2P landscapes from the divergence GFF."
          ) +
          enhanced_theme(base_size = 9.5) +
          ggplot2::theme(legend.position = "right")
        manifest <- rbind(
          manifest,
          save_plot(p, "enhanced_family_mean_kimura_landscape", 12.4, 8.5)
        )
      }

      # F4. Family mean-divergence distribution by subclass.
      dist_df <- div_valid[div_valid$mean_div <= max(0.8, max_kimura), , drop = FALSE]
      if (nrow(dist_df)) {
        p <- ggplot2::ggplot(
          dist_df,
          ggplot2::aes(x = class_base, y = mean_div, fill = class_base)
        ) +
          ggplot2::geom_violin(
            trim = TRUE, scale = "width", alpha = 0.72, linewidth = 0.25
          ) +
          ggplot2::geom_boxplot(
            width = 0.16, outlier.shape = NA, alpha = 0.9, linewidth = 0.3
          ) +
          ggplot2::coord_flip() +
          ggplot2::scale_fill_manual(values = class_palette, na.value = na_debug_color, drop = FALSE) +
          ggplot2::labs(
            title = paste("Repeat-family divergence distributions of", species_label),
            subtitle = "Distribution of family mean K2P values; family abundance is summarized separately",
            x = NULL, y = "Mean Kimura 2-parameter divergence (K2P)",
            fill = "Repeat subclass"
          ) +
          enhanced_theme() +
          ggplot2::theme(legend.position = "none")
        manifest <- rbind(
          manifest,
          save_plot(p, "enhanced_family_divergence_distribution", 9.3, 6.3)
        )
      }

      # F5. Top families: genomic coverage plus mean K2P and observed range.
      top_div <- div_valid[order(-div_valid$total_bp), , drop = FALSE]
      top_div <- utils::head(top_div, top_n)
      if (nrow(top_div)) {
        top_div$family_plot <- factor(top_div$family, levels = rev(top_div$family))
        p <- ggplot2::ggplot(
          top_div,
          ggplot2::aes(x = mean_div, y = family_plot, colour = class_base)
        ) +
          ggplot2::geom_segment(
            ggplot2::aes(
              x = ifelse(is.finite(min_div), min_div, mean_div),
              xend = ifelse(is.finite(max_div), max_div, mean_div),
              yend = family_plot
            ),
            colour = "grey72", linewidth = 0.65
          ) +
          ggplot2::geom_point(ggplot2::aes(size = total_bp), alpha = 0.92) +
          ggplot2::scale_colour_manual(values = class_palette, na.value = na_debug_color, drop = FALSE) +
          ggplot2::scale_size_continuous(
            range = c(2.2, 7.5),
            labels = function(x) paste0(format_mb(x), " Mb")
          ) +
          ggplot2::labs(
            title = paste("Top", top_n, "repeat families: coverage and divergence"),
            subtitle = paste("Families ranked by total genomic coverage in", species_label),
            x = "Mean K2P (segment = observed min–max K2P)",
            y = NULL,
            colour = "Repeat subclass",
            size = "Family coverage",
            caption = "Min–max segments are shown when Earl Grey reports finite family divergence bounds."
          ) +
          enhanced_theme(base_size = 9.2) +
          ggplot2::theme(legend.position = "right")
        manifest <- rbind(
          manifest,
          save_plot(p, "enhanced_top_family_divergence", 11.7, 8.7)
        )
      }

      # F6. Coverage-weighted subclass summary from the family-level table.
      family_subclass_summary <- do.call(
        rbind,
        lapply(split(div_valid, div_valid$class_base), function(d) {
          data.frame(
            subclass = as.character(d$class_base[1]),
            families = nrow(d),
            total_bp = sum(d$total_bp, na.rm = TRUE),
            mean_of_family_means = mean(d$mean_div, na.rm = TRUE),
            median_family_mean = stats::median(d$mean_div, na.rm = TRUE),
            coverage_weighted_mean_k2p = weighted_mean_safe(d$mean_div, d$total_bp),
            coverage_weighted_median_k2p = weighted_median_safe(d$mean_div, d$total_bp),
            stringsAsFactors = FALSE
          )
        })
      )
      family_subclass_summary <- family_subclass_summary[
        order(match(family_subclass_summary$subclass, subclass_label_levels)),
        , drop = FALSE
      ]
      utils::write.table(
        family_subclass_summary,
        file.path(outdir, paste0(species, "_family_divergence_subclass_summary.tsv")),
        sep = "\t", row.names = FALSE, quote = FALSE, na = "NA"
      )

      # F7. Superfamily × K2P heatmap, restricted to the most abundant superfamilies.
      heat <- div_valid[
        div_valid$mean_div <= max_kimura &
          !is.na(div_valid$superfamily) & nzchar(div_valid$superfamily),
        , drop = FALSE
      ]
      if (nrow(heat)) {
        sf_total <- stats::aggregate(total_bp ~ superfamily, data = heat, FUN = sum, na.rm = TRUE)
        sf_total <- sf_total[order(-sf_total$total_bp), , drop = FALSE]
        keep_sf <- utils::head(sf_total$superfamily, max(top_superfamilies, 8))
        heat <- heat[heat$superfamily %in% keep_sf, , drop = FALSE]
        heat$kimura_bin <- bin_numeric(heat$mean_div)
        heat_sum <- stats::aggregate(
          total_bp ~ superfamily + kimura_bin,
          data = heat, FUN = sum, na.rm = TRUE
        )
        heat_order <- sf_total$superfamily[sf_total$superfamily %in% keep_sf]
        heat_sum$superfamily <- factor(heat_sum$superfamily, levels = rev(heat_order))

        p <- ggplot2::ggplot(
          heat_sum,
          ggplot2::aes(x = kimura_bin, y = superfamily, fill = total_bp)
        ) +
          ggplot2::geom_tile(colour = "white", linewidth = 0.15) +
          kimura_scale() +
          ggplot2::scale_fill_viridis_c(
            option = "C", trans = "sqrt",
            labels = function(x) paste0(format_mb(x), " Mb")
          ) +
          ggplot2::labs(
            title = paste("Superfamily abundance-divergence matrix of", species_label),
            subtitle = paste("Top", length(keep_sf), "superfamilies by total family-level genomic coverage"),
            y = NULL, fill = "Coverage",
            caption = "Cell intensity is summed family coverage within each mean-K2P bin."
          ) +
          enhanced_theme(base_size = 9.3)
        manifest <- rbind(
          manifest,
          save_plot(p, "enhanced_superfamily_divergence_heatmap", 11.2, 7.2)
        )
      }
    }
  }
}

if (nrow(manifest) == 0) stop("No Earl Grey plots were created; required input files were missing, empty, or unreadable")

utils::write.table(manifest, file.path(outdir, "enhanced_plots_manifest.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(schema_audit, file.path(outdir, "enhanced_plots_schema.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

input_manifest <- data.frame(
  key = c(
    "species", "earlgrey_dir", "high_level_count", "summary_source", "summary_source_path",
    "family_level_count", "divergence_summary", "gff_with_kimura",
    "high_level_override", "family_level_override", "divergence_summary_override", "gff_override",
    "outdir", "dpi",
    "top_families", "top_superfamilies", "binwidth", "max_kimura", "y_scale", "axis_flip",
    "composition_tolerance", "debug_palette", "data_table_available", "scales_available", "patchwork_available"
  ),
  value = c(
    species, earlgrey_dir, high_file, summary_source, summary_source_path, family_file, div_file, gff_file,
    high_level_arg, family_level_arg, divergence_summary_arg, gff_arg,
    outdir, as.character(dpi), as.character(top_n), as.character(top_superfamilies), as.character(binwidth),
    as.character(max_kimura), y_scale, as.character(axis_flip), as.character(composition_tolerance),
    as.character(debug_palette), as.character(has_data_table), as.character(has_scales), as.character(has_patchwork)
  ), stringsAsFactors = FALSE
)
utils::write.table(input_manifest, file.path(outdir, "enhanced_plots_inputs.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

message("Created ", nrow(manifest), " plot files (", length(unique(manifest$plot)), " plot types) in: ", outdir)
message("Summary source: ", summary_source, if (!is.na(summary_source_path)) paste0(" [", summary_source_path, "]") else "")
message("Input manifest: ", file.path(outdir, "enhanced_plots_inputs.tsv"))
message("Plot manifest: ", file.path(outdir, "enhanced_plots_manifest.tsv"))
RSCRIPT

  chmod +x "$script_path"
}

# ============================================================
# VALIDATION
# ============================================================

validate_toggles

validate_runner ENV_RUNNER "$ENV_RUNNER"

for version_setting in TERRIER_REQUIRED_VERSION TETRIMMER_REQUIRED_VERSION EARLGREY_REQUIRED_VERSION EARLGREY_REQUIRED_REPEATMASKER EARLGREY_REQUIRED_REPEATMODELER EARLGREY_REQUIRED_DFAM; do
  version_value="${!version_setting}"
  [[ "$version_value" =~ ^[0-9]+(\.[0-9]+)+$ ]] || error "${version_setting} must be a dotted numeric version: ${version_value}"
done

validate_yes_no USE_LTRSTRUCT "$USE_LTRSTRUCT"
check_toggle RM_QUICK "$RM_QUICK"
check_toggle TERRIER_ENABLE "$TERRIER_ENABLE"
check_toggle TERRIER_UNKNOWN_ONLY "$TERRIER_UNKNOWN_ONLY"
check_toggle TETRIMMER_CONTINUE "$TETRIMMER_CONTINUE"
check_toggle TETRIMMER_ARCHIVE_INCOMPLETE "$TETRIMMER_ARCHIVE_INCOMPLETE"
check_toggle REPEAT_USE_REPBASE "$REPEAT_USE_REPBASE"
check_toggle TETRIMMER_DEDUP "$TETRIMMER_DEDUP"
check_toggle TETRIMMER_HMM "$TETRIMMER_HMM"
check_toggle TETRIMMER_GENOME_ANNO "$TETRIMMER_GENOME_ANNO"
check_toggle TETRIMMER_DEBUG "$TETRIMMER_DEBUG"
check_toggle TETRIMMER_COMPRESS_OUTPUT "$TETRIMMER_COMPRESS_OUTPUT"
check_toggle TETRIMMER_EXPORT_COVERAGE "$TETRIMMER_EXPORT_COVERAGE"
check_toggle TETRIMMER_SKIP_REST_SEQ "$TETRIMMER_SKIP_REST_SEQ"
check_toggle TETRIMMER_EXPERT_OVERRIDES "$TETRIMMER_EXPERT_OVERRIDES"
if [[ -n "$RM_SRAND" ]] && ! [[ "$RM_SRAND" =~ ^[0-9]+$ ]]; then
  echo "[ERROR] Invalid RM_SRAND=$RM_SRAND. Use an integer." >&2
  exit 1
fi
if [[ -n "$RM_GENOME_SAMPLE_SIZE_MAX" ]] && { ! [[ "$RM_GENOME_SAMPLE_SIZE_MAX" =~ ^[0-9]+$ ]] || [[ "$RM_GENOME_SAMPLE_SIZE_MAX" -lt 1 ]]; }; then
  echo "[ERROR] Invalid RM_GENOME_SAMPLE_SIZE_MAX=$RM_GENOME_SAMPLE_SIZE_MAX. Use a positive integer." >&2
  exit 1
fi
if ! [[ "$RM_ADDL_ROUNDS" =~ ^[0-9]+$ ]]; then
  echo "[ERROR] Invalid RM_ADDL_ROUNDS=$RM_ADDL_ROUNDS. Use a non-negative integer." >&2
  exit 1
fi
if [[ ! "$TETRIMMER_THREAD_CAP" =~ ^[0-9]+$ || "$TETRIMMER_THREAD_CAP" -lt 1 ]]; then
  echo "[ERROR] Invalid TETRIMMER_THREAD_CAP=$TETRIMMER_THREAD_CAP. Use a positive integer no greater than 16." >&2
  exit 1
fi
if [[ "$TETRIMMER_THREAD_CAP" -gt 16 ]]; then
  echo "[WARN] TETRIMMER_THREAD_CAP=$TETRIMMER_THREAD_CAP exceeds the hard safety limit of 16; using 16."
  TETRIMMER_THREAD_CAP=16
fi
if ! [[ "$TETRIMMER_MAX_THREAD_TIME" =~ ^[0-9]+$ ]] || [[ "$TETRIMMER_MAX_THREAD_TIME" -lt 1 ]]; then
  echo "[ERROR] Invalid TETRIMMER_MAX_THREAD_TIME=$TETRIMMER_MAX_THREAD_TIME. Use a positive number of seconds." >&2
  exit 1
fi
if ! [[ "$THREADS" =~ ^[0-9]+$ ]] || [[ "$THREADS" -lt 1 ]]; then
  echo "[ERROR] Invalid REPEATS_THREADS=$THREADS. Use a positive integer." >&2
  exit 1
fi
awk -v x="$FINAL_CDHIT_ID" 'BEGIN {exit !(x > 0 && x <= 1)}' ||
  error "FINAL_CDHIT_ID must be in (0,1]: $FINAL_CDHIT_ID"
[[ "$FINAL_CDHIT_WORD" =~ ^[0-9]+$ && "$FINAL_CDHIT_WORD" -ge 4 && "$FINAL_CDHIT_WORD" -le 11 ]] ||
  error "FINAL_CDHIT_WORD must be an integer in [4,11]: $FINAL_CDHIT_WORD"
awk -v x="$FINAL_CDHIT_AS" 'BEGIN {exit !(x >= 0 && x <= 1)}' ||
  error "FINAL_CDHIT_AS must be in [0,1]: $FINAL_CDHIT_AS"
awk -v x="$FINAL_CDHIT_AL" 'BEGIN {exit !(x >= 0 && x <= 1)}' ||
  error "FINAL_CDHIT_AL must be in [0,1]: $FINAL_CDHIT_AL"
[[ "$FINAL_CDHIT_MEMORY_MB" =~ ^[0-9]+$ && "$FINAL_CDHIT_MEMORY_MB" -ge 1 ]] ||
  error "FINAL_CDHIT_MEMORY_MB must be a positive integer: $FINAL_CDHIT_MEMORY_MB"
if is_true "$TERRIER_ENABLE"; then
  [[ "$TERRIER_THREADS" =~ ^[0-9]+$ && "$TERRIER_THREADS" -ge 1 ]] || error "TERRIER_THREADS must be a positive integer: $TERRIER_THREADS"
  [[ "$TERRIER_MIN_LENGTH" =~ ^[0-9]+$ && "$TERRIER_MIN_LENGTH" -ge 1 ]] || error "TERRIER_MIN_LENGTH must be a positive integer: $TERRIER_MIN_LENGTH"
  [[ "$TERRIER_MAX_LENGTH" =~ ^[0-9]+$ && "$TERRIER_MAX_LENGTH" -ge "$TERRIER_MIN_LENGTH" ]] || error "TERRIER_MAX_LENGTH must be >= TERRIER_MIN_LENGTH"
  [[ "$TERRIER_BATCH_SIZE" =~ ^[0-9]+$ && "$TERRIER_BATCH_SIZE" -ge 1 ]] || error "TERRIER_BATCH_SIZE must be a positive integer"
  awk -v x="$TERRIER_THRESHOLD" 'BEGIN{exit !(x >= 0 && x <= 1)}' || error "TERRIER_THRESHOLD must be in [0,1]: $TERRIER_THRESHOLD"
  awk -v x="$TERRIER_IMAGE_THRESHOLD" 'BEGIN{exit !(x >= 0 && x <= 1)}' || error "TERRIER_IMAGE_THRESHOLD must be in [0,1]: $TERRIER_IMAGE_THRESHOLD"
  case "$TERRIER_IMAGE_FORMAT" in png|pdf|svg|jpg|jpeg) ;; *) error "Unsupported TERRIER_IMAGE_FORMAT=$TERRIER_IMAGE_FORMAT" ;; esac
fi
if [[ "$THREADS" -lt 4 ]]; then
  echo "[ERROR] Earl Grey requires at least 4 threads because it derives RepeatMasker -pa as threads/4." >&2
  exit 1
fi
if [[ -z "$SPECIES" ]]; then
  echo "[ERROR] REPEATS_SPECIES resolved to an empty value." >&2
  exit 1
fi
if [[ "$SPECIES" =~ [[:space:]] ]]; then
  echo "[ERROR] REPEATS_SPECIES cannot contain whitespace: $SPECIES" >&2
  exit 1
fi

validate_yes_no EARLGREY_REMOVE_SHORT "$EARLGREY_REMOVE_SHORT"
validate_yes_no EARLGREY_SOFTMASK "$EARLGREY_SOFTMASK"
validate_yes_no EARLGREY_HELIANO "$EARLGREY_HELIANO"
validate_yes_no EARLGREY_VALIDATE_RESUME "$EARLGREY_VALIDATE_RESUME"
check_toggle EARLGREY_ENHANCED_PLOTS "$EARLGREY_ENHANCED_PLOTS"
check_toggle EARLGREY_ARCHIVE_INCOMPLETE "$EARLGREY_ARCHIVE_INCOMPLETE"
case "$EARLGREY_INCOMPLETE_POLICY" in
  resume|archive|error) ;;
  *)
    echo "[ERROR] Invalid EARLGREY_INCOMPLETE_POLICY=$EARLGREY_INCOMPLETE_POLICY. Use resume, archive, or error." >&2
    exit 1
    ;;
esac

if ! [[ "$EARLGREY_PLOT_DPI" =~ ^[0-9]+$ ]] || [[ "$EARLGREY_PLOT_DPI" -lt 1 ]]; then
  echo "[ERROR] Invalid EARLGREY_PLOT_DPI=$EARLGREY_PLOT_DPI. Use a positive integer." >&2
  exit 1
fi

if ! [[ "$EARLGREY_PLOT_TOP_FAMILIES" =~ ^[0-9]+$ ]] || [[ "$EARLGREY_PLOT_TOP_FAMILIES" -lt 1 ]]; then
  echo "[ERROR] Invalid EARLGREY_PLOT_TOP_FAMILIES=$EARLGREY_PLOT_TOP_FAMILIES. Use a positive integer." >&2
  exit 1
fi

for numeric_setting in EARLGREY_PLOT_BINWIDTH EARLGREY_PLOT_MAX_KIMURA EARLGREY_COMPOSITION_TOLERANCE; do
  numeric_value="${!numeric_setting}"
  if ! [[ "$numeric_value" =~ ^[0-9]*\.?[0-9]+$ ]]; then
    echo "[ERROR] ${numeric_setting} must be numeric: ${numeric_value}" >&2
    exit 1
  fi
done

awk -v x="$EARLGREY_PLOT_BINWIDTH" 'BEGIN {exit !(x > 0)}' || error "EARLGREY_PLOT_BINWIDTH must be greater than zero"
awk -v x="$EARLGREY_PLOT_MAX_KIMURA" 'BEGIN {exit !(x > 0)}' || error "EARLGREY_PLOT_MAX_KIMURA must be greater than zero"
awk -v x="$EARLGREY_COMPOSITION_TOLERANCE" 'BEGIN {exit !(x >= 0 && x < 1)}' || error "EARLGREY_COMPOSITION_TOLERANCE must be in [0,1)"

if is_true "$TETRIMMER_EXPERT_OVERRIDES"; then
  for x in TETRIMMER_CONS_THR TETRIMMER_EXT_THR TETRIMMER_GAP_THR TETRIMMER_GAP_NUL_THR TETRIMMER_CROP_END_DIV_THR TETRIMMER_CROP_END_GAP_THR; do
    v="${!x}"
    awk -v z="$v" 'BEGIN {exit !(z >= 0 && z <= 1)}' || error "$x must be in [0,1]: $v"
  done
  for x in TETRIMMER_MINI_ORF TETRIMMER_MAX_MSA_LINES TETRIMMER_TOP_MSA_LINES TETRIMMER_MIN_SEQ_NUM TETRIMMER_MIN_BLAST_LEN TETRIMMER_MAX_CLUSTER_NUM TETRIMMER_EXT_CHECK_WIN TETRIMMER_EXT_STEP TETRIMMER_MAX_EXT TETRIMMER_CROP_END_DIV_WIN TETRIMMER_CROP_END_GAP_WIN TETRIMMER_POLY_LEN; do
    v="${!x}"
    [[ "$v" =~ ^[0-9]+$ && "$v" -ge 1 ]] || error "$x must be a positive integer: $v"
  done
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
REPEATS_LIB_DIR="$(readlink -m "$REPEATS_LIB_DIR")"
DFAM_FAMDB_DIR="$(readlink -m "$DFAM_FAMDB_DIR")"
REPBASE_DIR="$(readlink -m "$REPBASE_DIR")"

migrate_earlgrey_v730_layout

DIR_GENOME="${PROJECT_ROOT}/00.genome"
DIR_KNOWN="${PROJECT_ROOT}/01.known_libraries"
DIR_RM="${PROJECT_ROOT}/02.repeatmodeler"
DIR_TERRIER="${PROJECT_ROOT}/02.repeatmodeler/terrier"
DIR_TETRIMMER="${PROJECT_ROOT}/03.tetrimmer"
TETRIMMER_PFAM_DIR="${TETRIMMER_PFAM_DIR:-${PROJECT_ROOT}/resources/pfam}"
TETRIMMER_RESULT_DIR=""
DIR_LIB="${PROJECT_ROOT}/04.final_library"
DIR_EARLGREY="${PROJECT_ROOT}/05.earlgrey_annotation_only"
DIR_FINAL="${PROJECT_ROOT}/06.final_outputs"
DIR_EARLGREY_ENHANCED_PLOTS="${DIR_FINAL}/earlgrey/enhanced_plots"
DIR_NON_REPEATS="${PROJECT_ROOT}/07.non_repetitive_regions"
DIR_LOGS="${PROJECT_ROOT}/logs"
EARLGREY_FAMDB_RUNTIME_DIR="${DIR_LOGS}/earlgrey_famdb_runtime"
EARLGREY_FAMDB_WRAPPER="${EARLGREY_FAMDB_RUNTIME_DIR}/famdb.py"
EARLGREY_FAMDB_REAL_PY=""
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
  "$DIR_TERRIER" \
  "$DIR_TETRIMMER" \
  "$DIR_LIB" \
  "$DIR_EARLGREY" \
  "$DIR_FINAL" \
  "$DIR_EARLGREY_ENHANCED_PLOTS" \
  "$DIR_NON_REPEATS" \
  "$DIR_LOGS" \
  "$REPEATS_CHECKPOINT_DIR"

mkdir -p "$(dirname "$TETRIMMER_PFAM_DIR")"
migrate_legacy_tetrimmer_pfam

GENOME_BASENAME="$(basename "$GENOME")"
GENOME_STEM="${GENOME_BASENAME%%.*}"
GENOME_LOCAL="${DIR_GENOME}/${GENOME_BASENAME}"

if is_true "$REPEATMODELER_SANITIZE_FASTA_HEADERS"; then
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

DFAM_TAXON_TAG="$(printf '%s' "$REPEAT_TAXON" | sed 's/[[:space:]]\+/_/g; s/[^A-Za-z0-9._-]/_/g')"
[[ -n "$DFAM_TAXON_TAG" ]] || fatal "Could not create a filesystem-safe repeat taxon tag from: $REPEAT_TAXON"
DFAM_FASTA="${DIR_KNOWN}/Dfam4.taxon_${DFAM_TAXON_TAG}.fa"
REPBASE_FASTA="${DIR_KNOWN}/RepBase.from_Libraries.fa"
KNOWN_RAW="${DIR_KNOWN}/known_repeats.Dfam4_RepBase.raw.fa"
KNOWN_NR="${DIR_KNOWN}/known_repeats.Dfam4_RepBase.nr.fa"

RM_LIBRARY="${DIR_RM}/${GENOME_STEM}-families.fa"

FINAL_RAW="${DIR_LIB}/${SPECIES}.TEtrimmer_Dfam4_RepBase.raw.fa"
FINAL_NR="${DIR_LIB}/${SPECIES}.TEtrimmer_Dfam4_RepBase.nr.fa"

# Detect outputs from previous runs before deciding which steps to execute.
auto_detect_outputs || true

if ! repeat_checkpoint_ready repeatmodeler || ! repeat_checkpoint_ready tetrimmer || ! repeat_checkpoint_ready earlgrey || known_library_needs_build; then
  prepare_dfam_famdb_dir
fi

if known_library_needs_build && is_true "$REPEAT_USE_REPBASE"; then
  resolve_repbase_embl_input
fi

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

  if ! repeat_checkpoint_ready repeatmodeler || known_library_needs_build; then
    require_dir "$DFAM_FAMDB_DIR" "Dfam FamDB directory"
  fi

  if known_library_needs_build && is_true "$REPEAT_USE_REPBASE"; then
    if [[ -d "$REPBASE_EMBL" ]]; then
      require_dir "$REPBASE_EMBL" "RepBase EMBL input"
    else
      require_file "$REPBASE_EMBL" "RepBase EMBL input"
    fi
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
    --env "REQUIRED_REPEATMASKER=$EARLGREY_REQUIRED_REPEATMASKER" \
    --env "REQUIRED_REPEATMODELER=$EARLGREY_REQUIRED_REPEATMODELER" \
    "$sif_repeatmasker" \
    bash -lc '
      set -Eeuo pipefail

      missing=""

      for tool in \
        perl \
        RepeatMasker \
        RepeatModeler \
        RepeatClassifier \
        BuildDatabase \
        cd-hit-est \
        gunzip
      do
        command -v "$tool" >/dev/null 2>&1 || \
          missing="$missing $tool"
      done

      command -v makeblastdb >/dev/null 2>&1 || \
        [[ -x /opt/rmblast/bin/makeblastdb ]] || \
        missing="$missing makeblastdb"

      [[ -x /opt/RepeatMasker/util/buildRMLibFromEMBL.pl ]] || \
        missing="$missing buildRMLibFromEMBL.pl"

      [[ -x /opt/FamDB/famdb.py ]] || \
        missing="$missing /opt/FamDB/famdb.py"

      [[ -s /opt/FamDB/famdb.conf ]] || \
        missing="$missing /opt/FamDB/famdb.conf"

      if [[ -n "$missing" ]]; then
        echo "[ERROR] Missing tools or configuration:$missing" >&2
        exit 1
      fi

      grep -Eq \
        "^[[:space:]]*FAMDB_DATA_DIR[[:space:]]*=[[:space:]]*/opt/RepeatMasker/Libraries/famdb[[:space:]]*$" \
        /opt/FamDB/famdb.conf || {
        echo "[ERROR] Unexpected FamDB database path" >&2
        cat /opt/FamDB/famdb.conf >&2
        exit 1
      }

      repeatmasker_version="$(RepeatMasker -v 2>&1)"
      repeatmodeler_version="$(RepeatModeler -version 2>&1)"
      printf "%s\n%s\n" "$repeatmasker_version" "$repeatmodeler_version"
      grep -Fq "$REQUIRED_REPEATMASKER" <<< "$repeatmasker_version" || {
        echo "[ERROR] Required RepeatMasker $REQUIRED_REPEATMASKER was not confirmed" >&2
        exit 1
      }
      grep -Fq "$REQUIRED_REPEATMODELER" <<< "$repeatmodeler_version" || {
        echo "[ERROR] Required RepeatModeler $REQUIRED_REPEATMODELER was not confirmed" >&2
        exit 1
      }
    '

  if known_library_needs_build || ! repeat_checkpoint_ready repeatmodeler; then
    echo "[INFO] Auditing external Dfam database"

    # shellcheck disable=SC2016
    singularity exec \
      --no-home \
      -B "$PROJECT_ROOT":"$PROJECT_ROOT" \
      -B "$DFAM_FAMDB_DIR":"$FAMDB_CONTAINER_DATA_DIR":ro \
      --env "FAMDB_DIR=$FAMDB_CONTAINER_DIR" \
      --env "FAMDB_DB_DIR=$FAMDB_CONTAINER_DATA_DIR" \
      --env "DFAM_RELEASE_TAG=$DFAM_RELEASE_TAG" \
      --env "DFAM_TAXON=$DFAM_TAXON" \
      --env "DFAM_FLAGS=$DFAM_FLAGS" \
      --env "DFAM_REQUIRED_PARTITIONS=$DFAM_REQUIRED_PARTITIONS" \
      --env "REQUIRED_DFAM_VERSION=$EARLGREY_REQUIRED_DFAM" \
      "$sif_repeatmasker" \
      bash -lc '
        set -Eeuo pipefail

        FAMDB_DIR="${FAMDB_DIR:?}"
        DFAM_RELEASE_TAG="${DFAM_RELEASE_TAG:?}"
        DFAM_TAXON="${DFAM_TAXON:?}"
        DFAM_FLAGS="${DFAM_FLAGS:-}"
        DFAM_REQUIRED_PARTITIONS="${DFAM_REQUIRED_PARTITIONS:?}"
        FAMDB_DB_DIR="${FAMDB_DB_DIR:?}"
        REQUIRED_DFAM_VERSION="${REQUIRED_DFAM_VERSION:?}"

        echo "FAMDB_DIR=${FAMDB_DIR}"
        echo "FAMDB_DB_DIR=${FAMDB_DB_DIR}"
        echo "DFAM_FLAGS=${DFAM_FLAGS}"
        echo "DFAM_TAXON=${DFAM_TAXON}"

        echo "[CHECK] Running: $FAMDB_DIR/famdb.py -i $FAMDB_DB_DIR info"
        # A conditional assignment preserves the command status under strict mode.
        if famdb_output="$("$FAMDB_DIR/famdb.py" -i "$FAMDB_DB_DIR" info 2>&1)"; then
          famdb_rc=0
        else
          famdb_rc=$?
        fi
        printf "%s\n" "$famdb_output"
        echo "[INFO] FamDB info exit status: $famdb_rc"
        if [[ "$famdb_rc" -ne 0 ]]; then
          echo "[ERROR] FamDB info failed with exit code $famdb_rc" >&2
          exit "$famdb_rc"
        fi

        grep -Eq "^Database[[:space:]]*:[[:space:]]*[^[:space:]]+" <<< "$famdb_output" || {
          echo "[ERROR] FamDB did not report a database name" >&2
          exit 1
        }

        grep -Eq "^Version[[:space:]]*:[[:space:]]*[^[:space:]]+" \
          <<< "$famdb_output" || {
          echo "[ERROR] FamDB did not report a version" >&2
          exit 1
        }

        # Compare the version literally; dots must not act as regex wildcards.
        observed_dfam_version="$(sed -nE "s/^Version[[:space:]]*:[[:space:]]*([^[:space:]]+).*/\1/p" <<< "$famdb_output")"
        [[ "$observed_dfam_version" == "$REQUIRED_DFAM_VERSION" ]] || {
          echo "[ERROR] Required Dfam version ${REQUIRED_DFAM_VERSION} was not confirmed (reported: ${observed_dfam_version})" >&2
          exit 1
        }

        echo "[INFO] FamDB partition report for: $DFAM_TAXON"

        if check_output="$("$FAMDB_DIR/famdb.py" -i "$FAMDB_DB_DIR" check "$DFAM_TAXON" 2>&1)"; then
          check_rc=0
        else
          check_rc=$?
        fi

        printf "%s\n" "$check_output"
        echo "[INFO] FamDB taxon check exit status: $check_rc"

        if [[ "$check_rc" -ne 0 ]]; then
          echo "[WARN] FamDB reports missing optional partitions."
          echo "[WARN] The configured 16-partition Dfam set is present; FamDB reports additional optional partitions."
        fi

        echo "[OK] Dfam database passed mode-aware validation."
      '
  fi

  {
    echo -e "Component\tPath"
    echo -e "RepeatMasker_SIF\t${sif_repeatmasker}"
    echo -e "Dfam_FamDB_host\t${DFAM_FAMDB_DIR}"
    echo -e "FamDB_software_container\t${FAMDB_CONTAINER_DIR}"
    echo -e "FamDB_config_container\t${FAMDB_CONTAINER_DIR}/famdb.conf"
    echo -e "FamDB_database_container\t${FAMDB_CONTAINER_DATA_DIR}"
    echo -e "RepeatModeler_container\t${REPEATMODELER_CONTAINER_DIR}"
    echo -e "Terrier_environment\t${ENV_TERRIER}"
    echo -e "Terrier_required_version\t${TERRIER_REQUIRED_VERSION}"
    echo -e "TEtrimmer_required_version\t${TETRIMMER_REQUIRED_VERSION}"
    echo -e "EarlGrey_required_version\t${EARLGREY_REQUIRED_VERSION}"
    echo -e "RepeatMasker_required_version\t${EARLGREY_REQUIRED_REPEATMASKER}"
    echo -e "RepeatModeler_required_version\t${EARLGREY_REQUIRED_REPEATMODELER}"
    echo -e "Dfam_required_version\t${EARLGREY_REQUIRED_DFAM}"
  } > "${DIR_LOGS}/repeat_tool_paths.tsv"

  if is_true "$TERRIER_ENABLE" && ! repeat_checkpoint_ready terrier; then
    echo "[INFO] Auditing Terrier environment"
    env_run "$ENV_TERRIER" env REQUIRED_TERRIER_VERSION="$TERRIER_REQUIRED_VERSION" TERRIER_UNKNOWN_ONLY_REQUIRED="$TERRIER_UNKNOWN_ONLY" bash -lc '
      set -Eeuo pipefail
      command -v terrier >/dev/null 2>&1 || { echo "[ERROR] terrier not found" >&2; exit 1; }
      help_output="$(terrier --help 2>&1)"
      for required_option in --input --output-csv --output-tips-csv --output-fasta --image-dir --image-format --image-threshold --threshold --min-length --max-length --batch-size --num-threads --unknown-only; do
        grep -Fq -- "$required_option" <<< "$help_output" || {
          echo "[ERROR] Terrier CLI is missing required option: $required_option" >&2
          exit 1
        }
      done
      if [[ "${TERRIER_UNKNOWN_ONLY_REQUIRED:-TRUE}" != "TRUE" ]]; then
        grep -Fq -- "--no-unknown-only" <<< "$help_output" || {
          echo "[ERROR] TERRIER_UNKNOWN_ONLY=FALSE requested, but --no-unknown-only is unavailable." >&2
          exit 1
        }
      fi
      python - <<"PYTERRIER"
import importlib.metadata as m, os
required=os.environ["REQUIRED_TERRIER_VERSION"]
for name in ("bio-terrier", "terrier"):
    try:
        version=m.version(name)
        break
    except m.PackageNotFoundError:
        version=None
if version is None:
    raise SystemExit("[ERROR] Terrier package metadata not found")
print(f"Terrier package version: {version}")
if not version.startswith(required):
    raise SystemExit(f"[ERROR] Terrier {required}.x required; found {version}")
PYTERRIER
    '
  fi

  if ! { repeat_checkpoint_ready tetrimmer; }; then
    echo "[INFO] Auditing TEtrimmer environment"
    resolve_tetrimmer_script

    # shellcheck disable=SC2016
    echo "[INFO] TEtrimmer preflight: checking Python entry point, RepeatMasker, version, and CLI options"
    if ! env_run "$ENV_TETRIMMER" env \
      REQUIRED_TETRIMMER_VERSION="$TETRIMMER_REQUIRED_VERSION" \
      TETRIMMER_SCRIPT="$TETRIMMER_SCRIPT" \
      TETRIMMER_CLASSIFICATION_REQUIRED="$TETRIMMER_CLASSIFICATION" \
      TETRIMMER_MODE_REQUIRED="$TETRIMMER_MODE" \
      TETRIMMER_HMM_REQUIRED="$TETRIMMER_HMM" \
      TETRIMMER_GENOME_ANNO_REQUIRED="$TETRIMMER_GENOME_ANNO" \
      TETRIMMER_DEBUG_REQUIRED="$TETRIMMER_DEBUG" \
      TETRIMMER_COMPRESS_REQUIRED="$TETRIMMER_COMPRESS_OUTPUT" \
      TETRIMMER_SKIP_REST_REQUIRED="$TETRIMMER_SKIP_REST_SEQ" \
      TETRIMMER_EXPORT_COVERAGE_REQUIRED="$TETRIMMER_EXPORT_COVERAGE" \
      bash -lc '
      set -Eeuo pipefail

      command -v python >/dev/null 2>&1 || {
        echo "[ERROR] python not found in the TEtrimmer environment." >&2
        exit 1
      }
      command -v timeout >/dev/null 2>&1 || {
        echo "[ERROR] GNU timeout not found; required for bounded TEtrimmer preflight checks." >&2
        exit 1
      }

      echo "[CHECK] Python: $(command -v python)"

      required_commands=(RepeatMasker)
      case "${TETRIMMER_CLASSIFICATION_REQUIRED:-none}" in
        none|classify_unknown) ;;
        classify_all) required_commands+=(RepeatClassifier) ;;
        *)
          echo "[ERROR] Invalid TETRIMMER_CLASSIFICATION=$TETRIMMER_CLASSIFICATION_REQUIRED" >&2
          exit 1
          ;;
      esac
      for required_command in "${required_commands[@]}"; do
        command -v "$required_command" >/dev/null 2>&1 || {
          echo "[ERROR] $required_command not found in the TEtrimmer environment." >&2
          exit 1
        }
        echo "[CHECK] $required_command: $(command -v "$required_command")"
      done

      [[ -s "$TETRIMMER_SCRIPT" ]] || {
        echo "[ERROR] TEtrimmer Python entry point not found or empty: $TETRIMMER_SCRIPT" >&2
        exit 1
      }
      [[ -r "$TETRIMMER_SCRIPT" ]] || {
        echo "[ERROR] TEtrimmer Python entry point is not readable: $TETRIMMER_SCRIPT" >&2
        exit 1
      }

      echo "[CHECK] TEtrimmer Python entry point: $TETRIMMER_SCRIPT"

      # Check the same Python entry point used for production execution.
      echo "[CHECK] Running TEtrimmer version"
      if version_output="$(timeout 60s python "$TETRIMMER_SCRIPT" --version 2>&1)"; then
        rc=0
      else
        rc=$?
      fi

      if [[ "$rc" -ne 0 ]]; then
        printf "%s\n" "$version_output" >&2
        if [[ "$rc" -eq 124 ]]; then
          echo "[ERROR] TEtrimmer --version exceeded the 60-second timeout." >&2
        else
          echo "[ERROR] TEtrimmer --version failed with exit code ${rc}." >&2
        fi
        exit "$rc"
      fi

      printf "%s\n" "$version_output"
      observed_version="$(grep -Eo "[0-9]+(\\.[0-9]+)+" <<< "$version_output" | head -n1 || true)"
      [[ -n "$observed_version" ]] || {
        echo "[ERROR] Could not parse TEtrimmer version from: $version_output" >&2
        exit 1
      }

      if [[ "$(printf "%s\n%s\n" "$REQUIRED_TETRIMMER_VERSION" "$observed_version" | sort -V | head -n1)" != "$REQUIRED_TETRIMMER_VERSION" ]]; then
        echo "[ERROR] TEtrimmer >=${REQUIRED_TETRIMMER_VERSION} required; found ${observed_version}." >&2
        exit 1
      fi

      echo "[CHECK] Running: python $TETRIMMER_SCRIPT --help"
      if help_output="$(timeout 120s python "$TETRIMMER_SCRIPT" --help 2>&1)"; then
        rc=0
      else
        rc=$?
      fi

      if [[ "$rc" -ne 0 ]]; then
        printf "%s\n" "$help_output" >&2
        if [[ "$rc" -eq 124 ]]; then
          echo "[ERROR] TEtrimmer --help exceeded the 120-second timeout." >&2
        else
          echo "[ERROR] TEtrimmer --help failed with exit code ${rc}." >&2
        fi
        exit "$rc"
      fi

      for required_option in --input_file --genome_file --output_dir --num_threads --preset --continue_analysis --dedup --pfam_dir --max_thread_time; do
        grep -Fq -- "$required_option" <<< "$help_output" || {
          echo "[ERROR] TEtrimmer CLI is missing required option: $required_option" >&2
          exit 1
        }
      done

      require_requested_option() {
        enabled="$1"
        option="$2"
        [[ "$enabled" == "TRUE" ]] || return 0
        grep -Fq -- "$option" <<< "$help_output" || {
          echo "[ERROR] TEtrimmer option requested by configuration but unavailable: $option" >&2
          exit 1
        }
      }

      if [[ "${TETRIMMER_MODE_REQUIRED:-curated}" == "fast" ]]; then
        grep -Fq -- "--fast_mode" <<< "$help_output" || {
          echo "[ERROR] TETRIMMER_MODE=fast requested, but --fast_mode is unavailable." >&2
          exit 1
        }
      fi

      if [[ "${TETRIMMER_CLASSIFICATION_REQUIRED:-none}" == "classify_all" ]]; then
        require_requested_option TRUE "--classify_all"
      fi

      require_requested_option "${TETRIMMER_HMM_REQUIRED:-FALSE}" "--hmm"
      require_requested_option "${TETRIMMER_GENOME_ANNO_REQUIRED:-FALSE}" "--genome_anno"
      require_requested_option "${TETRIMMER_DEBUG_REQUIRED:-FALSE}" "--debug"
      require_requested_option "${TETRIMMER_COMPRESS_REQUIRED:-FALSE}" "--compress_output"
      require_requested_option "${TETRIMMER_SKIP_REST_REQUIRED:-FALSE}" "--skip_rest_seq"
      require_requested_option "${TETRIMMER_EXPORT_COVERAGE_REQUIRED:-FALSE}" "--export_coverage"

      echo "[OK] TEtrimmer preflight completed successfully"
    ' > >(tee "${DIR_LOGS}/02_tetrimmer.preflight.txt") 2>&1; then
      fatal "TEtrimmer environment audit failed. See ${DIR_LOGS}/02_tetrimmer.preflight.txt"
    fi
  fi

  if is_true "$EARLGREY_ENHANCED_PLOTS" && ! { repeat_checkpoint_ready earlgrey_plots; }; then
    echo "[INFO] Auditing R environment for Earl Grey enhanced plots"

    env_run "$ENV_R" Rscript --vanilla -e 'if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Renv is missing required R package: ggplot2"); cat("Renv ggplot2 OK\n")' \
      > "${DIR_LOGS}/05_earlgrey_enhanced_plots.renv_preflight.log" \
      2>&1 \
      || error "Renv preflight failed for Earl Grey enhanced plots; install ggplot2 in ${ENV_R}"
  fi

# ============================================================
# REPEAT SUBSTAGE 0 — PREPARE DFAM + REPBASE KNOWN LIBRARY
# ============================================================

if repeat_checkpoint_ready known; then
    echo "[INFO] Step 0: existing known-repeat library found; skipping preparation: $KNOWN_NR"
  else
    echo "[INFO] Step 0: Preparing Dfam + RepBase known-repeat library"

  echo "[INFO] Exporting Dfam FamDB to FASTA"

  singularity exec \
    --no-home \
    -B "$PROJECT_ROOT":"$PROJECT_ROOT" \
    -B "$DFAM_FAMDB_DIR":"$FAMDB_CONTAINER_DATA_DIR":ro \
    --env "FAMDB_DIR=$FAMDB_CONTAINER_DIR" \
    --env "DFAM_FLAGS=$DFAM_FLAGS" \
    --env "DFAM_TAXON=$DFAM_TAXON" \
    --env "DFAM_FASTA=$DFAM_FASTA" \
    "$sif_repeatmasker" \
    bash -lc '
      set -Eeuo pipefail
      read -r -a dfam_flags <<< "${DFAM_FLAGS:-}"

      "$FAMDB_DIR/famdb.py" \
        families \
        --format fasta_name \
        --include-class-in-name \
        --ancestors \
        --descendants \
        "${dfam_flags[@]}" \
        "$DFAM_TAXON" \
        > "$DFAM_FASTA"
    ' \
    > "${DIR_LOGS}/00A_dfam_export.log" \
    2>&1

  require_file "$DFAM_FASTA" "Dfam FASTA export"
  DFAM_RECORDS=$(grep -c "^>" "$DFAM_FASTA" || true)
  if [[ "$DFAM_RECORDS" -lt 1 ]]; then
    echo "[ERROR] Dfam export produced no FASTA records: $DFAM_FASTA" >&2
    echo "[ERROR] Check taxon and partitions. See log: ${DIR_LOGS}/00A_dfam_export.log" >&2
    exit 1
  fi
  echo "[INFO] Dfam FASTA records: $DFAM_RECORDS"

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

  REPBASE_RECORDS=0
  if is_true "$REPEAT_USE_REPBASE"; then
    convert_repbase_to_fasta "$REPBASE_EMBL" "$REPBASE_FASTA"
    REPBASE_RECORDS=$(grep -c "^>" "$REPBASE_FASTA" || true)
    if [[ "$REPBASE_RECORDS" -lt 1 ]]; then
      echo "[ERROR] RepBase conversion produced no FASTA records: $REPBASE_FASTA" >&2
      echo "[ERROR] See log: ${DIR_LOGS}/00B_repbase_convert.log" >&2
      exit 1
    fi
    echo "[INFO] RepBase FASTA records: $REPBASE_RECORDS"
  else
    : > "$REPBASE_FASTA"
    echo "[INFO] RepBase integration disabled (REPEAT_USE_REPBASE=FALSE)"
  fi

  {
    echo -e "source\tpath"
    echo -e "Dfam_FamDB\t${DFAM_FAMDB_DIR}"
    echo -e "FamDB_software_container\t${FAMDB_CONTAINER_DIR}"
    echo -e "FamDB_config_container\t${FAMDB_CONTAINER_DIR}/famdb.conf"
    echo -e "FamDB_database_container\t${FAMDB_CONTAINER_DATA_DIR}"
    echo -e "Dfam_FASTA\t${DFAM_FASTA}"
    echo -e "RepBase_enabled\t${REPEAT_USE_REPBASE}"
    if is_true "$REPEAT_USE_REPBASE"; then
      echo -e "RepBase_input\t${REPBASE_EMBL}"
      echo -e "RepBase_FASTA\t${REPBASE_FASTA}"
    fi
    [[ -n "$ARTEFACTS_FASTA" ]] && echo -e "Artefacts_FASTA\t${ARTEFACTS_FASTA}"
  } > "${DIR_LOGS}/00_known_library_inputs.tsv"

  echo "[INFO] Building source-aware known-repeat library (no early clustering)"
  : > "$KNOWN_RAW"

  awk '/^>/ {sub(/^>/, ">DFAM|"); print; next} {print}' "$DFAM_FASTA" >> "$KNOWN_RAW"
  if is_true "$REPEAT_USE_REPBASE" && [[ -s "$REPBASE_FASTA" ]]; then
    awk '/^>/ {sub(/^>/, ">REPBASE|"); print; next} {print}' "$REPBASE_FASTA" >> "$KNOWN_RAW"
  fi
  if [[ -n "$ARTEFACTS_FASTA" && -s "$ARTEFACTS_FASTA" ]]; then
    awk '/^>/ {sub(/^>/, ">ARTEFACT|"); print; next} {print}' "$ARTEFACTS_FASTA" >> "$KNOWN_RAW"
  fi

  cp -f "$KNOWN_RAW" "$KNOWN_NR"
  require_file "$KNOWN_NR" "known-repeat library"
  repeat_mark_checkpoint known
fi

echo "[INFO] Known-repeat library: $KNOWN_NR"

# ============================================================
# REPEAT SUBSTAGE 1 — REPEATMODELER2
# ============================================================

if repeat_checkpoint_ready repeatmodeler; then
    echo "[INFO] Step 1: existing RepeatModeler library found; skipping RepeatModeler: $RM_LIBRARY"
  else
    echo "[INFO] Step 1: Running RepeatModeler2"

    archive_unrecoverable_repeatmodeler_dirs

  # Build RepeatModeler option string inside the host shell, then pass it to
  # the container. This allows recovery runs to use the exact same settings.
  RM_OPTS="-database '$GENOME_STEM' \
-threads '$THREADS' \
-famdb_dir '$FAMDB_CONTAINER_DIR'"

  if is_true "$USE_LTRSTRUCT"; then
    RM_OPTS="$RM_OPTS -LTRStruct"
  fi

  if is_true "$RM_QUICK"; then
    RM_OPTS="$RM_OPTS -quick"
  fi

  if [[ -n "$RM_SRAND" ]]; then
    RM_OPTS="$RM_OPTS -srand '$RM_SRAND'"
  fi

  if [[ -n "$RM_GENOME_SAMPLE_SIZE_MAX" ]]; then
    RM_OPTS="$RM_OPTS -genomeSampleSizeMax '$RM_GENOME_SAMPLE_SIZE_MAX'"
  fi

  if [[ "$RM_ADDL_ROUNDS" -gt 0 ]]; then
    RM_OPTS="$RM_OPTS -numAddlRounds '$RM_ADDL_ROUNDS'"
  fi

  if [[ -n "$REPEATMODELER_RECOVER_DIR" && ! -d "$REPEATMODELER_RECOVER_DIR" ]]; then
    fatal "REPEATMODELER_RECOVER_DIR is set but missing: $REPEATMODELER_RECOVER_DIR"
  fi

  {
    echo "# RepeatModeler command logic"
    echo "${REPEATMODELER_CONTAINER_DIR}/BuildDatabase -name '$GENOME_STEM' '$GENOME_LOCAL'"
    echo "${REPEATMODELER_CONTAINER_DIR}/RepeatModeler $RM_OPTS"
    if [[ -n "$REPEATMODELER_RECOVER_DIR" ]]; then
      echo "${REPEATMODELER_CONTAINER_DIR}/RepeatModeler $RM_OPTS -recoverDir '$REPEATMODELER_RECOVER_DIR'"
    elif is_true "$REPEATMODELER_RECOVER"; then
      echo "${REPEATMODELER_CONTAINER_DIR}/RepeatModeler $RM_OPTS -recoverDir '<newest RM_* under ${DIR_RM}>'"
    fi
    echo "DFAM_FAMDB_DIR_HOST=$DFAM_FAMDB_DIR"
    echo "FAMDB_DIR_CONTAINER=$FAMDB_CONTAINER_DIR"
    echo "FAMDB_CONFIG_CONTAINER=$FAMDB_CONTAINER_DIR/famdb.conf"
    echo "FAMDB_DATABASE_CONTAINER=$FAMDB_CONTAINER_DATA_DIR"
    echo "REPEATMODELER_DIR_CONTAINER=$REPEATMODELER_CONTAINER_DIR"
    echo "REPEATMODELER_RECOVER=$REPEATMODELER_RECOVER"
    echo "REPEATMODELER_RECOVER_DIR=$REPEATMODELER_RECOVER_DIR"
  } > "${DIR_LOGS}/01_repeatmodeler.command.sh"

  singularity exec \
    --no-home \
    -B "$PROJECT_ROOT":"$PROJECT_ROOT" \
    -B "$DFAM_FAMDB_DIR":"$FAMDB_CONTAINER_DATA_DIR":ro \
    --env "FAMDB_DIR=$FAMDB_CONTAINER_DIR" \
    --env "DFAM_REQUIRED_PARTITIONS=$DFAM_REQUIRED_PARTITIONS" \
    "$sif_repeatmasker" \
    bash -lc "
      set -Eeuo pipefail

      export FAMDB_DIR='$FAMDB_CONTAINER_DIR'
      FAMDB_DB_DIR='$FAMDB_CONTAINER_DATA_DIR'

      [[ -x \"\$FAMDB_DIR/famdb.py\" ]] || {
        echo '[ERROR] Missing FamDB executable:' \"\$FAMDB_DIR/famdb.py\" >&2
        exit 1
      }

      [[ -s \"\$FAMDB_DIR/famdb.conf\" ]] || {
        echo '[ERROR] Missing FamDB configuration:' \"\$FAMDB_DIR/famdb.conf\" >&2
        exit 1
      }

      read -r -a required_partitions <<< \"\$DFAM_REQUIRED_PARTITIONS\"
      [[ \"\${#required_partitions[@]}\" -gt 0 ]] || {
        echo '[ERROR] DFAM_REQUIRED_PARTITIONS is empty' >&2
        exit 1
      }
      for partition in \"\${required_partitions[@]}\"; do
        file=\"\$FAMDB_DB_DIR/\$partition.h5\"
        [[ -s \"\$file\" ]] || {
          echo '[ERROR] Missing required FamDB partition:' \"\$file\" >&2
          exit 1
        }
      done

      famdb_output=\$(\"\$FAMDB_DIR/famdb.py\" info 2>&1)
      printf '%s\n' \"\$famdb_output\"

      grep -Eq '^Database[[:space:]]*:' <<< \"\$famdb_output\" || {
        echo '[ERROR] FamDB database name could not be determined' >&2
        exit 1
      }

      grep -Eq '^Version[[:space:]]*:[[:space:]]*[^[:space:]]+' \
        <<< \"\$famdb_output\" || {
        echo '[ERROR] FamDB database version could not be determined' >&2
        exit 1
      }

      cd '$DIR_RM'

      '$REPEATMODELER_CONTAINER_DIR/BuildDatabase' \
        -name '$GENOME_STEM' \
        '$GENOME_LOCAL'

      newest_rm_dir=''

      if [[ -n '$REPEATMODELER_RECOVER_DIR' ]]; then
        newest_rm_dir='$REPEATMODELER_RECOVER_DIR'
      elif [[ '$REPEATMODELER_RECOVER' =~ ^(TRUE|True|true|1|YES|yes|Y|y)$ ]]; then
        newest_rm_dir=\$(find . -maxdepth 1 -type d -name 'RM_*' -printf '%T@\t%p\n' 2>/dev/null | sort -nr | head -n1 | cut -f2- || true)
      fi

      if [[ -n \"\$newest_rm_dir\" && -d \"\$newest_rm_dir\" ]]; then
        rm_consensi=\"\$newest_rm_dir/consensi.fa\"
        rm_stockholm=\"\$newest_rm_dir/families.stk\"
        rm_classified=\"\$newest_rm_dir/consensi.fa.classified\"
        rm_classified_stk=\"\$newest_rm_dir/families-classified.stk\"

        if [[ -s \"\$rm_classified\" ]]; then
          echo '[INFO] Existing classified RepeatModeler library:' \"\$rm_classified\"

          [[ -s \"\$rm_classified_stk\" ]] || {
            echo '[ERROR] Existing classified RepeatModeler directory lacks families-classified.stk:' \"\$newest_rm_dir\" >&2
            exit 1
          }

        elif [[ -s \"\$rm_consensi\" && -s \"\$rm_stockholm\" ]]; then
          echo '[INFO] Discovery completed; running RepeatClassifier directly.'

          (
            cd \"\$newest_rm_dir\"

            rm -f \
              consensi.fa.classified \
              families-classified.stk

            '$REPEATMODELER_CONTAINER_DIR/RepeatClassifier' \
              -consensi consensi.fa \
              -stockholm families.stk \
              -threads '$THREADS' \
              -famdb_dir '$FAMDB_CONTAINER_DIR'

            [[ -s consensi.fa.classified ]] || {
              echo '[ERROR] Missing consensi.fa.classified' >&2
              exit 1
            }

            [[ -s families-classified.stk ]] || {
              echo '[ERROR] Missing families-classified.stk' >&2
              exit 1
            }
          )

        else
          echo '[INFO] Attempting RepeatModeler recovery from:' \"\$newest_rm_dir\" >&2

          '$REPEATMODELER_CONTAINER_DIR/RepeatModeler' \
            $RM_OPTS \
            -recoverDir \"\$newest_rm_dir\"
        fi

        [[ -s \"\$rm_classified\" ]] || {
          echo '[ERROR] Classified RepeatModeler library was not produced' >&2
          exit 1
        }

        [[ -s \"\$rm_classified_stk\" ]] || {
          echo '[ERROR] Classified RepeatModeler Stockholm output was not produced' >&2
          exit 1
        }

        cp -f \"\$rm_classified\" '$RM_LIBRARY'
      else
        '$REPEATMODELER_CONTAINER_DIR/RepeatModeler' \
          $RM_OPTS

        newest_rm_dir=\$(find . \
          -maxdepth 1 \
          -type d \
          -name 'RM_*' \
          -printf '%T@\t%p\n' \
          | sort -nr \
          | head -n1 \
          | cut -f2-)

        [[ -s \"\$newest_rm_dir/consensi.fa.classified\" ]] || {
          echo '[ERROR] Fresh RepeatModeler run lacks classified output' >&2
          exit 1
        }

        [[ -s \"\$newest_rm_dir/families-classified.stk\" ]] || {
          echo '[ERROR] Fresh RepeatModeler run lacks classified Stockholm output' >&2
          exit 1
        }

        cp -f \
          \"\$newest_rm_dir/consensi.fa.classified\" \
          '$RM_LIBRARY'
      fi
    " \
    > "${DIR_LOGS}/01_repeatmodeler.log" \
    2>&1

  if [[ ! -s "$RM_LIBRARY" ]]; then
    echo "[ERROR] Published RepeatModeler library missing: $RM_LIBRARY" >&2
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
# REPEAT SUBSTAGE 2 — TERRIER
# ============================================================

# hierarchicalsoftmax checks a nonexistent score for deterministic one-child
# nodes. Repair only that guard in a pipeline-owned copy of the installed package.
prepare_terrier_inference_runtime() {
  local runtime_root="$1"
  local patcher="${DIR_LOGS}/terrier_inference_runtime_fix.py"
  cat > "$patcher" <<'PY_TERRIER_RUNTIME'
import importlib.util
import shutil
import sys
from pathlib import Path

runtime_root = Path(sys.argv[1])
if runtime_root.name != "terrier_inference_runtime" or runtime_root.is_symlink():
    raise SystemExit("Unsafe Terrier runtime directory")
spec = importlib.util.find_spec("hierarchicalsoftmax")
if spec is None or spec.origin is None:
    raise SystemExit("Terrier hierarchicalsoftmax package is missing")
source_dir = Path(spec.origin).parent
source_file = source_dir / "inference.py"
original = source_file.read_text()
old = "if threshold and predictions[node.softmax_start_index+prediction_child_index] < threshold:"
new = "if threshold and len(node.children) > 1 and predictions[node.softmax_start_index+prediction_child_index] < threshold:"
if original.count(old) != 1:
    raise SystemExit("Unsupported hierarchicalsoftmax greedy_predictions implementation")
patched = original.replace(old, new)
compile(patched, str(source_file), "exec")
if runtime_root.exists():
    shutil.rmtree(runtime_root)
runtime_root.mkdir(parents=True)
shutil.copytree(source_dir, runtime_root / "hierarchicalsoftmax",
                ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
(runtime_root / "hierarchicalsoftmax" / "inference.py").write_text(patched)
print("[INFO] Terrier runtime: skip threshold lookup for deterministic one-child nodes")
PY_TERRIER_RUNTIME
  env_run "$ENV_TERRIER" python "$patcher" "$runtime_root"

  env_run "$ENV_TERRIER" env "PYTHONPATH=${runtime_root}${PYTHONPATH:+:${PYTHONPATH}}" \
    python -c '
import sys
from pathlib import Path
from hierarchicalsoftmax import inference

expected = Path(sys.argv[1], "hierarchicalsoftmax", "inference.py").resolve()
if Path(inference.__file__).resolve() != expected:
    raise SystemExit("Terrier did not load the repaired hierarchicalsoftmax runtime")
' "$runtime_root"
}

TERRIER_LIBRARY="${TERRIER_LIBRARY_OVERRIDE:-${DIR_TERRIER}/repeatmodeler.terrier.fa}"
if ! is_true "$TERRIER_ENABLE"; then
  echo "[INFO] Step 2: Terrier disabled; TEtrimmer will use the RepeatModeler library directly."
  TERRIER_LIBRARY="$RM_LIBRARY"
elif repeat_checkpoint_ready terrier; then
  echo "[INFO] Step 2: existing Terrier outputs found; skipping: $TERRIER_LIBRARY"
else
  echo "[INFO] Step 2: Classifying RepeatModeler Unknown families with Terrier"
  mkdir -p "$DIR_TERRIER/images"

  TERRIER_ARGS=(
    --input "$RM_LIBRARY"
    --output-csv "$DIR_TERRIER/terrier_predictions.csv"
    --output-tips-csv "$DIR_TERRIER/terrier_tip_probabilities.csv"
    --output-fasta "$TERRIER_LIBRARY"
    --image-dir "$DIR_TERRIER/images"
    --image-format "$TERRIER_IMAGE_FORMAT"
    --image-threshold "$TERRIER_IMAGE_THRESHOLD"
    --threshold "$TERRIER_THRESHOLD"
    --min-length "$TERRIER_MIN_LENGTH"
    --max-length "$TERRIER_MAX_LENGTH"
    --batch-size "$TERRIER_BATCH_SIZE"
    --num-threads "$TERRIER_THREADS"
  )
  if is_true "$TERRIER_UNKNOWN_ONLY"; then
    TERRIER_ARGS+=(--unknown-only)
  else
    TERRIER_ARGS+=(--no-unknown-only)
  fi

  TERRIER_RUNTIME_ROOT="${DIR_LOGS}/terrier_inference_runtime"
  prepare_terrier_inference_runtime "$TERRIER_RUNTIME_ROOT"
  TERRIER_RUNTIME_PYTHONPATH="${TERRIER_RUNTIME_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"

  TERRIER_DOT_DIR=""
  if ! env_run "$ENV_TERRIER" bash -c 'command -v dot >/dev/null 2>&1'; then
    terrier_dot_fallback="${MINIFORGE_PREFIX}/envs/${ENV_R}/bin/dot"
    [[ -x "$terrier_dot_fallback" ]] || fatal \
      "Terrier probability images require Graphviz dot; missing from ${ENV_TERRIER} and ${terrier_dot_fallback}"
    TERRIER_DOT_DIR="$(dirname "$terrier_dot_fallback")"
    echo "[INFO] Terrier Graphviz: ${terrier_dot_fallback}"
  fi
  env_run "$ENV_TERRIER" bash -c '
    if [[ -n "$1" ]]; then export PATH="$1:$PATH"; fi
    printf "digraph test { a -> b }\n" | dot -T"$2" -o /dev/null
  ' _ "$TERRIER_DOT_DIR" "$TERRIER_IMAGE_FORMAT" \
    || fatal "Terrier Graphviz cannot render ${TERRIER_IMAGE_FORMAT} images"

  terrier_launch='if [[ -n "$1" ]]; then export PATH="$1:$PATH"; fi; shift; exec terrier "$@"'

  {
    echo "# Terrier command"
    printf '%q run -n %q env PYTHONPATH=%q bash -c %q _ %q' \
      "$ENV_RUNNER" "$ENV_TERRIER" "$TERRIER_RUNTIME_PYTHONPATH" \
      "$terrier_launch" "$TERRIER_DOT_DIR"
    printf ' %q' "${TERRIER_ARGS[@]}"
    echo
  } > "${DIR_LOGS}/02_terrier.command.sh"

  env_run "$ENV_TERRIER" env "PYTHONPATH=$TERRIER_RUNTIME_PYTHONPATH" \
    bash -c "$terrier_launch" _ "$TERRIER_DOT_DIR" "${TERRIER_ARGS[@]}" \
    2>&1 | tee "${DIR_LOGS}/02_terrier.log" \
    || fatal "Terrier failed; see ${DIR_LOGS}/02_terrier.log"

  require_file "$TERRIER_LIBRARY" "Terrier-classified RepeatModeler library"
  require_file "$DIR_TERRIER/terrier_predictions.csv" "Terrier prediction table"
  require_file "$DIR_TERRIER/terrier_tip_probabilities.csv" "Terrier tip-probability table"
  repeat_mark_checkpoint terrier
fi

echo "[INFO] Terrier library: $TERRIER_LIBRARY"

# ============================================================
# REPEAT SUBSTAGE 3 — TETRIMMER
# ============================================================

if repeat_checkpoint_ready tetrimmer; then
    echo "[INFO] Step 3: existing TEtrimmer library found; skipping TEtrimmer: $TETRIMMER_LIBRARY"
  else
    echo "[INFO] Step 3: Running TEtrimmer"
  echo "[INFO] TEtrimmer mode:   $TETRIMMER_MODE"
  echo "[INFO] TEtrimmer preset: $TETRIMMER_PRESET"

  TETRIMMER_THREADS="$THREADS"
  if [[ ! "$TETRIMMER_THREADS" =~ ^[0-9]+$ || "$TETRIMMER_THREADS" -lt 1 ]]; then
    echo "[ERROR] Invalid TEtrimmer thread count derived from THREADS=$THREADS." >&2
    exit 1
  fi
  if [[ "$TETRIMMER_THREADS" -gt "$TETRIMMER_THREAD_CAP" ]]; then
    echo "[WARN] Capping TEtrimmer threads from $TETRIMMER_THREADS to $TETRIMMER_THREAD_CAP to avoid RAM overflow."
    TETRIMMER_THREADS="$TETRIMMER_THREAD_CAP"
  fi

  TETRIMMER_ARGS=(
    --input_file "$TERRIER_LIBRARY"
    --genome_file "$GENOME_LOCAL"
    --output_dir "$DIR_TETRIMMER"
    --num_threads "$TETRIMMER_THREADS"
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

  TETRIMMER_CONTINUING="FALSE"
  case "$TETRIMMER_CONTINUE_MODE" in
    auto)
      if is_true "$TETRIMMER_CONTINUE" && tetrimmer_dir_resumable; then
        echo "[INFO] TEtrimmer continue mode is auto; using compatible interrupted output: $TETRIMMER_RESULT_DIR"
        # Replace --output_dir with the actual directory that contains the
        # interrupted TEtrimmer analysis.
        for ((i=0; i<${#TETRIMMER_ARGS[@]}; i++)); do
          if [[ "${TETRIMMER_ARGS[$i]}" == "--output_dir" ]]; then
            TETRIMMER_ARGS[$((i+1))]="$TETRIMMER_RESULT_DIR"
            break
          fi
        done
        TETRIMMER_ARGS+=(--continue_analysis)
        TETRIMMER_CONTINUING="TRUE"
      else
        echo "[INFO] TEtrimmer continue mode is auto; starting a fresh analysis."
      fi
      ;;
    always)
      if is_true "$TETRIMMER_CONTINUE"; then
        tetrimmer_dir_resumable || {
          echo "[ERROR] TETRIMMER_CONTINUE_MODE=always but no resumable TEtrimmer directory was found." >&2
          exit 1
        }
        echo "[INFO] TEtrimmer continue mode is always; using --continue_analysis in: $TETRIMMER_RESULT_DIR"
        for ((i=0; i<${#TETRIMMER_ARGS[@]}; i++)); do
          if [[ "${TETRIMMER_ARGS[$i]}" == "--output_dir" ]]; then
            TETRIMMER_ARGS[$((i+1))]="$TETRIMMER_RESULT_DIR"
            break
          fi
        done
        TETRIMMER_ARGS+=(--continue_analysis)
        TETRIMMER_CONTINUING="TRUE"
      fi
      ;;
    never)
      echo "[INFO] TEtrimmer continue mode is never; not using --continue_analysis"
      ;;
  esac

  EFFECTIVE_TETRIMMER_CLASSIFICATION="$TETRIMMER_CLASSIFICATION"

  case "$EFFECTIVE_TETRIMMER_CLASSIFICATION" in
    classify_unknown)
      # Legacy configuration value: classification remains disabled.
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

  mkdir -p "$TETRIMMER_PFAM_DIR"
  TETRIMMER_ARGS+=(--pfam_dir "$TETRIMMER_PFAM_DIR")

  # TEtrimmer 1.7.4 accepts an empty --pfam_dir and downloads/presses Pfam
  # automatically. Treat Pfam availability as an independent resource check,
  # not as a prerequisite for starting TEtrimmer.
  tetrimmer_pfam_complete=TRUE
  for pfam_file in     Pfam-A.hmm     Pfam-A.hmm.dat     Pfam-A.hmm.h3f     Pfam-A.hmm.h3i     Pfam-A.hmm.h3m     Pfam-A.hmm.h3p
  do
    if [[ ! -s "${TETRIMMER_PFAM_DIR}/${pfam_file}" ]]; then
      tetrimmer_pfam_complete=FALSE
      break
    fi
  done

  if is_true "$tetrimmer_pfam_complete"; then
    echo "[INFO] Reusing complete TEtrimmer-managed Pfam database: $TETRIMMER_PFAM_DIR"
  else
    echo "[INFO] Complete Pfam database is not present in $TETRIMMER_PFAM_DIR; TEtrimmer will initialize it automatically."
  fi

  TETRIMMER_ARGS+=(--max_thread_time "$TETRIMMER_MAX_THREAD_TIME")

  if is_true "$TETRIMMER_COMPRESS_OUTPUT"; then
    TETRIMMER_ARGS+=(--compress_output)
  fi
  if is_true "$TETRIMMER_EXPORT_COVERAGE"; then
    TETRIMMER_ARGS+=(--export_coverage)
  fi
  if is_true "$TETRIMMER_SKIP_REST_SEQ"; then
    if ! is_true "$TETRIMMER_CONTINUING"; then
      echo "[ERROR] TETRIMMER_SKIP_REST_SEQ=TRUE requires a resumable run with --continue_analysis." >&2
      exit 1
    fi
    TETRIMMER_ARGS+=(--skip_rest_seq)
  fi

  if is_true "$TETRIMMER_EXPERT_OVERRIDES"; then
    echo "[WARN] TEtrimmer expert overrides are enabled and take precedence over preset values."
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
  fi

  case "$TETRIMMER_MODE" in
    fast)
      TETRIMMER_ARGS+=(--fast_mode)
      ;;
    curated|debug)
      ;;
    *)
      echo "[ERROR] Invalid TETRIMMER_MODE=$TETRIMMER_MODE. Use fast, curated, or debug." >&2
      exit 1
      ;;
  esac

  resolve_tetrimmer_script

  {
    echo "# TEtrimmer command"
    echo "# TEtrimmer requested threads: $THREADS"
    echo "# TEtrimmer safety cap: $TETRIMMER_THREAD_CAP"
    echo "# TEtrimmer effective threads: $TETRIMMER_THREADS"
    echo "# TEtrimmer required version: $TETRIMMER_REQUIRED_VERSION"
    echo "# TEtrimmer Python entry point: $TETRIMMER_SCRIPT"
    echo "# TEtrimmer continuing: $TETRIMMER_CONTINUING"
    echo "# TEtrimmer expert overrides: $TETRIMMER_EXPERT_OVERRIDES"
    echo "# TEtrimmer persistent Pfam directory: $TETRIMMER_PFAM_DIR"
    echo "# TEtrimmer requested output root: $DIR_TETRIMMER"
    printf '%q run -n %q python %q' \
      "$ENV_RUNNER" "$ENV_TETRIMMER" "$TETRIMMER_SCRIPT"
    printf ' %q' "${TETRIMMER_ARGS[@]}"
    echo
  } > "${DIR_LOGS}/02_tetrimmer.command.sh"

  env_run "$ENV_TETRIMMER" \
    python "$TETRIMMER_SCRIPT" "${TETRIMMER_ARGS[@]}" \
    > "${DIR_LOGS}/02_tetrimmer.log" \
    2>&1

  # Validate TEtrimmer-managed Pfam independently from the refined-library
  # output. TEtrimmer 1.7.4 may populate --pfam_dir during this run.
  pfam_validation_failed=FALSE
  {
    printf 'field\tvalue\n'
    printf 'pfam_dir\t%s\n' "$TETRIMMER_PFAM_DIR"

    for pfam_file in       Pfam-A.hmm       Pfam-A.hmm.dat       Pfam-A.hmm.h3f       Pfam-A.hmm.h3i       Pfam-A.hmm.h3m       Pfam-A.hmm.h3p
    do
      pfam_path="${TETRIMMER_PFAM_DIR}/${pfam_file}"
      if [[ -s "$pfam_path" ]]; then
        printf '%s\t%s\n' "$pfam_file" "$pfam_path"
      else
        printf '%s\tMISSING\n' "$pfam_file"
        pfam_validation_failed=TRUE
      fi
    done

    if [[ -s "${TETRIMMER_PFAM_DIR}/Pfam-A.hmm" ]]; then
      printf 'pfam_hmm_sha256\t%s\n'         "$(sha256sum "${TETRIMMER_PFAM_DIR}/Pfam-A.hmm" | awk '{print $1}')"
    fi
  } > "${DIR_LOGS}/02_tetrimmer.pfam.tsv"

  if is_true "$pfam_validation_failed"; then
    echo "[ERROR] TEtrimmer Pfam initialization is incomplete: $TETRIMMER_PFAM_DIR" >&2
    echo "[ERROR] Expected non-empty files:" >&2
    for pfam_file in       Pfam-A.hmm       Pfam-A.hmm.dat       Pfam-A.hmm.h3f       Pfam-A.hmm.h3i       Pfam-A.hmm.h3m       Pfam-A.hmm.h3p
    do
      [[ -s "${TETRIMMER_PFAM_DIR}/${pfam_file}" ]] ||         echo "        ${TETRIMMER_PFAM_DIR}/${pfam_file}" >&2
    done
    echo "[ERROR] See ${DIR_LOGS}/02_tetrimmer.log and ${DIR_LOGS}/02_tetrimmer.pfam.tsv" >&2
    exit 1
  fi

  echo "[INFO] TEtrimmer Pfam database is complete: $TETRIMMER_PFAM_DIR"

  audit_tetrimmer_outputs

  if tetrimmer_famdb_errors_present; then
    fatal "TEtrimmer classification could not access Dfam. See ${DIR_LOGS}/02_tetrimmer.log and ${DIR_LOGS}/02_tetrimmer.preflight.txt"
  fi

  resolve_tetrimmer_result_dir || true
  if [[ -n "$TETRIMMER_RESULT_DIR" && -s "${TETRIMMER_RESULT_DIR}/TEtrimmer_consensus_merged.fasta" ]]; then
    TETRIMMER_LIBRARY="${TETRIMMER_RESULT_DIR}/TEtrimmer_consensus_merged.fasta"
  elif [[ -n "$TETRIMMER_RESULT_DIR" && -s "${TETRIMMER_RESULT_DIR}/TEtrimmer_consensus.fasta" ]]; then
    echo "[WARN] TEtrimmer_consensus_merged.fasta not found after TEtrimmer run."
    echo "[WARN] Falling back to TEtrimmer_consensus.fasta."
    TETRIMMER_LIBRARY="${TETRIMMER_RESULT_DIR}/TEtrimmer_consensus.fasta"
  else
    TETRIMMER_LIBRARY=""
  fi

  if [[ -z "$TETRIMMER_LIBRARY" || ! -s "$TETRIMMER_LIBRARY" ]]; then
    echo "[ERROR] TEtrimmer did not produce a refined consensus library." >&2
    echo "[ERROR] No complete direct or TEtrimmer_output_* result directory was detected under:" >&2
    echo "        $DIR_TETRIMMER" >&2
    echo "[ERROR] Pfam validation is independent and has already passed." >&2
    echo "[ERROR] Inspect:" >&2
    echo "        ${DIR_LOGS}/02_tetrimmer.log" >&2
    echo "        ${DIR_TETRIMMER}/TEtrimmer.log" >&2
    echo "[INFO] Candidate TEtrimmer result directories:" >&2
    find "$DIR_TETRIMMER" -mindepth 1 -maxdepth 1 -type d -name 'TEtrimmer_output_*' -printf '        %p\n' 2>/dev/null >&2 || true
    echo "[INFO] Files currently present under TEtrimmer output (max depth 3):" >&2
    find "$DIR_TETRIMMER" -maxdepth 3 -type f -printf '        %p\n' 2>/dev/null >&2 || true
    exit 1
  fi

  echo "[INFO] Resolved TEtrimmer result directory: $TETRIMMER_RESULT_DIR"

  repeat_mark_checkpoint tetrimmer
fi

echo "[INFO] TEtrimmer library: $TETRIMMER_LIBRARY"

# ============================================================
# REPEAT SUBSTAGE 4 — FINAL LIBRARY
# ============================================================

if repeat_checkpoint_ready final_library; then
    echo "[INFO] Step 4: existing final TE library found; skipping final-library construction: $FINAL_NR"
  else
    echo "[INFO] Step 4: Building final annotation library"

  : > "$FINAL_RAW"
  cat "$KNOWN_NR" >> "$FINAL_RAW"
  awk '/^>/ {sub(/^>/, ">DENOVO|"); print; next} {print}' "$TETRIMMER_LIBRARY" >> "$FINAL_RAW"

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
        -aS '$FINAL_CDHIT_AS' \
        -aL '$FINAL_CDHIT_AL' \
        -g 1 \
        -r 1 \
        -T '$THREADS' \
        -M '$FINAL_CDHIT_MEMORY_MB' \
        -d 0
    " \
    > "${DIR_LOGS}/03_final_library.log" \
    2>&1

  require_file "$FINAL_NR" "final TE annotation library"
  require_file "${FINAL_NR}.clstr" "CD-HIT final-library cluster file"

  "$python_bin" - "$FINAL_RAW" "$FINAL_NR" "${FINAL_NR}.clstr" "${DIR_LIB}/final_library.provenance.tsv" <<'PYPROV'
import re, sys
from pathlib import Path
raw, final, clstr, out = map(Path, sys.argv[1:])

def fasta_headers(p):
    hs=[]
    with p.open() as fh:
        for line in fh:
            if line.startswith('>'):
                hs.append(line[1:].strip())
    return hs

clusters=[]
current=[]
for line in clstr.read_text().splitlines():
    if line.startswith('>Cluster'):
        if current: clusters.append(current)
        current=[]
        continue
    m=re.search(r'>([^\.]+(?:\.[^\.]+)*?)\.\.\.', line)
    if m:
        seqid=m.group(1)
        current.append((seqid, line.rstrip().endswith('*')))
if current: clusters.append(current)

with out.open('w') as fh:
    fh.write('cluster_id\trepresentative_id\trepresentative_source\tmembers\tsources\n')
    for i, members in enumerate(clusters, 1):
        rep=next((x for x,isrep in members if isrep), members[0][0] if members else '')
        src=lambda x: x.split('|',1)[0] if '|' in x else 'UNKNOWN'
        fh.write(f'{i}\t{rep}\t{src(rep)}\t' + ','.join(x for x,_ in members) + '\t' + ','.join(sorted(set(src(x) for x,_ in members))) + '\n')
PYPROV

  repeat_mark_checkpoint final_library
fi

echo "[INFO] Final TE library: $FINAL_NR"

# ============================================================
# REPEAT SUBSTAGE 5 — EARL GREY ANNOTATION-ONLY
# ============================================================

if repeat_checkpoint_ready earlgrey; then
    echo "[INFO] Step 5: existing Earl Grey outputs found; skipping"
    audit_earlgrey_outputs
  else
    echo "[INFO] Step 5: Running Earl Grey >=${EARLGREY_REQUIRED_VERSION} annotation-only"

  earlgrey_prepare_library_ids "$FINAL_NR"
  echo "[INFO] Earl Grey uses only the consolidated custom library; Dfam/RepBase are not re-injected via -r."
  # Earl Grey 7.3.1 earlGreyAnnotationOnly supports -g/-s/-o/-t/-r/-l/-m/-d/-e only.
  # Do not pass the full-pipeline earlGrey -v resume-validation option here.

  require_file "$EARLGREY_LIBRARY" "Earl Grey consensus library"
  earlgrey_preflight
  prepare_earlgrey_output_dir
  earlgrey_archive_incompatible_masker_cache

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
    echo "# Earl Grey >=${EARLGREY_REQUIRED_VERSION} annotation-only command"
    echo "# FamDB real executable: $EARLGREY_FAMDB_REAL_PY"
    echo "# FamDB runtime wrapper: $EARLGREY_FAMDB_WRAPPER"
    printf '%q run -n %q env FAMDB_DIR=%q FAMDB_DATA_DIR=%q EARLGREY_SCRIPT_DIR=%q EARLGREY_DIVERGENCE_CALC=%q %q' \
      "$ENV_RUNNER" "$ENV_EARLGREY" "$EARLGREY_FAMDB_RUNTIME_DIR" "$DFAM_FAMDB_DIR" "$EARLGREY_SCRIPT_DIR" "$EARLGREY_DIVERGENCE_CALC" "$EARLGREY_ANNOTATION_CLI"
    printf ' %q' "${EARLGREY_ARGS[@]}"
    echo
  } > "${DIR_LOGS}/04_earlgrey_annotation_only.command.sh"

  if env_run_clean_r "$ENV_EARLGREY" \
      env \
      FAMDB_DIR="$EARLGREY_FAMDB_RUNTIME_DIR" \
      FAMDB_DATA_DIR="$DFAM_FAMDB_DIR" \
      EARLGREY_SCRIPT_DIR="$EARLGREY_SCRIPT_DIR" \
      EARLGREY_DIVERGENCE_CALC="$EARLGREY_DIVERGENCE_CALC" \
      "$EARLGREY_ANNOTATION_CLI" "${EARLGREY_ARGS[@]}" \
      > "${DIR_LOGS}/04_earlgrey_annotation_only.log" \
      2>&1; then
    :
  else
    earlgrey_rc=$?
    if earlgrey_repeatcraft_merge_failed; then
      fatal "Earl Grey RepeatCraft failed on a pipe-delimited repeat ID (rc=${earlgrey_rc}); check the EarlGrey-safe library and ID map. See ${DIR_LOGS}/04_earlgrey_annotation_only.log"
    fi
    if earlgrey_divergence_filename_errors_present; then
      fatal "Earl Grey divergence calculation failed because a repeat identifier was interpreted as a path (rc=${earlgrey_rc}). See ${DIR_LOGS}/04_earlgrey_annotation_only.log"
    fi
    fatal "Earl Grey annotation-only failed (rc=${earlgrey_rc}). See ${DIR_LOGS}/04_earlgrey_annotation_only.log"
  fi

  if earlgrey_famdb_errors_present; then
    fatal "Earl Grey could not access Dfam. See ${DIR_LOGS}/04_earlgrey_annotation_only.log and ${DIR_LOGS}/04_earlgrey.preflight.txt"
  fi

  earlgrey_core_outputs_complete || {
    echo "[ERROR] Earl Grey v${EARLGREY_REQUIRED_VERSION} finished without required final outputs." >&2
    echo "[ERROR] Required final outputs include:" >&2
    echo "        ${DIR_EARLGREY}/${SPECIES}_EarlGrey/${SPECIES}_summaryFiles/${SPECIES}.filteredRepeats.bed" >&2
    echo "        ${DIR_EARLGREY}/${SPECIES}_EarlGrey/${SPECIES}_summaryFiles/${SPECIES}.filteredRepeats.gff" >&2
    echo "        ${DIR_EARLGREY}/${SPECIES}_EarlGrey/${SPECIES}_summaryFiles/${SPECIES}.highLevelCount.txt" >&2
    echo "        ${DIR_EARLGREY}/${SPECIES}_EarlGrey/${SPECIES}_summaryFiles/${SPECIES}.familyLevelCount.txt" >&2
    echo "        ${DIR_EARLGREY}/${SPECIES}_EarlGrey/${SPECIES}_summaryFiles/${SPECIES}_divergence_summary_table.tsv" >&2
    echo "        ${DIR_EARLGREY}/${SPECIES}_EarlGrey/${SPECIES}_RepeatLandscape/${SPECIES}.filteredRepeats.withDivergence.gff" >&2
    echo "[ERROR] Check log: ${DIR_LOGS}/04_earlgrey_annotation_only.log" >&2
    exit 1
  }

  audit_earlgrey_outputs
  repeat_mark_checkpoint earlgrey
fi

# ============================================================
# REPEAT SUBSTAGE 6 — ENHANCED EARL GREY PLOTS
# ============================================================

if is_true "$EARLGREY_ENHANCED_PLOTS"; then
  if ! earlgrey_core_outputs_complete; then
    echo "[ERROR] Earl Grey final outputs are incomplete; refusing to create enhanced plots." >&2
    echo "[ERROR] Rerun repeat substage 5 after resolving stale ${DIR_EARLGREY} outputs." >&2
    exit 1
  fi

  EARLGREY_PLOTTER="${DIR_LOGS}/05_earlgrey_enhanced_plots.R"
  EARLGREY_PLOTTER_CANDIDATE="${EARLGREY_PLOTTER}.candidate.$$"
  write_earlgrey_enhanced_plotter "$EARLGREY_PLOTTER_CANDIDATE"
  if repeat_checkpoint_ready earlgrey_plots \
    && [[ -s "$EARLGREY_PLOTTER" ]] \
    && cmp -s "$EARLGREY_PLOTTER_CANDIDATE" "$EARLGREY_PLOTTER"; then
    rm -f "$EARLGREY_PLOTTER_CANDIDATE"
    echo "[INFO] Step 6: existing Earl Grey enhanced plots found; skipping"
  else
    echo "[INFO] Step 6: Creating Earl Grey enhanced ggplot2 plots"
    mv -f "$EARLGREY_PLOTTER_CANDIDATE" "$EARLGREY_PLOTTER"

    {
      echo "# Earl Grey enhanced ggplot2 plotting command"
      printf '%q run -n %q Rscript --vanilla %q' "$ENV_RUNNER" "$ENV_R" "$EARLGREY_PLOTTER"
      printf ' --species %q --earlgrey_dir %q --outdir %q --dpi %q --top_families %q --binwidth %q --max_kimura %q --composition_tolerance %q\n' \
        "$SPECIES" "$DIR_EARLGREY" "$DIR_EARLGREY_ENHANCED_PLOTS" "$EARLGREY_PLOT_DPI" "$EARLGREY_PLOT_TOP_FAMILIES" "$EARLGREY_PLOT_BINWIDTH" "$EARLGREY_PLOT_MAX_KIMURA" "$EARLGREY_COMPOSITION_TOLERANCE"
    } > "${DIR_LOGS}/05_earlgrey_enhanced_plots.command.sh"

    env_run "$ENV_R" Rscript --vanilla "$EARLGREY_PLOTTER" \
      --species "$SPECIES" \
      --earlgrey_dir "$DIR_EARLGREY" \
      --outdir "$DIR_EARLGREY_ENHANCED_PLOTS" \
      --dpi "$EARLGREY_PLOT_DPI" \
      --top_families "$EARLGREY_PLOT_TOP_FAMILIES" \
      --binwidth "$EARLGREY_PLOT_BINWIDTH" \
      --max_kimura "$EARLGREY_PLOT_MAX_KIMURA" \
      --composition_tolerance "$EARLGREY_COMPOSITION_TOLERANCE" \
      > "${DIR_LOGS}/05_earlgrey_enhanced_plots.log" \
      2>&1

    repeat_mark_checkpoint earlgrey_plots
  fi
else
  echo "[INFO] Step 6: Earl Grey enhanced plots disabled"
fi

# ============================================================
# REPEAT SUBSTAGE 7 — COLLECT OUTPUTS
# ============================================================

if repeat_checkpoint_ready collect_outputs; then
  echo "[INFO] Step 7: existing collected repeat outputs found; skipping collection"
else
  echo "[INFO] Step 7: Collecting final outputs"

  earlgrey_core_outputs_complete || {
    echo "[ERROR] Earl Grey final outputs are incomplete; refusing to collect final outputs." >&2
    echo "[ERROR] Rerun Step 4 after resolving stale ${DIR_EARLGREY} outputs." >&2
    exit 1
  }

  EG_ROOT="${DIR_EARLGREY}/${SPECIES}_EarlGrey"
  EG_RM_DIR="${EG_ROOT}/${SPECIES}_RepeatMasker_Against_Custom_Library"
  EG_SUMMARY_DIR="${EG_ROOT}/${SPECIES}_summaryFiles"
  EG_LANDSCAPE_DIR="${EG_ROOT}/${SPECIES}_RepeatLandscape"
  EG_MERGED_DIR="${EG_ROOT}/${SPECIES}_mergedRepeats"

  mkdir -p \
    "$DIR_FINAL/known_libraries" \
    "$DIR_FINAL/repeatmodeler" \
    "$DIR_FINAL/terrier" \
    "$DIR_FINAL/tetrimmer" \
    "$DIR_FINAL/final_library" \
    "$DIR_FINAL/earlgrey/native_plots" \
    "$DIR_FINAL/earlgrey/final_annotation" \
    "$DIR_FINAL/logs"

  collect_manifest="${DIR_FINAL}/logs/collect_outputs_manifest.tsv"
  collect_manifest_tmp="${collect_manifest}.tmp.$$"
  printf 'category\tsource\tdestination\tsize_bytes\tlink_or_copy\n' > "$collect_manifest_tmp"

  collect_record() {
    local category="$1"
    local source="$2"
    local destination="$3"
    local method="$4"
    local size=0
    if [[ -f "$source" ]]; then
      size="$(stat -c '%s' "$source")"
    elif [[ -d "$source" ]]; then
      size="$(du -sb "$source" 2>/dev/null | awk 'NR == 1 {print $1}')"
      size="${size:-0}"
    fi
    printf '%s\t%s\t%s\t%s\t%s\n' \
      "$category" "$source" "$destination" "$size" "$method" >> "$collect_manifest_tmp"
  }

  collect_copy() {
    local category="$1"
    local source="$2"
    local destination="$3"
    [[ -s "$source" ]] || return 0
    if [[ -d "$destination" ]]; then
      destination="${destination}/$(basename "$source")"
    fi
    cp -f "$source" "$destination"
    collect_record "$category" "$source" "$destination" "copy"
  }

  collect_link() {
    local category="$1"
    local source="$2"
    local destination="$3"
    [[ -d "$source" ]] || error "Missing Earl Grey directory for collection: $source"
    if [[ -e "$destination" && ! -L "$destination" ]]; then
      error "Collection link destination exists and is not a symlink: $destination"
    fi
    ln -sfn "$source" "$destination"
    collect_record "$category" "$source" "$destination" "link"
  }

  for f in "${DFAM_FASTA:-}" "${REPBASE_FASTA:-}" "${KNOWN_RAW:-}" "${KNOWN_NR:-}"; do
    [[ -n "$f" ]] && collect_copy "known_library" "$f" "$DIR_FINAL/known_libraries"
  done
  collect_copy "repeatmodeler" "${RM_LIBRARY:-}" "$DIR_FINAL/repeatmodeler"

  for f in \
    "${TERRIER_LIBRARY:-}" \
    "$DIR_TERRIER/terrier_predictions.csv" \
    "$DIR_TERRIER/terrier_tip_probabilities.csv" \
    "$DIR_LOGS/02_terrier.command.sh" \
    "$DIR_LOGS/02_terrier.log"
  do
    collect_copy "terrier" "$f" "$DIR_FINAL/terrier"
  done

  for f in \
    "${TETRIMMER_LIBRARY:-}" \
    "$DIR_TETRIMMER/TEtrimmer_consensus.fasta" \
    "$DIR_TETRIMMER/TEtrimmer_consensus_merged.fasta" \
    "$DIR_TETRIMMER/summary.txt" \
    "$DIR_TETRIMMER/Sequence_name_mapping.txt" \
    "$DIR_TETRIMMER/error_file.txt" \
    "$DIR_LOGS/02_tetrimmer.output_audit.txt"
  do
    collect_copy "tetrimmer" "$f" "$DIR_FINAL/tetrimmer"
  done

  for f in "${FINAL_RAW:-}" "${FINAL_NR:-}"; do
    [[ -n "$f" ]] && collect_copy "final_library" "$f" "$DIR_FINAL/final_library"
  done

  collect_link "earlgrey_repeatmasker" "$EG_RM_DIR" "$DIR_FINAL/earlgrey/repeatmasker_internal"
  collect_link "earlgrey_annotation" "$EG_SUMMARY_DIR" "$DIR_FINAL/earlgrey/summaryFiles"
  collect_link "earlgrey_annotation" "$EG_LANDSCAPE_DIR" "$DIR_FINAL/earlgrey/RepeatLandscape"
  collect_link "earlgrey_annotation" "$EG_MERGED_DIR" "$DIR_FINAL/earlgrey/mergedRepeats"

  for f in \
    "${EG_SUMMARY_DIR}/${SPECIES}.filteredRepeats.bed" \
    "${EG_SUMMARY_DIR}/${SPECIES}.filteredRepeats.gff" \
    "${EG_SUMMARY_DIR}/${SPECIES}.highLevelCount.txt" \
    "${EG_SUMMARY_DIR}/${SPECIES}.familyLevelCount.txt" \
    "${EG_SUMMARY_DIR}/${SPECIES}_divergence_summary_table.tsv" \
    "${EG_SUMMARY_DIR}/${SPECIES}.softmasked.fasta" \
    "${EG_LANDSCAPE_DIR}/${SPECIES}.filteredRepeats.withDivergence.gff"
  do
    collect_copy "earlgrey_annotation" "$f" "$DIR_FINAL/earlgrey/final_annotation"
  done

  for native_name in \
    "${SPECIES}_classification_landscape.pdf" \
    "${SPECIES}_split_class_landscape.pdf" \
    "${SPECIES}_superfamily_div_plot.pdf" \
    "${SPECIES}.summaryPie.pdf"
  do
    native_source="${EG_SUMMARY_DIR}/${native_name}"
    [[ -s "$native_source" ]] || native_source="${EG_LANDSCAPE_DIR}/${native_name}"
    collect_copy "earlgrey_native_plot" "$native_source" "$DIR_FINAL/earlgrey/native_plots"
  done

  while IFS= read -r -d '' f; do
    collect_record "earlgrey_enhanced_plot" "$f" "$f" "generated"
  done < <(find "$DIR_EARLGREY_ENHANCED_PLOTS" -maxdepth 1 -type f -size +0c -print0 2>/dev/null)

  while IFS= read -r -d '' f; do
    case "$f" in
      "$collect_manifest"|"$collect_manifest_tmp") continue ;;
    esac
    case "$(basename "$f")" in
      *.command.sh|*.log) category="log" ;;
      *) category="provenance" ;;
    esac
    collect_copy "$category" "$f" "$DIR_FINAL/logs"
  done < <(
    find "$DIR_LOGS" -maxdepth 1 -type f \
      \( -name "*.log" -o -name "*.command.sh" -o -name "*.txt" -o -name "*.tsv" -o -name "*.yml" -o -name "*.R" \) \
      -size +0c -print0 2>/dev/null
  )

  mv -f "$collect_manifest_tmp" "$collect_manifest"
  awk 'NR > 1 {found = 1} END {exit found ? 0 : 1}' "$collect_manifest" \
    || error "Collect outputs did not record any deliverables in ${collect_manifest}"

  repeat_mark_checkpoint collect_outputs
fi

	# STEP 7 — NON-REPETITIVE REGIONS BED
	# ============================================================

	if repeat_checkpoint_ready non_repeats; then
	    echo "[INFO] Step 7: existing non-repetitive regions BED found; skipping"
	  else
	  echo "[INFO] Step 7: Creating non-repetitive regions BED"

	  require_file "$GENOME_LOCAL" "genome for non-repetitive BED"
	  require "$bedtools_bin"

	  EARLGREY_FINAL_BED="${DIR_EARLGREY}/${SPECIES}_EarlGrey/${SPECIES}_summaryFiles/${SPECIES}.filteredRepeats.bed"
	  REPEATS_BED="${NON_REPEATS_REPEATS_BED:-$EARLGREY_FINAL_BED}"

	  require_file "$REPEATS_BED" "Earl Grey filtered repeats BED"

	  if [[ ! -s "${GENOME_LOCAL}.fai" ]]; then
	    "$samtools_bin" faidx "$GENOME_LOCAL" \
	      > "${DIR_LOGS}/07_non_repeats.samtools_faidx.log" \
	      2>&1 \
	      || error "samtools faidx failed for non-repetitive BED input: $GENOME_LOCAL"
	  fi

	  cut -f1,2 "${GENOME_LOCAL}.fai" > "$NON_REPEATS_GENOME_SIZES"

	  awk 'BEGIN{OFS="\t"} NF >= 3 && $2 ~ /^-?[0-9]+$/ && $3 ~ /^[0-9]+$/ {start=$2; end=$3; if (start < 0) start=0; if (end > start) print $1,start,end}' "$REPEATS_BED" \
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
	    printf 'awk %q %q > %q\n' 'BEGIN{OFS="\t"} NF >= 3 && $2 ~ /^-?[0-9]+$/ && $3 ~ /^[0-9]+$/ {start=$2; end=$3; if (start < 0) start=0; if (end > start) print $1,start,end}' "$REPEATS_BED" "$NON_REPEATS_RAW_BED"
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

	  repeat_stage_manifest="${DIR_FINAL}/repeat_annotation.summary.tsv"
  {
    printf 'metric\tvalue\n'
    printf 'dfam_taxon\t%s\n' "$DFAM_TAXON"
    printf 'dfam_release_tag\t%s\n' "$DFAM_RELEASE_TAG"
    printf 'dfam_records\t%s\n' "$(grep -c '^>' "$DFAM_FASTA" 2>/dev/null || true)"
    printf 'repbase_enabled\t%s\n' "$REPEAT_USE_REPBASE"
    printf 'repbase_records\t%s\n' "$(grep -c '^>' "$REPBASE_FASTA" 2>/dev/null || true)"
    printf 'repeatmodeler_records\t%s\n' "$(grep -c '^>' "$RM_LIBRARY" 2>/dev/null || true)"
    printf 'tetrimmer_records\t%s\n' "$(grep -c '^>' "$TETRIMMER_LIBRARY" 2>/dev/null || true)"
    printf 'pfam_dir\t%s\n' "$TETRIMMER_PFAM_DIR"
    printf 'precluster_records\t%s\n' "$(grep -c '^>' "$FINAL_RAW" 2>/dev/null || true)"
    printf 'final_nonredundant_records\t%s\n' "$(grep -c '^>' "$FINAL_NR" 2>/dev/null || true)"
    printf 'cdhit_clusters\t%s\n' "$(grep -c '^>Cluster' "${FINAL_NR}.clstr" 2>/dev/null || true)"
    printf 'final_library_sha256\t%s\n' "$(sha256sum "$FINAL_NR" | awk '{print $1}')"
    printf 'earlgrey_output\t%s\n' "${DIR_EARLGREY}/${SPECIES}_EarlGrey"
  } > "$repeat_stage_manifest"

  echo "[INFO] Step 8 repeat annotation complete. Outputs: ${DIR_FINAL}"

fi

# ============================================================================================================================
#
# STEP 9: Final Summary & Run Report — execution stats, resources, key metrics, artifacts
# Produces:
#   - ${dir7}/run_summary.txt
#   - ${dir7}/run_summary.tsv
#   - ${dir7}/run_summary.json
#   - ${dir7}/logs/summary_tail.log
# ============================================================================================================================

echo "[INFO] Step 9: Final Summary — starting"
mkdir -p "${dir7}" "${dir7}/logs" "${dir7}/checkpoints"
ck_summary="${dir7}/checkpoints/step9.summary.done"
summary_txt="${dir7}/run_summary.txt"
summary_tsv="${dir7}/run_summary.tsv"
summary_json="${dir7}/run_summary.json"
summary_tail="${dir7}/logs/summary_tail.log"

if [[ -f "${ck_summary}" && -s "${summary_txt}" && -s "${summary_tsv}" && -s "${summary_json}" && -s "${summary_tail}" ]]; then
  echo "[INFO] Step 9 checkpoint present — skipping final summary generation"
else
  rm -f "${ck_summary}"

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
# Step 7 writes these six BUSCO fields to its normalized metrics table.
busco_summary() {
  local metrics="${busco_dir}/busco.metrics.tsv"
  printf '%s %s %s %s %s %s\n' \
    "$(assessment_metric "$metrics" 'Complete percentage')" \
    "$(assessment_metric "$metrics" 'Single copy percentage')" \
    "$(assessment_metric "$metrics" 'Multi copy percentage')" \
    "$(assessment_metric "$metrics" 'Fragmented percentage')" \
    "$(assessment_metric "$metrics" 'Missing percentage')" \
    "$(assessment_metric "$metrics" 'n_markers')"
}

# Native assessment metric extractor from metric/value TSV.
assessment_metric() {
  local table="$1" key="$2"
  [[ -s "$table" ]] || { printf '%s\n' '-'; return 0; }
  awk -F'\t' -v k="$key" 'NR>1 && $1==k {gsub(/,/,"",$2); print $2; found=1; exit} END{if(!found) print "-"}' "$table"
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
assessment_dir="${dir7}/assessment"
busco_dir="${dir7}/busco"
btk_dir="${dir7}/blobtoolkit"

: "${species_fasta:=${dir6}/${species_name}.fasta}"      # canonical final assembly
assembly_metrics="${assessment_dir}/${base}.assembly_metrics.tsv"
busco_run_name="${busco_run_name:-final.assembly.busco}"
busco_run_dir_base="${busco_dir}/${busco_run_name}"
busco_run_dir="$(find "${busco_run_dir_base}" -maxdepth 1 -type d -name 'run_*' -print -quit 2>/dev/null || true)"
[[ -n "${busco_run_dir}" ]] || busco_run_dir="${busco_run_dir_base}"
blobdir_path="${btk_dir}/blobdir"

# Final summary requires the canonical published assembly.
final_target="${species_fasta}"
[[ -s "${final_target}" ]] || fatal "Step 9: canonical final assembly missing: ${final_target}"

final_base="$(basename "$final_target")"
final_base="${final_base%.fa}"

# Native assessment metrics.
N50=$(assessment_metric "$assembly_metrics" "n50")
L50=$(assessment_metric "$assembly_metrics" "l50")
CONTIGS=$(assessment_metric "$assembly_metrics" "scaffold_count")
GC=$(assessment_metric "$assembly_metrics" "gc_percent")
Ns100=$(assessment_metric "$assembly_metrics" "ns_per_100kbp")

# BUSCO headline stats
read -r BUSCO_C BUSCO_S BUSCO_D BUSCO_F BUSCO_M BUSCO_n < <(busco_summary "$busco_run_dir")

# Final assembly contig counts + size (human-readable)
read -r final_ctg_n final_ctg_bp final_ctg_size < <(contig_stats "$final_target")

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
repeat_dfam_library="${repeats_dir}/01.known_libraries/Dfam4.taxon_${DFAM_TAXON}.fa"
repeat_repbase_library="${repeats_dir}/01.known_libraries/RepBase.from_Libraries.fa"
repeat_known_library="${repeats_dir}/01.known_libraries/known_repeats.nr.fa"
repeat_rm_library="$(find "${repeats_dir}/02.repeatmodeler" -type f -name '*families.fa' -print -quit 2>/dev/null || true)"
repeat_tetrimmer_library="${repeats_dir}/03.tetrimmer/TEtrimmer_consensus_merged.fasta"
[[ -s "$repeat_tetrimmer_library" ]] || repeat_tetrimmer_library="${repeats_dir}/03.tetrimmer/TEtrimmer_consensus.fasta"
repeat_final_library="$(find "${repeats_dir}/04.final_library" -type f \( -name '*TEtrimmer*Dfam*RepBase*nr*.fa' -o -name '*final*nr*.fa' \) -print -quit 2>/dev/null || true)"
repeat_earlgrey_dir="${repeats_dir}/05.earlgrey_annotation_only"
repeat_species_name="${REPEATS_SPECIES:-${species_name}}"
repeat_earlgrey_rm_dir="${repeat_earlgrey_dir}/${repeat_species_name}_EarlGrey/${repeat_species_name}_RepeatMasker_Against_Custom_Library"
repeat_final_outputs_dir="${repeats_dir}/06.final_outputs"
repeat_non_repeats_bed="${repeats_dir}/07.non_repetitive_regions/non_repetitive_regions.bed"
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
  echo "Native assessment   : N50=${N50:-'-'} | L50=${L50:-'-'} | #contigs=${CONTIGS:-'-'} | GC%=${GC:-'-'} | Ns/100kb=${Ns100:-'-'}"
  echo "BUSCO (${BUSCO_LINEAGE:-'-'}): C=${BUSCO_C:-'-'}% [S=${BUSCO_S:-'-'}% | D=${BUSCO_D:-'-'}%] | F=${BUSCO_F:-'-'}% | M=${BUSCO_M:-'-'}% | n=${BUSCO_n:-'-'}"
  echo "BlobToolKit        : ${blobdir_path} (plots: ${btk_dir}/plots)"
  echo "Repeats            : ${repeat_summary_status} | root=${repeats_dir}"
  echo "Repeat libraries   : known=${repeat_known_library:-'-'} | de_novo=${repeat_rm_library:-'-'} | refined=${repeat_tetrimmer_library:-'-'} | final=${repeat_final_library:-'-'}"
  echo "Repeat annotation  : Earl Grey v${EARLGREY_REQUIRED_VERSION}=${repeat_earlgrey_dir}"
  echo "Internal RM output : ${repeat_earlgrey_rm_dir}"
  echo "Repeat toolchain   : RepeatModeler ${EARLGREY_REQUIRED_REPEATMODELER} | Terrier ${TERRIER_REQUIRED_VERSION} | TEtrimmer ${TETRIMMER_REQUIRED_VERSION} | Earl Grey ${EARLGREY_REQUIRED_VERSION} | RepeatMasker ${EARLGREY_REQUIRED_REPEATMASKER} | Dfam ${EARLGREY_REQUIRED_DFAM}"
  echo "Non-repetitive BED : ${repeat_non_repeats_bed}"
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
  echo "Assembly metrics    : ${assembly_metrics}"
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
  printf "assembly_N50\t%s\n" "$N50"
  printf "assembly_L50\t%s\n" "$L50"
  printf "assembly_scaffolds\t%s\n" "$CONTIGS"
  printf "assembly_GC_percent\t%s\n" "$GC"
  printf "assembly_Ns_per_100kb\t%s\n" "$Ns100"
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
  printf "repeats_earlgrey_repeatmasker_dir\t%s\n" "$repeat_earlgrey_rm_dir"
  printf "repeats_terrier_version\t%s\n" "$TERRIER_REQUIRED_VERSION"
  printf "repeats_tetrimmer_version\t%s\n" "$TETRIMMER_REQUIRED_VERSION"
  printf "repeats_earlgrey_version\t%s\n" "$EARLGREY_REQUIRED_VERSION"
  printf "repeats_repeatmasker_version\t%s\n" "$EARLGREY_REQUIRED_REPEATMASKER"
  printf "repeats_repeatmodeler_version\t%s\n" "$EARLGREY_REQUIRED_REPEATMODELER"
  printf "repeats_dfam_version\t%s\n" "$EARLGREY_REQUIRED_DFAM"
  printf "repeats_dfam_library\t%s\n" "$repeat_dfam_library"
  printf "repeats_repbase_library\t%s\n" "$repeat_repbase_library"
  printf "repeats_earlgrey_dir\t%s\n" "$repeat_earlgrey_dir"
  printf "repeats_final_outputs_dir\t%s\n" "$repeat_final_outputs_dir"
  printf "repeats_non_repetitive_bed\t%s\n" "$repeat_non_repeats_bed"
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
  printf '  "scaffolding": {"best_id": "%s", "k": "%s", "w": "%s", "score": "%s", "N50": "%s", "L50": "%s", "scaffolds": "%s", "Ns_per_100kb": "%s", "merqury_qv": "%s", "merqury_completeness": "%s", "yak_qv": "%s", "scores": "%s", "preselection_scores": "%s", "best_source": "%s"},\n' \
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
  printf '  "assembly_metrics": {"N50": "%s", "L50": "%s", "scaffolds": "%s", "GC_percent": "%s", "Ns_per_100kb": "%s"},\n' \
         "${N50:-}" "${L50:-}" "${CONTIGS:-}" "${GC:-}" "${Ns100:-}"
  printf '  "busco": {"lineage": "%s", "C": "%s", "S": "%s", "D": "%s", "F": "%s", "M": "%s", "n": "%s"},\n' \
         "${BUSCO_LINEAGE:-}" "${BUSCO_C:-}" "${BUSCO_S:-}" "${BUSCO_D:-}" "${BUSCO_F:-}" "${BUSCO_M:-}" "${BUSCO_n:-}"
  printf '  "blobtoolkit": {"blobdir": "%s", "plots_dir": "%s"},\n' \
         "${blobdir_path//\"/\\\"}" "${btk_dir}/plots"
  printf '  "repeats": {"status": "%s", "root": "%s", "dfam_library": "%s", "repbase_library": "%s", "known_library": "%s", "repeatmodeler_library": "%s", "tetrimmer_library": "%s", "final_library": "%s", "earlgrey_dir": "%s", "earlgrey_repeatmasker_dir": "%s", "final_outputs_dir": "%s", "non_repetitive_bed": "%s", "tetrimmer_version": "%s", "earlgrey_version": "%s", "repeatmasker_version": "%s", "repeatmodeler_version": "%s", "dfam_version": "%s", "checkpoints_dir": "%s", "logs_dir": "%s"},\n' \
         "${repeat_summary_status//\"/\\\"}" \
         "${repeats_dir//\"/\\\"}" \
         "${repeat_dfam_library//\"/\\\"}" \
         "${repeat_repbase_library//\"/\\\"}" \
         "${repeat_known_library//\"/\\\"}" \
         "${repeat_rm_library//\"/\\\"}" \
         "${repeat_tetrimmer_library//\"/\\\"}" \
         "${repeat_final_library//\"/\\\"}" \
         "${repeat_earlgrey_dir//\"/\\\"}" \
         "${repeat_earlgrey_rm_dir//\"/\\\"}" \
         "${repeat_final_outputs_dir//\"/\\\"}" \
         "${repeat_non_repeats_bed//\"/\\\"}" \
         "${TETRIMMER_REQUIRED_VERSION//\"/\\\"}" \
         "${EARLGREY_REQUIRED_VERSION//\"/\\\"}" \
         "${EARLGREY_REQUIRED_REPEATMASKER//\"/\\\"}" \
         "${EARLGREY_REQUIRED_REPEATMODELER//\"/\\\"}" \
         "${EARLGREY_REQUIRED_DFAM//\"/\\\"}" \
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
} > "${summary_tail}" || true

# --------------------------- stderr topline ---------------------------
{
  printf "RUN SUMMARY :: Final=%s | Size=%s | N50=%s | L50=%s | BUSCO_C=%s%% | Repeats=%s | Host=%s | Threads=%s | RAM=%sGB | Start=%s | End=%s | Duration=%s\n" \
    "$(basename "$final_target")" "$final_ctg_size" "${N50:-'-'}" "${L50:-'-'}" "${BUSCO_C:-'-'}" \
    "$repeat_summary_status" "$host" "$cpu_threads" "$mem_total_gb" "$start_h" "$end_h" "$runtime_h"
} >&2

[[ -s "${summary_txt}" && -s "${summary_tsv}" && -s "${summary_json}" && -s "${summary_tail}" ]]   || error "Step 9 summary generation incomplete"
touch "${ck_summary}"
fi

echo "[INFO] Pipeline complete. Summary: ${summary_txt} | ${summary_tsv} | ${summary_json}"
