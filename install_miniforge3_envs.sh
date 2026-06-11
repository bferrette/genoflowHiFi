#!/usr/bin/env bash
set -euo pipefail

################################################################################
# install_miniforge3_envs.sh
#
# Purpose:
#   Download and install Miniforge3, then create all genome-assembly and
#   repeat-annotation conda environments required by the pipeline.
#
# Default install prefix:
#   $HOME/miniforge3
#
# Usage:
#   bash install_miniforge3_envs.sh
#
# Optional:
#   MINIFORGE_PREFIX="$HOME/software/miniforge3" bash install_miniforge3_envs.sh
#   FORCE_RECREATE=1 bash install_miniforge3_envs.sh
#
# Notes:
#   - Miniforge3 is preferred over Anaconda here because it is conda-forge-native.
#   - Existing environments are skipped unless FORCE_RECREATE=1 is set.
#   - Channels are restricted to conda-forge + bioconda.
#   - Base is updated and receives compiler/compression/parallel tools before
#     any pipeline-specific environments are created.
#   - Non-interactive shells are configured explicitly by using full Miniforge3 paths,
#     conda init, PATH export, conda.sh, and mamba.sh.
################################################################################

#######################################
# 0. USER SETTINGS
#######################################

MINIFORGE_PREFIX="${MINIFORGE_PREFIX:-$HOME/miniforge3}"
MINIFORGE_BIN_DIR="${MINIFORGE_PREFIX}/bin"
MINIFORGE_PROFILE_DIR="${MINIFORGE_PREFIX}/etc/profile.d"
export MAMBA_EXE="${MINIFORGE_BIN_DIR}/mamba"
export MAMBA_ROOT_PREFIX="${MINIFORGE_PREFIX}"
FORCE_RECREATE="${FORCE_RECREATE:-0}"

# Run conda init for future interactive shells.
# Set RUN_CONDA_INIT=0 to skip this step.
RUN_CONDA_INIT="${RUN_CONDA_INIT:-1}"

CHANNEL_ARGS=(
  -c conda-forge
  -c bioconda
  --override-channels
  --strict-channel-priority
)

#######################################
# 1. LOGGING + CHECKS
#######################################

log() {
  printf '[INFO] %s\n' "$*"
}

warn() {
  printf '[WARN] %s\n' "$*" >&2
}

die() {
  printf '[ERROR] %s\n' "$*" >&2
  exit 1
}

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"
}

download_file() {
  local url="$1"
  local output="$2"

  if command -v curl >/dev/null 2>&1; then
    curl -C - -fL --retry 5 --retry-delay 5 --connect-timeout 30 -o "$output" "$url"
  elif command -v wget >/dev/null 2>&1; then
    wget --tries=5 --timeout=30 -O "$output" "$url"
  else
    die "Neither curl nor wget is available."
  fi
}

run_mamba() {
  export MAMBA_EXE="${MINIFORGE_BIN_DIR}/mamba"
  export MAMBA_ROOT_PREFIX="${MINIFORGE_PREFIX}"
  "${MINIFORGE_BIN_DIR}/mamba" --root-prefix "${MINIFORGE_PREFIX}" "$@"
}

#######################################
# 2. DETECT PLATFORM
#######################################

detect_miniforge_installer() {
  local os_name
  local arch_name
  local platform
  local arch

  os_name="$(uname -s)"
  arch_name="$(uname -m)"

  case "$os_name" in
    Linux)  platform="Linux" ;;
    Darwin) platform="MacOSX" ;;
    *) die "Unsupported OS: $os_name" ;;
  esac

  case "$arch_name" in
    x86_64|amd64) arch="x86_64" ;;
    aarch64|arm64) arch="aarch64" ;;
    ppc64le) arch="ppc64le" ;;
    *) die "Unsupported architecture: $arch_name" ;;
  esac

  # Miniforge uses MacOSX-arm64 rather than MacOSX-aarch64.
  if [[ "$platform" == "MacOSX" && "$arch" == "aarch64" ]]; then
    arch="arm64"
  fi

  printf 'Miniforge3-%s-%s.sh' "$platform" "$arch"
}

#######################################
# 3. INSTALL MINIFORGE3
#######################################

install_miniforge() {
  local installer_name
  local installer_url
  local installer_path

  if [[ -x "${MINIFORGE_PREFIX}/bin/conda" ]]; then
    log "Miniforge3 already exists at: ${MINIFORGE_PREFIX}"
    return 0
  fi

  installer_name="$(detect_miniforge_installer)"
  installer_url="https://github.com/conda-forge/miniforge/releases/latest/download/${installer_name}"
  installer_path="${TMPDIR:-/tmp}/${installer_name}"

  log "Installing Miniforge3 to: ${MINIFORGE_PREFIX}"
  log "Downloading installer: ${installer_url}"

  download_file "$installer_url" "$installer_path"

  bash "$installer_path" -b -p "$MINIFORGE_PREFIX"
  rm -f "$installer_path"

  log "Miniforge3 installation completed."
}

#######################################
# 4. LOAD CONDA/MAMBA
#######################################

load_conda() {
  # Miniforge3 command-line tools are in:
  #   ${MINIFORGE_PREFIX}/bin
  #
  # Conda/mamba shell initialization scripts are in:
  #   ${MINIFORGE_PROFILE_DIR}
  #
  # Both are added to PATH because non-interactive scripts do not reliably
  # inherit PATH settings from ~/.bashrc, ~/.profile, or conda init.
  export MAMBA_EXE="${MINIFORGE_BIN_DIR}/mamba"
  export MAMBA_ROOT_PREFIX="${MINIFORGE_PREFIX}"
  export PATH="${MINIFORGE_BIN_DIR}:${MINIFORGE_PREFIX}/condabin:${MINIFORGE_PROFILE_DIR}:${PATH}"

  [[ -f "${MINIFORGE_PROFILE_DIR}/conda.sh" ]] \
    || die "conda.sh not found in ${MINIFORGE_PROFILE_DIR}"

  # shellcheck source=/dev/null
  source "${MINIFORGE_PROFILE_DIR}/conda.sh"

  if [[ -f "${MINIFORGE_PROFILE_DIR}/mamba.sh" ]]; then
    # shellcheck source=/dev/null
    source "${MINIFORGE_PROFILE_DIR}/mamba.sh"
  else
    warn "mamba.sh not found in ${MINIFORGE_PROFILE_DIR}; continuing with conda.sh only."
  fi

  command -v conda >/dev/null 2>&1 \
    || die "conda is not available. Expected executable: ${MINIFORGE_BIN_DIR}/conda"

  conda activate base

  if ! command -v mamba >/dev/null 2>&1; then
    warn "mamba is not available on PATH after loading Miniforge3."
    warn "Installing mamba into the base environment."

    conda install -n base -c conda-forge --override-channels --strict-channel-priority -y mamba

    export MAMBA_EXE="${MINIFORGE_BIN_DIR}/mamba"
  export MAMBA_ROOT_PREFIX="${MINIFORGE_PREFIX}"
  export PATH="${MINIFORGE_BIN_DIR}:${MINIFORGE_PREFIX}/condabin:${MINIFORGE_PROFILE_DIR}:${PATH}"
    hash -r

    if [[ -f "${MINIFORGE_PROFILE_DIR}/mamba.sh" ]]; then
      # shellcheck source=/dev/null
      source "${MINIFORGE_PROFILE_DIR}/mamba.sh"
    fi
  fi

  command -v mamba >/dev/null 2>&1 \
    || die "mamba is not available. Expected executable: ${MINIFORGE_BIN_DIR}/mamba"

  log "Using Miniforge3 prefix: ${MINIFORGE_PREFIX}"
  log "Using profile directory: ${MINIFORGE_PROFILE_DIR}"
  log "Using conda: $(command -v conda)"
  log "Using mamba: $(command -v mamba)"
}


configure_shell_initialization() {
  # Kept intentionally minimal. Persistent shell setup is performed after
  # environment creation and cache cleanup by initialize_mamba_shell.
  :
}

initialize_mamba_shell() {
  log "Initializing mamba shell support for future Bash sessions."

  [[ -x "${MINIFORGE_BIN_DIR}/mamba" ]] \
    || die "mamba executable not found at ${MINIFORGE_BIN_DIR}/mamba"

  # Use the requested simple command:
  #   mamba shell init --shell bash --root-prefix=$HOME/miniforge3
  # The root prefix must be the Miniforge3 installation directory,
  # not the mamba executable path.
  "${MINIFORGE_BIN_DIR}/mamba" shell init --shell bash --root-prefix="${MINIFORGE_PREFIX}"

  log "mamba shell initialization completed."
  log "Restart the shell or run: source ~/.bashrc"
}

validate_miniforge_root_prefix() {
  log "Validating conda/mamba root-prefix configuration."

  export MAMBA_EXE="${MINIFORGE_BIN_DIR}/mamba"
  export MAMBA_ROOT_PREFIX="${MINIFORGE_PREFIX}"
  export PATH="${MINIFORGE_BIN_DIR}:${MINIFORGE_PREFIX}/condabin:${PATH}"

  [[ "${MAMBA_ROOT_PREFIX}" == "${MINIFORGE_PREFIX}" ]] || \
    die "MAMBA_ROOT_PREFIX mismatch: ${MAMBA_ROOT_PREFIX} != ${MINIFORGE_PREFIX}"

  if [[ -x "${MINIFORGE_BIN_DIR}/mamba" ]]; then
    "${MINIFORGE_BIN_DIR}/mamba" --root-prefix "${MINIFORGE_PREFIX}" env list >/dev/null || \
      die "mamba cannot list environments using root prefix: ${MINIFORGE_PREFIX}"
  fi

  "${MINIFORGE_BIN_DIR}/conda" info --base | grep -Fxq "${MINIFORGE_PREFIX}" || \
    warn "conda info --base differs from MINIFORGE_PREFIX; explicit paths will still be used."
}


#######################################
# 5. CONFIGURE CONDA/MAMBA
#######################################

configure_conda() {
  log "Configuring channels and solver."

  conda config --system --remove-key channels >/dev/null 2>&1 || true
  conda config --system --add channels conda-forge
  conda config --system --add channels bioconda
  conda config --system --set channel_priority strict
  conda config --system --set auto_activate false

  # New Miniforge installations generally use libmamba by default, but this is
  # harmless when supported and ignored if unavailable.
  conda config --system --set solver libmamba >/dev/null 2>&1 || true
}

#######################################
# 6. ENVIRONMENT HELPERS
#######################################

env_exists() {
  local env_name="$1"
  conda env list | awk '{print $1}' | grep -Fxq "$env_name"
}

remove_env_if_needed() {
  local env_name="$1"

  if env_exists "$env_name" && [[ "$FORCE_RECREATE" == "1" ]]; then
    log "Removing existing environment: ${env_name}"
    run_mamba env remove -n "$env_name" -y
  fi
}

create_env() {
  local env_name="$1"
  shift

  remove_env_if_needed "$env_name"

  if env_exists "$env_name"; then
    log "Skipping existing environment: ${env_name}"
    return 0
  fi

  log "Creating environment: ${env_name}"
  run_mamba create -n "$env_name" "${CHANNEL_ARGS[@]}" -y "$@"
}

update_base_env() {
  log "Updating Miniforge3 base environment."

  run_mamba update -n base "${CHANNEL_ARGS[@]}" -y --all
}

install_base_runtime_tools() {
  log "Installing compiler, compression, and parallel runtime tools into base."

  run_mamba install -n base "${CHANNEL_ARGS[@]}" -y \
  compilers \
  gcc_linux-64 \
  gxx_linux-64 \
  gfortran_linux-64 \
  gcc_linux-ppc64le \
  gxx_linux-ppc64le \
  gfortran_linux-ppc64le \
  autoconf \
  automake \
  libtool \
  pkg-config \
  make \
  cmake \
  binutils \
  zlib \
  bzip2 \
  xz \
  lz4-c \
  zstd \
  libcurl \
  libdeflate \
  ncurses \
  gsl \
  git \
  curl \
  wget \
  rsync \
  lftp \
  aria2 \
  pigz \
  parallel \
  perl \
  openjdk \
  bc \
  openmpi \
  dos2unix \
  tree
}

prepare_base_env() {
  update_base_env
  install_base_runtime_tools
}

#######################################
# 7. CREATE PIPELINE ENVIRONMENTS
#######################################

create_pipeline_envs() {
  create_env HIFIASMenv       hifiasm yak
  create_env NEXTPOLISH2env   nextpolish2
  create_env INSPECTORenv     inspector
  create_env LONGSTITCHenv    longstitch
  create_env QUASTenv         quast
  create_env BUSCOenv         busco=6.0.0 sepp=4.5.5
  create_env GEMOMAenv        gemoma=1.9
  create_env HICenv           chromap yahs
  create_env TETOOLSenv       repeatmasker=4.2.3 repeatmodeler=2.0.8
  create_env TETRIMMERenv     tetrimmer
  create_env EARLGREYenv      earlgrey bc repeatmodeler=2.0.7
  create_env QUALIMAPenv      qualimap
  create_env PAIRTOOLSenv     pairtools cooler
  create_env GETORGANELLEenv  getorganelle
  create_env MERQURYenv       merqury
  create_env GENOMESCOPE2env  genomescope2
}

#######################################
# 8. VALIDATION
#######################################

validate_envs() {
  log "Validating created environments."

  local envs=(
    HIFIASMenv
    NEXTPOLISH2env
    INSPECTORenv
    LONGSTITCHenv
    QUASTenv
    BUSCOenv
    GEMOMAenv
    HICenv
    TETOOLSenv
    TETRIMMERenv
    EARLGREYenv
    QUALIMAPenv
    PAIRTOOLSenv
    GETORGANELLEenv
    MERQURYenv
    GENOMESCOPE2env
  )

  for env_name in "${envs[@]}"; do
    if env_exists "$env_name"; then
      log "Environment available: ${env_name}"
    else
      die "Environment missing after installation: ${env_name}"
    fi
  done
}

#######################################
# 9. CLEANUP
#######################################

clean_mamba_cache() {
  log "Cleaning conda/mamba package cache."

  if ! run_mamba clean --all -y; then
    warn "mamba clean failed; continuing because environment validation already passed."
  fi
}

#######################################
# 10. MAIN
#######################################

main() {
  need_cmd uname
  need_cmd bash
  need_cmd awk
  need_cmd grep

  install_miniforge
  load_conda
  configure_shell_initialization
  validate_miniforge_root_prefix
  configure_conda

  prepare_base_env

  create_pipeline_envs
  validate_envs
  clean_mamba_cache

  initialize_mamba_shell

  cat <<EOF

[INFO] All environments were created successfully.

Miniforge3 prefix:
  ${MINIFORGE_PREFIX}

Miniforge3 profile directory:
  ${MINIFORGE_PROFILE_DIR}

For this script, PATH was set using:
  ${MINIFORGE_BIN_DIR}:${MINIFORGE_PROFILE_DIR}

Mamba shell initialization was configured with:
  ${MINIFORGE_BIN_DIR}/mamba shell init --shell bash --root-prefix="${MINIFORGE_PREFIX}"

Restart the shell or run:
  source ~/.bashrc

To activate an environment:
  mamba activate HIFIASMenv

Conda activation also remains available:
  source "${MINIFORGE_PROFILE_DIR}/conda.sh"
  conda activate HIFIASMenv

To run without shell activation:
  ${MINIFORGE_BIN_DIR}/mamba --root-prefix "${MINIFORGE_PREFIX}" run -n HIFIASMenv <command>

To recreate all environments from scratch:
  FORCE_RECREATE=1 bash $0
EOF
}

main "$@"