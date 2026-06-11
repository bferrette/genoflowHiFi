#!/usr/bin/env bash
set -euo pipefail

# Containers:
# Braker 3
singularity build "${HOME}/software/braker3.sif" docker://teambraker/braker3:latest

# Blobtoolkit
singularity build "${HOME}/software/blobtoolkit.sif" docker://genomehubs/blobtoolkit:latest

# GATK 4
singularity build "${HOME}/software/gatk4.sif" docker://broadinstitute/gatk:latest

# GATK 3
singularity build "${HOME}/software/gatk3.sif" docker://broadinstitute/gatk3:3.8-1

# Google DeepVariant
singularity build "${HOME}/software/deepvariant.sif" docker://google/deepvariant:1.9.0
singularity build "${HOME}/software/deepvariant.sif" docker://google/deepvariant:latest

# Dfam TE Tools Container
singularity build "${HOME}/software/dfam-tetools-latest.sif" docker://dfam/tetools:latest

# MitoHiFi
singularity build "${HOME}/software/mitohifi.sif" docker://ghcr.io/marcelauliano/mitohifi:latest

# Helixer
singularity build "${HOME}/software/helixer.sif" docker://gglyptodon/helixer-docker:latest

# Juicer
singularity build "${HOME}/software/juicer.sif" docker://aidenlab/juicer-kasm:latest

# FCS-adaptor and FCS-GX
singularity build "${HOME}/software/fcs-adaptor.sif" docker://ncbi/fcs-adaptor:latest
singularity build "${HOME}/software/fcs-gx.sif" docker://ncbi/fcs-gx:latest

# polishCLR
singularity build "${HOME}/software/polishclr.sif" docker://csiva2022/polishclr:latest

# NCBI eGAPx
singularity build "${HOME}/software/egapx.sif" docker://ncbi/egapx:latest
singularity build "${HOME}/software/egapx.sif" docker://ncbi/egapx:0.5.0

# Tetrimmer
singularity build "${HOME}/software/tetrimmer.sif" docker://quay.io/biocontainers/tetrimmer:1.7.1--hdfd78af_0
