<img src="genoflowHiFi_badge.png" align="right" width="250"/>

# genoflowHiFi

## PacBio HiFi reference-genome assembly workflow

`genoflowHiFi` is a Bash workflow for end-to-end vertebrate reference-genome assembly from PacBio HiFi long-read data.

The pipeline coordinates primary assembly, optional NCBI FCS decontamination, iterative DeepVariant polishing, Inspector correction, LongStitch scaffold selection, final polishing, assembly evaluation, repeat discovery and annotation, and consolidated reporting. It is designed for Linux systems where explicit paths, restart-aware checkpoints, controlled conda environments, Apptainer or Singularity containers, and structured output directories are required.

The repository contains:

```text
genoflowHiFi_v0.1.0.sh        Main PacBio HiFi assembly pipeline
install_miniforge3_envs.sh    Miniforge3 and conda environment installer
genoflowHiFi_badge.png        Repository badge
containers/                   Apptainer/Singularity definition files and build notes
README.md                     Repository documentation
REPOSITORY_MAINTENANCE.md     Repository maintenance guide
LICENSE                       GPL-3.0 license
.gitignore                    Output and raw-data exclusion rules
```

## What does `genoflowHiFi` do?

`genoflowHiFi` implements a restart-aware genome assembly workflow:

```text
PacBio HiFi reads
        |
Hifiasm primary assembly
        |
FCS-Adaptor and FCS-GX decontamination
        |
DeepVariant pre-polishing
        |
Inspector correction
        |
LongStitch scaffold grid search and robust selection
        |
DeepVariant final polishing and optional final FCS cleanup
        |
QUAST, BUSCO, BlobToolKit, Qualimap, and MultiQC evaluation
        |
RepeatModeler2, TEtrimmer, RepeatMasker, and Earl Grey repeat annotation
        |
Final assembly, repeat libraries, statistics, logs, and run summary
```

The primary deliverable is:

```text
06.final_asm/<Scientific_Name>.fasta
```

## Main workflow stages

### 1. Hifiasm assembly

Runs `hifiasm` on PacBio HiFi reads and produces the primary haplotype-aware draft assembly used by downstream refinement.

### 2. FCS decontamination

Runs optional NCBI FCS-Adaptor and FCS-GX screening to remove sequencing adaptors, foreign contaminants, and taxonomically inconsistent sequence before polishing.

### 3. DeepVariant pre-polishing

Maps HiFi reads back to the draft assembly, calls high-confidence homozygous corrections, filters them, and applies iterative consensus corrections.

### 4. Inspector correction

Runs reference-free long-read assembly evaluation and correction to detect and repair local assembly errors before scaffolding.

### 5. Scaffolding selection

Runs a LongStitch parameter grid and ranks scaffold candidates using QUAST metrics, GenomeScope2 genome profiling, optional yak screening, Merqury k-mer validation, and robust scoring tables.

### 6. Final assembly polishing

Runs post-scaffold DeepVariant polishing, optional final FCS cleanup, and canonical FASTA export with species-derived headers.

### 7. Assembly evaluation

Runs QUAST, BUSCO v6, BlobToolKit, Qualimap, and MultiQC to generate structural, completeness, taxonomic, coverage, and mapping-quality summaries.

### 8. Repeat annotation

Builds and consolidates repeat evidence with RepeatModeler2, Dfam, RepBase, TEtrimmer, RepeatMasker, and Earl Grey annotation-only reporting.

### 9. Final summary

Writes run summaries, key metrics, selected artifacts, software provenance, resource use, and log tails.

## System requirements

`genoflowHiFi` is intended for modern Linux systems with high-capacity storage and enough memory for long-read assembly, polishing, scaffold search, and repeat annotation.

Core runtime requirements:

```text
Bash 4+
GNU Parallel
Miniforge3, conda, or mamba
Apptainer or Singularity
Python 3
R
pigz
curl or wget
rsync
```

Core workflow tools:

```text
hifiasm
yak
minimap2
SAMtools 1.23.1+
BCFtools 1.23.1+
BEDTools
DeepVariant
GATK
Inspector
LongStitch
GenomeScope2
QUAST
Merqury
BUSCO v6
BlobToolKit
Qualimap
MultiQC
RepeatModeler2
TEtrimmer
RepeatMasker
Earl Grey
NCBI FCS-Adaptor
NCBI FCS-GX
```

The repository includes `install_miniforge3_envs.sh` to create the conda environments used by the workflow. Container definition files are stored under `containers/`.

## Input requirements

Required inputs:

```text
PacBio HiFi reads in FASTQ, FASTQ.GZ, FASTA, or FASTA.GZ format
NCBI Taxonomy ID or explicit species name
NCBI FCS-GX database
BUSCO lineage
```

Optional inputs:

```text
Pre-converted HiFi FASTA for LongStitch scaffolding
PacBio CLR reads for LongStitch scaffolding
Illumina reads or k-mer resources for GenomeScope2, yak, and Merqury validation
Dfam and RepBase repeat-library resources
Prebuilt Apptainer/Singularity containers
```

## Basic usage

Clone the repository:

```bash
git clone https://github.com/bferrette/genoflowHiFi.git
cd genoflowHiFi
```

Make the scripts executable:

```bash
chmod +x genoflowHiFi_v0.1.0.sh install_miniforge3_envs.sh
```

Install Miniforge3 and pipeline conda environments:

```bash
bash install_miniforge3_envs.sh
```

Build or provide required containers, then run the pipeline with explicit project variables:

```bash
working_dir="/path/to/project" \
pacbio_hifi="/path/to/hifi_reads.fq.gz" \
pacbio_hifi_fasta="/path/to/hifi_reads.fasta" \
TAXID="13604" \
BUSCO_LINEAGE="actinopterygii_odb12" \
GX_ROOT="/path/to/ncbi/refseq/fcs-gx/gxdb" \
BLAST_DB_ROOT="/path/to/ncbi/refseq/core_nt/core_nt" \
nthreads=32 \
njobs=4 \
bash genoflowHiFi_v0.1.0.sh 2>&1 | tee -a genoflowHiFi_v0.1.0.$(date +%Y%m%d_%H%M%S).log
```

## Configuration

Project variables:

```text
working_dir
pacbio_hifi
pacbio_hifi_fasta
pacbio_clr
TAXID
SPECIES_NAME
SPECIES_HINT
nthreads
njobs
```

Step toggles:

```text
run_hifiasm
run_fcs_early
run_deepvariant
run_inspector
run_scaffolding
run_final_asm
run_stats
run_repeats
```

Environment variables:

```text
CONDA_SH
MAMBA_SH
ENV_RUNNER
ENV_HIFIASM
ENV_INSPECTOR
ENV_LONGSTITCH
ENV_MERQURY
ENV_GENOMESCOPE2
ENV_R
ENV_QUALIMAP
ENV_QUAST
ENV_BUSCO
ENV_TETRIMMER
ENV_EARLGREY
```

Tool and container variables:

```text
bcftools_bin
bedtools_bin
blastn_bin
minimap2_bin
samtools_bin
sif_blobtoolkit
sif_deepvariant
sif_fcs_adaptor
sif_fcs_gx
sif_gatk4
sif_repeatmasker
ADAPTOR_RUNNER
GX_ROOT
GX_DB
GX_MANIFEST
BLAST_DB_ROOT
BUSCO_LINEAGE
```

Scaffolding and validation variables:

```text
GENOME_SIZE
GENOME_SIZE_EST
LONGSTITCH_K_VALUES
LONGSTITCH_W_VALUES
LONGSTITCH_ROUNDS
LONGSTITCH_MODE
QUAST_MIN_CONTIG
QUAST_KMER_SIZE
MERQURY_K
MERQURY_TOP_N
YAK_K
YAK_BLOOM_BITS
```

Repeat annotation variables:

```text
REPEATS_PROJECT_ROOT
REPEATS_GENOME
REPEATS_SPECIES
REPEATS_THREADS
DFAM_FAMDB_DIR
DFAM_TAXON
REPBASE_EMBL
TETRIMMER_MODE
TETRIMMER_PRESET
RM_ENGINE
EARLGREY_REPEATMASKER_SPECIES
```

Defaults are defined near the top of `genoflowHiFi_v0.1.0.sh`.

## Output layout

```text
01.hifiasm/          Primary hifiasm assembly
02.fcs_decontam/    FCS-Adaptor and FCS-GX decontamination outputs
03.deepvariant/     Pre-scaffold DeepVariant polishing outputs
04.inspector/       Inspector correction outputs
05.scaffolding/     LongStitch candidates, scoring tables, and selected scaffold set
06.final_asm/       Final polished and decontaminated assembly
07.asm_stats/       QUAST, BUSCO, BlobToolKit, Qualimap, MultiQC, and run summaries
08.repeats/         Repeat libraries, RepeatMasker outputs, Earl Grey reports, and BED masks
09.summary/         Reserved summary directory
```

## Containers

Container support files are stored in `containers/`:

```text
blobtoolkit.def          BlobToolKit image definition
gatk4.def                GATK 4 image definition
tetools.def              RepeatMasker and RepeatModeler image definition
containers_creation.sh   Reference build commands for additional images
```

Build commands depend on the local container runtime and writable software prefix. Review each definition file before building.

## Validation

Run shell syntax checks before committing changes:

```bash
bash -n genoflowHiFi_v0.1.0.sh
bash -n install_miniforge3_envs.sh
bash -n containers/containers_creation.sh
```

Optional static linting:

```bash
shellcheck genoflowHiFi_v0.1.0.sh install_miniforge3_envs.sh containers/containers_creation.sh
```

Validation scope:

```text
bash -n checks shell syntax only.
shellcheck checks static shell risks.
Neither command proves full runtime correctness on genome assembly data.
```

## License

`genoflowHiFi` is distributed under the GNU General Public License v3.0. See `LICENSE`.
