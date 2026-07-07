<img src="genoflowHiFi_badge.png" align="right" width="500"/>

# genoflowHiFi

## PacBio HiFi and optional Hi-C reference-genome assembly workflows

`genoflowHiFi` is a Bash workflow collection for end-to-end vertebrate reference-genome assembly from PacBio HiFi long-read data. It supports two assembly modes:

- HiFi-only assembly and scaffolding with `genoflowHiFi_v0.1.0.sh` or `genoflowHiFi_v0.2.0.sh`.
- HiFi plus Hi-C chromosome-level scaffolding with `genoflowHiFi_HiC_v0.1.0.sh` or `genoflowHiFi_HiC_v0.2.0.sh`.

The pipelines coordinate primary assembly, optional NCBI FCS decontamination, iterative DeepVariant polishing, Inspector correction, scaffolding, final polishing, assembly evaluation, repeat discovery and annotation, and consolidated reporting. They are designed for Linux systems where explicit paths, restart-aware checkpoints, controlled conda environments, Apptainer or Singularity containers, and structured output directories are required.

The repository contains:

```text
genoflowHiFi_v0.1.0.sh          PacBio HiFi assembly pipeline
genoflowHiFi_v0.2.0.sh          PacBio HiFi assembly pipeline with enhanced repeats
genoflowHiFi_HiC_v0.1.0.sh      PacBio HiFi + Hi-C chromosome-level assembly pipeline
genoflowHiFi_HiC_v0.2.0.sh      PacBio HiFi + Hi-C pipeline with enhanced repeats
pacbio_hifi_deepvariant_v0.1.0.sh
install_miniforge3_envs.sh      Miniforge3 and conda environment installer
genoflowHiFi_badge.png          Repository badge
containers/                     Apptainer/Singularity definition files and build notes
README.md                       Repository documentation
REPOSITORY_MAINTENANCE.md       Repository maintenance guide
LICENSE                         GPL-3.0 license
```

## Workflow modes

### HiFi-only mode

Use `genoflowHiFi_v0.1.0.sh` or `genoflowHiFi_v0.2.0.sh` when only PacBio HiFi reads are available, or when long-read scaffolding is preferred.

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

Primary deliverable:

```text
06.final_asm/<Scientific_Name>.fasta
```

### HiFi plus Hi-C mode

Use `genoflowHiFi_HiC_v0.1.0.sh` or `genoflowHiFi_HiC_v0.2.0.sh` when paired Hi-C reads are available and chromosome-level scaffolding is required. These scripts replace the LongStitch step with a YaHS Hi-C scaffolding workflow and add a required manual curation boundary after Step 5.

```text
PacBio HiFi reads + paired Hi-C reads
        |
Hifiasm primary assembly
        |
FCS-Adaptor and FCS-GX decontamination
        |
DeepVariant pre-polishing
        |
Inspector correction
        |
Chromap Hi-C alignment and YaHS scaffolding
        |
Manual YaHS/JBAT curation
        |
TGS-GapCloser2 gap filling
        |
DeepVariant final polishing and optional final FCS cleanup
        |
QUAST, BUSCO, BlobToolKit, Qualimap, and MultiQC evaluation
        |
RepeatModeler2, TEtrimmer, RepeatMasker, and Earl Grey repeat annotation
        |
Final assembly, repeat libraries, statistics, logs, and run summary
```

Primary deliverable:

```text
07.final_asm/<Scientific_Name>.fasta
```

## Hi-C two-pass execution requirement

The Hi-C scripts must be run in two consecutive passes because YaHS/JBAT curation is manual.

Pass 1 runs through Step 5 and creates the YaHS scaffolds and JBAT curation bundle:

```bash
working_dir="/path/to/project" \
pacbio_hifi="/path/to/hifi_reads.fq.gz" \
pacbio_hifi_fasta="/path/to/hifi_reads.fa.gz" \
hic_r1="/path/to/hic_R1.fastq.gz" \
hic_r2="/path/to/hic_R2.fastq.gz" \
TAXID="13604" \
BUSCO_LINEAGE="actinopterygii_odb12" \
run_hifiasm=TRUE \
run_fcs_early=TRUE \
run_deepvariant=TRUE \
run_inspector=TRUE \
run_yahs=TRUE \
run_tgsgapcloser2=FALSE \
run_final_asm=FALSE \
run_stats=FALSE \
run_repeats=FALSE \
bash genoflowHiFi_HiC_v0.2.0.sh 2>&1 | tee -a genoflowHiFi_HiC_v0.2.0.pass1.$(date +%Y%m%d_%H%M%S).log
```

After Pass 1, open the generated YaHS/JBAT files, manually curate the assembly, and save the curated FASTA. By default, the scripts expect:

```text
05.yahs/curated_assembly/<base>.curated.fa
```

Pass 2 resumes from the curated assembly and runs the downstream pipeline:

```bash
working_dir="/path/to/project" \
pacbio_hifi="/path/to/hifi_reads.fq.gz" \
pacbio_hifi_fasta="/path/to/hifi_reads.fa.gz" \
hic_r1="/path/to/hic_R1.fastq.gz" \
hic_r2="/path/to/hic_R2.fastq.gz" \
curated_assembly="/path/to/project/05.yahs/curated_assembly/<base>.curated.fa" \
TAXID="13604" \
BUSCO_LINEAGE="actinopterygii_odb12" \
run_hifiasm=FALSE \
run_fcs_early=FALSE \
run_deepvariant=FALSE \
run_inspector=FALSE \
run_yahs=FALSE \
run_tgsgapcloser2=TRUE \
run_final_asm=TRUE \
run_stats=TRUE \
run_repeats=TRUE \
bash genoflowHiFi_HiC_v0.2.0.sh 2>&1 | tee -a genoflowHiFi_HiC_v0.2.0.pass2.$(date +%Y%m%d_%H%M%S).log
```

## Main workflow stages

### 1. Hifiasm assembly

Runs `hifiasm` on PacBio HiFi reads and produces the primary haplotype-aware draft assembly used by downstream refinement. In the Hi-C scripts, paired Hi-C reads can also be passed to hifiasm through `hic_r1` and `hic_r2`.

### 2. FCS decontamination

Runs optional NCBI FCS-Adaptor and FCS-GX screening to remove sequencing adaptors, foreign contaminants, and taxonomically inconsistent sequence before polishing.

### 3. DeepVariant pre-polishing

Maps HiFi reads back to the draft assembly, calls high-confidence homozygous corrections, filters them, and applies iterative consensus corrections.

### 4. Inspector correction

Runs reference-free long-read assembly evaluation and correction to detect and repair local assembly errors before scaffolding.

### 5. Scaffolding

HiFi-only scripts run a LongStitch parameter grid and rank scaffold candidates using QUAST metrics, GenomeScope2 genome profiling, optional yak screening, Merqury k-mer validation, and robust scoring tables.

Hi-C scripts run Chromap paired-end Hi-C alignment, YaHS scaffolding, and JBAT file generation. This stage intentionally stops at the manual curation boundary before downstream gap closing and polishing.

### 6. Gap closing or final assembly polishing

HiFi-only scripts run post-scaffold DeepVariant polishing, optional final FCS cleanup, and canonical FASTA export with species-derived headers.

Hi-C scripts run TGS-GapCloser2 on the manually curated YaHS assembly, then use the gap-closed assembly as preferred input for final DeepVariant polishing and optional final FCS cleanup.

### 7. Assembly evaluation

Runs QUAST, BUSCO v6, BlobToolKit, Qualimap, and MultiQC to generate structural, completeness, taxonomic, coverage, and mapping-quality summaries.

### 8 or 9. Repeat annotation

Builds and consolidates repeat evidence with RepeatModeler2, Dfam, RepBase, TEtrimmer, RepeatMasker, and Earl Grey annotation-only reporting.

The v0.2 scripts include the enhanced repeat-library contract:

```text
Dfam 4 FamDB archive validation and optional decompression
RepBase input discovery from EMBL, EMBL.GZ, EMBL archives, REF, REF.GZ, or extracted directories
Dfam4_RepBase known-library naming
```

### Final summary

Writes run summaries, key metrics, selected artifacts, software provenance, resource use, and log tails.

## System requirements

`genoflowHiFi` is intended for modern Linux systems with high-capacity storage and enough memory for long-read assembly, polishing, scaffold search, Hi-C scaffolding, and repeat annotation.

Core runtime requirements:

```text
Bash 4+
GNU Parallel
Miniforge3, conda, or mamba
Apptainer or Singularity
Python 3
R
Java for juicer_tools in Hi-C mode
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
Chromap
YaHS
juicer_tools
TGS-GapCloser2
RepeatModeler2
TEtrimmer
RepeatMasker
Earl Grey
NCBI FCS-Adaptor
NCBI FCS-GX
```

The repository includes `install_miniforge3_envs.sh` to create the conda environments used by the workflow. The installer creates the Hi-C environment as `HICenv` with `chromap` and `yahs`; it also creates `LONGSTITCHenv` for the HiFi-only scripts and the repeat/evaluation environments used by both modes. Container definition files are stored under `containers/`.

## Input requirements

Required inputs for all modes:

```text
PacBio HiFi reads in FASTQ, FASTQ.GZ, FASTA, or FASTA.GZ format
NCBI Taxonomy ID or explicit species name
NCBI FCS-GX database when FCS-GX is enabled
BUSCO lineage
```

Additional required inputs for Hi-C mode:

```text
Paired Hi-C reads: hic_r1 and hic_r2
Manual curated YaHS FASTA for the Step 6-10 pass: curated_assembly
juicer_tools JAR for JBAT/HIC generation
```

Optional inputs:

```text
Pre-converted HiFi FASTA for LongStitch or TGS-GapCloser2
PacBio CLR reads for LongStitch scaffolding in HiFi-only mode
Illumina reads or k-mer resources for GenomeScope2, yak, and Merqury validation in HiFi-only mode
scaff_map_tsv for optional final scaffold renaming in Hi-C mode
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
chmod +x genoflowHiFi_v0.1.0.sh genoflowHiFi_v0.2.0.sh
chmod +x genoflowHiFi_HiC_v0.1.0.sh genoflowHiFi_HiC_v0.2.0.sh
chmod +x install_miniforge3_envs.sh
```

Install Miniforge3 and pipeline conda environments:

```bash
bash install_miniforge3_envs.sh
```

Build or provide required containers, then run the selected pipeline with explicit project variables.

### HiFi-only example

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
bash genoflowHiFi_v0.2.0.sh 2>&1 | tee -a genoflowHiFi_v0.2.0.$(date +%Y%m%d_%H%M%S).log
```

### Hi-C chromosome-level example

Use the two-pass command pattern shown in the Hi-C two-pass execution section.

## Configuration

Shared project variables:

```text
working_dir
pacbio_hifi
pacbio_hifi_fasta
TAXID
SPECIES_NAME
SPECIES_HINT
nthreads
njobs
```

HiFi-only variables:

```text
pacbio_clr
run_scaffolding
GENOME_SIZE
GENOME_SIZE_EST
LONGSTITCH_K_VALUES
LONGSTITCH_W_VALUES
LONGSTITCH_ROUNDS
LONGSTITCH_MODE
MERQURY_K
MERQURY_TOP_N
YAK_K
YAK_BLOOM_BITS
```

Hi-C variables:

```text
hic_r1
hic_r2
scaff_map_tsv
curated_assembly
tgsgapcloser_reads
tgsgapcloser2_bin
juicer_tools_jar
run_yahs
run_tgsgapcloser2
ENV_HIC
```

Shared step toggles:

```text
run_hifiasm
run_fcs_early
run_deepvariant
run_inspector
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
ENV_HIC
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

Repeat annotation variables:

```text
REPEATS_PROJECT_ROOT
REPEATS_GENOME
REPEATS_SPECIES
REPEATS_THREADS
REPEATS_LIB_DIR
DFAM_FAMDB_DIR
DFAM_TAXON
DFAM_FLAGS
DFAM_AUTODECOMPRESS_GZ
DFAM_VALIDATE_MD5
DFAM_RELEASE_TAG
REPBASE_DIR
REPBASE_EMBL
REPBASE_EMBL_PATTERN
TETRIMMER_MODE
TETRIMMER_PRESET
RM_ENGINE
EARLGREY_REPEATMASKER_SPECIES
```

Defaults are defined near the top of each pipeline script.

## Output layout

### HiFi-only scripts

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

### Hi-C scripts

```text
01.hifiasm/          Primary hifiasm assembly
02.fcs_decontam/    FCS-Adaptor and FCS-GX decontamination outputs
03.deepvariant/     Pre-Hi-C DeepVariant polishing outputs
04.inspector/       Inspector correction outputs
05.yahs/            Chromap, YaHS, JBAT/HIC files, and curated assembly handoff
06.tgsgapcloser2/   TGS-GapCloser2 gap-closing outputs
07.final_asm/       Final polished and decontaminated assembly
08.asm_stats/       QUAST, BUSCO, BlobToolKit, Qualimap, and MultiQC outputs
09.repeats/         Repeat libraries, RepeatMasker outputs, Earl Grey reports, and BED masks
10.summary/         Final run summary files
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
bash -n genoflowHiFi_v0.2.0.sh
bash -n genoflowHiFi_HiC_v0.1.0.sh
bash -n genoflowHiFi_HiC_v0.2.0.sh
bash -n install_miniforge3_envs.sh
bash -n containers/containers_creation.sh
```

Optional static linting:

```bash
shellcheck -x genoflowHiFi_v0.1.0.sh genoflowHiFi_v0.2.0.sh
shellcheck -x genoflowHiFi_HiC_v0.1.0.sh genoflowHiFi_HiC_v0.2.0.sh
shellcheck install_miniforge3_envs.sh containers/containers_creation.sh
```

Validation scope:

```text
bash -n checks shell syntax only.
shellcheck checks static shell risks.
Neither command proves full runtime correctness on genome assembly data.
```

## License

`genoflowHiFi` is distributed under the GNU General Public License v3.0. See `LICENSE`.
