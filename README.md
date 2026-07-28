<img src="genoflowHiFi_badge.png" align="right" width="420" alt="genoflowHiFi badge"/>

# genoflowHiFi

Reproducible Bash workflows for chromosome-level genome assembly and functional annotation.

## Public release scope

This repository releases exactly two Bash workflows:

| Script | Purpose |
|---|---|
| `genoflowHiFi_HiC_v0.1.0.sh` | PacBio HiFi assembly with Hi-C scaffolding, manual curation, polishing, evaluation, repeat annotation, and reporting |
| `functional_annotation_v2.sh` | RNA-seq processing, protein functional annotation, consensus tables, and descriptive expression and genome-annotation graphics |

Other development Bash workflows are not part of this public distribution.

`functional_annotation_v2_R/` contains the three R helpers used by the public functional-annotation workflow. They are supporting source files, not additional Bash workflows.

## 1. HiFi plus Hi-C genome assembly

`genoflowHiFi_HiC_v0.1.0.sh` implements the following stages:

```text
PacBio HiFi reads + paired Hi-C reads
        |
Hifiasm primary assembly
        |
Optional FCS-Adaptor and FCS-GX decontamination
        |
DeepVariant pre-polishing
        |
Inspector correction
        |
Chromap alignment and YaHS scaffolding
        |
Required manual YaHS/JBAT curation
        |
TGS-GapCloser2 gap filling
        |
Final DeepVariant polishing and optional FCS-GX cleanup
        |
QUAST, BUSCO, BlobToolKit, Qualimap, and MultiQC evaluation
        |
RepeatModeler2, TEtrimmer, RepeatMasker, and Earl Grey
        |
Final assembly and run summary
```

The primary assembly is written to:

```text
07.final_asm/<Scientific_Name>.fasta
```

### Required two-pass execution

The workflow has a hard manual-curation boundary after YaHS Step 5.

Pass 1 creates the YaHS scaffolds and JBAT curation bundle:

```bash
working_dir="/path/to/project" \
pacbio_hifi="/path/to/hifi_reads.fastq.gz" \
pacbio_hifi_fasta="/path/to/hifi_reads.fasta.gz" \
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
bash genoflowHiFi_HiC_v0.1.0.sh
```

Curate the JBAT assembly and save the curated FASTA. The default handoff path is:

```text
05.yahs/curated_assembly/<base>.curated.fa
```

Pass 2 resumes from that curated assembly:

```bash
working_dir="/path/to/project" \
pacbio_hifi="/path/to/hifi_reads.fastq.gz" \
pacbio_hifi_fasta="/path/to/hifi_reads.fasta.gz" \
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
bash genoflowHiFi_HiC_v0.1.0.sh
```

Important configuration variables are defined near the top of the script. Local executable, database, container, taxonomy, and input paths must be reviewed before execution.

### Assembly output layout

```text
01.hifiasm/          Primary assembly
02.fcs_decontam/    FCS screening and cleanup
03.deepvariant/     Pre-scaffold polishing
04.inspector/       Assembly correction
05.yahs/            Hi-C alignment, YaHS scaffolds, and curation handoff
06.tgsgapcloser2/   Gap-closing outputs
07.final_asm/       Final assembly
08.asm_stats/       Assembly evaluation
09.repeats/         Repeat libraries and annotations
10.summary/         Final run summaries
```

## 2. Functional annotation

`functional_annotation_v2.sh` provides modular stages for:

- SRA download and FASTQ conversion.
- `fastp` quality control and trimming.
- HISAT2 alignment and indexed BAM generation.
- Gene-level counting with `Rsubread::featureCounts`.
- eggNOG-mapper annotation.
- InterProScan domain and protein-family annotation.
- DIAMOND searches against curated Swiss-Prot proteins.
- KofamScan KEGG Orthology assignment.
- Optional dbCAN CAZyme annotation.
- Consensus tables and descriptive RNA-seq and genome-annotation graphics.

The expression component supports descriptive fold-change plots for designs with one sample per organ. These outputs are not formal differential-expression tests because biological replication is absent.

### Run modes

| `RUN_MODE` | Behavior |
|---|---|
| `annotation_only` | Runs the five protein-annotation methods; this is the default |
| `all` | Enables all processing, annotation, and plotting stages |
| `plots_only` | Generates the R-based tables and plots from existing inputs |
| `preflight` | Validates all enabled inputs, databases, environments, and tools without running stages |
| `startup_check` | Exercises initialization with all stages disabled |
| `custom` | Uses only explicitly supplied `RUN_*` stage switches |

Basic invocation:

```bash
PROJECT_DIR="/path/to/project" \
EGAPX_OUT_DIR="/path/to/egapx/out" \
RUN_MODE="annotation_only" \
bash functional_annotation_v2.sh 8 4
```

The two positional arguments are threads per job and parallel jobs. Every study-specific path and stage switch can be overridden with environment variables.

The default sample metadata contains the four SRA accessions used in the associated study. Supply `01_metadata/sra_metadata.tsv` or set `RUN_CREATE_METADATA=TRUE` with an edited accession table for another study.

### Required R helpers

The following files must remain beside the Bash script under `functional_annotation_v2_R/`:

```text
run_featurecounts_rsubread.R
rnaseq_single_sample_plots.R
egapx_genome_annotation_graphics_karyoploter.R
```

Alternatively, set `R_SOURCE_DIR`, `FEATURECOUNTS_R_SOURCE`, `RNASEQ_R_SOURCE`, and `GENOME_GRAPHICS_R_SOURCE` to external copies.

### Functional-annotation output layout

```text
00_logs/                    Logs and stage manifests
01_metadata/                SRA and sample metadata
02_sra/                     SRA downloads
03_fastq/                   Raw FASTQ files
04_fastp/                   Trimmed reads and fastp reports
05_multiqc/                 MultiQC reports
06_hisat2/                  HISAT2 index, BAM files, and summaries
07_featurecounts/           Gene-count tables
08_eggnog_mapper/           eggNOG-mapper outputs
09_interproscan/            InterProScan outputs
10_diamond_annotation/      DIAMOND/Swiss-Prot outputs
11_kofamscan/               KofamScan outputs
12_dbcan/                   dbCAN outputs
13_consensus_annotation/    Combined tables, scripts, and plots
```

## Installation and dependencies

Clone the repository and make the two public Bash workflows executable:

```bash
git clone https://github.com/bferrette/genoflowHiFi.git
cd genoflowHiFi
chmod 755 genoflowHiFi_HiC_v0.1.0.sh functional_annotation_v2.sh
```

The workflows coordinate external software; they do not install it. Major assembly dependencies include Hifiasm, Chromap, YaHS, DeepVariant, Inspector, TGS-GapCloser2, QUAST, BUSCO, BlobToolKit, RepeatModeler2, TEtrimmer, RepeatMasker, Earl Grey, GNU Parallel, and Apptainer or Singularity.

Major functional-annotation dependencies include SRA Toolkit, `fastp`, HISAT2, SAMtools, GNU Parallel, MultiQC, eggNOG-mapper, InterProScan, DIAMOND, KofamScan, dbCAN, mamba, and R with the packages declared in `R_REQUIRED_PACKAGES`.

Exact versions, environment names, paths, databases, and optional stages are recorded in the scripts. Full production runs require substantial storage and compute resources.

The optional definitions under `containers/` document BlobToolKit, GATK 4, and
TE-tools container builds. They do not install or build containers automatically.

## Validation

Syntax validation:

```bash
bash -n genoflowHiFi_HiC_v0.1.0.sh
bash -n functional_annotation_v2.sh
```

Static analysis, when ShellCheck is installed:

```bash
shellcheck -x genoflowHiFi_HiC_v0.1.0.sh
shellcheck -x functional_annotation_v2.sh
```

Low-impact initialization check for the annotation workflow:

```bash
project_dir="$(mktemp -d)"
PROJECT_DIR="${project_dir}" RUN_MODE=startup_check \
  bash functional_annotation_v2.sh 1 1
```

These checks do not replace validation on representative sequencing and annotation data.

## Citation and Zenodo

Create a versioned GitHub release only after the public branch contains the intended release files and no private content in its reachable Git history. Enable the repository in Zenodo before creating the GitHub release. Zenodo will archive the release and assign a version DOI.

For a scientific publication, cite the version DOI for the exact archived release. Add the DOI badge and formatted citation here after Zenodo has created the record.

```text
DOI: pending Zenodo release
Citation: pending Zenodo metadata
```

The Zenodo DOI applies to the archived software release containing both public Bash workflows and their supporting files. Separate DOIs require separate Zenodo records.

## License

This software is distributed under the GNU General Public License v3.0. See `LICENSE`.
