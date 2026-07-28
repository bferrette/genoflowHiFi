# genoflowHiFi: PacBio HiFi–Hi-C Genome Assembly and Functional Annotation Workflows

## Release information

- Version: `0.1.0`
- Release date: `2026-07-28`
- Version DOI: <https://doi.org/10.5281/zenodo.21650438>
- Concept DOI: <https://doi.org/10.5281/zenodo.21650437>

## Overview

This release provides reproducible Bash workflows for chromosome-level genome
assembly from PacBio HiFi and Hi-C sequencing data and for downstream RNA-seq
processing and functional annotation.

## Included workflows

### `genoflowHiFi_HiC_v0.1.0.sh`

End-to-end genome assembly workflow covering:

- Hifiasm primary assembly.
- Optional FCS-Adaptor and FCS-GX decontamination.
- DeepVariant consensus polishing.
- Inspector assembly evaluation and correction.
- Chromap alignment and YaHS Hi-C scaffolding.
- Manual JBAT curation.
- TGS-GapCloser2 gap closing.
- QUAST, BUSCO, BlobToolKit, Qualimap, and MultiQC evaluation.
- RepeatModeler2, TEtrimmer, RepeatMasker, and Earl Grey repeat analysis.
- Consolidated provenance, artifact, and run-summary reporting.

### `functional_annotation_v2.sh`

Modular RNA-seq and functional-annotation workflow supporting:

- SRA retrieval and FASTQ conversion.
- `fastp` quality control and trimming.
- HISAT2 alignment.
- Gene-level counting with `Rsubread::featureCounts`.
- eggNOG-mapper annotation.
- InterProScan protein-family and domain annotation.
- DIAMOND searches against curated Swiss-Prot proteins.
- KofamScan KEGG Orthology assignment.
- Optional dbCAN CAZyme annotation.
- Consensus annotation tables.
- Descriptive RNA-seq and genome-annotation graphics.

### Supporting R programs

The `functional_annotation_v2_R/` directory contains the R programs required
for feature counting, descriptive RNA-seq analysis, and EGAPx
genome-annotation graphics:

```text
functional_annotation_v2_R/
├── egapx_genome_annotation_graphics_karyoploter.R
├── rnaseq_single_sample_plots.R
└── run_featurecounts_rsubread.R
```

### `install_miniforge3_envs.sh`

When included in the final archive, this optional utility installs Miniforge
and provisions the conda environments used by the workflows.

## Hi-C assembly execution model

The assembly workflow uses a mandatory two-pass design.

1. Run through Step 5 to generate the YaHS scaffolds and JBAT curation files.
2. Manually curate the assembly.
3. Resume Steps 6–10 using the curated FASTA supplied through
   `curated_assembly`.

The default curated-assembly handoff is:

```text
05.yahs/curated_assembly/<base>.curated.fa
```

The principal assembly output is:

```text
07.final_asm/<Scientific_Name>.fasta
```

## Functional-annotation execution profiles

| `RUN_MODE` | Behavior |
|---|---|
| `annotation_only` | Run the five protein functional-annotation methods |
| `all` | Run RNA-seq processing, annotation, and plotting stages |
| `plots_only` | Generate tables and graphics from existing results |
| `preflight` | Validate enabled inputs, databases, environments, and tools |
| `startup_check` | Test initialization with computational stages disabled |
| `custom` | Enable individual stages through explicit `RUN_*` variables |

Expression fold-change plots generated from one sample per biological group are
descriptive. They are not formal differential-expression tests because
biological replication is absent.

## Validation

The released Bash workflows passed:

- Bash syntax validation with `bash -n`.
- Static analysis with ShellCheck.
- Low-impact initialization testing for `functional_annotation_v2.sh`.

The supporting R programs passed R parsing validation.

These checks do not constitute full production validation. Complete execution
requires representative sequencing data, external software, reference
databases, containers, substantial storage, and study-specific configuration.
Local paths, taxonomy information, input files, database locations, and
computational resources must be reviewed before use.

## License

The software is distributed under the GNU General Public License v3.0. See
`LICENSE`.

## Repository

<https://github.com/bferrette/genoflowHiFi>

## Citation

Use the version-specific Zenodo DOI when citing the exact software release used
in an analysis:

```text
[AUTHORS]. (2026). genoflowHiFi: PacBio HiFi–Hi-C Genome Assembly and
Functional Annotation Workflows (version 0.1.0). Zenodo.
https://doi.org/10.5281/zenodo.21650438
```
