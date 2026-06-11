# genoflowHiFi Repository Maintenance Guide

This guide covers configuration, validation, commit workflow, release workflow, and access control for:

```text
genoflowHiFi -> Bash workflow for PacBio HiFi reference-genome assembly
```

Repository root on the current workstation:

```bash
/datadisk1/research/scripts/genoflow/genoflowHiFi
```

GitHub repository:

```text
https://github.com/bferrette/genoflowHiFi
```

## 1. Repository Role

`genoflowHiFi` contains the PacBio HiFi genome-assembly workflow and supporting environment/container setup files.

Tracked files should include:

```text
genoflowHiFi_v0.1.0.sh
install_miniforge3_envs.sh
genoflowHiFi_badge.png
containers/blobtoolkit.def
containers/gatk4.def
containers/tetools.def
containers/containers_creation.sh
README.md
REPOSITORY_MAINTENANCE.md
LICENSE
.gitignore
```

The repository provides:

```text
Hifiasm primary assembly
NCBI FCS-Adaptor and FCS-GX decontamination
DeepVariant iterative polishing
Inspector correction
LongStitch scaffold-grid search
GenomeScope2, QUAST, yak, and Merqury candidate evaluation
Final polished assembly export
QUAST, BUSCO, BlobToolKit, Qualimap, and MultiQC evaluation
RepeatModeler2, TEtrimmer, RepeatMasker, and Earl Grey repeat annotation
Consolidated run summaries and artifact manifests
```

## 2. Files Not To Commit

Do not commit generated pipeline outputs, logs, caches, temporary directories, raw sequencing data, BAM files, BCF files, VCF files, FASTA assemblies, FASTQ reads, or built container images.

Excluded by `.gitignore`:

```text
*.log
*.tmp
tmp.*/
01.hifiasm/
02.fcs_decontam/
03.deepvariant/
04.inspector/
05.scaffolding/
06.final_asm/
07.asm_stats/
08.repeats/
09.summary/
*.fq
*.fq.gz
*.fastq
*.fastq.gz
*.fa
*.fa.gz
*.fasta
*.fasta.gz
*.fai
*.gfa
*.bam
*.bam.bai
*.cram
*.cram.crai
*.sam
*.bcf
*.bcf.csi
*.vcf
*.vcf.gz
*.vcf.gz.tbi
*.sif
*.img
*.sqsh
.DS_Store
__pycache__/
*.py[cod]
```

Before committing, verify that only source files, documentation, badge files, and repository metadata are staged.

## 3. Check Repository Status

```bash
cd /datadisk1/research/scripts/genoflow/genoflowHiFi
git status --short --branch
git log --oneline --max-count=5
git remote -v
```

Expected remote:

```text
origin -> git@github.com:bferrette/genoflowHiFi.git
```

If the repository has not been initialized:

```bash
git init
git branch -M main
git remote add origin git@github.com:bferrette/genoflowHiFi.git
```

If the remote is missing or uses HTTPS, set the SSH remote:

```bash
git remote set-url origin git@github.com:bferrette/genoflowHiFi.git
```

## 4. SSH Configuration

The repository can use the default GitHub SSH host:

```bash
cd /datadisk1/research/scripts/genoflow/genoflowHiFi
git remote set-url origin git@github.com:bferrette/genoflowHiFi.git
```

Test access:

```bash
ssh -T git@github.com
```

Expected output:

```text
Hi bferrette! You've successfully authenticated, but GitHub does not provide shell access.
```

## 5. Validate Before Commit

Run Bash syntax validation:

```bash
cd /datadisk1/research/scripts/genoflow/genoflowHiFi
bash -n genoflowHiFi_v0.1.0.sh
bash -n install_miniforge3_envs.sh
bash -n containers/containers_creation.sh
```

Optional shell linting:

```bash
shellcheck genoflowHiFi_v0.1.0.sh
shellcheck install_miniforge3_envs.sh
shellcheck containers/containers_creation.sh
```

Optional container-definition inspection:

```bash
apptainer build --dry-run /tmp/blobtoolkit.sif containers/blobtoolkit.def
apptainer build --dry-run /tmp/gatk4.sif containers/gatk4.def
apptainer build --dry-run /tmp/tetools.sif containers/tetools.def
```

Validation scope:

```text
bash -n checks shell syntax only.
shellcheck checks static shell risks.
Apptainer dry-run checks definition parsing only.
None of these commands proves full runtime correctness on genome assembly data.
```

## 6. Commit Updates

Use targeted staging.

Documentation-only update:

```bash
cd /datadisk1/research/scripts/genoflow/genoflowHiFi
git status --short
git add README.md REPOSITORY_MAINTENANCE.md
git commit -m "Update genoflowHiFi documentation"
```

Pipeline script update:

```bash
cd /datadisk1/research/scripts/genoflow/genoflowHiFi
bash -n genoflowHiFi_v0.1.0.sh
bash -n install_miniforge3_envs.sh
git status --short
git add genoflowHiFi_v0.1.0.sh install_miniforge3_envs.sh README.md REPOSITORY_MAINTENANCE.md
git commit -m "Update genoflowHiFi pipeline"
```

Repository metadata update:

```bash
cd /datadisk1/research/scripts/genoflow/genoflowHiFi
git status --short
git add .gitignore LICENSE README.md REPOSITORY_MAINTENANCE.md genoflowHiFi_badge.png containers
git commit -m "Add genoflowHiFi repository metadata"
```

Initial repository population:

```bash
cd /datadisk1/research/scripts/genoflow/genoflowHiFi
bash -n genoflowHiFi_v0.1.0.sh
bash -n install_miniforge3_envs.sh
bash -n containers/containers_creation.sh
git status --short
git add .gitignore LICENSE README.md REPOSITORY_MAINTENANCE.md genoflowHiFi_badge.png genoflowHiFi_v0.1.0.sh install_miniforge3_envs.sh containers
git commit -m "Initial genoflowHiFi repository population"
```

## 7. Pull Before Push

Fetch and integrate remote changes before pushing:

```bash
cd /datadisk1/research/scripts/genoflow/genoflowHiFi
git fetch origin
git pull --rebase origin main
git push
```

If conflicts occur:

```bash
git status
# Edit conflicted files.
git add <resolved-files>
git rebase --continue
git push
```

Abort a bad rebase:

```bash
git rebase --abort
```

## 8. Release Checklist

Before tagging a release:

```bash
cd /datadisk1/research/scripts/genoflow/genoflowHiFi
git status --short
bash -n genoflowHiFi_v0.1.0.sh
bash -n install_miniforge3_envs.sh
bash -n containers/containers_creation.sh
git log --oneline --max-count=5
```

Confirm:

```text
README usage matches the current script variables.
REPOSITORY_MAINTENANCE.md lists all tracked repository files.
.gitignore excludes generated assembly outputs and raw sequencing data.
The badge image is present.
Container definition files are current.
No generated FASTA, FASTQ, BAM, BCF, VCF, SIF, log, cache, or output directory is staged.
```

Tag after validation:

```bash
git tag -a v0.1.0 -m "genoflowHiFi v0.1.0"
git push origin main
git push origin v0.1.0
```

## 9. Access Control

Use SSH remotes for write access:

```text
git@github.com:bferrette/genoflowHiFi.git
```

Use HTTPS only for read-only clone examples:

```text
https://github.com/bferrette/genoflowHiFi.git
```
