# genoflowHiFi Public-Release Maintenance

## Publication boundary

The public repository may contain only these Bash workflows:

```text
genoflowHiFi_HiC_v0.1.0.sh
functional_annotation_v2.sh
```

The R files under `functional_annotation_v2_R/` are required supporting source
files for `functional_annotation_v2.sh`. Documentation, licensing metadata, the
badge, and non-Bash container definitions may also be public.

All other Bash files in the local working directory are private. Exact ignore
rules in `.gitignore` reduce accidental staging, but ignore rules do not conceal
files that are already tracked or present in Git history.

## Required clean-history publication

Do not change the existing GitHub repository to public until its reachable
history has been checked. This checkout has historical commits containing
private Bash workflows. Removing those files in a later commit is insufficient:
older commits remain downloadable.

Use one of these publication models:

1. Preferred: create a new public repository from a clean export containing
   only the approved release files.
2. Rewrite every branch and tag in the existing private repository, verify the
   rewritten object graph, then force-push. This is appropriate only when
   retaining the repository identity is essential.

Do not add the existing private repository as a remote of the clean public
export. Do not push historical tags from the private repository.

## Intended public file set

```text
.gitignore
LICENSE
README.md
REPOSITORY_MAINTENANCE.md
genoflowHiFi_badge.png
genoflowHiFi_HiC_v0.1.0.sh
functional_annotation_v2.sh
functional_annotation_v2_R/egapx_genome_annotation_graphics_karyoploter.R
functional_annotation_v2_R/rnaseq_single_sample_plots.R
functional_annotation_v2_R/run_featurecounts_rsubread.R
containers/blobtoolkit.def
containers/gatk4.def
containers/tetools.def
```

Remove container definitions from the export if they are not required by the
published workflow contract.

## Pre-release checks

Run from the clean public export:

```bash
bash -n genoflowHiFi_HiC_v0.1.0.sh
bash -n functional_annotation_v2.sh
shellcheck -x genoflowHiFi_HiC_v0.1.0.sh
shellcheck -x functional_annotation_v2.sh

project_dir="$(mktemp -d)"
PROJECT_DIR="${project_dir}" RUN_MODE=startup_check \
  bash functional_annotation_v2.sh 1 1

git status --short
git ls-files
git log --all --name-only --pretty=format: | sort -u
git rev-list --all |
  xargs -r git grep -nE \
    'ghp_|github_pat_|AKIA|BEGIN .* PRIVATE KEY|password[[:space:]]*='
```

Expected Bash files from both the index and all reachable history:

```text
functional_annotation_v2.sh
genoflowHiFi_HiC_v0.1.0.sh
```

Verify file modes:

```text
755  public Bash and R scripts
644  README, license, metadata, images, and container definitions
```

## Zenodo release

Zenodo's GitHub integration archives a GitHub release. Therefore, every file in
the release repository is part of the archived software record.

1. Make the clean repository public.
2. Add valid `CITATION.cff` or `.zenodo.json` metadata using the final author,
   ORCID, affiliation, title, version, keywords, and related publication data.
3. Enable the public repository in Zenodo.
4. Create a versioned GitHub release.
5. Wait for Zenodo ingestion and verify the archived file list and metadata.
6. Add the resulting version DOI to `README.md` and cite that version DOI in the
   scientific publication.

A single repository release produces one software record covering both public
workflows. Use separate Zenodo deposits if each script requires an independent
DOI.
