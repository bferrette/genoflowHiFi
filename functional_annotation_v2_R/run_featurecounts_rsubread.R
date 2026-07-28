#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Rsubread)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) return(default)
  if (idx == length(args)) stop("Missing value for argument: ", flag)
  args[[idx + 1]]
}

as_logical_arg <- function(x) {
  toupper(as.character(x)) %in% c("TRUE", "T", "1", "YES", "Y")
}

bam_dir        <- get_arg("--bam_dir")
gtf            <- get_arg("--gtf")
out_counts     <- get_arg("--out_counts")
out_matrix     <- get_arg("--out_matrix")
out_summary    <- get_arg("--out_summary", paste0(out_counts, ".summary"))
threads        <- as.integer(get_arg("--threads", "1"))
feature_type   <- get_arg("--feature_type", "exon")
attr_type      <- get_arg("--attr_type", "gene_id")
strandness     <- as.integer(get_arg("--strandness", "0"))
is_paired      <- as_logical_arg(get_arg("--paired", "TRUE"))
count_pairs    <- as_logical_arg(get_arg("--count_read_pairs", "TRUE"))
allow_multi    <- as_logical_arg(get_arg("--allow_multi_overlap", "FALSE"))
count_multi    <- as_logical_arg(get_arg("--count_multimappers", "FALSE"))

if (!dir.exists(bam_dir)) stop("BAM directory not found: ", bam_dir)
if (!file.exists(gtf)) stop("GTF file not found: ", gtf)
if (is.na(threads) || threads < 1) stop("Invalid --threads value")
if (!strandness %in% c(0L, 1L, 2L)) stop("--strandness must be 0, 1, or 2")

bam_files <- list.files(
  bam_dir,
  pattern = "\\.hisat2\\.sorted\\.bam$",
  full.names = TRUE
)

bam_files <- sort(bam_files)

if (length(bam_files) == 0) {
  stop("No sorted BAM files found in: ", bam_dir)
}

dir.create(dirname(out_counts), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_matrix), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_summary), recursive = TRUE, showWarnings = FALSE)

message("Running Rsubread::featureCounts")
message("BAM files: ", length(bam_files))
message("GTF: ", gtf)
message("Threads: ", threads)
message("Feature type: ", feature_type)
message("Attribute: ", attr_type)
message("Strandness: ", strandness)
message("Paired-end: ", is_paired)
message("Count read pairs: ", count_pairs)
message("Count multimappers: ", count_multi)
message("Allow multi-overlap: ", allow_multi)

fc <- Rsubread::featureCounts(
  files = bam_files,
  annot.ext = gtf,
  isGTFAnnotationFile = TRUE,
  GTF.featureType = feature_type,
  GTF.attrType = attr_type,
  useMetaFeatures = TRUE,
  allowMultiOverlap = allow_multi,
  countMultiMappingReads = count_multi,
  strandSpecific = strandness,
  isPairedEnd = is_paired,
  countReadPairs = count_pairs,
  nthreads = threads
)

counts <- as.data.frame(fc$counts, check.names = FALSE)
annotation <- as.data.frame(fc$annotation, check.names = FALSE)
stat <- as.data.frame(fc$stat, check.names = FALSE)

clean_names <- basename(colnames(counts))
clean_names <- sub("\\.hisat2\\.sorted\\.bam$", "", clean_names)
colnames(counts) <- clean_names

matrix_out <- cbind(gene_id = rownames(counts), counts)
write.table(
  matrix_out,
  file = out_matrix,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

full_out <- cbind(annotation, counts)
write.table(
  full_out,
  file = out_counts,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  stat,
  file = out_summary,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Wrote full Rsubread featureCounts table: ", out_counts)
message("Wrote clean count matrix: ", out_matrix)
message("Wrote assignment summary table: ", out_summary)
