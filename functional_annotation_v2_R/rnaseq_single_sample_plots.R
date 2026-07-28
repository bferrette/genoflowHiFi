#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(gt)
  library(gtsummary)
  library(scales)
  library(ggrepel)
})

tidyplots_available <- requireNamespace("tidyplots", quietly = TRUE)

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) return(default)
  args[[idx + 1]]
}

counts_file <- get_arg("--counts")
metadata_file <- get_arg("--metadata")
eggnog_file <- get_arg("--eggnog", "")
interpro_file <- get_arg("--interpro", "")
diamond_file <- get_arg("--diamond", "")
kofam_file <- get_arg("--kofam", "")
dbcan_dir <- get_arg("--dbcan_dir", "")
hisat2_summary_dir <- get_arg("--hisat2_summary_dir", "")
outdir <- get_arg("--outdir")
plotdir <- get_arg("--plotdir")
tabledir <- get_arg("--tabledir")
annotation_dir <- get_arg("--annotation_dir")
run_descriptive_volcano <- get_arg("--run_descriptive_volcano", "TRUE")
run_tidyplots <- get_arg("--run_tidyplots", "TRUE")
denominator <- get_arg("--denominator", "Muscle")
top_n <- as.integer(get_arg("--top_n", "30"))
pseudocount <- as.numeric(get_arg("--pseudocount", "1"))
lfc_threshold <- as.numeric(get_arg("--lfc_threshold", "1"))
plot_dpi <- as.integer(get_arg("--dpi", "600"))
plot_width <- as.numeric(get_arg("--width", "10"))
plot_height <- as.numeric(get_arg("--height", "7"))
plot_formats <- unique(tolower(strsplit(get_arg("--formats", "png,pdf,svg"), ",", fixed = TRUE)[[1]]))
plot_formats <- trimws(plot_formats)
if (!length(plot_formats) || any(!plot_formats %in% c("png", "pdf", "svg"))) {
  stop("--formats must be a comma-separated subset of png,pdf,svg")
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(plotdir, recursive = TRUE, showWarnings = FALSE)
dir.create(tabledir, recursive = TRUE, showWarnings = FALSE)
dir.create(annotation_dir, recursive = TRUE, showWarnings = FALSE)

save_plot_all <- function(plot, filename, width = plot_width, height = plot_height, dpi = plot_dpi) {
  if ("png" %in% plot_formats) {
    ggsave(file.path(plotdir, paste0(filename, ".png")), plot, width = width, height = height, dpi = dpi, bg = "white", limitsize = FALSE)
  }
  if ("pdf" %in% plot_formats) {
    ggsave(file.path(plotdir, paste0(filename, ".pdf")), plot, width = width, height = height, device = cairo_pdf, bg = "white", limitsize = FALSE)
  }
  if ("svg" %in% plot_formats) {
    ggsave(file.path(plotdir, paste0(filename, ".svg")), plot, width = width, height = height, device = "svg", bg = "white", limitsize = FALSE)
  }
}

wrap_plot_text <- function(x, width = 80) {
  stringr::str_wrap(x, width = width)
}

plot_theme <- function(base_size = 12, legend_position = "right") {
  theme_classic(base_size = base_size) +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(face = "bold", lineheight = 1.05, margin = margin(b = 4)),
      plot.subtitle = element_text(size = rel(0.82), lineheight = 1.05, margin = margin(b = 10)),
      plot.margin = margin(t = 18, r = 30, b = 18, l = 18),
      axis.title = element_text(size = rel(0.92)),
      axis.text = element_text(size = rel(0.82)),
      legend.position = legend_position,
      legend.title = element_text(size = rel(0.9)),
      legend.text = element_text(size = rel(0.82))
    )
}

if (run_tidyplots == "TRUE" && !tidyplots_available) {
  message("RUN_TIDYPLOTS=TRUE but tidyplots is not installed; tidyplots outputs will be skipped.")
}

read_counts <- function(path) {
  readr::read_tsv(path, show_col_types = FALSE) %>%
    mutate(gene_id = as.character(gene_id))
}

read_metadata <- function(path) {
  readr::read_tsv(path, show_col_types = FALSE) %>%
    mutate(
      sample = as.character(sample),
      accession = as.character(accession),
      organ = as.character(organ)
    )
}

read_eggnog <- function(path) {
  if (!nzchar(path) || !file.exists(path)) {
    return(tibble(protein_id = character(), eggnog_description = character()))
  }

  raw <- readr::read_tsv(
    path,
    comment = "##",
    show_col_types = FALSE,
    col_names = TRUE,
    col_types = readr::cols(.default = readr::col_character()),
    name_repair = "unique_quiet"
  )

  # eggNOG-mapper files often use "#query" as first column.
  names(raw)[1] <- "protein_id"

  raw %>%
    transmute(
      protein_id = as.character(protein_id),
      eggnog_seed_ortholog = if ("seed_ortholog" %in% names(raw)) seed_ortholog else NA_character_,
      eggnog_description = if ("Description" %in% names(raw)) Description else NA_character_,
      eggnog_preferred_name = if ("Preferred_name" %in% names(raw)) Preferred_name else NA_character_,
      eggnog_gos = if ("GOs" %in% names(raw)) GOs else NA_character_,
      eggnog_kegg_ko = if ("KEGG_ko" %in% names(raw)) KEGG_ko else NA_character_,
      eggnog_kegg_pathway = if ("KEGG_Pathway" %in% names(raw)) KEGG_Pathway else NA_character_,
      eggnog_pfam = if ("PFAMs" %in% names(raw)) PFAMs else NA_character_
    )
}

read_interpro <- function(path) {
  if (!nzchar(path) || !file.exists(path)) {
    return(tibble(protein_id = character()))
  }

  # InterProScan TSV columns:
  # 1 protein accession, 4 analysis, 5 signature accession, 6 signature description,
  # 12 InterPro accession, 13 InterPro description, 14 GO, 15 pathways.
  raw <- readr::read_tsv(
    path,
    col_names = FALSE,
    show_col_types = FALSE,
    col_types = readr::cols(.default = readr::col_character())
  )

  if (ncol(raw) < 13) {
    return(tibble(protein_id = character()))
  }

  raw %>%
    transmute(
      protein_id = as.character(X1),
      analysis = as.character(X4),
      signature = as.character(X5),
      signature_description = as.character(X6),
      interpro_id = as.character(X12),
      interpro_description = as.character(X13),
      go_terms = if (ncol(raw) >= 14) as.character(X14) else NA_character_,
      pathways = if (ncol(raw) >= 15) as.character(X15) else NA_character_
    ) %>%
    group_by(protein_id) %>%
    summarise(
      interpro_ids = paste(unique(na.omit(interpro_id[interpro_id != "-"])), collapse = ";"),
      interpro_descriptions = paste(unique(na.omit(interpro_description[interpro_description != "-"])), collapse = ";"),
      pfam_ids = paste(unique(na.omit(signature[analysis == "Pfam"])), collapse = ";"),
      pfam_descriptions = paste(unique(na.omit(signature_description[analysis == "Pfam"])), collapse = ";"),
      interpro_go_terms = paste(unique(na.omit(go_terms[go_terms != "-"])), collapse = ";"),
      interpro_pathways = paste(unique(na.omit(pathways[pathways != "-"])), collapse = ";"),
      .groups = "drop"
    )
}

read_diamond <- function(path) {
  if (!nzchar(path) || !file.exists(path)) {
    return(tibble(protein_id = character()))
  }

  raw <- readr::read_tsv(
    path,
    comment = "#",
    col_names = FALSE,
    show_col_types = FALSE,
    col_types = readr::cols(.default = readr::col_character())
  )
  if (nrow(raw) == 0) return(tibble(protein_id = character()))

  if (ncol(raw) >= 17) {
    names(raw)[seq_len(17)] <- c(
      "qseqid", "qlen", "sseqid", "slen", "pident", "length", "mismatch", "gapopen",
      "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovhsp", "scovhsp", "stitle"
    )
  }

  needed <- c("qseqid", "sseqid", "pident", "evalue", "bitscore", "qcovhsp", "scovhsp", "stitle")
  for (nm in needed) {
    if (!nm %in% names(raw)) raw[[nm]] <- NA
  }

  raw %>%
    filter(qseqid != "qseqid") %>%
    mutate(
      qseqid = as.character(qseqid),
      sseqid = as.character(sseqid),
      stitle = as.character(stitle),
      pident = suppressWarnings(as.numeric(pident)),
      evalue = suppressWarnings(as.numeric(evalue)),
      bitscore = suppressWarnings(as.numeric(bitscore)),
      qcovhsp = suppressWarnings(as.numeric(qcovhsp)),
      scovhsp = suppressWarnings(as.numeric(scovhsp))
    ) %>%
    arrange(qseqid, evalue, desc(bitscore)) %>%
    group_by(qseqid) %>%
    slice(1) %>%
    ungroup() %>%
    transmute(
      protein_id = qseqid,
      diamond_sseqid = sseqid,
      diamond_pident = pident,
      diamond_evalue = evalue,
      diamond_bitscore = bitscore,
      diamond_qcovhsp = qcovhsp,
      diamond_scovhsp = scovhsp,
      diamond_title = stitle
    )
}

read_kofam <- function(path) {
  if (!nzchar(path) || !file.exists(path)) {
    return(tibble(protein_id = character()))
  }

  raw <- readr::read_tsv(
    path,
    show_col_types = FALSE,
    col_names = FALSE,
    comment = "#",
    col_types = readr::cols(.default = readr::col_character())
  )
  if (nrow(raw) == 0 || ncol(raw) < 2) return(tibble(protein_id = character()))

  marker_like <- ncol(raw) >= 3 && all(as.character(raw[[1]]) %in% c("*", "", " ", "-", NA_character_))

  if (marker_like) {
    protein <- raw[[2]]
    ko <- raw[[3]]
    score <- if (ncol(raw) >= 5) raw[[5]] else NA
    evalue <- if (ncol(raw) >= 6) raw[[6]] else NA
  } else {
    protein <- raw[[1]]
    ko <- raw[[2]]
    score <- if (ncol(raw) >= 4) raw[[4]] else NA
    evalue <- if (ncol(raw) >= 5) raw[[5]] else NA
  }

  tibble(
    protein_id = as.character(protein),
    ko = as.character(ko),
    score = suppressWarnings(as.numeric(score)),
    evalue = suppressWarnings(as.numeric(evalue))
  ) %>%
    filter(!is.na(protein_id), protein_id != "", !is.na(ko), ko != "") %>%
    group_by(protein_id) %>%
    summarise(
      kofam_ko = paste(unique(na.omit(ko)), collapse = ";"),
      kofam_best_score = suppressWarnings(max(score, na.rm = TRUE)),
      kofam_best_evalue = suppressWarnings(min(evalue, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(
      kofam_best_score = ifelse(is.infinite(kofam_best_score), NA_real_, kofam_best_score),
      kofam_best_evalue = ifelse(is.infinite(kofam_best_evalue), NA_real_, kofam_best_evalue)
    )
}

read_dbcan_manifest <- function(path) {
  if (!nzchar(path) || !dir.exists(path)) {
    return(tibble(dbcan_output_dir = character(), dbcan_files = character()))
  }

  files <- list.files(path, recursive = TRUE, full.names = FALSE)
  tibble(
    dbcan_output_dir = path,
    dbcan_files = paste(files, collapse = ";")
  )
}

read_hisat2_alignment_qc <- function(summary_dir) {
  if (!nzchar(summary_dir) || !dir.exists(summary_dir)) {
    return(tibble(sample = character(), category = character(), reads = numeric()))
  }

  files <- list.files(
    summary_dir,
    pattern = "\\.hisat2\\.summary\\.txt$",
    full.names = TRUE
  )

  if (length(files) == 0) {
    return(tibble(sample = character(), category = character(), reads = numeric()))
  }

  parse_num <- function(pattern, x) {
    hit <- stringr::str_match(x, pattern)[, 2]
    hit <- hit[!is.na(hit)]
    if (length(hit) == 0) return(0)
    sum(suppressWarnings(as.numeric(gsub(",", "", hit))), na.rm = TRUE)
  }

  purrr::map_dfr(files, function(f) {
    x <- readLines(f, warn = FALSE)
    sample <- basename(f) %>% sub("\\.hisat2\\.summary\\.txt$", "", .)

    total_pairs <- parse_num("^\\s*Total pairs:\\s*([0-9,]+)", x)
    total_reads <- parse_num("^\\s*Total reads:\\s*([0-9,]+)", x)
    total <- ifelse(total_pairs > 0, total_pairs, total_reads)

    conc_0 <- parse_num("^\\s*Aligned concordantly 0 time:\\s*([0-9,]+)", x)
    conc_1 <- parse_num("^\\s*Aligned concordantly 1 time:\\s*([0-9,]+)", x)
    conc_gt1 <- parse_num("^\\s*Aligned concordantly >1 times:\\s*([0-9,]+)", x)
    disc_1 <- parse_num("^\\s*Aligned discordantly 1 time:\\s*([0-9,]+)", x)

    unpaired_0 <- parse_num("^\\s*Aligned 0 time:\\s*([0-9,]+)", x)
    unpaired_1 <- parse_num("^\\s*Aligned 1 time:\\s*([0-9,]+)", x)
    unpaired_gt1 <- parse_num("^\\s*Aligned >1 times:\\s*([0-9,]+)", x)

    if (total_pairs > 0) {
      tibble(
        sample = sample,
        category = c(
          "Concordant unique",
          "Concordant multi",
          "Discordant unique",
          "Mate-rescue unique",
          "Mate-rescue multi",
          "Unmapped"
        ),
        reads = c(
          conc_1,
          conc_gt1,
          disc_1,
          unpaired_1 / 2,
          unpaired_gt1 / 2,
          unpaired_0 / 2
        )
      )
    } else if (total_reads > 0) {
      tibble(
        sample = sample,
        category = c("Unique", "Multi", "Unmapped"),
        reads = c(unpaired_1, unpaired_gt1, unpaired_0)
      )
    } else {
      tibble(sample = character(), category = character(), reads = numeric())
    }
  }) %>%
    filter(reads > 0)
}

counts <- read_counts(counts_file)
metadata <- read_metadata(metadata_file)

# Align sample IDs in metadata with count matrix column names.
count_samples <- setdiff(names(counts), "gene_id")

metadata <- metadata %>%
  mutate(
    expected_sample = paste(accession, organ, sep = "."),
    sample = if_else(expected_sample %in% count_samples, expected_sample, sample)
  ) %>%
  filter(sample %in% count_samples)

if (nrow(metadata) == 0) {
  stop("No metadata sample IDs match count matrix columns.")
}

counts <- counts %>%
  select(gene_id, all_of(metadata$sample))

count_mat <- counts %>%
  column_to_rownames("gene_id") %>%
  as.matrix()

storage.mode(count_mat) <- "numeric"

# Basic expression summaries ----------------------------------------------------

library_summary <- tibble(
  sample = colnames(count_mat),
  total_counts = colSums(count_mat),
  detected_genes_count_ge_1 = colSums(count_mat >= 1),
  detected_genes_count_ge_10 = colSums(count_mat >= 10)
) %>%
  left_join(metadata, by = "sample") %>%
  arrange(organ)

readr::write_tsv(library_summary, file.path(outdir, "library_summary.tsv"))

library_summary %>%
  select(sample, accession, organ, total_counts, detected_genes_count_ge_1, detected_genes_count_ge_10) %>%
  gt() %>%
  fmt_number(columns = c(total_counts, detected_genes_count_ge_1, detected_genes_count_ge_10), decimals = 0) %>%
  tab_header(
    title = "RNA-seq library summary",
    subtitle = "Descriptive summary from one RNA-seq library per organ"
  ) %>%
  gtsave(file.path(tabledir, "library_summary.html"))

# CPM normalization for descriptive plots.
lib_sizes <- colSums(count_mat)
cpm_mat <- sweep(count_mat, 2, lib_sizes, "/") * 1e6
log2cpm_mat <- log2(cpm_mat + pseudocount)

cpm_tbl <- cpm_mat %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  as_tibble()

log2cpm_tbl <- log2cpm_mat %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  as_tibble()

readr::write_tsv(cpm_tbl, file.path(outdir, "cpm_matrix.tsv"))
readr::write_tsv(log2cpm_tbl, file.path(outdir, "log2cpm_matrix.tsv"))

# Sample Spearman correlation heatmap ------------------------------------------

sample_cor_tbl <- cor(log2cpm_mat, method = "spearman") %>%
  as.data.frame(check.names = FALSE) %>%
  rownames_to_column("sample_x") %>%
  as_tibble() %>%
  pivot_longer(
    -sample_x,
    names_to = "sample_y",
    values_to = "spearman_correlation"
  ) %>%
  mutate(spearman_correlation = replace_na(spearman_correlation, 0))

readr::write_tsv(
  sample_cor_tbl,
  file.path(outdir, "sample_spearman_correlation_matrix_long.tsv")
)

sample_cor_min <- min(sample_cor_tbl$spearman_correlation, na.rm = TRUE)
if (!is.finite(sample_cor_min) || sample_cor_min >= 1) {
  sample_cor_min <- 0.95
}

sample_cor_height <- max(6, 0.55 * n_distinct(sample_cor_tbl$sample_x) + 3)

p_sample_cor <- sample_cor_tbl %>%
  ggplot(aes(x = sample_x, y = sample_y, fill = spearman_correlation)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = number(spearman_correlation, accuracy = 0.01)), size = 2.7) +
  plot_theme(base_size = 10) +
  labs(
    title = wrap_plot_text("Sample Spearman correlation heatmap", 78),
    subtitle = wrap_plot_text("Calculated from log2(CPM + pseudocount) expression profiles.", 92),
    x = NULL,
    y = NULL,
    fill = "Spearman correlation"
  ) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 18)) +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 18)) +
  scale_fill_gradient(
    low = "#440154",
    high = "#FDE725",
    limits = c(sample_cor_min, 1),
    oob = scales::squish
  ) +
  coord_fixed(clip = "off") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot_all(
  p_sample_cor,
  "sample_spearman_correlation_heatmap",
  width = max(7.5, 0.75 * n_distinct(sample_cor_tbl$sample_x) + 4),
  height = sample_cor_height
)

# Optional tidyplots: sample Spearman correlation heatmap -----------------------

if (run_tidyplots == "TRUE" && tidyplots_available) {
  library(tidyplots)

  p_sample_cor_tidyplots <- sample_cor_tbl %>%
    tidyplot(x = sample_x, y = sample_y, color = spearman_correlation) %>%
    add_heatmap() %>%
    adjust_colors(colors_continuous_inferno) %>%
    remove_x_axis_title() %>%
    remove_y_axis_title() %>%
    adjust_x_axis(rotate_labels = 45) %>%
    adjust_legend_title("Spearman correlation") %>%
    adjust_legend_position("right") %>%
    adjust_theme_details(legend.key.height = grid::unit(1, "null"))

  save_plot_all(
    p_sample_cor_tidyplots,
    "sample_spearman_correlation_heatmap_tidyplots",
    width = max(7.5, 0.75 * n_distinct(sample_cor_tbl$sample_x) + 4),
    height = sample_cor_height
  )
}

# Annotation merge --------------------------------------------------------------

eggnog <- read_eggnog(eggnog_file)
interpro <- read_interpro(interpro_file)
diamond <- read_diamond(diamond_file)
kofam <- read_kofam(kofam_file)
dbcan_manifest <- read_dbcan_manifest(dbcan_dir)

# Protein-level merge. This is intentionally kept separate from the gene-level
# count matrix because EGAPx gene IDs and protein IDs are not guaranteed to be
# identical. Build an explicit EGAPx gene/transcript/protein map before using
# this table for final gene-level product naming.
annotation_by_protein <- list(eggnog, interpro, diamond, kofam) %>%
  purrr::reduce(full_join, by = "protein_id")

annotation_by_gene_like_id <- annotation_by_protein %>%
  rename(gene_id = protein_id)

readr::write_tsv(annotation_by_protein, file.path(annotation_dir, "combined_functional_annotation_by_protein_id.tsv"))
readr::write_tsv(dbcan_manifest, file.path(annotation_dir, "dbcan_output_manifest.tsv"))

annotated_log2cpm <- log2cpm_tbl %>%
  left_join(annotation_by_gene_like_id, by = "gene_id")

readr::write_tsv(annotated_log2cpm, file.path(outdir, "annotated_log2cpm_matrix.tsv"))

# Plot 1: library size ----------------------------------------------------------

p_library <- library_summary %>%
  ggplot(aes(x = organ, y = total_counts)) +
  geom_col(width = 0.62, fill = "#4B5563") +
  geom_text(aes(label = comma(total_counts)), vjust = -0.45, size = 3.1) +
  plot_theme(base_size = 12, legend_position = "none") +
  labs(
    title = wrap_plot_text("RNA-seq library size by organ", 78),
    subtitle = wrap_plot_text("One RNA-seq library per organ; values are total assigned featureCounts reads.", 92),
    x = "Organ",
    y = "Total assigned count depth"
  ) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 16)) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.18))) +
  coord_cartesian(clip = "off")

save_plot_all(p_library, "library_size_by_organ")

# Plot 2: detected genes --------------------------------------------------------

detected_long <- library_summary %>%
  select(sample, organ, detected_genes_count_ge_1, detected_genes_count_ge_10) %>%
  pivot_longer(
    starts_with("detected_genes"),
    names_to = "threshold",
    values_to = "detected_genes"
  ) %>%
  mutate(
    threshold = recode(
      threshold,
      detected_genes_count_ge_1 = "Count ≥ 1",
      detected_genes_count_ge_10 = "Count ≥ 10"
    )
  )

p_detected <- detected_long %>%
  ggplot(aes(x = organ, y = detected_genes, fill = threshold)) +
  geom_col(position = position_dodge(width = 0.68), width = 0.58) +
  plot_theme(base_size = 12) +
  labs(
    title = wrap_plot_text("Detected genes by organ", 78),
    subtitle = wrap_plot_text("Detection thresholds based on raw gene-level counts.", 92),
    x = "Organ",
    y = "Detected genes",
    fill = "Threshold"
  ) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 16)) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.12))) +
  scale_fill_manual(values = c("Count ≥ 1" = "#E76F51", "Count ≥ 10" = "#2A9D8F")) +
  coord_cartesian(clip = "off")

save_plot_all(p_detected, "detected_genes_by_organ")

# Plot 3: PCA ------------------------------------------------------------------

# PCA is descriptive here; with four samples it shows global expression structure.
pca <- prcomp(t(log2cpm_mat), center = TRUE, scale. = FALSE)
pca_var <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)

pca_tbl <- as_tibble(pca$x, rownames = "sample") %>%
  left_join(metadata, by = "sample")

readr::write_tsv(pca_tbl, file.path(outdir, "pca_coordinates.tsv"))

p_pca <- pca_tbl %>%
  ggplot(aes(x = PC1, y = PC2, label = organ)) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "grey85") +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "grey85") +
  geom_point(size = 3.5) +
  geom_text_repel(
    size = 3.2,
    box.padding = 0.45,
    point.padding = 0.25,
    min.segment.length = 0,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  plot_theme(base_size = 12, legend_position = "none") +
  labs(
    title = wrap_plot_text("PCA of log2(CPM + pseudocount) expression profiles", 78),
    subtitle = wrap_plot_text("Descriptive PCA; one library per organ.", 92),
    x = paste0("PC1: ", pca_var[1], "% variance"),
    y = paste0("PC2: ", pca_var[2], "% variance")
  ) +
  scale_x_continuous(expand = expansion(mult = 0.12)) +
  scale_y_continuous(expand = expansion(mult = 0.12)) +
  coord_cartesian(clip = "off")

save_plot_all(p_pca, "PCA_log2CPM_expression")

# Optional tidyplots: featureCounts assignment QC -------------------------------

if (run_tidyplots == "TRUE" && tidyplots_available) {
  library(tidyplots)

  fc_summary_file <- paste0(counts_file, ".summary")

  if (file.exists(fc_summary_file)) {
    featurecounts_qc <- readr::read_tsv(
      fc_summary_file,
      show_col_types = FALSE
    ) %>%
      pivot_longer(
        -Status,
        names_to = "sample",
        values_to = "reads"
      ) %>%
      rename(category = Status) %>%
      mutate(reads = suppressWarnings(as.numeric(reads))) %>%
      filter(!is.na(reads), reads > 0) %>%
      group_by(sample) %>%
      mutate(sample_total = sum(reads, na.rm = TRUE)) %>%
      ungroup() %>%
      mutate(sample = fct_reorder(sample, sample_total))

    readr::write_tsv(
      featurecounts_qc,
      file.path(outdir, "featurecounts_assignment_summary_long.tsv")
    )

    featurecounts_qc_height <- max(4.5, 0.5 * n_distinct(featurecounts_qc$sample) + 2)

    p_featurecounts_abs <- featurecounts_qc %>%
      tidyplot(x = reads, y = sample, color = category) %>%
      add_barstack_absolute(reverse = TRUE) %>%
      theme_minimal_x() %>%
      adjust_x_axis(title = "Number of reads", cut_short_scale = TRUE) %>%
      remove_y_axis_title() %>%
      remove_legend_title() %>%
      adjust_legend_position("right")

    save_plot_all(
      p_featurecounts_abs,
      "featurecounts_assignment_absolute_tidyplots",
      width = 8.5,
      height = featurecounts_qc_height
    )

    p_featurecounts_rel <- featurecounts_qc %>%
      tidyplot(x = reads, y = sample, color = category) %>%
      add_barstack_relative(reverse = TRUE) %>%
      theme_minimal_x() %>%
      adjust_x_axis(title = "Percentage of reads", labels = scales::percent) %>%
      remove_y_axis_title() %>%
      remove_legend_title() %>%
      adjust_legend_position("right")

    save_plot_all(
      p_featurecounts_rel,
      "featurecounts_assignment_relative_tidyplots",
      width = 8.5,
      height = featurecounts_qc_height
    )
  }
}

# Optional tidyplots: HISAT2 alignment QC --------------------------------------

if (run_tidyplots == "TRUE" && tidyplots_available) {
  library(tidyplots)

  alignment_qc <- read_hisat2_alignment_qc(hisat2_summary_dir) %>%
    group_by(sample) %>%
    mutate(sample_total = sum(reads, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(sample = fct_reorder(sample, sample_total))

  if (nrow(alignment_qc) > 0) {
    readr::write_tsv(
      alignment_qc,
      file.path(outdir, "hisat2_alignment_summary_long.tsv")
    )

    alignment_qc_height <- max(4.5, 0.5 * n_distinct(alignment_qc$sample) + 2)

    p_alignment_abs <- alignment_qc %>%
      tidyplot(x = reads, y = sample, color = category) %>%
      add_barstack_absolute(reverse = TRUE) %>%
      theme_minimal_x() %>%
      adjust_x_axis(title = "Number of read pairs", cut_short_scale = TRUE) %>%
      remove_y_axis_title() %>%
      remove_legend_title() %>%
      adjust_legend_position("right")

    save_plot_all(
      p_alignment_abs,
      "hisat2_alignment_absolute_tidyplots",
      width = 8.5,
      height = alignment_qc_height
    )

    p_alignment_rel <- alignment_qc %>%
      tidyplot(x = reads, y = sample, color = category) %>%
      add_barstack_relative(reverse = TRUE) %>%
      theme_minimal_x() %>%
      adjust_x_axis(title = "Percentage of read pairs", labels = scales::percent) %>%
      remove_y_axis_title() %>%
      remove_legend_title() %>%
      adjust_legend_position("right")

    save_plot_all(
      p_alignment_rel,
      "hisat2_alignment_relative_tidyplots",
      width = 8.5,
      height = alignment_qc_height
    )
  }
}

# Top expressed genes per organ -------------------------------------------------

long_expr <- log2cpm_tbl %>%
  pivot_longer(-gene_id, names_to = "sample", values_to = "log2cpm") %>%
  left_join(metadata, by = "sample") %>%
  mutate(cpm = 2^log2cpm - pseudocount)

top_expr <- long_expr %>%
  group_by(organ) %>%
  slice_max(order_by = log2cpm, n = top_n, with_ties = FALSE) %>%
  ungroup() %>%
  left_join(annotation_by_gene_like_id, by = "gene_id")

readr::write_tsv(top_expr, file.path(outdir, "top_expressed_genes_per_organ.tsv"))

top_expr %>%
  select(any_of(c("organ", "gene_id", "log2cpm", "cpm", "eggnog_description", "interpro_descriptions", "pfam_ids", "diamond_title", "kofam_ko"))) %>%
  arrange(organ, desc(log2cpm)) %>%
  gt(groupname_col = "organ") %>%
  fmt_number(columns = c(log2cpm, cpm), decimals = 2) %>%
  tab_header(
    title = "Top expressed genes per organ",
    subtitle = "Descriptive ranking by log2(CPM + pseudocount)"
  ) %>%
  gtsave(file.path(tabledir, "top_expressed_genes_per_organ.html"))

# Plot 4: top expressed heatmap-like tile plot ---------------------------------

top_gene_ids <- top_expr %>%
  distinct(gene_id) %>%
  pull(gene_id)

heat_tbl <- long_expr %>%
  filter(gene_id %in% top_gene_ids) %>%
  group_by(gene_id) %>%
  mutate(z_log2cpm = as.numeric(scale(log2cpm))) %>%
  ungroup() %>%
  mutate(
    z_log2cpm = replace_na(z_log2cpm, 0),
    z_plot = pmin(pmax(z_log2cpm, -2.5), 2.5)
  )

heat_gene_count <- n_distinct(heat_tbl$gene_id)
heatmap_show_gene_labels <- heat_gene_count <= 50
heatmap_height <- max(6.5, min(18, 0.13 * heat_gene_count + 3))

p_heat <- heat_tbl %>%
  ggplot(aes(x = organ, y = fct_reorder(gene_id, log2cpm, .fun = max), fill = z_plot)) +
  geom_tile() +
  plot_theme(base_size = 10) +
  labs(
    title = wrap_plot_text("Relative expression of top organ-expressed genes", 78),
    subtitle = wrap_plot_text("Row-scaled log2(CPM + pseudocount); values clipped at ±2.5 only for color stability.", 92),
    x = "Organ",
    y = if (heatmap_show_gene_labels) "Gene" else "Gene labels hidden; see top_expressed_genes_per_organ.tsv",
    fill = "Row z-score"
  ) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 16)) +
  scale_y_discrete(labels = if (heatmap_show_gene_labels) waiver() else function(x) rep("", length(x))) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-2.5, 2.5),
    oob = scales::squish
  ) +
  coord_cartesian(clip = "off")

save_plot_all(p_heat, "top_expressed_genes_heatmap", width = 10, height = heatmap_height)

# Compatibility output for the previous tidyplots heatmap filename. The ggplot
# version is used here because tidyplots renders dense gene heatmaps as a small
# panel with overlapping gene labels and excessive whitespace.
save_plot_all(
  p_heat,
  "top_expressed_genes_heatmap_tidyplots",
  width = 10,
  height = heatmap_height
)

# Descriptive fold-change and volcano-like plots -------------------------------

if (run_descriptive_volcano == "TRUE") {
  organs <- metadata$organ
  names(organs) <- metadata$sample

  organ_samples <- metadata %>%
    select(sample, organ)

  if (!(denominator %in% organ_samples$organ)) {
    denominator <- organ_samples$organ[[1]]
  }

  denominator_sample <- organ_samples %>%
    filter(organ == denominator) %>%
    slice(1) %>%
    pull(sample)

  for (i in seq_len(nrow(organ_samples))) {
    numerator <- organ_samples$organ[[i]]
    numerator_sample <- organ_samples$sample[[i]]

    if (numerator == denominator) next

    fc_tbl <- tibble(
      gene_id = rownames(cpm_mat),
      numerator = numerator,
      denominator = denominator,
      numerator_cpm = cpm_mat[, numerator_sample],
      denominator_cpm = cpm_mat[, denominator_sample],
      mean_cpm = (numerator_cpm + denominator_cpm) / 2,
      log2_fold_change = log2((numerator_cpm + pseudocount) / (denominator_cpm + pseudocount)),
      abs_log2_fold_change = abs(log2_fold_change),
      direction = case_when(
        log2_fold_change >= lfc_threshold ~ paste0("Higher in ", numerator),
        log2_fold_change <= -lfc_threshold ~ paste0("Higher in ", denominator),
        TRUE ~ "Low fold-change"
      )
    ) %>%
      left_join(annotation_by_gene_like_id, by = "gene_id") %>%
      arrange(desc(abs_log2_fold_change))

    out_name <- paste0("descriptive_fold_change_", numerator, "_vs_", denominator)

    readr::write_tsv(fc_tbl, file.path(outdir, paste0(out_name, ".tsv")))

    top_labels <- fc_tbl %>%
      filter(abs_log2_fold_change >= lfc_threshold) %>%
      group_by(direction) %>%
      slice_max(abs_log2_fold_change, n = 5, with_ties = FALSE) %>%
      ungroup()

    p_fc <- fc_tbl %>%
      ggplot(aes(x = log2_fold_change, y = log10(mean_cpm + pseudocount), color = direction)) +
      geom_point(alpha = 0.35, size = 0.75) +
      geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = "dashed", linewidth = 0.35) +
      geom_text_repel(
        data = top_labels,
        aes(label = gene_id),
        size = 2.2,
        box.padding = 0.35,
        point.padding = 0.15,
        min.segment.length = 0,
        max.overlaps = 30,
        show.legend = FALSE
      ) +
      plot_theme(base_size = 12) +
      labs(
        title = wrap_plot_text(paste0("Descriptive fold-change plot: ", numerator, " vs ", denominator), 78),
        subtitle = wrap_plot_text("Descriptive only: one RNA-seq library per organ; no P-values estimated.", 92),
        x = expression(log[2]~fold~change),
        y = expression(log[10]~mean~CPM),
        color = "Class"
      ) +
      scale_x_continuous(expand = expansion(mult = 0.08)) +
      scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
      coord_cartesian(clip = "off")

    save_plot_all(p_fc, out_name, width = 10, height = 7)
  }
}

# Summary table using gtsummary -------------------------------------------------

metadata_summary <- metadata %>%
  select(sample, accession, organ) %>%
  tbl_summary(
    by = organ,
    statistic = list(all_categorical() ~ "{n}"),
    missing = "no"
  )

as_gt(metadata_summary) %>%
  tab_header(
    title = "RNA-seq sample metadata summary",
    subtitle = "One SRA run per organ"
  ) %>%
  gtsave(file.path(tabledir, "sample_metadata_summary_gtsummary.html"))

message("RNA-seq descriptive plots and tables completed.")
