#!/usr/bin/env Rscript

# egapx_genome_annotation_graphics_karyoploter.R
#
# Complete graphical summary of NCBI EGAPx genome annotation and functional
# annotation support. Produces ggplot2 circular/linear summaries, summary plots,
# and larger karyoploteR density and multiple-data-type plots with regions, markers, coverage, ribbons, and functional support layers. It also creates a tidyverse-derived EGAPx gene/transcript/protein ID map and per-method descriptive summaries/plots for eggNOG, InterProScan, DIAMOND/Swiss-Prot, KofamScan, and dbCAN.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(Rsamtools)
  library(GenomicRanges)
  library(Biostrings)
  library(ggplot2)
  library(scales)
  library(RColorBrewer)
  library(karyoploteR)
})

# Optional table/grammar packages are used when available.
# The script still writes TSV outputs if gt/gtsummary are not installed.

# ==============================================================================
# CLI
# ==============================================================================

opt_list <- list(
  make_option(c("--gff"), type = "character", default = "complete.genomic.gff",
              help = "EGAPx GFF/GFF3 file (default: %default)"),
  make_option(c("--fasta"), type = "character", default = "complete.genomic.fna",
              help = "Genome FASTA used by EGAPx (default: %default)"),
  make_option(c("--out"), type = "character", default = "egapx_genome_annotation",
              help = "Output prefix (default: %default)"),

  make_option(c("--top_n"), type = "integer", default = 15,
              help = "Plot only top N contigs by length; 0 = all after min_len filter (default: %default)"),
  make_option(c("--min_len"), type = "integer", default = 0,
              help = "Minimum contig length in bp (default: %default)"),
  make_option(c("--window"), type = "integer", default = 1000000,
              help = "Window size for density/GC tracks in bp (default: %default)"),
  make_option(c("--step"), type = "integer", default = 0,
              help = "Step size for windows in bp; 0 = non-overlapping windows (default: %default)"),

  make_option(c("--feature_counts"), type = "character", default = "",
              help = "Optional EGAPx feature_counts.txt summary"),
  make_option(c("--official_count_source"), type = "character", default = "egapx",
              help = "Official annotation-count source for non-spatial summary labels: egapx or gff (default: %default)"),
  make_option(c("--plot_count_source"), type = "character", default = "gff",
              help = "Spatial count source for genome tracks. Keep as gff because genomic coordinates come from GFF features (default: %default)"),
  make_option(c("--plot_count_comparison"), type = "logical", default = TRUE,
              help = "Export side-by-side EGAPx feature_counts.txt versus GFF-derived count comparison plots/tables (default: %default)"),
  make_option(c("--id_map"), type = "character", default = "",
              help = "Optional gene/transcript/protein ID map with columns gene_id and protein_id"),
  make_option(c("--eggnog"), type = "character", default = "",
              help = "Optional eggNOG-mapper .emapper.annotations file"),
  make_option(c("--interpro"), type = "character", default = "",
              help = "Optional InterProScan TSV file"),
  make_option(c("--diamond"), type = "character", default = "",
              help = "Optional DIAMOND/Swiss-Prot TSV file"),
  make_option(c("--kofam"), type = "character", default = "",
              help = "Optional KofamScan detail TSV file"),
  make_option(c("--dbcan_dir"), type = "character", default = "",
              help = "Optional dbCAN output directory"),
  make_option(c("--create_id_map"), type = "logical", default = TRUE,
              help = "Create EGAPx gene/transcript/protein map from the GFF if missing (default: %default)"),
  make_option(c("--force_id_map"), type = "logical", default = FALSE,
              help = "Recreate the EGAPx gene/transcript/protein map even if it already exists (default: %default)"),
  make_option(c("--id_map_out"), type = "character", default = "",
              help = "Output path for created gene/transcript/protein map; default uses --id_map if supplied or <out>.egapx_gene_transcript_protein_map.tsv"),
  make_option(c("--method_summary_plots"), type = "logical", default = TRUE,
              help = "Export descriptive plots and tables summarizing genes identified by each functional annotation method (default: %default)"),
  make_option(c("--method_karyoplots"), type = "logical", default = TRUE,
              help = "Export one karyoploteR multiple-data-type plot per functional annotation method (default: %default)"),
  make_option(c("--method_gene_table_top_n"), type = "integer", default = 100,
              help = "Number of per-method gene examples to include in optional gt HTML tables (default: %default)"),

  make_option(c("--gc"), type = "logical", default = TRUE,
              help = "Compute and plot GC fraction per window (default: %default)"),
  make_option(c("--cds_metric"), type = "character", default = "count",
              help = "CDS metric per window: count or bp (bp = covered bp fraction) (default: %default)"),
  make_option(c("--include_functional_tracks"), type = "logical", default = TRUE,
              help = "Add optional functional-support tracks when ID map and annotation files are available (default: %default)"),
  make_option(c("--include_problem_track"), type = "logical", default = TRUE,
              help = "Add corrected/problematic CDS density track (default: %default)"),

  make_option(c("--circular"), type = "logical", default = TRUE,
              help = "Export circular ggplot2 plot (default: %default)"),
  make_option(c("--linear"), type = "logical", default = TRUE,
              help = "Export linear faceted ggplot2 genome-track plot (default: %default)"),
  make_option(c("--summary_plots"), type = "logical", default = TRUE,
              help = "Export annotation summary and density-distribution plots (default: %default)"),
  make_option(c("--karyoplot_density"), type = "logical", default = TRUE,
              help = "Export karyoploteR density-track plot (default: %default)"),
  make_option(c("--karyoplot_multitype"), type = "logical", default = TRUE,
              help = "Export karyoploteR multiple-data-type plot (default: %default)"),

  make_option(c("--gene_rich_z"), type = "double", default = 2,
              help = "Absolute z-score threshold for gene-rich/gene-poor windows (default: %default)"),
  make_option(c("--clip_quant"), type = "double", default = 0.995,
              help = "Upper quantile used to clip color intensity (default: %default)"),
  make_option(c("--alpha_min"), type = "double", default = 0.55,
              help = "Minimum alpha for circular track bars (default: %default)"),
  make_option(c("--alpha_max"), type = "double", default = 1.00,
              help = "Maximum alpha for circular track bars (default: %default)"),

  make_option(c("--label_style"), type = "character", default = "chr",
              help = "Contig label style: chr or seqid (default: %default)"),
  make_option(c("--label_top_k"), type = "integer", default = 15,
              help = "Label only top K contigs in circular plot (default: %default)"),
  make_option(c("--major_tick_bp"), type = "integer", default = 50000000,
              help = "Major tick interval per contig in bp (default: %default)"),
  make_option(c("--minor_tick_bp"), type = "integer", default = 10000000,
              help = "Minor tick interval per contig in bp (default: %default)"),
  make_option(c("--show_tick_labels"), type = "logical", default = FALSE,
              help = "Show major tick labels around the circular ring (default: %default)"),
  make_option(c("--show_boundaries"), type = "logical", default = TRUE,
              help = "Draw contig boundary lines in circular plot (default: %default)"),
  make_option(c("--show_legend"), type = "logical", default = TRUE,
              help = "Add legend/caption annotations where supported (default: %default)"),

  make_option(c("--theme"), type = "character", default = "light",
              help = "Theme for circular plot: light or dark (default: %default)"),
  make_option(c("--palette"), type = "character", default = "Spectral",
              help = "Default palette: Spectral|RdYlBu|BrBG|YlOrRd|PuBuGn|viridis|turbo|plasma|magma|inferno|cividis (default: %default)"),
  make_option(c("--gene_palette"), type = "character", default = "Spectral",
              help = "Palette for all-gene track (default: %default)"),
  make_option(c("--coding_palette"), type = "character", default = "YlOrRd",
              help = "Palette for protein-coding track (default: %default)"),
  make_option(c("--noncoding_palette"), type = "character", default = "PuBuGn",
              help = "Palette for non-coding track (default: %default)"),
  make_option(c("--pseudo_palette"), type = "character", default = "Greys",
              help = "Palette for pseudogene track (default: %default)"),
  make_option(c("--cds_palette"), type = "character", default = "YlGnBu",
              help = "Palette for CDS track (default: %default)"),
  make_option(c("--problem_palette"), type = "character", default = "OrRd",
              help = "Palette for problem-CDS track (default: %default)"),
  make_option(c("--gc_palette"), type = "character", default = "BrBG",
              help = "Palette for GC track (default: %default)"),

  make_option(c("--width"), type = "double", default = 18,
              help = "Circular figure width in inches (default: %default)"),
  make_option(c("--height"), type = "double", default = 18,
              help = "Circular figure height in inches (default: %default)"),
  make_option(c("--linear_width"), type = "double", default = 18,
              help = "Linear figure width in inches (default: %default)"),
  make_option(c("--linear_height"), type = "double", default = 14,
              help = "Linear figure height in inches (default: %default)"),
  make_option(c("--karyo_width"), type = "double", default = 34,
              help = "karyoploteR figure width in inches (default: %default)"),
  make_option(c("--karyo_height"), type = "double", default = 26,
              help = "karyoploteR figure height in inches (default: %default)"),
  make_option(c("--karyo_base_cex"), type = "double", default = 0.55,
              help = "Base-number label size for karyoploteR plots (default: %default)"),
  make_option(c("--karyo_label_cex"), type = "double", default = 0.68,
              help = "Track-label size for karyoploteR plots (default: %default)"),
  make_option(c("--karyo_legend_cex"), type = "double", default = 0.85,
              help = "Legend size for karyoploteR plots (default: %default)"),
  make_option(c("--karyo_legend_inset"), type = "double", default = 0.030,
              help = "Positive inset for karyoploteR bottom legends. Larger values move legends upward and prevent bottom clipping (default: %default)"),
  make_option(c("--karyo_legend_x_intersp"), type = "double", default = 0.70,
              help = "Horizontal spacing multiplier for karyoploteR legends (default: %default)"),
  make_option(c("--karyo_legend_y_intersp"), type = "double", default = 0.85,
              help = "Vertical spacing multiplier for karyoploteR legends (default: %default)"),
  make_option(c("--karyo_track_bg"), type = "character", default = "#FAFAFA",
              help = "Background color for karyoploteR track panels (default: %default)"),
  make_option(c("--karyo_top_oma"), type = "double", default = 6.5,
              help = "Outer top margin for karyoploteR plots; increases space for titles/subtitles (default: %default)"),
  make_option(c("--karyo_bottom_oma"), type = "double", default = 8.0,
              help = "Outer bottom margin for karyoploteR plots; increases space for legends (default: %default)"),
  make_option(c("--karyo_left_oma"), type = "double", default = 1.2,
              help = "Outer left margin for karyoploteR plots (default: %default)"),
  make_option(c("--karyo_right_oma"), type = "double", default = 1.2,
              help = "Outer right margin for karyoploteR plots (default: %default)"),
  make_option(c("--karyo_title_cex"), type = "double", default = 1.25,
              help = "Outer title size for karyoploteR plots (default: %default)"),
  make_option(c("--karyo_subtitle_cex"), type = "double", default = 0.95,
              help = "Outer subtitle size for karyoploteR plots (default: %default)"),
  make_option(c("--karyo_title_line"), type = "double", default = 2.9,
              help = "Outer-margin line used for karyoploteR titles (default: %default)"),
  make_option(c("--karyo_subtitle_line"), type = "double", default = 1.55,
              help = "Outer-margin line used for karyoploteR subtitles (default: %default)"),
  make_option(c("--gc_ymin"), type = "double", default = 0.30,
              help = "Lower GC axis limit in karyoploteR GC tracks (default: %default)"),
  make_option(c("--gc_ymax"), type = "double", default = 0.55,
              help = "Upper GC axis limit in karyoploteR GC tracks (default: %default)"),
  make_option(c("--gc_ribbon_k"), type = "integer", default = 7,
              help = "Odd-number smoothing window for GC ribbon in karyoploteR multitype plot; 0 disables ribbon (default: %default)"),
  make_option(c("--top_marker_n"), type = "integer", default = 40,
              help = "Number of strongest gene-rich windows to label with kpPlotMarkers in multitype plot (default: %default)"),
  make_option(c("--plot_functional_coverage"), type = "logical", default = TRUE,
              help = "Plot combined functional-support coverage in karyoploteR multitype plot (default: %default)"),
  make_option(c("--plot_functional_subtracks"), type = "logical", default = TRUE,
              help = "Plot separate functional-support density subtracks in karyoploteR multitype plot (default: %default)"),
  make_option(c("--dpi"), type = "integer", default = 600,
              help = "PNG DPI (default: %default)"),

  make_option(c("--no_svg"), action = "store_true", default = FALSE,
              help = "Disable SVG export"),
  make_option(c("--no_pdf"), action = "store_true", default = FALSE,
              help = "Disable PDF export"),
  make_option(c("--no_png"), action = "store_true", default = FALSE,
              help = "Disable PNG export")
)

opt <- parse_args(OptionParser(option_list = opt_list))

if (!file.exists(opt$gff)) stop("GFF not found: ", opt$gff)
if (!file.exists(opt$fasta)) stop("FASTA not found: ", opt$fasta)
if (!(tolower(opt$cds_metric) %in% c("count", "bp"))) stop("--cds_metric must be count or bp")
if (!(tolower(opt$label_style) %in% c("chr", "seqid"))) stop("--label_style must be chr or seqid")
if (!(tolower(opt$official_count_source) %in% c("egapx", "gff"))) stop("--official_count_source must be egapx or gff")
if (tolower(opt$plot_count_source) != "gff") {
  warning("Spatial genome tracks require GFF-derived coordinates. Forcing plot_count_source='gff'.")
  opt$plot_count_source <- "gff"
}
if (!(tolower(opt$theme) %in% c("light", "dark"))) stop("--theme must be light or dark")

window <- as.integer(opt$window)
step <- as.integer(if (opt$step <= 0) opt$window else opt$step)

# ==============================================================================
# Helpers
# ==============================================================================

msg <- function(...) message(sprintf(...))

safe_file <- function(x) {
  is.character(x) && length(x) == 1L && nzchar(x) && file.exists(x) && file.info(x)$size > 0
}

safe_dir <- function(x) {
  is.character(x) && length(x) == 1L && nzchar(x) && dir.exists(x)
}

get_palette <- function(name, n = 256, reverse = FALSE) {
  nm <- tolower(name)
  brewer_names <- rownames(RColorBrewer::brewer.pal.info)
  brewer_match <- brewer_names[tolower(brewer_names) == nm]
  if (length(brewer_match) == 1L) {
    max_n <- RColorBrewer::brewer.pal.info[brewer_match, "maxcolors"]
    pal <- RColorBrewer::brewer.pal(max_n, brewer_match)
    pal <- grDevices::colorRampPalette(pal)(n)
    if (reverse) pal <- rev(pal)
    return(pal)
  }
  if (requireNamespace("viridisLite", quietly = TRUE)) {
    if (nm %in% c("viridis", "magma", "inferno", "plasma", "cividis", "turbo")) {
      pal <- viridisLite::viridis(n, option = nm)
      if (reverse) pal <- rev(pal)
      return(pal)
    }
  }
  if (nm %in% c("mako")) {
    pal <- grDevices::colorRampPalette(c("#0B0405", "#0F2B5B", "#1669A3", "#2FB7B0", "#F2F2F2"))(n)
    if (reverse) pal <- rev(pal)
    return(pal)
  }
  pal <- grDevices::colorRampPalette(c("#132B43", "#56B1F7", "#F7F7F7", "#F0766D", "#7F0000"))(n)
  if (reverse) pal <- rev(pal)
  pal
}

extract_attr <- function(x, key) {
  pattern <- paste0("(^|;)", key, "=([^;]+)")
  hit <- grepl(pattern, x)
  out <- rep(NA_character_, length(x))
  out[hit] <- sub(pattern, "\\2", x[hit])
  out[hit] <- sub(";.*$", "", out[hit])
  URLdecode(out)
}

clean_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("^\\*+", "", x)
  x <- sub("\\s+.*$", "", x)
  x <- sub("^>", "", x)
  x <- sub("^gene-", "", x)
  x <- sub("^rna-", "", x)
  x <- sub("^cds-", "", x)
  x <- sub("^protein-", "", x)
  x <- sub("^lcl\\|", "", x)
  x <- sub("^gnl\\|[^|]+\\|", "", x)
  x <- sub(".*\\|", "", x)
  x
}

id_variants <- function(x) {
  x <- unique(stats::na.omit(as.character(x)))
  x <- trimws(x)
  x <- x[nzchar(x)]
  if (length(x) == 0L) return(character())

  # Split common compound fields without destroying accession-like IDs.
  x_split <- unlist(strsplit(x, "[;,]"), use.names = FALSE)
  x_split <- trimws(x_split)
  x_split <- x_split[nzchar(x_split)]

  v <- unique(c(
    x,
    x_split,
    clean_id(x),
    clean_id(x_split),
    sub("\\.\\d+$", "", x),
    sub("\\.\\d+$", "", clean_id(x)),
    sub("^.*\\|", "", x),
    sub("^.*\\|", "", clean_id(x)),
    sub("^cds-", "", x),
    sub("^rna-", "", x),
    sub("^gene-", "", x),
    sub("^protein-", "", x)
  ))
  v <- trimws(v)
  v <- v[!is.na(v) & nzchar(v)]
  unique(v)
}

first_existing_col <- function(dt, candidates) {
  hit <- candidates[candidates %in% names(dt)]
  if (length(hit) == 0L) NA_character_ else hit[1L]
}

read_table_auto <- function(path) {
  fread(path, sep = "\t", header = TRUE, quote = "", fill = TRUE, data.table = TRUE)
}

# ------------------------------------------------------------------------------
# Tidyverse-based EGAPx gene/transcript/protein map creation
# ------------------------------------------------------------------------------

create_egapx_id_map_tidyverse <- function(gff_path, out_path) {
  if (!requireNamespace("readr", quietly = TRUE) ||
      !requireNamespace("dplyr", quietly = TRUE) ||
      !requireNamespace("stringr", quietly = TRUE)) {
    stop("Creating the ID map requires the tidyverse packages readr, dplyr and stringr in Renv.")
  }

  msg("[ID map] Creating tidyverse gene/transcript/protein map: %s", out_path)

  attr_extract_tv <- function(x, key) {
    pattern <- paste0("(?:^|;)", key, "=([^;]+)")
    out <- stringr::str_match(x, pattern)[, 2]
    utils::URLdecode(out)
  }

  transcript_types <- c(
    "mRNA", "transcript", "lnc_RNA", "ncRNA", "rRNA", "tRNA",
    "snRNA", "snoRNA", "miRNA", "pseudogenic_transcript"
  )

  gff_tbl <- readr::read_tsv(
    gff_path,
    comment = "#",
    col_names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attr"),
    col_types = readr::cols(.default = readr::col_character(), start = readr::col_integer(), end = readr::col_integer()),
    progress = FALSE
  )

  genes_tbl <- gff_tbl |>
    dplyr::filter(.data$type == "gene") |>
    dplyr::mutate(
      gene_id = clean_id(attr_extract_tv(.data$attr, "ID")),
      gene_name = attr_extract_tv(.data$attr, "Name"),
      locus_tag = attr_extract_tv(.data$attr, "locus_tag"),
      gene_biotype = dplyr::coalesce(
        attr_extract_tv(.data$attr, "gene_biotype"),
        attr_extract_tv(.data$attr, "biotype")
      )
    ) |>
    dplyr::select(.data$gene_id, .data$gene_name, .data$locus_tag, .data$seqid, .data$start, .data$end, .data$strand, .data$gene_biotype)

  tx_tbl <- gff_tbl |>
    dplyr::filter(.data$type %in% transcript_types) |>
    dplyr::mutate(
      transcript_id = clean_id(attr_extract_tv(.data$attr, "ID")),
      gene_id = clean_id(attr_extract_tv(.data$attr, "Parent")),
      transcript_biotype = dplyr::coalesce(
        attr_extract_tv(.data$attr, "transcript_biotype"),
        attr_extract_tv(.data$attr, "biotype")
      ),
      transcript_product = attr_extract_tv(.data$attr, "product")
    ) |>
    dplyr::select(.data$gene_id, .data$transcript_id, .data$transcript_biotype, .data$transcript_product)

  cds_tbl <- gff_tbl |>
    dplyr::filter(.data$type == "CDS") |>
    dplyr::mutate(
      transcript_id = clean_id(attr_extract_tv(.data$attr, "Parent")),
      protein_id = clean_id(dplyr::coalesce(
        attr_extract_tv(.data$attr, "protein_id"),
        attr_extract_tv(.data$attr, "ID"),
        attr_extract_tv(.data$attr, "Name")
      )),
      product = attr_extract_tv(.data$attr, "product")
    ) |>
    dplyr::filter(!is.na(.data$transcript_id), !is.na(.data$protein_id), .data$protein_id != "") |>
    dplyr::group_by(.data$transcript_id, .data$protein_id) |>
    dplyr::summarise(product = paste(unique(stats::na.omit(.data$product)), collapse = ";"), .groups = "drop")

  id_map_tbl <- tx_tbl |>
    dplyr::left_join(cds_tbl, by = "transcript_id") |>
    dplyr::left_join(genes_tbl, by = "gene_id") |>
    dplyr::mutate(product = dplyr::coalesce(.data$product, .data$transcript_product)) |>
    dplyr::select(
      .data$gene_id, .data$gene_name, .data$locus_tag, .data$transcript_id, .data$protein_id,
      .data$seqid, .data$start, .data$end, .data$strand, .data$gene_biotype,
      .data$transcript_biotype, .data$product
    ) |>
    dplyr::arrange(.data$seqid, .data$start, .data$gene_id, .data$transcript_id, .data$protein_id)

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(id_map_tbl, out_path)
  invisible(out_path)
}

make_gene_coord_all_from_gff <- function(gff_path) {
  dt <- suppressWarnings(fread(
    gff_path, sep = "\t", header = FALSE, quote = "", comment.char = "#",
    col.names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attr"),
    colClasses = c("character", "character", "character", "integer", "integer",
                   "character", "character", "character", "character")
  ))
  dt <- dt[type == "gene"]
  if (nrow(dt) == 0L) return(data.table(gene_id = character(), seqid = character(), start = integer(), end = integer()))
  dt[, ID := extract_attr(attr, "ID")]
  dt[, gene_id_attr := extract_attr(attr, "gene_id")]
  dt[, gene_id := clean_id(fifelse(!is.na(gene_id_attr), gene_id_attr, ID))]
  dt <- dt[!is.na(gene_id) & gene_id != "", .(gene_id, seqid, start, end)]
  dt <- merge(dt, seq_df[, .(seqid, plot_label)], by = "seqid", all.x = TRUE)
  unique(dt)
}

save_plot_all <- function(plot, prefix, width, height, dpi = 600, bg = "white") {
  if (!isTRUE(opt$no_png)) {
    ggsave(
      paste0(prefix, ".png"), plot,
      width = width, height = height, units = "in", dpi = dpi, bg = bg,
      device = if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
    )
  }
  if (!isTRUE(opt$no_pdf)) {
    ggsave(paste0(prefix, ".pdf"), plot, width = width, height = height, units = "in", bg = bg)
  }
  if (!isTRUE(opt$no_svg)) {
    ggsave(
      paste0(prefix, ".svg"), plot,
      width = width, height = height, units = "in", bg = bg,
      device = if (requireNamespace("svglite", quietly = TRUE)) svglite::svglite else "svg"
    )
  }
}

open_device <- function(filename, device, width, height, dpi = 600, bg = "white") {
  if (device == "png") {
    if (requireNamespace("ragg", quietly = TRUE)) {
      ragg::agg_png(filename, width = width, height = height, units = "in", res = dpi, background = bg)
    } else {
      png(filename, width = width, height = height, units = "in", res = dpi, bg = bg)
    }
  } else if (device == "pdf") {
    pdf(filename, width = width, height = height, bg = bg)
  } else if (device == "svg") {
    if (requireNamespace("svglite", quietly = TRUE)) {
      svglite::svglite(filename, width = width, height = height, bg = bg)
    } else {
      svg(filename, width = width, height = height, bg = bg)
    }
  } else {
    stop("Unsupported device: ", device)
  }
}

add_karyo_outer_title <- function(title, subtitle = "") {
  # karyoploteR uses base graphics. Drawing titles in the outer margin avoids
  # clipping when many chromosomes/tracks are plotted or when PDF/SVG devices
  # use tight bounding boxes. Increase --karyo_top_oma if long titles are used.
  graphics::mtext(
    title,
    side = 3,
    outer = TRUE,
    line = opt$karyo_title_line,
    cex = opt$karyo_title_cex,
    font = 2
  )

  if (is.character(subtitle) && length(subtitle) == 1L && nzchar(subtitle)) {
    graphics::mtext(
      subtitle,
      side = 3,
      outer = TRUE,
      line = opt$karyo_subtitle_line,
      cex = opt$karyo_subtitle_cex
    )
  }
}

run_karyo_with_margins <- function(plot_fun) {
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar), add = TRUE)

  graphics::par(
    oma = c(
      opt$karyo_bottom_oma,
      opt$karyo_left_oma,
      opt$karyo_top_oma,
      opt$karyo_right_oma
    ),
    xpd = NA
  )

  plot_fun()
}

save_karyo_all <- function(plot_fun, prefix, width, height, dpi = 600, bg = "white") {
  if (!isTRUE(opt$no_png)) {
    open_device(paste0(prefix, ".png"), "png", width, height, dpi, bg)
    run_karyo_with_margins(plot_fun)
    dev.off()
  }
  if (!isTRUE(opt$no_pdf)) {
    open_device(paste0(prefix, ".pdf"), "pdf", width, height, dpi, bg)
    run_karyo_with_margins(plot_fun)
    dev.off()
  }
  if (!isTRUE(opt$no_svg)) {
    open_device(paste0(prefix, ".svg"), "svg", width, height, dpi, bg)
    run_karyo_with_margins(plot_fun)
    dev.off()
  }
}

clip_upper <- function(x, q) {
  if (all(!is.finite(x))) return(x)
  u <- as.numeric(stats::quantile(x, probs = q, na.rm = TRUE, names = FALSE))
  pmin(x, u)
}

scale_to_alpha <- function(x, alpha_min, alpha_max) {
  rng <- range(x, finite = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) {
    return(rep((alpha_min + alpha_max) / 2, length(x)))
  }
  alpha_min + (x - rng[1]) / (rng[2] - rng[1]) * (alpha_max - alpha_min)
}

rescale_to_ring <- function(x, y0, y1) {
  rng <- range(x, finite = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) {
    return(rep((y0 + y1) / 2, length(x)))
  }
  y0 + (x - rng[1]) / (rng[2] - rng[1]) * (y1 - y0)
}

dt_to_gr <- function(dt, seq_col = "seqid", start_col = "start", end_col = "end") {
  if (nrow(dt) == 0L) return(GRanges())
  GRanges(
    seqnames = dt[[seq_col]],
    ranges = IRanges(start = as.integer(dt[[start_col]]), end = as.integer(dt[[end_col]]))
  )
}

# ==============================================================================
# FASTA index and contig lengths
# ==============================================================================

fai_path <- paste0(opt$fasta, ".fai")
if (!file.exists(fai_path)) {
  msg("[FASTA] Index not found; creating .fai with Rsamtools::indexFa()")
  Rsamtools::indexFa(opt$fasta)
}

fa <- FaFile(opt$fasta)
open(fa)
on.exit({ try(close(fa), silent = TRUE) }, add = TRUE)

fai <- scanFaIndex(fa)
seq_df <- data.table(
  seqid = as.character(seqnames(fai)),
  length_int = as.integer(seqlengths(fai))
)
seq_df[, length := as.numeric(length_int)]

seq_df <- seq_df[length_int >= opt$min_len]
setorder(seq_df, -length_int, seqid)
if (nrow(seq_df) == 0L) stop("No contigs passed filtering.")

if (opt$top_n > 0 && nrow(seq_df) > opt$top_n) {
  seq_df <- seq_df[seq_len(opt$top_n)]
}

seq_df[, chr_label := paste0("Chr", seq_len(.N))]
seq_df[, plot_label := if (tolower(opt$label_style) == "chr") chr_label else seqid]
keep_seq <- seq_df$seqid

msg("[Genome] Contigs plotted: %d", length(keep_seq))
msg("[Genome] Total plotted bp: %s", format(sum(seq_df$length), scientific = FALSE, big.mark = ","))

seq_df[, offset := c(0, head(cumsum(length), -1L))]
seq_df[, start_global := offset]
seq_df[, end_global := offset + length]
global_len <- sum(seq_df$length)

fwrite(seq_df[, .(chr_label, plot_label, seqid, length, start_global, end_global)],
       paste0(opt$out, ".contig_label_map.tsv"), sep = "\t")

karyo_genome <- toGRanges(seq_df[, .(chr = plot_label, start = 1L, end = length_int)])
chr_map <- seq_df[, .(original_seqid = seqid, chr_label = plot_label, length = length_int)]

# ==============================================================================
# Genome-wide GC and GFF
# ==============================================================================

msg("[Genome] Computing genome-wide GC%% from plotted contigs...")
chunk_bp <- 5e6

compute_gc_streamed <- function(seqname, seqlen_int) {
  starts <- seq.int(1L, seqlen_int, by = chunk_bp)
  ends <- pmin(starts + chunk_bp - 1L, seqlen_int)
  gr <- GRanges(seqnames = seqname, ranges = IRanges(starts, ends))
  s <- getSeq(fa, gr)
  gc <- letterFrequency(s, letters = c("G", "C"), as.prob = FALSE)
  list(gc = sum(gc[, "G"] + gc[, "C"]), n = sum(width(s)))
}

gc_tot <- 0
n_tot <- 0
for (i in seq_len(nrow(seq_df))) {
  r <- compute_gc_streamed(seq_df$seqid[i], seq_df$length_int[i])
  gc_tot <- gc_tot + r$gc
  n_tot <- n_tot + r$n
}
genome_gc <- if (n_tot > 0) gc_tot / n_tot else NA_real_

msg("[GFF] Reading: %s", opt$gff)
gff <- suppressWarnings(
  fread(
    opt$gff,
    sep = "\t",
    header = FALSE,
    quote = "",
    comment.char = "#",
    col.names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attr"),
    colClasses = c(
      "character", "character", "character", "integer", "integer",
      "character", "character", "character", "character"
    )
  )
)

gff <- gff[seqid %in% keep_seq]
gff <- merge(gff, seq_df[, .(seqid, plot_label)], by = "seqid", all.x = TRUE)

gff[, `:=`(
  ID = extract_attr(attr, "ID"),
  Parent = extract_attr(attr, "Parent"),
  gene_id_attr = extract_attr(attr, "gene_id"),
  gene_biotype = extract_attr(attr, "gene_biotype"),
  biotype = extract_attr(attr, "biotype"),
  product = extract_attr(attr, "product"),
  exception = extract_attr(attr, "exception"),
  protein_id_attr = extract_attr(attr, "protein_id"),
  name_attr = extract_attr(attr, "Name")
)]

gff[, ID_clean := clean_id(ID)]
gff[, Parent_clean := clean_id(Parent)]
gff[, gene_id_clean := clean_id(fifelse(!is.na(gene_id_attr), gene_id_attr, ID))]

attr_all <- tolower(gff$attr)
gff[, pseudo := grepl("pseudo=true|pseudogene|pseudo", attr_all)]
gff[, partial := grepl("partial=true|partial", attr_all)]
gff[, frameshift := grepl("frameshift", attr_all)]
gff[, premature_stop := grepl("premature.stop|premature stop|stop codon|internal stop", attr_all)]
gff[, correction := grepl("correction|corrected", attr_all)]

# ==============================================================================
# Feature classes
# ==============================================================================

genes_all <- gff[type == "gene", .(
  seqid, plot_label, start, end, ID, ID_clean, gene_id_clean, gene_biotype, biotype, pseudo, attr
)]
genes_all[, class_text := tolower(paste(gene_biotype, biotype, attr))]

genes_pc <- genes_all[
  grepl("protein_coding|protein-coding", class_text, ignore.case = TRUE) & !pseudo,
  .(seqid, plot_label, start, end, gene_id = gene_id_clean)
]

genes_ncrna <- genes_all[
  grepl("lncrna|ncrna|rrna|trna|snrna|snorna|mirna|antisense|non.coding|non-coding|misc_rna",
        class_text, ignore.case = TRUE) & !pseudo,
  .(seqid, plot_label, start, end, gene_id = gene_id_clean)
]

genes_pseudo <- genes_all[
  pseudo | grepl("pseudogene|pseudo", class_text, ignore.case = TRUE),
  .(seqid, plot_label, start, end, gene_id = gene_id_clean)
]

if (nrow(genes_pc) == 0L && nrow(genes_ncrna) == 0L && nrow(genes_pseudo) == 0L) {
  warning("No explicit gene classes detected from attributes. Broad gene track will still be available.")
}

genes_broad <- genes_all[, .(seqid, plot_label, start, end, gene_id = gene_id_clean)]

cds_all <- gff[type == "CDS", .(
  seqid, plot_label, start, end, ID, ID_clean, Parent, Parent_clean,
  protein_id_attr = clean_id(protein_id_attr),
  name_attr = clean_id(name_attr),
  partial, frameshift, premature_stop, correction, attr
)]

cds_problem <- cds_all[
  partial | frameshift | premature_stop | correction,
  .(seqid, plot_label, start, end)
]

msg("[GFF] Genes: %s | protein-coding: %s | ncRNA genes: %s | pseudogenes: %s | CDS: %s | problem CDS: %s",
    format(nrow(genes_broad), big.mark = ","),
    format(nrow(genes_pc), big.mark = ","),
    format(nrow(genes_ncrna), big.mark = ","),
    format(nrow(genes_pseudo), big.mark = ","),
    format(nrow(cds_all), big.mark = ","),
    format(nrow(cds_problem), big.mark = ","))

gene_coord <- genes_broad[, .(gene_id, seqid, plot_label, start, end)]
gene_coord <- gene_coord[!is.na(gene_id) & gene_id != ""]
gene_coord <- unique(gene_coord)

rna_to_gene <- gff[type %in% c("mRNA", "transcript", "lnc_RNA", "ncRNA", "tRNA", "rRNA", "miRNA", "snoRNA", "snRNA"),
                   .(transcript_id = ID_clean, gene_id = Parent_clean)]
rna_to_gene <- rna_to_gene[!is.na(transcript_id) & !is.na(gene_id)]
rna_to_gene <- unique(rna_to_gene)

cds_to_parent <- rbindlist(list(
  cds_all[!is.na(ID_clean) & ID_clean != "" & !is.na(Parent_clean),
          .(protein_id = ID_clean, transcript_id = Parent_clean)],
  cds_all[!is.na(protein_id_attr) & protein_id_attr != "" & !is.na(Parent_clean),
          .(protein_id = protein_id_attr, transcript_id = Parent_clean)],
  cds_all[!is.na(name_attr) & name_attr != "" & !is.na(Parent_clean),
          .(protein_id = name_attr, transcript_id = Parent_clean)]
), use.names = TRUE, fill = TRUE)
cds_to_parent <- unique(cds_to_parent[!is.na(protein_id) & protein_id != "" & !is.na(transcript_id)])
cds_to_gene <- merge(cds_to_parent, rna_to_gene, by = "transcript_id", all.x = TRUE)
cds_to_gene[, gene_id := fifelse(is.na(gene_id), transcript_id, gene_id)]
cds_to_gene <- unique(cds_to_gene[, .(protein_id, transcript_id, gene_id)])

genes_broad_gr <- dt_to_gr(genes_broad, seq_col = "plot_label")
genes_pc_gr <- dt_to_gr(genes_pc, seq_col = "plot_label")
genes_ncrna_gr <- dt_to_gr(genes_ncrna, seq_col = "plot_label")
genes_pseudo_gr <- dt_to_gr(genes_pseudo, seq_col = "plot_label")
cds_all_gr <- dt_to_gr(cds_all, seq_col = "plot_label")
cds_problem_gr <- dt_to_gr(cds_problem, seq_col = "plot_label")

# ==============================================================================
# Windows and track calculations
# ==============================================================================

make_windows_dt <- function(seqid, seqlen_int, win, step) {
  starts <- seq.int(1L, seqlen_int, by = step)
  ends <- pmin(starts + win - 1L, seqlen_int)
  data.table(seqid = seqid, start = starts, end = ends)
}

wins <- rbindlist(lapply(seq_len(nrow(seq_df)), function(i) {
  make_windows_dt(seq_df$seqid[i], seq_df$length_int[i], window, step)
}), use.names = TRUE)

wins[, win_id := .I]
wins <- merge(wins, seq_df[, .(seqid, plot_label)], by = "seqid", all.x = TRUE)

calc_counts_by_overlap <- function(feat_dt, wins_dt) {
  if (nrow(feat_dt) == 0L) return(wins_dt[, .(win_id, count = 0L)])
  feat <- feat_dt[, .(seqid, start, end)]
  feat <- feat[!is.na(seqid) & !is.na(start) & !is.na(end)]
  if (nrow(feat) == 0L) return(wins_dt[, .(win_id, count = 0L)])
  setkey(feat, seqid, start, end)
  w <- wins_dt[, .(seqid, start, end, win_id)]
  setkey(w, seqid, start, end)
  hits <- foverlaps(feat, w, nomatch = 0L)
  out <- hits[, .N, by = win_id]
  setnames(out, "N", "count")
  out <- merge(wins_dt[, .(win_id)], out, by = "win_id", all.x = TRUE)
  out[is.na(count), count := 0L]
  out[]
}

calc_cds_bp_fraction <- function(cds_dt, wins_dt) {
  if (nrow(cds_dt) == 0L) return(wins_dt[, .(win_id, seqid, plot_label, start, end, value = 0)])
  wins_gr <- GRanges(seqnames = wins_dt$seqid, ranges = IRanges(wins_dt$start, wins_dt$end), win_id = wins_dt$win_id)
  cds_gr <- GRanges(seqnames = cds_dt$seqid, ranges = IRanges(cds_dt$start, cds_dt$end))
  cds_gr <- reduce(cds_gr)
  ol <- findOverlaps(wins_gr, cds_gr, ignore.strand = TRUE)
  if (length(ol) == 0L) return(wins_dt[, .(win_id, seqid, plot_label, start, end, value = 0)])
  wi <- queryHits(ol)
  ci <- subjectHits(ol)
  inter <- pintersect(ranges(wins_gr)[wi], ranges(cds_gr)[ci])
  bp <- width(inter)
  dt <- data.table(win_id = mcols(wins_gr)$win_id[wi], bp = bp)
  dt <- dt[, .(bp = sum(bp)), by = win_id]
  out <- merge(wins_dt[, .(win_id, seqid, plot_label, start, end)], dt, by = "win_id", all.x = TRUE)
  out[is.na(bp), bp := 0L]
  out[, value := bp / (end - start + 1L)]
  out[]
}

make_count_track <- function(feat_dt, label, group, per_mb = TRUE) {
  counts <- calc_counts_by_overlap(feat_dt, wins)
  out <- merge(wins[, .(win_id, seqid, plot_label, start, end)], counts, by = "win_id")
  out[, value := if (per_mb) count / (window / 1e6) else count]
  out[, `:=`(track = label, group = group)]
  out[]
}

msg("[Tracks] Computing structural tracks...")
track_list <- list(
  make_count_track(genes_broad, "All genes / Mb", "Structural"),
  make_count_track(genes_pc, "Protein-coding genes / Mb", "Structural"),
  make_count_track(genes_ncrna, "Non-coding genes / Mb", "Structural"),
  make_count_track(genes_pseudo, "Pseudogenes / Mb", "Structural")
)

if (tolower(opt$cds_metric) == "count") {
  track_list[["CDS / Mb"]] <- make_count_track(cds_all[, .(seqid, start, end)], "CDS / Mb", "Structural")
} else {
  cds_bp <- calc_cds_bp_fraction(cds_all, wins)
  cds_bp[, `:=`(count = NA_integer_, track = "CDS bp fraction", group = "Structural")]
  track_list[["CDS bp fraction"]] <- cds_bp
}

if (isTRUE(opt$include_problem_track)) {
  track_list[["Problem CDS / Mb"]] <- make_count_track(cds_problem, "Corrected/problem CDS / Mb", "Quality")
}

if (isTRUE(opt$gc)) {
  msg("[Tracks] Computing GC fraction per window...")
  wins_gr <- GRanges(seqnames = wins$seqid, ranges = IRanges(wins$start, wins$end), win_id = wins$win_id)
  seqs <- getSeq(fa, wins_gr)
  gc_counts <- letterFrequency(seqs, letters = c("G", "C"), as.prob = FALSE)
  gc_frac <- (gc_counts[, "G"] + gc_counts[, "C"]) / width(seqs)
  track_gc <- data.table(
    win_id = mcols(wins_gr)$win_id,
    seqid = as.character(seqnames(wins_gr)),
    start = start(wins_gr),
    end = end(wins_gr),
    count = NA_integer_,
    value = as.numeric(gc_frac),
    track = "GC fraction",
    group = "Composition"
  )
  track_gc <- merge(track_gc, seq_df[, .(seqid, plot_label)], by = "seqid", all.x = TRUE)
  track_list[["GC fraction"]] <- track_gc
}

# ==============================================================================
# Optional functional annotation tracks and method-level summaries
# ==============================================================================

# Create the ID map from the same GFF when requested. This standardizes all
# method summaries and spatial plots to the GFF-derived gene coordinates.
id_map_path <- opt$id_map
if (!nzchar(id_map_path)) {
  id_map_path <- if (nzchar(opt$id_map_out)) opt$id_map_out else paste0(opt$out, ".egapx_gene_transcript_protein_map.tsv")
}
if (isTRUE(opt$create_id_map) && (!safe_file(id_map_path) || isTRUE(opt$force_id_map))) {
  id_map_path <- create_egapx_id_map_tidyverse(opt$gff, id_map_path)
}
opt$id_map <- id_map_path

gene_coord_all <- make_gene_coord_all_from_gff(opt$gff)
fwrite(gene_coord_all, paste0(opt$out, ".gff_gene_coordinates_all.tsv"), sep = "\t")

read_id_map <- function(path) {
  if (!safe_file(path)) return(data.table())
  dt <- read_table_auto(path)
  gene_col <- first_existing_col(dt, c("gene_id", "GeneID", "gene", "locus_tag"))
  prot_col <- first_existing_col(dt, c("protein_id", "ProteinID", "protein", "cds_id", "CDS_ID"))
  tx_col <- first_existing_col(dt, c("transcript_id", "rna_id", "mrna_id", "transcript"))
  out <- data.table()
  if (!is.na(gene_col) && !is.na(prot_col)) {
    out <- dt[, .(gene_id = clean_id(get(gene_col)), protein_id = clean_id(get(prot_col)))]
  } else if (!is.na(gene_col) && !is.na(tx_col)) {
    out <- dt[, .(gene_id = clean_id(get(gene_col)), transcript_id = clean_id(get(tx_col)))]
  }
  unique(out)
}

id_map_user <- read_id_map(opt$id_map)
id_map <- unique(rbindlist(list(
  if (nrow(cds_to_gene) > 0L) cds_to_gene[, .(gene_id, protein_id)] else data.table(),
  if (nrow(id_map_user) > 0L && "protein_id" %in% names(id_map_user)) id_map_user[, .(gene_id, protein_id)] else data.table()
), use.names = TRUE, fill = TRUE))
id_map <- id_map[!is.na(gene_id) & !is.na(protein_id) & gene_id != "" & protein_id != ""]
fwrite(id_map, paste0(opt$out, ".protein_to_gene_id_map.used.tsv"), sep = "\t")

# Build a tolerant lookup table because functional tools may report protein IDs
# in slightly different forms: raw FASTA ID, cds-* ID, accession.version,
# accession without version, lcl|*, gnl|*, or pipe-delimited IDs.
make_id_lookup <- function(map_dt) {
  if (nrow(map_dt) == 0L) return(data.table(protein_key = character(), protein_id = character(), gene_id = character()))
  out <- rbindlist(lapply(seq_len(nrow(map_dt)), function(i) {
    data.table(
      protein_key = id_variants(map_dt$protein_id[i]),
      protein_id = map_dt$protein_id[i],
      gene_id = map_dt$gene_id[i]
    )
  }), use.names = TRUE, fill = TRUE)
  out <- out[!is.na(protein_key) & protein_key != "" & !is.na(gene_id) & gene_id != ""]
  unique(out)
}

id_lookup <- make_id_lookup(id_map)
fwrite(id_lookup, paste0(opt$out, ".protein_to_gene_id_lookup.expanded.tsv"), sep = "\t")

protein_ids_to_gene_coords_all <- function(protein_ids) {
  keys <- id_variants(protein_ids)
  if (length(keys) == 0L) {
    return(data.table(gene_id = character(), seqid = character(), plot_label = character(), start = integer(), end = integer()))
  }

  # Primary route: method protein/query IDs -> expanded protein lookup -> GFF genes.
  genes_from_protein <- if (nrow(id_lookup) > 0L) unique(id_lookup[protein_key %in% keys, gene_id]) else character()

  # Fallback route: some tools may report gene IDs instead of protein IDs.
  gene_keys <- unique(clean_id(keys))
  genes_direct <- unique(gene_coord_all[gene_id %in% gene_keys, gene_id])

  genes <- unique(c(genes_from_protein, genes_direct))
  out <- gene_coord_all[gene_id %in% genes, .(gene_id, seqid, plot_label, start, end)]
  unique(out)
}

protein_ids_to_gene_coords <- function(protein_ids) {
  out <- protein_ids_to_gene_coords_all(protein_ids)
  out <- out[!is.na(plot_label) & plot_label != "", .(seqid, plot_label, start, end)]
  unique(out)
}

make_method_gene_table <- function(method_name, protein_ids) {
  keys <- id_variants(protein_ids)
  if (length(keys) == 0L) {
    return(data.table(method = method_name, protein_id = character(), gene_id = character(),
                      seqid = character(), plot_label = character(), start = integer(), end = integer(),
                      on_plotted_contig = logical(), mapping_source = character()))
  }

  map_prot <- if (nrow(id_lookup) > 0L) {
    unique(id_lookup[protein_key %in% keys, .(protein_id, gene_id)])
  } else {
    data.table(protein_id = character(), gene_id = character())
  }
  if (nrow(map_prot) > 0L) map_prot[, mapping_source := "protein_id_map"]

  map_gene <- unique(gene_coord_all[gene_id %in% clean_id(keys), .(gene_id)])
  if (nrow(map_gene) > 0L) {
    map_gene[, `:=`(protein_id = gene_id, mapping_source = "direct_gene_id")]
    setcolorder(map_gene, c("protein_id", "gene_id", "mapping_source"))
  } else {
    map_gene <- data.table(protein_id = character(), gene_id = character(), mapping_source = character())
  }

  map <- unique(rbindlist(list(map_prot, map_gene), use.names = TRUE, fill = TRUE))
  if (nrow(map) == 0L) {
    return(data.table(method = method_name, protein_id = character(), gene_id = character(),
                      seqid = character(), plot_label = character(), start = integer(), end = integer(),
                      on_plotted_contig = logical(), mapping_source = character()))
  }

  out <- merge(map, gene_coord_all, by = "gene_id", all.x = TRUE, allow.cartesian = TRUE)
  out[, method := method_name]
  out[, on_plotted_contig := !is.na(plot_label) & plot_label != ""]
  setcolorder(out, c("method", "gene_id", "protein_id", "mapping_source", "seqid", "plot_label", "start", "end", "on_plotted_contig"))
  unique(out)
}

read_tsv_keep_hash_header <- function(path, expected_header_regex = NULL) {
  if (!safe_file(path)) return(data.table())
  ln <- readLines(path, warn = FALSE)
  ln <- ln[nzchar(ln)]
  # Keep single-# header lines such as eggNOG '#query', but remove metadata lines.
  ln <- ln[!grepl("^##", ln)]
  if (!is.null(expected_header_regex)) {
    idx <- grep(expected_header_regex, ln, ignore.case = TRUE)
    if (length(idx) > 0L) ln <- ln[idx[1L]:length(ln)]
  }
  if (length(ln) == 0L) return(data.table())
  ln[1L] <- sub("^#", "", ln[1L])
  fread(text = paste(ln, collapse = "\n"), sep = "\t", header = TRUE, quote = "", fill = TRUE, data.table = TRUE)
}

read_eggnog_proteins <- function(path) {
  dt <- read_tsv_keep_hash_header(path, expected_header_regex = "^#?query\\t|^#?query_name\\t")
  if (nrow(dt) == 0L) return(character())
  qcol <- first_existing_col(dt, c("query", "#query", "query_name", "qseqid"))
  if (is.na(qcol)) qcol <- names(dt)[1L]
  unique(dt[[qcol]])
}

read_interpro_proteins <- function(path) {
  if (!safe_file(path)) return(character())
  dt <- fread(path, sep = "\t", header = FALSE, quote = "", fill = TRUE, data.table = TRUE)
  if (ncol(dt) < 1L || nrow(dt) == 0L) return(character())
  ids <- dt[[1L]]
  ids <- ids[!is.na(ids) & ids != "" & !grepl("^#", ids)]
  unique(ids)
}

read_diamond_proteins <- function(path) {
  if (!safe_file(path)) return(character())
  # DIAMOND with --header simple may start with a commented header. Without
  # --header, the first column is still qseqid. Try both safely.
  dt <- tryCatch(read_tsv_keep_hash_header(path, expected_header_regex = "^#?qseqid\\t|^#?query\\t"),
                 error = function(e) data.table())
  if (nrow(dt) > 0L) {
    qcol <- first_existing_col(dt, c("qseqid", "query", "#query"))
    if (is.na(qcol)) qcol <- names(dt)[1L]
    return(unique(dt[[qcol]]))
  }
  dt <- fread(path, sep = "\t", header = FALSE, quote = "", fill = TRUE, data.table = TRUE)
  if (ncol(dt) < 1L) return(character())
  unique(dt[[1L]])
}

read_kofam_proteins <- function(path) {
  if (!safe_file(path)) return(character())
  dt <- fread(path, sep = "\t", header = FALSE, quote = "", fill = TRUE, data.table = TRUE)
  if (nrow(dt) == 0L || ncol(dt) < 1L) return(character())

  # detail-tsv often has '*' in column 1 for above-threshold hits and the query
  # in column 2. Other formats may put the query in column 1.
  first <- as.character(dt[[1L]])
  if (ncol(dt) >= 2L && any(first %in% c("*", ""), na.rm = TRUE)) {
    candidate <- dt[[2L]]
  } else {
    candidate <- dt[[1L]]
  }
  candidate <- as.character(candidate)
  candidate <- candidate[!is.na(candidate) & candidate != "" & !grepl("^#", candidate)]
  unique(candidate)
}

read_dbcan_proteins <- function(dbcan_dir) {
  if (!safe_dir(dbcan_dir)) return(character())
  files <- list.files(dbcan_dir, pattern = "overview|cazy|CAZy|dbcan|dbCAN", full.names = TRUE, ignore.case = TRUE)
  files <- files[file.info(files)$size > 0]
  if (length(files) == 0L) return(character())
  ids <- character()
  for (f in files) {
    dt <- tryCatch(fread(f, sep = "\t", header = TRUE, quote = "", fill = TRUE, data.table = TRUE),
                   error = function(e) data.table())
    if (nrow(dt) == 0L) next
    qcol <- first_existing_col(dt, c("Gene ID", "gene_id", "GeneID", "protein_id", "Protein ID", "seqid", "ID", "name"))
    if (!is.na(qcol)) ids <- c(ids, dt[[qcol]])
  }
  unique(ids)
}

functional_gr_list <- list()
method_gr_list <- list()
method_gene_list_dt <- data.table()
method_summary_dt <- data.table()

method_protein_ids <- list(
  "EggNOG" = read_eggnog_proteins(opt$eggnog),
  "InterProScan" = read_interpro_proteins(opt$interpro),
  "DIAMOND Swiss-Prot" = read_diamond_proteins(opt$diamond),
  "KofamScan KEGG KO" = read_kofam_proteins(opt$kofam),
  "dbCAN CAZyme" = read_dbcan_proteins(opt$dbcan_dir)
)

method_input_diagnostics <- rbindlist(lapply(names(method_protein_ids), function(nm) {
  ids <- unique(as.character(method_protein_ids[[nm]]))
  ids <- ids[!is.na(ids) & ids != ""]
  keys <- id_variants(ids)
  data.table(
    method = nm,
    input_ids = length(ids),
    expanded_match_keys = length(keys),
    keys_matching_protein_lookup = if (exists("id_lookup") && nrow(id_lookup) > 0L) uniqueN(id_lookup[protein_key %in% keys, protein_key]) else 0L,
    example_input_ids = paste(utils::head(ids, 5), collapse = ";")
  )
}), use.names = TRUE, fill = TRUE)
fwrite(method_input_diagnostics, paste0(opt$out, ".functional_method_input_id_diagnostics.tsv"), sep = "\t")

if (nrow(id_map) == 0L) {
  warning("No usable protein_id-to-gene_id map was found. Method summaries and functional tracks will be skipped.")
} else {
  msg("[Methods] Mapping functional annotation proteins to GFF genes...")
  method_gene_tables <- lapply(names(method_protein_ids), function(nm) {
    make_method_gene_table(nm, method_protein_ids[[nm]])
  })
  method_gene_list_dt <- rbindlist(method_gene_tables, use.names = TRUE, fill = TRUE)
  if (nrow(method_gene_list_dt) > 0L) {
    fwrite(method_gene_list_dt, paste0(opt$out, ".functional_method_gene_list.tsv"), sep = "\t")
  }

  method_summary_dt <- rbindlist(lapply(names(method_protein_ids), function(nm) {
    tab <- method_gene_list_dt[method == nm]
    data.table(
      method = nm,
      input_proteins_or_queries = length(unique(clean_id(method_protein_ids[[nm]]))),
      genes_mapped_to_gff = uniqueN(tab[!is.na(gene_id) & gene_id != "", gene_id]),
      genes_on_plotted_contigs = uniqueN(tab[on_plotted_contig == TRUE & !is.na(gene_id) & gene_id != "", gene_id]),
      plotted_fraction = ifelse(uniqueN(tab[!is.na(gene_id) & gene_id != "", gene_id]) > 0,
                                uniqueN(tab[on_plotted_contig == TRUE & !is.na(gene_id) & gene_id != "", gene_id]) /
                                  uniqueN(tab[!is.na(gene_id) & gene_id != "", gene_id]),
                                NA_real_)
    )
  }), use.names = TRUE, fill = TRUE)
  fwrite(method_summary_dt, paste0(opt$out, ".functional_method_summary.tsv"), sep = "\t")

  if (isTRUE(opt$include_functional_tracks)) {
    msg("[Tracks] Computing optional functional annotation tracks...")
    functional_inputs <- list(
      "eggNOG-supported genes / Mb" = protein_ids_to_gene_coords(method_protein_ids[["EggNOG"]]),
      "InterPro-supported genes / Mb" = protein_ids_to_gene_coords(method_protein_ids[["InterProScan"]]),
      "Swiss-Prot-supported genes / Mb" = protein_ids_to_gene_coords(method_protein_ids[["DIAMOND Swiss-Prot"]]),
      "KEGG KO-supported genes / Mb" = protein_ids_to_gene_coords(method_protein_ids[["KofamScan KEGG KO"]]),
      "dbCAN CAZyme genes / Mb" = protein_ids_to_gene_coords(method_protein_ids[["dbCAN CAZyme"]])
    )
    method_to_track <- c(
      "EggNOG" = "eggNOG-supported genes / Mb",
      "InterProScan" = "InterPro-supported genes / Mb",
      "DIAMOND Swiss-Prot" = "Swiss-Prot-supported genes / Mb",
      "KofamScan KEGG KO" = "KEGG KO-supported genes / Mb",
      "dbCAN CAZyme" = "dbCAN CAZyme genes / Mb"
    )
    for (nm in names(functional_inputs)) {
      if (nrow(functional_inputs[[nm]]) > 0L) {
        track_list[[nm]] <- make_count_track(functional_inputs[[nm]], nm, "Functional")
        functional_gr_list[[nm]] <- dt_to_gr(functional_inputs[[nm]], seq_col = "plot_label")
      }
    }
    for (method_nm in names(method_to_track)) {
      tr <- method_to_track[[method_nm]]
      if (tr %in% names(functional_gr_list)) method_gr_list[[method_nm]] <- functional_gr_list[[tr]]
    }
  }
}

# ==============================================================================
# Map tracks to global coordinates
# ==============================================================================

map_to_global <- function(dt) {
  dt2 <- merge(dt, seq_df[, .(seqid, offset, plot_label)], by = "seqid", all.x = TRUE, suffixes = c("", ".seq"))
  if ("plot_label.seq" %in% names(dt2)) {
    dt2[, plot_label := fifelse(is.na(plot_label), plot_label.seq, plot_label)]
    dt2[, plot_label.seq := NULL]
  }
  if (any(is.na(dt2$offset))) {
    bad <- unique(dt2[is.na(offset), seqid])
    stop("Some seqids are missing offsets: ", paste(bad, collapse = ", "))
  }
  dt2[, `:=`(
    xmin = offset + (as.numeric(start) - 1),
    xmax = offset + as.numeric(end),
    midpoint = (as.numeric(start) + as.numeric(end)) / 2,
    midpoint_mb = ((as.numeric(start) + as.numeric(end)) / 2) / 1e6
  )]
  dt2[]
}

track_long <- rbindlist(track_list, use.names = TRUE, fill = TRUE)
track_long_g <- map_to_global(track_long)

nonzero_tracks <- track_long_g[, .(max_value = max(value, na.rm = TRUE)), by = track]
nonzero_tracks <- nonzero_tracks[is.finite(max_value) & (max_value > 0 | track == "GC fraction"), track]
track_long_g <- track_long_g[track %in% nonzero_tracks]

preferred_order <- c(
  "All genes / Mb",
  "Protein-coding genes / Mb",
  "Non-coding genes / Mb",
  "Pseudogenes / Mb",
  "CDS / Mb",
  "CDS bp fraction",
  "Corrected/problem CDS / Mb",
  "eggNOG-supported genes / Mb",
  "InterPro-supported genes / Mb",
  "Swiss-Prot-supported genes / Mb",
  "KEGG KO-supported genes / Mb",
  "dbCAN CAZyme genes / Mb",
  "GC fraction"
)
track_order <- c(preferred_order[preferred_order %in% unique(track_long_g$track)],
                 setdiff(unique(track_long_g$track), preferred_order))
track_long_g[, track := factor(track, levels = track_order)]

fwrite(track_long_g[
  , .(seqid, plot_label, start, end, xmin, xmax, midpoint_mb, track, group, count, value)
], paste0(opt$out, ".window_tracks.tsv"), sep = "\t")

gene_track <- track_long_g[as.character(track) == "All genes / Mb"]
gene_rich_gr <- GRanges()
gene_poor_gr <- GRanges()

if (nrow(gene_track) > 0L) {
  gene_track[, gene_density_z := as.numeric(scale(value))]
  gene_track[, gene_density_class := fifelse(
    gene_density_z >= opt$gene_rich_z, "gene-rich",
    fifelse(gene_density_z <= -opt$gene_rich_z, "gene-poor", "average")
  )]
  outliers <- gene_track[
    gene_density_class != "average",
    .(seqid, plot_label, start, end, value, gene_density_z, gene_density_class)
  ]
  fwrite(outliers, paste0(opt$out, ".gene_density_outlier_windows.tsv"), sep = "\t")

  if (nrow(outliers[gene_density_class == "gene-rich"]) > 0L) {
    gene_rich_gr <- dt_to_gr(outliers[gene_density_class == "gene-rich"], seq_col = "plot_label")
  }
  if (nrow(outliers[gene_density_class == "gene-poor"]) > 0L) {
    gene_poor_gr <- dt_to_gr(outliers[gene_density_class == "gene-poor"], seq_col = "plot_label")
  }
}

gc_track <- track_long_g[as.character(track) == "GC fraction"]
gc_gr <- if (nrow(gc_track) > 0L) {
  GRanges(seqnames = gc_track$plot_label, ranges = IRanges(gc_track$start, gc_track$end), gc = gc_track$value)
} else {
  GRanges()
}

problem_windows_dt <- track_long_g[as.character(track) == "Corrected/problem CDS / Mb" & value > 0]
problem_windows_gr <- if (nrow(problem_windows_dt) > 0L) {
  GRanges(seqnames = problem_windows_dt$plot_label,
          ranges = IRanges(problem_windows_dt$start, problem_windows_dt$end),
          score = problem_windows_dt$value)
} else {
  GRanges()
}

# ==============================================================================
# Feature-count summary
# ==============================================================================

parse_feature_counts <- function(path) {
  if (!safe_file(path)) return(data.table())
  txt <- readLines(path, warn = FALSE)
  txt <- trimws(txt)
  txt <- txt[nzchar(txt)]
  txt <- txt[grepl(":", txt)]
  txt <- txt[!grepl("^Overall Counts", txt)]
  key <- sub(":.*$", "", txt)
  val <- sub("^.*:\\s*", "", txt)
  val <- gsub(",", "", val)
  suppressWarnings(num <- as.numeric(val))
  dt <- data.table(feature = key, count = num)
  dt <- dt[!is.na(count)]
  dt[, feature := trimws(feature)]
  dt[]
}

get_feature_count_value <- function(dt, feature_name) {
  if (is.null(dt) || nrow(dt) == 0L) return(NA_real_)
  x <- dt[feature == feature_name, count]
  if (length(x) == 0L) NA_real_ else x[1L]
}

make_gff_feature_counts <- function() {
  transcript_types <- c("mRNA", "transcript")
  ncrna_types <- c("lnc_RNA", "ncRNA", "tRNA", "rRNA", "miRNA", "snoRNA", "snRNA")

  classified_gene_ids <- unique(c(genes_pc$gene_id, genes_ncrna$gene_id, genes_pseudo$gene_id))
  classified_gene_ids <- classified_gene_ids[!is.na(classified_gene_ids) & classified_gene_ids != ""]
  other_genes <- genes_broad[!(gene_id %in% classified_gene_ids)]

  out <- data.table(
    feature = c(
      "genes",
      "genes (protein coding)",
      "genes (non coding)",
      "genes (pseudogene-like)",
      "genes (other/unclassified)",
      "mRNAs/transcripts",
      "non-coding RNA features",
      "CDSs",
      "CDSs (corrected/problematic)"
    ),
    count = c(
      nrow(genes_broad),
      nrow(genes_pc),
      nrow(genes_ncrna),
      nrow(genes_pseudo),
      nrow(other_genes),
      nrow(gff[type %in% transcript_types]),
      nrow(gff[type %in% ncrna_types]),
      nrow(cds_all),
      nrow(cds_problem)
    ),
    source = "GFF plotted contigs"
  )

  out[, feature_group := fifelse(
    grepl("^genes", feature), "Genes",
    fifelse(grepl("mRNAs|transcripts", feature), "mRNAs/transcripts",
    fifelse(grepl("^CDSs", feature), "CDSs",
    fifelse(grepl("non-coding RNA features", feature), "ncRNA features", "Other")))
  )]
  out[]
}

make_annotation_count_diagnostics <- function(gff_counts_dt, egapx_counts_dt) {
  data.table(
    metric = c(
      "GFF genes on plotted contigs",
      "GFF protein-coding genes on plotted contigs",
      "GFF non-coding genes on plotted contigs",
      "GFF pseudogene-like genes on plotted contigs",
      "GFF CDS features on plotted contigs",
      "GFF corrected/problematic CDS features on plotted contigs",
      "EGAPx feature_counts.txt genes",
      "EGAPx feature_counts.txt protein-coding genes",
      "EGAPx feature_counts.txt non-coding genes",
      "EGAPx feature_counts.txt non-transcribed pseudogenes",
      "EGAPx feature_counts.txt CDSs"
    ),
    count = c(
      get_feature_count_value(gff_counts_dt, "genes"),
      get_feature_count_value(gff_counts_dt, "genes (protein coding)"),
      get_feature_count_value(gff_counts_dt, "genes (non coding)"),
      get_feature_count_value(gff_counts_dt, "genes (pseudogene-like)"),
      get_feature_count_value(gff_counts_dt, "CDSs"),
      get_feature_count_value(gff_counts_dt, "CDSs (corrected/problematic)"),
      get_feature_count_value(egapx_counts_dt, "genes"),
      get_feature_count_value(egapx_counts_dt, "genes (protein coding)"),
      get_feature_count_value(egapx_counts_dt, "genes (non coding)"),
      get_feature_count_value(egapx_counts_dt, "genes (non-transcribed pseudo)"),
      get_feature_count_value(egapx_counts_dt, "CDSs")
    ),
    source = c(rep("GFF plotted contigs", 6), rep("EGAPx feature_counts.txt", 5))
  )
}

feature_counts_dt <- parse_feature_counts(opt$feature_counts)
if (nrow(feature_counts_dt) > 0L) {
  fwrite(feature_counts_dt, paste0(opt$out, ".feature_counts.parsed.tsv"), sep = "\t")
}

gff_feature_counts_dt <- make_gff_feature_counts()
fwrite(gff_feature_counts_dt, paste0(opt$out, ".gff_feature_counts.tsv"), sep = "\t")

annotation_count_diagnostics <- make_annotation_count_diagnostics(gff_feature_counts_dt, feature_counts_dt)
fwrite(annotation_count_diagnostics, paste0(opt$out, ".annotation_count_diagnostics.tsv"), sep = "\t")

make_standardized_count_comparison <- function(gff_counts_dt, egapx_counts_dt) {
  # This table intentionally separates count source, scope, and counting unit.
  # EGAPx feature_counts.txt is the official genome-wide summary.
  # GFF counts are coordinate-backed and are therefore used for spatial plots.
  egapx <- data.table(
    feature_class = c("Total genes", "Protein-coding genes", "Non-coding genes",
                      "Pseudogenes", "Transcripts/mRNAs", "CDSs"),
    count = c(
      get_feature_count_value(egapx_counts_dt, "genes"),
      get_feature_count_value(egapx_counts_dt, "genes (protein coding)"),
      get_feature_count_value(egapx_counts_dt, "genes (non coding)"),
      get_feature_count_value(egapx_counts_dt, "genes (non-transcribed pseudo)"),
      get_feature_count_value(egapx_counts_dt, "mRNAs"),
      get_feature_count_value(egapx_counts_dt, "CDSs")
    ),
    source = "EGAPx feature_counts.txt",
    scope = "genome-wide",
    counting_unit = c("gene models", "gene models", "gene models",
                      "non-transcribed pseudogene models", "mRNA models", "CDS models")
  )

  gff <- data.table(
    feature_class = c("Total genes", "Protein-coding genes", "Non-coding genes",
                      "Pseudogene-like features", "Transcript-like features", "CDS intervals"),
    count = c(
      get_feature_count_value(gff_counts_dt, "genes"),
      get_feature_count_value(gff_counts_dt, "genes (protein coding)"),
      get_feature_count_value(gff_counts_dt, "genes (non coding)"),
      get_feature_count_value(gff_counts_dt, "genes (pseudogene-like)"),
      get_feature_count_value(gff_counts_dt, "mRNAs/transcripts"),
      get_feature_count_value(gff_counts_dt, "CDSs")
    ),
    source = "GFF-derived plotted-contig features",
    scope = sprintf("top %d contigs; %.2f Gb", nrow(seq_df), sum(seq_df$length) / 1e9),
    counting_unit = c("GFF gene feature rows", "GFF gene feature rows", "GFF gene feature rows",
                      "GFF pseudogene-like rows detected by parser", "GFF transcript-like rows", "GFF CDS interval rows")
  )

  out <- rbindlist(list(egapx, gff), use.names = TRUE, fill = TRUE)
  out <- out[!is.na(count)]
  out[, feature_class := factor(feature_class, levels = c(
    "Total genes", "Protein-coding genes", "Non-coding genes", "Pseudogenes",
    "Pseudogene-like features", "Transcripts/mRNAs", "Transcript-like features",
    "CDSs", "CDS intervals"
  ))]
  out[]
}

standardized_count_comparison_dt <- make_standardized_count_comparison(gff_feature_counts_dt, feature_counts_dt)
fwrite(standardized_count_comparison_dt, paste0(opt$out, ".annotation_count_comparison_standardized.tsv"), sep = "\t")

# Official count values for plot labels. These are non-spatial summaries; spatial tracks still use GFF coordinates.
if (tolower(opt$official_count_source) == "egapx" && nrow(feature_counts_dt) > 0L) {
  official_count_source_label <- "EGAPx feature_counts.txt (genome-wide)"
  official_genes_total <- get_feature_count_value(feature_counts_dt, "genes")
  official_pc_total <- get_feature_count_value(feature_counts_dt, "genes (protein coding)")
  official_nc_total <- get_feature_count_value(feature_counts_dt, "genes (non coding)")
  official_pseudo_total <- get_feature_count_value(feature_counts_dt, "genes (non-transcribed pseudo)")
  official_mrna_total <- get_feature_count_value(feature_counts_dt, "mRNAs")
  official_cds_total <- get_feature_count_value(feature_counts_dt, "CDSs")
} else {
  official_count_source_label <- "GFF-derived plotted-contig features"
  official_genes_total <- get_feature_count_value(gff_feature_counts_dt, "genes")
  official_pc_total <- get_feature_count_value(gff_feature_counts_dt, "genes (protein coding)")
  official_nc_total <- get_feature_count_value(gff_feature_counts_dt, "genes (non coding)")
  official_pseudo_total <- get_feature_count_value(gff_feature_counts_dt, "genes (pseudogene-like)")
  official_mrna_total <- get_feature_count_value(gff_feature_counts_dt, "mRNAs/transcripts")
  official_cds_total <- get_feature_count_value(gff_feature_counts_dt, "CDSs")
}

fmt_count <- function(x) {
  ifelse(is.na(x), "NA", format(x, scientific = FALSE, big.mark = ","))
}

# ==============================================================================
# Color setup
# ==============================================================================

track_colors <- c(
  "All genes / Mb" = "#3288BD",
  "Protein-coding genes / Mb" = "#D53E4F",
  "Non-coding genes / Mb" = "#1B9E77",
  "Pseudogenes / Mb" = "#666666",
  "CDS / Mb" = "#F46D43",
  "CDS bp fraction" = "#F46D43",
  "Corrected/problem CDS / Mb" = "#7F0000",
  "eggNOG-supported genes / Mb" = "#1B9E77",
  "InterPro-supported genes / Mb" = "#7570B3",
  "Swiss-Prot-supported genes / Mb" = "#E7298A",
  "KEGG KO-supported genes / Mb" = "#E6AB02",
  "dbCAN CAZyme genes / Mb" = "#A6761D",
  "GC fraction" = "#5E4FA2"
)

method_colors <- c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#A6761D", "#666666")


track_palette_map <- c(
  "All genes / Mb" = opt$gene_palette,
  "Protein-coding genes / Mb" = opt$coding_palette,
  "Non-coding genes / Mb" = opt$noncoding_palette,
  "Pseudogenes / Mb" = opt$pseudo_palette,
  "CDS / Mb" = opt$cds_palette,
  "CDS bp fraction" = opt$cds_palette,
  "Corrected/problem CDS / Mb" = opt$problem_palette,
  "eggNOG-supported genes / Mb" = "Greens",
  "InterPro-supported genes / Mb" = "Purples",
  "Swiss-Prot-supported genes / Mb" = "RdPu",
  "KEGG KO-supported genes / Mb" = "YlOrBr",
  "dbCAN CAZyme genes / Mb" = "Oranges",
  "GC fraction" = opt$gc_palette
)

# ==============================================================================
# Circular plot
# ==============================================================================

if (isTRUE(opt$circular)) {
  is_dark <- tolower(opt$theme) == "dark"
  col_bg <- if (is_dark) "#0B0F14" else "white"
  col_fg <- if (is_dark) "#E7EEF7" else "#111111"
  col_grid <- if (is_dark) "#2A3440" else "#D0D0D0"

  ideo <- seq_df[, .(seqid, plot_label, xmin = start_global, xmax = end_global, length = length)]
  ideo[, y0 := 3.80]
  ideo[, y1 := 4.35]
  ideo[, grp := (.I %% 2)]
  ideo[, fill := ifelse(grp == 0, "#BFD7EA", "#F7E1AE")]

  bounds <- seq_df[, .(x = start_global)]
  lab_k <- min(opt$label_top_k, nrow(seq_df))
  labels <- copy(seq_df[seq_len(lab_k)])
  labels[, x := as.numeric((start_global + end_global) / 2)]

  make_ticks <- function(seqid, seqlen_int, offset_num, major, minor) {
    major_at <- seq.int(0L, seqlen_int, by = major)
    minor_at <- seq.int(0L, seqlen_int, by = minor)
    major_dt <- data.table(seqid = seqid, pos = major_at, kind = "major")
    minor_dt <- data.table(seqid = seqid, pos = setdiff(minor_at, major_at), kind = "minor")
    dt <- rbindlist(list(major_dt, minor_dt), use.names = TRUE)
    dt[, x := as.numeric(offset_num) + as.numeric(pos)]
    dt
  }

  ticks <- rbindlist(lapply(seq_len(nrow(seq_df)), function(i) {
    make_ticks(seq_df$seqid[i], seq_df$length_int[i], seq_df$offset[i],
               as.integer(opt$major_tick_bp), as.integer(opt$minor_tick_bp))
  }), use.names = TRUE)

  sep <- 0.14
  track_height <- 0.78
  r_start <- 4.70
  track_radii <- data.table(track = factor(track_order, levels = track_order))
  track_radii[, idx := .I]
  track_radii[, y0 := r_start + (idx - 1) * (track_height + sep)]
  track_radii[, y1 := y0 + track_height]
  r_label <- max(track_radii$y1) + 0.60

  track_plot <- merge(track_long_g, track_radii[, .(track, y0, y1)], by = "track", all.x = TRUE)
  track_plot[, value_clip := clip_upper(value, opt$clip_quant), by = track]
  track_plot[, ybot := y0]
  track_plot[, ytop := rescale_to_ring(value, y0[1], y1[1]), by = track]

  track_plot[, fill := {
    pal_name <- track_palette_map[as.character(track[1])]
    if (is.na(pal_name)) pal_name <- opt$palette
    pal <- get_palette(pal_name, 256)
    domain <- range(value_clip, finite = TRUE)
    if (!is.finite(domain[1]) || !is.finite(domain[2]) || domain[1] == domain[2]) {
      rep(pal[128], .N)
    } else {
      col_fun <- scales::col_numeric(palette = pal, domain = domain)
      col_fun(value_clip)
    }
  }, by = track]

  track_plot[, alpha := scale_to_alpha(value_clip, opt$alpha_min, opt$alpha_max), by = track]

  ticks[, y0 := 4.38]
  ticks[kind == "major", y1 := 4.58]
  ticks[kind == "minor", y1 := 4.50]

  tick_labels <- ticks[kind == "major"]
  tick_labels[, label := paste0(format(pos / 1e6, trim = TRUE), " Mb")]
  tick_labels[, y := 4.70]
  labels[, y := r_label]

  center_txt <- sprintf(
    paste0(
      "Official genes: %s\n",
      "Official protein-coding: %s\n",
      "Count source: %s\n",
      "Plotted contigs: %d\n",
      "Plotted bp: %s\n",
      "GFF genes plotted: %s\n",
      "GFF CDS intervals: %s\n",
      "GC: %s%%\n",
      "Window: %s bp"
    ),
    fmt_count(official_genes_total),
    fmt_count(official_pc_total),
    official_count_source_label,
    length(keep_seq),
    format(sum(seq_df$length), scientific = FALSE, big.mark = ","),
    format(nrow(genes_broad), scientific = FALSE, big.mark = ","),
    format(nrow(cds_all), scientific = FALSE, big.mark = ","),
    ifelse(is.na(genome_gc), "NA", sprintf("%.2f", 100 * genome_gc)),
    format(window, scientific = FALSE, big.mark = ",")
  )

  msg("[Plot] Building circular plot...")

  p_circular <- ggplot() +
    geom_rect(data = ideo, aes(xmin = xmin, xmax = xmax, ymin = y0, ymax = y1, fill = fill),
              color = if (is_dark) "#1B2430" else "grey25", linewidth = 0.18) +
    {if (isTRUE(opt$show_boundaries)) geom_segment(
      data = bounds, aes(x = x, xend = x, y = 3.8, yend = r_label),
      linewidth = 0.12, color = col_grid, alpha = 0.35
    ) else NULL} +
    geom_rect(data = track_plot,
              aes(xmin = xmin, xmax = xmax, ymin = ybot, ymax = ytop, fill = fill, alpha = alpha),
              color = NA) +
    geom_segment(data = ticks, aes(x = x, xend = x, y = y0, yend = y1),
                 linewidth = 0.20, color = col_fg, alpha = 0.80) +
    {if (isTRUE(opt$show_tick_labels) && nrow(tick_labels) > 0L) geom_text(
      data = tick_labels, aes(x = x, y = y, label = label), size = 2.0, color = col_fg, alpha = 0.90
    ) else NULL} +
    geom_text(data = labels, aes(x = x, y = y, label = plot_label), size = 2.6, angle = 90, vjust = 0, color = col_fg) +
    geom_text(data = track_radii, aes(x = global_len * 0.01, y = (y0 + y1) / 2, label = as.character(track)),
              inherit.aes = FALSE, hjust = 0, size = 2.5, fontface = "bold", color = col_fg) +
    annotate("text", x = global_len / 2, y = 0, label = center_txt,
             size = 4.0, fontface = "bold", lineheight = 1.08, color = col_fg) +
    coord_polar(start = 0, direction = 1, clip = "off") +
    scale_x_continuous(limits = c(0, global_len), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, r_label + 0.9), expand = c(0, 0)) +
    scale_fill_identity() +
    scale_alpha_identity() +
    labs(title = "EGAPx genome annotation landscape",
         subtitle = sprintf(
           "Official summary: %s; spatial tracks: GFF-derived features on top %d contigs (%.2f Gb); window = %s bp",
           official_count_source_label, nrow(seq_df), sum(seq_df$length) / 1e9, comma(window)
         )) +
    theme_void(base_size = 12) +
    theme(plot.background = element_rect(fill = col_bg, color = NA),
          panel.background = element_rect(fill = col_bg, color = NA),
          plot.title = element_text(hjust = 0.5, face = "bold", color = col_fg),
          plot.subtitle = element_text(hjust = 0.5, color = col_fg),
          plot.margin = margin(18, 18, 18, 18))

  if (isTRUE(opt$show_legend)) {
    p_circular <- p_circular + annotate(
      "text", x = global_len * 0.50, y = r_label + 0.50,
      label = "Darker / taller bars indicate higher window-level values within each track; colors are scaled independently by track.",
      size = 3.1, color = col_fg
    )
  }

  save_plot_all(p_circular, opt$out, width = opt$width, height = opt$height, dpi = opt$dpi, bg = col_bg)
}

# ==============================================================================
# Linear genome-track plot
# ==============================================================================

if (isTRUE(opt$linear)) {
  msg("[Plot] Building linear faceted plot...")

  linear_dt <- copy(track_long_g)
  linear_dt[, track_chr := as.character(track)]
  linear_dt[, value_plot := fifelse(
    track_chr %in% c("CDS / Mb", "Corrected/problem CDS / Mb", "Pseudogenes / Mb"),
    log10(value + 1), value)]
  linear_dt[, track_label := fifelse(
    track_chr %in% c("CDS / Mb", "Corrected/problem CDS / Mb", "Pseudogenes / Mb"),
    paste0(track_chr, " [log10(x + 1)]"), track_chr)]
  linear_dt[, plot_label := factor(plot_label, levels = seq_df$plot_label)]

  p_linear <- ggplot(linear_dt, aes(x = midpoint_mb, y = value_plot, fill = track_chr, color = track_chr)) +
    geom_area(alpha = 0.80, linewidth = 0.10) +
    facet_grid(track_label ~ plot_label, scales = "free", space = "free_x") +
    scale_fill_manual(values = track_colors, guide = "none") +
    scale_color_manual(values = track_colors, guide = "none") +
    labs(title = "Linear EGAPx annotation tracks across plotted contigs",
         subtitle = "Window-level values; selected skewed tracks are log10-transformed for readability",
         x = "Position within contig (Mb)", y = "Window value") +
    theme_bw(base_size = 10) +
    theme(strip.text.x = element_text(face = "bold", size = 8),
          strip.text.y = element_text(face = "bold", size = 8),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          axis.text.y = element_text(size = 6),
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))

  save_plot_all(p_linear, paste0(opt$out, ".linear_tracks"),
                width = opt$linear_width, height = opt$linear_height, dpi = opt$dpi, bg = "white")
}

# ==============================================================================
# Summary plots
# ==============================================================================

if (isTRUE(opt$summary_plots)) {
  msg("[Plot] Building summary plots...")

  density_dt <- copy(track_long_g)
  density_dt <- density_dt[is.finite(value)]
  density_dt[, track_chr := as.character(track)]
  density_dt[, value_density := fifelse(
    track_chr %in% c("CDS / Mb", "Corrected/problem CDS / Mb", "Pseudogenes / Mb"),
    log10(value + 1), value)]
  density_dt[, x_label := fifelse(
    track_chr %in% c("CDS / Mb", "Corrected/problem CDS / Mb", "Pseudogenes / Mb"),
    paste0(track_chr, " [log10(x + 1)]"), track_chr)]

  p_density <- ggplot(density_dt, aes(x = value_density, fill = track_chr)) +
    geom_histogram(bins = 60, color = "white", linewidth = 0.15, alpha = 0.90) +
    facet_wrap(~ x_label, scales = "free", ncol = 2) +
    scale_fill_manual(values = track_colors, guide = "none") +
    labs(title = "Distribution of window-level EGAPx annotation metrics",
         subtitle = sprintf("Window size = %s bp; plotted contigs = %d; total plotted length = %.2f Gb",
                            comma(window), nrow(seq_df), sum(seq_df$length) / 1e9),
         x = "Window value", y = "Number of windows") +
    theme_bw(base_size = 11) +
    theme(strip.text = element_text(face = "bold"), panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

  save_plot_all(p_density, paste0(opt$out, ".window_density_distributions"),
                width = 12, height = max(8, 1.3 * length(track_order)), dpi = opt$dpi, bg = "white")

  if (nrow(gff_feature_counts_dt) > 0L) {
    gff_plot <- copy(gff_feature_counts_dt)
    gff_feature_order <- c(
      "genes",
      "genes (protein coding)",
      "genes (non coding)",
      "genes (pseudogene-like)",
      "genes (other/unclassified)",
      "mRNAs/transcripts",
      "non-coding RNA features",
      "CDSs",
      "CDSs (corrected/problematic)"
    )
    gff_plot[, feature := factor(feature, levels = rev(gff_feature_order[gff_feature_order %in% feature]))]
    gff_plot[, feature_group := factor(feature_group, levels = c("Genes", "mRNAs/transcripts", "CDSs", "ncRNA features", "Other"))]

    gff_group_cols <- c(
      "Genes" = "#3288BD",
      "mRNAs/transcripts" = "#66C2A5",
      "CDSs" = "#FEE08B",
      "ncRNA features" = "#ABDDA4",
      "Other" = "#999999"
    )

    p_gff_fc <- ggplot(gff_plot, aes(x = feature, y = count, fill = feature_group)) +
      geom_col(width = 0.78, color = "grey20", linewidth = 0.15) +
      geom_text(aes(label = scales::comma(count)), hjust = -0.08, size = 3.0) +
      coord_flip(clip = "off") +
      scale_fill_manual(values = gff_group_cols, name = "Feature class") +
      scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.16))) +
      labs(
        title = "GFF-derived annotation feature summary",
        subtitle = sprintf("Primary count source for plots; counts are from the GFF after contig filtering: top %d contigs, %.2f Gb plotted", nrow(seq_df), sum(seq_df$length) / 1e9),
        x = NULL,
        y = "Feature count"
      ) +
      theme_bw(base_size = 11) +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "right"
      )

    save_plot_all(p_gff_fc, paste0(opt$out, ".annotation_feature_counts_gff"),
                  width = 10.5, height = 6.8, dpi = opt$dpi, bg = "white")
  }

  if (nrow(feature_counts_dt) > 0L) {
    keep_features <- c(
      "genes", "genes (protein coding)", "genes (non coding)", "genes (non-transcribed pseudo)",
      "genes (Ig TCR segment)", "mRNAs", "mRNAs (fully supported)", "non-coding RNAs",
      "non-coding RNAs (fully supported)", "pseudo transcripts", "CDSs", "CDSs (fully supported)",
      "CDSs (major correction)", "CDSs (premature stop)", "CDSs (has frameshifts)"
    )
    fc_plot <- feature_counts_dt[feature %in% keep_features]
    fc_plot[, feature_group := fifelse(
      grepl("^genes", feature), "Genes",
      fifelse(grepl("^mRNAs", feature), "mRNAs",
      fifelse(grepl("^CDSs", feature), "CDSs",
      fifelse(grepl("non-coding RNAs", feature), "ncRNAs",
      fifelse(grepl("pseudo", feature), "Pseudogenes", "Other"))))
    )]
    fc_plot[, feature_group := factor(feature_group, levels = c("Genes", "mRNAs", "CDSs", "ncRNAs", "Pseudogenes", "Other"))]
    fc_plot[, feature := factor(feature, levels = rev(keep_features[keep_features %in% feature]))]

    feature_group_cols <- c("Genes" = "#3288BD", "mRNAs" = "#66C2A5", "CDSs" = "#FEE08B",
                            "ncRNAs" = "#ABDDA4", "Pseudogenes" = "#D53E4F", "Other" = "#999999")

    p_fc <- ggplot(fc_plot, aes(x = feature, y = count, fill = feature_group)) +
      geom_col(width = 0.78, color = "grey20", linewidth = 0.15) +
      geom_text(aes(label = scales::comma(count)), hjust = -0.08, size = 3.0) +
      coord_flip(clip = "off") +
      scale_fill_manual(values = feature_group_cols, name = "Feature class") +
      scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.14))) +
      labs(title = "EGAPx feature_counts.txt annotation summary",
           subtitle = "Comparison-only summary parsed from EGAPx feature_counts.txt; all spatial plots use GFF-derived counts",
           x = NULL, y = "Feature count") +
      theme_bw(base_size = 11) +
      theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
            plot.title = element_text(face = "bold", hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5), legend.position = "right")

    save_plot_all(p_fc, paste0(opt$out, ".annotation_feature_counts_egapx"),
                  width = 10, height = 7, dpi = opt$dpi, bg = "white")
  }

  if (isTRUE(opt$plot_count_comparison) && nrow(standardized_count_comparison_dt) > 0L) {
    comp_plot <- copy(standardized_count_comparison_dt)
    comp_plot <- comp_plot[feature_class %in% c(
      "Total genes", "Protein-coding genes", "Non-coding genes",
      "Pseudogenes", "Pseudogene-like features", "Transcripts/mRNAs",
      "Transcript-like features", "CDSs", "CDS intervals"
    )]
    comp_plot[, feature_class := factor(as.character(feature_class), levels = rev(c(
      "Total genes", "Protein-coding genes", "Non-coding genes",
      "Pseudogenes", "Pseudogene-like features", "Transcripts/mRNAs",
      "Transcript-like features", "CDSs", "CDS intervals"
    )))]

    p_comp <- ggplot(comp_plot, aes(x = feature_class, y = count, fill = source)) +
      geom_col(position = position_dodge2(width = 0.78, preserve = "single"),
               width = 0.70, color = "grey25", linewidth = 0.12) +
      geom_text(aes(label = comma(count)),
                position = position_dodge2(width = 0.78, preserve = "single"),
                hjust = -0.08, size = 2.8) +
      coord_flip(clip = "off") +
      scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.18))) +
      scale_fill_manual(values = c(
        "EGAPx feature_counts.txt" = "#3288BD",
        "GFF-derived plotted-contig features" = "#D53E4F"
      ), name = "Count source") +
      labs(
        title = "Standardized comparison of EGAPx and GFF-derived annotation counts",
        subtitle = "EGAPx values are official genome-wide summaries; GFF values are coordinate-backed counts used for spatial plots",
        x = NULL,
        y = "Feature count"
      ) +
      theme_bw(base_size = 11) +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "bottom"
      )

    save_plot_all(p_comp, paste0(opt$out, ".annotation_feature_counts_comparison_standardized"),
                  width = 12, height = 8, dpi = opt$dpi, bg = "white")
  }
}



# ==============================================================================
# Functional method descriptive summaries
# ==============================================================================

if (isTRUE(opt$method_summary_plots) && nrow(method_summary_dt) > 0L) {
  msg("[Plot] Building functional method summary plots and tables...")

  method_plot_dt <- copy(method_summary_dt)
  method_plot_dt[, method := factor(method, levels = method[order(genes_on_plotted_contigs)])]

  p_method_counts <- ggplot(method_plot_dt, aes(x = method, y = genes_on_plotted_contigs, fill = method)) +
    geom_col(width = 0.75, color = "grey25", linewidth = 0.15, show.legend = FALSE) +
    geom_text(aes(label = comma(genes_on_plotted_contigs)), hjust = -0.08, size = 3.2) +
    coord_flip(clip = "off") +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.15))) +
    scale_fill_manual(values = setNames(method_colors[seq_len(nrow(method_plot_dt))], as.character(method_plot_dt$method))) +
    labs(
      title = "Genes identified by each functional annotation method",
      subtitle = "Counts are mapped to GFF-derived gene coordinates on the plotted contigs",
      x = NULL,
      y = "GFF genes on plotted contigs"
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))

  save_plot_all(p_method_counts, paste0(opt$out, ".functional_method_gene_counts"),
                width = 9, height = 5.5, dpi = opt$dpi, bg = "white")

  method_long <- method_summary_dt[, .(
    method,
    `All GFF genes` = genes_mapped_to_gff,
    `Plotted-contig GFF genes` = genes_on_plotted_contigs
  )]
  method_long <- melt(method_long, id.vars = "method", variable.name = "count_scope", value.name = "genes")
  method_long[, method := factor(method, levels = method_summary_dt$method[order(method_summary_dt$genes_mapped_to_gff)])]

  p_method_compare <- ggplot(method_long, aes(x = method, y = genes, fill = count_scope)) +
    geom_col(position = position_dodge(width = 0.78), width = 0.70, color = "grey25", linewidth = 0.12) +
    coord_flip() +
    scale_y_continuous(labels = comma) +
    scale_fill_brewer(palette = "Set2", name = "Count source") +
    labs(
      title = "Genome-wide versus plotted-contig functional annotation support",
      subtitle = "All counts are derived by mapping method protein IDs to GFF gene coordinates",
      x = NULL,
      y = "Number of GFF genes"
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5), legend.position = "bottom")

  save_plot_all(p_method_compare, paste0(opt$out, ".functional_method_gene_counts_comparison"),
                width = 10, height = 6, dpi = opt$dpi, bg = "white")

  if (requireNamespace("gt", quietly = TRUE)) {
    gt_tbl <- gt::gt(method_summary_dt) |>
      gt::tab_header(
        title = "Functional annotation support by method",
        subtitle = "Protein-level method outputs mapped to GFF-derived gene coordinates"
      ) |>
      gt::fmt_number(columns = c("input_proteins_or_queries", "genes_mapped_to_gff", "genes_on_plotted_contigs"), decimals = 0) |>
      gt::fmt_percent(columns = "plotted_fraction", decimals = 1)
    gt::gtsave(gt_tbl, paste0(opt$out, ".functional_method_summary.gt.html"))

    if (nrow(method_gene_list_dt) > 0L) {
      examples <- method_gene_list_dt[on_plotted_contig == TRUE][
        , head(.SD, opt$method_gene_table_top_n), by = method
      ][, .(method, gene_id, protein_id, seqid, plot_label, start, end)]
      gt_examples <- gt::gt(examples) |>
        gt::tab_header(
          title = "Example GFF genes identified by each functional annotation method",
          subtitle = sprintf("First %d plotted-contig genes per method; complete list is in the TSV output", opt$method_gene_table_top_n)
        ) |>
        gt::fmt_number(columns = c("start", "end"), decimals = 0)
      gt::gtsave(gt_examples, paste0(opt$out, ".functional_method_gene_examples.gt.html"))
    }
  }

  if (requireNamespace("gtsummary", quietly = TRUE) && nrow(method_gene_list_dt) > 0L) {
    gtsum_tbl <- method_gene_list_dt[, .(
      method,
      on_plotted_contig,
      contig = fifelse(is.na(plot_label), "Not plotted", plot_label)
    )] |>
      gtsummary::tbl_summary(by = method, include = c(on_plotted_contig, contig))
    if (requireNamespace("gt", quietly = TRUE)) {
      gt::gtsave(gtsummary::as_gt(gtsum_tbl), paste0(opt$out, ".functional_method_gene_list.gtsummary.html"))
    }
  }
}

# ==============================================================================
# karyoploteR helpers for multi-data-type tracks
# ==============================================================================

# Build simple moving mean/sd tracks for GC fraction. This is used with
# kpPlotRibbon(), following the multiple-data-type tutorial idea of overlaying a
# smoothed trend and its local variability on raw points/lines.
make_gc_smooth_dt <- function(gc_dt, k = 7) {
  if (k <= 1L || nrow(gc_dt) == 0L) return(data.table())
  if (k %% 2 == 0L) k <- k + 1L
  gc_dt <- copy(gc_dt)
  gc_dt <- gc_dt[is.finite(value)]
  if (nrow(gc_dt) == 0L) return(data.table())
  setorder(gc_dt, plot_label, start)
  out <- gc_dt[, {
    n <- .N
    half <- floor(k / 2)
    vals <- value
    means <- rep(NA_real_, n)
    sds <- rep(NA_real_, n)
    for (ii in seq_len(n)) {
      lo <- max(1L, ii - half)
      hi <- min(n, ii + half)
      vv <- vals[lo:hi]
      means[ii] <- mean(vv, na.rm = TRUE)
      sds[ii] <- stats::sd(vv, na.rm = TRUE)
    }
    data.table(
      start = start,
      end = end,
      midpoint = as.integer((start + end) / 2),
      gc_mean = means,
      gc_low = pmax(opt$gc_ymin, means - sds),
      gc_high = pmin(opt$gc_ymax, means + sds)
    )
  }, by = plot_label]
  out[is.finite(gc_mean) & is.finite(gc_low) & is.finite(gc_high)]
}

gc_smooth_dt <- make_gc_smooth_dt(gc_track, opt$gc_ribbon_k)
gc_smooth_gr <- if (nrow(gc_smooth_dt) > 0L) {
  GRanges(
    seqnames = gc_smooth_dt$plot_label,
    ranges = IRanges(gc_smooth_dt$midpoint, gc_smooth_dt$midpoint),
    gc_mean = gc_smooth_dt$gc_mean,
    gc_low = gc_smooth_dt$gc_low,
    gc_high = gc_smooth_dt$gc_high
  )
} else {
  GRanges()
}

# Strongest gene-rich windows are labelled as markers. These are not meant as
# biological gene names; they identify high-density windows that can be reported
# in supplementary tables.
gene_rich_marker_gr <- GRanges()
if (exists("gene_track") && nrow(gene_track) > 0L && "gene_density_z" %in% names(gene_track)) {
  marker_dt <- gene_track[gene_density_z >= opt$gene_rich_z]
  if (nrow(marker_dt) > 0L) {
    setorder(marker_dt, -gene_density_z)
    marker_dt <- marker_dt[seq_len(min(opt$top_marker_n, nrow(marker_dt)))]
    marker_dt[, marker_label := paste0(plot_label, ":", round(start / 1e6, 1), "Mb")]
    gene_rich_marker_gr <- GRanges(
      seqnames = marker_dt$plot_label,
      ranges = IRanges(as.integer((marker_dt$start + marker_dt$end) / 2),
                       as.integer((marker_dt$start + marker_dt$end) / 2)),
      label = marker_dt$marker_label,
      score = marker_dt$gene_density_z
    )
  }
}

all_functional_gr <- if (length(functional_gr_list) > 0L) {
  Reduce(c, functional_gr_list)
} else {
  GRanges()
}

method_colors <- c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#A6761D", "#666666")

# ==============================================================================
# karyoploteR plotting
# ==============================================================================

karyo_cols <- list(
  all_genes = track_colors[["All genes / Mb"]],
  protein_coding = track_colors[["Protein-coding genes / Mb"]],
  noncoding = track_colors[["Non-coding genes / Mb"]],
  pseudogene = track_colors[["Pseudogenes / Mb"]],
  cds = track_colors[["CDS / Mb"]],
  problem = track_colors[["Corrected/problem CDS / Mb"]],
  gc = track_colors[["GC fraction"]],
  gene_rich = "#D7191C",
  gene_poor = "#2C7BB6"
)

plot_karyoploteR_density_body <- function() {
  pp <- getDefaultPlotParams(plot.type = 2)
  pp$leftmargin <- 0.12
  pp$rightmargin <- 0.03
  pp$topmargin <- 55
  pp$bottommargin <- 65
  pp$data1height <- 900
  pp$ideogramheight <- 45

  karyo_title <- "EGAPx karyoploteR density tracks"
  karyo_subtitle <- sprintf("Top %d contigs; window = %s bp; plotted length = %.2f Gb",
                            nrow(chr_map), comma(opt$window), sum(seq_df$length) / 1e9)

  kp <- plotKaryotype(
    genome = karyo_genome,
    plot.type = 2,
    chromosomes = chr_map$chr_label,
    plot.params = pp,
    main = ""
  )
  add_karyo_outer_title(karyo_title, karyo_subtitle)
  kpAddBaseNumbers(kp, tick.dist = opt$major_tick_bp, minor.tick.dist = opt$minor_tick_bp, cex = opt$karyo_base_cex)

  tracks <- data.table(
    name = c("All genes", "Protein-coding", "Non-coding", "Pseudogenes", "CDS", "Problem CDS", "GC"),
    r0 = c(0.86, 0.72, 0.58, 0.44, 0.30, 0.16, 0.02),
    r1 = c(0.98, 0.84, 0.70, 0.56, 0.42, 0.28, 0.14)
  )

  for (i in seq_len(nrow(tracks))) {
    kpDataBackground(kp, data.panel = 1, r0 = tracks$r0[i], r1 = tracks$r1[i], col = opt$karyo_track_bg)
    kpAddLabels(kp, labels = tracks$name[i], data.panel = 1, r0 = tracks$r0[i], r1 = tracks$r1[i],
                cex = opt$karyo_label_cex, label.margin = 0.018)
  }

  kpPlotDensity(kp, data = genes_broad_gr, window.size = opt$window, col = adjustcolor(karyo_cols$all_genes, 0.75), border = NA, data.panel = 1, r0 = tracks[name == "All genes"]$r0, r1 = tracks[name == "All genes"]$r1)
  kpPlotDensity(kp, data = genes_pc_gr, window.size = opt$window, col = adjustcolor(karyo_cols$protein_coding, 0.75), border = NA, data.panel = 1, r0 = tracks[name == "Protein-coding"]$r0, r1 = tracks[name == "Protein-coding"]$r1)
  kpPlotDensity(kp, data = genes_ncrna_gr, window.size = opt$window, col = adjustcolor(karyo_cols$noncoding, 0.75), border = NA, data.panel = 1, r0 = tracks[name == "Non-coding"]$r0, r1 = tracks[name == "Non-coding"]$r1)
  kpPlotDensity(kp, data = genes_pseudo_gr, window.size = opt$window, col = adjustcolor(karyo_cols$pseudogene, 0.75), border = NA, data.panel = 1, r0 = tracks[name == "Pseudogenes"]$r0, r1 = tracks[name == "Pseudogenes"]$r1)
  kpPlotDensity(kp, data = cds_all_gr, window.size = opt$window, col = adjustcolor(karyo_cols$cds, 0.70), border = NA, data.panel = 1, r0 = tracks[name == "CDS"]$r0, r1 = tracks[name == "CDS"]$r1)
  kpPlotDensity(kp, data = cds_problem_gr, window.size = opt$window, col = adjustcolor(karyo_cols$problem, 0.80), border = NA, data.panel = 1, r0 = tracks[name == "Problem CDS"]$r0, r1 = tracks[name == "Problem CDS"]$r1)

  if (length(gc_gr) > 0L) {
    kpLines(kp, data = gc_gr, y = mcols(gc_gr)$gc, ymin = opt$gc_ymin, ymax = opt$gc_ymax,
            col = karyo_cols$gc, lwd = 1.2, data.panel = 1,
            r0 = tracks[name == "GC"]$r0, r1 = tracks[name == "GC"]$r1)
    kpAbline(kp, h = genome_gc, ymin = opt$gc_ymin, ymax = opt$gc_ymax,
             data.panel = 1, r0 = tracks[name == "GC"]$r0, r1 = tracks[name == "GC"]$r1,
             col = adjustcolor("#222222", 0.55), lty = 2)
  }

  if (isTRUE(opt$show_legend)) {
    legend("bottom", horiz = TRUE, inset = opt$karyo_legend_inset, xpd = NA, bty = "n",
           cex = opt$karyo_legend_cex,
           x.intersp = opt$karyo_legend_x_intersp,
           y.intersp = opt$karyo_legend_y_intersp,
           fill = unlist(karyo_cols[c("all_genes", "protein_coding", "noncoding", "pseudogene", "cds", "problem", "gc")]),
           legend = c("All genes", "Protein-coding", "Non-coding", "Pseudogenes", "CDS", "Problem CDS", "GC"))
  }
}

plot_karyoploteR_multitype_body <- function() {
  pp <- getDefaultPlotParams(plot.type = 2)
  pp$leftmargin <- 0.13
  pp$rightmargin <- 0.03
  pp$topmargin <- 55
  pp$bottommargin <- 70
  pp$data1height <- 760
  pp$data2height <- 480
  pp$ideogramheight <- 45

  karyo_title <- "EGAPx multiple-data-type annotation landscape"
  karyo_subtitle <- sprintf("Density, GC, gene-rich/poor regions, problem-CDS markers, and functional support; window = %s bp",
                            comma(opt$window))

  kp <- plotKaryotype(
    genome = karyo_genome,
    plot.type = 2,
    chromosomes = chr_map$chr_label,
    plot.params = pp,
    main = ""
  )
  add_karyo_outer_title(karyo_title, karyo_subtitle)
  kpAddBaseNumbers(kp, tick.dist = opt$major_tick_bp, minor.tick.dist = opt$minor_tick_bp, cex = opt$karyo_base_cex)

  # Panel 1 follows the tutorial pattern of combining rectangles, points/lines,
  # rolling summaries, and markers above the ideograms. Panel 2 shows regions
  # and coverage below the ideograms.
  tracks1 <- data.table(
    name = c("All genes density", "Protein-coding density", "CDS density", "GC fraction"),
    r0 = c(0.76, 0.55, 0.34, 0.07),
    r1 = c(0.95, 0.74, 0.53, 0.27)
  )
  tracks2 <- data.table(
    name = c("Gene-rich / gene-poor regions", "Problem CDS markers", "Functional coverage", "Functional support subtracks"),
    r0 = c(0.72, 0.48, 0.23, 0.03),
    r1 = c(0.95, 0.68, 0.43, 0.18)
  )

  for (i in seq_len(nrow(tracks1))) {
    kpDataBackground(kp, data.panel = 1, r0 = tracks1$r0[i], r1 = tracks1$r1[i], col = opt$karyo_track_bg)
    kpAddLabels(kp, labels = tracks1$name[i], data.panel = 1, r0 = tracks1$r0[i], r1 = tracks1$r1[i],
                cex = opt$karyo_label_cex, label.margin = 0.018)
  }
  for (i in seq_len(nrow(tracks2))) {
    kpDataBackground(kp, data.panel = 2, r0 = tracks2$r0[i], r1 = tracks2$r1[i], col = opt$karyo_track_bg)
    kpAddLabels(kp, labels = tracks2$name[i], data.panel = 2, r0 = tracks2$r0[i], r1 = tracks2$r1[i],
                cex = opt$karyo_label_cex, label.margin = 0.018)
  }

  # Density tracks.
  kpPlotDensity(kp, data = genes_broad_gr, window.size = opt$window,
                col = adjustcolor(karyo_cols$all_genes, 0.70), border = NA,
                data.panel = 1, r0 = tracks1[name == "All genes density"]$r0,
                r1 = tracks1[name == "All genes density"]$r1)
  kpPlotDensity(kp, data = genes_pc_gr, window.size = opt$window,
                col = adjustcolor(karyo_cols$protein_coding, 0.72), border = NA,
                data.panel = 1, r0 = tracks1[name == "Protein-coding density"]$r0,
                r1 = tracks1[name == "Protein-coding density"]$r1)
  kpPlotDensity(kp, data = cds_all_gr, window.size = opt$window,
                col = adjustcolor(karyo_cols$cds, 0.62), border = NA,
                data.panel = 1, r0 = tracks1[name == "CDS density"]$r0,
                r1 = tracks1[name == "CDS density"]$r1)

  # GC line, genome-average reference, and optional ribbon representing local mean ± sd.
  if (length(gc_gr) > 0L) {
    if (length(gc_smooth_gr) > 0L && opt$gc_ribbon_k > 1L) {
      try(kpPlotRibbon(kp, data = gc_smooth_gr,
                       y0 = mcols(gc_smooth_gr)$gc_low,
                       y1 = mcols(gc_smooth_gr)$gc_high,
                       ymin = opt$gc_ymin, ymax = opt$gc_ymax,
                       col = adjustcolor(karyo_cols$gc, 0.18), border = NA,
                       data.panel = 1,
                       r0 = tracks1[name == "GC fraction"]$r0,
                       r1 = tracks1[name == "GC fraction"]$r1), silent = TRUE)
      try(kpLines(kp, data = gc_smooth_gr,
                  y = mcols(gc_smooth_gr)$gc_mean,
                  ymin = opt$gc_ymin, ymax = opt$gc_ymax,
                  col = adjustcolor(karyo_cols$gc, 0.80), lwd = 1.4,
                  data.panel = 1,
                  r0 = tracks1[name == "GC fraction"]$r0,
                  r1 = tracks1[name == "GC fraction"]$r1), silent = TRUE)
    }
    kpLines(kp, data = gc_gr, y = mcols(gc_gr)$gc, ymin = opt$gc_ymin, ymax = opt$gc_ymax,
            col = adjustcolor(karyo_cols$gc, 0.45), lwd = 0.75,
            data.panel = 1, r0 = tracks1[name == "GC fraction"]$r0,
            r1 = tracks1[name == "GC fraction"]$r1)
    kpAbline(kp, h = genome_gc, ymin = opt$gc_ymin, ymax = opt$gc_ymax,
             data.panel = 1, r0 = tracks1[name == "GC fraction"]$r0,
             r1 = tracks1[name == "GC fraction"]$r1,
             col = adjustcolor("#333333", 0.55), lty = 2)
  }

  # Region tracks: gene-rich windows at the top of the panel, gene-poor windows at the bottom.
  if (length(gene_rich_gr) > 0L) {
    kpRect(kp, data = gene_rich_gr, y0 = 0, y1 = 1,
           col = adjustcolor(karyo_cols$gene_rich, 0.65), border = NA,
           data.panel = 2, r0 = tracks2[name == "Gene-rich / gene-poor regions"]$r0,
           r1 = tracks2[name == "Gene-rich / gene-poor regions"]$r0 + 0.10)
  }
  if (length(gene_poor_gr) > 0L) {
    kpRect(kp, data = gene_poor_gr, y0 = 0, y1 = 1,
           col = adjustcolor(karyo_cols$gene_poor, 0.65), border = NA,
           data.panel = 2, r0 = tracks2[name == "Gene-rich / gene-poor regions"]$r1 - 0.10,
           r1 = tracks2[name == "Gene-rich / gene-poor regions"]$r1)
  }

  # Label a limited number of the strongest gene-rich windows as markers, following
  # the tutorial's use of kpPlotMarkers for selected labelled positions.
  if (length(gene_rich_marker_gr) > 0L) {
    try(kpPlotMarkers(kp, data = gene_rich_marker_gr,
                      labels = mcols(gene_rich_marker_gr)$label,
                      label.color = "#333333", line.color = karyo_cols$gene_rich,
                      r0 = tracks2[name == "Gene-rich / gene-poor regions"]$r0,
                      r1 = tracks2[name == "Gene-rich / gene-poor regions"]$r1 + 0.08,
                      data.panel = 2, cex = 0.45, label.margin = 4), silent = TRUE)
  }

  # Problem-CDS windows as both point markers and vertical segments, making rare
  # windows visible even when density is low.
  if (length(problem_windows_gr) > 0L) {
    problem_mid <- start(problem_windows_gr) + width(problem_windows_gr) / 2
    problem_y <- pmin(mcols(problem_windows_gr)$score, 10)
    kpPoints(kp, chr = as.character(seqnames(problem_windows_gr)), x = problem_mid,
             y = problem_y, ymin = 0, ymax = 10,
             col = karyo_cols$problem, pch = 16, cex = 0.38,
             data.panel = 2, r0 = tracks2[name == "Problem CDS markers"]$r0,
             r1 = tracks2[name == "Problem CDS markers"]$r1)
    try(kpSegments(kp, chr = as.character(seqnames(problem_windows_gr)),
                   x0 = problem_mid, x1 = problem_mid,
                   y0 = 0, y1 = problem_y, ymin = 0, ymax = 10,
                   col = adjustcolor(karyo_cols$problem, 0.40), lwd = 0.5,
                   data.panel = 2, r0 = tracks2[name == "Problem CDS markers"]$r0,
                   r1 = tracks2[name == "Problem CDS markers"]$r1), silent = TRUE)
  }

  # Functional-support coverage and stacked subtracks.
  if (isTRUE(opt$plot_functional_coverage) && length(all_functional_gr) > 0L) {
    kpPlotCoverage(kp, data = all_functional_gr,
                   col = adjustcolor("#1B9E77", 0.55), border = NA,
                   data.panel = 2,
                   r0 = tracks2[name == "Functional coverage"]$r0,
                   r1 = tracks2[name == "Functional coverage"]$r1)
    kpAbline(kp, h = 0.5, col = adjustcolor("#333333", 0.45), lty = 2,
             data.panel = 2,
             r0 = tracks2[name == "Functional coverage"]$r0,
             r1 = tracks2[name == "Functional coverage"]$r1)
  }

  if (isTRUE(opt$plot_functional_subtracks) && length(functional_gr_list) > 0L) {
    fun_cols <- c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#A6761D")
    nfun <- length(functional_gr_list)
    idx <- 0
    for (nm in names(functional_gr_list)) {
      idx <- idx + 1
      if (length(functional_gr_list[[nm]]) == 0L) next
      sub_r0 <- tracks2[name == "Functional support subtracks"]$r0 + (idx - 1) / nfun * (tracks2[name == "Functional support subtracks"]$r1 - tracks2[name == "Functional support subtracks"]$r0)
      sub_r1 <- tracks2[name == "Functional support subtracks"]$r0 + idx / nfun * (tracks2[name == "Functional support subtracks"]$r1 - tracks2[name == "Functional support subtracks"]$r0)
      kpPlotDensity(kp, data = functional_gr_list[[nm]], window.size = opt$window,
                    col = adjustcolor(fun_cols[((idx - 1) %% length(fun_cols)) + 1], 0.75), border = NA,
                    data.panel = 2, r0 = sub_r0, r1 = sub_r1)
    }
  }

  if (isTRUE(opt$show_legend)) {
    legend("bottom", horiz = TRUE, inset = opt$karyo_legend_inset, xpd = NA, bty = "n",
           cex = opt$karyo_legend_cex,
           x.intersp = opt$karyo_legend_x_intersp,
           y.intersp = opt$karyo_legend_y_intersp,
           fill = c(karyo_cols$all_genes, karyo_cols$protein_coding, karyo_cols$cds,
                    karyo_cols$gc, karyo_cols$gene_rich, karyo_cols$gene_poor, karyo_cols$problem,
                    "#1B9E77"),
           legend = c("All genes", "Protein-coding", "CDS", "GC", "Gene-rich", "Gene-poor", "Problem CDS", "Functional support"))
  }
}


plot_karyoploteR_method_multitype_body <- function(method_name, method_gr, method_col = "#1B9E77") {
  pp <- getDefaultPlotParams(plot.type = 2)
  pp$leftmargin <- 0.13
  pp$rightmargin <- 0.03
  pp$topmargin <- 55
  pp$bottommargin <- 70
  pp$data1height <- 760
  pp$data2height <- 480
  pp$ideogramheight <- 45

  method_count <- length(unique(as.character(seqnames(method_gr)) |> paste(start(method_gr), end(method_gr))))

  karyo_title <- sprintf("%s-supported genes", method_name)
  karyo_subtitle <- sprintf("Multiple-data-type EGAPx landscape; method-supported genes, density, GC, regions, and markers; window = %s bp",
                            comma(opt$window))

  kp <- plotKaryotype(
    genome = karyo_genome,
    plot.type = 2,
    chromosomes = chr_map$chr_label,
    plot.params = pp,
    main = ""
  )
  add_karyo_outer_title(karyo_title, karyo_subtitle)
  kpAddBaseNumbers(kp, tick.dist = opt$major_tick_bp, minor.tick.dist = opt$minor_tick_bp, cex = opt$karyo_base_cex)

  tracks1 <- data.table(
    name = c("All genes density", paste0(method_name, " gene density"), "CDS density", "GC fraction"),
    r0 = c(0.76, 0.55, 0.34, 0.07),
    r1 = c(0.95, 0.74, 0.53, 0.27)
  )
  tracks2 <- data.table(
    name = c(paste0(method_name, " supported regions"), paste0(method_name, " coverage"), "Gene-rich / gene-poor regions", "Problem CDS markers"),
    r0 = c(0.72, 0.48, 0.23, 0.03),
    r1 = c(0.95, 0.68, 0.43, 0.18)
  )

  for (i in seq_len(nrow(tracks1))) {
    kpDataBackground(kp, data.panel = 1, r0 = tracks1$r0[i], r1 = tracks1$r1[i], col = opt$karyo_track_bg)
    kpAddLabels(kp, labels = tracks1$name[i], data.panel = 1, r0 = tracks1$r0[i], r1 = tracks1$r1[i],
                cex = opt$karyo_label_cex, label.margin = 0.018)
  }
  for (i in seq_len(nrow(tracks2))) {
    kpDataBackground(kp, data.panel = 2, r0 = tracks2$r0[i], r1 = tracks2$r1[i], col = opt$karyo_track_bg)
    kpAddLabels(kp, labels = tracks2$name[i], data.panel = 2, r0 = tracks2$r0[i], r1 = tracks2$r1[i],
                cex = opt$karyo_label_cex, label.margin = 0.018)
  }

  kpPlotDensity(kp, data = genes_broad_gr, window.size = opt$window,
                col = adjustcolor(karyo_cols$all_genes, 0.60), border = NA,
                data.panel = 1, r0 = tracks1[1]$r0, r1 = tracks1[1]$r1)
  kpPlotDensity(kp, data = method_gr, window.size = opt$window,
                col = adjustcolor(method_col, 0.78), border = NA,
                data.panel = 1, r0 = tracks1[2]$r0, r1 = tracks1[2]$r1)
  kpPlotDensity(kp, data = cds_all_gr, window.size = opt$window,
                col = adjustcolor(karyo_cols$cds, 0.55), border = NA,
                data.panel = 1, r0 = tracks1[3]$r0, r1 = tracks1[3]$r1)

  if (length(gc_gr) > 0L) {
    if (length(gc_smooth_gr) > 0L && opt$gc_ribbon_k > 1L) {
      try(kpPlotRibbon(kp, data = gc_smooth_gr,
                       y0 = mcols(gc_smooth_gr)$gc_low, y1 = mcols(gc_smooth_gr)$gc_high,
                       ymin = opt$gc_ymin, ymax = opt$gc_ymax,
                       col = adjustcolor(karyo_cols$gc, 0.16), border = NA,
                       data.panel = 1, r0 = tracks1[4]$r0, r1 = tracks1[4]$r1), silent = TRUE)
      try(kpLines(kp, data = gc_smooth_gr, y = mcols(gc_smooth_gr)$gc_mean,
                  ymin = opt$gc_ymin, ymax = opt$gc_ymax,
                  col = adjustcolor(karyo_cols$gc, 0.80), lwd = 1.4,
                  data.panel = 1, r0 = tracks1[4]$r0, r1 = tracks1[4]$r1), silent = TRUE)
    }
    kpLines(kp, data = gc_gr, y = mcols(gc_gr)$gc, ymin = opt$gc_ymin, ymax = opt$gc_ymax,
            col = adjustcolor(karyo_cols$gc, 0.45), lwd = 0.75,
            data.panel = 1, r0 = tracks1[4]$r0, r1 = tracks1[4]$r1)
    kpAbline(kp, h = genome_gc, ymin = opt$gc_ymin, ymax = opt$gc_ymax,
             data.panel = 1, r0 = tracks1[4]$r0, r1 = tracks1[4]$r1,
             col = adjustcolor("#333333", 0.55), lty = 2)
  }

  # Method-specific annotations: regions, points, and coverage.
  kpPlotRegions(kp, data = method_gr, col = adjustcolor(method_col, 0.45), border = NA,
                data.panel = 2, r0 = tracks2[1]$r0, r1 = tracks2[1]$r1)
  kpPlotCoverage(kp, data = method_gr, col = adjustcolor(method_col, 0.55), border = NA,
                 data.panel = 2, r0 = tracks2[2]$r0, r1 = tracks2[2]$r1)
  method_mid <- start(method_gr) + width(method_gr) / 2
  kpPoints(kp, chr = as.character(seqnames(method_gr)), x = method_mid,
           y = rep(1, length(method_gr)), ymin = 0, ymax = 1,
           col = adjustcolor(method_col, 0.75), pch = 16, cex = 0.18,
           data.panel = 2, r0 = tracks2[1]$r0, r1 = tracks2[1]$r1)

  if (length(gene_rich_gr) > 0L) {
    kpRect(kp, data = gene_rich_gr, y0 = 0, y1 = 1,
           col = adjustcolor(karyo_cols$gene_rich, 0.50), border = NA,
           data.panel = 2, r0 = tracks2[3]$r0, r1 = tracks2[3]$r0 + 0.09)
  }
  if (length(gene_poor_gr) > 0L) {
    kpRect(kp, data = gene_poor_gr, y0 = 0, y1 = 1,
           col = adjustcolor(karyo_cols$gene_poor, 0.50), border = NA,
           data.panel = 2, r0 = tracks2[3]$r1 - 0.09, r1 = tracks2[3]$r1)
  }
  if (length(problem_windows_gr) > 0L) {
    problem_mid <- start(problem_windows_gr) + width(problem_windows_gr) / 2
    problem_y <- pmin(mcols(problem_windows_gr)$score, 10)
    kpPoints(kp, chr = as.character(seqnames(problem_windows_gr)), x = problem_mid,
             y = problem_y, ymin = 0, ymax = 10,
             col = karyo_cols$problem, pch = 16, cex = 0.32,
             data.panel = 2, r0 = tracks2[4]$r0, r1 = tracks2[4]$r1)
  }

  if (isTRUE(opt$show_legend)) {
    legend("bottom", horiz = TRUE, inset = opt$karyo_legend_inset, xpd = NA, bty = "n",
           cex = opt$karyo_legend_cex,
           x.intersp = opt$karyo_legend_x_intersp,
           y.intersp = opt$karyo_legend_y_intersp,
           fill = c(karyo_cols$all_genes, method_col, karyo_cols$cds, karyo_cols$gc,
                    karyo_cols$gene_rich, karyo_cols$gene_poor, karyo_cols$problem),
           legend = c("All genes", method_name, "CDS", "GC", "Gene-rich", "Gene-poor", "Problem CDS"))
  }
}

if (isTRUE(opt$method_karyoplots) && length(method_gr_list) > 0L) {
  msg("[Plot] Building per-method karyoploteR multiple-data-type plots...")
  mi <- 0L
  for (method_nm in names(method_gr_list)) {
    mi <- mi + 1L
    gr <- method_gr_list[[method_nm]]
    if (length(gr) == 0L) next
    safe_method <- gsub("[^A-Za-z0-9]+", "_", tolower(method_nm))
    col <- method_colors[((mi - 1L) %% length(method_colors)) + 1L]
    save_karyo_all(function() plot_karyoploteR_method_multitype_body(method_nm, gr, col),
                   paste0(opt$out, ".karyoploter_method_", safe_method),
                   width = opt$karyo_width, height = opt$karyo_height, dpi = opt$dpi, bg = "white")
  }
}

if (isTRUE(opt$karyoplot_density)) {
  msg("[Plot] Building karyoploteR density-track plot...")
  save_karyo_all(plot_karyoploteR_density_body, paste0(opt$out, ".karyoploter_density"),
                 width = opt$karyo_width, height = opt$karyo_height, dpi = opt$dpi, bg = "white")
}

if (isTRUE(opt$karyoplot_multitype)) {
  msg("[Plot] Building karyoploteR multiple-data-type plot...")
  save_karyo_all(plot_karyoploteR_multitype_body, paste0(opt$out, ".karyoploter_multitype"),
                 width = opt$karyo_width, height = opt$karyo_height, dpi = opt$dpi, bg = "white")
}

if (length(warnings()) > 0) {
  print(warnings())
}

msg("[Done] Created enhanced EGAPx genome annotation graphics with prefix: %s", opt$out)
