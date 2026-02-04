#!/usr/bin/env Rscript

# 02_qc: Quality control and filtering

suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(SummarizedExperiment)
  library(Matrix)
})

pipeline_root <- Sys.getenv("PIPELINE_ROOT", getwd())
source(file.path(pipeline_root, "scripts", "helpers", "cli_utils.R"))

args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
sample_id <- args$sample_id
in_rds <- args$in_rds
out_processed <- if (!is.null(args$out_processed)) args$out_processed else file.path(getwd(), "data", "processed")
out_outputs <- if (!is.null(args$out_outputs)) args$out_outputs else file.path(getwd(), "outputs")
log_file <- args$log_file

if (!is.null(sample_id) && sample_id != "" && !is.null(in_rds) && in_rds != "") {
  stop("Provide either --sample_id or --in_rds, not both.")
}

start_time <- Sys.time()

qc_dir <- file.path(out_outputs, "qc")
plot_dir <- file.path(qc_dir, "plots")

dir.create(out_processed, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

qc_params <- list(
  min_counts = 500,
  min_genes = 200,
  max_mito = 20,
  keep_tissue_only = TRUE
)

input_files <- list.files(out_processed, pattern = "\\.01_ingest\\.rds$", full.names = TRUE)
if (!is.null(sample_id) && sample_id != "") {
  input_files <- file.path(out_processed, paste0(sample_id, ".01_ingest.rds"))
  if (!file.exists(input_files)) {
    legacy <- file.path(out_processed, paste0(sample_id, ".rds"))
    if (file.exists(legacy)) input_files <- legacy
  }
}
if (!is.null(in_rds) && in_rds != "") {
  input_files <- in_rds
}

input_files <- input_files[file.exists(input_files)]
if (length(input_files) == 0) {
  legacy <- list.files(out_processed, pattern = "^[0-9]+\\.rds$", full.names = TRUE)
  input_files <- legacy
}
if (length(input_files) == 0) stop("No input .rds files found for 02_qc.")

get_sample_id <- function(path) sub("(\\.01_ingest|)\\.rds$", "", basename(path))

get_gene_symbols <- function(spe) {
  rd <- SummarizedExperiment::rowData(spe)
  for (col in c("gene_name", "symbol", "gene_symbol", "gene_id")) {
    if (col %in% colnames(rd)) {
      vals <- rd[[col]]
      if (is.character(vals)) return(vals)
    }
  }
  rn <- rownames(spe)
  if (is.character(rn)) rn else NULL
}

compute_qc_metrics <- function(spe) {
  counts <- SummarizedExperiment::assay(spe, "counts")
  if (!inherits(counts, "dgCMatrix")) {
    counts <- Matrix::Matrix(counts, sparse = TRUE)
    SummarizedExperiment::assay(spe, "counts") <- counts
  }

  total_counts <- Matrix::colSums(counts)
  detected_genes <- Matrix::colSums(counts > 0)

  gene_symbols <- get_gene_symbols(spe)
  percent_mito <- NULL
  if (!is.null(gene_symbols)) {
    mito_idx <- startsWith(toupper(gene_symbols), "MT-")
    if (any(mito_idx)) {
      mito_counts <- Matrix::colSums(counts[mito_idx, , drop = FALSE])
      percent_mito <- (mito_counts / pmax(total_counts, 1)) * 100
    }
  }

  tissue_flag <- NULL
  cd <- SummarizedExperiment::colData(spe)
  if ("in_tissue" %in% colnames(cd)) tissue_flag <- cd$in_tissue
  else if ("tissue" %in% colnames(cd)) tissue_flag <- cd$tissue

  list(
    spe = spe,
    total_counts = total_counts,
    detected_genes = detected_genes,
    percent_mito = percent_mito,
    tissue_flag = tissue_flag
  )
}

plot_qc <- function(path, pre, post) {
  png(path, width = 1400, height = 900)
  par(mfrow = c(2, 3))

  hist(log10(pre$total_counts + 1), breaks = 50, main = "Pre log10(counts+1)", xlab = "log10(counts+1)")
  hist(pre$detected_genes, breaks = 50, main = "Pre detected genes", xlab = "detected genes")
  if (!is.null(pre$percent_mito)) {
    hist(pre$percent_mito, breaks = 50, main = "Pre % mito", xlab = "percent mito")
  } else {
    plot.new(); text(0.5, 0.5, "Pre % mito: NA")
  }

  hist(log10(post$total_counts + 1), breaks = 50, main = "Post log10(counts+1)", xlab = "log10(counts+1)")
  hist(post$detected_genes, breaks = 50, main = "Post detected genes", xlab = "detected genes")
  if (!is.null(post$percent_mito)) {
    hist(post$percent_mito, breaks = 50, main = "Post % mito", xlab = "percent mito")
  } else {
    plot.new(); text(0.5, 0.5, "Post % mito: NA")
  }
  dev.off()
}

summary_rows <- vector("list", length(input_files))

for (i in seq_along(input_files)) {
  path <- input_files[[i]]
  sid <- if (!is.null(sample_id) && sample_id != "") sample_id else get_sample_id(path)
  sample_start <- Sys.time()
  log_line(log_file, paste0("02_qc: processing ", sid))

  spe <- readRDS(path)
  qc <- compute_qc_metrics(spe)
  spe <- qc$spe

  keep <- qc$total_counts >= qc_params$min_counts & qc$detected_genes >= qc_params$min_genes
  if (!is.null(qc$percent_mito) && !is.null(qc_params$max_mito)) {
    keep <- keep & (qc$percent_mito <= qc_params$max_mito)
  }
  if (!is.null(qc$tissue_flag) && isTRUE(qc_params$keep_tissue_only)) {
    keep <- keep & (as.numeric(qc$tissue_flag) == 1)
  }

  spe_filtered <- spe[, keep]
  if (nrow(SpatialExperiment::imgData(spe_filtered)) > 0) {
    SpatialExperiment::imgData(spe_filtered) <- SpatialExperiment::imgData(spe_filtered)[0, ]
  }

  out_new <- file.path(out_processed, paste0(sid, ".02_qc.rds"))
  out_legacy <- file.path(out_processed, paste0(sid, "_qc.rds"))
  saveRDS(spe_filtered, out_new)
  if (!file.exists(out_legacy)) saveRDS(spe_filtered, out_legacy)

  qc_after <- compute_qc_metrics(spe_filtered)
  plot_qc(file.path(plot_dir, paste0(sid, "_qc.png")), qc, qc_after)

  pct_removed <- (1 - (ncol(spe_filtered) / ncol(spe))) * 100
  summary_rows[[i]] <- data.frame(
    sample_id = sid,
    spots_before = ncol(spe),
    spots_after = ncol(spe_filtered),
    pct_removed = pct_removed,
    median_counts_before = median(qc$total_counts),
    median_counts_after = median(qc_after$total_counts),
    median_genes_before = median(qc$detected_genes),
    median_genes_after = median(qc_after$detected_genes),
    median_mito_before = if (is.null(qc$percent_mito)) NA_real_ else median(qc$percent_mito),
    median_mito_after = if (is.null(qc_after$percent_mito)) NA_real_ else median(qc_after$percent_mito),
    stringsAsFactors = FALSE
  )

  runtime_sec <- as.numeric(difftime(Sys.time(), sample_start, units = "secs"))
  write_summary("02_qc", sid, ncol(spe), ncol(spe_filtered), runtime_sec, "success", out_outputs)

  rm(spe, spe_filtered)
  gc()
}

thresholds_path <- file.path(qc_dir, "qc_thresholds.csv")
write.csv(
  data.frame(
    min_counts = qc_params$min_counts,
    min_genes = qc_params$min_genes,
    max_mito = qc_params$max_mito,
    keep_tissue_only = qc_params$keep_tissue_only,
    stringsAsFactors = FALSE
  ),
  thresholds_path,
  row.names = FALSE
)

summary_path <- file.path(qc_dir, "qc_summary_postfilter.csv")
write.csv(do.call(rbind, summary_rows), summary_path, row.names = FALSE)

if (!file.exists(summary_path)) stop("Missing outputs/qc/qc_summary_postfilter.csv")
log_line(log_file, "02_qc complete")
message("02_qc complete")
