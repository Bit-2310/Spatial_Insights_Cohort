#!/usr/bin/env Rscript

# 03_normalize: Normalization and per-sample features

suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(SummarizedExperiment)
  library(SingleCellExperiment)
  library(scater)
  library(scran)
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

out_dir <- file.path(out_outputs, "integration")
plot_dir <- file.path(out_dir, "plots")

dir.create(out_processed, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

input_files <- list.files(out_processed, pattern = "\\.02_qc\\.rds$", full.names = TRUE)
if (!is.null(sample_id) && sample_id != "") {
  input_files <- file.path(out_processed, paste0(sample_id, ".02_qc.rds"))
  if (!file.exists(input_files)) {
    legacy <- file.path(out_processed, paste0(sample_id, "_qc.rds"))
    if (file.exists(legacy)) input_files <- legacy
  }
}
if (!is.null(in_rds) && in_rds != "") {
  input_files <- in_rds
}
input_files <- input_files[file.exists(input_files)]
if (length(input_files) == 0) {
  legacy <- list.files(out_processed, pattern = "_qc\\.rds$", full.names = TRUE)
  input_files <- legacy
}
if (length(input_files) == 0) stop("No input *_qc.rds files found for 03_normalize.")

get_sample_id <- function(path) sub("(\\.02_qc|_qc)\\.rds$", "", basename(path))

summary_rows <- vector("list", length(input_files))

for (i in seq_along(input_files)) {
  path <- input_files[[i]]
  sid <- if (!is.null(sample_id) && sample_id != "") sample_id else get_sample_id(path)
  sample_start <- Sys.time()
  log_line(log_file, paste0("03_normalize: processing ", sid))

  spe <- readRDS(path)
  spots_in <- ncol(spe)

  spe <- scater::logNormCounts(spe)
  var_fit <- scran::modelGeneVar(spe)
  var_fit <- var_fit[order(var_fit$bio, decreasing = TRUE), , drop = FALSE]
  top_n <- min(2000, nrow(var_fit))
  top_genes <- rownames(var_fit)[seq_len(top_n)]

  spe <- scater::runPCA(spe, subset_row = top_genes, ncomponents = 30, name = "PCA")

  if (nrow(SpatialExperiment::imgData(spe)) > 0) {
    SpatialExperiment::imgData(spe) <- SpatialExperiment::imgData(spe)[0, ]
  }

  out_new <- file.path(out_processed, paste0(sid, ".03_normalize.rds"))
  out_legacy <- file.path(out_processed, paste0(sid, "_norm.rds"))
  saveRDS(spe, out_new)
  if (!file.exists(out_legacy)) saveRDS(spe, out_legacy)

  pca <- SingleCellExperiment::reducedDim(spe, "PCA")
  percent_var <- attr(pca, "percentVar")
  if (is.null(percent_var)) percent_var <- rep(NA_real_, ncol(pca))

  png(file.path(plot_dir, paste0(sid, "_pca_var.png")), width = 1200, height = 800)
  plot(
    seq_along(percent_var),
    percent_var,
    type = "b",
    xlab = "PC",
    ylab = "Percent variance explained",
    main = paste0(sid, " PCA variance")
  )
  dev.off()

  summary_rows[[i]] <- data.frame(
    sample_id = sid,
    spots = ncol(spe),
    genes = nrow(spe),
    n_pcs = ncol(pca),
    top_variable_genes_count = length(top_genes),
    stringsAsFactors = FALSE
  )

  runtime_sec <- as.numeric(difftime(Sys.time(), sample_start, units = "secs"))
  write_summary("03_normalize", sid, spots_in, ncol(spe), runtime_sec, "success", out_outputs)

  rm(spe, var_fit, top_genes, pca)
  gc()
}

summary_path <- file.path(out_dir, "per_sample_norm_summary.csv")
write.csv(do.call(rbind, summary_rows), summary_path, row.names = FALSE)

if (!file.exists(summary_path)) stop("Missing outputs/integration/per_sample_norm_summary.csv")
log_line(log_file, "03_normalize complete")
message("03_normalize complete")
