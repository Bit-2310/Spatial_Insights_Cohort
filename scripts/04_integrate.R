#!/usr/bin/env Rscript

# 04_integrate: Cohort integration

suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(SummarizedExperiment)
  library(SingleCellExperiment)
  library(scater)
  library(scran)
  library(batchelor)
})

pipeline_root <- Sys.getenv("PIPELINE_ROOT", getwd())
source(file.path(pipeline_root, "scripts", "helpers", "cli_utils.R"))

args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
in_files <- args$in_files
out_processed <- if (!is.null(args$out_processed)) args$out_processed else file.path(getwd(), "data", "processed")
out_outputs <- if (!is.null(args$out_outputs)) args$out_outputs else file.path(getwd(), "outputs")
log_file <- args$log_file

start_time <- Sys.time()

out_dir <- file.path(out_outputs, "integration")
plot_dir <- file.path(out_dir, "plots")

dir.create(out_processed, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

input_files <- character(0)
if (!is.null(in_files) && in_files != "") {
  input_files <- unlist(strsplit(in_files, ",", fixed = TRUE))
  input_files <- trimws(input_files)
  input_files <- input_files[nzchar(input_files)]
  input_files <- input_files[file.exists(input_files)]
  if (length(input_files) == 0) {
    stop("No valid files found in --in_files.")
  }
}
if (length(input_files) == 0) {
  input_files <- list.files(out_processed, pattern = "\\.03_normalize\\.rds$", full.names = TRUE)
}
if (length(input_files) == 0) {
  input_files <- list.files(out_processed, pattern = "_norm\\.rds$", full.names = TRUE)
}
if (length(input_files) == 0) stop("No input *_norm.rds files found for 04_integrate.")

get_sample_id <- function(path) sub("(\\.03_normalize|_norm)\\.rds$", "", basename(path))

sce_list <- vector("list", length(input_files))
names(sce_list) <- vapply(input_files, get_sample_id, character(1))
for (i in seq_along(input_files)) sce_list[[i]] <- readRDS(input_files[[i]])

hvg_list <- lapply(sce_list, function(sce) {
  var_fit <- scran::modelGeneVar(sce)
  var_fit <- var_fit[order(var_fit$bio, decreasing = TRUE), , drop = FALSE]
  top_n <- min(2000, nrow(var_fit))
  rownames(var_fit)[seq_len(top_n)]
})

hvg_union <- unique(unlist(hvg_list))
if (length(hvg_union) > 2000) {
  bio_mat <- vapply(sce_list, function(sce) {
    vf <- scran::modelGeneVar(sce)
    vf$bio
  }, numeric(nrow(sce_list[[1]])))
  rownames(bio_mat) <- rownames(sce_list[[1]])
  mean_bio <- rowMeans(bio_mat[hvg_union, , drop = FALSE])
  hvg_union <- names(sort(mean_bio, decreasing = TRUE))[seq_len(2000)]
}

pcs_used <- 1:20
method_name <- "fastMNN"

integrated <- batchelor::fastMNN(
  sce_list,
  subset.row = hvg_union,
  d = max(pcs_used)
)

has_umap <- FALSE
if (requireNamespace("uwot", quietly = TRUE)) {
  integrated <- scater::runUMAP(integrated, dimred = "corrected", name = "UMAP")
  has_umap <- TRUE
}

out_path_new <- file.path(out_processed, "cohort.04_integrate.rds")
out_path_legacy <- file.path(out_processed, "dlpfc_cohort_integrated.rds")
saveRDS(integrated, out_path_new)
if (!file.exists(out_path_legacy)) saveRDS(integrated, out_path_legacy)

summary_path <- file.path(out_dir, "cohort_integration_summary.csv")
write.csv(
  data.frame(
    n_samples = length(sce_list),
    total_spots = ncol(integrated),
    hvgs_used = length(hvg_union),
    pcs_used = paste(pcs_used, collapse = ","),
    method_name = method_name,
    stringsAsFactors = FALSE
  ),
  summary_path,
  row.names = FALSE
)

pc <- SingleCellExperiment::reducedDim(integrated, "corrected")
sample_id <- SummarizedExperiment::colData(integrated)$sample_id
if (is.null(sample_id) && "batch" %in% colnames(SummarizedExperiment::colData(integrated))) {
  sample_id <- SummarizedExperiment::colData(integrated)$batch
  SummarizedExperiment::colData(integrated)$sample_id <- sample_id
}

png(file.path(plot_dir, "cohort_pca_by_sample.png"), width = 1200, height = 900)
plot(
  pc[, 1], pc[, 2],
  col = as.integer(factor(sample_id)),
  pch = 16, cex = 0.6,
  xlab = "PC1 (corrected)",
  ylab = "PC2 (corrected)",
  main = "Cohort corrected PCs by sample"
)
legend(
  "topright",
  legend = levels(factor(sample_id)),
  col = seq_along(levels(factor(sample_id))),
  pch = 16,
  cex = 0.7
)

dev.off()

if (has_umap) {
  umap <- SingleCellExperiment::reducedDim(integrated, "UMAP")
  png(file.path(plot_dir, "cohort_umap_by_sample.png"), width = 1200, height = 900)
  plot(
    umap[, 1], umap[, 2],
    col = as.integer(factor(sample_id)),
    pch = 16, cex = 0.6,
    xlab = "UMAP1",
    ylab = "UMAP2",
    main = "Cohort UMAP by sample"
  )
  legend(
    "topright",
    legend = levels(factor(sample_id)),
    col = seq_along(levels(factor(sample_id))),
    pch = 16,
    cex = 0.7
  )
  dev.off()
}

runtime_sec <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
write_summary("04_integrate", "cohort", ncol(integrated), ncol(integrated), runtime_sec, "success", out_outputs)

if (!file.exists(out_path_new)) stop("Missing data/processed/cohort.04_integrate.rds")
if (!file.exists(summary_path)) stop("Missing outputs/integration/cohort_integration_summary.csv")
if (!file.exists(file.path(plot_dir, "cohort_pca_by_sample.png"))) stop("Missing cohort_pca_by_sample.png")
if (has_umap && !file.exists(file.path(plot_dir, "cohort_umap_by_sample.png"))) stop("Missing cohort_umap_by_sample.png")

log_line(log_file, "04_integrate complete")
message("04_integrate complete")
