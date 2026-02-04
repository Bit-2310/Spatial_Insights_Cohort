#!/usr/bin/env Rscript

# 06_svg: Spatially variable genes

suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(SummarizedExperiment)
  library(SingleCellExperiment)
  library(Matrix)
  library(scran)
  library(RANN)
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

svg_dir <- file.path(out_outputs, "svg")
svg_plot_dir <- file.path(svg_dir, "plots")

dir.create(out_processed, recursive = TRUE, showWarnings = FALSE)
dir.create(svg_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(svg_plot_dir, recursive = TRUE, showWarnings = FALSE)

if (!is.null(in_rds) && in_rds != "") {
  domain_files <- in_rds
} else if (!is.null(sample_id) && sample_id != "") {
  domain_files <- file.path(out_processed, paste0(sample_id, ".05_domains.rds"))
  if (!file.exists(domain_files)) {
    legacy <- file.path(out_processed, paste0(sample_id, "_domains.rds"))
    if (file.exists(legacy)) domain_files <- legacy
  }
} else {
  domain_files <- list.files(out_processed, pattern = "\\.05_domains\\.rds$", full.names = TRUE)
}

domain_files <- domain_files[file.exists(domain_files)]
if (length(domain_files) == 0) {
  domain_files <- list.files(out_processed, pattern = "_domains\\.rds$", full.names = TRUE)
}
if (length(domain_files) == 0) stop("No input *_domains.rds files found for 06_svg.")

get_sample_id <- function(path) sub("(\\.05_domains|_domains)\\.rds$", "", basename(path))

knn_graph <- function(coords, k = 6) {
  nn <- RANN::nn2(coords, k = k + 1)
  idx <- nn$nn.idx[, -1, drop = FALSE]
  n <- nrow(coords)
  i <- rep(seq_len(n), each = k)
  j <- as.vector(t(idx))
  Matrix::sparseMatrix(i = i, j = j, x = 1, dims = c(n, n))
}

moran_scores <- function(expr, W) {
  n <- ncol(expr)
  S0 <- sum(W)
  X <- t(expr)
  X <- scale(X, center = TRUE, scale = FALSE)
  WX <- W %*% X
  num <- Matrix::colSums(X * WX)
  denom <- Matrix::colSums(X * X)
  (n / S0) * (num / denom)
}

cohort_svg_rows <- list()

for (i in seq_along(domain_files)) {
  domain_path <- domain_files[[i]]
  sid <- if (!is.null(sample_id) && sample_id != "") sample_id else get_sample_id(domain_path)
  sample_start <- Sys.time()
  norm_path <- file.path(out_processed, paste0(sid, ".03_normalize.rds"))
  if (!file.exists(norm_path)) {
    legacy_norm <- file.path(out_processed, paste0(sid, "_norm.rds"))
    if (file.exists(legacy_norm)) norm_path <- legacy_norm
  }

  log_line(log_file, paste0("06_svg: processing ", sid))

  spe_dom <- readRDS(domain_path)
  spe_expr <- spe_dom
  if (!"logcounts" %in% SummarizedExperiment::assayNames(spe_expr)) {
    if (!file.exists(norm_path)) stop("Missing ", norm_path, " for sample ", sid)
    spe_expr <- readRDS(norm_path)
  }

  if (!"logcounts" %in% SummarizedExperiment::assayNames(spe_expr)) stop("logcounts assay missing for sample ", sid)

  coords <- SpatialExperiment::spatialCoords(spe_dom)
  if (is.null(coords)) stop("Spatial coordinates missing for sample ", sid)

  W <- knn_graph(coords, k = 6)
  W <- (W + Matrix::t(W))
  W@x[W@x > 0] <- 1

  var_fit <- scran::modelGeneVar(spe_expr)
  var_fit <- var_fit[order(var_fit$bio, decreasing = TRUE), , drop = FALSE]
  top_n <- min(500, nrow(var_fit))
  top_genes <- rownames(var_fit)[seq_len(top_n)]

  expr <- SummarizedExperiment::assay(spe_expr, "logcounts")[top_genes, , drop = FALSE]
  scores <- moran_scores(expr, W)

  svg_df <- data.frame(
    gene = top_genes,
    score = as.numeric(scores),
    stringsAsFactors = FALSE
  )
  svg_df <- svg_df[order(svg_df$score, decreasing = TRUE), , drop = FALSE]
  svg_df$rank <- seq_len(nrow(svg_df))

  svg_ranked_path <- file.path(svg_dir, paste0("svg_ranked_", sid, ".csv"))
  write.csv(svg_df[, c("gene", "score", "rank")], svg_ranked_path, row.names = FALSE)

  top12 <- head(svg_df, 12)
  png(file.path(svg_plot_dir, paste0(sid, "_top_svg.png")), width = 1200, height = 800)
  barplot(
    height = top12$score,
    names.arg = top12$gene,
    las = 2,
    cex.names = 0.7,
    ylab = "Moran's I",
    main = paste0(sid, " top 12 SVGs")
  )
  dev.off()

  cohort_svg_rows[[sid]] <- head(svg_df, 50)

  out_cache <- file.path(out_processed, paste0(sid, ".06_svg.rds"))
  saveRDS(list(svg_ranked = svg_df), out_cache)

  runtime_sec <- as.numeric(difftime(Sys.time(), sample_start, units = "secs"))
  write_summary("06_svg", sid, ncol(spe_dom), ncol(spe_dom), runtime_sec, "success", out_outputs)

  rm(spe_dom, spe_expr, expr, W)
  gc()
}

all_top50 <- do.call(rbind, lapply(names(cohort_svg_rows), function(sid) {
  df <- cohort_svg_rows[[sid]]
  df$sample_id <- sid
  df
}))

if (is.null(all_top50) || nrow(all_top50) == 0) {
  stop("No SVG rows found to summarize.")
}

n_samples_in_top50 <- tapply(all_top50$sample_id, all_top50$gene, function(x) length(unique(x)))
mean_rank <- tapply(all_top50$rank, all_top50$gene, mean)
median_rank <- tapply(all_top50$rank, all_top50$gene, median)
mean_score <- tapply(all_top50$score, all_top50$gene, mean)

genes <- names(n_samples_in_top50)
svg_summary <- data.frame(
  gene = genes,
  n_samples_in_top50 = as.integer(n_samples_in_top50[genes]),
  mean_rank = as.numeric(mean_rank[genes]),
  median_rank = as.numeric(median_rank[genes]),
  mean_score = as.numeric(mean_score[genes]),
  stringsAsFactors = FALSE
)

svg_summary <- svg_summary[order(-svg_summary$n_samples_in_top50, svg_summary$mean_rank), , drop = FALSE]
write.csv(svg_summary, file.path(svg_dir, "svg_cohort_summary.csv"), row.names = FALSE)

consensus_top50 <- head(svg_summary, 50)
write.csv(consensus_top50, file.path(svg_dir, "svg_consensus_top50.csv"), row.names = FALSE)

expected_ids <- vapply(domain_files, get_sample_id, character(1))
missing_svg <- expected_ids[
  !file.exists(file.path(svg_dir, paste0("svg_ranked_", expected_ids, ".csv")))
]

if (length(missing_svg) > 0) stop("Missing SVG ranked CSVs for samples: ", paste(missing_svg, collapse = ", "))
if (!file.exists(file.path(svg_dir, "svg_cohort_summary.csv"))) stop("Missing outputs/svg/svg_cohort_summary.csv")
if (!file.exists(file.path(svg_dir, "svg_consensus_top50.csv"))) stop("Missing outputs/svg/svg_consensus_top50.csv")

log_line(log_file, "06_svg complete")
message("06_svg complete")
