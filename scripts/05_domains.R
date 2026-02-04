#!/usr/bin/env Rscript

# 05_domains: Spatial domains and stability

suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(SummarizedExperiment)
  library(SingleCellExperiment)
  library(mclust)
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

domains_dir <- file.path(out_outputs, "domains")
plot_dir <- file.path(domains_dir, "plots")

dir.create(out_processed, recursive = TRUE, showWarnings = FALSE)
dir.create(domains_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

input_files <- list.files(out_processed, pattern = "\\.03_normalize\\.rds$", full.names = TRUE)
if (!is.null(sample_id) && sample_id != "") {
  input_files <- file.path(out_processed, paste0(sample_id, ".03_normalize.rds"))
  if (!file.exists(input_files)) {
    legacy <- file.path(out_processed, paste0(sample_id, "_norm.rds"))
    if (file.exists(legacy)) input_files <- legacy
  }
}
if (!is.null(in_rds) && in_rds != "") {
  input_files <- in_rds
}
input_files <- input_files[file.exists(input_files)]
if (length(input_files) == 0) {
  legacy <- list.files(out_processed, pattern = "_norm\\.rds$", full.names = TRUE)
  input_files <- legacy
}
if (length(input_files) == 0) stop("No input *_norm.rds files found for 05_domains.")

get_sample_id <- function(path) sub("(\\.03_normalize|_norm)\\.rds$", "", basename(path))

settings <- c(6, 7, 8)
main_setting <- settings[2]
use_bayesspace <- requireNamespace("BayesSpace", quietly = TRUE)
method_name <- if (use_bayesspace) "BayesSpace" else "kmeans"

safe_drop_images <- function(spe) {
  img <- try(SpatialExperiment::imgData(spe), silent = TRUE)
  if (!inherits(img, "try-error") && !is.null(img) && nrow(img) > 0) {
    SpatialExperiment::imgData(spe) <- img[0, ]
  }
  spe
}

ensure_col_row <- function(spe_pre, spe_raw) {
  cd <- SummarizedExperiment::colData(spe_pre)
  if (all(c("col", "row") %in% colnames(cd))) return(spe_pre)
  if (all(c("array_col", "array_row") %in% colnames(cd))) {
    SummarizedExperiment::colData(spe_pre)$col <- as.integer(cd$array_col)
    SummarizedExperiment::colData(spe_pre)$row <- as.integer(cd$array_row)
    return(spe_pre)
  }
  coords <- SpatialExperiment::spatialCoords(spe_raw)
  if (is.null(coords) || ncol(coords) < 2) stop("Spatial coordinates missing for BayesSpace.")
  SummarizedExperiment::colData(spe_pre)$col <- as.integer(round(coords[, 1]))
  SummarizedExperiment::colData(spe_pre)$row <- as.integer(round(coords[, 2]))
  spe_pre
}

plot_domains <- function(path, coords, labels, title) {
  png(path, width = 1200, height = 900)
  cols <- as.integer(factor(labels))
  plot(
    coords[, 1], coords[, 2],
    col = cols, pch = 16, cex = 0.6,
    xlab = "x", ylab = "y",
    main = title
  )
  dev.off()
}

run_domains_bayesspace <- function(spe_pre, q) {
  set.seed(100 + q)
  spe_q <- BayesSpace::spatialCluster(
    spe_pre,
    q = q,
    platform = "Visium",
    init.method = "mclust",
    nrep = 50,
    burn.in = 20
  )
  as.integer(SummarizedExperiment::colData(spe_q)$spatial.cluster)
}

run_domains_kmeans <- function(spe, q) {
  pcs <- SingleCellExperiment::reducedDim(spe, "PCA")
  if (is.null(pcs)) stop("PCA not found in normalized object.")
  coords <- SpatialExperiment::spatialCoords(spe)
  if (is.null(coords)) stop("Spatial coordinates missing from object.")
  x <- cbind(pcs[, 1:20, drop = FALSE], scale(coords))
  set.seed(100 + q)
  kmeans(x, centers = q, nstart = 10)$cluster
}

cohort_rows <- vector("list", length(input_files))

for (i in seq_along(input_files)) {
  path <- input_files[[i]]
  sid <- if (!is.null(sample_id) && sample_id != "") sample_id else get_sample_id(path)
  sample_start <- Sys.time()
  log_line(log_file, paste0("05_domains: processing ", sid, " using ", method_name))

  spe <- readRDS(path)
  out_path_new <- file.path(out_processed, paste0(sid, ".05_domains.rds"))
  out_path_legacy <- file.path(out_processed, paste0(sid, "_domains.rds"))

  spe_pre <- NULL
  if (use_bayesspace) {
    spe_pre <- BayesSpace::spatialPreprocess(
      spe,
      platform = "Visium",
      n.PCs = 20,
      log.normalize = TRUE
    )
    spe_pre <- ensure_col_row(spe_pre, spe)
  }

  labels_list <- vector("list", length(settings))
  for (j in seq_along(settings)) {
    q <- settings[[j]]
    if (use_bayesspace) {
      labels_list[[j]] <- run_domains_bayesspace(spe_pre, q)
    } else {
      labels_list[[j]] <- run_domains_kmeans(spe, q)
    }
    if (length(labels_list[[j]]) != ncol(spe)) stop("Domain labels length mismatch for q=", q, " in sample ", sid)
  }

  pairs <- combn(seq_along(settings), 2)
  pair_agree <- apply(pairs, 2, function(idx) {
    mclust::adjustedRandIndex(labels_list[[idx[1]]], labels_list[[idx[2]]])
  })
  stability_score <- mean(pair_agree)

  SummarizedExperiment::colData(spe)$domain_q6 <- labels_list[[1]]
  SummarizedExperiment::colData(spe)$domain_q7 <- labels_list[[2]]
  SummarizedExperiment::colData(spe)$domain_q8 <- labels_list[[3]]
  SummarizedExperiment::colData(spe)$domain_main <- labels_list[[2]]

  spe <- safe_drop_images(spe)

  saveRDS(spe, out_path_new)
  if (!file.exists(out_path_legacy)) saveRDS(spe, out_path_legacy)

  coords <- SpatialExperiment::spatialCoords(spe)
  plot_domains(
    path = file.path(plot_dir, paste0(sid, "_domains.png")),
    coords = coords,
    labels = SummarizedExperiment::colData(spe)$domain_main,
    title = paste0(sid, " domains (q=", main_setting, ", ", method_name, ")")
  )

  counts <- table(SummarizedExperiment::colData(spe)$domain_main)
  write.csv(
    data.frame(domain = names(counts), spot_count = as.integer(counts)),
    file.path(domains_dir, paste0(sid, "_domain_counts.csv")),
    row.names = FALSE
  )

  stability_path <- file.path(domains_dir, paste0(sid, "_stability.csv"))
  stability_df <- data.frame(
    setting_1 = settings[pairs[1, ]],
    setting_2 = settings[pairs[2, ]],
    agreement = pair_agree,
    stability_score = stability_score,
    stringsAsFactors = FALSE
  )
  write.csv(stability_df, stability_path, row.names = FALSE)

  top_counts <- sort(counts, decreasing = TRUE)
  top_sizes <- paste0(
    names(top_counts)[seq_len(min(3, length(top_counts)))],
    ":",
    as.integer(top_counts[seq_len(min(3, length(top_counts)))]),
    collapse = ";"
  )

  cohort_rows[[i]] <- data.frame(
    sample_id = sid,
    n_domains = length(counts),
    stability_score = stability_score,
    top_domain_sizes = top_sizes,
    stringsAsFactors = FALSE
  )

  runtime_sec <- as.numeric(difftime(Sys.time(), sample_start, units = "secs"))
  write_summary("05_domains", sid, ncol(spe), ncol(spe), runtime_sec, "success", out_outputs)

  rm(spe, spe_pre, labels_list)
  gc()
}

cohort_summary_path <- file.path(domains_dir, "domains_cohort_summary.csv")
write.csv(do.call(rbind, cohort_rows), cohort_summary_path, row.names = FALSE)

png(file.path(plot_dir, "domains_stability_by_sample.png"), width = 1200, height = 800)
barplot(
  vapply(cohort_rows, function(x) x$stability_score, numeric(1)),
  names.arg = vapply(cohort_rows, function(x) x$sample_id, character(1)),
  las = 2,
  ylab = "Stability score (mean ARI)",
  main = "Domain stability by sample"
)

dev.off()

expected_ids <- vapply(input_files, get_sample_id, character(1))
missing_domains <- expected_ids[!file.exists(file.path(out_processed, paste0(expected_ids, ".05_domains.rds")))]
if (length(missing_domains) > 0) stop("Missing domain .rds files for samples: ", paste(missing_domains, collapse = ", "))
if (!file.exists(cohort_summary_path)) stop("Missing outputs/domains/domains_cohort_summary.csv")
if (!file.exists(file.path(plot_dir, "domains_stability_by_sample.png"))) stop("Missing domains_stability_by_sample.png")

log_line(log_file, "05_domains complete")
message("05_domains complete")
