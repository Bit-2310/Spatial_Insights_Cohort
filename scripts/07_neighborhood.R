#!/usr/bin/env Rscript

# 07_neighborhood: Neighborhood adjacency

suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(SummarizedExperiment)
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

nb_dir <- file.path(out_outputs, "neighborhood")

dir.create(out_processed, recursive = TRUE, showWarnings = FALSE)
dir.create(nb_dir, recursive = TRUE, showWarnings = FALSE)

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
if (length(domain_files) == 0) stop("No input *_domains.rds files found for 07_neighborhood.")

get_sample_id <- function(path) sub("(\\.05_domains|_domains)\\.rds$", "", basename(path))

for (i in seq_along(domain_files)) {
  domain_path <- domain_files[[i]]
  sid <- if (!is.null(sample_id) && sample_id != "") sample_id else get_sample_id(domain_path)
  sample_start <- Sys.time()

  log_line(log_file, paste0("07_neighborhood: processing ", sid))

  spe_dom <- readRDS(domain_path)
  coords <- SpatialExperiment::spatialCoords(spe_dom)
  if (is.null(coords)) stop("Spatial coordinates missing for sample ", sid)

  domains <- SummarizedExperiment::colData(spe_dom)$domain_main
  if (is.null(domains)) stop("domain_main missing for sample ", sid)

  n <- length(domains)
  nn <- RANN::nn2(coords, k = 7)
  idx <- nn$nn.idx[, -1, drop = FALSE]
  i_idx <- rep(seq_len(n), each = 6)
  j_idx <- as.vector(t(idx))

  dom_i <- as.character(domains[i_idx])
  dom_j <- as.character(domains[j_idx])

  pairs <- paste(pmin(dom_i, dom_j), pmax(dom_i, dom_j), sep = "|")
  pair_counts <- sort(table(pairs), decreasing = TRUE)

  dom_levels <- sort(unique(as.character(domains)))
  mat <- matrix(0, nrow = length(dom_levels), ncol = length(dom_levels),
                dimnames = list(dom_levels, dom_levels))
  for (p in names(pair_counts)) {
    parts <- strsplit(p, "\\|", fixed = FALSE)[[1]]
    a <- parts[1]
    b <- parts[2]
    count <- as.integer(pair_counts[[p]])
    mat[a, b] <- mat[a, b] + count
    mat[b, a] <- mat[b, a] + count
  }

  nb_mat_path <- file.path(nb_dir, paste0("neighborhood_matrix_", sid, ".csv"))
  write.csv(mat, nb_mat_path, row.names = TRUE)

  top_pairs <- head(pair_counts, 25)
  top_pairs_df <- data.frame(
    domain_pair = names(top_pairs),
    score = as.integer(top_pairs),
    stringsAsFactors = FALSE
  )
  nb_pairs_path <- file.path(nb_dir, paste0("neighborhood_top_pairs_", sid, ".csv"))
  write.csv(top_pairs_df, nb_pairs_path, row.names = FALSE)

  out_cache <- file.path(out_processed, paste0(sid, ".07_neighborhood.rds"))
  saveRDS(list(matrix = mat, top_pairs = top_pairs_df), out_cache)

  runtime_sec <- as.numeric(difftime(Sys.time(), sample_start, units = "secs"))
  write_summary("07_neighborhood", sid, ncol(spe_dom), ncol(spe_dom), runtime_sec, "success", out_outputs)

  rm(spe_dom)
  gc()
}

nb_pair_files <- list.files(nb_dir, pattern = "neighborhood_top_pairs_.*\\.csv$", full.names = TRUE)
pair_rows <- lapply(nb_pair_files, function(path) {
  df <- read.csv(path, stringsAsFactors = FALSE)
  df$sample_id <- sub("^neighborhood_top_pairs_(.*)\\.csv$", "\\1", basename(path))
  df
})
pair_rows <- do.call(rbind, pair_rows)

if (is.null(pair_rows) || nrow(pair_rows) == 0) {
  stop("No neighborhood pair rows found to summarize.")
}

n_samples_top25 <- tapply(pair_rows$sample_id, pair_rows$domain_pair, function(x) length(unique(x)))
mean_score <- tapply(pair_rows$score, pair_rows$domain_pair, mean)
median_score <- tapply(pair_rows$score, pair_rows$domain_pair, median)

pairs <- names(n_samples_top25)
nb_summary <- data.frame(
  domain_pair = pairs,
  n_samples_top25 = as.integer(n_samples_top25[pairs]),
  mean_score = as.numeric(mean_score[pairs]),
  median_score = as.numeric(median_score[pairs]),
  stringsAsFactors = FALSE
)

nb_summary <- nb_summary[order(-nb_summary$n_samples_top25, -nb_summary$mean_score), , drop = FALSE]
write.csv(nb_summary, file.path(nb_dir, "neighborhood_cohort_summary.csv"), row.names = FALSE)

consensus_pairs <- head(nb_summary, 50)
write.csv(consensus_pairs, file.path(nb_dir, "neighborhood_consensus_pairs.csv"), row.names = FALSE)

expected_ids <- vapply(domain_files, get_sample_id, character(1))
missing_nb_mat <- expected_ids[
  !file.exists(file.path(nb_dir, paste0("neighborhood_matrix_", expected_ids, ".csv")))
]
missing_nb_pairs <- expected_ids[
  !file.exists(file.path(nb_dir, paste0("neighborhood_top_pairs_", expected_ids, ".csv")))
]

if (length(missing_nb_mat) > 0) stop("Missing neighborhood matrices for samples: ", paste(missing_nb_mat, collapse = ", "))
if (length(missing_nb_pairs) > 0) stop("Missing neighborhood top pairs for samples: ", paste(missing_nb_pairs, collapse = ", "))
if (!file.exists(file.path(nb_dir, "neighborhood_cohort_summary.csv"))) stop("Missing outputs/neighborhood/neighborhood_cohort_summary.csv")
if (!file.exists(file.path(nb_dir, "neighborhood_consensus_pairs.csv"))) stop("Missing outputs/neighborhood/neighborhood_consensus_pairs.csv")

log_line(log_file, "07_neighborhood complete")
message("07_neighborhood complete")
