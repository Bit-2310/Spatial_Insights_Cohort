#!/usr/bin/env Rscript

# 01_ingest: Data ingest and caching only

suppressPackageStartupMessages({
  library(spatialLIBD)
  library(SpatialExperiment)
  library(Matrix)
  library(dplyr)
  library(data.table)
})

pipeline_root <- Sys.getenv("PIPELINE_ROOT", getwd())
source(file.path(pipeline_root, "scripts", "helpers", "cli_utils.R"))

args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
sample_id <- args$sample_id
dataset <- if (!is.null(args$dataset) && args$dataset != "") args$dataset else "dlpfc"
input_spe_rds <- args$input_spe_rds
spe_source <- args$spe_source
out_processed <- if (!is.null(args$out_processed)) args$out_processed else file.path(getwd(), "data", "processed")
out_outputs <- if (!is.null(args$out_outputs)) args$out_outputs else file.path(getwd(), "outputs")
log_file <- args$log_file

start_time <- Sys.time()

raw_dir <- file.path(pipeline_root, "data", "raw")

# Directory setup
for (d in c(raw_dir, out_processed, file.path(out_outputs, "qc"))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

source_path <- NULL
if (!is.null(spe_source) && spe_source != "") {
  source_path <- spe_source
} else if (!is.null(input_spe_rds) && input_spe_rds != "") {
  source_path <- input_spe_rds
}

if (!is.null(source_path)) {
  if (!file.exists(source_path)) {
    alt_path <- file.path(pipeline_root, source_path)
    if (file.exists(alt_path)) {
      source_path <- alt_path
    } else {
      stop("input_spe_rds not found: ", source_path)
    }
  }
  log_line(log_file, paste0("01_ingest: reading input spe from ", source_path))
  spe_all <- readRDS(source_path)
} else if (dataset == "dlpfc") {
  log_line(log_file, "01_ingest: fetching DLPFC spatial dataset")
  spe_all <- spatialLIBD::fetch_data(type = "spe")
} else {
  stop("Unsupported dataset: ", dataset, ". Supported options: dlpfc or --input_spe_rds.")
}
if (!"sample_id" %in% colnames(colData(spe_all))) {
  stop("sample_id missing in fetched object")
}

sample_ids <- NULL
if ("sample_id" %in% colnames(colData(spe_all))) {
  sample_ids <- unique(colData(spe_all)$sample_id)
}
if (!is.null(sample_id) && sample_id != "") {
  if (!is.null(sample_ids) && sample_id %in% sample_ids) {
    sample_ids <- sample_id
  } else {
    # For single-sample inputs without matching sample_id, treat as one sample
    sample_ids <- sample_id
    SummarizedExperiment::colData(spe_all)$sample_id <- sample_id
  }
}

qc_summary <- list()

for (sid in sample_ids) {
  sample_start <- Sys.time()
  log_line(log_file, paste0("01_ingest: processing ", sid))

  out_file_new <- file.path(out_processed, paste0(sid, ".01_ingest.rds"))
  out_file_legacy <- file.path(out_processed, paste0(sid, ".rds"))

  if (file.exists(out_file_new)) {
    spe <- readRDS(out_file_new)
  } else {
    spe <- spe_all[, spe_all$sample_id == sid]

    if ("imgData" %in% slotNames(spe)) {
      imgData(spe) <- imgData(spe)[0, ]
    }

    saveRDS(spe, out_file_new)
    if (!file.exists(out_file_legacy)) {
      saveRDS(spe, out_file_legacy)
    }
  }

  counts <- counts(spe)
  n_spots <- ncol(counts)
  n_genes <- nrow(counts)

  cd <- colData(spe)
  tissue_col <- intersect(c("in_tissue", "tissue"), colnames(cd))

  pct_tissue <- NA_real_
  if (length(tissue_col) == 1) {
    vals <- cd[[tissue_col]]
    if (is.logical(vals)) {
      pct_tissue <- mean(vals, na.rm = TRUE) * 100
    } else {
      pct_tissue <- mean(as.numeric(vals) == 1, na.rm = TRUE) * 100
    }
  }

  qc_summary[[sid]] <- data.frame(
    sample_id = sid,
    number_of_spots = n_spots,
    number_of_genes = n_genes,
    percent_tissue_spots = pct_tissue
  )

  runtime_sec <- as.numeric(difftime(Sys.time(), sample_start, units = "secs"))
  write_summary("01_ingest", sid, n_spots, n_spots, runtime_sec, "success", out_outputs)
}

qc_df <- bind_rows(qc_summary)
qc_out <- file.path(out_outputs, "qc", "initial_qc_summary.csv")
fwrite(qc_df, qc_out)

missing <- setdiff(
  paste0(sample_ids, ".01_ingest.rds"),
  list.files(out_processed)
)
if (length(missing) > 0) {
  stop("Missing cached samples: ", paste(missing, collapse = ", "))
}
if (!file.exists(qc_out)) stop("QC summary file not created")

log_line(log_file, "01_ingest complete")
message("01_ingest complete: ingest + caching successful")
