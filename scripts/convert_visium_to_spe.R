#!/usr/bin/env Rscript

# Convert raw 10x Visium folders into SpatialExperiment RDS files.
# This does not alter the pipeline; it prepares inputs for --input_spe_rds or spe_source.

suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(parallel)
  library(DropletUtils)
  library(SummarizedExperiment)
  library(SingleCellExperiment)
})

root_dir <- file.path("data", "raw", "breast_cancer_visium")
if (!dir.exists(root_dir)) stop("Missing data directory: ", root_dir)

get_cores <- function() {
  opt <- getOption("N_CORES")
  if (!is.null(opt) && !is.na(opt)) return(as.integer(opt))
  cores <- parallel::detectCores()
  if (is.na(cores)) cores <- 1L
  min(8L, cores)
}

sample_dirs <- list.dirs(root_dir, recursive = FALSE, full.names = TRUE)
sample_ids <- basename(sample_dirs)
sample_ids <- sample_ids[grepl("^BC_[0-9]{2}$", sample_ids)]
sample_dirs <- file.path(root_dir, sample_ids)

if (length(sample_ids) == 0) stop("No BC_XX directories found under ", root_dir)

convert_one <- function(sample_id, sample_dir) {
  filt_dir <- file.path(sample_dir, "filtered_feature_bc_matrix")
  spatial_dir <- file.path(sample_dir, "spatial")
  out_path <- file.path(sample_dir, paste0(sample_id, ".spe.rds"))

  if (file.exists(out_path)) {
    return(data.frame(sample_id = sample_id, spe_source = out_path,
                      status = "skipped", message = "RDS exists", stringsAsFactors = FALSE))
  }

  if (!dir.exists(filt_dir) || !dir.exists(spatial_dir)) {
    return(data.frame(sample_id = sample_id, spe_source = out_path,
                      status = "skipped", message = "missing filtered_feature_bc_matrix or spatial",
                      stringsAsFactors = FALSE))
  }

  spe <- tryCatch({
    SpatialExperiment::read10xVisium(sample_dir)
  }, error = function(e) e)

  if (inherits(spe, "error")) {
    message("read10xVisium failed for ", sample_id, ": ", conditionMessage(spe))
    positions_path <- file.path(spatial_dir, "tissue_positions_list.csv")
    if (!file.exists(positions_path)) {
      positions_path <- file.path(spatial_dir, "tissue_positions.csv")
    }
    if (!file.exists(positions_path)) {
      return(data.frame(sample_id = sample_id, spe_source = out_path,
                        status = "failed", message = "missing tissue_positions_list.csv",
                        stringsAsFactors = FALSE))
    }

    pos <- read.csv(positions_path, header = FALSE, stringsAsFactors = FALSE)
    colnames(pos) <- c(
      "barcode", "in_tissue", "array_row", "array_col",
      "pxl_row_in_fullres", "pxl_col_in_fullres"
    )
    rownames(pos) <- pos$barcode

    sce <- DropletUtils::read10xCounts(filt_dir)
    barcodes <- colnames(sce)
    if (length(barcodes) == 0) {
      bc_path <- file.path(filt_dir, "barcodes.tsv.gz")
      if (!file.exists(bc_path)) bc_path <- file.path(filt_dir, "barcodes.tsv")
      if (!file.exists(bc_path)) {
        return(data.frame(sample_id = sample_id, spe_source = out_path,
                          status = "failed", message = "missing barcodes.tsv",
                          stringsAsFactors = FALSE))
      }
      barcodes <- readLines(gzfile(bc_path))
      colnames(sce) <- barcodes
    }
    pos <- pos[barcodes, , drop = FALSE]
    if (nrow(pos) == 0) {
      return(data.frame(sample_id = sample_id, spe_source = out_path,
                        status = "failed", message = "no matching barcodes",
                        stringsAsFactors = FALSE))
    }

    coords <- cbind(pos$pxl_col_in_fullres, pos$pxl_row_in_fullres)
    colnames(coords) <- c("x", "y")

    spe <- SpatialExperiment::SpatialExperiment(
      assays = SummarizedExperiment::assays(sce),
      rowData = SummarizedExperiment::rowData(sce),
      colData = cbind(SummarizedExperiment::colData(sce), pos[barcodes, , drop = FALSE]),
      spatialCoords = coords
    )
  }

  if (!"sample_id" %in% colnames(SummarizedExperiment::colData(spe))) {
    SummarizedExperiment::colData(spe)$sample_id <- sample_id
  }

  saveRDS(spe, out_path)
  data.frame(sample_id = sample_id, spe_source = out_path,
             status = "success", message = "converted", stringsAsFactors = FALSE)
}

cores <- get_cores()
message("Using cores: ", cores)

is_unix <- .Platform$OS.type == "unix"
results <- if (is_unix) {
  parallel::mclapply(seq_along(sample_ids), function(i) {
    convert_one(sample_ids[[i]], sample_dirs[[i]])
  }, mc.cores = cores)
} else {
  lapply(seq_along(sample_ids), function(i) {
    convert_one(sample_ids[[i]], sample_dirs[[i]])
  })
}

summary_df <- do.call(rbind, results)
summary_path <- file.path(root_dir, "conversion_summary.csv")
write.table(summary_df, summary_path, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)

manifest_path <- file.path(root_dir, "breast_cancer_manifest.csv")
manifest_df <- summary_df[
  summary_df$status == "success" | (summary_df$status == "skipped" & summary_df$message == "RDS exists"),
  c("sample_id", "spe_source")
]
write.table(manifest_df, manifest_path, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)

message("Summary written to: ", summary_path)
message("Manifest written to: ", manifest_path)
message("Done.")
