#!/usr/bin/env Rscript

# Collect step summary TSVs into wide and long CSVs.

suppressPackageStartupMessages({
  library(utils)
})

pipeline_root <- Sys.getenv("PIPELINE_ROOT", getwd())
source(file.path(pipeline_root, "scripts", "helpers", "cli_utils.R"))

args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
outdir <- if (!is.null(args$outdir) && args$outdir != "") args$outdir else getwd()

outputs_dir <- file.path(outdir, "outputs")
summary_dirs <- list.files(outputs_dir, pattern = "^0[0-9]_.*", full.names = TRUE)

summary_files <- unlist(lapply(summary_dirs, function(d) {
  list.files(file.path(d, "summary"), pattern = "\\.tsv$", full.names = TRUE)
}))

if (length(summary_files) == 0) {
  long_out_dir <- file.path(outputs_dir, "logs")
  dir.create(long_out_dir, recursive = TRUE, showWarnings = FALSE)
  long_path <- file.path(long_out_dir, "step_summaries_long.csv")
  wide_path <- file.path(long_out_dir, "step_summaries_wide.csv")
  write.csv(data.frame(), long_path, row.names = FALSE)
  write.csv(data.frame(), wide_path, row.names = FALSE)
  message("No step summaries found. Wrote empty summary files.")
  quit(save = "no", status = 0)
}

read_one <- function(path) {
  df <- read.delim(path, stringsAsFactors = FALSE)
  step <- basename(dirname(dirname(path)))
  df$step <- step
  df
}

long_df <- do.call(rbind, lapply(summary_files, read_one))

required <- c("sample_id", "spots_in", "spots_out", "runtime_sec", "status", "step")
missing_cols <- setdiff(required, colnames(long_df))
if (length(missing_cols) > 0) {
  stop("Missing required columns in summaries: ", paste(missing_cols, collapse = ", "))
}

long_out_dir <- file.path(outputs_dir, "logs")
dir.create(long_out_dir, recursive = TRUE, showWarnings = FALSE)

long_path <- file.path(long_out_dir, "step_summaries_long.csv")
write.csv(long_df, long_path, row.names = FALSE)

wide_by_step <- split(long_df, long_df$step)
wide_list <- lapply(names(wide_by_step), function(step) {
  df <- wide_by_step[[step]]
  df <- df[, c("sample_id", "spots_in", "spots_out", "runtime_sec", "status")]
  colnames(df)[-1] <- paste0(step, "_", colnames(df)[-1])
  df
})

wide_df <- Reduce(function(x, y) merge(x, y, by = "sample_id", all = TRUE), wide_list)
wide_path <- file.path(long_out_dir, "step_summaries_wide.csv")
write.csv(wide_df, wide_path, row.names = FALSE)

message("Wrote ", long_path, " and ", wide_path)
