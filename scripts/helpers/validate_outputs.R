#!/usr/bin/env Rscript

# Validate expected outputs for selected pipeline steps.

suppressPackageStartupMessages({
  library(utils)
})

pipeline_root <- Sys.getenv("PIPELINE_ROOT", getwd())
source(file.path(pipeline_root, "scripts", "helpers", "cli_utils.R"))

args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
outdir <- if (!is.null(args$outdir) && args$outdir != "") args$outdir else getwd()
manifest <- if (!is.null(args$manifest) && args$manifest != "") args$manifest else NULL
steps_arg <- if (!is.null(args$steps) && args$steps != "") args$steps else "01,02,03,04,05,06,07"
steps <- unique(trimws(strsplit(steps_arg, ",", fixed = TRUE)[[1]]))
strict <- if (!is.null(args$strict) && args$strict != "") tolower(args$strict) %in% c("true", "1", "yes") else TRUE
sample_ids_arg <- args$sample_ids
sample_ids_file <- args$sample_ids_file

processed_dir <- file.path(outdir, "data", "processed")
outputs_dir <- file.path(outdir, "outputs")

collect_sample_ids <- function() {
  if (!is.null(sample_ids_arg) && sample_ids_arg != "") {
    ids <- trimws(strsplit(sample_ids_arg, ",", fixed = TRUE)[[1]])
    ids <- ids[nzchar(ids)]
    return(unique(ids))
  }
  if (!is.null(sample_ids_file) && sample_ids_file != "" && file.exists(sample_ids_file)) {
    ids <- readLines(sample_ids_file, warn = FALSE)
    ids <- ids[nzchar(ids)]
    return(unique(ids))
  }
  if (!is.null(manifest) && file.exists(manifest)) {
    mf <- read.csv(manifest, stringsAsFactors = FALSE)
    if (!"sample_id" %in% colnames(mf)) stop("manifest missing sample_id column")
    ids <- mf$sample_id
    ids <- ids[!is.na(ids) & ids != ""]
    if (length(ids) == 0) stop("manifest has no sample_id values")
    return(unique(ids))
  }

  patterns <- c(
    "\\.01_ingest\\.rds$",
    "\\.02_qc\\.rds$",
    "\\.03_normalize\\.rds$",
    "\\.05_domains\\.rds$",
    "\\.06_svg\\.rds$",
    "\\.07_neighborhood\\.rds$"
  )
  files <- unlist(lapply(patterns, function(pat) {
    list.files(processed_dir, pattern = pat, full.names = FALSE)
  }))
  if (length(files) == 0) return(character(0))
  ids <- sub("\\.(01_ingest|02_qc|03_normalize|05_domains|06_svg|07_neighborhood)\\.rds$", "", files)
  unique(ids)
}

sample_ids <- collect_sample_ids()
if (length(sample_ids) == 0) {
  stop("No sample IDs found. Provide --manifest or ensure data/processed contains cached outputs.")
}

missing <- character(0)
warnings <- character(0)

check_files <- function(paths, label) {
  for (p in paths) {
    if (!file.exists(p)) missing <<- c(missing, paste(label, p, sep = ": "))
  }
}

if ("01" %in% steps) {
  paths <- file.path(processed_dir, paste0(sample_ids, ".01_ingest.rds"))
  check_files(paths, "01_ingest")
  check_files(file.path(outputs_dir, "qc", "initial_qc_summary.csv"), "01_ingest")
}

if ("02" %in% steps) {
  paths <- file.path(processed_dir, paste0(sample_ids, ".02_qc.rds"))
  check_files(paths, "02_qc")
  check_files(file.path(outputs_dir, "qc", "qc_summary_postfilter.csv"), "02_qc")
  check_files(file.path(outputs_dir, "qc", "qc_thresholds.csv"), "02_qc")
}

if ("03" %in% steps) {
  paths <- file.path(processed_dir, paste0(sample_ids, ".03_normalize.rds"))
  check_files(paths, "03_normalize")
  check_files(file.path(outputs_dir, "integration", "per_sample_norm_summary.csv"), "03_normalize")
}

if ("04" %in% steps) {
  check_files(file.path(processed_dir, "cohort.04_integrate.rds"), "04_integrate")
  check_files(file.path(outputs_dir, "integration", "cohort_integration_summary.csv"), "04_integrate")
  check_files(file.path(outputs_dir, "integration", "plots", "cohort_pca_by_sample.png"), "04_integrate")
  umap_path <- file.path(outputs_dir, "integration", "plots", "cohort_umap_by_sample.png")
  if (!file.exists(umap_path)) {
    warnings <- c(warnings, paste("04_integrate: missing optional", umap_path))
  }
}

if ("05" %in% steps) {
  paths <- file.path(processed_dir, paste0(sample_ids, ".05_domains.rds"))
  check_files(paths, "05_domains")
  check_files(file.path(outputs_dir, "domains", "domains_cohort_summary.csv"), "05_domains")
}

if ("06" %in% steps) {
  paths <- file.path(outputs_dir, "svg", paste0("svg_ranked_", sample_ids, ".csv"))
  check_files(paths, "06_svg")
  check_files(file.path(outputs_dir, "svg", "svg_cohort_summary.csv"), "06_svg")
  check_files(file.path(outputs_dir, "svg", "svg_consensus_top50.csv"), "06_svg")
}

if ("07" %in% steps) {
  mat_paths <- file.path(outputs_dir, "neighborhood", paste0("neighborhood_matrix_", sample_ids, ".csv"))
  pair_paths <- file.path(outputs_dir, "neighborhood", paste0("neighborhood_top_pairs_", sample_ids, ".csv"))
  check_files(mat_paths, "07_neighborhood")
  check_files(pair_paths, "07_neighborhood")
  check_files(file.path(outputs_dir, "neighborhood", "neighborhood_cohort_summary.csv"), "07_neighborhood")
  check_files(file.path(outputs_dir, "neighborhood", "neighborhood_consensus_pairs.csv"), "07_neighborhood")
}

validation_dir <- file.path(outputs_dir, "validation")
dir.create(validation_dir, recursive = TRUE, showWarnings = FALSE)

run_timestamp <- format(as.POSIXct(Sys.time(), tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ")
status <- if (length(missing) == 0) "pass" else "fail"

json_escape <- function(x) {
  x <- gsub("\\\\", "\\\\\\\\", x)
  x <- gsub("\"", "\\\\\"", x)
  x <- gsub("\n", "\\\\n", x)
  x
}

to_json <- function(x) {
  if (is.null(x)) return("null")
  if (is.logical(x)) return(if (length(x) == 1) ifelse(x, "true", "false") else paste0("[", paste(ifelse(x, "true", "false"), collapse = ","), "]"))
  if (is.numeric(x)) {
    if (length(x) == 1) return(ifelse(is.na(x), "null", as.character(x)))
    vals <- vapply(x, function(v) ifelse(is.na(v), "null", as.character(v)), character(1))
    return(paste0("[", paste(vals, collapse = ","), "]"))
  }
  if (is.character(x)) {
    if (length(x) == 1) return(paste0("\"", json_escape(x), "\""))
    vals <- vapply(x, function(v) paste0("\"", json_escape(v), "\""), character(1))
    return(paste0("[", paste(vals, collapse = ","), "]"))
  }
  if (is.list(x)) {
    is_named <- !is.null(names(x))
    if (!is_named) {
      vals <- vapply(x, to_json, character(1))
      return(paste0("[", paste(vals, collapse = ","), "]"))
    }
    parts <- mapply(function(k, v) {
      paste0("\"", json_escape(k), "\":", to_json(v))
    }, names(x), x, SIMPLIFY = TRUE, USE.NAMES = FALSE)
    return(paste0("{", paste(parts, collapse = ","), "}"))
  }
  "\"unknown\""
}

report <- list(
  run_timestamp = run_timestamp,
  outdir = normalizePath(outdir),
  steps = steps,
  sample_ids = sample_ids,
  missing = missing,
  warnings = warnings,
  status = status
)

report_json <- file.path(validation_dir, "validation_report.json")
report_txt <- file.path(validation_dir, "validation_report.txt")

writeLines(to_json(report), report_json)
writeLines(c(
  paste("run_timestamp:", run_timestamp),
  paste("outdir:", normalizePath(outdir)),
  paste("steps:", paste(steps, collapse = ",")),
  paste("sample_ids:", paste(sample_ids, collapse = ",")),
  paste("status:", status),
  if (length(missing) == 0) "missing: none" else "missing:",
  if (length(missing) == 0) "" else paste("- ", missing, collapse = "\n"),
  if (length(warnings) == 0) "warnings: none" else "warnings:",
  if (length(warnings) == 0) "" else paste("- ", warnings, collapse = "\n")
), report_txt)

if (length(missing) > 0 && strict) {
  message("Missing outputs:")
  message(paste("- ", missing, collapse = "\n"))
  quit(save = "no", status = 1)
}

if (length(warnings) > 0) {
  message("Warnings:")
  message(paste("- ", warnings, collapse = "\n"))
}

message("Validation complete: ", status)
