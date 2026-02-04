#!/usr/bin/env Rscript

# Write run metadata and provenance to outputs/run_metadata.*

suppressPackageStartupMessages({
  library(utils)
})

pipeline_root <- Sys.getenv("PIPELINE_ROOT", getwd())
source(file.path(pipeline_root, "scripts", "helpers", "cli_utils.R"))

args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
out_json <- if (!is.null(args$out_json) && args$out_json != "") args$out_json else "run_metadata.json"
out_txt <- if (!is.null(args$out_txt) && args$out_txt != "") args$out_txt else "run_metadata.txt"
manifest <- args$manifest
max_samples <- args$max_samples
steps <- args$steps
outdir <- args$outdir
dataset <- args$dataset
input_spe_rds <- args$input_spe_rds
sample_ids <- args$sample_ids

qc_min_counts <- args$qc_min_counts
qc_min_genes <- args$qc_min_genes
qc_max_mito <- args$qc_max_mito
qc_keep_tissue_only <- args$qc_keep_tissue_only
q_sweep <- args$q_sweep
hvg_count_norm <- args$hvg_count_norm
hvg_count_svg <- args$hvg_count_svg
knn_k_svg <- args$knn_k_svg
knn_k_nb <- args$knn_k_nb

run_timestamp <- format(as.POSIXct(Sys.time(), tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ")

get_git_commit <- function(root) {
  cmd <- sprintf("git -C %s rev-parse HEAD", shQuote(root))
  out <- suppressWarnings(tryCatch(system(cmd, intern = TRUE), error = function(e) NULL))
  if (is.null(out) || length(out) == 0) return("unknown")
  out[[1]]
}

get_nextflow_version <- function() {
  ver <- Sys.getenv("NXF_VER")
  if (!is.null(ver) && ver != "") return(ver)
  "unknown"
}

sha256_file <- function(path) {
  if (!file.exists(path)) return("unknown")
  cmd <- sprintf("sha256sum %s", shQuote(path))
  out <- suppressWarnings(tryCatch(system(cmd, intern = TRUE), error = function(e) NULL))
  if (is.null(out) || length(out) == 0) return("unknown")
  strsplit(out[[1]], "\\s+")[[1]][1]
}

get_pkg_version <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) return("unknown")
  as.character(utils::packageVersion(pkg))
}

pkg_versions <- list(
  spatialLIBD = get_pkg_version("spatialLIBD"),
  SpatialExperiment = get_pkg_version("SpatialExperiment"),
  SingleCellExperiment = get_pkg_version("SingleCellExperiment"),
  scater = get_pkg_version("scater"),
  scran = get_pkg_version("scran"),
  batchelor = get_pkg_version("batchelor"),
  BayesSpace = get_pkg_version("BayesSpace")
)

sample_ids_vec <- character(0)
if (!is.null(sample_ids) && sample_ids != "") {
  sample_ids_vec <- trimws(strsplit(sample_ids, ",", fixed = TRUE)[[1]])
  sample_ids_vec <- sample_ids_vec[nzchar(sample_ids_vec)]
}

metadata <- list(
  run_timestamp = run_timestamp,
  repo_root = normalizePath(pipeline_root),
  git_commit = get_git_commit(pipeline_root),
  nextflow_version = get_nextflow_version(),
  conda_env_file_checksum = sha256_file(file.path(pipeline_root, "environment.yml")),
  params = list(
    manifest = manifest,
    max_samples = max_samples,
    steps = steps,
    outdir = outdir,
    dataset = dataset,
    input_spe_rds = input_spe_rds,
    qc_thresholds = list(
      min_counts = qc_min_counts,
      min_genes = qc_min_genes,
      max_mito = qc_max_mito,
      keep_tissue_only = qc_keep_tissue_only
    ),
    q_sweep = q_sweep,
    hvg_count_norm = hvg_count_norm,
    hvg_count_svg = hvg_count_svg,
    knn_k_svg = knn_k_svg,
    knn_k_nb = knn_k_nb
  ),
  sample_ids = sample_ids_vec,
  n_samples = length(sample_ids_vec),
  r_session = list(
    r_version = R.version.string,
    package_versions = pkg_versions
  )
)

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

writeLines(to_json(metadata), out_json)

txt_lines <- c(
  paste("run_timestamp:", run_timestamp),
  paste("repo_root:", normalizePath(pipeline_root)),
  paste("git_commit:", metadata$git_commit),
  paste("nextflow_version:", metadata$nextflow_version),
  paste("conda_env_file_checksum:", metadata$conda_env_file_checksum),
  paste("dataset:", dataset),
  paste("input_spe_rds:", ifelse(is.null(input_spe_rds) || input_spe_rds == "", "none", input_spe_rds)),
  paste("manifest:", ifelse(is.null(manifest), "none", manifest)),
  paste("max_samples:", ifelse(is.null(max_samples), "none", max_samples)),
  paste("steps:", steps),
  paste("outdir:", outdir),
  paste("sample_ids:", ifelse(length(sample_ids_vec) == 0, "none", paste(sample_ids_vec, collapse = ","))),
  "qc_thresholds:",
  paste("  min_counts:", qc_min_counts),
  paste("  min_genes:", qc_min_genes),
  paste("  max_mito:", qc_max_mito),
  paste("  keep_tissue_only:", qc_keep_tissue_only),
  paste("q_sweep:", q_sweep),
  paste("hvg_count_norm:", hvg_count_norm),
  paste("hvg_count_svg:", hvg_count_svg),
  paste("knn_k_svg:", knn_k_svg),
  paste("knn_k_nb:", knn_k_nb),
  paste("r_version:", metadata$r_session$r_version)
)
writeLines(txt_lines, out_txt)
