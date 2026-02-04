load_csv_safe <- function(path) {
  if (!file.exists(path)) return(NULL)
  read.csv(path, stringsAsFactors = FALSE)
}

load_rds_safe <- function(path) {
  if (!file.exists(path)) return(NULL)
  readRDS(path)
}

list_samples_from_domains <- function(processed_dir) {
  files_new <- list.files(processed_dir, pattern = "\\.05_domains\\.rds$", full.names = TRUE)
  files_legacy <- list.files(processed_dir, pattern = "_domains\\.rds$", full.names = TRUE)
  files <- unique(c(files_new, files_legacy))
  if (length(files) == 0) return(character(0))
  base <- basename(files)
  base <- sub("\\.05_domains\\.rds$", "", base)
  sub("_domains\\.rds$", "", base)
}

list_samples_from_svg <- function(svg_dir) {
  files <- list.files(svg_dir, pattern = "^svg_ranked_.*\\.csv$", full.names = TRUE)
  if (length(files) == 0) return(character(0))
  sub("^svg_ranked_(.*)\\.csv$", "\\1", basename(files))
}

get_sample_ids <- function(processed_dir, svg_dir) {
  s1 <- list_samples_from_domains(processed_dir)
  s2 <- list_samples_from_svg(svg_dir)
  samples <- unique(c(s1, s2))
  samples[order(samples)]
}

resolve_domains_rds <- function(processed_dir, sample_id) {
  new_path <- file.path(processed_dir, paste0(sample_id, ".05_domains.rds"))
  legacy_path <- file.path(processed_dir, paste0(sample_id, "_domains.rds"))
  if (file.exists(new_path)) return(new_path)
  if (file.exists(legacy_path)) return(legacy_path)
  new_path
}

resolve_norm_rds <- function(processed_dir, sample_id) {
  new_path <- file.path(processed_dir, paste0(sample_id, ".03_normalize.rds"))
  legacy_path <- file.path(processed_dir, paste0(sample_id, "_norm.rds"))
  if (file.exists(new_path)) return(new_path)
  if (file.exists(legacy_path)) return(legacy_path)
  new_path
}

get_svg_genes <- function(svg_dir, sample_id) {
  path <- file.path(svg_dir, paste0("svg_ranked_", sample_id, ".csv"))
  df <- load_csv_safe(path)
  if (is.null(df) || nrow(df) == 0) return(character(0))
  df$gene
}

get_default_gene <- function(sample_id, svg_dir) {
  path <- file.path(svg_dir, paste0("svg_ranked_", sample_id, ".csv"))
  df <- load_csv_safe(path)
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df$gene[1]
}

get_default_pair <- function(sample_id, nb_dir) {
  path <- file.path(nb_dir, paste0("neighborhood_top_pairs_", sample_id, ".csv"))
  df <- load_csv_safe(path)
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df$domain_pair[1]
}

load_domain_frame <- function(domain_rds) {
  spe <- readRDS(domain_rds)
  coords <- SpatialExperiment::spatialCoords(spe)
  if (is.null(coords)) return(NULL)
  data.frame(
    x = coords[, 1],
    y = coords[, 2],
    domain = as.factor(SummarizedExperiment::colData(spe)$domain_main),
    barcode = colnames(spe),
    sample_id = as.character(SummarizedExperiment::colData(spe)$sample_id),
    stringsAsFactors = FALSE
  )
}

load_gene_values <- function(sample_id, processed_dir, gene) {
  norm_path <- resolve_norm_rds(processed_dir, sample_id)
  dom_path <- resolve_domains_rds(processed_dir, sample_id)

  spe <- NULL
  if (file.exists(dom_path)) spe <- readRDS(dom_path)
  if (is.null(spe) || !"logcounts" %in% SummarizedExperiment::assayNames(spe)) {
    spe <- readRDS(norm_path)
  }
  if (is.null(spe) || !"logcounts" %in% SummarizedExperiment::assayNames(spe)) return(NULL)
  if (!(gene %in% rownames(spe))) return(NULL)

  coords <- SpatialExperiment::spatialCoords(spe)
  if (is.null(coords)) return(NULL)

  expr <- SummarizedExperiment::assay(spe, "logcounts")[gene, ]
  data.frame(
    x = coords[, 1],
    y = coords[, 2],
    value = as.numeric(expr),
    barcode = colnames(spe),
    stringsAsFactors = FALSE
  )
}

load_svg_ranked <- function(sample_id, svg_dir) {
  path <- file.path(svg_dir, paste0("svg_ranked_", sample_id, ".csv"))
  df <- load_csv_safe(path)
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df
}

load_svg_score_frame <- function(sample_id, processed_dir, svg_dir, top_n) {
  df <- load_svg_ranked(sample_id, svg_dir)
  if (is.null(df)) return(NULL)
  top_genes <- df$gene[seq_len(min(top_n, nrow(df)))]

  norm_path <- resolve_norm_rds(processed_dir, sample_id)
  spe <- readRDS(norm_path)
  if (is.null(spe) || !"logcounts" %in% SummarizedExperiment::assayNames(spe)) return(NULL)

  coords <- SpatialExperiment::spatialCoords(spe)
  if (is.null(coords)) return(NULL)

  expr <- SummarizedExperiment::assay(spe, "logcounts")[top_genes, , drop = FALSE]
  expr <- as.matrix(expr)
  z <- t(scale(t(expr)))
  score <- colMeans(z, na.rm = TRUE)
  cutoff <- quantile(score, 0.9, na.rm = TRUE)

  data.frame(
    x = coords[, 1],
    y = coords[, 2],
    score = as.numeric(score),
    hotspot = score >= cutoff,
    barcode = colnames(spe),
    stringsAsFactors = FALSE
  )
}

load_neighborhood_matrix <- function(sample_id, nb_dir) {
  path <- file.path(nb_dir, paste0("neighborhood_matrix_", sample_id, ".csv"))
  if (!file.exists(path)) return(NULL)
  mat <- read.csv(path, row.names = 1, check.names = FALSE)
  as.matrix(mat)
}

read_run_metadata <- function(project_root) {
  json_path <- file.path(project_root, "outputs", "run_metadata.json")
  if (!file.exists(json_path)) return(NULL)
  if (!requireNamespace("jsonlite", quietly = TRUE)) return(NULL)
  tryCatch(jsonlite::fromJSON(json_path), error = function(e) NULL)
}
