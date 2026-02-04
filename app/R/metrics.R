metrics_env <- new.env(parent = emptyenv())

set_metrics_context <- function(ctx) {
  metrics_env$ctx <- ctx
}

get_metrics_context <- function() {
  ctx <- metrics_env$ctx
  if (is.null(ctx)) stop("Metrics context not set.")
  ctx
}

safe_read_csv <- function(path) {
  if (!file.exists(path)) return(NULL)
  read.csv(path, stringsAsFactors = FALSE)
}

safe_read_rds <- function(path) {
  if (!file.exists(path)) return(NULL)
  readRDS(path)
}

resolve_domains_rds <- function(processed_dir, sample_id) {
  new_path <- file.path(processed_dir, paste0(sample_id, ".05_domains.rds"))
  legacy_path <- file.path(processed_dir, paste0(sample_id, "_domains.rds"))
  if (file.exists(new_path)) return(new_path)
  if (file.exists(legacy_path)) return(legacy_path)
  new_path
}

validate_domains_object <- function(spe) {
  errors <- character(0)
  if (is.null(spe)) errors <- c(errors, "Domain object not loaded.")
  if (!"domain_main" %in% colnames(SummarizedExperiment::colData(spe))) {
    errors <- c(errors, "domain_main missing in colData.")
  }
  coords <- SpatialExperiment::spatialCoords(spe)
  if (is.null(coords) || nrow(coords) == 0) {
    errors <- c(errors, "Spatial coordinates missing.")
  }
  if (!"sample_id" %in% colnames(SummarizedExperiment::colData(spe))) {
    errors <- c(errors, "sample_id missing in colData.")
  }
  errors
}

validate_svg_table <- function(df) {
  if (is.null(df)) return("SVG table not found.")
  required <- c("gene", "rank")
  if (!all(required %in% colnames(df))) return("SVG table missing required columns.")
  if (!("score" %in% colnames(df) || "moran_I" %in% colnames(df))) {
    return("SVG table missing score or moran_I column.")
  }
  NULL
}

validate_neighborhood_matrix <- function(mat, tol = 1e-6) {
  if (is.null(mat)) return("Neighborhood matrix not found.")
  if (nrow(mat) != ncol(mat)) return("Neighborhood matrix is not square.")
  if (max(abs(mat - t(mat)), na.rm = TRUE) > tol) return("Neighborhood matrix is not symmetric.")
  NULL
}

compute_domains_metrics <- function(sample_id) {
  ctx <- get_metrics_context()
  dom_path <- resolve_domains_rds(ctx$processed_dir, sample_id)
  stab_path <- file.path(ctx$domains_dir, paste0(sample_id, "_stability.csv"))

  if (!file.exists(dom_path)) {
    message("Metrics error: missing ", dom_path)
    return(list(ok = FALSE, error = paste("Missing", dom_path), command = "Rscript scripts/05_domains.R", metrics = NULL))
  }

  spe <- safe_read_rds(dom_path)
  errs <- validate_domains_object(spe)
  if (length(errs) > 0) {
    message("Metrics error: ", paste(errs, collapse = " "))
    return(list(ok = FALSE, error = paste(errs, collapse = " "), command = "Rscript scripts/05_domains.R", metrics = NULL))
  }

  domains <- SummarizedExperiment::colData(spe)$domain_main
  counts <- table(domains)
  n_spots <- ncol(spe)
  n_domains <- length(counts)
  largest_domain_pct <- round(max(counts) / n_spots * 100, 2)

  stability_score <- NA_real_
  if (file.exists(stab_path)) {
    stab <- safe_read_csv(stab_path)
    if (!is.null(stab) && "stability_score" %in% colnames(stab)) {
      stability_score <- mean(stab$stability_score, na.rm = TRUE)
    }
  }

  list(ok = TRUE, metrics = list(
    n_spots = n_spots,
    n_domains = n_domains,
    largest_domain_pct = largest_domain_pct,
    stability_score = stability_score
  ))
}

compute_svg_metrics <- function(sample_id) {
  ctx <- get_metrics_context()
  svg_path <- file.path(ctx$svg_dir, paste0("svg_ranked_", sample_id, ".csv"))
  cohort_path <- file.path(ctx$svg_dir, "svg_cohort_summary.csv")

  df <- safe_read_csv(svg_path)
  err <- validate_svg_table(df)
  if (!is.null(err)) {
    message("Metrics error: ", err)
    return(list(ok = FALSE, error = paste("Missing or invalid SVG table at", svg_path, ".", err), command = "Rscript scripts/06_svg.R", metrics = NULL))
  }

  score_col <- if ("score" %in% colnames(df)) "score" else "moran_I"
  top_gene <- df$gene[1]
  top_score <- df[[score_col]][1]
  n_genes_tested <- nrow(df)

  consensus_hit <- NA_integer_
  cohort <- safe_read_csv(cohort_path)
  if (!is.null(cohort) && "gene" %in% colnames(cohort) && "n_samples_in_top50" %in% colnames(cohort)) {
    idx <- match(top_gene, cohort$gene)
    if (!is.na(idx)) consensus_hit <- cohort$n_samples_in_top50[idx]
  }

  list(ok = TRUE, metrics = list(
    top_gene = top_gene,
    top_score = round(top_score, 4),
    n_genes_tested = n_genes_tested,
    consensus_hit = consensus_hit
  ))
}

compute_neighborhood_metrics <- function(sample_id) {
  ctx <- get_metrics_context()
  mat_path <- file.path(ctx$nb_dir, paste0("neighborhood_matrix_", sample_id, ".csv"))
  top_path <- file.path(ctx$nb_dir, paste0("neighborhood_top_pairs_", sample_id, ".csv"))
  cohort_path <- file.path(ctx$nb_dir, "neighborhood_cohort_summary.csv")

  mat <- NULL
  if (file.exists(mat_path)) {
    mat <- as.matrix(read.csv(mat_path, row.names = 1, check.names = FALSE))
  }
  err <- validate_neighborhood_matrix(mat)
  if (!is.null(err)) {
    message("Metrics error: ", err)
    return(list(ok = FALSE, error = paste("Missing or invalid neighborhood matrix at", mat_path, ".", err), command = "Rscript scripts/07_neighborhood.R", metrics = NULL))
  }

  density <- round(sum(mat > 0) / (nrow(mat) * ncol(mat)), 4)

  top_df <- safe_read_csv(top_path)
  top_pair <- if (!is.null(top_df) && nrow(top_df) > 0) top_df$domain_pair[1] else NA_character_
  top_pair_score <- if (!is.null(top_df) && nrow(top_df) > 0) top_df$score[1] else NA_real_

  consensus_count <- NA_integer_
  cohort <- safe_read_csv(cohort_path)
  if (!is.null(cohort) && "domain_pair" %in% colnames(cohort) && "n_samples_top25" %in% colnames(cohort)) {
    idx <- match(top_pair, cohort$domain_pair)
    if (!is.na(idx)) consensus_count <- cohort$n_samples_top25[idx]
  }

  list(ok = TRUE, metrics = list(
    top_pair = top_pair,
    top_pair_score = top_pair_score,
    matrix_density = density,
    consensus_count = consensus_count
  ))
}

compute_cohort_metrics <- function() {
  ctx <- get_metrics_context()
  dom_sum <- safe_read_csv(file.path(ctx$domains_dir, "domains_cohort_summary.csv"))
  svg_sum <- safe_read_csv(file.path(ctx$svg_dir, "svg_cohort_summary.csv"))
  nb_sum <- safe_read_csv(file.path(ctx$nb_dir, "neighborhood_cohort_summary.csv"))

  n_samples <- if (!is.null(dom_sum) && "sample_id" %in% colnames(dom_sum)) length(unique(dom_sum$sample_id)) else 0
  mean_stability <- if (!is.null(dom_sum) && "stability_score" %in% colnames(dom_sum)) {
    round(mean(dom_sum$stability_score, na.rm = TRUE), 3)
  } else {
    NA_real_
  }
  stability_range <- if (!is.null(dom_sum) && "stability_score" %in% colnames(dom_sum)) {
    rng <- range(dom_sum$stability_score, na.rm = TRUE)
    paste0(round(rng[1], 3), " to ", round(rng[2], 3))
  } else {
    NA_character_
  }
  top_consensus_pair <- if (!is.null(nb_sum) && nrow(nb_sum) > 0) nb_sum$domain_pair[1] else NA_character_
  top_consensus_gene <- if (!is.null(svg_sum) && nrow(svg_sum) > 0) svg_sum$gene[1] else NA_character_

  domains_median <- if (!is.null(dom_sum) && "n_domains" %in% colnames(dom_sum)) {
    median(dom_sum$n_domains, na.rm = TRUE)
  } else {
    NA_real_
  }

  list(ok = TRUE, metrics = list(
    n_samples = n_samples,
    mean_stability = mean_stability,
    stability_range = stability_range,
    top_consensus_pair = top_consensus_pair,
    top_consensus_gene = top_consensus_gene,
    domains_median = domains_median
  ))
}
