#!/usr/bin/env Rscript

# Download raw 10x Visium breast cancer datasets into data/raw/breast_cancer_visium/
# Base R only + utils + parallel

suppressPackageStartupMessages({
  library(utils)
  library(parallel)
})

root_dir <- file.path("data", "raw", "breast_cancer_visium")
dir.create(root_dir, recursive = TRUE, showWarnings = FALSE)

get_cores <- function() {
  opt <- getOption("N_CORES")
  if (!is.null(opt) && !is.na(opt)) {
    return(as.integer(opt))
  }
  cores <- parallel::detectCores()
  if (is.na(cores)) cores <- 1L
  min(8L, cores)
}

samples <- data.frame(
  sample_id = c("BC_01", "BC_02"),
  dataset_name = c(
    "V1_Breast_Cancer_Block_A_Section_1",
    "V1_Breast_Cancer_Block_A_Section_2"
  ),
  description = c(
    "Human breast cancer, Block A Section 1, Visium Spatial Gene Expression.",
    "Human breast cancer, Block A Section 2, Visium Spatial Gene Expression."
  ),
  stringsAsFactors = FALSE
)

base_url <- "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0"
samples$urls <- I(lapply(samples$dataset_name, function(ds) {
  c(
    file.path(base_url, ds, paste0(ds, "_filtered_feature_bc_matrix.tar.gz")),
    file.path(base_url, ds, paste0(ds, "_spatial.tar.gz"))
  )
}))
samples$expected_contents <- I(lapply(seq_len(nrow(samples)), function(i) {
  c("filtered_feature_bc_matrix", "spatial")
}))

expand_to_40 <- function(df) {
  if (nrow(df) >= 40) return(df)
  extra_ids <- sprintf("BC_%02d", (nrow(df) + 1):40)
  extra <- data.frame(
    sample_id = extra_ids,
    dataset_name = rep("TBD", length(extra_ids)),
    description = rep("Reserved for future dataset.", length(extra_ids)),
    stringsAsFactors = FALSE
  )
  extra$urls <- I(lapply(seq_len(nrow(extra)), function(i) character(0)))
  extra$expected_contents <- I(lapply(seq_len(nrow(extra)), function(i) c("filtered_feature_bc_matrix", "spatial")))
  rbind(df, extra)
}

samples <- expand_to_40(samples)

download_and_extract <- function(sample_row) {
  sample_id <- sample_row$sample_id
  urls <- sample_row$urls[[1]]
  expected <- sample_row$expected_contents[[1]]
  out_dir <- file.path(root_dir, sample_id)

  if (dir.exists(out_dir) && length(list.files(out_dir)) > 0) {
    return(data.frame(sample_id = sample_id, source_url = paste(urls, collapse = ";"),
                      status = "skipped", message = "directory exists", stringsAsFactors = FALSE))
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if (length(urls) == 0) {
    return(data.frame(sample_id = sample_id, source_url = "", status = "skipped",
                      message = "no URL assigned", stringsAsFactors = FALSE))
  }

  for (u in urls) {
    tmp <- tempfile(fileext = if (grepl("\\.zip$", u)) ".zip" else ".tar.gz")
    message("Downloading ", sample_id, ": ", u)
    ok <- tryCatch({
      utils::download.file(u, destfile = tmp, mode = "wb", quiet = TRUE)
      TRUE
    }, error = function(e) {
      message("Download failed for ", sample_id, ": ", conditionMessage(e))
      FALSE
    })
    if (!ok) {
      return(data.frame(sample_id = sample_id, source_url = u, status = "failed",
                        message = "download failed", stringsAsFactors = FALSE))
    }

    message("Extracting ", sample_id, ": ", basename(tmp))
    if (grepl("\\.zip$", tmp)) {
      utils::unzip(tmp, exdir = out_dir)
    } else {
      utils::untar(tmp, exdir = out_dir)
    }
  }

  # If expected directories are nested, move them to top-level
  for (name in expected) {
    target <- file.path(out_dir, name)
    if (!dir.exists(target)) {
      matches <- list.files(out_dir, pattern = paste0("^", name, "$"), recursive = TRUE, full.names = TRUE)
      if (length(matches) > 0) {
        dir.create(target, showWarnings = FALSE)
        file.rename(matches[1], target)
      }
    }
  }

  missing <- expected[!dir.exists(file.path(out_dir, expected))]
  if (length(missing) > 0) {
    return(data.frame(sample_id = sample_id, source_url = paste(urls, collapse = ";"),
                      status = "failed", message = paste("missing:", paste(missing, collapse = ",")),
                      stringsAsFactors = FALSE))
  }

  data.frame(sample_id = sample_id, source_url = paste(urls, collapse = ";"),
             status = "success", message = "downloaded", stringsAsFactors = FALSE)
}

cores <- get_cores()
message("Using cores: ", cores)

is_unix <- .Platform$OS.type == "unix"
results <- if (is_unix) {
  parallel::mclapply(split(samples, seq_len(nrow(samples))), download_and_extract, mc.cores = cores)
} else {
  lapply(split(samples, seq_len(nrow(samples))), download_and_extract)
}

summary_df <- do.call(rbind, results)
summary_path <- file.path(root_dir, "download_summary.csv")
write.csv(summary_df, summary_path, row.names = FALSE)

readme_path <- file.path(root_dir, "README.md")
readme_lines <- c(
  "# Breast cancer Visium raw data",
  "",
  "Mapping of internal IDs to source datasets:",
  ""
)

for (i in seq_len(nrow(samples))) {
  sid <- samples$sample_id[i]
  name <- samples$dataset_name[i]
  desc <- samples$description[i]
  urls <- samples$urls[[i]]
  if (length(urls) == 0) {
    readme_lines <- c(readme_lines, paste0("- ", sid, " → ", name, " | ", desc))
  } else {
    readme_lines <- c(readme_lines, paste0("- ", sid, " → ", name, " | ", desc),
                      paste0("  URL: ", paste(urls, collapse = " ; ")))
  }
}

writeLines(readme_lines, readme_path)

message("Summary written to: ", summary_path)
message("README written to: ", readme_path)
message("Done.")
