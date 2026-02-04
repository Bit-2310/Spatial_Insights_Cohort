parse_cli_args <- function(args) {
  out <- list()
  if (length(args) == 0) return(out)
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (grepl("^--", key)) {
      val <- if (i + 1 <= length(args)) args[[i + 1]] else NA_character_
      out[[sub("^--", "", key)]] <- val
      i <- i + 2
    } else {
      i <- i + 1
    }
  }
  out
}

log_line <- function(log_file, msg) {
  if (is.null(log_file) || is.na(log_file) || log_file == "") return(invisible(NULL))
  dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
  cat(paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ", msg, "\n"),
      file = log_file, append = TRUE)
}

write_summary <- function(step, sample_id, spots_in, spots_out, runtime_sec, status, out_outputs) {
  summary_dir <- file.path(out_outputs, step, "summary")
  dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
  out_path <- file.path(summary_dir, paste0(sample_id, ".tsv"))
  df <- data.frame(
    sample_id = sample_id,
    spots_in = spots_in,
    spots_out = spots_out,
    runtime_sec = round(runtime_sec, 2),
    status = status,
    stringsAsFactors = FALSE
  )
  write.table(df, out_path, sep = "\t", row.names = FALSE, quote = FALSE)
}
