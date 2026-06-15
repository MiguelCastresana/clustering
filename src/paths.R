# Shared paths for clustering benchmark scripts.

script_args <- commandArgs(trailingOnly = FALSE)
script_file_arg <- grep("^--file=", script_args, value = TRUE)
script_dir <- if (length(script_file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_file_arg[1]), mustWork = TRUE))
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

repo_root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
benchmark_root <- Sys.getenv("CLUSTERING_BENCHMARK_DIR", unset = repo_root)

input_file <- function(...) {
  file.path(repo_root, "input", ...)
}

benchmark_file <- function(...) {
  file.path(benchmark_root, ...)
}

benchmark_output <- function(...) {
  path <- benchmark_file(...)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  path
}
