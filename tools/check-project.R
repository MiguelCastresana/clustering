#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
project_root <- if (length(file_arg) > 0) {
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."), mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

required_paths <- c(
  "README.md",
  "src",
  "input",
  "results",
  "cluster_algorithms",
  "src/clustering_genesets.R",
  "src/roc_curves.R"
)

missing <- required_paths[!file.exists(file.path(project_root, required_paths))]
if (length(missing) > 0) {
  stop("Missing required paths:\n", paste(missing, collapse = "\n"), call. = FALSE)
}

r_files <- list.files(file.path(project_root, "src"), pattern = "\\.R$", recursive = TRUE, full.names = TRUE)

parse_one <- function(path) {
  tryCatch(
    {
      parse(path)
      TRUE
    },
    error = function(err) {
      message("Parse failed: ", path)
      message(conditionMessage(err))
      FALSE
    }
  )
}

ok <- vapply(r_files, parse_one, logical(1))
if (!all(ok)) {
  stop("One or more R files failed to parse.", call. = FALSE)
}

message("Project structure and R syntax checks passed.")
