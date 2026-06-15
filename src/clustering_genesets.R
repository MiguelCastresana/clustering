# Cluster gene sets with MGclus, MCL, and Infomap.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
project_root <- if (length(file_arg) > 0) {
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."), mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

project_file <- function(...) {
  file.path(project_root, ...)
}

input_dir <- project_file("cluster_algorithms")
output_dir <- project_file("results", "clusters")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

input_rdata <- project_file("input", "input_TP_genesets.RData")
if (!file.exists(input_rdata)) {
  stop("Missing input file: ", input_rdata, call. = FALSE)
}
load(input_rdata)

if (!exists("kegg_kegg")) {
  stop("Expected object `kegg_kegg` in ", input_rdata, call. = FALSE)
}

write_modules <- function(clusters, source_sets, output_path) {
  names(clusters) <- names(source_sets)
  modules <- data.frame(gene = character(), module = character())

  for (i in seq_along(source_sets)) {
    if (length(clusters[[i]]) == 0) {
      next
    }

    module_index <- 1
    for (cluster in clusters[[i]]) {
      genes <- unlist(cluster, use.names = FALSE)
      if (length(genes) == 0) {
        next
      }

      pathway_name <- names(source_sets)[i]
      pathway_name <- gsub("([[:punct:]])|\\s+", "_", pathway_name)
      pathway_name <- gsub("[[:blank:]]", "", pathway_name)
      module_id <- toupper(paste0("g_", module_index, "_#", pathway_name))

      modules <- rbind(
        modules,
        data.frame(gene = genes, module = module_id, stringsAsFactors = FALSE)
      )
      module_index <- module_index + 1
    }
  }

  modules <- modules[modules$module %in% names(which(table(modules$module) > 2)), ]
  write.table(modules, output_path, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  invisible(modules)
}

run_mgclus <- function(source_sets) {
  mgclus_jar <- project_file("cluster_algorithms", "mgclusjar.jar")
  if (!file.exists(mgclus_jar)) {
    stop("Missing MGclus jar: ", mgclus_jar, call. = FALSE)
  }

  one_kegg_path <- file.path(input_dir, "onekegg.tsv")
  cluster_path <- file.path(input_dir, "onekeggcluster")
  clusters <- vector("list", length(source_sets))

  for (i in seq_along(source_sets)) {
    write.table(source_sets[[i]], one_kegg_path, sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE)
    status <- system2(
      "java",
      c("-jar", mgclus_jar, "-f", one_kegg_path, "-w", "T", "-o", cluster_path)
    )
    if (status != 0) {
      stop("MGclus failed for gene set index ", i, call. = FALSE)
    }
    clusters[[i]] <- strsplit(scan(cluster_path, what = "", sep = "\n", quiet = TRUE), "[[:space:]]+")
  }

  write_modules(clusters, source_sets, file.path(output_dir, "mgclus.tsv"))
}

run_mcl <- function(source_sets) {
  mcl_executable <- project_file("cluster_algorithms", "mcl", "mcl-14-137", "src", "shmcl", "mcl")
  if (!file.exists(mcl_executable)) {
    stop("Missing MCL executable: ", mcl_executable, call. = FALSE)
  }

  one_kegg_path <- file.path(input_dir, "onekegg.tsv")
  cluster_path <- file.path(input_dir, "onekeggcluster.mcl")
  clusters <- vector("list", length(source_sets))

  for (i in seq_along(source_sets)) {
    write.table(source_sets[[i]], one_kegg_path, sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE)
    status <- system2(mcl_executable, c(one_kegg_path, "--abc", "-o", cluster_path))
    if (status != 0) {
      stop("MCL failed for gene set index ", i, call. = FALSE)
    }
    clusters[[i]] <- strsplit(scan(cluster_path, what = "", sep = "\n", quiet = TRUE), "[[:space:]]+")
  }

  write_modules(clusters, source_sets, file.path(output_dir, "mcl.tsv"))
}

run_infomap <- function(source_sets) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("The `igraph` package is required for Infomap clustering.", call. = FALSE)
  }

  clusters <- vector("list", length(source_sets))
  for (i in seq_along(source_sets)) {
    graph <- igraph::graph_from_data_frame(source_sets[[i]][, 1:2], directed = FALSE)
    weights <- if (ncol(source_sets[[i]]) >= 3) as.vector(unlist(source_sets[[i]][, 3])) else NULL
    communities <- igraph::cluster_infomap(graph, e.weights = weights)
    clusters[[i]] <- igraph::communities(communities)
  }

  write_modules(clusters, source_sets, file.path(output_dir, "infomap.tsv"))
}

run_mgclus(kegg_kegg)
run_mcl(kegg_kegg)
run_infomap(kegg_kegg)
