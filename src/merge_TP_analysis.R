

library(dplyr)
library(stringr)
library(parallel)
library(ANUBIX)

# --- 1. Load & preprocess input ---
load_network <- function(path, cutoff = 0.75) {
  read.delim(path, header = TRUE, stringsAsFactors = FALSE) %>%
    setNames(c("score", "dummy1", "from", "to", tail(names(.), -4))) %>%
    filter(score >= cutoff) %>%
    select(from, to, score)
}

load_pathways <- function(path) {
  read.delim(path, header = TRUE, stringsAsFactors = FALSE)
}

load_genesets <- function(path) {
  read.delim(path, header = FALSE, stringsAsFactors = FALSE) %>%
    setNames(c("gene", "geneset"))
}

# --- 2. Run constrained ANUBIX enrichment ---

run_anubix_fp <- function(network, pathways, genesets,
                          sampling     = 2000,
                          cutoff       = 0.75,
                          cores        = 2,
                          network_type = "weighted") {
  # build link matrix
  links <- anubix_links(
    network      = network,
    pathways     = pathways,
    cutoff       = cutoff,
    network_type = network_type
  )
  
  # perform constrained enrichment and return
  result <- anubix(
    network       = network,
    links_matrix  = links,
    genesets      = genesets,
    pathways      = pathways,
    sampling      = sampling,
    cores         = cores,
    cutoff        = cutoff,
    network_type  = network_type
  )
  return(result)
}

# --- 3. Helper to clean & min‑p summarize ---
clean_and_minp <- function(df, set_col, path_col, p_col, extra_cols = NULL) {
  df %>%
    mutate(
      set_clean  = str_remove_all({{set_col}},  "[^[:alnum:]]"),
      path_clean = str_remove_all({{path_col}}, "[^[:alnum:]]")
    ) %>%
    filter(str_detect(set_clean, path_clean)) %>%
    group_by(across(all_of(c(deparse(substitute(set_col)), deparse(substitute(path_col)))))) %>%
    slice_min(order_by = {{p_col}}, n = 1) %>%
    ungroup() %>%
    select(-set_clean, -path_clean)
}

# --- 4. Process other methods ---
process_binox <- function(path) {
  df <- read.delim(path, header = TRUE, stringsAsFactors = FALSE)
  df %>%
    filter(rel == "+") %>%
    mutate(p_adj = p.adjust(p.value, method = "BH")) %>%
    filter(p_adj <= 0.05) %>%
    clean_and_minp(df = ., set_col = geneset, path_col = pathway, p_col = p.value) %>%
    transmute(
      geneset,
      pathway,
      CB.fdr = p_adj,
      sizeA,
      sizeB,
      rel
    )
}

process_neat <- function(path) {
  df <- read.delim(path, header = TRUE, stringsAsFactors = FALSE)
  df %>%
    filter(nab > expected_nab) %>%
    mutate(p_adj = p.adjust(nab.pval, method = "BH")) %>%
    filter(p_adj <= 0.05) %>%
    clean_and_minp(df = ., set_col = geneset, path_col = pathway, p_col = nab.pval) %>%
    transmute(
      geneset,
      pathway,
      NEAT.fdr = p_adj,
      nab,
      expected_nab
    )
}

process_gea <- function(path, universe = 20000) {
  df <- read.delim(path, header = TRUE, stringsAsFactors = FALSE)
  df <- df %>%
    mutate(
      pval = map2_dbl(
        overlap, setSize,
        ~ fisher.test(matrix(c(.x, .y - .x, universe, universe), 2, 2))$p.value
      ),
      p_adj = p.adjust(pval, method = "BH")
    ) %>%
    filter(p_adj <= 0.05) %>%
    clean_and_minp(df = ., set_col = geneset, path_col = pathway, p_col = pval) %>%
    transmute(
      geneset,
      pathway,
      GEA.fdr = p_adj,
      GEA.k = overlap
    )
}

# --- 5. Merge all results ---
merge_results <- function(tp, binox = NULL, neat = NULL, gea = NULL) {
  result <- tp
  for (tbl in list(binox, neat, gea)) {
    if (!is.null(tbl)) {
      result <- full_join(result, tbl, by = c("geneset", "pathway"))
    }
  }
  result
}

# --- 6. Main execution for TP (no clustering) ---
network   <- load_network("/scratch/2020_clustering/fc4.1")
pathways  <- load_pathways("/scratch/2020_clustering/KEGG_h_sapiens")
genesets  <- load_genesets("/scratch/2020_clustering/TP/TP_genesets")

tp_noclust <- run_anubix_tp(network, pathways, genesets)
binox_nocl <- process_binox("/scratch/2020_clustering/benchmark/TP/binox_TP_nocluster")
neat_nocl  <- process_neat("/scratch/2020_clustering/benchmark/TP/neat_TP_nocluster")
gea_nocl   <- process_gea("/scratch/2020_clustering/benchmark/TP/gea_TP_nocluster")

merged_nocl <- merge_results(tp_noclust, binox_nocl, neat_nocl, gea_nocl)
write.table(merged_nocl,
            "/scratch/2020_clustering/FINAL/all_merge_TP_NOcluster.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- 7. Repeat for clustering methods ---
cluster_methods <- c(infomap = "infomap_tp",
                     mcl     = "mcl_tp",
                     mgclus  = "mgclus_tp")

for (method in names(cluster_methods)) {
  gs_path  <- file.path("/scratch/2020_clustering/benchmark/TP/", paste0(cluster_methods[[method]], ".tsv"))
  gs       <- load_genesets(gs_path)
  tp_res   <- run_anubix_tp(network, pathways, gs)
  bx       <- process_binox(sprintf("/scratch/2020_clustering/benchmark/TP/results/binox_%s.tsv",      method))
  nt       <- process_neat(sprintf("/scratch/2020_clustering/benchmark/TP/results/neat_%s.tsv",        method))
  ga       <- process_gea(sprintf("/scratch/2020_clustering/benchmark/TP/results/gea_%s.tsv",          method))
  merged   <- merge_results(tp_res, bx, nt, ga)
  out_path <- sprintf("/scratch/2020_clustering/benchmark/TP/results/all_merge_TP_cluster_%s.tsv", method)
  write.table(merged, out_path, sep = "\t", row.names = FALSE, quote = FALSE)
}
