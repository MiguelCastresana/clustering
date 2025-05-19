
# Load libraries
library(dplyr)
library(ANUBIX) 

# 1. Load network (filter by cutoff)
load_network <- function(net_path, cutoff = 0.75) {
  net <- read.delim(net_path, header = TRUE, stringsAsFactors = FALSE)
  colnames(net)[1:3] <- c("from", "to", "score")
  net %>%
    filter(score >= cutoff) %>%
    select(from, to, score)
}

# 2. Load pathways file
load_pathways <- function(pathway_path) {
  read.delim(pathway_path, header = TRUE, stringsAsFactors = FALSE)
}
# Load genesets file
load_genesets <- function(path) {
  df <- read.delim(path, header = FALSE, stringsAsFactors = FALSE)
  colnames(df) <- c("gene", "geneset")
  df
}

# 3. Run ANUBIX constrained enrichment
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
  
  # perform constrained enrichment
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

# 4. Clean and merge ANUBIX, BinoX, NEAT, GEA results
merge_enrichment_results <- function(annx_res, binox_path = NULL,
                                     neat_path = NULL, gea_path = NULL) {
  # ANUBIX
  anubix_df <- as.data.frame(annx_res)

  
  # BinoX
  if (!is.null(binox_path)) {
    cb <- read.delim(binox_path, header = TRUE, stringsAsFactors = FALSE)
    cb <- cb %>%
      filter(rel == "+") %>%
      mutate(CB.fdr = p.adjust(p.value, method = "BH")) %>%
      select(geneset, pathway, CB.fdr, sizeA, sizeB, rel)
  } else {
    cb <- NULL
  }
  
  # NEAT
  if (!is.null(neat_path)) {
    neat <- read.delim(neat_path, header = TRUE, stringsAsFactors = FALSE)
    neat <- neat %>%
      filter(nab > expected_nab) %>%
      mutate(NEAT.fdr = p.adjust(nab.pval, method = "BH")) %>%
      select(geneset, pathway, NEAT.fdr, nab, expected_nab)
  } else {
    neat <- NULL
  }
  
  # GEA
  if (!is.null(gea_path)) {
    gea_raw <- read.delim(gea_path, header = TRUE, stringsAsFactors = FALSE)
    pvals <- apply(gea_raw, 1, function(row) {
      m <- as.numeric(row[6]) - 1
      n <- as.numeric(row[4]) - as.numeric(row[6])
      k <- as.numeric(row[5])
      total <- 20000 - k
      fisher.test(matrix(c(m, n, k, total), nrow = 2))$p.value
    })
    gea <- gea_raw %>%
      mutate(GEA.fdr = p.adjust(pvals, method = "BH")) %>%
      select(geneset, pathway, GEA.fdr, k)
  } else {
    gea <- NULL
  }
  
  # Merge all results
  result <- anubix_df
  for (tab in list(cb, neat, gea)) {
    if (!is.null(tab)) {
      result <- full_join(result, tab, by = c("geneset", "pathway"))
    }
  }
  result
}

# ===== Example Execution =====
# File paths
network_path <- "/scratch/2020_clustering/fc4.1"
pathway_path <- "/scratch/2020_clustering/KEGG_h_sapiens"
genesets_path <- "/scratch/2020_clustering/FP/FP_genesets"
# 4.1 No clustering
net        <- load_network(network_path)
paths      <- load_pathways(pathway_path)
genesets      <- load_genesets(genesets_path)
no_annx    <- run_anubix_fp(net, paths, genesets)
no_merged  <- merge_enrichment_results(
  annx_res   = no_annx,
  binox_path = "/scratch/2020_clustering/benchmark/FP/results/binox_FP_nocluster_biased_50",
  neat_path  = "/scratch/2020_clustering/benchmark/FP/results/170520FP_neat_nocluster_50.tsv",
  gea_path   = "/scratch/2020_clustering/benchmark/FP/results/gea_fp_nocluster"
)
write.table(no_merged, "/scratch/2020_clustering/benchmark/FP/results/all_merge_nocluster_FP_50",
            sep = "\t", row.names = FALSE, quote = FALSE)

# 4.2 Clustering methods
cluster_files <- list(
  infomap = "/scratch/2020_clustering/benchmark/FP/FINAL_infomap_FP_50_allclusters",
  mcl     = "/scratch/2020_clustering/benchmark/FP/FINAL_mcl_FP_50_allclusters",
  mgclus  = "/scratch/2020_clustering/benchmark/FP/FINAL_mgclus_FP_50_allclusters"
)

for (method in names(cluster_files)) {
  gs     <- read.delim(cluster_files[[method]], header = TRUE, stringsAsFactors = FALSE)
  annx   <- run_anubix_fp(net, paths, gs)
  merged <- merge_enrichment_results(
    annx_res   = annx,
    binox_path = sprintf("/scratch/2020_clustering/benchmark/FP/results/FINAL_binox_FP_cluster_50_%s", method),
    neat_path  = sprintf("/scratch/2020_clustering/benchmark/FP/results/280520FINAL_NEAT_FP_cluster_50_%s.tsv", method),
    gea_path   = sprintf("/scratch/2020_clustering/benchmark/FP/results/FINAL_gea_fp_cluster_50_%s", method)
  )
  out_file <- sprintf("/scratch/2020_clustering/benchmark/FP/results/FINAL_all_merge_cluster_FP_50_%s", method)
  write.table(merged, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
}
