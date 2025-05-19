
# Load required libraries
library(dplyr)

# 1. Load a geneset file (FP or cluster output)
load_genesets <- function(path, header = TRUE) {
  df <- read.delim(path, header = header, stringsAsFactors = FALSE)
  names(df)[1:2] <- c("gene", "geneset")
  return(df)
}

# 2. Compute the size of each geneset
compute_geneset_sizes <- function(geneset_df) {
  geneset_df %>%
    count(geneset, name = "size") %>%
    arrange(geneset)
}

# 3. Load network and build degree list
load_degree_list <- function(net_path, score_col = 3, cutoff = 0.75) {
  net <- read.delim(net_path, header = TRUE, stringsAsFactors = FALSE)
  colnames(net)[1:3] <- c("from", "to", "score")
  net_filt <- net %>% filter(score >= cutoff)
  degree_list <- data.frame(
    gene = c(net_filt$from, net_filt$to),
    stringsAsFactors = FALSE
  ) %>%
    count(gene, name = "degree")
  return(degree_list)
}

# 4. Load MSigDB frequency/probability table
load_msigdb_probs <- function(msigdb_path) {
  msig <- read.delim(msigdb_path, header = TRUE, stringsAsFactors = FALSE)
  freq_tbl <- msig %>%
    count(gene) %>%
    mutate(prob = n / sum(n)) %>%
    select(gene, prob)
  return(freq_tbl)
}

# 5. Generate random FP genesets matching given sizes
sample_random_fp <- function(size_tbl, msig_prob, n_iter = 100) {
  replicate(n_iter, {
    lapply(size_tbl$size, function(sz) {
      sample(msig_prob$gene, sz, prob = msig_prob$prob, replace = FALSE)
    })
  }, simplify = FALSE)
}

# 6. Compute mean and median degree for each geneset
compute_degree_stats <- function(geneset_list, degree_list) {
  stats <- lapply(geneset_list, function(genes) {
    degrees <- degree_list %>% filter(gene %in% genes) %>% pull(degree)
    degrees[is.na(degrees)] <- 0
    c(mean = mean(degrees), median = median(degrees))
  })
  stats_df <- do.call(rbind, stats)
  rownames(stats_df) <- names(geneset_list)
  as.data.frame(stats_df)
}

# 7. Summarize degree stats across iterations
summarize_iteration_stats <- function(stats_list) {
  df <- bind_rows(lapply(stats_list, as.data.frame), .id = "iter")
  df %>%
    group_by(iter) %>%
    summarize(
      mean_of_means     = mean(mean),
      median_of_medians = median(median)
    )
}

# 8. Analyze existing FP or cluster-based genesets
analyze_existing_sets <- function(geneset_df, degree_list) {
  genes_by_set <- split(geneset_df$gene, geneset_df$geneset)
  stats <- compute_degree_stats(genes_by_set, degree_list)
  summary_stats <- stats %>%
    summarize(
      overall_mean   = mean(mean),
      overall_median = median(median)
    )
  list(per_set = stats, summary = summary_stats)
}

# ===== Example Usage =====
# Define file paths
fp_path       <- "/scratch/2020_clustering/benchmark/FP/FINAL_FP_50"
msigdb_path   <- "/scratch/2020_clustering/benchmark/msigdb_v7"
network_path  <- "/scratch/2020_clustering/benchmark/fc4.1"
cluster_paths <- list(
  infomap = "/scratch/2020_clustering/benchmark/FP/FINAL_infomap_FP_50_allclusters",
  mcl     = "/scratch/2020_clustering/benchmark/FP/FINAL_mcl_FP_50_allclusters",
  mgclus  = "/scratch/2020_clustering/benchmark/FP/FINAL_mgclus_FP_50_allclusters"
)

# 1. Load inputs
fp_df     <- load_genesets(fp_path)
deg_list  <- load_degree_list(network_path)
msig_prob <- load_msigdb_probs(msigdb_path)

# 2. Analyze original FP sets
fp_analysis <- analyze_existing_sets(fp_df, deg_list)
print("Original FP sets degree summary:")
print(fp_analysis$summary)

# 3. Generate and summarize random FP iterations
table_sizes    <- compute_geneset_sizes(fp_df)
random_fps      <- sample_random_fp(table_sizes, msig_prob, n_iter = 100)
iter_stats      <- lapply(random_fps, compute_degree_stats, degree_list = deg_list)
summary_df      <- summarize_iteration_stats(iter_stats)
random_medians  <- summary_df$median_of_medians
print("Random FP iterations summary (medians):")
print(summary_df)

# 4. Compare cluster methods vs random using Wilcoxon rank-sum test
cluster_results <- lapply(names(cluster_paths), function(name) {
  df            <- load_genesets(cluster_paths[[name]])
  analysis      <- analyze_existing_sets(df, deg_list)
  cluster_meds  <- analysis$per_set$median
  test_rr       <- wilcox.test(cluster_meds, random_medians)
  list(
    method           = name,
    summary_stats    = analysis$summary,
    wilcox_statistic = unname(test_rr$statistic),
    p_value          = test_rr$p.value
  )
})
print("Cluster vs Random Wilcoxon results:")
print(cluster_results)

# 5. Pairwise comparisons among cluster methods
cluster_medians <- lapply(names(cluster_paths), function(name) {
  df <- load_genesets(cluster_paths[[name]])
  analyze_existing_sets(df, deg_list)$per_set$median
})
names(cluster_medians) <- names(cluster_paths)
pairwise_results <- combn(names(cluster_medians), 2, function(pair) {
  x    <- cluster_medians[[pair[1]]]
  y    <- cluster_medians[[pair[2]]]
  test <- wilcox.test(x, y)
  data.frame(
    comparison       = paste(pair, collapse = " vs "),
    wilcox_statistic = unname(test$statistic),
    p_value          = test$p.value,
    stringsAsFactors = FALSE
  )
}, simplify = FALSE) %>% bind_rows()
print("Pairwise cluster comparisons (Wilcoxon):")
print(pairwise_results)
