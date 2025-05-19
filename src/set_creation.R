

# Load libraries
library(dplyr)
library(fastmatch)
library(parallel)

# 1. Load the functional network and compute degree list
load_network <- function(net_path, score_col = "score", cutoff = 0.75) {
  net <- read.delim(net_path, header = TRUE, stringsAsFactors = FALSE)
  net <- net %>% filter(.data[[score_col]] >= cutoff)
  # net columns: from, to, score
  degree_list <- data.frame(
    gene = c(net[[3]], net[[4]]),
    stringsAsFactors = FALSE
  ) %>%
    count(gene, name = "degree")
  return(list(net = net, degree_list = degree_list))
}

# 2. Load KEGG (or Reactome) overlaps
load_pathways <- function(path_file, header = FALSE, skip = 1) {
  kegg <- read.delim(path_file, header = header, skip = skip,
                     stringsAsFactors = FALSE)
  names(kegg)[1:2] <- c("gene", "pathway")
  return(kegg)
}

# 3. Compute pathway sizes and split into bins of ~45
compute_path_bins <- function(kegg_df, bin_size = 45) {
  sizes <- kegg_df %>%
    group_by(pathway) %>%
    summarise(size = n(), .groups = 'drop') %>%
    arrange(size)
  bins <- split(sizes$pathway, ceiling(seq_along(sizes$size) / bin_size))
  return(bins)
}

# 4. Sample one pathway from each bin to form a geneset
sample_genesets <- function(path_bins, num_sets = 100) {
  lapply(seq_len(num_sets), function(i) {
    sapply(path_bins, function(bin) sample(bin, 1), USE.NAMES = FALSE)
  })
}

# 5. Annotate each gene–pathway pair with network degree
annotate_degrees <- function(kegg_df, degree_list) {
  kegg_df %>%
    left_join(degree_list, by = c("gene" = "gene")) %>%
    mutate(degree = ifelse(is.na(degree), 0, degree))
}

# 6. Build long DataFrame of geneset membership
build_geneset_df <- function(annotated_kegg, genesets) {
  df_list <- lapply(seq_along(genesets), function(i) {
    paths <- genesets[[i]]
    subset <- annotated_kegg %>% filter(pathway %in% paths)
    subset$geneset <- paste0("GENESET_", i, "_#", paste(paths, collapse = "/"))
    subset[, c("gene", "pathway", "geneset", "degree")]
  })
  bind_rows(df_list)
}

# 7. Generate false-positive (noise) half by degree matching
generate_false_half <- function(geneset_df, full_kegg_df, degree_list) {
  all_genes <- unique(full_kegg_df$gene)
  output <- lapply(seq_len(nrow(geneset_df)), function(i) {
    target_deg <- geneset_df$degree[i]
    geneset_id <- geneset_df$geneset[i]
    # exclude genes already in this geneset
    existing <- geneset_df %>% filter(geneset == geneset_id) %>% pull(gene)
    candidates <- degree_list %>%
      filter(!gene %in% existing, degree == target_deg) %>%
      pull(gene)
    if (length(candidates) == 0) {
      sample(setdiff(full_kegg_df$gene, degree_list$gene), 1)  # random from outside
    } else {
      sample(candidates, 1)
    }
  })
  data.frame(
    gene = unlist(output),
    geneset = geneset_df$geneset,
    stringsAsFactors = FALSE
  )
}

# 8. Save KEGG sub-networks for each geneset
save_geneset_networks <- function(net_df, geneset_df, out_file) {
  geneset_nets <- lapply(unique(geneset_df$geneset), function(gs) {
    genes <- unique(geneset_df$gene[geneset_df$geneset == gs])
    sub_net <- net_df %>%
      filter(from %in% genes & to %in% genes) %>%
      select(from, to, score)
    sub_net
  })
  names(geneset_nets) <- unique(geneset_df$geneset)
  save(geneset_nets, file = out_file)
}

# 9. Generate biased false-positive sets from MSigDB
generate_biased_fps <- function(msigdb_file, output_file, num_sets = 100,
                                min_size = 10, max_size = 500) {
  msig <- read.delim(msigdb_file, header = TRUE, stringsAsFactors = FALSE)
  freq_tbl <- msig %>% count(gene) %>%
    mutate(prob = n / sum(n))
  fps <- data.frame()
  for (i in seq_len(num_sets)) {
    size <- sample(min_size:max_size, 1)
    genes <- sample(freq_tbl$gene, size, prob = freq_tbl$prob)
    fps <- bind_rows(fps,
                     data.frame(gene = genes,
                                geneset = paste0("geneset", i),
                                stringsAsFactors = FALSE))
  }
  write.table(fps, output_file, sep = "\t",
              col.names = TRUE, row.names = FALSE, quote = FALSE)
  return(fps)
}

# 10. Correlate network degree with MSigDB occurrence
degree_vs_msigdb <- function(network_file, msigdb_file) {
  net <- read.delim(network_file, header = TRUE, stringsAsFactors = FALSE)
  net <- filter(net, score >= 0.75)
  deg_tbl <- data.frame(gene = c(net[[3]], net[[4]]),
                        stringsAsFactors = FALSE) %>%
    count(gene, name = "degree")
  msig <- read.delim(msigdb_file, header = TRUE, stringsAsFactors = FALSE)
  msig_tbl <- msig %>% count(gene, name = "times")
  res <- deg_tbl %>%
    inner_join(msig_tbl, by = "gene") %>%
    summarize(correlation = cor(degree, times))
  return(res$correlation)
}

# ===== Example Execution =====
# user_paths <- "/scratch/2020_clustering/KEGGB_overlap_2020"
# user_net   <- "/scratch/2020_clustering/fc4.1"
# msig_file  <- "/scratch/2020_clustering/msigdb/msigdb_v7"

# 1. Load data
net_data <- load_network(user_net)
kegg_df <- load_pathways(user_paths)

# 2. Annotate with degrees
kegg_annot <- annotate_degrees(kegg_df, net_data$degree_list)

# 3. Create and save true genesets
path_bins  <- compute_path_bins(kegg_df)
genesets   <- sample_genesets(path_bins)
geneset_df <- build_geneset_df(kegg_annot, genesets)

# 4. Generate false-half noise and combine
false_df   <- generate_false_half(geneset_df, kegg_df, net_data$degree_list)
combined_df <- bind_rows(
  geneset_df %>% select(gene, geneset),
  false_df
)
write.table(combined_df, "/scratch/2020_clustering/TP/TP_genesets",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# 5. Save sub-networks
save_geneset_networks(net_data$net %>% rename(from = V3, to = V4, score = V1),
                      combined_df, "/scratch/2020_clustering/benchmark/TP/TP_genesets.RData")

# 6. Generate biased FPs
biased_fps <- generate_biased_fps(msig_file,
                                  "/scratch/2020_clustering/benchmark/FP/FP_genesets")

# 7. Correlation analysis
corr_value <- degree_vs_msigdb(user_net, msig_file)
print(paste("Degree–MSigDB correlation: ", round(corr_value, 4)))
