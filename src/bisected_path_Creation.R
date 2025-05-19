
# Benchmark overlap and split KEGG pathways on FunCoup network
# Load libraries

library(dplyr)
library(igraph)

# --- 1. Load and filter network ---
net <- read.delim("/scratch2/pathbix_website/FC5.0_H.sapiens_compact", header = TRUE)
net_filt <- net %>%
  filter(V1 >= 0.8) %>%
  select(gene1 = V3, gene2 = V4, weight = V1)

# --- 2. Load MSigDB and KEGG ---
msigdb <- read.delim("/scratch/NAR_june_ANUBIX/msigdb/msigdb_v7", header = TRUE)
KEGG_raw <- read.delim(
  "/scratch/PathBIX_all/pathbix/data/Rdata/Pathway_databases/ReactomePathways/REACTOME_h_sapiens",
  header = TRUE
)
paths <- unique(KEGG_raw$V2)

# --- 3. Compute Jaccard overlap metrics ---
# Mean overlap between MSigDB and each KEGG pathway
msig_overlap <- function(msig, kegg_list) {
  scores <- numeric()
  groups <- unique(msig$V2)
  for (grp in groups) {
    gs <- msig$V1[msig$V2 == grp]
    for (path in paths) {
      kp <- kegg_list[[path]]
      l <- min(length(gs), length(kp))
      scores <- c(scores, length(intersect(gs, kp)) / l)
    }
  }
  list(mean = mean(scores), median = median(scores))
}

kegg_list <- split(KEGG_raw$V1, KEGG_raw$V2)
overlap_stats <- msig_overlap(msigdb, kegg_list)
overlap_stats

# --- 4. Determine per-pathway overlap counts ---
g <- table(KEGG_raw$V2)
avg_rate <- overlap_stats$mean
# number of overlapping genes per pathway
ov_genes <- round(avg_rate * as.numeric(g))

# --- 5. Split KEGG pathways into two balanced halves with overlap ---
# Compute node degrees
deg_df <- table(c(net_filt$gene1, net_filt$gene2))

split_pathway <- function(genes, degrees) {
  ord <- order(degrees[genes], decreasing = TRUE)
  gs <- genes[ord]
  groupA <- gs[seq(1, length(gs), by = 2)]
  groupB <- setdiff(genes, groupA)
  list(A = groupA, B = groupB)
}

# First bisection
dfA <- data.frame(); dfB <- data.frame()
for (i in seq_along(paths)) {
  path <- paths[i]
  genes <- KEGG_raw$V1[KEGG_raw$V2 == path]
  sp <- split_pathway(genes, deg_df)
  dfA <- bind_rows(dfA, data.frame(gene = sp$A, pathway = path))
  dfB <- bind_rows(dfB, data.frame(gene = sp$B, pathway = path))
}

# Add overlapping genes back to one half randomly
finalA <- data.frame(); finalB <- data.frame()
for (i in seq_along(paths)) {
  path <- paths[i]
  a <- dfA$gene[dfA$pathway == path]
  b <- dfB$gene[dfB$pathway == path]
  n_ov <- ov_genes[i]
  if (n_ov > 0) {
    if (sample(c(TRUE, FALSE), 1)) {
      add <- sample(b, n_ov)
      a2 <- c(a, add); b2 <- b
    } else {
      add <- sample(a, n_ov)
      b2 <- c(b, add); a2 <- a
    }
  } else {
    a2 <- a; b2 <- b
  }
  finalA <- bind_rows(finalA, data.frame(gene = a2, pathway = path))
  finalB <- bind_rows(finalB, data.frame(gene = b2, pathway = path))
}

# --- 6. Save split sets ---
write.table(finalA, 
            file = "/scratch2/Clustering_2022/ReactomeA_overlap",
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(finalB, 
            file = "/scratch2/Clustering_2022/ReactomeB_overlap",
            sep = "\t", row.names = FALSE, quote = FALSE)


