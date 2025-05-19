# Load required libraries
library(ggplot2)
library(knitr)
library(data.table)
library(reshape2)
library(sets)
library(ggpubr)

# Define Jaccard similarity function
jaccard <- function(x, y) {
  gset_similarity(as.set(x), as.set(y), "Jaccard")
}

# Define function to calculate maximum similarity
mediaf <- function(x) {
  max(x)
}

# Prepare true data and cluster modules
prepare_data <- function(true_path, module_path) {
  trueclus <- read.delim(true_path, header = TRUE)
  groups <- unique(as.vector(trueclus[, 3]))
  
  trueclus[, 4] <- unlist(lapply(groups, function(g) rep(which(groups == g), sum(trueclus[, 3] == g))))
  trueclus <- trueclus[, c(1, 4, 2)]
  
  modules <- read.delim(module_path, header = TRUE)
  sub <- gsub('(.*)#', '', as.vector(modules[, 2]))
  unique_groups <- unique(sub)
  
  modules[, 3] <- unlist(lapply(unique_groups, function(g) rep(which(unique_groups == g), sum(sub == g))))
  
  names(trueclus) <- c("gene", "geneset", "path")
  names(modules) <- c("gene", "module", "geneset")
  
  merged <- merge(trueclus, modules, by = c("gene", "geneset"))
  unique(merged)
}

# Compute Jaccard index per geneset
compute_jaccard_scores <- function(test) {
  media <- numeric()
  for (i in seq_along(unique(test$geneset))) {
    set <- test[test$geneset == i, ]
    path_genes <- by(set, set$path, function(x) x$gene)
    module_genes <- by(set, set$module, function(x) x$gene)
    output <- sapply(path_genes, function(x) sapply(module_genes, function(y) jaccard(x, y)))
    
    maximums <- numeric()
    if (!is.null(dim(output))) {
      df_output <- as.data.frame(output)
      for (j in seq_along(unique(set$path))) {
        max1 <- max(df_output[, j])
        pos <- which(df_output == max1, arr.ind = TRUE)
        max2 <- max(df_output[pos[1], ])
        maximums[j] <- max2
        df_output <- df_output[-pos[1], -pos[2], drop = FALSE]
      }
    }
    media[i] <- mean(maximums)
  }
  media
}

# Process all clustering methods
mcl_test <- prepare_data("100_3paths_paper_pathwayrecovery", "100_3paths_paper_pathwayrecovery_mcl")
mcl <- compute_jaccard_scores(mcl_test)

mgclus_test <- prepare_data("100_3paths_paper_pathwayrecovery", "100_3paths_paper_pathwayrecovery_mgclus")
mgclus <- compute_jaccard_scores(mgclus_test)

infomap_test <- prepare_data("100_3paths_paper_pathwayrecovery", "100_3paths_paper_pathwayrecovery_infomap")
infomap <- compute_jaccard_scores(infomap_test)

# Combine and analyze results
d1 <- c(mgclus, mcl, infomap)
d2 <- rep(c("MGclus", "MCL", "Infomap"), each = length(mgclus))
dat <- data.frame(jaccard_index = as.numeric(d1), Method = d2)

# Statistical tests
res.aov <- aov(jaccard_index ~ Method, data = dat)
summary(res.aov)
TukeyHSD(res.aov)

# Plot
my_comparisons <- list(c("MGclus", "MCL"), c("MCL", "Infomap"), c("MGclus", "Infomap"))
ggboxplot(dat, x = "Method", y = "jaccard_index", 
          order = c("MGclus", "MCL", "Infomap"), add = "jitter",
          ylab = "Jaccard index", xlab = "Method") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 13),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15)) +
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method = "wilcox.test") +
  stat_compare_means(label.y = 1.14, method = "kruskal.test")



