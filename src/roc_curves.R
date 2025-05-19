
library(tidyverse)
library(cowplot)

# ----------------------------------------------------------------------------- 
## CONFIG
base_dir <- "~/FINAL RESULTS"

tp_files <- list(
  cluster    = file.path(base_dir, "INFOMAP_allresults_TP"),
  nocluster  = file.path(base_dir, "noclustering_allresults_TP")
)

fp_files <- list(
  cluster    = file.path(base_dir, "INFOMAP_allresults_FP"),
  nocluster  = file.path(base_dir, "noclustering_allresults_FP")
)

methods <- c("BinoX", "NEAT", "GEA", "ANUBIX")

# ==================== LOAD TP AND FP FILES ====================================
# Load TP geneset data for pathway extraction
tp_pathway_df <- read.delim(tp_files$cluster, header = FALSE, stringsAsFactors = FALSE)
colnames(tp_pathway_df)[1:2] <- c("Geneset", "Pathway")

# Extract all unique HSA pathways from geneset column
tp_pathway_ids <- tp_pathway_df$Geneset |>
  strsplit("__") |>
  sapply(`[`, 2) |>
  strsplit("_") |>
  unlist() |>
  unique() |>
  .[grepl("^HSA\\d{5}$", .)]

# Load FP data
fp_pathway_df <- read.delim(fp_files$cluster, header = FALSE, stringsAsFactors = FALSE)
colnames(fp_pathway_df)[1:2] <- c("Geneset", "Pathway")

# Parameters
alpha_threshold    <- 0.05
n_true_positives   <- length(tp_pathway_ids)  # total unique pathways in TP gene sets
n_true_negatives   <- length(unique(fp_pathway_df$Geneset)) * length(unique(fp_pathway_df$Pathway))

# ==================== GENERIC HELPERS ========================================
read_metrics <- function(path, keep_cols) {
  read_delim(path, delim = "\t", show_col_types = FALSE) %>% 
    select(all_of(keep_cols))
}

melt_pvals <- function(df, cols, tp_flag) {
  map2_dfr(cols, methods, \(col, m) {
    tibble(value = df[[col]],
           tp_flag = tp_flag,
           method  = m) |>
      filter(value < alpha_threshold)
  })
}

make_roc <- function(tp_df, fp_df) {
  bind_rows(tp_df, fp_df) |>
    arrange(value) |>
    group_by(method) |>
    mutate(
      tp = cumsum(tp_flag == 0) / n_true_positives,
      fp = cumsum(tp_flag == 1) / n_true_negatives
    ) |>
    ungroup()
}

process_variant <- function(tag) {
  tp <- read_metrics(tp_files[[tag]], c(1, 2, 5, 12, 15, 17))
  fp <- read_metrics(fp_files[[tag]], c(1, 2, 3, 7, 10, 12))
  
  tp_long <- melt_pvals(tp, c(3:6), tp_flag = 0)
  fp_long <- melt_pvals(fp, c(3:6), tp_flag = 1)
  
  make_roc(tp_long, fp_long) |>
    mutate(clustering = if_else(tag == "cluster", "YES", "NO"))
}

# ============================= RUN PIPELINE ===================================
roc_cluster    <- process_variant("cluster")
roc_nocluster  <- process_variant("nocluster")

roc_all <- bind_rows(roc_cluster, roc_nocluster) %>%
  mutate(
    method     = factor(method, levels = methods),
    clustering = factor(clustering, levels = c("NO", "YES"))
  )

# ============================= PLOT ===========================================
ggplot(roc_all, aes(x = fp, y = tp, color = method)) +
  geom_line(aes(linetype = clustering, size = clustering, alpha = clustering)) +
  scale_color_manual(values = c(
    ANUBIX = "darkred",
    BinoX  = "green3",
    NEAT   = "darkorange1",
    GEA    = "blue"
  )) +
  scale_size_manual(values = c(NO = 0.7, YES = 1.5)) +
  scale_alpha_manual(values = c(NO = 0.5, YES = 1)) +
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 1)) +
  labs(x = "False Positive Rate (FPR)", y = "True Positive Rate (TPR)") +
  theme_minimal(base_family = "Arial") +
  theme(
    legend.position = c(0.55, 0.40),
    axis.text = element_text(size = 14, colour = "black"),
    axis.title = element_text(size = 16),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA)
  )
