# Load necessary libraries
library(tidyverse)
library(parallel)
library(fastmatch)

# Function to load and process KEGGA data
load_and_process_kegga <- function(file_path) {
  KEGGA <- read_delim(file_path, delim = "\t", col_names = FALSE, skip = 1)
  colnames(KEGGA) <- c("Gene", "Pathway")
  KEGGA
}

# Load KEGGA data
KEGGA <- load_and_process_kegga("D:/Clustering_2022/KEGGA_overlap_2020")

# Load and filter network data function
load_and_filter_network <- function(file_path, threshold = 0.8) {
  # Load the network data
  net <- read_delim(file_path, delim = "\t", col_names = TRUE)
  
  # Rename columns to meaningful names if they don't match expectations
  colnames(net) <- c("PFC", "FBS_max", "Gene1", "Gene2")
  
  # Filter based on threshold
  net <- net %>%
    filter(PFC >= threshold) %>%
    select(Gene1, Gene2, PFC)
  
  # Create degree list
  degree_list <- table(c(net$Gene1, net$Gene2)) %>%
    as.data.frame() %>%
    rename(Gene = Var1, Degree = Freq)
  
  list(net = net, degree_list = degree_list)
}

network_data <- load_and_filter_network("D:/Clustering_2022/FC5.0_H.sapiens_compact")
net <- network_data$net
degree_list_total <- network_data$degree_list

# Filter KEGGA data to keep only genes present in the degree list
KEGGA <- KEGGA %>%
  filter(Gene %in% degree_list_total$Gene)

# Split KEGGA by pathways
kegga_split <- split(KEGGA, KEGGA$Pathway)

# Calculate the size of each pathway group
path_sizes <- KEGGA %>%
  group_by(Pathway) %>%
  summarize(Size = n()) %>%
  ungroup() %>%
  arrange(Size)

# Filter and split pathway sizes into groups
filtered_sizes <- path_sizes %>%
  filter(Size > 5) %>%
  pull(Size)

split_sizes <- split(filtered_sizes, ceiling(seq_along(filtered_sizes)/60))

# Create data frames for each size split
path_size_splits <- map(split_sizes, function(size) {
  path_sizes %>%
    filter(Size %in% size)
})

# Function to generate random pathways
generate_random_paths <- function(path_size_splits, n = 50) {
  replicate(n, {
    map_chr(path_size_splits, ~ sample(.x$Pathway, 1))
  }, simplify = FALSE)
}

# Generate group list with random pathways
group_list <- generate_random_paths(path_size_splits)

# Calculate degree values for KEGGA genes
KEGGA <- KEGGA %>%
  rowwise() %>%
  mutate(Degree = degree_list_total$Degree[match(Gene, degree_list_total$Gene)] %>% replace_na(0))

# Generate gene sets from the group list
generate_gene_sets <- function(group_list, KEGGA) {
  map_dfr(seq_along(group_list), function(i) {
    dat <- map_dfr(group_list[[i]], function(group) {
      genes <- KEGGA %>%
        filter(Pathway == group) %>%
        select(Gene) %>%
        mutate(Pathway = group)
      
      genes
    })
    dat <- dat %>%
      mutate(Geneset = paste("geneset", i, sep = ""))
    
    dat
  })
}

dat1 <- generate_gene_sets(group_list, KEGGA)

# Generate names for each geneset
generate_names <- function(group_list) {
  map_chr(seq_along(group_list), function(i) {
    paste("GENESET", "_", i, "_#", paste(group_list[[i]], collapse = "/"), sep = "")
  })
}

nombres <- generate_names(group_list)

# Ensure that each gene set in dat1 is assigned the correct group name
dat1 <- map2_dfr(group_list, nombres, function(groups, name) {
  dat1_subset <- dat1 %>%
    filter(Pathway %in% groups) %>%
    mutate(Group = name)
  return(dat1_subset)
})

# Check the structure to ensure it worked correctly
str(dat1)

# Combine geneset data with names (if necessary, but should already be done)
if (!"Group" %in% colnames(dat1)) {
  dat1 <- map2_dfr(group_list, nombres, function(groups, name) {
    dat1_subset <- dat1 %>%
      filter(Pathway %in% groups) %>%
      mutate(Group = name)
    return(dat1_subset)
  })
}

# Check the structure of dat1 to ensure it's correctly formed
str(dat1)

# If no further modifications are needed, simply proceed with your next steps
no_cluster <- dat1


# Combine geneset data with names
no_cluster <- dat1 %>%
  mutate(Group = rep(nombres, each = nrow(dat1)))

# Clean and save final no_cluster data
no_cluster1 <- no_cluster %>%
  select(Gene, Group) %>%
  distinct()

write_tsv(no_cluster1, "D:/Clustering_2022/6paths_fc5")
