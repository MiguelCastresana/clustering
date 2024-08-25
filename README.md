# Clustered pathway enrichment analysis

This repository contains scripts for clustering genes based on functional pathway data using various algorithms, including MGclus, MCL, and Infomap.

## Overview

The scripts provided here are designed to generate the data for the benchmark and to cluster genes into functional modules using different clustering methods. These methods are applied to pathway data to identify groups of genes that may share functional relationships. The clustering results can be used for further analysis, such as identifying key pathways involved in biological processes.

For a detailed explanation of the advantages of clustering for pathway enrichment analysis tools, check this paper: [Clustered pathway enrichment analysis](https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2022.855766/full).

The method was implemented in the R package ANUBIX under the function anubix_clustering, see the repository [here](https://github.com/MiguelCastresana/anubix) for more info.

## Repository Contents

This repository provides scripts for the following clustering methods:

- **MGclus**: A method for clustering employing shared neighbors
- **MCL (Markov Cluster Algorithm)**: A popular algorithm for clustering data based on flow simulations.
- **Infomap**: A method that uses information theory to cluster genes into modules.


## Requirements

The following R packages are required to run the scripts:

- **readr**
- **dplyr**
- **purrr**
- **tibble**
- **igraph**

- A functional association network is required. For instance, FunCoup human network, download [here](https://funcoup.org/downloads/download.action?type=network&instanceID=24480085&fileName=FC5.0_H.sapiens_compact.gz)

### 1. `data_generation_clustering.R`

**Purpose:**  
This script is responsible for generating data used for the clustering benchmark. It includes functions to load and process KEGG pathway database data and to filter network data.


### 2. `clustering.R`

**Purpose:**  
This script performs the clustering of the input genesets with MGclus, MCL, and Infomap.

**Inputs:**  
- `input_TP_genesets.RData`: This file contains the gene set data required for the analysis.




**Contact**
Miguel Castresana Aguirre (miguel.castresana.aguirre@ki.se)