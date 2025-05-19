# Clustered pathway enrichment analysis

Gene sets derived from experiments are often heterogeneous, meaning that they represent multiple pathways, see visual example:
<p align="center">
  <img src="Figure.png" width="400">
</p>

A way to counteract this is to cluster the gene set into more homogenous parts before performing pathway analysis on each module. We explored whether network-based pre-clustering of a query gene set can improve pathway analysis:

<p align="center">
  <img src="Figure1.png" width="400">
</p>


The methods MCL, Infomap, and MGclus were used to cluster the gene set projected onto the FunCoup network. We characterized how well these methods are able to detect individual pathways in multi-pathway gene sets, and applied each of the clustering methods in combination with four pathway analysis methods: Gene Enrichment Analysis, BinoX, NEAT, and ANUBIX.

## Overview

The scripts provided here are designed to generate the data for the benchmark and to cluster genes into functional modules using different clustering methods. These methods are applied to pathway data to identify groups of genes that may share functional relationships. The clustering results can be used for further analysis, such as identifying key pathways involved in biological processes.

For a detailed explanation of the advantages of clustering for pathway enrichment analysis tools, check this paper: [Clustered pathway enrichment analysis](https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2022.855766/full).

The method was implemented in the R package ANUBIX under the function anubix_clustering, see the repository [here](https://github.com/MiguelCastresana/anubix) for more info.

Additionally it is also implemented in our Web-tool [PathBIX](https://pathbix.sbc.su.se/).

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

## Scripts Overview

### R Scripts:

**Analysis:**

1. **bisected_path_Creation.R**: Bisects pathways keeping an expected overlap. We do this to create our own sets of ground truth gene sets.

2. **set_creation.R**: Creates TP and FP gene sets.

3. **clustering_genesets.R**: Clusters gene sets using MGclus, MCL and Infomap

4. **cluster_jaccard_comparison.R**: Benchmarking clustering methods for pathway recovery using Jaccard similarity between known pathways and predicted modules

5. **Cluster_Degree_Analysis.R**: Compare network degree distributions between real and random gene sets to evaluate bias in false positives and clustering methods.

6. **merge_FP_analysis.R**: Merge FP results

7. **merge_TP_analysis.R**: Merge TP results

8. **roc_curves.R**: Generate the ROC curves combining the clustering and non clustering results




**Contact**
Miguel Castresana Aguirre (miguel.castresana.aguirre@ki.se)