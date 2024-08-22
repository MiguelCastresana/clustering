# Clustered pathway enrichment analysis

This repository contains scripts for clustering genes based on functional pathway data using various algorithms, including MGclus, MCL, and Infomap.

## Overview

The scripts provided here are designed to generate the data for the benchmark and to cluster genes into functional modules using different clustering methods. These methods are applied to pathway data to identify groups of genes that may share functional relationships. The clustering results can be used for further analysis, such as identifying key pathways involved in biological processes.

For a detailed explanation of the advantages of clustering for pathway enrichment analysis tools, check this paper: [Clustered pathway enrichment analysis](https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2022.855766/full).

## Repository Contents

This repository provides scripts for the following clustering methods:

- **MGclus**: A method for clustering employing shared neighbors
- **MCL (Markov Cluster Algorithm)**: A popular algorithm for clustering data based on flow simulations.
- **Infomap**: A method that uses information theory to cluster genes into modules.


