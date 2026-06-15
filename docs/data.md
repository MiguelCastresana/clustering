# Data Notes

This repository includes the benchmark inputs and result tables used by the original analysis.

## Committed Inputs

`input/` contains:

- True-positive and false-positive gene sets.
- KEGG overlap resources.
- Precomputed pathway-recovery inputs for MCL, MGclus, Infomap, and the non-clustered baseline.
- `input_TP_genesets.RData`, which stores the `kegg_kegg` object used by `src/clustering_genesets.R`.

## Committed Results

`results/` contains benchmark outputs for clustered and non-clustered methods. These files are used by merge and plotting scripts.

## External Network

New benchmark runs require a functional association network, for example the FunCoup human compact network. Download current networks from:

<https://funcoup.org/downloads/>

Keep large downloaded networks outside Git unless they are intentionally curated as small test fixtures.

## Bundled Algorithms

`cluster_algorithms/` contains original benchmark assets for MGclus and MCL. They are kept for reproducibility, but future locally generated cluster files are ignored by `.gitignore`.
