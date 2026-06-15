# Workflow

Run scripts from the repository root.

## 1. Build Benchmark Gene Sets

```bash
Rscript src/bisected_path_Creation.R
Rscript src/set_creation.R
```

These scripts construct true-positive and false-positive gene sets from pathway resources.

## 2. Cluster Projected Gene Sets

```bash
Rscript src/clustering_genesets.R
```

This runs MGclus, MCL, and Infomap and writes module assignments to `results/clusters/`.

## 3. Compare Cluster Recovery

```bash
Rscript src/cluster_jaccard_comparison.R
Rscript src/Cluster_Degree_Analysis.R
```

These scripts evaluate module/pathway overlap and network-degree behavior.

## 4. Merge Benchmark Outputs

```bash
Rscript src/merge_TP_analysis.R
Rscript src/merge_FP_analysis.R
```

The merge scripts combine pathway-analysis outputs across methods.

## 5. Plot ROC Curves

```bash
Rscript src/roc_curves.R
```

The plotting script reads from `results/` by default and compares clustered against non-clustered workflows.

## Notes

- Full benchmark reruns require the original network/pathway resources and external tools.
- `tools/check-project.R` is a quick repository-health check, not a scientific validation run.
