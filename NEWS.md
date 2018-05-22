# pathdindR 1.2

## Major Changes
- Added the new function `auto_choose_clusters`, which calculates the pairwise distance values, hierarchically clusters the pathways using these distances and chooses the optimal number of clusters `k` automatically as the value which maximizes the average silhouette width. The function returns a data frame with the cluster assignments and the representative/member statuses of each pathway. `auto_choose_clusters` is an alternative to `choose_clusters`, which needs manual selection of the number of clusters.

- Added the `Fold_Enrichment` column to the resulting data frame of `enrichment`, and as a corollary to the resulting data frame of `run_pathfindR`.

- Added the option to plot a bubble chart displaying the enrichment results in `run_pathfindR`. The x-axis corresponds to fold enrichment values while the y-axis indicates the enriched pathways. Size of the bubble indicates the number of DEGs in the given pathway. Color indicates the -log10(lowest-p) value.

## Minor changes and bug fixes

# pathfindR 1.1

## Major changes
- Added the `gene_sets` option in `run_pathfindR` to chose between different gene sets. Available gene sets are `KEGG`, `Reactome`, `BioCarta` and Gene Ontology gene sets (`GO-BP`, `GO-CC` and `GO-MF`).
- `cluster_pathways` automatically recognizes the ID type and chooses the gene sets accordingly.

## Minor changes and bug fixes
- Fixed issue regarding p values < 1e-13. No active subnetworks were found when there were p values < 1e-13. These are now changed to 1e-13 in the function `input_processing`.
- In `input_processing`, genes for which no interactions are found in the PIN are now removed before active subnetwork search
- Duplicated gene symbols no longer raise an error. If there are duplicated symbols, the lowest p value is chosen for each gene symbol in the function `input_processing`.
- To prevent the formation of nested folders, by default and on errors, the function `run_pathfindR` returns to the user's working directory.
- Citation information are now provided for our BioRxiv pre-print.