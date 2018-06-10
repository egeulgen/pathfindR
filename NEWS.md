# pathfindR 1.2

## Major Changes
- Implemented the (per subject) pathway scoring function `calculate_pw_scores` and the function to plot the heatmap of pathway scores per subject `plot_scores`.

- Added the `auto` parameter to `choose_clusters`. When `auto == TRUE` (default), the function chooses the optimal number of clusters `k` automatically, as the value which maximizes the average silhouette width. It then returns a data frame with the cluster assignments and the representative/member statuses of each pathway.

- Added the `Fold_Enrichment` column to the resulting data frame of `enrichment`, and as a corollary to the resulting data frame of `run_pathfindR`.

- Added the option `bubble` to plot a bubble chart displaying the enrichment results in `run_pathfindR` using the helper function `enrichment_chart`. To plot the bubble chart set `bubble = TRUE` in `run_pathfindR` or use `enrichment_chart(your_result_df)`. 

## Minor changes and bug fixes
- Add the paramater `silent_option` to `run_pathfindR`. When `silent_option == TRUE` (default), the console outputs during active subnetwork search are printed to a file named "console_out.txt". If `silent_option == FALSE`, the output is printed on the screen. Default was set to `TRUE` because multiple console outputs are simultaneously printed when runnning in parallel.

- Added the `list_active_snw_genes` parameter to `run_pathfindR`. When `list_active_snw_genes == TRUE`, the function adds the column `non_DEG_Active_Snw_Genes`, which reports the non-DEG active subnetwork genes for the active subnetwork which was enriched for the given pathway with the lowest p value.

- Added the data `RA_clustered`, which is the example output of the clustering workflow.

- In the function, `run_pathfindR` added the option to specify the argument `output_dir` which specifies the directory to be created under the current working directory for storing the result HTML files. `output_dir` is "pathfindR_Results" by default.

- `run_pathfindR` now checks whether the output directory (`output_dir`) already exists and if it exists, stops and displays an error message. This was implemented to prevent writing over existing results.

- `genes_table.html` now contains a second table displaying the input gene symbols for which there were no interactions in the PIN.

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
