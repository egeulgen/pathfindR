
<!-- README.md is generated from README.Rmd. Please edit that file -->

# <img src="https://github.com/egeulgen/pathfindR/blob/master/inst/extdata/logo.png?raw=true" align="left" height=150/> pathfindR: An R Package for Enrichment Analysis Utilizing Active Subnetworks

<!-- badges: start -->

[![Travis-CI Build
Status](https://travis-ci.org/egeulgen/pathfindR.svg?branch=master)](https://travis-ci.org/egeulgen/pathfindR)
[![Codecov test
coverage](https://codecov.io/gh/egeulgen/pathfindR/branch/master/graph/badge.svg)](https://codecov.io/gh/egeulgen/pathfindR)
[![CRAN
version](http://www.r-pkg.org/badges/version-ago/pathfindR)](https://cran.r-project.org/package=pathfindR)
[![CRAN total
downloads](https://cranlogs.r-pkg.org/badges/grand-total/pathfindR)](https://cran.r-project.org/package=pathfindR)
[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

# Overview

`pathfindR` is a tool for enrichment analysis via active subnetworks.
The package also offers functionalities to cluster the enriched terms
and identify representative terms in each cluster, to score the enriched
terms per sample and to visualize analysis results.

The functionalities of pathfindR is described in detail in *Ulgen E,
Ozisik O, Sezerman OU. 2019. pathfindR: An R Package for Comprehensive
Identification of Enriched Pathways in Omics Data Through Active
Subnetworks. Front. Genet. <https://doi.org/10.3389/fgene.2019.00858>*

For detailed documentation, see [pathfindR’s
website](https://egeulgen.github.io/pathfindR/) or [the pathfindR
wiki](https://github.com/egeulgen/pathfindR/wiki).

# Installation

You can install the released version of pathfindR from CRAN via:

``` r
install.packages("pathfindR")
```

or via [pak](https://pak.r-lib.org/) (this might be preferable given
`pathfindR`’s Bioconductor dependencies):

``` r
install.packages("pak") # if you have not installed "pak"
pak::pkg_install("pathfindR")
```

And the development version from GitHub via
[pak](https://pak.r-lib.org/):

``` r
install.packages("pak") # if you have not installed "pak"
pak::pkg_install("egeulgen/pathfindR")
```

> **IMPORTANT NOTE** For the active subnetwork search component to work,
> the user must have [Java (\>= 8.0)](https://www.java.com/en/)
> installed and path/to/java must be in the PATH environment variable.

We also have docker images available on [Docker
Hub](https://hub.docker.com/repository/docker/egeulgen/pathfindr) and on
[GitHub](https://github.com/egeulgen/pathfindR/packages):

``` bash
# pull image for latest release
docker pull egeulgen/pathfindr:latest

# pull image for specific version (e.g. 1.4.1)
docker pull egeulgen/pathfindr:1.4.1

# pull image for latest development version
docker pull egeulgen/pathfindr:dev
```

# Enrichment Analysis with pathfindR

![pathfindR Enrichment
Workflow](https://github.com/egeulgen/pathfindR/blob/master/vignettes/pathfindr.png?raw=true
"pathfindr Enrichment Workflow")

This workflow takes in a data frame consisting of “gene symbols”,
“change values” (optional) and “associated p values”:

| Gene\_symbol | logFC  | FDR\_p  |
| :----------- | :----: | :-----: |
| FAM110A      | \-0.69 | 3.4e-06 |
| RNASE2       |  1.35  | 1.0e-05 |
| S100A8       |  1.54  | 3.5e-05 |
| S100A9       |  1.03  | 2.3e-04 |

After input testing, any gene symbol that is not in the chosen
protein-protein interaction network (PIN) is converted to an alias
symbol if there is an alias that is in the PIN. After mapping the input
genes with the associated p values onto the PIN, active subnetwork
search is performed. The resulting active subnetworks are then filtered
based on their scores and the number of significant genes they contain.

> An active subnetwork can be defined as a group of interconnected genes
> in a protein-protein interaction network (PIN) that predominantly
> consists of significantly altered genes. In other words, active
> subnetworks define distinct disease-associated sets of interacting
> genes, whether discovered through the original analysis or discovered
> because of being in interaction with a significant gene.

These filtered list of active subnetworks are then used for enrichment
analyses, i.e. using the genes in each of the active subnetworks, the
significantly enriched terms (pathways/gene sets) are identified.
Enriched terms with adjusted p values larger than the given threshold
are discarded and the lowest adjusted p value (over all active
subnetworks) for each term is kept. This process of `active subnetwork
search + enrichment analyses` is repeated for a selected number of
iterations, performed in parallel. Over all iterations, the lowest and
the highest adjusted-p values, as well as number of occurrences over all
iterations are reported for each significantly enriched term.

This workflow can be run using the function `run_pathfindR()`:

``` r
library(pathfindR)
output_df <- run_pathfindR(input_df)
```

This wrapper function performs the active-subnetwork-oriented enrichment
analysis and returns a data frame of enriched terms (as well as
visualization of enriched terms and an HTML report):

![pathfindR Enrichment
Chart](https://github.com/egeulgen/pathfindR/blob/master/vignettes/enrichment_chart.png?raw=true
"Enrichment Chart")

Some useful arguments are:

``` r
# change the output directory
output_df <- run_pathfindR(input_df, output_dir = "/top/secret/results")

# change the gene sets used for analysis (default = "KEGG")
output_df <- run_pathfindR(input_df, gene_sets = "GO-MF")

# change the PIN for active subnetwork search (default = Biogrid)
output_df <- run_pathfindR(input_df, pin_name_path = "IntAct")
# or use an external PIN of your choice
output_df <- run_pathfindR(input_df, pin_name_path = "/path/to/myPIN.sif")

# change the number of iterations (default = 10)
output_df <- run_pathfindR(input_df, iterations = 25) 

# report the non-significant active subnetwork genes (for later analyses)
output_df <- run_pathfindR(input_df, list_active_snw_genes = TRUE)
```

The available PINs are “Biogrid”, “STRING”, “GeneMania”, “IntAct”,
“KEGG” and “mmu\_STRING”. The available gene sets are “KEGG”,
“Reactome”, “BioCarta”, “GO-All”, “GO-BP”, “GO-CC”, “GO-MF”, and
“mmu\_KEGG”. You also use a custom PIN (see `?return_pin_path`) or a
custom gene set (see `?fetch_gene_set`)

> As of the latest development version, pathfindR offers utility
> functions for obtaining organism-specific PIN data (for now, only
> BioGRID PINs) and organism-specific gene sets (KEGG and Reactome) data
> via `get_pin_file()` and `get_gene_sets_list()`, respectively.

# Clustering of the Enriched Terms

![Enriched Terms Clustering
Workflow](https://github.com/egeulgen/pathfindR/blob/master/vignettes/term_clustering.png?raw=true
"Enriched Terms Clustering Workflow") The wrapper function for this
workflow is `cluster_enriched_terms()`.

This workflow first calculates the pairwise kappa statistics between the
enriched terms. The function then performs hierarchical clustering (by
default), automatically determines the optimal number of clusters by
maximizing the average silhouette width and returns a data frame with
cluster assignments.

``` r
# default settings
clustered_df <- cluster_enriched_terms(output_df)

# display the heatmap of hierarchical clustering
clustered_df <- cluster_enriched_terms(output_df, plot_hmap = TRUE)

# display the dendrogram and automatically-determined clusters
clustered_df <- cluster_enriched_terms(output_df, plot_dend = TRUE)

# change agglomeration method (default = "average") for hierarchical clustering
clustered_df <- cluster_enriched_terms(output_df, clu_method = "centroid")
```

Alternatively, the `fuzzy` clustering method (as described in Huang DW,
Sherman BT, Tan Q, et al. The DAVID Gene Functional Classification Tool:
a novel biological module-centric algorithm to functionally analyze
large gene lists. Genome Biol. 2007;8(9):R183.) can be used:

``` r
clustered_df_fuzzy <- cluster_enriched_terms(output_df, method = "fuzzy")
```

# Visualization of Enrichment Results

## Term-Gene Heatmap

The function `term_gene_heatmap()` can be utilized to visualize the
heatmap of enriched terms by the involved input genes. This heatmap
allows visual identification of the input genes involved in the enriched
terms, as well as the common or distinct genes between different terms.
If the input data frame (same as in `run_pathfindR()`) is supplied, the
tile colors indicate the change values.

![Term-Gene
Heatmap](https://github.com/egeulgen/pathfindR/blob/master/vignettes/hmap.png?raw=true
"Term-Gene Heatmap")

## Term-Gene Graph

The function `term_gene_graph()` (adapted from the Gene-Concept network
visualization by the R package `enrichplot`) can be utilized to
visualize which significant genes are involved in the enriched terms.
The function creates the term-gene graph, displaying the connections
between genes and biological terms (enriched pathways or gene sets).
This allows for the investigation of multiple terms to which significant
genes are related. The graph also enables determination of the degree of
overlap between the enriched terms by identifying shared and/or distinct
significant genes.

![Term-Gene
Graph](https://github.com/egeulgen/pathfindR/blob/master/vignettes/term_gene.png?raw=true
"Term-Gene Graph")

## UpSet Plot

UpSet plots are plots of the intersections of sets as a matrix. This
function creates a ggplot object of an UpSet plot where the x-axis is
the UpSet plot of intersections of enriched terms. By default (i.e.,
`method = "heatmap"`), the main plot is a heatmap of genes at the
corresponding intersections, colored by up/down regulation (if
`genes_df` is provided, colored by change values). If `method =
"barplot"`, the main plot is bar plots of the number of genes at the
corresponding intersections. Finally, if `method = "boxplot"` and
`genes_df` is provided, then the main plot displays the boxplots of
change values of the genes at the corresponding intersections.

![UpSet
plot](https://github.com/egeulgen/pathfindR/blob/master/vignettes/upset.png?raw=true
"UpSet Plot")

# Per Sample Enriched Term Scores

![Agglomerated Scores for all Enriched Terms per
Sample](https://github.com/egeulgen/pathfindR/blob/master/vignettes/score_hmap.png?raw=true
"Scoring per Sample")

The function `score_terms()` can be used to calculate the agglomerated z
score of each enriched term per sample. This allows the user to
individually examine the scores and infer how a term is overall altered
(activated or repressed) in a given sample or a group of samples.

# Comparison of 2 pathfindR Results

The function `combine_pathfindR_results()` allows combination of two
pathfindR active-subnetwork-oriented enrichment analysis results for
investigating common and distinct terms between the groups. Below is an
example for comparing two different results using rheumatoid
arthritis-related data.

``` r
combined_df <- combine_pathfindR_results(result_A = RA_output, 
                                         result_B = RA_comparison_output)
```

By default, `combine_pathfindR_results()` plots the term-gene graph for
the common terms in the combined results. The function
`combined_results_graph()` can be used to create this graph (using only
selected terms etc.) later on.

``` r
combined_results_graph(combined_df, selected_terms = c("hsa04144", "hsa04141", "hsa04140"))
```

![Combined Results
Graph](https://github.com/egeulgen/pathfindR/blob/master/vignettes/combined_graph.png?raw=true
"Combined Results Graph")
