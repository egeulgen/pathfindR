# <img src="https://static.wixstatic.com/media/69bf8f_9dd2c46cfa8041a69c6be5a18cd74af8~mv2_d_4218_3060_s_4_2.png/v1/fill/w_537,h_280,al_c,usm_0.66_1.00_0.01/69bf8f_9dd2c46cfa8041a69c6be5a18cd74af8~mv2_d_4218_3060_s_4_2.png" align="left" height=120/> pathfindR : An R Package for Pathway Enrichment Analysis Utilizing Active Subnetworks

[![Travis-CI Build Status](https://travis-ci.org/egeulgen/pathfindR.svg?branch=master)](https://travis-ci.org/egeulgen/pathfindR) [![CRAN version](http://www.r-pkg.org/badges/version-ago/pathfindR)](https://cran.r-project.org/package=pathfindR) [![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/pathfindR)](https://cran.r-project.org/package=pathfindR) [![Rdoc](http://www.rdocumentation.org/badges/version/pathfindR)](http://www.rdocumentation.org/packages/pathfindR)

## Overview

`pathfindR` is a tool for pathway enrichment analysis via active subnetworks. The package also offers the option to cluster the enriched pathways and choose representative pathways. The method is described in detail in _Ulgen E, Ozisik O, Sezerman OU. 2018. pathfindR: An R Package for Pathway Enrichment Analysis Utilizing Active Subnetworks. bioRxiv. [https://doi.org/10.1101/272450](https://doi.org/10.1101/272450)_

## Wiki
See [the pathfindR wiki](https://github.com/egeulgen/pathfindR/wiki) for detailed documentation.

## Installation

From CRAN (release):
```r
install.packages("pathfindR")
```

From GitHub (development):
```r
install.packages("devtools") # if you have not installed "devtools" package
devtools::install_github("egeulgen/pathfindR")
```

See the [wiki page](https://github.com/egeulgen/pathfindR/wiki/Installation) for more details.

## Overview of the Enrichment Workflow

![pathfindR Enrichment Workflow](./vignettes/pathfindr.png?raw=true "pathfindr Enrichment Workflow")
This workflow takes in a data frame consisting of Gene Symbol, log-fold-change and adjusted-p values. After input testing, any gene symbol that is not in the protein-protein interaction network (PIN) is converted to an alias symbol if there is an alias that is in the PIN. Next, active subnetwork search is performed. Pathway enrichment analyses are then performed using the genes in each of the identified active subnetworks. Pathways with adjusted p values larger than a given threshold are discarded. The lowest adjusted p value (over all active subnetworks) for each pathway is kept. This process of active subnetwork search and enrichment analyses is repeated for a selected number of iterations, which is done in parallel. Over all iterations, the lowest and the highest adjusted-p values, as well as number of occurrences are reported for each enriched pathway.

This workflow can be run using the function `run_pathfindR`:

```r
RA_output <- run_pathfindR(RA_input)

# to change the output directory
RA_output <- run_pathfindR(RA_input, output = "new_directory")

# to change the PIN (default = Biogrid)
RA_output <- run_pathfindR(RA_input, pin_name = "IntAct")
# to use an external PIN of user's choice
RA_output <- run_pathfindR(RA_input, pin_name = "/path/to/myPIN.sif")

# to change the active subnetwork search algorithm (default = "GR", i.e. greedy algorithm)
# for simulated annealing:
RA_output <- run_pathfindR(RA_input, search_method = "SA")

# to change the number of iterations (default = 10)
RA_output <- run_pathfindR(RA_input, iterations = 5) 

# to manually specify the number processes used during parallel loop by foreach
# defaults to the number of detected cores 
RA_output <- run_pathfindR(RA_input, n_processes = 2)

# to report the non-DEG active subnetwork genes
RA_output <- run_pathfindR(RA_input, list_active_snw_genes = TRUE)
```

See the [wiki page](https://github.com/egeulgen/pathfindR/wiki/Enrichment%20Documentation) for more details.

## Overview of the Clustering Workflow

![Pathway Clustering Workflow](./vignettes/pw_clustering.png?raw=true "Pathway Clustering Workflow")
The wrapper function for this workflow is `choose_clusters()`.

This workflow first calculates the pairwise distances between the pathways in the resulting data frame. By default, the function automatically determines the optimal number of clusters, by maximizing the average silhouette width and returns a data frame with cluster assignments.

```r
# default settings
RA_clustered <- choose_clusters(RA_output)

# to display the heatmap of pathway clustering
RA_clustered <- choose_clusters(RA_output, plot_heatmap = TRUE)

# to display the dendrogram and clusters
RA_clustered <- choose_clusters(RA_output, plot_dend = TRUE)

# to change agglomeration method (default = "average")
RA_clustered <- choose_clusters(RA_output, agg_method = "centroid", plot_dend = TRUE)
```

Alternatively, manual selection of the height at which to cut the dendrogram can be performed. For this, the user should set the `auto` parameter to `FALSE`. Via a shiny app, presented as an HTML document, the hierarchical clustering dendrogram is visualized. In this HTML document, the user can select the agglomeration method and the distance value at which to cut the tree. The dendrogram with the cut-off value marked with a red line is dynamically visualized and the resulting cluster assignments of the pathways along with annotation of representative pathways (chosen by smallest lowest p value) are presented as a table. This table can be saved as a csv file via pressing the button `Get Pathways w\ Cluster Info`. Example usage:

```r
choose_clusters(RA_output, auto = FALSE)
```

See the [wiki page](https://github.com/egeulgen/pathfindR/wiki/Clustering%20Documentation) for more details.

## Overview of the Pathway Scoring Functionality

![Pathway Scoring per Sample](./vignettes/pw_score_hmap.png?raw=true "Pathway Scoring per Sample")
 
The function `calculate_pw_scores` can be used to calculate the pathway scores per sample. This allows the user to individually examine the scores and infer whether a pathway is activated or repressed in a given sample.

For a set of pathways <img src="https://latex.codecogs.com/gif.latex?\inline&space;P&space;=&space;\{P_1,&space;P_2,&space;...&space;,&space;P_k\}" title="P = \{P_1, P_2, ... , P_k\}" />, where each <img src="https://latex.codecogs.com/gif.latex?\inline&space;P_i" title="P_i" /> contains a set of genes, i.e. <img src="https://latex.codecogs.com/gif.latex?\inline&space;P_i&space;=&space;\{g_1,&space;g_2,&space;...\}" title="P_i = \{g_1, g_2, ...\}" />, the pathway score matrix _PS_ is defined as:

<img src="https://latex.codecogs.com/gif.latex?\inline&space;PS_{p,s}&space;=&space;\frac{1}{k}&space;\sum_{g&space;\in&space;P_p}&space;GS_{g,s}" title="PS_{p,s} = \frac{1}{k} \sum_{g \in P_p} GS_{g,s}" /> for each pathway _p_ and for each sample _s_.

_GS_ is the gene score per sample matrix and is defined as:
<img src="https://latex.codecogs.com/gif.latex?\inline&space;GS_{g,s}&space;=&space;(EM_{g,s}&space;-&space;\bar{x}_g)&space;/&space;sd_g" title="GS_{g,s} = (EM_{g,s} - \bar{x}_g) / sd_g" /> where _EM_ is the expression matrix (columns are samples, rows are genes), <img src="https://latex.codecogs.com/gif.latex?\inline&space;\bar{x}_g" title="\bar{x}_g" /> is the mean expression value of the gene and <img src="https://latex.codecogs.com/gif.latex?\inline&space;sd_g" title="sd_g" /> is the standard deviaton of the expression values for the gene.

See the [wiki page](https://github.com/egeulgen/pathfindR/wiki/Pathway-Scoring) for more details.

## Dependencies
For the active subnetwork search component to work, the user must have [JAVA](https://www.java.com/en/download/manual.jsp) installed and path/to/java must be in the PATH environment variable.

## Resources
The PINs were gathered from various resources:
- [Biogrid](https://downloads.thebiogrid.org/BioGRID)
- [GeneMania](http://genemania.org/data/): only interactions with weights >= 0.0006 were kept.
- [IntAct](https://www.ebi.ac.uk/intact/)
- KEGG PIN - created via an in-house script.

The gene sets for enrichment analyses were gathered from:
- KEGG - created using the R package KEGGREST (_Dan Tenenbaum (2018). KEGGREST: Client-side REST access to KEGG. R package version 1.18.1._)
- [Reactome](https://reactome.org/download-data)
- [BioCarta](http://software.broadinstitute.org/gsea/msigdb/genesets.jsp?collection=CP:BIOCARTA)
- [Gene Ontology gene sets](http://www.go2msig.org/cgi-bin/prebuilt.cgi?taxid=9606) - High quality GO annotations only
