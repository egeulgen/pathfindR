## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
suppressPackageStartupMessages(library(pathfindR))
data("RA_input")
knitr::kable(head(RA_input))

## ----eval=FALSE----------------------------------------------------------
#  RA_output <- run_pathfindR(RA_input)

## ----eval=FALSE----------------------------------------------------------
#  # to change the PIN (default = Biogrid)
#  RA_output <- run_pathfindR(RA_input, pin_name = "IntAct")
#  # to use an external PIN of user's choice
#  RA_output <- run_pathfindR(RA_input, pin_name = "/path/to/myPIN.sif")
#  
#  # available gene sets are KEGG, Reactome, BioCarta, GO-BP, GO-CC and GO-MF
#  # default is KEGG
#  # to change the gene sets used for enrichment analysis
#  RA_output <- run_pathfindR(RA_input, gene_sets = "BioCarta")
#  
#  # to change the active subnetwork search algorithm (default = "GR", i.e. greedy algorithm)
#  # for simulated annealing:
#  RA_output <- run_pathfindR(RA_input, search_method = "SA")
#  
#  # to change the number of iterations (default = 10)
#  RA_output <- run_pathfindR(RA_input, iterations = 5)
#  
#  # to manually specify the number processes used during parallel loop by foreach
#  # defaults to the number of detected cores
#  RA_output <- run_pathfindR(RA_input, n_processes = 2)

## ------------------------------------------------------------------------
data("RA_output")
knitr::kable(head(RA_output, 2))

## ---- fig.height=4, fig.width=8------------------------------------------
data("RA_output")
RA_clustered <- choose_clusters(RA_output)
## First 2 rows of clustered pathways data frame
knitr::kable(head(RA_clustered, 2))
## The 16 representative pathways
knitr::kable(RA_clustered[RA_clustered$Status == "Representative", ])

# to display the heatmap of pathway clustering
RA_clustered <- choose_clusters(RA_output, plot_heatmap = TRUE)

# to display the dendrogram and clusters
RA_clustered <- choose_clusters(RA_output, plot_dend = TRUE)

# to change agglomeration method (default = "average")
RA_clustered <- choose_clusters(RA_output, agg_method = "centroid", plot_dend = TRUE)

## ----eval=FALSE----------------------------------------------------------
#  choose_clusters(RA_output, auto = FALSE)

