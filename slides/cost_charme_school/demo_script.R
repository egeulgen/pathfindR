##################################################
## Project: pathfindR
## Script purpose: COST CHARME Summer Training
## School, Istanbul - pathfindR hands-on demonstration
## Date: Sep 1, 2019
## Author: Ege Ulgen
##################################################

# Installation ------------------------------------------------------------
# For the active subnetwork search component to work(i.e., in order to
# run pathfindR), the user must have Java installed and
# path/to/java must be in the PATH environment variable.
# For Windows users, to configure the PATH environment variable see:
# https://github.com/egeulgen/pathfindR/wiki/Installation#configuration-of-java-on-windows

install.packages("devtools") # if you have not installed "devtools" package
devtools::install_github("egeulgen/pathfindR")

library(pathfindR)

# Enrichment Analysis -----------------------------------------------------
## demo input file = RA_input
?RA_input
dim(RA_input)
head(RA_input)

## demo runs
?run_pathfindR

# takes a while (use `visualize_pathways = FALSE` for faster runs)
RA_demo_out1 <- run_pathfindR(RA_input,
                         iterations = 1) # change number of iter.s to 1 for demo
dim(RA_demo_out1)
head(RA_demo_out1)

# faster non-default run
RA_demo_out2 <- run_pathfindR(RA_input,
                              iterations = 1, # change number of iter.s to 1 for demo
                              gene_sets = "BioCarta", # change from default ("KEGG")
                              pin_name_path = "GeneMania", # change from default ("Biogrid")
                              output = "DEMO_OUTPUT") # change output directory

# Pathway Clustering ------------------------------------------------------
## demo enrichment result file = RA_output
?RA_output
dim(RA_output)
head(RA_output)

?cluster_pathways
RA_demo_clu1 <- cluster_pathways(RA_output) # hierarchical (default)
RA_demo_clu2 <- cluster_pathways(RA_output,
                                 method = "fuzzy")

head(RA_demo_clu1)
head(RA_demo_clu2)

## Plot enrichment chart grouped by clusters
enrichment_chart(RA_demo_clu1,
                 plot_by_cluster = TRUE)

## Example Output for the pathfindR Clustering Workflow
?RA_clustered

# Term-gene graph ---------------------------------------------------------
?term_gene_graph
### `options(stringsAsFactors = TRUE)` if `stringsAsFactors` is set as FALSE in .Rprofile

term_gene_graph(RA_output) # top 10 terms(default)

## Graph using representative pathways
RA_representative <- RA_demo_clu1[RA_demo_clu1$Status == "Representative", ]
term_gene_graph(RA_representative,
                num_terms = NULL, # to plot using all terms
                use_names = TRUE) # use pw names instead of IDs

# Pathway Scoring ---------------------------------------------------------
## Expression matrix = RA_exp_mat
?RA_exp_mat

## Vector of "Case" IDs
cases <- c("GSM389703", "GSM389704", "GSM389706", "GSM389708", "GSM389711",
           "GSM389714", "GSM389716", "GSM389717", "GSM389719", "GSM389721",
           "GSM389722", "GSM389724", "GSM389726", "GSM389727", "GSM389730",
           "GSM389731", "GSM389733", "GSM389735")

?calculate_pw_scores

## Calculate pathway scores and plot heatmap for top 10 enriched pathways
score_matrix <- calculate_pw_scores(RA_output[1:10, ],
                                    RA_exp_mat,
                                    cases)

## Calculate pathway scores and plot heatmap for representative patways
score_matrix <- calculate_pw_scores(RA_representative,
                                    RA_exp_mat,
                                    cases)

# works if cases are not supplied as well
score_matrix <- calculate_pw_scores(RA_representative,
                                    RA_exp_mat)
