#'Wrapper Function for pathfindR Workflow
#'
#'\code{run_pathfindR} is the wrapper function for the pathfindR workflow
#'
#'This function takes in a data frame consisting of Gene Symbol, log-fold-change
#'and adjusted-p values. After input testing, any gene symbols that are not in
#'the PIN are converted to alias symbols if the alias is in the PIN. Next,
#'active subnetwork search is performed. Pathway enrichment analysis is
#'performed using the genes in each of the active subnetworks. Pathways with
#'adjusted-p values lower than \code{enrichment_threshold} are discarded. The
#'lowest adjusted-p value (over all subnetworks) for each pathway is kept. This
#'process of active subnetwork search and enrichment is repeated  for a selected
#'number of \code{iterations}, which is done in parallel. Over all iterations,
#'the lowest and the highest adjusted-p values, as well as number of occurrences
#'are reported for each enriched pathway.
#'
#'@inheritParams input_testing
#'@param enrichment_threshold threshold used when filtering individual pathway
#'  enrichment results
#'@param adj_method correction method to be used for adjusting p-values of
#'  pathway enrichment results (Default: 'bonferroni')
#'@param search_method algorithm to use when performing active subnetwork
#'  search. Options are greedy search (GR), simulated annealing (SA) or genetic
#'  algorithm (GA) for the search (Default:GR. Can be one of c("GR", "SA",
#'  "GA"))
#'@param use_all_positives if TRUE: in GA, adds an individual with all positive
#'  nodes. In SA, initializes candidate solution with all positive nodes.
#'  (Default = FALSE)
#'@param saTemp0 initial temperature for SA (Default: 1.0)
#'@param saTemp1 final temperature for SA (Default: 0.01)
#'@param saIter iteration number for SA (Default: 10000)
#'@param gaPop population size for GA (Default: 400)
#'@param gaIter iteration number for GA (Default: 10000)
#'@param gaThread number of threads to be used in GA (Default: 5)
#'@param gaMut the mutation rate for GA (Default: 0)
#'@param grMaxDepth sets max depth in greedy search. set to 0 for no limit
#'  (Default: 1)
#'@param grSearchDepth sets search depth in greedy search (Default: 1)
#'@param grOverlap sets overlap threshold for results of greedy search (Default:
#'  0.5)
#'@param grSubNum sets number of subnetworks to be presented in the results
#'  (Default: 1000)
#'@param iterations number of iterations for active subnetwork search and
#'  enrichment analyses (Default = 10. Gets set to 1 for Genetic Algorithm)
#'@param n_processes optional argument for specifying the number of processes
#'  used by foreach. If not specified, the function determines this
#'  automatically (Default == NULL. Gets set to 1 for Genetic Algorithm)
#'@inheritParams return_pin_path
#'@param score_thr active subnetwork score threshold (Default = 3)
#'@param sig_gene_thr threshold for minimum number of significant genes (Default = 2)
#'@param gene_sets the gene sets to be used for enrichment analysis. Available gene sets
#'  are KEGG, Reactome, BioCarta, GO-BP, GO-CC, GO-MF or Custom. If "Custom", the arguments
#'  custom_genes and custom pathways must be specified. (Default = "KEGG")
#'@param custom_genes a list containing the genes involved in each custom pathway. Each element
#' is a vector of gene symbols located in the given pathway. Names correspond to
#' the ID of the pathway.
#'@param custom_pathways A list containing the descriptions for each custom pathway. Names of the
#' list correspond to the ID of the pathway.
#'@param bubble boolean value. If TRUE, a bubble chart displaying the enrichment
#' results is plotted. (default = TRUE)
#'@param output_dir the directory to be created under the current working
#' directory where the output and intermediate files are saved (default: "pathfindR_Results")
#'@param list_active_snw_genes boolean value indicating whether or not to report
#' the non-DEG active subnetwork genes for the active subnetwork which was enriched for
#' the given pathway with the lowest p value (default = FALSE)
#'@param silent_option boolean value indicating whether or not to print to the console (FALSE)
#'or print to a file (TRUE) during active subnetwork search (default = TRUE)
#'
#'@return Data frame of pathfindR enrichment results. Columns are: \describe{
#'   \item{ID}{KEGG ID of the enriched pathway}
#'   \item{Pathway}{Description of the enriched pathway}
#'   \item{Fold_Enrichment}{Fold enrichment value for the enriched pathway}
#'   \item{occurrence}{the number of iterations that the given pathway was found to enriched over all iterations}
#'   \item{lowest_p}{the lowest adjusted-p value of the given pathway over all iterations}
#'   \item{highest_p}{the highest adjusted-p value of the given pathway over all iterations}
#'   \item{non_DEG_Active_Snw_Genes (OPTIONAL)}{the non-DEG active subnetwork genes, comma-separated}
#'   \item{Up_regulated}{the up-regulated genes in the input involved in the given pathway, comma-separated}
#'   \item{Down_regulated}{the down-regulated genes in the input involved in the given pathway, comma-separated}
#' }
#'  The function also creates an HTML report with the pathfindR enrichment
#'  results linked to the visualizations of the pathways in addition to
#'  the table of converted gene symbols. This report can be found in
#'  "`output_dir`/results.html" under the current working directory.
#'
#'  Optionally, a bubble chart of enrichment results are plotted. The x-axis
#'  corresponds to fold enrichment values while the y-axis indicates the enriched
#'  pathways. Size of the bubble indicates the number of DEGs in the given pathway.
#'  Color indicates the -log10(lowest-p) value; the more red it gets, the more significant
#'  the pathway is.
#'
#'@import knitr
#'@import rmarkdown
#'@import parallel
#'@import doParallel
#'@import foreach
#'@import graphics
#'
#'@export
#'
#'@section Warning: Depending on the protein interaction network of your choice,
#'  active subnetwork finding component of pathfindR may take a very long time
#'  to finish.
#'
#'@seealso \code{\link{input_testing}} for input testing,
#'  \code{\link{input_processing}} for input processing,
#'  \code{\link{parseActiveSnwSearch}} for parsing an active subnetwork search
#'  output, \code{\link{enrichment}} for pathway enrichment analysis and
#'  \code{\link{pathmap}} for annotation of involved genes and visualization of
#'  pathways. See \code{\link[foreach]{foreach}} for details on parallel
#'  execution of looping constructs. See \code{\link{choose_clusters}} for
#'  clustering the resulting enriched pathways and partitioning into
#'  clusters.
#'
#' @examples
#' \dontrun{
#' run_pathfindR(RA_input)
#' }
run_pathfindR <- function(input, p_val_threshold = 5e-2,
                          enrichment_threshold = 5e-2,
                          adj_method = "bonferroni",
                          search_method = "GR",
                          use_all_positives = FALSE,
                          saTemp0 = 1, saTemp1 = 0.01, saIter = 10000,
                          gaPop = 400, gaIter = 10000, gaThread = 5, gaMut = 0,
                          grMaxDepth = 1, grSearchDepth = 1,
                          grOverlap = 0.5, grSubNum = 1000,
                          iterations = 10, n_processes = NULL,
                          pin_name_path = "Biogrid",
                          score_thr = 3, sig_gene_thr = 2,
                          gene_sets = "KEGG",
                          custom_genes = NULL, custom_pathways = NULL,
                          bubble = TRUE,
                          output_dir = "pathfindR_Results",
                          list_active_snw_genes = FALSE,
                          silent_option = TRUE) {
  ## Argument checks
  if (!search_method %in% c("GR", "SA", "GA"))
    stop("`search_method` must be one of \"GR\", \"SA\", \"GA\"")

  if (!gene_sets %in% c("KEGG", "Reactome", "BioCarta",
                        "GO-BP", "GO-CC", "GO-MF", "Custom"))
    stop("`gene_sets` must be one of KEGG, Reactome, BioCarta, GO-BP, GO-CC, GO-MF or Custom")

  if (gene_sets == "Custom" & (is.null(custom_genes) | is.null(custom_pathways)))
    stop("You must provide both `custom_genes` and `custom_pathways` if `gene_sets` is `Custom`!")

  if (!is.logical(use_all_positives))
    stop("the argument `use_all_positives` must be either TRUE or FALSE")

  if (!is.logical(bubble))
    stop("the argument `bubble` must be either TRUE or FALSE")

  if (!is.logical(silent_option))
    stop("the argument `silent_option` must be either TRUE or FALSE")

  ## create output dir
  if (dir.exists(output_dir)) {
    warning(paste0("There already is a directoy named \"", output_dir,
                   "\". Changing to \"", output_dir, "(1)\" not to overwrite the previous results."))
    output_dir <- paste0(output_dir, "(1)")
  }

  dir.create(output_dir, recursive = TRUE)
  setwd(output_dir)

  ## turn silent_option into an argument
  silent_option <- ifelse(silent_option, " > console_out.txt", "")

  ## If search_method is GA, set iterations as 1
  if (search_method == "GA")
    iterations <- n_processes <- 1

  ## If iterations == 1, set n_processes to 1
  if (iterations == 1)
    n_processes <- 1

  ## turn use_all_positives into the java argument
  use_all_positives <- ifelse(use_all_positives, " -useAllPositives", "")

  ## absolute paths for cytoscape and pin
  active_search_path <- normalizePath(
    system.file("java/ActiveSubnetworkSearch.jar",
                package = "pathfindR"))
  pin_path <- return_pin_path(pin_name_path)

  ## Check input
  message("## Testing input\n\n")
  input_testing(input, p_val_threshold)

  ## Process input
  message("## Processing input. Converting gene symbols, if necessary\n\n")
  input_processed <- input_processing(input, p_val_threshold, pin_path)

  dir.create("active_snw_search")
  utils::write.table(input_processed[, c("GENE", "P_VALUE")],
                     "./active_snw_search/input_for_search.txt",
                     row.names = FALSE, quote = FALSE, sep = "\t")

  ## Prep for parallel run
  message("## Performing Active Subnetwork Search and Enrichment\n")
  # calculate the number of cores, if necessary
  if (is.null(n_processes))
    n_processes <- parallel::detectCores() - 1
  # Initiate the clusters
  cl <- parallel::makeCluster(n_processes)
  doParallel::registerDoParallel(cl)

  dirs <- rep("", iterations)
  for (i in 1:iterations) {
    dir.create(paste0("./active_snw_search/search", i))
    dirs[i] <- normalizePath(paste0("./active_snw_search/search", i))
  }

  geneInitProbs <- seq(from = 0.01, to = 0.2, length.out = iterations)

  org_dir <- getwd()

  `%dopar%` <- foreach::`%dopar%`
  final_res <- foreach::foreach(i = 1:iterations, .combine = rbind) %dopar% {
    setwd(dirs[i])

    # running Active Subnetwork Search
    system(paste0("java -Xss4m -jar \"", active_search_path, "\"",
                  " -sif=\"", pin_path,"\"",
                  " -sig=../input_for_search.txt",
                  " -method=", search_method,
                  use_all_positives,
                  " -saTemp0=", saTemp0,
                  " -saTemp1=", saTemp1,
                  " -saIter=", format(saIter, scientific = F),
                  " -geneInitProb=", geneInitProbs[i],
                  " -gaPop=", gaPop,
                  " -gaIter=", gaIter,
                  " -gaThread=", gaThread,
                  " -gaMut=", gaMut,
                  " -grMaxDepth=", grMaxDepth,
                  " -grSearchDepth=", grSearchDepth,
                  " -grOverlap=", grOverlap,
                  " -grSubNum=", grSubNum, silent_option))

    # parse
    snws <- pathfindR::parseActiveSnwSearch(
      "resultActiveSubnetworkSearch.txt", input_processed$GENE)

    message(paste0("Found ", length(snws), " active subnetworks\n\n"))

    if (gene_sets == "KEGG") {
      genes_by_pathway <- pathfindR::kegg_genes
      pathways_list <- pathfindR::kegg_pathways
    } else if (gene_sets == "Reactome") {
      genes_by_pathway <- pathfindR::reactome_genes
      pathways_list <- pathfindR::reactome_pathways
    } else if (gene_sets == "BioCarta") {
      genes_by_pathway <- pathfindR::biocarta_genes
      pathways_list <- pathfindR::biocarta_pathways
    } else if (gene_sets == "GO-BP") {
      genes_by_pathway <- pathfindR::go_bp_genes
      pathways_list <- pathfindR::go_bp_pathways
    } else if (gene_sets == "GO-CC") {
      genes_by_pathway <- pathfindR::go_cc_genes
      pathways_list <- pathfindR::go_cc_pathways
    } else if (gene_sets == "GO-MF") {
      genes_by_pathway <- pathfindR::go_mf_genes
      pathways_list <- pathfindR::go_mf_pathways
    } else if (gene_sets == "Custom") {
      genes_by_pathway <- custom_genes
      pathways_list <- custom_pathways
    }

    ## enrichment per subnetwork
    enrichment_res <- lapply(snws, function(x)
      pathfindR::enrichment(genes_by_pathway, x, pathways_list,
                            adj_method, enrichment_threshold, pin_path,
                            DEG_vec = input_processed$GENE))
    enrichment_res <- Reduce(rbind, enrichment_res)

    ## delete non_DEG_Active_Snw_Genes if list_active_snw_genes == FALSE
    if (!list_active_snw_genes)
      enrichment_res$non_DEG_Active_Snw_Genes <- NULL

    if (!is.null(enrichment_res)) {
      ## keep lowest p for each pathway
      idx <- order(enrichment_res$adj_p)
      enrichment_res <- enrichment_res[idx, ]
      enrichment_res <- enrichment_res[!duplicated(enrichment_res$ID), ]
    }

    enrichment_res
  }
  parallel::stopCluster(cl)
  setwd(org_dir)

  if (is.null(final_res)) {
    setwd("..")
    warning("Did not find any enriched pathways!")
    return(data.frame())
  }

  ## Annotate lowest p, highest p and occurrence
  message("## Processing the enrichment results over all iterations \n\n")

  lowest_p <- tapply(final_res$adj_p, final_res$ID, min)
  highest_p <- tapply(final_res$adj_p, final_res$ID, max)
  occurrence <- tapply(final_res$adj_p, final_res$ID, length)

  idx <- match(final_res$ID, names(lowest_p))
  final_res$lowest_p <- as.numeric(lowest_p[idx])

  idx <- match(final_res$ID, names(highest_p))
  final_res$highest_p <- as.numeric(highest_p[idx])

  idx <- match(final_res$ID, names(occurrence))
  final_res$occurrence <- as.numeric(occurrence[idx])

  ## reformat data frame
  keep <- c("ID", "Pathway", "Fold_Enrichment", "occurrence", "lowest_p", "highest_p")
  if (list_active_snw_genes)
    keep <- c(keep, "non_DEG_Active_Snw_Genes")

  final_res <- final_res[, keep]
  final_res <- final_res[order(final_res$lowest_p), ]
  final_res <- final_res[!duplicated(final_res$ID), ]
  rownames(final_res) <- NULL

  message("## Annotating involved genes and visualizing pathways\n\n")
  if (gene_sets == "KEGG") {
    ## Annotate involved genes and generate pathway maps
    genes_df <- input_processed[, c("GENE", "CHANGE")]
    rownames(genes_df) <- genes_df$GENE
    genes_df <- genes_df[, -1, drop = FALSE]
    final_res <- pathmap(final_res, genes_df)
  } else {

    if (gene_sets == "Reactome") {
      genes_by_pathway <- pathfindR::reactome_genes
    } else if (gene_sets == "BioCarta") {
      genes_by_pathway <- pathfindR::biocarta_genes
    } else if (gene_sets == "GO-BP") {
      genes_by_pathway <- pathfindR::go_bp_genes
    } else if (gene_sets == "GO-CC") {
      genes_by_pathway <- pathfindR::go_cc_genes
    } else if (gene_sets == "GO-MF") {
      genes_by_pathway <- pathfindR::go_mf_genes
    } else if (gene_sets == "Custom") {
      genes_by_pathway <- custom_genes
    }

    upreg <- input_processed$GENE[input_processed$CHANGE >= 0]
    downreg <- input_processed$GENE[input_processed$CHANGE < 0]

    final_res$Down_regulated <- final_res$Up_regulated <- NA

    for (i in 1:nrow(final_res)) {
      idx <- which(names(genes_by_pathway) == final_res$ID[i])
      temp <- genes_by_pathway[[idx]]
      final_res$Up_regulated[i] <- paste(temp[temp %in% upreg], collapse = ", ")
      final_res$Down_regulated[i] <- paste(temp[temp %in% downreg], collapse = ", ")
    }

  }


  message("## Creating HTML report\n\n")
  ## Create report
  rmarkdown::render(system.file("rmd/results.Rmd", package = "pathfindR"),
                    output_dir = ".")
  rmarkdown::render(system.file("rmd/all_pathways.Rmd", package = "pathfindR"),
                    params = list(df = final_res, gset = gene_sets), output_dir = ".")
  rmarkdown::render(system.file("rmd/genes_table.Rmd", package = "pathfindR"),
                    params = list(df = input_processed, original_df = input), output_dir = ".")

  setwd("..")

  ## Bubble Chart
  if (bubble) {
    message("Plotting the enrichment bubble chart\n\n")
    graphics::plot(enrichment_chart(final_res))
  }

  message("Pathway enrichment results and converted genes ")
  message("can be found in \"results.html\" ")
  message(paste0("in the folder \"", output_dir, "\"\n\n"))
  message("Run choose_clusters() for clustering pathways\n\n")

  return(final_res)
}

#' Plot the Bubble Chart of Enrichment Results
#'
#' This function is used to plot a bubble chart displaying the enrichment
#' results.
#'
#' @param result_df a data frame that must contain the following columns:\describe{
#'   \item{Pathway}{Description of the enriched pathway}
#'   \item{Fold_Enrichment}{Fold enrichment value for the enriched pathway}
#'   \item{lowest_p}{the lowest adjusted-p value of the given pathway over all iterations}
#'   \item{Up_regulated}{the up-regulated genes in the input involved in the given pathway, comma-separated}
#'   \item{Down_regulated}{the down-regulated genes in the input involved in the given pathway, comma-separated}
#'   \item{Cluster(OPTIONAL)}{the cluster to which the pathway is assigned}
#' }
#' @param plot_by_cluster boolean value indicating whether or not to group the
#' pathways by cluster (works if "Cluster" is a column of `result_df`).
#'
#' @return a `ggplot2` object containing the bubble chart. The x-axis corresponds to
#' fold enrichment values while the y-axis indicates the enriched pathways. Size of
#' the bubble indicates the number of DEGs in the given pathway. Color indicates
#' the -log10(lowest-p) value. The closer the color is to red, the more significant
#' the enrichment is. Optionally, if "Cluster" is a column of `result_df` and
#' plot_by_cluster == TRUE, the pathways are grouped by clusters.
#'
#' @import ggplot2
#' @export
#'
#' @examples
#' g <- enrichment_chart(RA_output)
enrichment_chart <- function(result_df, plot_by_cluster = FALSE) {
  necessary <- c("Pathway", "Fold_Enrichment", "lowest_p", "Up_regulated", "Down_regulated")
  if (!all(necessary %in% colnames(result_df)))
    stop("The input data frame must have the columns: Pathway, Fold_Enrichment, lowest_p, Up_regulated, Down_regulated")

  if (!is.logical(plot_by_cluster))
    stop("plot_by_cluster must be either TRUE or FALSE")

  # sort by lowest adj.p
  result_df <- result_df[order(result_df$lowest_p), ]

  n <- sapply(result_df$Up_regulated, function(x) length(unlist(strsplit(x, ", "))))
  n <- n + sapply(result_df$Down_regulated, function(x) length(unlist(strsplit(x, ", "))))

  result_df$Pathway <- factor(result_df$Pathway, levels = rev(result_df$Pathway))

  g <- ggplot2::ggplot(result_df, ggplot2::aes_(x = ~Fold_Enrichment, y = ~Pathway))
  g <- g + ggplot2::geom_point(ggplot2::aes(color = -log10(result_df$lowest_p),
                                            size = n), na.rm = TRUE)
  g <- g + ggplot2::theme_bw()
  g <- g + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10),
                          axis.text.y = ggplot2::element_text(size = 10),
                          plot.title = ggplot2::element_blank())
  g <- g + ggplot2::xlab("Fold Enrichment") + ggplot2::ylab('')
  g <- g + ggplot2::labs(size = "# of DEGs", color = "-log10(lowest-p)")
  g <- g + ggplot2::scale_color_continuous(low = "#f5efef", high = "red")

  if (plot_by_cluster & "Cluster" %in% colnames(result_df)) {
    g <- g + ggplot2::facet_grid(result_df$Cluster~., scales = "free_y", space = "free", drop = TRUE)
  } else if (plot_by_cluster) {
    warning("For plotting by cluster, there must a column named `Cluster` in the input data frame!")
  }

  return(g)
}

#' Cluster Pathways and Partition the Dendrogram
#'
#' This function first calculates the pairwise distances between the
#' pathways in the \code{result_df} data frame. Next, using this distance
#' matrix, the pathways are clustered via hierarchical clustering. By default,
#' the average silhouette width for each possible number of clusters is
#' calculated. The optimal number of clusters is selected as the one with the
#' highest average silhouette width. The dendrogram is cut into this optimal
#' number of clusters, and the pathways with the lowest p value within each
#' cluster are chosen as representative pathways. If 'auto == FALSE", the user
#' can manually select at which height to cut the dendrogram via a shiny application.
#' See "Chen, Y. A. et al. Integrated pathway clusters with coherent biological
#' themes for target prioritisation. PLoS One 9, e99030,
#' doi:10.1371/journal.pone.0099030 (2014)." for details on the method of
#' pathway clustering.
#'
#' @param result_df data frame of enriched pathways. Must-have columns are: \enumerate{
#'   \item{ID}{KEGG ID of the enriched pathway}
#'   \item{lowest_p}{the lowest adjusted-p value of the given pathway over all iterations}
#'   \item{highest_p}{the highest adjusted-p value of the given pathway over all iterations}
#'   }
#' @param auto boolean value indicating whether to select the optimal number of clusters
#' automatically. If FALSE, a shiny application is displayed, where the user can manually
#' partition the clustering dendrogram (default: TRUE).
#' @param agg_method the agglomeration method to be used if plotting heatmap. Must be one of "ward.D", "ward.D2",
#' "single", "complete", "average", "mcquitty", "median" or "centroid" (default: "average").
#' @param plot_heatmap boolean value indicating whether or not to plot the heat
#'   map of pathway clustering (default: FALSE).
#' @param plot_dend boolean value indicating whether or not to plot the dendrogram
#'   partitioned into the optimal number of clusters, shown by red rectangles (default: FALSE)
#' @param use_names boolean value indicating whether to use gene set names instead of gene set ids (default: FALSE)
#' @param custom_genes a list containing the genes involved in each custom pathway. Each element
#' is a vector of gene symbols located in the given pathway. Names correspond to
#' the ID of the pathway. Must be provided if `result_df` was generated using custom
#' gene sets.
#'
#' @return  If 'auto' is FALSE, manual partitioning can be performed. Via a shiny HTML document, the
#'   hierarchical clustering dendrogram is visualized. In this HTML document,
#'   the user can select the agglomeration method and the distance value at
#'   which to cut the tree. The resulting cluster assignments of the pathways
#'   along with annotation of representative pathways (chosen by smallest lowest
#'   p value) are presented as a table and this table can be saved as a csv
#'   file.
#'   If 'auto' is TRUE, automatic partitioning of clusters is performed. The function
#'   adds 2 additional columns to the input data frame and returns it: \describe{
#'   \item{Cluster}{the cluster to which the pathway is assigned}
#'   \item{Status}{whether the pathway is the "Representative" pathway in its cluster or only a "Member"}
#' }
#'
#' @import fpc
#' @import knitr
#' @import shiny
#' @import rmarkdown
#' @import stats
#' @export
#' @seealso See \code{\link{calculate_pwd}} for calculation of pairwise
#'   distances between enriched pathways. See \code{\link[stats]{hclust}}
#'   for more information on hierarchical clustering. See \code{\link{run_pathfindR}}
#'   for the wrapper function of the pathfindR enrichment workflow.
#'
#' @examples
#' choose_clusters(RA_output)
choose_clusters <- function(result_df, auto = TRUE, agg_method = "average",
                            plot_heatmap = FALSE, plot_dend = FALSE, use_names = FALSE, custom_genes = NULL) {
  ## argument checks
  if (!is.logical(auto))
    stop("The argument `auto` must be either TRUE or FALSE!")

  if (!is.logical(plot_heatmap))
    stop("The argument `plot_heatmap` must be either TRUE or FALSE!")

  valid <- c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")
  if (!agg_method %in% valid)
    stop("`agg_method` must be one of ward.D, ward.D2, single, complete, average, mcquitty, median or centroid!")

  if (!is.logical(plot_heatmap))
    stop("The argument `plot_dend` must be either TRUE or FALSE!")

  ## Check if clustering should be performed

  if (nrow(result_df) < 3)
  {
    warning("There are less than 3 pathways in result_df so clustering is not performed!")
    result_df$Cluster <- 1:nrow(result_df)
    result_df$Status <- "Representative"
    return(result_df)
  }

  ## Create PWD matrix
  message("Calculating pairwise distances between pathways\n\n")
  PWD_mat <- pathfindR::calculate_pwd(result_df$ID, result_df$Pathway,
                                      agg_method = agg_method,
                                      plot_heatmap = plot_heatmap,
                                      use_names,
                                      custom_genes)
  if (!auto) {
    message("Creating the shiny app\n\n")
    parameters <- list(df = result_df, mat = PWD_mat, use_names = use_names)
    rmarkdown::run(system.file("rmd/clustering.Rmd", package = "pathfindR"),
                   render_args = list(output_dir = ".", params = parameters))
  } else {
    ### Calculate PWDs and Cluster
    message("Clustering pathways\n\n")
    hclu <- stats::hclust(as.dist(PWD_mat), method = agg_method)

    ### Optimal k
    message("Calculating the optimal number of clusters (based on average silhouette width)\n\n")
    kmax <- nrow(PWD_mat) - 1
    avg_sils <- c()
    for (i in 2:kmax)
      avg_sils <- c(avg_sils, fpc::cluster.stats(stats::as.dist(PWD_mat),
                                                 stats::cutree(hclu, k = i),
                                                 silhouette = TRUE)$avg.silwidth)

    k_opt <- (2:kmax)[which.max(avg_sils)]
    if (plot_dend) {
      to_plot <- hclu
      if (use_names)
        to_plot$labels <- result_df$Pathway[match(to_plot$labels, result_df$ID)]
      graphics::plot(to_plot, hang = -1)
      stats::rect.hclust(to_plot, k = k_opt)
    }
    message(paste("The maximum average silhouette width was", round(max(avg_sils), 2),
              "for k =", k_opt, "\n\n"))

    ### Return Optimal Clusters
    clusters <- cutree(hclu, k = k_opt)

    result_df$Cluster <- clusters[match(result_df$ID, names(clusters))]
    tmp <- result_df$lowest_p
    names(tmp) <- result_df$ID
    tmp <- tapply(tmp, result_df$Cluster, function(x) names(x)[which.min(x)])
    result_df$Status <- ifelse(result_df$ID %in% tmp, "Representative", "Member")

    result_df <- result_df[order(result_df$Status, decreasing = TRUE), ]
    result_df <- result_df[order(result_df$Cluster), ]

    message("Returning the resulting data frame\n\n")
    return(result_df)
    }
}

#' Return The Path to Given Protein-Protein Interaction Network (PIN)
#'
#' This function returns the path/to/PIN.sif. While the default PINs are
#' Biogrid, GeneMania, IntAct and KEGG, the user can choose to use any other PIN
#' by specifying the path/to/PIN.sif. All PINs to be used in this workflow must
#' have 3 columns with no header and be tab-separated. Columns 1 and 3 must be
#' interacting proteins' HGNC gene symbols, column 2 must be a column with all
#' rows consisting of "pp".
#'
#' @param pin_name_path Name of the chosen PIN or path/to/PIN.sif. If PIN name,
#'   must be one of c("Biogrid", "GeneMania", "IntAct", "KEGG"). If
#'   path/to/PIN.sif, the file must comply with the PIN specifications. Defaults
#'   to "Biogrid".
#'
#' @return A character value that contains the path to chosen PIN.
#'
#' @export
#' @seealso See \code{\link{run_pathfindR}} for the wrapper function of the
#'   pathfindR workflow
#' @examples
#' pin_path <- return_pin_path("Biogrid")
#' pin_path <- return_pin_path("KEGG")

return_pin_path <- function(pin_name_path = "Biogrid") {
  if (pin_name_path %in% c("Biogrid", "GeneMania",
                           "IntAct", "KEGG"))
    path <- normalizePath(system.file(paste0("extdata/", pin_name_path, ".sif"),
                                      package = "pathfindR"))
  else if (file.exists(normalizePath(pin_name_path))) {
    path <- normalizePath(pin_name_path)
    pin <- utils::read.delim(file = path,
                             header = FALSE, stringsAsFactors = FALSE)
    if (ncol(pin) != 3) {
      setwd("..")
      stop("The PIN file must have 3 columns and be tab-separated")
    }

    if (any(pin[, 2] != "pp")) {
      setwd("..")
      stop("The second column of the PIN file must all be \"pp\" ")
    }
  } else {
    setwd("..")
    stop(paste0("The chosen PIN must be one of:\n",
                "Biogrid, GeneMania, IntAct, KEGG or a valid /path/to/SIF"))

  }
  return(path)
}
