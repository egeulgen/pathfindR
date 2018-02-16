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
#'the lowest and the highest adjusted-p values, as well as number of occurences
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
#'@param use_all_positives if TRUE in GA, adds an individual with all positive
#'  nodes. In SA, initializes candidate solution with all positive nodes.
#'  (Default = FALSE)
#'@param saTemp0 initial temperature for SA (Default: 1.0)
#'@param saTemp1 final temperature for SA (Default: 0.01)
#'@param saIter iteration number for SA (Default: 10000)
#'@param gaPop population size for GA (Default: 400)
#'@param gaIter iteration number for GA (Default: 10000)
#'@param gaThread number of threads to be used in GA (Default: 5)
#'@param gaMut the mutation rate for GA (Default: 0)
#'@param grMaxDepth sets max depth in greedy search, 0 for no limit (Default: 2)
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
#'@inheritParams current_KEGG
#'@inheritParams return_pin_path
#'
#'@import knitr
#'
#'@return Data frame of pathview enrichment results. Columns are: "ID",
#'  "Pathway", "occurence", "lowest_p", "highest_p". "ID" is the KEGG ID for a
#'  given pathway and "Pathway" is the name. "occurence" indicates the number of
#'  times the given pathway is encountered over all iterations. "lowest_p" and
#'  "highest_p" indicate the lowest and highest adjusted p values over all
#'  iterations. The function also creates an HTML report with the pathview
#'  enrichment results linked to the visualizations of the pathways in addition
#'  to the table of converted gene symbols. This report can be found in
#'  "results.html".
#'
#'@export
#'
#'@section Warning: Depending on the protein interaction network of your choice,
#'  active subnetwork finding component of pathfindR may take a very long time
#'  to finish. Therefore, overnight runs are recommended.
#'
#'@seealso \code{\link{input_testing}} for input testing,
#'  \code{\link{input_processing}} for input processing,
#'  \code{\link{current_KEGG}} for KEGG pathway genes retrieval,
#'  \code{\link{parseActiveSnwSearch}} for parsing a jActive modules output,
#'  \code{\link{enrichment}} for pathway enrichment analysis and
#'  \code{\link{pathmap}} for annotation of involved genes and visualization of
#'  pathways. See \code{\link[foreach]{foreach}} for details on parallel
#'  execution of looping constructs. See \code{\link{choose_clusters}} for
#'  clustering the resulting enriched pathways.
#'
#' @examples
#' \dontrun{
#' run_pathfindR(RA_input)
#' }
run_pathfindR <- function(input, p_val_threshold = 5e-2,
                          enrichment_threshold = 0.05,
                          adj_method = "bonferroni",
                          search_method = "GR",
                          use_all_positives = FALSE,
                          saTemp0 = 1, saTemp1 = 0.01, saIter = 10000,
                          gaPop = 400, gaIter = 10000, gaThread = 5, gaMut = 0,
                          grMaxDepth = 1, grSearchDepth = 1,
                          grOverlap = 0.5, grSubNum = 1000,
                          iterations = 10, n_processes = NULL,
                          kegg_update = FALSE, pin_name_path = "Biogrid") {

  if (!search_method %in% c("GR", "SA", "GA"))
    stop("search_method must be one of \"GR\", \"SA\", \"GA\"")
  if (!is.logical(use_all_positives))
    stop("search_method must be logical")

  if (search_method != "GR")
    iterations <- n_processes <- 1

  use_all_positives <- ifelse(use_all_positives, " -useAllPositives", "")

  ## absolute paths for cytoscape and pin
  active_search_path <- normalizePath(
    system.file("java/ActiveSubnetworkSearch.jar",
                package = "pathfindR"))
  pin_path <- return_pin_path(pin_name_path)

  ## Check input
  cat("## Testing input\n\n")
  input_testing(input, p_val_threshold) # perform input testing

  ## Process input
  cat("## Processing input. Converting gene symbols, if necessary\n\n")
  input_processed <- input_processing(input, p_val_threshold, pin_path)

  dir.create("active_snw_search")
  utils::write.table(input_processed[, c("GENE", "P_VALUE")],
                     "./active_snw_search/input_for_search.txt",
                     row.names = FALSE, quote = FALSE, sep = "\t")

  ## get current KEGG pathways and kegg id, pathway names
  cat("## Retreiving most current KEGG pathway genes\n\n")
  pw_genes <- current_KEGG(kegg_update)
  pathways_list <- KEGGREST::keggList("pathway", "hsa")

  ## Prep for parallel run
  cat("## Performing Active Subnetwork Search and Enrichment\n")

  if (is.null(n_processes))
    n_processes <- parallel::detectCores()
  cl <- snow::makeCluster(n_processes)
  doSNOW::registerDoSNOW(cl)

  dirs <- rep("", iterations)
  for (i in 1:iterations) {
    dir.create(paste0("./active_snw_search/search", i))
    dirs[i] <- normalizePath(paste0("./active_snw_search/search", i))
  }

  `%dopar%` <- foreach::`%dopar%`
  final_res <- foreach::foreach(i = 1:iterations, .combine = rbind) %dopar% {
    setwd(dirs[i])

    # running jactivemodules
    system(paste0("java -Xss4m -jar ", active_search_path,
                  " -sif=", pin_path,
                  " -sig=../input_for_search.txt",
                  " -method=", search_method,
                  use_all_positives,
                  " -saTemp0=", saTemp0,
                  " -saTemp1=", saTemp1,
                  " -saIter=", format(saIter, scientific = F),
                  " -gaPop=", gaPop,
                  " -gaIter=", gaIter,
                  " -gaThread=", gaThread,
                  " -gaMut=", gaMut,
                  " -grMaxDepth=", grMaxDepth,
                  " -grSearchDepth=", grSearchDepth,
                  " -grOverlap=", grOverlap,
                  " -grSubNum=", grSubNum))

    # parse
    snws <- pathfindR::parseActiveSnwSearch(
      "resultActiveSubnetworkSearch.txt", input_processed$GENE)

    cat(paste0("Found ", length(snws), " active subnetworks\n\n"))

    ## enrichment per subnetwork
    enrichment_res <- lapply(snws, function(x)
      pathfindR::enrichment(pw_genes, x, pathways_list,
                            adj_method, enrichment_threshold, pin_path))
    enrichment_res <- Reduce(rbind, enrichment_res)

    if (!is.null(enrichment_res)) {
      ## keep lowest p for each pathway
      idx <- order(enrichment_res$adj_p)
      enrichment_res <- enrichment_res[idx, ]
      enrichment_res <- enrichment_res[!duplicated(enrichment_res$ID), ]
    }

    enrichment_res
  }
  snow::stopCluster(cl)

  if (is.null(final_res))
    stop("Did not find any enriched pathways!")

  ## Annotate lowest p, highest p and occurence
  cat("## Processing the enrichment results over all iterations \n\n")

  lowest_p <- tapply(final_res$adj_p, final_res$ID, min)
  highest_p <- tapply(final_res$adj_p, final_res$ID, max)
  occurence <- tapply(final_res$adj_p, final_res$ID, length)

  idx <- match(final_res$ID, names(lowest_p))
  final_res$lowest_p <- as.numeric(lowest_p[idx])

  idx <- match(final_res$ID, names(highest_p))
  final_res$highest_p <- as.numeric(highest_p[idx])

  idx <- match(final_res$ID, names(occurence))
  final_res$occurence <- as.numeric(occurence[idx])

  ## reformat data frame
  keep <- c("ID", "Pathway", "occurence", "lowest_p", "highest_p")
  final_res <- final_res[, keep]
  final_res <- final_res[order(final_res$lowest_p), ]
  final_res <- final_res[!duplicated(final_res$ID), ]
  rownames(final_res) <- NULL

  cat("## Annotating involved genes and visualizing pathways\n\n")
  ## Annotate involved genes and generate pathway maps
  genes_df <- input_processed[, c("GENE", "CHANGE")]
  rownames(genes_df) <- genes_df$GENE
  genes_df <- genes_df[, -1, drop = FALSE]
  final_res <- pathmap(final_res, genes_df)

  cat("## Creating HTML report\n\n")
  ## Create report
  rmarkdown::render(system.file("rmd/results.Rmd", package = "pathfindR"),
                    output_dir = ".")
  rmarkdown::render(system.file("rmd/all_pathways.Rmd", package = "pathfindR"),
                    params = list(df = final_res), output_dir = ".")

  rmarkdown::render(system.file("rmd/genes_table.Rmd", package = "pathfindR"),
                    params = list(df = input_processed), output_dir = ".")

  cat("Pathway enrichment results and converted genes ")
  cat("can be found in \"results.html\"\n\n")
  cat("Run choose_clusters() for clustering pathways\n\n")

  return(final_res)
}

#' Cluster Pathways and Dynamically Cut the Dendrogram
#'
#' See "Chen, Y. A. et al. Integrated pathway clusters with coherent biological
#' themes for target prioritisation. PLoS One 9, e99030,
#' doi:10.1371/journal.pone.0099030 (2014)." for details on the method of
#' pathway clustering.
#'
#' @param result_df resulting data frame of the pathfindR main workflow.
#' @param ... optional arguments for \code{cluster_pathways}
#'
#' @return This function first calculates the pairwise distances between the
#'   pathways in the \code{result_df} data frame. Via a shiny HTML document, the
#'   hierarchical clustering dendrogram is visualized. In this HTML document,
#'   the user can select the agglomeration method and the distance value at
#'   which to cut the tree. The resulting cluster assignments of the pathways
#'   along with annotation of representative pathways (chosen by smallest lowest
#'   p value) are presented as a table and this table can be saved as a csv
#'   file.
#'
#' @export
#' @import knitr
#' @import shiny
#'
#' @seealso See \code{\link{cluster_pathways}} for calculation of pairwise distances
#' between enriched pathways.
#'
#' @examples
#' \dontrun{
#' choose_clusters(RA_output)
#' }
choose_clusters <- function(result_df, ...) {
  cat("Calculating pairwise distances between pathways\n\n")
  PWD_mat <- cluster_pathways(result_df$ID, ...)

  cat("Creating shiny app\n\n")
  parameters <- list(df = result_df, mat = PWD_mat)
  rmarkdown::run(system.file("rmd/clustering.Rmd", package = "pathfindR"),
                 render_args = list(output_dir = ".", params = parameters))
}

#' Return The Path to Given Protein Interaction Network (PIN)
#'
#' This function returns the path/to/PIN.sif. While the default PINs are
#' Biogrid, GeneMania, IntAct and KEGG, the user can choose to use any other PIN
#' by specifying the path/to/PPI.sif. All PINs to be used in this workflow must
#' have 3 columns with no header and be tab-separated. Columns 1 and 3 must be
#' interactor proteins' HGNC gene symbols, column 2 must be a column with all
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
    if (ncol(pin) != 3)
      stop("The PIN file must have 3 columns and be tab-separated")
    if (any(pin[, 2] != "pp"))
      stop("The second column of the PIN file must all be \"pp\" ")
  }

  else
    stop(paste0("The chosen PIN must be one of:\n",
                "Biogrid, GeneMania, IntAct or KEGG"))

  return(path)
}
