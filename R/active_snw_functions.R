#' Parse Active Subnetwork Search Output File and Filter the Subnetworks
#'
#' @param active_snw_path path to the output of an Active Subnetwork Search.
#' @param signif_genes the vector of significant genes.
#' @param score_quan_thr active subnetwork score quantile threshold (Default = 0.80)
#' @param sig_gene_thr threshold for minimum number of affected genes (Default = 10)
#'
#' @return A list of genes in every active subnetwork that has a score greater than
#' the `score_quan_thr`th quantile and that has at least `sig_gene_thr` affected genes.
#' @export
#'
#' @seealso See \code{\link{run_pathfindR}} for the wrapper function of the
#'   pathfindR workflow
#'
#' @examples
#' \dontshow{
#' filterActiveSnws(normalizePath(system.file("extdata/resultActiveSubnetworkSearch.txt",
#' package = "pathfindR")), pathfindR::RA_input$Gene.symbol)
#' }
#' \dontrun{
#' filterActiveSnws("path/to/output", significant_genes)
#' }
filterActiveSnws <- function(active_snw_path, signif_genes,
                             score_quan_thr = 0.80, sig_gene_thr = 10) {

  output <- readLines(active_snw_path)

  if (length(output) == 0)
    return(NULL)

  score_vec <- c()
  subnetworks <- list()
  for (i in 1:length(output)) {
    snw <- output[[i]]

    snw <- unlist(strsplit(snw, " "))

    score_vec <- c(score_vec, as.numeric(snw[1]))
    subnetworks[[i]] <- snw[-1]
  }

  # keep subnetworks with score over the "score_quan_thr"th quantile
  score_thr <- stats::quantile(score_vec, score_quan_thr)
  subnetworks <- subnetworks[as.numeric(score_vec) > as.numeric(score_thr)]
  # select subnetworks with at least sig_gene_thr significant genes
  snw_sig_counts <- sapply(subnetworks, function(snw) sum(snw %in% signif_genes))
  cond <- (snw_sig_counts >= sig_gene_thr)
  subnetworks <- subnetworks[cond]

  return(subnetworks)
}

#' Perform Active Subnetwork Search
#'
#' @param input_for_search input the input data that active subnetwork search uses. The input
#' must be a data frame containing at least these three columns: \describe{
#'   \item{describe}{HGNC Gene Symbol}
#'   \item{P_VALUE}{p value obtained through a test, e.g. differential expression/methylation}
#' }
#' @param pin_path path to the Protein Interaction Network (PIN) file used in
#'   the analysis
#' @param snws_file name for active subnetwork search output data
#' @param dir_for_parallel_run directory for parallel run iteration.
#' Only used in the wrapper function (see ?run_pathfindR) (Default = NULL)
#' @param score_quan_thr active subnetwork score quantile threshold (Default = 0.80)
#' @param sig_gene_thr threshold for minimum number of significant genes (Default = 10)
#' @param search_method algorithm to use when performing active subnetwork
#'  search. Options are greedy search (GR), simulated annealing (SA) or genetic
#'  algorithm (GA) for the search (Default = GR).
#' @param silent_option boolean value indicating whether to print the messages to the console (FALSE)
#' or print to a file (TRUE) during active subnetwork search (default = TRUE). This option was added
#' because during parallel runs, the console messages get mixed up.
#' @param use_all_positives if TRUE: in GA, adds an individual with all positive
#'  nodes. In SA, initializes candidate solution with all positive nodes. (Default = FALSE)
#' @param geneInitProbs For SA and GA, probability of adding a gene in initial solution (Default = 0.1)
#' @param saTemp0 Initial temperature for SA (Default = 1.0)
#' @param saTemp1 Final temperature for SA (Default = 0.01)
#' @param saIter Iteration number for SA (Default = 10000)
#' @param gaPop Population size for GA (Default = 400)
#' @param gaIter Iteration number for GA (Default = 200)
#' @param gaThread Number of threads to be used in GA (Default = 5)
#' @param gaMut For GA, applies mutation with given mutation rate (Default = 0, i.e. mutation off)
#' @param grMaxDepth Sets max depth in greedy search, 0 for no limit (Default = 1)
#' @param grSearchDepth Search depth in greedy search (Default = 1)
#' @param grOverlap Overlap threshold for results of greedy search (Default = 0.5)
#' @param grSubNum Number of subnetworks to be presented in the results (Default = 1000)
#'
#' @return A list of genes in every identified active subnetwork that has a score greater than
#' the `score_quan_thr`th quantile and that has at least `sig_gene_thr` affected genes.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' active_snw_search(input_for_search, pin_path = "path/to/PIN", search_method = "GR")
#' active_snw_search(input_for_search, pin_path = "path/to/PIN",
#' search_method = "SA", saTemp0 = 2, saTemp1 = 0.05)
#' }
active_snw_search <- function(input_for_search, pin_path,
                              snws_file = "active_snws",
                              dir_for_parallel_run = NULL,
                              score_quan_thr = 0.80, sig_gene_thr = 10,
                              search_method = "GR",
                              silent_option = TRUE,
                              use_all_positives = FALSE,
                              geneInitProbs = 0.1,
                              saTemp0 = 1, saTemp1 = 0.01, saIter = 10000,
                              gaPop = 400, gaIter = 10000, gaThread = 5, gaMut = 0,
                              grMaxDepth = 1, grSearchDepth = 1,
                              grOverlap = 0.5, grSubNum = 1000) {
  ############ Argument checks
  if (!search_method %in% c("GR", "SA", "GA"))
    stop("`search_method` must be one of \"GR\", \"SA\", \"GA\"")

  if (!is.logical(use_all_positives))
    stop("the argument `use_all_positives` must be either TRUE or FALSE")

  if (!is.logical(silent_option))
    stop("the argument `silent_option` must be either TRUE or FALSE")

  ############ Initial Steps
  ## If dir_for_parallel_run is provided, change working dir to dir_for_parallel_run
  if(!is.null(dir_for_parallel_run))
    setwd(dir_for_parallel_run)

  ## turn silent_option into argument
  silent_option <- ifelse(silent_option, paste0(" > ./active_snw_search/console_out_", snws_file, ".txt"), "")

  ## turn use_all_positives into the java argument
  use_all_positives <- ifelse(use_all_positives, " -useAllPositives", "")

  ## absolute path for active snw search jar
  active_search_path <- normalizePath(
    system.file("java/ActiveSubnetworkSearch.jar",
                package = "pathfindR"))

  ## create directory for active subnetworks
  if (!dir.exists("active_snw_search"))
    dir.create("active_snw_search")

  if(!file.exists("active_snw_search/input_for_search.txt")) {
    utils::write.table(input_for_search[, c("GENE", "P_VALUE")],
                       "active_snw_search/input_for_search.txt",
                       col.names = FALSE,
                       row.names = FALSE,
                       quote = FALSE, sep = "\t")
  }

  input_path <- normalizePath("active_snw_search/input_for_search.txt")

  ############ Run active Subnetwork Search
  # running Active Subnetwork Search
  system(paste0("java -Xss4m -jar \"", active_search_path, "\"",
                " -sif=\"", pin_path,"\"",
                " -sig=\"", input_path, "\"",
                " -method=", search_method,
                use_all_positives,
                " -saTemp0=", saTemp0,
                " -saTemp1=", saTemp1,
                " -saIter=", format(saIter, scientific = FALSE),
                " -geneInitProb=", geneInitProbs,
                " -gaPop=", gaPop,
                " -gaIter=", gaIter,
                " -gaThread=", gaThread,
                " -gaMut=", gaMut,
                " -grMaxDepth=", grMaxDepth,
                " -grSearchDepth=", grSearchDepth,
                " -grOverlap=", grOverlap,
                " -grSubNum=", grSubNum, silent_option))

  snws_file <- paste0("active_snw_search/", snws_file, ".txt")
  file.rename(from = "resultActiveSubnetworkSearch.txt",
              to = snws_file)

  ############ Parse and filter active subnetworks
  snws <- pathfindR::filterActiveSnws(
    active_snw_path = snws_file,
    signif_genes = input_for_search$GENE,
    score_quan_thr = score_quan_thr,
    sig_gene_thr = sig_gene_thr)

  message(paste0("Found ", length(snws), " active subnetworks\n\n"))

  ## If dir_for_parallel_run is provided, change working dir back to analysis output dir
  if(!is.null(dir_for_parallel_run))
    setwd("../..")
  return(snws)
}
