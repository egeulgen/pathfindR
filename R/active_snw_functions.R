#' Parse Active Subnetwork Search Output File and Filter the Subnetworks
#'
#' @param active_snw_path path to the output of an Active Subnetwork Search.
#' @param sig_genes_vec vector of significant gene symbols. In the scope of this
#'   package, these are the input genes that were used for active subnetwork search
#' @param score_quan_thr active subnetwork score quantile threshold (Default = 0.80)
#' Must be between 0 and 1 or set to -1 for not filtering
#' @param sig_gene_thr threshold for the minimum proportion of significant genes in
#' the subnetwork (Default = 0.02) If the number of genes to use as threshold is
#' calculated to be < 2 (e.g. 50 signif. genes x 0.01 = 0.5), the threshold number
#' is set to 2
#'
#' @return A list of genes in every active subnetwork that has a score greater than
#' the `score_quan_thr`th quantile and that contains at least `sig_gene_thr` of significant genes.
#' @export
#'
#' @seealso See \code{\link{run_pathfindR}} for the wrapper function of the
#'   pathfindR enrichment workflow
#'
#' @examples
#' path2snw_list <- system.file("extdata/resultActiveSubnetworkSearch.txt",
#'                               package = "pathfindR")
#' filtered <- filterActiveSnws(active_snw_path = path2snw_list,
#'                              sig_genes_vec = pathfindR::RA_input$Gene.symbol)
#' unlink("active_snw_search", recursive = TRUE)
filterActiveSnws <- function(active_snw_path, sig_genes_vec,
                             score_quan_thr = 0.80, sig_gene_thr = 0.02) {
  ## Arg. checks
  active_snw_path <- suppressWarnings(normalizePath(active_snw_path))

  if (!file.exists(active_snw_path))
    stop("The active subnetwork file does not exist! Check the `active_snw_path` argument")

  if (!is.atomic(sig_genes_vec))
    stop("`sig_genes_vec` should be a vector")

  if (!is.numeric(score_quan_thr))
    stop("`score_quan_thr` should be numeric")
  if (score_quan_thr != -1 & (score_quan_thr > 1 | score_quan_thr < 0))
    stop("`score_quan_thr` should be in [0, 1] or -1 (if not filtering)")

  if(!is.numeric(sig_gene_thr))
    stop("`sig_gene_thr` should be numeric")
  if(sig_gene_thr < 0 | sig_gene_thr > 1)
    stop("`sig_gene_thr` should be in [0, 1]")

  output <- readLines(active_snw_path)

  if (length(output) == 0) {
    return(NULL)
  }

  score_vec <- c()
  subnetworks <- list()
  for (i in base::seq_len(length(output))) {
    snw <- output[[i]]

    snw <- unlist(strsplit(snw, "\\s"))

    score_vec <- c(score_vec, as.numeric(snw[1]))
    subnetworks[[i]] <- snw[-1]
  }

  # keep subnetworks with score over the "score_quan_thr"th quantile
  if (score_quan_thr == -1) {
    score_thr <- min(score_vec) - 1
  } else {
    score_thr <- stats::quantile(score_vec, score_quan_thr)
  }
  subnetworks <- subnetworks[as.numeric(score_vec) > as.numeric(score_thr)]

  # select subnetworks containing at least "sig_gene_thr" of significant genes
  snw_sig_counts <- vapply(subnetworks, function(snw_genes)
    sum(base::toupper(snw_genes) %in% base::toupper(sig_genes_vec)), 1)
  sig_gene_num_thr <- sig_gene_thr * length(sig_genes_vec)
  sig_gene_num_thr <- max(2, sig_gene_num_thr)
  cond <- (snw_sig_counts >= sig_gene_num_thr)
  subnetworks <- subnetworks[cond]

  return(subnetworks)
}

#' Perform Active Subnetwork Search
#'
#' @param input_for_search input the input data that active subnetwork search uses. The input
#' must be a data frame containing at least these 2 columns: \describe{
#'   \item{GENE}{Gene Symbol}
#'   \item{P_VALUE}{p value obtained through a test, e.g. differential expression/methylation}
#' }
#' @inheritParams return_pin_path
#' @param snws_file name for active subnetwork search output data
#' \strong{without file extension} (default = "active_snws")
#' @param dir_for_parallel_run (previously created) directory for a parallel run iteration.
#' Used in the wrapper function (see ?run_pathfindR) (Default = NULL)
#' @inheritParams filterActiveSnws
#' @param search_method algorithm to use when performing active subnetwork
#'  search. Options are greedy search (GR), simulated annealing (SA) or genetic
#'  algorithm (GA) for the search (default = "GR").
#' @param silent_option boolean value indicating whether to print the messages
#' to the console (FALSE) or not (TRUE, this will print to a temp. file) during
#' active subnetwork search (default = TRUE). This option was added because
#' during parallel runs, the console messages get disorderly printed.
#' @param use_all_positives if TRUE: in GA, adds an individual with all positive
#'  nodes. In SA, initializes candidate solution with all positive nodes. (default = FALSE)
#' @param geneInitProbs For SA and GA, probability of adding a gene in initial solution (default = 0.1)
#' @param saTemp0 Initial temperature for SA (default = 1.0)
#' @param saTemp1 Final temperature for SA (default = 0.01)
#' @param saIter Iteration number for SA (default = 10000)
#' @param gaPop Population size for GA (default = 400)
#' @param gaIter Iteration number for GA (default = 200)
#' @param gaThread Number of threads to be used in GA (default = 5)
#' @param gaCrossover Applies crossover with the given probability in GA (default = 1, i.e. always perform crossover)
#' @param gaMut For GA, applies mutation with given mutation rate (default = 0, i.e. mutation off)
#' @param grMaxDepth Sets max depth in greedy search, 0 for no limit (default = 1)
#' @param grSearchDepth Search depth in greedy search (default = 1)
#' @param grOverlap Overlap threshold for results of greedy search (default = 0.5)
#' @param grSubNum Number of subnetworks to be presented in the results (default = 1000)
#'
#' @return A list of genes in every identified active subnetwork that has a score greater than
#' the `score_quan_thr`th quantile and that has at least `sig_gene_thr` affected genes.
#'
#' @export
#'
#' @examples
#' processed_df <- RA_input[1:15, -2]
#' colnames(processed_df) <- c("GENE", "P_VALUE")
#' GR_snws <- active_snw_search(input_for_search = processed_df,
#'                              pin_name_path = "KEGG",
#'                              search_method = "GR")
active_snw_search <- function(input_for_search,
                              pin_name_path = "Biogrid",
                              snws_file = "active_snws",
                              dir_for_parallel_run = NULL,
                              score_quan_thr = 0.80, sig_gene_thr = 0.02,
                              search_method = "GR",
                              silent_option = TRUE,
                              use_all_positives = FALSE,
                              geneInitProbs = 0.1,
                              saTemp0 = 1, saTemp1 = 0.01, saIter = 10000,
                              gaPop = 400, gaIter = 10000,
                              gaThread = 5, gaCrossover = 1, gaMut = 0,
                              grMaxDepth = 1, grSearchDepth = 1,
                              grOverlap = 0.5, grSubNum = 1000) {

  ############ Argument checks
  # input_for_search
  if (!is.data.frame(input_for_search)) {
    stop("`input_for_search` should be data frame")
  }
  cnames <- c("GENE", "P_VALUE")
  if (any(!cnames %in% colnames(input_for_search))) {
    stop("`input_for_search` should contain the columns ",
         paste(dQuote(cnames), collapse = ","))
  }

  # pin_name_path (fetch pin path)
  pin_path <- return_pin_path(pin_name_path)

  # snws_file
  if (!suppressWarnings(file.create(file.path(tempdir(check = TRUE),
                                              snws_file)))) {
    stop("`snws_file` may be containing forbidden characters. Please change and try again")
  }

  # search_method
  valid_mets <- c("GR", "SA", "GA")
  if (!search_method %in% valid_mets) {
    stop("`search_method` should be one of ",
         paste(dQuote(valid_mets), collapse = ", "))
  }

  # silent_option
  if (!is.logical(silent_option)) {
    stop("`silent_option` should be either TRUE or FALSE")
  }

  # use_all_positives
  if (!is.logical(use_all_positives)) {
    stop("`use_all_positives` should be either TRUE or FALSE")
  }

  ############ Initial Steps
  ## If dir_for_parallel_run is provided,
  # change working dir to dir_for_parallel_run
  if (!is.null(dir_for_parallel_run)) {
    org_dir <- getwd()
    on.exit(setwd(org_dir))
    setwd(dir_for_parallel_run)
  }

  ## turn silent_option into shell argument
  tmp_out <- file.path(tempdir(check = TRUE), paste0("console_out_", snws_file, ".txt"))
  silent_option <- ifelse(silent_option, paste0(" > ", tmp_out), "")

  ## turn use_all_positives into the java argument
  use_all_positives <- ifelse(use_all_positives, " -useAllPositives", "")

  ## absolute path for active snw search jar
  active_search_jar_path <- system.file("java/ActiveSubnetworkSearch.jar",
                                        package = "pathfindR")

  ## create directory for active subnetworks
  if (!dir.exists("active_snw_search")) {
    dir.create("active_snw_search")
  }

  if (!file.exists("active_snw_search/input_for_search.txt")) {
    input_for_search$GENE <- base::toupper(input_for_search$GENE)
    utils::write.table(input_for_search[, c("GENE", "P_VALUE")],
                       "active_snw_search/input_for_search.txt",
                       col.names = FALSE,
                       row.names = FALSE,
                       quote = FALSE,
                       sep = "\t")
  }

  input_path <- normalizePath("active_snw_search/input_for_search.txt")

  ############ Run active Subnetwork Search
  # running Active Subnetwork Search
  system(paste0(
    "java -Xss4m -jar \"", active_search_jar_path, "\"",
    " -sif=\"", pin_path, "\"",
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
    " -gaCrossover=", gaCrossover,
    " -gaMut=", gaMut,
    " -grMaxDepth=", grMaxDepth,
    " -grSearchDepth=", grSearchDepth,
    " -grOverlap=", grOverlap,
    " -grSubNum=", grSubNum, silent_option
  ))

  snws_file <- file.path("active_snw_search",
                         paste0(snws_file, ".txt"))
  file.rename(from = "resultActiveSubnetworkSearch.txt",
              to = snws_file)

  ############ Parse and filter active subnetworks
  snws <- pathfindR::filterActiveSnws(active_snw_path  = snws_file,
                                      sig_genes_vec = input_for_search$GENE,
                                      score_quan_thr = score_quan_thr,
                                      sig_gene_thr = sig_gene_thr)

  message(paste0("Found ", length(snws), " active subnetworks\n\n"))

  return(snws)
}
