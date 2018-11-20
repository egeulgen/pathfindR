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
#'@param visualize_pathways Boolean value to indicate whether or not to create pathway diagrams.
#'@param human_genes boolean to indicate whether the input genes are
#'  human gene symbols or not (default = TRUE)
#'@param enrichment_threshold threshold used when filtering individual iterations' pathway
#'  enrichment results
#'@param adj_method correction method to be used for adjusting p-values of
#'  pathway enrichment results (Default: 'bonferroni', see ?p.adjust)
#'@param search_method algorithm to use when performing active subnetwork
#'  search. Options are greedy search (GR), simulated annealing (SA) or genetic
#'  algorithm (GA) for the search (Default:GR. Can be one of c("GR", "SA",
#'  "GA"))
#'@param use_all_positives if TRUE: in GA, adds an individual with all positive
#'  nodes. In SA, initializes candidate solution with all positive nodes.
#'  (Default = FALSE)
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
#'@param iterations number of iterations for active subnetwork search and
#'  enrichment analyses (Default = 10. Gets set to 1 for Genetic Algorithm)
#'@param n_processes optional argument for specifying the number of processes
#'  used by foreach. If not specified, the function determines this
#'  automatically (Default == NULL. Gets set to 1 for Genetic Algorithm)
#'@inheritParams return_pin_path
#'@param score_quan_thr active subnetwork score quantile threshold (Default = 0.80)
#'@param sig_gene_thr threshold for minimum number of significant genes (Default = 10)
#'@param gene_sets the gene sets to be used for enrichment analysis. Available gene sets
#'  are KEGG, Reactome, BioCarta, GO-All, GO-BP, GO-CC, GO-MF or Custom. If "Custom", the arguments
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
#'@param silent_option boolean value indicating whether to print the messages to the console (FALSE)
#' or print to a file (TRUE) during active subnetwork search (default = TRUE). This option was added
#' because during parallel runs, the console messages get mixed up.
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
#'@section Warning: Especially depending on the protein interaction network,
#'  the algorithm and the number of iterations you choose, active subnetwork
#'  search component of pathfindR may take a very long time to finish.
#'
#'@seealso
#' \code{\link{input_testing}} for input testing, \code{\link{input_processing}} for input processing,
#' \code{\link{active_snw_search}} for active subnetwork search and subnetwork filtering,
#' \code{\link{enrichment_analyses}} for enrichment analysis (using the active subnetworks),
#' \code{\link{summarize_enrichment_results}} for summarizing the active-subnetwork-oriented enrichment results,
#' \code{\link{annotate_pathway_DEGs}} for annotation of affected genes in the given gene sets,
#' \code{\link{visualize_pws}} for visualization of pathway diagrams,
#' \code{\link{enrichment_chart}} for a visual summary of the pathfindR enrichment result,
#' \code{\link[foreach]{foreach}} for details on parallel execution of looping constructs,
#' \code{\link{cluster_pathways}} for clustering the resulting enriched pathways and partitioning into clusters.
#'
#' @examples
#' \dontrun{
#' run_pathfindR(RA_input)
#' }
run_pathfindR <- function(input, p_val_threshold = 5e-2,
                          visualize_pathways = TRUE,
                          human_genes = TRUE,
                          enrichment_threshold = 5e-2,
                          adj_method = "bonferroni",
                          search_method = "GR",
                          use_all_positives = FALSE,
                          saTemp0 = 1, saTemp1 = 0.01, saIter = 10000,
                          gaPop = 400, gaIter = 200, gaThread = 5, gaMut = 0,
                          grMaxDepth = 1, grSearchDepth = 1,
                          grOverlap = 0.5, grSubNum = 1000,
                          iterations = 10, n_processes = NULL,
                          pin_name_path = "Biogrid",
                          score_quan_thr = 0.80, sig_gene_thr = 10,
                          gene_sets = "KEGG",
                          custom_genes = NULL, custom_pathways = NULL,
                          bubble = TRUE,
                          output_dir = "pathfindR_Results",
                          list_active_snw_genes = FALSE,
                          silent_option = TRUE) {
  ############ Argument checks
  # Active Subnetwork Search
  if (!search_method %in% c("GR", "SA", "GA"))
    stop("`search_method` must be one of \"GR\", \"SA\", \"GA\"")

  if (!is.logical(use_all_positives))
    stop("the argument `use_all_positives` must be either TRUE or FALSE")

  if (!is.logical(silent_option))
    stop("the argument `silent_option` must be either TRUE or FALSE")

  # Gene Sets
  if (!gene_sets %in% c("KEGG", "Reactome", "BioCarta",
                        "GO-All", "GO-BP", "GO-CC", "GO-MF", "Custom"))
    stop("`gene_sets` must be one of KEGG, Reactome, BioCarta, GO-All, GO-BP, GO-CC, GO-MF or Custom")

  if (gene_sets == "Custom" & (is.null(custom_genes) | is.null(custom_pathways)))
    stop("You must provide both `custom_genes` and `custom_pathways` if `gene_sets` is `Custom`!")

  # Enrichment chart option
  if (!is.logical(bubble))
    stop("the argument `bubble` must be either TRUE or FALSE")

  ############ Initial Steps
  ## absolute path to PIN
  pin_path <- return_pin_path(pin_name_path)

  ## create output dir
  dir_changed <- FALSE
  output_dir_init <- output_dir
  while(dir.exists(output_dir)) {
    if (grepl("\\(\\d+\\)$" , output_dir)) {
      output_dir <- unlist(strsplit(output_dir, "\\(" ))
      suffix <- as.numeric(sub("\\)", "", output_dir[2])) + 1
      output_dir <- paste0(output_dir[1], "(", suffix, ")")
    } else {
      output_dir <- paste0(output_dir, "(1)")
    }
    dir_changed <- TRUE
  }

  if (dir_changed) {
    warning(paste0("There already is a directory named \"", output_dir_init,
                   "\".\nWriting the result to \"", output_dir,
                   "\" not to overwrite the previous results."))
  }

  org_dir <- getwd()
  dir.create(output_dir, recursive = TRUE)
  output_dir <- normalizePath(output_dir)
  setwd(output_dir)
  output_dir <- getwd()

  ## If search_method is GA, set iterations as 1
  if (search_method == "GA")
    iterations <- n_processes <- 1

  ## If iterations == 1, set n_processes to 1
  if (iterations == 1)
    n_processes <- 1

  ## Set initial probabilities
  geneInitProbs <- seq(from = 0.01, to = 0.2, length.out = iterations)

  ############ Input testing and Processing
  ## Check input
  message("## Testing input\n\n")
  pathfindR::input_testing(input, p_val_threshold, org_dir)

  ## Process input
  message("## Processing input. Converting gene symbols, if necessary (and if human gene symbols provied)\n\n")
  input_processed <- pathfindR::input_processing(input, p_val_threshold, pin_path, org_dir, human_genes)

  ############ Active Subnetwork Search and Enrichment
  ## Prep for parallel run
  message("## Performing Active Subnetwork Search and Enrichment\n")
  # calculate the number of cores, if necessary
  if (is.null(n_processes))
    n_processes <- parallel::detectCores() - 1
  # Initiate the clusters
  cl <- parallel::makeCluster(n_processes)
  doParallel::registerDoParallel(cl)

  dirs <- c()
  for(i in 1:iterations) {
    dir_i <- paste0("active_snw_searches/Iteration_", i)
    dir.create(dir_i, recursive = TRUE)
    dirs <- c(dirs, dir_i)
  }


  `%dopar%` <- foreach::`%dopar%`
  combined_res <- foreach::foreach(i = 1:iterations, .combine = rbind) %dopar% {

    ## Active Subnetwork Search
    snws <- pathfindR::active_snw_search(input_processed, pin_path,
                                         snws_file = paste0("active_snws_", i),
                                         dir_for_parallel_run = dirs[i],
                                         score_quan_thr, sig_gene_thr,
                                         search_method,
                                         silent_option,
                                         use_all_positives,
                                         geneInitProbs[i],
                                         saTemp0, saTemp1, saIter,
                                         gaPop, gaIter, gaThread, gaMut,
                                         grMaxDepth, grSearchDepth,
                                         grOverlap, grSubNum)

    enrichment_res <- pathfindR::enrichment_analyses(snws = snws,
                                                     input_genes = input_processed$GENE,
                                                     gene_sets = gene_sets,
                                                     custom_genes = custom_genes,
                                                     custom_pathways = custom_pathways,
                                                     pin_path = pin_path,
                                                     adj_method = adj_method,
                                                     enrichment_threshold = enrichment_threshold,
                                                     list_active_snw_genes = list_active_snw_genes)

    enrichment_res
  }
  parallel::stopCluster(cl)
  setwd(output_dir)

  ## In case no enrichment was found
  if (is.null(combined_res)) {
    setwd(org_dir)
    warning("Did not find any enriched pathways!")
    return(data.frame())
  }

  ############ Process Enrichment Results of All Iterations
  message("## Processing the enrichment results over all iterations \n\n")
  final_res <- pathfindR::summarize_enrichment_results(combined_res, list_active_snw_genes)

  ############ Annotation of Involved DEGs and Visualization
  message("## Annotating involved genes and visualizing pathways\n\n")

  ##### Annotate Involved DEGs by up/down-regulation status
  final_res <- pathfindR::annotate_pathway_DEGs(final_res, input_processed, gene_sets, custom_genes)

  ##### Visualize the Pathways (If KEGG human, KEGG diagram. Otherwise, Interactions of Genes in the PIN)
  if (visualize_pathways)
    pathfindR::visualize_pws(final_res, input_processed, gene_sets, pin_name_path)

  ############ Create HTML Report
  message("## Creating HTML report\n\n")
  ## Create report
  rmarkdown::render(system.file("rmd/results.Rmd", package = "pathfindR"),
                    output_dir = ".")
  rmarkdown::render(system.file("rmd/enriched_pathways.Rmd", package = "pathfindR"),
                    params = list(df = final_res, gset = gene_sets, vis_cond = visualize_pathways, out_dir = output_dir), output_dir = ".")
  rmarkdown::render(system.file("rmd/conversion_table.Rmd", package = "pathfindR"),
                    params = list(df = input_processed, original_df = input), output_dir = ".")

  setwd(org_dir)

  ############ Bubble Chart
  if (bubble) {
    message("Plotting the enrichment bubble chart\n\n")
    graphics::plot(pathfindR::enrichment_chart(final_res))
  }

  message(paste0("Found ", nrow(final_res), " enriched pathways\n\n"))

  message("Pathway enrichment results and table of converted genes ")
  message("can be found in \"results.html\" ")
  message(paste0("in the folder \"", output_dir, "\"\n\n"))
  message("Run choose_clusters() for clustering pathways\n\n")

  return(final_res)
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
#' @param org_dir path/to/original/directory, supplied by run_pathfindR (default = NULL)
#'
#' @return A character value that contains the path to chosen PIN.
#'
#' @export
#' @seealso See \code{\link{run_pathfindR}} for the wrapper function of the
#'   pathfindR workflow
#' @examples
#' pin_path <- return_pin_path("Biogrid")
#' pin_path <- return_pin_path("KEGG")

return_pin_path <- function(pin_name_path = "Biogrid", org_dir = NULL) {
  if (is.null(org_dir))
    org_dir <- getwd()

  if (pin_name_path %in% c("Biogrid", "GeneMania",
                           "IntAct", "KEGG"))
    path <- normalizePath(system.file(paste0("extdata/", pin_name_path, ".sif"),
                                      package = "pathfindR"))
  else if (file.exists(normalizePath(pin_name_path))) {
    path <- normalizePath(pin_name_path)
    pin <- utils::read.delim(file = path,
                             header = FALSE, stringsAsFactors = FALSE)
    if (ncol(pin) != 3) {
      setwd(org_dir)
      stop("The PIN file must have 3 columns and be tab-separated")
    }

    if (any(pin[, 2] != "pp")) {
      setwd(org_dir)
      stop("The second column of the PIN file must all be \"pp\" ")
    }
  } else {
    setwd(org_dir)
    stop(paste0("The chosen PIN must be one of:\n",
                "Biogrid, GeneMania, IntAct, KEGG or a valid /path/to/SIF"))

  }
  return(path)
}

#' Input Testing
#'
#' @param input the input data that pathfindR uses. The input must be a data
#'   frame with three columns: \enumerate{
#'   \item Gene Symbol (HGNC Gene Symbol)
#'   \item Change value, e.g. log(fold change) (Not obligatory)
#'   \item adjusted p value associated with test, e.g. differential expression/methylation
#' }
#' @param p_val_threshold the adjusted-p value threshold to use when filtering
#'   the input data frame. Must a numeric value between 0 and 1.
#' @param org_dir path/to/original/directory, supplied by run_pathfindR (default = NULL)
#'
#' @return Only checks if the input and the threshold follows the required
#'   specifications.
#' @export
#' @seealso See \code{\link{run_pathfindR}} for the wrapper function of the
#'   pathfindR workflow
#' @examples
#' input_testing(RA_input, 0.05)
input_testing <- function(input, p_val_threshold, org_dir = NULL){
  if (is.null(org_dir))
    org_dir <- getwd()

  if (!is.data.frame(input)) {
    setwd(org_dir)
    stop("the input is not a data frame")
  }

  if (ncol(input) < 2){
    setwd(org_dir)
    stop("There must be at least 2 columns in the input data frame")
  }

  if (!is.numeric(p_val_threshold)){
    setwd(org_dir)
    stop("`p_val_threshold` must be a numeric value between 0 and 1")
  }

  if (p_val_threshold > 1 | p_val_threshold < 0){
    setwd(org_dir)
    stop("`p_val_threshold` must be between 0 and 1")
  }

  p_column <- ifelse(ncol(input) == 3, 3, 2)
  if (!all(is.numeric(input[, p_column]))) {
    setwd(org_dir)
    stop("p values, provided in the third column, must all be numeric")
  }

  if (any(input[, p_column] > 1 | input[, p_column] < 0)) {
    setwd(org_dir)
    stop("p values, provided in the third column, must all be between 0 and 1")
  }

  message("The input looks OK\n\n")
}

#' Process Input
#'
#' @param input the input data that pathfindR uses. The input must be a data
#'   frame with three columns: \enumerate{
#'   \item Gene Symbol (HGNC Gene Symbol)
#'   \item Change value, e.g. log(fold change) (Not obligatory)
#'   \item adjusted p value associated with test, e.g. differential expression/methylation
#' }
#' @param p_val_threshold the adjusted-p value threshold to use when filtering
#'   the input data frame
#' @param pin_path path to the Protein Interaction Network (PIN) file used in
#'   the analysis
#' @param org_dir path/to/original/directory, supplied by run_pathfindR (default = NULL)
#' @param human_genes boolean to indicate whether the input genes are human gene symbols or not (default = TRUE)
#'
#' @return This function first filters the input so that all p values are less
#'   than or equal to the threshold. Next, gene symbols that are not found in
#'   the PIN are identified. If aliases of these gene symbols are found in the
#'   PIN, the symbols are converted to the corresponding aliases. The
#'   resulting data frame containing the original gene symbols, the updated
#'   symbols, change values and p values is then returned.
#' @export
#'
#' @seealso See \code{\link{run_pathfindR}} for the wrapper function of the
#'   pathfindR workflow
#'
#' @examples
#' \dontshow{
#' input_processing(RA_input[1:20,], 0.05, return_pin_path("KEGG"))
#' input_processing(RA_input[1:20,], 0.05, return_pin_path("KEGG"), human_genes = FALSE)
#' }
#' \dontrun{
#' input_processing(RA_input, 0.05, return_pin_path("KEGG"))
#' }
#'
input_processing <- function(input, p_val_threshold, pin_path, org_dir = NULL, human_genes = TRUE) {
  if (is.null(org_dir))
    org_dir <- getwd()

  input <- as.data.frame(input)
  if (ncol(input) == 2)
    input <- data.frame(GENE = input[, 1],
                        CHANGE = rep(100, nrow(input)),
                        P_VALUE = input[, 2])

  colnames(input) <- c("GENE", "CHANGE", "P_VALUE")

  ## Turn GENE into character
  if (is.factor(input$GENE)) {
    warning("The gene column was turned into character from factor.")
    input$GENE <- as.character(input$GENE)
  }

  ## Discard larger than p-value threshold
  input <- input[input$P_VALUE <= p_val_threshold, ]

  ## Choose lowest p for each gene
  if (anyDuplicated(input$GENE)) {
    warning("Duplicated genes found!\nChoosing the lowest p value for each gene")
    input <- input[order(input$P_VALUE, decreasing = FALSE), ]
    input <- input[!duplicated(input$GENE), ]
  }

  ## Fix p < 1e-13
  if (any(input$P_VALUE < 1e-13)) {
    warning("pathfindR cannot handle p values < 1e-13\nThese were changed to 1e-13")
    input$P_VALUE <- ifelse(input$P_VALUE < 1e-13, 1e-13, input$P_VALUE)
  }

  ## load and prep pin
  pin <- utils::read.delim(file = pin_path,
                           header = FALSE, stringsAsFactors = FALSE)
  pin$V2 <- NULL

  ## Genes not in pin
  missing <- input$GENE[!input$GENE %in% c(pin[, 1], pin[, 2])]

  if (human_genes & length(missing) != 0) {
    ## use sql to get alias table and gene_info table (contains the symbols)
    ## first open the database connection
    db_con <- org.Hs.eg.db::org.Hs.eg_dbconn()
    ## write the SQL query
    sql_query <-
      "SELECT * FROM alias, gene_info WHERE alias._id == gene_info._id;"
    ## execute the query on the database
    alias_symbol <- DBI::dbGetQuery(db_con, sql_query)

    select_alias <- function(result, converted, idx) {
      if (idx == 0)
        return("NOT_FOUND")
      else if (result[idx] %in% converted[, 2])
        return(result[idx - 1])
      else
        return(result[idx])
    }

    ## loop for getting all symbols
    converted <- c()
    for (i in 1:length(missing)) {
      result <- alias_symbol[alias_symbol$alias_symbol == missing[i],
                             c("alias_symbol", "symbol")]
      result <- alias_symbol[alias_symbol$symbol %in% result$symbol,
                             c("alias_symbol", "symbol")]
      result <- result$alias_symbol[result$alias_symbol %in%
                                      c(pin[, 1], pin[, 2])]
      ## avoid duplicate entries
      to_add <- select_alias(result, converted, length(result))
      converted <- rbind(converted, c(missing[i], to_add))
    }

    ## Convert to appropriate symbol
    input$new_gene <- input$GENE
    input$new_gene[match(converted[, 1], input$new_gene)] <- converted[, 2]
  } else {
    input$new_gene <- ifelse(input$GENE %in% missing, "NOT_FOUND", input$GENE)
  }

  ## number and percent still missing
  n <- sum(input$new_gene == "NOT_FOUND")
  perc <- n / nrow(input) * 100

  if (n == nrow(input)) {
    setwd(org_dir)
    stop("None of the genes were in the PIN\nPlease check your gene symbols")
  }

  ## Give out warning indicating the number of still missing
  if (n != 0) {
    message(paste0("Could not find any interactions for ",
                   n,
                   " (", round(perc, 2), "%) genes in the PIN\n\n"))
  } else {
    message(paste0("Found interactions for all genes in the PIN\n\n"))
  }

  ## reorder columns
  input <- input[, c(1, 4, 2, 3)]
  colnames(input) <- c("old_GENE", "GENE", "CHANGE", "P_VALUE")

  input <- input[input$GENE != "NOT_FOUND", ]

  ## Keep lowest p value for duplicated genes
  input <- input[order(input$P_VALUE), ]
  input <- input[!duplicated(input$GENE), ]

  return(input)
}

#' Annotate the Affected Genes in the Provided Pathways
#'
#' @param result_df data frame of enrichment results. The only must-have column is "ID".
#' @param input_processed input data processed via `input_processing`
#' @param gene_sets the gene sets used for enrichment analysis. Possible gene sets
#'  are KEGG, Reactome, BioCarta, GO-All, GO-BP, GO-CC, GO-MF or Custom (Default = "KEGG"). If "Custom", the
#'  custom_genes argument must be specified.
#' @param custom_genes a list containing the genes involved in each custom gene set. Each element
#' is a vector of gene symbols located in the given pathway. Names correspond to
#' the ID of the pathway.
#'
#' @return The original data frame with two additional columns:  \describe{
#'   \item{Up_regulated}{Up-regulated input genes in the given pathway}
#'   \item{Down_regulated}{Down-regulated input genes in the given pathway}
#' }
#' @export
#'
#' @examples
#' example_gene_data <- RA_input
#' colnames(example_gene_data) <- c("GENE", "CHANGE", "P_VALUE")
#'
#' annotated_result <- annotate_pathway_DEGs(RA_output, example_gene_data)
#'
annotate_pathway_DEGs <- function(result_df, input_processed, gene_sets = "KEGG", custom_genes = NULL) {
  ############ Load Gene Set Data
  gene_sets_df <- data.frame("Gene Set Name" = c("KEGG", "Reactome", "BioCarta",
                                                 "GO-All", "GO-BP", "GO-CC", "GO-MF"),
                             "genes_by_pathway" = c("kegg_genes", "reactome_genes", "biocarta_genes",
                                                    "go_all_genes", "go_bp_genes", "go_cc_genes", "go_mf_genes"),
                             "pathways_list" = c("kegg_pathways", "reactome_pathways", "biocarta_pathways",
                                                 "go_all_pathways", "go_bp_pathways", "go_cc_pathways", "go_mf_pathways"))

  if (gene_sets %in% gene_sets_df$Gene.Set.Name) {
    idx <- which(gene_sets_df$Gene.Set.Name == gene_sets)
    genes_name <- gene_sets_df$genes_by_pathway[idx]
    genes_by_pathway <- base::eval(parse(text = paste0("pathfindR::", genes_name)))
  } else if (gene_sets == "Custom") {
    genes_by_pathway <- custom_genes
  }

  ############ Annotate up/down-regulated pathway genes

  ## Up/Down-regulated genes
  upreg <- input_processed$GENE[input_processed$CHANGE >= 0]
  downreg <- input_processed$GENE[input_processed$CHANGE < 0]

  ## Annotation
  annotated_df <- result_df
  annotated_df$Down_regulated <- annotated_df$Up_regulated <- NA
  for (i in 1:nrow(annotated_df)) {
    idx <- which(names(genes_by_pathway) == annotated_df$ID[i])
    temp <- genes_by_pathway[[idx]]
    annotated_df$Up_regulated[i] <- paste(temp[temp %in% upreg], collapse = ", ")
    annotated_df$Down_regulated[i] <- paste(temp[temp %in% downreg], collapse = ", ")
  }

  return(annotated_df)
}
