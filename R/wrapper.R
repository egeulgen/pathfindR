.onAttach <- function(libname, pkgname) {
  packageStartupMessage("##############################################################################
                            Welcome to pathfindr
##############################################################################")
}

#' Wrapper Function for pathfindr Workflow
#'
#' \code{run_pathfindr} is the wrapper function for the pathfindr workflow
#'
#' This function takes in a data frame consisting of Gene Symbol,
#' log-fold-change and adjusted-p values. After input testing, any gene symbols
#' that are not in the PIN are converted to alias symbols if the alias is in the
#' PIN. Next, active subnetwork search is performed. Pathway enrichment analysis
#' is performed using the genes in each of the active subnetworks. Pathways with
#' adjusted-p values lower than \code{enrichment_threshold} are discarded. The
#' lowest adjusted-p value (over all subnetworks) for each pathway is kept. This
#' process of active subnetwork search and enrichment is repeated  for a
#' selected number of \code{iterations}, which is done in parallel. Over all
#' iterations, the lowest and the highest adjusted-p values, as well as number
#' of occurences are reported for each enriched pathway.
#'
#' @inheritParams input_testing
#' @param enrichment_threshold threshold used when filtering individual pathway
#'   enrichment results
#' @param adj_method correction method to be used for adjusting p-values of
#'   pathway enrichment results
#' @param iterations number of iterations for active subnetwork search and
#'   enrichment analyses
#' @param n_processes optional argument for specifying the number of processes
#'   used by foreach. The function determines this automatically
#' @param n_snw optional argument for specifying the number of active
#'   subnetworks when running jActive modules. (Default = 5)
#' @param overlap_threshold optional argument for specifying the overlap
#'   thresholds when running jActive modules. (Default = 0.8)
#' @inheritParams current_KEGG
#' @inheritParams return_pin_path
#'
#' @return Data frame of pathview enrichment results. Columns are: "ID",
#'   "Pathway", "occurence", "lowest_p", "highest_p". "ID" is the KEGG ID for a
#'   given pathway and "Pathway" is the name. "occurence" indicates the number
#'   of times the given pathway is encountered over all iterations. "lowest_p"
#'   and "highest_p" indicate the lowest and highest adjusted p values over all
#'   iterations. The function also creates an HTML report with the pathview
#'   enrichment results linked to the visualizations of the pathways in addition
#'   to the table of converted gene symbols. This report can be found in
#'   "results.html".
#'
#' @export
#'
#' @section Warning: Depending on the protein interaction network of your
#'   choice, active subnetwork finding component of pathfindr may take a very
#'   long time to finish. Therefore, overnight runs are recommended.
#'
#' @seealso \code{\link{input_testing}} for input testing,
#'   \code{\link{input_processing}} for input processing,
#'   \code{\link{current_KEGG}} for KEGG pathway genes retrieval,
#'   \code{\link{parsejActive}} for parsing a jActive modules output,
#'   \code{\link{enrichment}} for pathway enrichment analysis and
#'   \code{\link{pathmap}} for annotation of involved genes and visualization of
#'   pathways. See \code{\link[foreach]{foreach}} for details on parallel
#'   execution of looping constructs. See \code{\link{choose_clusters}} for
#'   clustering the resulting enriched pathways.
#'
#' @examples
#' run_pathfindr(input_data_frame)
run_pathfindr <- function(input, p_val_threshold = 5e-2,
                      enrichment_threshold = 1e-4,
                      adj_method = "bonferroni",
                      iterations = 10, n_processes = NULL,
                      n_snw = 5, overlap_threshold = 0.8,
                      kegg_update = F, pin_name = "KEGG") {
  ## absolute paths for cytoscape and pin
  jactive_path <- normalizePath(system.file("java/myCytoscape.jar",
                                            package = "pathfindr"))
  pin_path <- return_pin_path(pin_name)

  ## Check input
  cat("Testing input\n\n")
  input_testing(input, p_val_threshold) # perform input testing

  ## Process input
  cat("Processing input. Converting gene symbols, if necessary\n\n")
  input_processed <- input_processing(input, p_val_threshold, pin_path)
  write.table(input_processed[, c("GENE", "SPOTPvalue")],
              "./input_for_jactive.txt", row.names = F, quote = F, sep = "\t")

  ## get current KEGG pathways and kegg id, pathway names
  cat("Retreiving most current KEGG pathway genes\n\n")
  pw_genes <- current_KEGG(kegg_update)
  pathways_list <- KEGGREST::keggList("pathway", "hsa")

  ## Prep for parallel run
  cat("Running jActive modules and enrichment\n")
  cat("Any java window that opens will close once the task is finished\n")
  cat("DO NOT close the java window(s)!\n\n")

  if (is.null(n_processes))
    n_processes <- parallel::detectCores()
  cl <- snow::makeCluster(n_processes)
  doSNOW::registerDoSNOW(cl)

  dir.create("jActive")
  dirs <- rep("", iterations)
  for (i in 1:iterations) {
    dir.create(paste0("./jActive/jActive", i))
    dirs[i] <- normalizePath(paste0("./jActive/jActive", i))
  }

  `%dopar%` <- foreach::`%dopar%`
  final_res <- foreach::foreach(i = 1:iterations, .combine = rbind) %dopar% {
    setwd(dirs[i])

    # running jactivemodules
    system(paste0("java -jar ", jactive_path,
                  " -N ", pin_path,
                  " -m ../../input_for_jactive.txt",
                  " -so ./jActive.txt -np ", n_snw,
                  " -ot ", overlap_threshold))

    # parse
    jactive_output <- read.table("jActive.txt", stringsAsFactors = F)
    snws <- pathfindr::parsejActive(jactive_output, input_processed$GENE)

    cat(paste0("Found ", length(snws), " active subnetworks\n\n"))

    ## enrichment per subnetwork
    enrichment_res <- lapply(snws, function(x)
        pathfindr::enrichment(pw_genes, x, pathways_list, adj_method, enrichment_threshold, pin_path))
    enrichment_res <- Reduce(rbind, enrichment_res)

    ## keep lowest p for each pathway
    idx <- order(enrichment_res$adj_p)
    enrichment_res <- enrichment_res[idx, ]
    enrichment_res <- enrichment_res[!duplicated(enrichment_res$ID), ]

    enrichment_res
  }
  snow::stopCluster(cl)

  ## Annotate lowest p, highest p and occurence
  cat("Processing the enrichment results over all iterations \n\n")

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

  cat("Annotating involved genes and generating pathway diagrams\n\n")
  ## Annotate involved genes and generate pathway maps
  genes_df <- input_processed[, c("GENE", "CHANGE")]
  rownames(genes_df) <- genes_df$GENE
  genes_df <- genes_df[, -1, drop = F]
  final_res <- pathmap(final_res, genes_df)

  cat("Creating HTML report\n\n")
  ## Create report
  rmarkdown::render(system.file("rmd/results.Rmd", package = "pathfindr"),
                    output_dir = ".")
  rmarkdown::render(system.file("rmd/all_pathways.Rmd", package = "pathfindr"),
                    params = list(df = final_res), output_dir = ".")

  rmarkdown::render(system.file("rmd/genes_table.Rmd", package = "pathfindr"),
                    params = list(df = input_processed), output_dir = ".")

  cat("Pathway enrichment results and converted genes
can be found in \"results.html\"\n\n")
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
#' @param result_df resulting data frame of the pathfindr main workflow.
#' @inheritParams cluster_pathways
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
#'
#' @examples
#' choose_clusters(result_df)
choose_clusters <- function(result_df, ...) {
  cat("Calculating pairwise distances between pathways\n\n")
  PWD_mat <- cluster_pathways(result_df$ID, ...)

  cat("Creating shiny app\n\n")
  parameters <- list(df = result_df, mat = PWD_mat)
  rmarkdown::run(system.file("rmd/clustering.Rmd", package = "pathfindr"),
                 render_args = list(output_dir = ".", params = parameters))
}

#' Return The Path to Given Protein Interaction Network (PIN)
#'
#' @param pin_name Name of the chosen PIN. Must be one of c("Biogrid", "STRING",
#'   "GeneMania", "BioPlex", "HitPredict", "IntAct", "KEGG"). Defaults to "KEGG".
#'
#' @return A character value that contains the path to chosen PIN.
#'
#' @export
#' @seealso See \code{\link{run_pathfindr}} for the wrapper function of the
#'   pathfindr workflow
#' @examples
#' pin_path <- return_pin_path("Biogrid")

return_pin_path <- function(pin_name = "KEGG") {
  if (!pin_name %in% c("Biogrid", "STRING", "GeneMania", "BioPlex",
                       "HitPredict", "IntAct", "KEGG"))
    stop(paste0("The chosen PIN must be one of:\n",
                "Biogrid, STRING, GeneMania, BioPlex, HitPredict, IntAct or KEGG"))

  path <- normalizePath(system.file(paste0("data/", pin_name, ".sif"),
                                    package = "pathfindr"))
  return(path)
}
