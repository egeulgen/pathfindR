.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to pathfindr")
}

#' Wrapper Function for pathfindr
#'
#' @inheritParams input_testing
#' @param enrichment_threshold threshold used when filtering individual pathway
#'   enrichment results
#' @param adj_method correction method to be used for adjusting p-values of
#'   pathway enrichment results
#' @param iterations number of iterations for active subnetwork search and
#'   enrichment analyses
#' @param ncores optional argument for specifying the number of cores. The
#'   function determines this automatically
#' @param n_snw optional argument for specifying the number of active
#'   subnetworks when running jActive modules. Defaults to 1000.
#' @param overlap_threshold optional argument for specifying the overlap
#'   thresholds when running jActive modules. Defaults to 0.5.
#'
#' @return Data frame of pathview enrichment results. Columns are: "ID",
#'   "Pathway", "occurence", "lowest_p", "highest_p". "ID" is the KEGG ID for a
#'   given pathway and "Pathway" is the name. "occurence" indicates the number
#'   of times the given pathway is encountered over all iterations. "lowest_p"
#'   and "highest_p" indicate the lowest and highest adjusted p values over all
#'   iterations.
#' @export
#'
#' @examples
#' pathfindr(input_data_frame)
pathfindr <- function(input, p_val_threshold = 0.05,
                      enrichment_threshold = 1e-4,
                      adj_method = "bonferroni",
                      iterations = 10, ncores = NULL,
                      n_snw = 1000, overlap_threshold = 0.5) {
  ## absolute paths for cytoscape and ppi
  jactive_path <- system.file("java/myCytoscape.jar", package = "pathfindr")
  ppi_path <- system.file("data/humanPPI.sif", package = "pathfindr")
  package_dir <- system.file(package = "pathfindr")

  ## Check input
  cat("Testing input\n\n")
  input_testing(input, p_val_threshold) # perform input testing

  ## Process input
  cat("Processing input. Converting gene symbols, if necessary\n\n")
  input_processed <- input_processing(input, p_val_threshold, ppi_path)
  write.table(input_processed[, c("GENE", "SPOTPvalue")],
              "./input_processed.txt", row.names = F, quote = F, sep = "\t")

  ## get current KEGG pathways and kegg id, pathway names
  cat("Retreiving most current KEGG pathway genes\n\n")
  pw_genes <- current_KEGG()
  pathways_list <- KEGGREST::keggList("pathway", "hsa")

  ## Prep for parallel run
  cat("Running jActive modules and enrichment\n")
  cat("Any java window that opens will close once the task is finished\n")
  cat("DO NOT close the java window(s)!\n\n")

  if (is.null(ncores))
    ncores <- parallel::detectCores()
  cl <- snow::makeCluster(ncores)
  doSNOW::registerDoSNOW(cl)

  dir.create("jActive")
  dirs <- rep("", iterations)
  for (i in 1:iterations) {
    dir.create(paste0("./jActive/jActive",i))
    dirs[i] <- normalizePath(paste0("./jActive/jActive", i))
  }

  `%dopar%` <- foreach::`%dopar%`
  final_res <- foreach::foreach(i = 1:iterations, .combine = rbind) %dopar% {
    devtools::load_all(package_dir)
    setwd(dirs[i])

    # running jactivemodules
    system(paste0("java -jar ", jactive_path,
                  " -N ", ppi_path,
                  " -m ../../input_processed.txt",
                  " -so ./jActive.txt -np ", n_snw,
                  " -ot ", overlap_threshold))

    # parse
    jactive_output <- read.table("jActive.txt", stringsAsFactors = F)
    snws <- pathfindr::parsejActive(jactive_output, input_processed$GENE)

    cat(paste0("Found ", length(snws), " active subnetworks\n\n"))

    ## enrichment per subnetwork
    enrichment_res <- lapply(snws, function(x)
        pathfindr::enrichment(pw_genes, x, pathways_list, adj_method))
    enrichment_res <- Reduce(rbind, enrichment_res)

    ## remove p larger than enrichment_threshold and keep lowest p for each pathway
    cond <- enrichment_res$adj_p <= enrichment_threshold
    enrichment_res <- enrichment_res[cond, ]
    idx <- order(enrichment_res$adj_p, enrichment_res$ID)
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

  return(final_res)
}

#' Cluster Pathways and Dynamically Cut the Dendrogram
#'
#' @param result_df resulting data frame of the pathfindr main workflow.
#'
#' @return This function first calculates the pairwise distances between the
#'   pathways in the \code{result_df} data frame. Via a shiny HTML document, the
#'   hierarchical clustering dendrogram is visualized. In this HTML document,
#'   the user can select the value at which to cut the tree and the resulting
#'   representative pathways are presented as a table and pathways with cluster
#'   assignments are saved as a csv file to the current directory
#' @export
#'
#' @examples
#' choose_clusters(result_df)
choose_clusters <- function(result_df) {
  cat("Calculating pairwise distances between pathways\n\n")
  PWD_mat <- cluster_pathways(result_df$ID)

  cat("Creating shiny app\n\n")
  parameters <- list(df = result_df, mat = PWD_mat)
  rmarkdown::run(system.file("rmd/clustering.Rmd", package = "pathfindr"),
                 render_args = list(output_dir = ".", params = parameters))
}
