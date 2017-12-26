pathfindr <- function(input, p_val_threshold = 0.05,
                      enrichment_threshold = 1e-4,
                      adj_method = "bonferroni",
                      iterations = 10, ncores = NULL,
                      n_snw = 1000, overlap_threshold = 0.5) {
  ## absolute paths for cytoscape and ppi
  cytoscape_path <- system.file("java/myCytoscape.jar", package = "pathfindr")
  ppi_path <- system.file("data/humanPPI.sif", package = "pathfindr")
  package_dir <- system.file(package = "pathfindr")

  ## Process inputs
  cat("Processing input\n\n")
  input_processed <- input_processing(input, p_val_threshold, ppi_path)
  write.table(input_processed[, c("GENE", "SPOTPvalue")],
              "./input_processed.txt", row.names = F, quote = F, sep = "\t")

  ## get current KEGG pathways
  cat("Retreiving most current KEGG pathway genes\n\n")
  pw_genes <- current_KEGG()
  pathways_list <- KEGGREST::keggList("pathway", "hsa")

  ## Prep for parallel run
  cat("Running jActive modules and enrichment\n")
  cat("Any java window that opens will close once the task is finished\n")
  cat("DO NOT close the java window(s)!\n\n")

  suppressPackageStartupMessages(library(foreach))

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

  final_res <- foreach(i = 1:iterations, .combine = rbind) %dopar% {
    devtools::load_all(package_dir)
    setwd(dirs[i])

    # running jactivemodules
    system(paste0("java -jar ", cytoscape_path,
                  " -N ", ppi_path,
                  " -m ../../input_processed.txt",
                  " -so ./jActive.txt -np ", n_snw,
                  " -ot ", overlap_threshold))

    # parse
    jactive_output <- read.table("jActive.txt", stringsAsFactors = F)
    cat(paste0(input_processed$GENE[1], "WORKING\n")) # sanity check
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
  cat("Processing the enrichment results\n\n")

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

  cat("Annotating involved genes and generating pathway diagrams\n")
  ## Create report
  rmarkdown::render(system.file("rmd/results.Rmd", package = "pathfindr"),
                    output_dir = ".")
  rmarkdown::render(system.file("rmd/all_pathways.Rmd", package = "pathfindr"),
                    params = list(df = final_res), output_dir = ".")

  rmarkdown::render(system.file("rmd/genes_table.Rmd", package = "pathfindr"),
                    params = list(df = input_processed), output_dir = ".")

  return(final_res)
}
