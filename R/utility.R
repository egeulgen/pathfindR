#' Wrapper for Active Subnetwork Search + Enrichment over Single/Multiple Iteration(s)
#'
#' @param input_processed processed input data frame
#' @param pin_path path/to/PIN/file
#' @param gset_list list for gene sets
#' @inheritParams run_pathfindR
#' @inheritParams active_snw_search
#' @inheritParams enrichment_analyses
#' @param iterations number of iterations for active subnetwork search and
#'  enrichment analyses (Default = 10)
#' @param n_processes optional argument for specifying the number of processes
#'  used by foreach. If not specified, the function determines this
#'  automatically (Default == NULL. Gets set to 1 for Genetic Algorithm)
#'
#' @return Data frame of combined pathfindR enrichment results
active_snw_enrichment_wrapper <- function(input_processed,
                                          pin_path,
                                          gset_list,
                                          enrichment_threshold,
                                          list_active_snw_genes,
                                          adj_method = "bonferroni",
                                          search_method = "GR",
                                          use_all_positives = FALSE,
                                          iterations = 10, n_processes = NULL,
                                          score_quan_thr = 0.8, sig_gene_thr = 0.02,
                                          saTemp0 = 1, saTemp1 = 0.01, saIter = 10000,
                                          gaPop = 400, gaIter = 200, gaThread = 5,
                                          gaCrossover = 1, gaMut = 0,
                                          grMaxDepth = 1, grSearchDepth = 1,
                                          grOverlap = 0.5, grSubNum = 1000,
                                          silent_option = TRUE) {
  message("## Performing Active Subnetwork Search and Enrichment")
  ############ Argument checks
  # Active Subnetwork Search Method
  valid_mets <- c("GR", "SA", "GA")
  if (!search_method %in% valid_mets) {
    stop(
      "`search_method` should be one of ",
      paste(dQuote(valid_mets), collapse = ", ")
    )
  }

  ## If search_method is GA, set iterations as 1
  if (search_method == "GA") {
    warning("`iterations` is set to 1 because `search_method = \"GA\"`",
      call. = FALSE
    )
    iterations <- 1
  }

  if (!is.null(n_processes)) {
    if (!is.numeric(n_processes)) {
      stop("`n_processes` should be either NULL or a positive integer")
    }
    if (n_processes < 1) {
      stop("`n_processes` should be > 1")
    }
  }

  # calculate the number of processes, if necessary
  if (is.null(n_processes)) {
    n_processes <- parallel::detectCores() - 1
  }

  ## If iterations < n_processes, set n_processes to iterations
  if (iterations < n_processes & iterations != 1) {
    message("`n_processes` is set to `iterations` because `iterations` < `n_processes`")
    n_processes <- iterations
  }

  if (!is.logical(use_all_positives)) {
    stop("`use_all_positives` should be either TRUE or FALSE")
  }

  if (!is.logical(silent_option)) {
    stop("`silent_option` should be either TRUE or FALSE")
  }

  if (!is.numeric(iterations)) {
    stop("`iterations` should be a positive integer")
  }
  if (iterations < 1) {
    stop("`iterations` should be >= 1")
  }

  geneInitProbs <- 0.1
  dirs <- c()
  if (iterations > 1) {
    geneInitProbs <- seq(from = 0.01, to = 0.2, length.out = iterations)

    for (i in base::seq_len(iterations)) {
      dir_i <- file.path("active_snw_searches", paste0("Iteration_", i))
      dir.create(dir_i, recursive = TRUE)
      dirs <- c(dirs, dir_i)
    }
  }

  single_iter_wrapper <- function(i = NULL) {
    snws_file <- "active_snws"
    dir_for_parallel_run <- NULL
    if (!is.null(i)) {
      snws_file <- paste0("active_snws_", i)
      dir_for_parallel_run <- dirs[i]
    }
    snws <- active_snw_search(
      input_for_search = input_processed,
      pin_name_path = pin_path,
      snws_file = snws_file,
      dir_for_parallel_run = dir_for_parallel_run,
      score_quan_thr = score_quan_thr,
      sig_gene_thr = sig_gene_thr,
      search_method = search_method,
      seedForRandom = ifelse(is.null(i), 1234, i),
      silent_option = silent_option,
      use_all_positives = use_all_positives,
      geneInitProbs = ifelse(!is.null(i), geneInitProbs[i], geneInitProbs),
      saTemp0 = saTemp0, saTemp1 = saTemp1, saIter = saIter,
      gaPop = gaPop, gaIter = gaIter,
      gaThread = gaThread,
      gaCrossover = gaCrossover, gaMut = gaMut,
      grMaxDepth = grMaxDepth, grSearchDepth = grSearchDepth,
      grOverlap = grOverlap, grSubNum = grSubNum
    )
    enrichment_res <- enrichment_analyses(
      snws = snws,
      sig_genes_vec = input_processed$GENE,
      pin_name_path = pin_path,
      genes_by_term = gset_list$genes_by_term,
      term_descriptions = gset_list$term_descriptions,
      adj_method = adj_method,
      enrichment_threshold = enrichment_threshold,
      list_active_snw_genes = list_active_snw_genes
    )
    return(enrichment_res)
  }

  if (iterations == 1) {
    combined_res <- single_iter_wrapper()
  } else {
    cl <- parallel::makeCluster(n_processes, setup_strategy = "sequential")
    doParallel::registerDoParallel(cl)
    `%dopar%` <- foreach::`%dopar%`
    combined_res <- foreach::foreach(i = 1:iterations, .combine = rbind) %dopar% {
      single_iter_wrapper(i)
    }
    parallel::stopCluster(cl)
  }
  return(combined_res)
}



#' Configure Output Directory Name
#'
#' @inheritParams run_pathfindR
#'
#' @return /path/to/output/dir
configure_output_dir <- function(output_dir = NULL) {
  output_dir_init <- output_dir
  output_dir <- ifelse(
    is.null(output_dir),
    file.path(tempdir(check = TRUE), "pathfindR_results"),
    output_dir
  )
  dir_changed <- FALSE
  while (dir.exists(output_dir)) {
    output_dir <- sub("/$", "", output_dir)
    if (grepl("\\(\\d+\\)$", output_dir)) {
      output_dir <- unlist(strsplit(output_dir, "\\("))
      suffix <- as.numeric(sub("\\)", "", output_dir[2])) + 1
      output_dir <- paste0(output_dir[1], "(", suffix, ")")
    } else {
      output_dir <- paste0(output_dir, "(1)")
    }
    dir_changed <- TRUE
  }

  if (dir_changed & !is.null(output_dir_init)) {
    message(paste0(
      "There is already a directory named \"", output_dir_init,
      "\".\nWriting the result to \"", output_dir,
      "\" not to overwrite any previous results."
    ))
  }
  return(output_dir)
}

#' Create HTML Report of pathfindR Results
#'
#' @inheritParams run_pathfindR
#' @param input_processed processed input data frame
#' @param final_res final pathfindR result data frame
create_HTML_report <- function(input, input_processed, final_res) {
  message("## Creating HTML report")
  rmarkdown::render(
    input = system.file("rmd", "results.Rmd",
      package = "pathfindR"
    ),
    output_dir = "."
  )
  rmarkdown::render(
    input = system.file("rmd", "enriched_terms.Rmd",
      package = "pathfindR"
    ),
    params = list(df = final_res),
    output_dir = "."
  )
  rmarkdown::render(
    input = system.file("rmd", "conversion_table.Rmd",
      package = "pathfindR"
    ),
    params = list(
      df = input_processed,
      original_df = input
    ),
    output_dir = "."
  )
}

#' Input Testing
#'
#' @param input the input data that pathfindR uses. The input must be a data
#'   frame with three columns: \enumerate{
#'   \item Gene Symbol (Gene Symbol)
#'   \item Change value, e.g. log(fold change) (OPTIONAL)
#'   \item p value, e.g. adjusted p value associated with differential expression
#' }
#' @param p_val_threshold the p value threshold to use when filtering
#'   the input data frame. Must a numeric value between 0 and 1. (default = 0.05)
#'
#' @return Only checks if the input and the threshold follows the required
#'   specifications.
#' @export
#' @seealso See \code{\link{run_pathfindR}} for the wrapper function of the
#'   pathfindR workflow
#' @examples
#' input_testing(example_pathfindR_input, 0.05)
input_testing <- function(input, p_val_threshold = 0.05) {
  message("## Testing input")

  if (!is.data.frame(input)) {
    stop("the input is not a data frame")
  }

  if (ncol(input) != 2 & ncol(input) != 3) {
    stop("the input should have 2 or 3 columns")
  }

  if (nrow(input) < 2) {
    stop("There must be at least 2 rows (genes) in the input data frame")
  }

  if (!is.numeric(p_val_threshold)) {
    stop("`p_val_threshold` must be a numeric value between 0 and 1")
  }

  if (p_val_threshold > 1 | p_val_threshold < 0) {
    stop("`p_val_threshold` must be between 0 and 1")
  }

  # if changes are provided, p vals are in col. 3, else in col. 2
  p_column <- ifelse(ncol(input) == 3, 3, 2)

  if (any(is.na(input[, p_column]))) {
    stop("p values cannot contain NA values")
  }

  if (!all(is.numeric(input[, p_column]))) {
    stop("p values must all be numeric")
  }

  if (any(input[, p_column] > 1 | input[, p_column] < 0)) {
    stop("p values must all be between 0 and 1")
  }

  message("The input looks OK")
}

#' Process Input
#' @inheritParams input_testing
#' @inheritParams active_snw_search
#' @inheritParams return_pin_path
#' @param convert2alias boolean to indicate whether or not to convert gene symbols
#' in the input that are not found in the PIN to an alias symbol found in the PIN
#' (default = TRUE) IMPORTANT NOTE: the conversion uses human gene symbols/alias symbols.
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
#' processed_df <- input_processing(
#'   input = example_pathfindR_input[1:5, ],
#'   pin_name_path = "KEGG"
#' )
#' processed_df <- input_processing(
#'   input = example_pathfindR_input[1:10, ],
#'   pin_name_path = "KEGG",
#'   convert2alias = FALSE
#' )
input_processing <- function(input, p_val_threshold = 0.05,
                             pin_name_path = "Biogrid", convert2alias = TRUE) {
  message("## Processing input. Converting gene symbols,
          if necessary (and if human gene symbols provided)")

  if (!is.logical(convert2alias)) {
    stop("`convert2alias` should be either TRUE or FALSE")
  }

  pin_path <- return_pin_path(pin_name_path)

  if (ncol(input) == 2) {
    input <- data.frame(
      GENE = input[, 1],
      CHANGE = rep(1e6, nrow(input)),
      P_VALUE = input[, 2]
    )
  }

  colnames(input) <- c("GENE", "CHANGE", "P_VALUE")

  ## Turn GENE into character
  if (is.factor(input$GENE)) {
    warning("The gene column was turned into character from factor.", call. = FALSE)
    input$GENE <- as.character(input$GENE)
  }

  message("Number of genes provided in input: ", nrow(input))
  ## Discard larger than p-value threshold
  if (sum(input$P_VALUE <= p_val_threshold) == 0) {
    stop(
      "No input p value is lower than the provided threshold (",
      p_val_threshold, ")"
    )
  }
  input <- input[input$P_VALUE <= p_val_threshold, ]
  message("Number of genes in input after p-value filtering: ", nrow(input))

  ## Choose lowest p for each gene
  if (anyDuplicated(input$GENE)) {
    warning("Duplicated genes found! The lowest p value for each gene was selected", call. = FALSE)

    input <- input[order(input$P_VALUE, decreasing = FALSE), ]
    input <- input[!duplicated(input$GENE), ]
  }

  ## Fix p < 1e-13
  if (any(input$P_VALUE < 1e-13)) {
    message("pathfindR cannot handle p values < 1e-13. These were changed to 1e-13")
    input$P_VALUE <- ifelse(input$P_VALUE < 1e-13, 1e-13, input$P_VALUE)
  }

  ## load and prep pin
  pin <- utils::read.delim(
    file = pin_path,
    header = FALSE, stringsAsFactors = FALSE
  )

  ## Genes not in pin
  PIN_genes <- c(base::toupper(pin[, 1]), base::toupper(pin[, 3]))
  missing_symbols <- input$GENE[!base::toupper(input$GENE) %in% PIN_genes]
  non_missing_symbols <- input$GENE[base::toupper(input$GENE) %in% PIN_genes]

  if (convert2alias & length(missing_symbols) != 0) {
    ## use SQL to get alias table and gene_info table (contains the symbols)
    ## first open the database connection
    db_con <- org.Hs.eg.db::org.Hs.eg_dbconn()
    ## the SQL query
    sql_query <-
      "SELECT * FROM alias, gene_info WHERE alias._id == gene_info._id;"
    ## execute the query on the database
    hsa_alias_df <- DBI::dbGetQuery(db_con, sql_query)

    select_alias <- function(result, converted, idx) {
      while (idx > 0) {
        if (!result[idx] %in% c(converted[, 2], non_missing_symbols)) {
          return(result[idx])
        }
        idx <- idx - 1
      }
      return("NOT_FOUND")
    }

    ## loop for getting all symbols
    converted <- c()
    for (i in base::seq_len(length(missing_symbols))) {
      result <- hsa_alias_df[
        hsa_alias_df$alias_symbol == missing_symbols[i],
        c("alias_symbol", "symbol")
      ]
      result <- hsa_alias_df[
        hsa_alias_df$symbol %in% result$symbol,
        c("alias_symbol", "symbol")
      ]
      result <- result$alias_symbol[base::toupper(result$alias_symbol) %in% PIN_genes]
      ## avoid duplicate entries
      to_add <- select_alias(result, converted, length(result))
      converted <- rbind(converted, c(missing_symbols[i], to_add))
    }

    ## Convert to appropriate symbol
    input$new_gene <- input$GENE
    input$new_gene[match(converted[, 1], input$new_gene)] <- converted[, 2]
  } else {
    input$new_gene <- ifelse(input$GENE %in% missing_symbols, "NOT_FOUND", input$GENE)
  }

  ## number and percent still missing
  n <- sum(input$new_gene == "NOT_FOUND")
  perc <- n / nrow(input) * 100

  if (n == nrow(input)) {
    stop("None of the genes were in the PIN\nPlease check your gene symbols")
  }

  ## Give out warning indicating the number of still missing
  if (n != 0) {
    message(paste0(
      "Could not find any interactions for ",
      n, " (", round(perc, 2), "%) genes in the PIN"
    ))
  } else {
    message(paste0("Found interactions for all genes in the PIN"))
  }

  ## reorder columns
  input <- input[, c(1, 4, 2, 3)]
  colnames(input) <- c("old_GENE", "GENE", "CHANGE", "P_VALUE")

  input <- input[input$GENE != "NOT_FOUND", ]

  ## Keep lowest p value for duplicated genes
  input <- input[order(input$P_VALUE), ]
  input <- input[!duplicated(input$GENE), ]

  ## Check that at least two genes remain
  if (nrow(input) < 2) {
    stop("After processing, 1 gene (or no genes) could be mapped to the PIN")
  }

  message("Final number of genes in input: ", nrow(input))

  return(input)
}

#' Annotate the Affected Genes in the Provided Enriched Terms
#'
#' Function to annotate the involved affected (input) genes in each term.
#'
#' @param result_df data frame of enrichment results.
#'  The only must-have column is "ID".
#' @param input_processed input data processed via \code{\link{input_processing}}
#' @param genes_by_term List that contains genes for each gene set. Names of
#'   this list are gene set IDs (default = kegg_genes)
#'
#' @return The original data frame with two additional columns:  \describe{
#'   \item{Up_regulated}{the up-regulated genes in the input involved in the given term's gene set, comma-separated}
#'   \item{Down_regulated}{the down-regulated genes in the input involved in the given term's gene set, comma-separated}
#' }
#' @export
#'
#' @examples
#' example_gene_data <- example_pathfindR_input
#' colnames(example_gene_data) <- c("GENE", "CHANGE", "P_VALUE")
#'
#' annotated_result <- annotate_term_genes(
#'   result_df = example_pathfindR_output,
#'   input_processed = example_gene_data
#' )
annotate_term_genes <- function(result_df,
                                input_processed,
                                genes_by_term = pathfindR.data::kegg_genes) {
  message("## Annotating involved genes and visualizing enriched terms")
  ### Argument checks
  if (!is.data.frame(result_df)) {
    stop("`result_df` should be a data frame")
  }
  if (!"ID" %in% colnames(result_df)) {
    stop("`result_df` should contain an \"ID\" column")
  }

  if (!is.data.frame(input_processed)) {
    stop("`input_processed` should be a data frame")
  }
  if (!all(c("GENE", "CHANGE") %in% colnames(input_processed))) {
    stop("`input_processed` should contain the columns \"GENE\" and \"CHANGE\"")
  }

  if (!is.list(genes_by_term)) {
    stop("`genes_by_term` should be a list of term gene sets")
  }
  if (is.null(names(genes_by_term))) {
    stop("`genes_by_term` should be a named list (names are gene set IDs)")
  }

  ### Annotate up/down-regulated term-related genes
  ## Up/Down-regulated genes
  upreg <- base::toupper(input_processed$GENE[input_processed$CHANGE >= 0])
  downreg <- base::toupper(input_processed$GENE[input_processed$CHANGE < 0])

  ## Annotation
  annotated_df <- result_df
  annotated_df$Down_regulated <- annotated_df$Up_regulated <- NA
  for (i in base::seq_len(nrow(annotated_df))) {
    idx <- which(names(genes_by_term) == annotated_df$ID[i])
    temp <- genes_by_term[[idx]]
    annotated_df$Up_regulated[i] <- paste(temp[base::toupper(temp) %in% upreg],
      collapse = ", "
    )
    annotated_df$Down_regulated[i] <- paste(temp[base::toupper(temp) %in% downreg],
      collapse = ", "
    )
  }

  return(annotated_df)
}


#' Fetch Gene Set Objects
#'
#' Function for obtaining the gene sets per term and the term descriptions to
#' be used for enrichment analysis.
#'
#' @param gene_sets Name of the gene sets to be used for enrichment analysis.
#'  Available gene sets are "KEGG", "Reactome", "BioCarta", "GO-All",
#'  "GO-BP", "GO-CC", "GO-MF", "cell_markers", "mmu_KEGG" or "Custom".
#'  If "Custom", the arguments \code{custom_genes} and \code{custom_descriptions}
#'  must be specified. (Default = "KEGG")
#' @param min_gset_size minimum number of genes a term must contain (default = 10)
#' @param max_gset_size maximum number of genes a term must contain (default = 300)
#' @param custom_genes a list containing the genes involved in each custom
#'  term. Each element is a vector of gene symbols located in the given custom
#'  term. Names should correspond to the IDs of the custom terms.
#' @param custom_descriptions A vector containing the descriptions for each
#'  custom  term. Names of the vector should correspond to the IDs of the custom
#'  terms.
#'
#' @return a list containing 2 elements \describe{
#'   \item{genes_by_term}{list of vectors of genes contained in each term}
#'   \item{term_descriptions}{vector of descriptions per each term}
#' }
#'
#' @export
#'
#' @examples
#' KEGG_gset <- fetch_gene_set()
#' GO_MF_gset <- fetch_gene_set("GO-MF", min_gset_size = 20, max_gset_size = 100)
fetch_gene_set <- function(gene_sets = "KEGG",
                           min_gset_size = 10,
                           max_gset_size = 300,
                           custom_genes = NULL,
                           custom_descriptions = NULL) {
  ### Argument checks
  all_gs_opts <- c(
    "KEGG", "Reactome", "BioCarta",
    "GO-All", "GO-BP", "GO-CC", "GO-MF",
    "cell_markers",
    "mmu_KEGG", "Custom"
  )
  if (!gene_sets %in% all_gs_opts) {
    stop("`gene_sets` should be one of ", paste(dQuote(all_gs_opts), collapse = ", "))
  }

  if (!is.numeric(min_gset_size)) {
    stop("`min_gset_size` should be numeric")
  }
  if (!is.numeric(max_gset_size)) {
    stop("`max_gset_size` should be numeric")
  }


  ### Custom Gene Sets
  if (gene_sets == "Custom") {
    if (is.null(custom_genes) | is.null(custom_descriptions)) {
      stop("`custom_genes` and `custom_descriptions` must be provided if `gene_sets = \"Custom\"`")
    }

    if (!is.list(custom_genes)) {
      stop("`custom_genes` should be a list of term gene sets")
    }
    if (is.null(names(custom_genes))) {
      stop("`custom_genes` should be a named list (names are gene set IDs)")
    }

    if (!is.atomic(custom_descriptions)) {
      stop("`custom_descriptions` should be a vector of term gene descriptions")
    }
    if (is.null(names(custom_descriptions))) {
      stop("`custom_descriptions` should be a named vector (names are gene set IDs)")
    }

    # filter by size
    gset_lens <- vapply(custom_genes, length, 1)
    keep <- which(gset_lens >= min_gset_size & gset_lens <= max_gset_size)
    custom_genes <- custom_genes[keep]
    custom_descriptions <- custom_descriptions[names(custom_genes)]

    return(list(
      genes_by_term = custom_genes,
      term_descriptions = custom_descriptions
    ))
  }

  ### Built-in Gene Sets
  ## GO gene sets
  if (grepl("^GO", gene_sets)) {
    genes_by_term <- pathfindR.data::go_all_genes

    GO_df <- pathfindR.data:::GO_all_terms_df
    term_descriptions <- GO_df$GO_term
    names(term_descriptions) <- GO_df$GO_ID

    if (gene_sets == "GO-BP") {
      tmp <- GO_df$GO_ID[GO_df$Category == "Process"]
      genes_by_term <- genes_by_term[tmp]
      term_descriptions <- term_descriptions[tmp]
    } else if (gene_sets == "GO-CC") {
      tmp <- GO_df$GO_ID[GO_df$Category == "Component"]
      genes_by_term <- genes_by_term[tmp]
      term_descriptions <- term_descriptions[tmp]
    } else if (gene_sets == "GO-MF") {
      tmp <- GO_df$GO_ID[GO_df$Category == "Function"]
      genes_by_term <- genes_by_term[tmp]
      term_descriptions <- term_descriptions[tmp]
    }

    ## non-GO (KEGG, Reactome, BioCarta, mmu_KEGG)
  } else {
    if (gene_sets == "KEGG") {
      genes_by_term <- pathfindR.data::kegg_genes
      term_descriptions <- pathfindR.data::kegg_descriptions
    } else if (gene_sets == "Reactome") {
      genes_by_term <- pathfindR.data::reactome_genes
      term_descriptions <- pathfindR.data::reactome_descriptions
    } else if (gene_sets == "BioCarta") {
      genes_by_term <- pathfindR.data::biocarta_genes
      term_descriptions <- pathfindR.data::biocarta_descriptions
    } else if (gene_sets == "mmu_KEGG") {
      genes_by_term <- pathfindR.data::mmu_kegg_genes
      term_descriptions <- pathfindR.data::mmu_kegg_descriptions
    } else {
      genes_by_term <- pathfindR.data::cell_markers_gsets
      term_descriptions <- pathfindR.data::cell_markers_descriptions
    }
  }

  # filter by size
  term_lens <- vapply(genes_by_term, length, 1)
  keep <- which(term_lens >= min_gset_size & term_lens <= max_gset_size)
  genes_by_term <- genes_by_term[keep]
  term_descriptions <- term_descriptions[names(genes_by_term)]

  return(list(
    genes_by_term = genes_by_term,
    term_descriptions = term_descriptions
  ))
}

#' Return The Path to Given Protein-Protein Interaction Network (PIN)
#'
#' This function returns the absolute path/to/PIN.sif. While the default PINs are
#' "Biogrid", "STRING", "GeneMania", "IntAct", "KEGG" and "mmu_STRING". The user can also
#' use any other PIN by specifying the "path/to/PIN.sif". All PINs to be used
#' in this package must formatted as SIF files: i.e. have 3 columns with no
#' header, no row names and be tab-separated. Columns 1 and 3 must be
#' interactors' gene symbols, column 2 must be a column with all
#' rows consisting of "pp".
#'
#' @param pin_name_path Name of the chosen PIN or absolute/path/to/PIN.sif. If PIN name,
#'   must be one of c("Biogrid", "STRING", "GeneMania", "IntAct", "KEGG", "mmu_STRING"). If
#'   path/to/PIN.sif, the file must comply with the PIN specifications. (Default = "Biogrid")
#'
#' @return The absolute path to chosen PIN.
#'
#' @export
#' @seealso See \code{\link{run_pathfindR}} for the wrapper function of the
#'   pathfindR workflow
#' @examples
#' \dontrun{
#' pin_path <- return_pin_path("GeneMania")
#' }
return_pin_path <- function(pin_name_path = "Biogrid") {
  ## Default PINs
  valid_opts <- c(
    "Biogrid", "STRING", "GeneMania", "IntAct", "KEGG",
    "mmu_STRING", "/path/to/custom/SIF"
  )
  if (pin_name_path %in% valid_opts[-length(valid_opts)]) {
    path <- file.path(tempdir(check = TRUE), paste0(pin_name_path, ".sif"))
    if (!file.exists(path)) {
      adj_list <- utils::getFromNamespace(paste0(tolower(pin_name_path), "_adj_list"),
        ns = "pathfindR.data"
      )

      pin_df <- lapply(seq_along(adj_list),
        function(i, nm, val) {
          data.frame(base::toupper(nm[[i]]),
            "pp",
            base::toupper(val[[i]]),
            stringsAsFactors = FALSE
          )
        },
        val = adj_list, nm = names(adj_list)
      )
      pin_df <- base::do.call("rbind", pin_df)
      utils::write.table(pin_df,
        path,
        sep = "\t",
        row.names = FALSE, col.names = FALSE, quote = FALSE
      )
    }
    path <- normalizePath(path)

    ## Custom PIN
  } else if (file.exists(suppressWarnings(normalizePath(pin_name_path)))) {
    path <- normalizePath(pin_name_path)
    pin <- utils::read.delim(
      file = path, quote = "",
      header = FALSE, stringsAsFactors = FALSE
    )
    if (ncol(pin) != 3) {
      stop("The PIN file must have 3 columns and be tab-separated")
    }

    if (any(pin[, 2] != "pp")) {
      stop("The second column of the PIN file must all be \"pp\" ")
    }

    if (any(grepl("[a-z]", pin[, 1])) | any(grepl("[a-z]", pin[, 3]))) {
      pin[, 1] <- base::toupper(pin[, 1])
      pin[, 3] <- base::toupper(pin[, 3])

      path <- file.path(tempdir(check = TRUE), "custom_PIN.sif")
      utils::write.table(pin,
        path,
        sep = "\t",
        row.names = FALSE, col.names = FALSE, quote = FALSE
      )
      path <- normalizePath(path)
    }
  } else {
    stop(
      "The chosen PIN must be one of:\n",
      paste(dQuote(valid_opts), collapse = ", ")
    )
  }
  return(path)
}
