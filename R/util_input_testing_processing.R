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
#'   input = example_pathfindR_input[1:5, ],
#'   pin_name_path = "KEGG",
#'   convert2alias = FALSE
#' )
input_processing <- function(input, p_val_threshold = 0.05, pin_name_path = "Biogrid",
                             convert2alias = TRUE) {
  message("## Processing input. Converting gene symbols,
          if necessary (and if human gene symbols provided)")

  if (!is.logical(convert2alias)) {
    stop("`convert2alias` should be either TRUE or FALSE")
  }

  pin_path <- return_pin_path(pin_name_path)

  if (ncol(input) == 2) {
    input <- data.frame(
      GENE = input[, 1], CHANGE = rep(1e+06, nrow(input)),
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
      "No input p value is lower than the provided threshold (", p_val_threshold,
      ")"
    )
  }
  input <- input[input$P_VALUE <= p_val_threshold, ]
  message("Number of genes in input after p-value filtering: ", nrow(input))

  ## Choose lowest p for each gene
  if (anyDuplicated(input$GENE)) {
    warning("Duplicated genes found! The lowest p value for each gene was selected",
      call. = FALSE
    )

    input <- input[order(input$P_VALUE, decreasing = FALSE), ]
    input <- input[!duplicated(input$GENE), ]
  }

  ## Fix p < 1e-13
  if (any(input$P_VALUE < 1e-13)) {
    message("pathfindR cannot handle p values < 1e-13. These were changed to 1e-13")
    input$P_VALUE <- ifelse(input$P_VALUE < 1e-13, 1e-13, input$P_VALUE)
  }

  ## load and prep pin
  pin <- utils::read.delim(file = pin_path, header = FALSE)

  ## Genes not in pin
  PIN_genes <- c(base::toupper(pin[, 1]), base::toupper(pin[, 3]))
  missing_symbols <- input$GENE[!base::toupper(input$GENE) %in% PIN_genes]
  non_missing_symbols <- input$GENE[base::toupper(input$GENE) %in% PIN_genes]


  if (convert2alias & !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    message(
      "Package 'org.Hs.eg.db' is not installed; returning input genes unchanged.\n",
      "Install it with:\n",
      "  if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')\n",
      "  BiocManager::install('org.Hs.eg.db')"
    )
    convert2alias <- FALSE
  }


  if (convert2alias & length(missing_symbols) != 0) {
    ## use SQL to get alias table and gene_info table (contains the
    ## symbols) first open the database connection
    db_con <- org.Hs.eg.db::org.Hs.eg_dbconn()
    ## the SQL query
    sql_query <- "SELECT * FROM alias, gene_info WHERE alias._id == gene_info._id;"
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
      result <- hsa_alias_df[hsa_alias_df$symbol %in% result$symbol, c(
        "alias_symbol",
        "symbol"
      )]
      result <- result$alias_symbol[base::toupper(result$alias_symbol) %in%
        PIN_genes]
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
    message(paste0("Could not find any interactions for ", n, " (", round(
      perc,
      2
    ), "%) genes in the PIN"))
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
