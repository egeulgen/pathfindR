#' Input Tesing
#'
#' @param input the input data that pathfindr uses. The input must be a data
#'   frame with three columns: 1. Gene Symbol (HGNC Gene Symbol) 2. Change
#'   value, e.g. log(fold change) 3. adjusted p value associated with test, e.g.
#'   differential expression/methylation
#' @param p_val_threshold the adjusted-p value threshold to use when filtering
#'   the input data frame. Must a numeric value between 0 and 1.
#'
#' @return Only checks if the input and the threshold follows the required
#'   specifications
#' @export
#'
#' @examples
#' input_testing()
input_testing <- function(input, p_val_threshold){
  if (!is.data.frame(input))
    stop("the input is not a data frame")
  if (ncol(input) > 3)
    stop("there are more than 3 columns in the input data frame")
  if (ncol(input) < 3)
    stop("there are less than 3 columns in the input data frame")

  if (!is.numeric(p_val_threshold))
    stop("p_val_threshold must be a numeric value between 0 and 1")
  if (p_val_threshold > 1 | p_val_threshold < 0)
    stop("p_val_threshold must be between 0 and 1")

  if (!all(is.numeric(input[, 3])))
    stop("p values, provided in the third column, must all be numeric")
  if (any(input[, 3] > 1 | input[, 3] < 0))
    stop("p values, provided in the third column, must all be between 0 and 1")

  if (anyDuplicated(input[, 1]))
    stop("Duplicated genes found, please choose ones with lowest p-values among replicates")

  cat("The input looks OK\n\n")
}

#' Process Input
#'
#' @param input the input data that pathfindr uses
#' @param p_val_threshold the adjusted-p value threshold to use when filtering
#'   the input data frame
#' @param ppi_path path to the PPI network file used in the analysis
#'
#' @return This function first filters the input so that all p values are less
#'   than or equal to the threshold. Next, gene symbols that are not found in
#'   the PPI are identified. If aliases of these gene symbols are found in the
#'   PPI network, the symbols are converted to the corresponding aliases. The
#'   resulting data frame containing the original gene symbols, the updated
#'   symbols, change values and p values is then returned.
#' @export
#'
#' @examples
#' input_processing(input, p_val_threshold = 0.05, "path/to/ppi")
input_processing <- function(input, p_val_threshold, ppi_path) {
  colnames(input) <- c("GENE", "CHANGE", "P_VALUE")

  ## Discard larger than p-value threshold
  input <- input[input[, 2] <= p_val_threshold, ]

  ## load and prep ppi
  ppi <- read.delim(file = ppi_path, header = F, stringsAsFactors = F)
  ppi <- apply(ppi, 1, function(x) strsplit(x, " pp "))
  ppi <- lapply(ppi, function(x) x$V1)
  ppi <- ppi[sapply(ppi, length) == 2]
  ppi <- ppi[sapply(ppi, function(x) x[1] != x[2])] ## removing loops
  ppi <- Reduce(rbind, ppi)

  ## Genes not in ppi
  missing <- input$GENE[!input$GENE %in% c(ppi[, 1], ppi[, 2])]

  ## use sql to get alias table and gene_info table (contains the symbols)
  ## first open the database connection
  db_con <- org.Hs.eg.db::org.Hs.eg_dbconn()
  ## write the SQL query
  sql_query <-
    "SELECT * FROM alias, gene_info WHERE alias._id == gene_info._id;"
  ## execute the query on the database
  alias_symbol <- DBI::dbGetQuery(db_con, sql_query)

  ## loop for getting all symbols
  converted <- c()
  for (i in 1:length(missing)) {
    result <- alias_symbol[alias_symbol$alias_symbol == missing[i], c(2, 5)]
    result <- alias_symbol[alias_symbol$symbol %in% result$symbol, c(2, 5)]
    result <- result$alias_symbol[result$alias_symbol %in%
                                    c(ppi[, 1], ppi[, 2])]

    converted <- rbind(converted, c(missing[i],
                                    ifelse(length(result) == 0,
                                           "NOT_FOUND", result[1])))
  }

  ## Give out warning indicating the number of still missing
  perc <- sum(converted[, 2] == "NOT_FOUND") / nrow(input) * 100
  if (sum(converted[, 2] == "NOT_FOUND") != 0)
    cat(paste0("Could not find ",
               sum(converted[, 2] == "NOT_FOUND"),
               " (", round(perc), "%) genes in the PPI network\n\n"))

  ## Convert to appropriate symbol
  converted <- converted[converted[, 2] != "NOT_FOUND", ]
  input$UpdatedGene <- input$GENE
  input$UpdatedGene[match(converted[, 1], input$UpdatedGene)] <- converted[, 2]

  input <- input[, c(1, 4, 2, 3)]
  colnames(input) <- c("old_GENE", "GENE", "CHANGE", "SPOTPvalue")

  return(input)
}
