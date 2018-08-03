#' Input Testing
#'
#' @param input the input data that pathfindR uses. The input must be a data
#'   frame with three columns: \enumerate{
#'   \item Gene Symbol (HGNC Gene Symbol)
#'   \item Change value, e.g. log(fold change)
#'   \item adjusted p value associated with test, e.g. differential expression/methylation
#' }
#' @param p_val_threshold the adjusted-p value threshold to use when filtering
#'   the input data frame. Must a numeric value between 0 and 1.
#'
#' @return Only checks if the input and the threshold follows the required
#'   specifications.
#' @export
#' @seealso See \code{\link{run_pathfindR}} for the wrapper function of the
#'   pathfindR workflow
#' @examples
#' input_testing(RA_input, 0.05)
input_testing <- function(input, p_val_threshold){
  if (!is.data.frame(input)) {
    setwd("..")
    stop("the input is not a data frame")
  }

  if (ncol(input) != 3){
    setwd("..")
    stop("There must be exactly 3 columns in the input data frame")
  }

  if (!is.numeric(p_val_threshold)){
    setwd("..")
    stop("`p_val_threshold` must be a numeric value between 0 and 1")
  }

  if (p_val_threshold > 1 | p_val_threshold < 0){
    setwd("..")
    stop("`p_val_threshold` must be between 0 and 1")
  }

  if (!all(is.numeric(input[, 3]))) {
    setwd("..")
    stop("p values, provided in the third column, must all be numeric")
  }

  if (any(input[, 3] > 1 | input[, 3] < 0)) {
    setwd("..")
    stop("p values, provided in the third column, must all be between 0 and 1")
  }

  message("The input looks OK\n\n")
}

#' Process Input
#'
#' @param input the input data that pathfindR uses. The input must be a data
#'   frame with three columns: \enumerate{
#'   \item Gene Symbol (HGNC Gene Symbol)
#'   \item Change value, e.g. log(fold change)
#'   \item adjusted p value associated with test, e.g. differential expression/methylation
#' }
#' @param p_val_threshold the adjusted-p value threshold to use when filtering
#'   the input data frame
#' @param pin_path path to the Protein Interaction Network (PIN) file used in
#'   the analysis
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
#' input_processing(RA_input[1,], 0.05, return_pin_path("KEGG"))
#' }
#' \dontrun{
#' input_processing(RA_input, 0.05, return_pin_path("KEGG"))
#' }
input_processing <- function(input, p_val_threshold, pin_path) {
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

  if (length(missing) != 0) {
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

    ## Give out warning indicating the number of still missing
    n <- sum(converted[, 2] == "NOT_FOUND")
    perc <- n / nrow(input) * 100
    if (sum(converted[, 2] == "NOT_FOUND") != 0)
      message(paste0("Could not find any interactions for ",
                     n,
                     " (", round(perc, 2), "%) genes in the PIN\n\n"))

    ## Convert to appropriate symbol
    input$new_gene <- input$GENE
    input$new_gene[match(converted[, 1], input$new_gene)] <- converted[, 2]

    input <- input[, c(1, 4, 2, 3)]
    colnames(input) <- c("old_GENE", "GENE", "CHANGE", "P_VALUE")

    input <- input[input$GENE != "NOT_FOUND", ]

    if (nrow(input) == 0) {
      setwd("..")
      stop("None of the genes were in the PIN\nPlease check your gene symbols")
    }

    input <- input[order(input$P_VALUE), ]
    input <- input[!duplicated(input$GENE), ]
  } else {
    message(paste0("Found interactions for all genes in the PIN\n\n"))
    input$old_GENE <- input$GENE
    input <- input[, c(4, 1, 2, 3)]
  }
  return(input)
}
