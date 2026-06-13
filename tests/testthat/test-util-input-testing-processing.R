set.seed(123)

test_that("`input_testing()` -- works as expected", {
  expect_message(
    input_testing(input = example_pathfindR_input, p_val_threshold = 0.05),
    "The input looks OK"
  )

  expect_error(input_testing(input = matrix(), p_val_threshold = 0.05), "the input is not a data frame")

  expect_error(input_testing(
    input = example_pathfindR_input[, 1, drop = FALSE],
    p_val_threshold = 0.05
  ), "the input should have 2 or 3 columns")

  expect_error(
    input_testing(input = example_pathfindR_input[1, ], p_val_threshold = 0.05),
    "There must be at least 2 rows \\(genes\\) in the input data frame"
  )

  expect_error(
    input_testing(input = example_pathfindR_input, p_val_threshold = "INVALID"),
    "`p_val_threshold` must be a numeric value between 0 and 1"
  )

  expect_error(
    input_testing(input = example_pathfindR_input, p_val_threshold = -1),
    "`p_val_threshold` must be between 0 and 1"
  )

  tmp <- example_pathfindR_input
  tmp$adj.P.Val <- NA
  expect_error(input_testing(input = tmp, p_val_threshold = 0.05), "p values cannot contain NA values")

  tmp <- example_pathfindR_input
  tmp$adj.P.Val <- "INVALID"
  expect_error(input_testing(input = tmp, p_val_threshold = 0.05), "p values must all be numeric")

  tmp <- example_pathfindR_input
  tmp$adj.P.Val[1] <- -1
  expect_error(input_testing(input = tmp, p_val_threshold = 0.05), "p values must all be between 0 and 1")
})

test_that("`input_processing()` -- works as expected", {
  input_df <- example_pathfindR_input[1:10, ]
  toy_PIN <- data.frame(
    V1 = sample(example_pathfindR_input$Gene.symbol, 100),
    V2 = "pp", V3 = sample(example_pathfindR_input$Gene.symbol, 100)
  )
  mockery::stub(input_processing, "return_pin_path", NULL)
  mockery::stub(input_processing, "utils::read.delim", toy_PIN)

  expect_is(processed_df <- input_processing(input_df), "data.frame")
  expect_true(ncol(processed_df) == 4)
  expect_true(nrow(processed_df) <= nrow(example_pathfindR_input))

  # no change values provided
  input_df2 <- input_df[, -2]
  expect_is(processed_df2 <- suppressWarnings(input_processing(input_df2)), "data.frame")
  expect_true(ncol(processed_df2) == 4)
  expect_true(all(processed_df2$CHANGE == 1e+06))

  toy_PIN2 <- rbind(toy_PIN, data.frame(
    V1 = c("SERPINA3", "ARHGAP17"), V2 = "pp",
    V3 = c("ACT", "GIG25")
  ))
  mockery::stub(input_processing, "utils::read.delim", toy_PIN2)

  # multiple mapping
  input_multimap <- input_df
  input_multimap$Gene.symbol[1] <- "GIG24"
  input_multimap$Gene.symbol[2] <- "ACT"
  input_multimap$Gene.symbol[3] <- "AACT"
  input_multimap$Gene.symbol[4] <- "GIG25"
  expect_is(processed_df3 <- input_processing(input_multimap), "data.frame")
})

test_that("`input_processing()` -- errors and warnings work", {
  input_df <- example_pathfindR_input[1:10, ]

  toy_PIN <- data.frame(V1 = sample(input_df$Gene.symbol, 7), V2 = "pp", V3 = sample(
    input_df$Gene.symbol,
    7
  ))
  mockery::stub(input_processing, "return_pin_path", NULL)
  mockery::stub(input_processing, "utils::read.delim", toy_PIN)

  input_df$Gene.symbol <- as.factor(input_df$Gene.symbol)
  expect_warning(input_processing(input_df,
    p_val_threshold = 0.05, pin_name_path = "Biogrid",
    convert2alias = TRUE
  ), "The gene column was turned into character from factor.")

  expect_error(input_processing(example_pathfindR_input,
    p_val_threshold = 1e-100,
    pin_name_path = "Biogrid"
  ), "No input p value is lower than the provided threshold \\(1e-100\\)")

  input_dup <- example_pathfindR_input[1:3, ]
  input_dup <- rbind(input_dup, input_dup[1, ])
  expect_warning(
    input_processing(input_dup, p_val_threshold = 0.05, pin_name_path = "Biogrid"),
    "Duplicated genes found! The lowest p value for each gene was selected"
  )

  low_sig_input <- example_pathfindR_input[1:3, ]
  low_sig_input$adj.P.Val <- 1e-15
  expect_message(res <- input_processing(low_sig_input,
    p_val_threshold = 0.05,
    pin_name_path = "Biogrid"
  ), "pathfindR cannot handle p values < 1e-13. These were changed to 1e-13")
  expect_true(all(res$P_VALUE == 1e-13))

  invalid_genes_input <- low_sig_input
  invalid_genes_input$Gene.symbol <- paste0(
    LETTERS[seq_len(nrow(invalid_genes_input))],
    "INVALID"
  )
  expect_error(
    input_processing(invalid_genes_input, p_val_threshold = 0.05, pin_name_path = "Biogrid"),
    "None of the genes were in the PIN\nPlease check your gene symbols"
  )

  low_sig_input$Gene.symbol[1] <- "INVALID_A"
  low_sig_input$Gene.symbol[2] <- "INVALID_B"
  low_sig_input$Gene.symbol[3] <- toy_PIN$V1[1]
  expect_error(
    input_processing(low_sig_input, p_val_threshold = 0.05, pin_name_path = "Biogrid"),
    "After processing, 1 gene \\(or no genes\\) could be mapped to the PIN"
  )

  expect_error(input_processing(low_sig_input,
    p_val_threshold = 0.05, pin_name_path = "Biogrid",
    convert2alias = "INVALID"
  ), "`convert2alias` should be either TRUE or FALSE")
})
