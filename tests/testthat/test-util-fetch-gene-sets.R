test_that("`fetch_gene_sets()` -- can fetch all gene set objects", {
  skip_on_cran()
  for (gset_name in c(
    "KEGG", "mmu_KEGG", "Reactome", "BioCarta", "cell_markers",
    "GO-All", "GO-BP", "GO-CC", "GO-MF"
  )) {
    expect_is(gset_obj <- fetch_gene_sets(
      gene_sets = gset_name, min_gset_size = 10,
      max_gset_size = 300
    ), "list")
    expect_is(gset_obj$genes_by_term, "list")
    expect_is(gset_obj$term_descriptions, "character")
    expect_true(length(gset_obj$genes_by_term) == length(gset_obj$term_descriptions))
    tmp <- vapply(gset_obj$genes_by_term, length, 1L)
    expect_true(min(tmp) >= 10 & max(tmp) <= 300)
  }
  # Custom
  gset_obj <- fetch_gene_sets(
    gene_sets = "Custom", min_gset_size = 20, max_gset_size = 200,
    custom_genes = kegg_genes, custom_descriptions = kegg_descriptions
  )
  expect_is(gset_obj$genes_by_term, "list")
  expect_is(gset_obj$term_descriptions, "character")
  expect_true(length(gset_obj$genes_by_term) == length(gset_obj$term_descriptions))
  tmp <- vapply(gset_obj$genes_by_term, length, 1L)
  expect_true(min(tmp) >= 20 & max(tmp) <= 200)
})


test_that("`fetch_gene_sets()` -- min/max_gset_size args correctly filter gene sets", {
  skip_on_cran()
  min_max_pairs <- list(c(min = 10, max = 300), c(min = 50, max = 200))
  num_of_terms_after_size_filtering <- c()
  for (idx in seq_along(min_max_pairs)) {
    cur_vals <- min_max_pairs[[idx]]
    expect_is(gset_obj <- fetch_gene_sets(
      gene_sets = "KEGG", min_gset_size = cur_vals["min"],
      max_gset_size = cur_vals["max"]
    ), "list")
    sizes_of_terms <- vapply(gset_obj$genes_by_term, length, 1L)
    expect_true(min(sizes_of_terms) >= cur_vals["min"] & max(sizes_of_terms) <=
      cur_vals["max"])
    num_of_terms_after_size_filtering <- c(
      num_of_terms_after_size_filtering,
      length(gset_obj$genes_by_term)
    )
  }

  expect_true(num_of_terms_after_size_filtering[2] < num_of_terms_after_size_filtering[1])
})

test_that("`fetch_gene_sets()` -- for 'Custom' gene set, check if the custom objects are provided", {
  expect_error(fetch_gene_sets(gene_sets = "Custom"), "`custom_genes` and `custom_descriptions` must be provided if `gene_sets = \"Custom\"`")
  expect_error(
    fetch_gene_sets(gene_sets = "Custom", custom_genes = kegg_genes),
    "`custom_genes` and `custom_descriptions` must be provided if `gene_sets = \"Custom\"`"
  )
  expect_error(
    fetch_gene_sets(gene_sets = "Custom", custom_descriptions = kegg_descriptions),
    "`custom_genes` and `custom_descriptions` must be provided if `gene_sets = \"Custom\"`"
  )
})

test_that("`fetch_gene_sets()` -- argument checks work", {
  all_gs_opts <- c(
    "KEGG", "Reactome", "BioCarta", "GO-All", "GO-BP", "GO-CC",
    "GO-MF", "cell_markers", "mmu_KEGG", "Custom"
  )
  expect_error(fetch_gene_sets(gene_sets = "INVALID"), paste0(
    "`gene_sets` should be one of ",
    paste(dQuote(all_gs_opts), collapse = ", ")
  ))

  expect_error(fetch_gene_sets(min_gset_size = "INVALID"), "`min_gset_size` should be numeric")

  expect_error(fetch_gene_sets(max_gset_size = "INVALID"), "`max_gset_size` should be numeric")

  expect_error(
    fetch_gene_sets(gene_sets = "Custom", custom_genes = "INVALID", custom_descriptions = ""),
    "`custom_genes` should be a list of term gene sets"
  )
  expect_error(
    fetch_gene_sets(gene_sets = "Custom", custom_genes = list(), custom_descriptions = ""),
    "`custom_genes` should be a named list \\(names are gene set IDs\\)"
  )

  expect_error(fetch_gene_sets(
    gene_sets = "Custom", custom_genes = kegg_genes,
    custom_descriptions = list()
  ), "`custom_descriptions` should be a vector of term gene descriptions")
  expect_error(fetch_gene_sets(
    gene_sets = "Custom", custom_genes = kegg_genes,
    custom_descriptions = 1:3
  ), "`custom_descriptions` should be a named vector \\(names are gene set IDs\\)")
})
