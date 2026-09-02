test_that("test Format_BaHZING", {
  # Load data
  data("iHMP_Reduced")
  # Format microbiome data
  formatted_data <- Format_BaHZING(iHMP_Reduced)
  # Test format data
  testthat::expect_equal(object = length(formatted_data), expected = 8)
  testthat::expect_true(is.list(formatted_data))
  testthat::expect_true("Table" %in% names(formatted_data))
  testthat::expect_equal(formatted_data$taxa_levels,
                         c("Phylum", "Class", "Order", "Family", "Genus", "Species"))
  testthat::expect_equal(length(formatted_data$taxon_columns),
                         ncol(formatted_data$Table[[1]]) -
                           ncol(data.frame(phyloseq::sample_data(iHMP_Reduced))))
})

test_that("custom taxa_levels drive hierarchy matrices", {
  data("iHMP_Reduced")

  levels <- c("Phylum", "Class", "Order", "Genus")
  formatted_data <- Format_BaHZING(iHMP_Reduced, taxa_levels = levels)

  testthat::expect_equal(formatted_data$taxa_levels, levels)
  testthat::expect_true(all(c(
    "Class.Phylum.Matrix",
    "Order.Class.Matrix",
    "Genus.Order.Matrix"
  ) %in% names(formatted_data)))
  testthat::expect_false("Species.Genus.Matrix" %in% names(formatted_data))
})

test_that("invalid or unavailable taxa_levels fail clearly", {
  data("iHMP_Reduced")

  testthat::expect_error(
    Format_BaHZING(iHMP_Reduced, taxa_levels = c("Phylum", "Family", "Foo")),
    "distinct first letters"
  )
  testthat::expect_error(
    Format_BaHZING(iHMP_Reduced, taxa_levels = c("Phylum", "Cohort", "Species")),
    "missing required taxa_levels: Cohort"
  )
})

test_that("a nonconsecutive four-level hierarchy is formatted generically", {
  data("iHMP_Reduced")

  custom_levels <- c("Phylum", "Family", "Genus", "Species")
  formatted_data <- Format_BaHZING(iHMP_Reduced, taxa_levels = custom_levels)

  testthat::expect_equal(length(formatted_data), 6)
  testthat::expect_equal(formatted_data$taxa_levels, custom_levels)
  testthat::expect_true(all(c(
    "Family.Phylum.Matrix",
    "Genus.Family.Matrix",
    "Species.Genus.Matrix"
  ) %in% names(formatted_data)))
  testthat::expect_false(any(grepl("Class|Order", names(formatted_data))))
})

test_that("a Kingdom column containing d__ values does not stack prefixes", {
  data("iHMP_Reduced")

  PS <- iHMP_Reduced
  tt <- tax_table(PS)
  colnames(tt)[colnames(tt) == "Domain"] <- "Kingdom"
  tax_table(PS) <- tt

  formatted_data <- Format_BaHZING(PS)

  testthat::expect_true(all(grepl(
    "^k__(?!d__)", formatted_data$taxon_columns, perl = TRUE
  )))
  testthat::expect_false(any(grepl("^k__d__", formatted_data$taxon_columns)))
})



test_that("Error thrown when taxonomic levels are less than 2", {

  data("iHMP_Reduced")

  PS <- iHMP_Reduced
  # Modify the taxonomic table to have only one level
  tax_table(PS) <- tax_table(tax_table(PS)[,1])

  # Expect that the function stops with an error message
  testthat::expect_error(Format_BaHZING(PS), "Need > 1 taxonomic level")
  })


test_that("If species level not present, create species column", {

  data("iHMP_Reduced")

  PS <- iHMP_Reduced
  # phyloseq object without a 'Species' column
  tax_table(PS) <- tax_table(PS)[ , !colnames(tax_table(PS)) %in% "Species"]

  formatted_data <- Format_BaHZING(PS)

  # Extract the taxonomic table from the result
  taxa_table_result <- formatted_data$Table[[1]]

  # Check if 'Species' column is created
  testthat::expect_true(grepl("s__unclassified",
                              colnames(taxa_table_result)[ncol(taxa_table_result)]))
  testthat::expect_true(grepl("s__unclassified",
                              colnames(taxa_table_result)[ncol(taxa_table_result)-10]))
  testthat::expect_true(grepl("s__unclassified",
                              colnames(taxa_table_result)[ncol(taxa_table_result)-20]))
})
