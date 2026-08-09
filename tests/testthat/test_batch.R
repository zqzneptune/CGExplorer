library(testthat)
library(CGExplorer)

test_that("batch operations and plot_replicate_growth_curves work", {
  reg <- new_registry(project_name = "Batch Test")
  b1 <- new_batch("Batch_1", c("uuid_a", "uuid_b"))
  
  reg <- add_batch(reg, b1)
  expect_equal(length(reg@batches), 1)
  expect_true("Batch_1" %in% names(reg@batches))
  
  reg <- remove_batch(reg, "Batch_1")
  expect_equal(length(reg@batches), 0)
  
  # Test plot_replicate_growth_curves
  dummy_rep_df <- data.frame(
    Time_Hours = c(0, 1, 2, 0, 1, 2),
    OD_Raw = c(0.1, 0.2, 0.4, 0.12, 0.22, 0.45),
    Row = c("A", "A", "A", "A", "A", "A"),
    Column = c(1, 1, 1, 1, 1, 1),
    Gene = c("geneA", "geneA", "geneA", "geneA", "geneA", "geneA"),
    Rep_Suffix = c("P1", "P1", "P1", "P2", "P2", "P2")
  )
  
  plt <- plot_replicate_growth_curves(dummy_rep_df, stain_hr = 2)
  expect_s3_class(plt, "plotly")
})
