if (!suppressWarnings(require(CGExplorer, quietly = TRUE))) {
  r_files <- list.files("../../R", full.names = TRUE)
  if (length(r_files) > 0) lapply(r_files, source)
}

test_that("batch operations and plot_replicate_growth_curves work", {
  reg <- new_registry(project_name = "Batch Test")
  b1 <- new_batch("Batch_1", c("uuid_a", "uuid_b"))
  
  reg <- add_batch(reg, b1)
  expect_equal(length(reg@batches), 1)
  expect_true("Batch_1" %in% names(reg@batches))
  
  reg <- remove_batch(reg, "Batch_1")
  expect_equal(length(reg@batches), 0)
  
  # Test Plate Usability Flag and Layout Uniformity Validation
  layout1 <- data.frame(Row = c("A", "A"), Column = c(1, 2), Type = c("Blank", "WT"), Gene = c("Blank", "WT"), stringsAsFactors = FALSE)
  layout2 <- data.frame(Row = c("A", "A"), Column = c(1, 2), Type = c("Blank", "Mutant"), Gene = c("Blank", "GeneX"), stringsAsFactors = FALSE)
  
  growth_df <- data.frame(
    DateTime = c(as.POSIXct("2026-08-11 10:00:00"), as.POSIXct("2026-08-11 10:00:00"),
                 as.POSIXct("2026-08-11 11:00:00"), as.POSIXct("2026-08-11 11:00:00")),
    Row = c("A", "A", "A", "A"),
    Column = c(1, 2, 1, 2),
    OD = c(0.1, 0.1, 0.1, 0.5), # Row A Col 1 Blank (OD stays 0.1), Col 2 WT (grows 0.1 -> 0.5)
    stringsAsFactors = FALSE
  )
  
  p1 <- new_plate("S1L1", growth_df, layout1)
  p2 <- new_plate("S1L2", growth_df, layout1)
  p3_diff <- new_plate("S1L3", growth_df, layout2)
  
  p1 <- compute_metrics(p1)
  p2 <- compute_metrics(p2)
  p3_diff <- compute_metrics(p3_diff)
  
  # Default auto QC pass
  expect_true(isTRUE(p1@qc_flags$qc_passed))
  
  # Add to registry
  reg <- add_plate(reg, p1)
  reg <- add_plate(reg, p2)
  reg <- add_plate(reg, p3_diff)
  
  # Valid batch plates (same layout, passed QC)
  expect_true(validate_batch_plates(reg, c(p1@uuid, p2@uuid)))
  
  # Mismatched layout batch creation should error
  expect_error(validate_batch_plates(reg, c(p1@uuid, p3_diff@uuid)), "Layout mismatch")
  
  # FAILED plate quarantine check
  p1@qc_flags$user_override <- "FAILED"
  p1 <- compute_metrics(p1)
  expect_false(p1@qc_flags$qc_passed)
  reg <- update_plate(reg, p1)
  
  expect_error(validate_batch_plates(reg, c(p1@uuid, p2@uuid)), "is marked FAILED")
  
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
