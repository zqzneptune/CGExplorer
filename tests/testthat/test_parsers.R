library(testthat)
library(CGExplorer)

test_that("parse_plate_layout auto-detects 96 and 384 well formats", {
  # Create a mock 96-well layout
  layout_96 <- data.frame(
    Row = LETTERS[1:8]
  )
  for (i in 1:12) {
    layout_96[[as.character(i)]] <- "WT"
  }
  
  tmp_96 <- tempfile(fileext = ".tsv")
  readr::write_tsv(layout_96, tmp_96)
  
  res_96 <- parse_plate_layout(tmp_96)
  expect_equal(res_96$format, "96")
  expect_equal(res_96$n_wells, 96L)
  expect_equal(res_96$n_rows, 8L)
  expect_equal(res_96$n_cols, 12L)
  expect_equal(nrow(res_96$data), 96L)
  
  unlink(tmp_96)
  
  # Create a mock 384-well layout
  layout_384 <- data.frame(
    Row = LETTERS[1:16]
  )
  for (i in 1:24) {
    layout_384[[as.character(i)]] <- "Mutant"
  }
  
  tmp_384 <- tempfile(fileext = ".tsv")
  readr::write_tsv(layout_384, tmp_384)
  
  res_384 <- parse_plate_layout(tmp_384)
  expect_equal(res_384$format, "384")
  expect_equal(res_384$n_wells, 384L)
  expect_equal(nrow(res_384$data), 384L)
  
  unlink(tmp_384)
})

test_that("parse_plate_layout fails on invalid sizes", {
  layout_invalid <- data.frame(
    Row = LETTERS[1:10]
  )
  for (i in 1:15) {
    layout_invalid[[as.character(i)]] <- "Mutant"
  }
  
  tmp_inv <- tempfile(fileext = ".tsv")
  readr::write_tsv(layout_invalid, tmp_inv)
  
  expect_error(parse_plate_layout(tmp_inv), "Expected 96 or 384 wells; found 150")
  
  unlink(tmp_inv)
})
