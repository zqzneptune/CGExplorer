library(testthat)
library(CGExplorer)

test_that("app launches and S4 classes work", {
  # Test PlateRegistry constructor
  reg <- new_registry(project_name = "Test Project")
  expect_s4_class(reg, "PlateRegistry")
  expect_equal(reg@project_name, "Test Project")
  
  # Test Batch constructor
  batch <- new_batch("Test Batch", c("uuid1", "uuid2"))
  expect_s4_class(batch, "Batch")
  expect_equal(batch@name, "Test Batch")
  
  # Test Plate constructor
  dummy_growth <- data.frame(
    DateTime = c(Sys.time(), Sys.time() + 3600),
    Row = c("A", "A"),
    Column = c(1, 2),
    OD = c(0.1, 0.5)
  )
  dummy_layout <- data.frame(
    Row = "A", Column = 1, Well_ID = "A1", Gene = "WT"
  )
  
  plate <- new_plate(slot_id = "S1L1", growth_data = dummy_growth, layout = dummy_layout)
  expect_s4_class(plate, "Plate")
  expect_equal(plate@slot_id, "S1L1")
})
