if (!suppressWarnings(require(CGExplorer, quietly = TRUE))) {
  r_files <- list.files("../../R", full.names = TRUE)
  if (length(r_files) > 0) lapply(r_files, source)
}

make_pdf_test_plate <- function(slot_id, treatment, rep) {
  layout <- data.frame(
    Row = c("A", "A", "A", "B"), Column = c("01", "02", "03", "01"),
    Well_ID = c("A1", "A2", "A3", "B1"),
    Type = c("Blank", "WT_Control", "Mutant", "Mutant"),
    Gene = c("Blank", "WT", "geneB", "geneA (dup)"),
    stringsAsFactors = FALSE
  )
  times <- as.POSIXct("2026-08-11 10:00:00", tz = "America/Halifax") + 3600 * 0:3
  growth <- expand.grid(Row_Col = seq_len(nrow(layout)), t = seq_along(times))
  growth <- data.frame(
    DateTime = times[growth$t],
    Row = layout$Row[growth$Row_Col],
    Column = layout$Column[growth$Row_Col],
    OD = dplyr::case_when(
      layout$Type[growth$Row_Col] == "Blank" ~ 0.1,
      layout$Gene[growth$Row_Col] == "geneB" ~ 0.05,          # below blank -> negative corrected OD
      TRUE ~ 0.1 + 0.1 * (growth$t - 1)
    ),
    stringsAsFactors = FALSE
  )
  p <- new_plate(slot_id, growth, layout)
  p@treatment <- treatment
  p@replicate <- as.integer(rep)
  p@media <- "SCFM"
  p
}

test_that("growth curve data puts one panel per plate and orders strains", {
  plates <- list(make_pdf_test_plate("S1L1", "Drug1", 1), make_pdf_test_plate("S2L1", "NoDrug", 1))
  df <- build_growth_curve_data(plates, mode = "raw")
  expect_equal(nlevels(df$Panel_Title), 2)
  # NoDrug panels sort first
  expect_match(levels(df$Panel_Title)[1], "^NoDrug R1")
  # Parenthetical annotations are stripped from strain names
  expect_true("geneA" %in% df$Clean_Gene)
  expect_equal(order_growth_curve_strains(unique(df$Clean_Gene)),
               c("WT", "Blank", "geneA", "geneB"))
})

test_that("PDF data uses unmodified raw OD and unclipped blank-subtracted OD", {
  p <- make_pdf_test_plate("S1L1", "NoDrug", 1)
  raw <- build_growth_curve_data(list(p), mode = "raw")
  expect_equal(sort(raw$OD_Plot), sort(p@assays$growth$OD))

  corr <- build_growth_curve_data(list(p), mode = "corrected")
  # Blank median is 0.1 at every timepoint
  expect_equal(unique(round(corr$OD_Plot[corr$Clean_Gene == "geneB"], 6)), -0.05)
  expect_equal(unique(round(corr$OD_Plot[corr$Clean_Gene == "Blank"], 6)), 0)
  wt <- corr[corr$Clean_Gene == "WT", ]
  expect_equal(round(wt$OD_Plot[order(wt$Time_Hours)], 6), c(0, 0.1, 0.2, 0.3))
})

test_that("export_batch_growth_curves_pdf writes one page per strain", {
  reg <- new_registry(project_name = "PDF Test")
  p1 <- make_pdf_test_plate("S1L1", "NoDrug", 1)
  p2 <- make_pdf_test_plate("S1L2", "NoDrug", 2)
  reg <- add_plate(reg, p1)
  reg <- add_plate(reg, p2)
  reg <- add_batch(reg, new_batch("B1", c(p1@uuid, p2@uuid)))

  f <- tempfile(fileext = ".pdf")
  calls <- 0
  n <- export_batch_growth_curves_pdf(reg, "B1", f, mode = "raw",
                                      progress = function(i, n, s) calls <<- calls + 1)
  expect_equal(n, 4)
  expect_equal(calls, 4)
  expect_true(file.exists(f) && file.size(f) > 0)

  # Corrected mode keeps every strain, including blanks
  f2 <- tempfile(fileext = ".pdf")
  expect_equal(export_batch_growth_curves_pdf(reg, "B1", f2, mode = "corrected"), 4)

  expect_error(export_batch_growth_curves_pdf(reg, "nope", tempfile()), "not found")
})
