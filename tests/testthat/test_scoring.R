if (!suppressWarnings(require(CGExplorer, quietly = TRUE))) {
  r_files <- list.files("../../R", full.names = TRUE)
  if (length(r_files) > 0) lapply(r_files, source)
}

make_synthetic_pheno_df <- function() {
  set.seed(42)
  genes <- c("WT", paste0("gene", sprintf("%02d", 1:20)))
  rows <- list()
  for (trt in c("NoDrug", "DrugA", "DrugB")) {
    for (rep in 1:3) {
      for (g in genes) {
        is_wt <- g == "WT"
        base_auc <- if (is_wt) 12.0 else stats::rnorm(1, 11.5, 0.8)
        if (g == "gene01" && trt == "DrugA") base_auc <- 3.0   # strong sensitive hit
        if (g == "gene02" && trt == "DrugA") base_auc <- 20.0  # strong resistant hit
        if (g == "gene20") base_auc <- 0.2                     # low-growth mutant
        drug_eff <- if (trt == "NoDrug") 1.0 else 0.85
        val_auc <- max(0.05, base_auc * drug_eff + stats::rnorm(1, 0, 0.25))
        val_end <- val_auc / 20
        rows[[length(rows) + 1L]] <- data.frame(
          Batch = "TestBatch",
          Plate_UUID = paste("P", trt, rep, sep = "_"),
          Media = "SCFM",
          Treatment = trt,
          Replicate = rep,
          Gene = g,
          Type = if (is_wt) "WT_Control" else "Mutant",
          Exclude_From_Scoring = "N",
          Is_NoDrug_Control = if (trt == "NoDrug") "Y" else "N",
          Growth_Initial_OD = 0.05,
          Growth_Max_OD = 0.05 + val_end,
          Growth_Endpoint_OD = val_end,
          Growth_AUC = val_auc,
          Biofilm_Stain_OD = val_auc * 0.1,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  dplyr::bind_rows(rows)
}

test_that("huber_pheno_stats handles normal, outlier, single, and empty vectors", {
  empty_res <- huber_pheno_stats(c(NA_real_, NA_real_))
  expect_true(is.na(empty_res$mu_p))
  expect_equal(empty_res$k_eff, 0)

  single_res <- huber_pheno_stats(5.2)
  expect_equal(single_res$mu_p, 5.2)
  expect_equal(single_res$var_p, 0)
  expect_equal(single_res$k_eff, 1)

  # Outlier should be down-weighted relative to arithmetic mean
  vals <- c(10.0, 10.1, 9.9, 10.2, 50.0)
  rob <- huber_pheno_stats(vals)
  expect_lt(rob$mu_p, mean(vals))
  expect_gt(rob$mu_p, 9.5)
  expect_lt(rob$k_eff, length(vals))
})

test_that("score_batch_metric_combined computes all 4 strategies cleanly", {
  pheno <- make_synthetic_pheno_df()

  res_unfiltered <- score_batch_dataset(
    pheno,
    strategy_type = "Metric_Combined",
    metrics = c("Growth_AUC", "Biofilm_Stain_OD"),
    control_drug = "NoDrug",
    exclude_low_growth = FALSE
  )

  expect_s3_class(res_unfiltered, "data.frame")
  expect_setequal(unique(res_unfiltered$Combine_Strategy), c("Min", "Max", "Average", "Robust_Variance_Aware"))
  # 20 mutants x 2 drugs x 2 metrics x 4 strategies = 320
  expect_equal(nrow(res_unfiltered), 320L)
  expect_true(all(c("cg_score", "log_cg_score", "z_score", "p_value", "adj_p", "se_log2_cg", "shrinkage_factor") %in% names(res_unfiltered)))

  # Strong sensitive (gene01) and resistant (gene02) signals in DrugA Growth_AUC
  avg_a <- res_unfiltered %>%
    dplyr::filter(Combine_Strategy == "Average", Treatment == "DrugA", Metric == "Growth_AUC")
  expect_lt(avg_a$log_cg_score[avg_a$Gene == "gene01"], -1.0)
  expect_gt(avg_a$log_cg_score[avg_a$Gene == "gene02"], 0.5)

  # Robust Variance-Aware populates SE and shrinkage_factor in (0, 1]
  rob_a <- res_unfiltered %>%
    dplyr::filter(Combine_Strategy == "Robust_Variance_Aware", Treatment == "DrugA", Metric == "Growth_AUC")
  expect_true(all(!is.na(rob_a$se_log2_cg) & rob_a$se_log2_cg >= 0))
  expect_true(all(rob_a$shrinkage_factor > 0 & rob_a$shrinkage_factor <= 1))

  # Optional low-growth QC excludes gene20
  res_qc <- score_batch_metric_combined(
    pheno,
    metrics = "Growth_AUC",
    control_drug = "NoDrug",
    exclude_low_growth = TRUE
  )
  expect_true("gene20" %in% attr(res_qc, "excluded_low_growth_genes"))
  expect_false("gene20" %in% res_qc$Gene)
})

test_that("S4 score_batch and mod_scoring_server work end-to-end for One-vs-One and Merge", {
  # Build a minimal valid PlateRegistry with 4 plates (2 NoDrug, 2 DrugA)
  layout_df <- data.frame(
    Row = rep("A", 12),
    Column = sprintf("%02d", 1:12),
    Type = c("Blank", "Blank", "WT_Control", "WT_Control", rep("Mutant", 8)),
    Gene = c("Blank", "Blank", "WT", "WT", paste0("mut", 1:8)),
    Rep_Suffix = rep("R1", 12),
    stringsAsFactors = FALSE
  )

  make_plate_obj <- function(slot, trt, rep_num, mult = 1.0) {
    t_vals <- c(0, 5, 10)
    recs <- lapply(seq_len(nrow(layout_df)), function(i) {
      g <- layout_df$Gene[i]
      tp <- layout_df$Type[i]
      inc <- if (tp == "Blank") 0.0 else if (g == "mut1" && trt == "DrugA") 0.05 else 0.4 * mult
      data.frame(
        DateTime = as.POSIXct("2026-10-01 00:00:00") + t_vals * 3600,
        Row = layout_df$Row[i],
        Column = layout_df$Column[i],
        OD = c(0.08, 0.08 + inc * 0.5, 0.08 + inc),
        stringsAsFactors = FALSE
      )
    })
    growth_df <- dplyr::bind_rows(recs)
    p <- new_plate(slot, growth_df, layout_df)
    p@media <- "SCFM"
    p@treatment <- trt
    p@replicate <- as.integer(rep_num)
    compute_metrics(p)
  }

  p1 <- make_plate_obj("S1L1", "NoDrug", 1, 1.0)
  p2 <- make_plate_obj("S1L2", "NoDrug", 2, 1.02)
  p3 <- make_plate_obj("S1L3", "DrugA",  1, 0.85)
  p4 <- make_plate_obj("S1L4", "DrugA",  2, 0.82)

  reg <- new_registry(project_name = "Scoring Test", project_dir = tempdir())
  for (p in list(p1, p2, p3, p4)) reg <- add_plate(reg, p)
  b <- new_batch("B1", c(p1@uuid, p2@uuid, p3@uuid, p4@uuid))
  reg <- add_batch(reg, b)

  # S4 score_batch with Metric_Combined
  b_mc <- score_batch(b, reg, strategies = "Metric_Combined", metrics = "Growth_AUC", control_drug = "NoDrug")
  sc_mc <- get_scores(b_mc, "Metric_Combined", "Growth_AUC")
  expect_true(nrow(sc_mc) == 32L) # 8 mutants x 4 strategies
  expect_setequal(unique(sc_mc$Combine_Strategy), c("Min", "Max", "Average", "Robust_Variance_Aware"))

  # Shiny mod_scoring_server end-to-end
  rv_mock <- shiny::reactiveValues(registry = reg)
  shiny::testServer(mod_scoring_server, args = list(rv = rv_mock), {
    # 1. Select batch & One-vs-One mode, then click Run Scoring
    session$setInputs(
      selected_batch = "B1",
      metric = "Growth_AUC",
      control_drug = "NoDrug",
      score_mode = "one_vs_one",
      strategy = "NoDrug_Plates",
      replicate_pair = "Control_Rep_1_vs_Drug_Rep_1",
      score_alpha = 0.05
    )
    session$setInputs(btn_score = 1)
    expect_equal(nrow(local_rv$scoring_run), 8L)
    invisible(output$scoring_plot)
    invisible(output$dt_scores)

    # 2. Merge mode (Variance-Aware)
    session$setInputs(
      score_mode = "merge",
      combine_strategy = "Robust_Variance_Aware",
      exclude_low_growth = FALSE,
      btn_score = 2
    )
    expect_equal(nrow(local_rv$scoring_run), 8L)
    expect_equal(nrow(local_rv$scoring_run_all4), 32L)
    invisible(output$scoring_plot)
    invisible(output$dt_scores)
    invisible(output$strategy_comparison_plot)
    invisible(output$strategy_comparison_table)
  })
})
