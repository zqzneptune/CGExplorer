#' @include classes.R generics.R
NULL

#' Compute metrics for a Plate
#'
#' @param plate A Plate object
#' @param blank_thr Numeric threshold for flagging blanks
#' @param growth_thr Numeric threshold for flagging non-growing wells
#' @return A Plate object with @metrics and @qc_flags populated
#' @noRd
setMethod("compute_metrics", "Plate", function(plate, blank_thr = 0.05, growth_thr = 0.02, ...) {
  # 1. Combine Growth, Staining, and Final assays with Layout
  growth_df <- plate@assays$growth
  if (is.null(growth_df) || nrow(growth_df) == 0) stop("No growth data available on this plate.")
  
  growth_df$Is_Stained <- "N"
  
  stain_df <- plate@assays$staining
  if (!is.null(stain_df) && nrow(stain_df) > 0) {
    stain_df$Is_Stained <- "Y"
  } else {
    stain_df <- NULL
  }
  
  combined <- dplyr::bind_rows(growth_df, stain_df)
  if (!"OD_Raw" %in% names(combined) && "OD" %in% names(combined)) {
    combined$OD_Raw <- combined$OD
  }
  
  # Join layout
  combined <- combined %>% dplyr::left_join(plate@layout, by = c("Row", "Column"))
  
  # Ensure Time_Hours is calculated relative to t0
  t0 <- plate@t0
  combined$Time_Hours <- round(as.numeric(difftime(combined$DateTime, t0, units = "hours")), digits = 2)
  
  # 2. Subtract Blanks
  res <- subtract_plate_blanks(combined, blank_threshold = blank_thr)
  corr_df <- res$data
  
  # 3. Extract Metrics
  growth_m <- extract_growth_metrics(corr_df, growth_threshold = growth_thr)
  biofilm_m <- extract_biofilm_metrics(corr_df)
  
  if (!is.null(biofilm_m)) {
    rel_biofilm <- compute_relative_biofilm(growth_m, biofilm_m)
    # Merge into a single biofilm metrics dataframe if desired, but PRD asks for them separate or combined?
    # PRD says: "corrected", "growth", "biofilm"
  } else {
    rel_biofilm <- NULL
  }
  
  # WT Stats
  wt_stats <- growth_m %>%
    dplyr::filter(Type == "WT" | Type == "WT_Control")
  
  # Automated Failure Checks
  all_blanks_bad <- isTRUE(res$all_blanks_contaminated)
  all_wt_bad <- (nrow(wt_stats) > 0 && all(wt_stats$Is_Not_Growing == "Y"))
  auto_qc_passed <- !(all_blanks_bad || all_wt_bad)
  
  # Check for existing user override flag
  existing_user_status <- if (!is.null(plate@qc_flags$user_override)) plate@qc_flags$user_override else NULL
  final_qc_passed <- if (!is.null(existing_user_status)) {
    identical(existing_user_status, "PASSED")
  } else {
    auto_qc_passed
  }
  
  failure_reasons <- character(0)
  if (all_blanks_bad) failure_reasons <- c(failure_reasons, "All blank wells contaminated")
  if (all_wt_bad) failure_reasons <- c(failure_reasons, "All WT control wells failed to grow")
  if (!is.null(existing_user_status) && existing_user_status == "FAILED") {
    failure_reasons <- c(failure_reasons, "Manually set to FAILED by user")
  }

  # 4. Populate Plate
  plate@metrics <- list(
    corrected = corr_df,
    growth = growth_m,
    biofilm = rel_biofilm
  )
  
  plate@qc_flags <- list(
    qc_passed = final_qc_passed,
    auto_qc_passed = auto_qc_passed,
    user_override = existing_user_status,
    failure_reasons = failure_reasons,
    blank_stats = res$blank_stats,
    removed_blanks = res$removed_blanks,
    wt_stats = wt_stats
  )
  
  return(plate)
})
