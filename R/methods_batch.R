#' @include classes.R generics.R
NULL

#' Validate plate eligibility and layout uniformity for a batch
#'
#' @param registry PlateRegistry
#' @param uuids character vector of plate UUIDs
#' @return TRUE if valid, otherwise throws an error
#' @export
validate_batch_plates <- function(registry, uuids) {
  if (length(uuids) < 2) stop("A Batch must reference at least 2 plates.")
  
  plates <- registry@plates[uuids]
  for (u in uuids) {
    p <- plates[[u]]
    if (is.null(p)) stop(sprintf("Plate %s not found in registry.", u))
    
    # 1. Usability / QC Passed Check
    is_passed <- if (!is.null(p@qc_flags$qc_passed)) isTRUE(p@qc_flags$qc_passed) else TRUE
    if (!is_passed) {
      plate_name <- ifelse(is.na(p@label) || p@label == "", p@slot_id, p@label)
      stop(sprintf("Plate '%s' is marked FAILED and cannot be included in a batch.", plate_name))
    }
  }
  
  # 2. Strict Layout Uniformity Check
  ref_layout <- plates[[1]]@layout
  if (nrow(ref_layout) > 0) {
    for (i in 2:length(plates)) {
      curr_layout <- plates[[i]]@layout
      if (!isTRUE(all.equal(ref_layout, curr_layout, check.attributes = FALSE))) {
        stop("Layout mismatch detected: All plates in a batch must share the exact same 384-well layout sheet.")
      }
    }
  }
  
  return(TRUE)
}

#' Assemble pheno_df for scoring
#'
#' @param registry PlateRegistry
#' @param uuids character vector of plate UUIDs
#' @return data.frame ready for scoring
#' @noRd
assemble_pheno_df <- function(registry, uuids, batch_name) {
  validate_batch_plates(registry, uuids)
  
  list_dfs <- list()
  for (u in uuids) {
    p <- registry@plates[[u]]
    if (is.null(p)) next
    if (is.null(p@metrics) || length(p@metrics) == 0) {
      tryCatch({
        p <- compute_metrics(p)
        registry@plates[[u]] <- p
      }, error = function(e) NULL)
    }
    
    growth_m <- p@metrics$growth
    if (is.null(growth_m)) {
      warning(sprintf("Plate %s lacks computed metrics, skipping.", u))
      next
    }
    
    # Extract phenotype columns
    df <- growth_m %>% 
      dplyr::mutate(
        Plate_UUID = u,
        Batch = batch_name,
        Treatment = p@treatment,
        Media = p@media,
        Replicate = p@replicate,
        # Default flags
        Exclude_From_Scoring = "N"
      )
    
    # Check removed blanks flag
    if (!is.null(p@qc_flags$removed_blanks) && nrow(p@qc_flags$removed_blanks) > 0) {
      bad_wells <- paste0(p@qc_flags$removed_blanks$Row, p@qc_flags$removed_blanks$Column)
      df <- df %>% dplyr::mutate(
        Is_Contaminated_Blank = ifelse(paste0(Row, Column) %in% bad_wells, "Y", "N")
      )
    } else {
      df$Is_Contaminated_Blank <- "N"
    }
    
    has_nodrug_col <- "Is_NoDrug_Control" %in% names(df)
    is_nodrug_vec <- if (has_nodrug_col) (df$Is_NoDrug_Control == TRUE | df$Is_NoDrug_Control == "Y") else FALSE
    
    df <- df %>% dplyr::mutate(
      Is_NonGrowing_Control = ifelse(is_nodrug_vec & Is_Not_Growing == "Y", "Y", "N"),
      Exclude_From_Scoring = ifelse(Is_Contaminated_Blank == "Y" | Is_NonGrowing_Control == "Y", "Y", "N")
    )
    
    # Join biofilm metrics if they exist
    biofilm_m <- p@metrics$biofilm
    if (!is.null(biofilm_m)) {
      grp_cols <- intersect(names(df), names(biofilm_m))
      df <- df %>% dplyr::left_join(biofilm_m, by = grp_cols)
    }
    
    list_dfs[[u]] <- df
  }
  dplyr::bind_rows(list_dfs)
}

#' @noRd
setMethod("score_batch", "Batch", function(batch, registry, strategies, metrics, control_drug = "NoDrug", ...) {
  if (length(batch@plate_uuids) == 0) stop("Batch has no plates.")
  
  pheno_df <- assemble_pheno_df(registry, batch@plate_uuids, batch@name)
  if (nrow(pheno_df) == 0) stop("No valid metrics found across batch plates.")
  
  if (is.null(batch@scores)) batch@scores <- list()
  
  for (strat in strategies) {
    if (is.null(batch@scores[[strat]])) batch@scores[[strat]] <- list()
    
    res <- score_batch_dataset(pheno_df, strategy_type = strat, metrics = metrics, control_drug = control_drug)
    if (nrow(res) > 0) {
      for (met in unique(res$Metric)) {
        batch@scores[[strat]][[met]] <- res %>% dplyr::filter(Metric == met)
      }
    }
  }
  
  batch@scoring_params <- list(
    strategies = strategies,
    metrics = metrics,
    control_drug = control_drug
  )
  
  batch
})

#' @noRd
setMethod("get_scores", "Batch", function(batch, strategy, metric, ...) {
  if (is.null(batch@scores)) return(NULL)
  if (is.null(batch@scores[[strategy]])) return(NULL)
  return(batch@scores[[strategy]][[metric]])
})
