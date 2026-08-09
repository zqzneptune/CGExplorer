#' @import dplyr tidyr
NULL

#' Subtract blanks at the plate level
#'
#' @param df data.frame with Row, Column, Time_Hours, OD_Raw, Type, Is_Stained
#' @param blank_threshold numeric, threshold for bad blanks
#' @return list with corrected `data` and `removed_blanks`
#' @export
subtract_plate_blanks <- function(df, blank_threshold = 0.05) {
  if (!"OD_Raw" %in% names(df) && "OD" %in% names(df)) {
    df$OD_Raw <- df$OD
  }
  
  growth_blanks <- df %>%
    dplyr::filter(Is_Stained == "N" & Type == "Blank") %>%
    dplyr::group_by(Row, Column) %>%
    dplyr::arrange(Time_Hours, .by_group = TRUE) %>%
    dplyr::summarise(
      Initial_OD = dplyr::first(OD_Raw),
      Final_OD = dplyr::last(OD_Raw),
      Growth_Diff = Final_OD - Initial_OD,
      .groups = "drop"
    )
  
  bad_blanks <- growth_blanks %>% dplyr::filter(Growth_Diff > blank_threshold)
  
  clean_blanks <- df %>%
    dplyr::filter(Type == "Blank") %>%
    dplyr::anti_join(bad_blanks, by = c("Row", "Column"))
    
  # If all blanks are contaminated, fallback to all blanks
  if (nrow(clean_blanks) == 0 && nrow(growth_blanks) > 0) {
    message("All blanks are contaminated. Falling back to using all blanks.")
    clean_blanks <- df %>% dplyr::filter(Type == "Blank")
    bad_blanks <- data.frame(Row = character(), Column = character())
  }
  
  plate_blank_medians <- clean_blanks %>%
    dplyr::group_by(Is_Stained, Time_Hours) %>%
    dplyr::summarise(Blank_Median_OD = median(OD_Raw, na.rm = TRUE), .groups = "drop")
  
  df_corr <- df %>%
    dplyr::left_join(plate_blank_medians, by = c("Is_Stained", "Time_Hours")) %>%
    dplyr::mutate(
      Blank_Median_OD = tidyr::replace_na(Blank_Median_OD, 0),
      OD_Corrected = pmax(OD_Raw - Blank_Median_OD, 0.01)
    )
  
  return(list(data = df_corr, removed_blanks = bad_blanks, blank_stats = growth_blanks))
}

#' Extract growth metrics
#'
#' @param df data.frame corrected
#' @param growth_threshold numeric, threshold for valid growth
#' @return data.frame of metrics
#' @export
extract_growth_metrics <- function(df, growth_threshold = 0.02) {
  grp_cols <- intersect(names(df), c("Row", "Column", "Well_ID", "Type", "Gene", "Rep_Suffix", "Is_Edge", "Is_NoDrug_Control"))
  
  df %>%
    dplyr::filter(Is_Stained == "N") %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(grp_cols))) %>%
    dplyr::arrange(Time_Hours, .by_group = TRUE) %>%
    dplyr::summarise(
      Growth_Initial_OD = dplyr::first(OD_Corrected),
      Growth_Endpoint_OD = dplyr::last(OD_Corrected),
      Growth_Endpoint_OD_Raw = dplyr::last(OD_Raw),
      Growth_Max_OD = max(OD_Corrected, na.rm = TRUE),
      Growth_Time_Duration = max(Time_Hours) - min(Time_Hours),
      Growth_AUC = sum(diff(Time_Hours) * (head(OD_Corrected, -1) + tail(OD_Corrected, -1)) / 2),
      Growth_AUC_Rate = ifelse(Growth_Time_Duration > 0, Growth_AUC / Growth_Time_Duration, 0),
      Growth_Max_Rate = if (dplyr::n() > 1) max(diff(OD_Corrected) / diff(Time_Hours), na.rm = TRUE) else 0,
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      Is_Not_Growing = ifelse(
        Type != "Blank" & (Growth_Max_OD - Growth_Initial_OD) < growth_threshold, "Y", "N"
      )
    )
}

#' Extract biofilm metrics
#'
#' @param df data.frame corrected
#' @return data.frame of biofilm metrics
#' @export
extract_biofilm_metrics <- function(df) {
  if (!any(df$Is_Stained == "Y")) return(NULL)
  
  grp_cols <- intersect(names(df), c("Row", "Column", "Well_ID", "Type", "Gene", "Rep_Suffix", "Is_Edge", "Is_NoDrug_Control"))
  
  df %>%
    dplyr::filter(Is_Stained == "Y") %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(grp_cols))) %>%
    dplyr::summarise(
      Biofilm_Mean_OD = mean(OD_Corrected, na.rm = TRUE),
      Biofilm_Stain_OD = dplyr::last(OD_Corrected),
      Biofilm_SD_OD = sd(OD_Corrected, na.rm = TRUE),
      .groups = "drop"
    )
}

#' Compute relative biofilm
#'
#' @param growth_m data.frame growth metrics
#' @param biofilm_m data.frame biofilm metrics
#' @return data.frame joined and scored
#' @export
compute_relative_biofilm <- function(growth_m, biofilm_m) {
  if (is.null(biofilm_m)) return(NULL)
  
  grp_cols <- intersect(names(growth_m), c("Row", "Column", "Well_ID", "Type", "Gene", "Rep_Suffix", "Is_Edge", "Is_NoDrug_Control"))
  
  dplyr::inner_join(growth_m, biofilm_m, by = grp_cols) %>%
    dplyr::mutate(Biofilm_Score_RB = Biofilm_Mean_OD / pmax(Growth_Endpoint_OD, 0.001))
}
