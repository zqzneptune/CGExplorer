#' Detect staining addition from a kinetic read
#'
#' Computes per-well OD delta between consecutive cycles. If the median plate-wide
#' delta exceeds the threshold, flags that cycle as the start of the staining assay.
#'
#' @param df data.frame with Plate, Row, Column, DateTime, Time_Hours, OD_Raw
#' @param threshold numeric, OD jump threshold (default 0.3)
#' @return A list with `has_staining`, `jump_cycle`, `jump_time`, `staining_hr`
#' @export
detect_staining_jump <- function(df, threshold = 0.3) {
  if (!"OD_Raw" %in% colnames(df) && "OD" %in% colnames(df)) {
    df$OD_Raw <- df$OD
  }
  if (!"Time_Hours" %in% colnames(df) && "DateTime" %in% colnames(df)) {
    t0 <- min(df$DateTime, na.rm = TRUE)
    df$Time_Hours <- round(as.numeric(difftime(df$DateTime, t0, units = "hours")), digits = 2)
  }
  
  # Sort by well and time
  df <- df %>% 
    dplyr::arrange(Plate, Row, Column, DateTime) %>%
    dplyr::group_by(Plate, Row, Column) %>%
    dplyr::mutate(OD_Delta = OD_Raw - dplyr::lag(OD_Raw),
                  Cycle = dplyr::row_number()) %>%
    dplyr::ungroup()
  
  # Calculate plate-wide median delta per cycle
  cycle_stats <- df %>%
    dplyr::group_by(Cycle, DateTime, Time_Hours) %>%
    dplyr::summarise(Median_Delta = median(OD_Delta, na.rm = TRUE), .groups = "drop") %>%
    dplyr::filter(!is.na(Median_Delta))
  
  jump <- cycle_stats %>% dplyr::filter(Median_Delta > threshold)
  
  if (nrow(jump) > 0) {
    first_jump <- jump[1, ]
    return(list(
      has_staining = TRUE,
      jump_cycle = first_jump$Cycle,
      jump_time = first_jump$DateTime,
      staining_hr = first_jump$Time_Hours
    ))
  } else {
    return(list(has_staining = FALSE))
  }
}
