#' Parse Tecan Plate Reader output (.txt)
#'
#' @param filepath character.
#' @return A data.frame with Plate, DateTime, Row, Column, OD
#' @export
parse_tecan_txt <- function(filepath) {
  if (!file.exists(filepath)) stop(paste("File not found:", filepath))
  text_content <- readr::read_file(filepath)
  lines <- stringr::str_split(text_content, "\r?\n")[[1]]
  data_list <- list()
  i <- 1
  
  while (i <= length(lines)) {
    line <- stringr::str_trim(lines[i])
    if (stringr::str_detect(line, "(?i)^Plate:")) {
      plate_id <- stringr::str_extract(line, "(?<=Plate: )\\w+|(?<=Plate:)\\w+")
      date_str <- stringr::str_extract(line, "(?<=Date: )[^;]+|(?<=Date:)[^;]+")
      time_str <- stringr::str_extract(line, "(?<=Time: )[^ -]+|(?<=Time:)[^ -]+")
      
      dt_str <- paste(stringr::str_trim(date_str), stringr::str_trim(time_str))
      dt_object <- tryCatch(
        {
          lubridate::parse_date_time(
            dt_str,
            orders = c("ymd HMS", "mdy HMS", "dmy HMS", "ymd HM", "mdy HM", "dmy HM", "Ymd HMS", "Ymd HM"),
            quiet = TRUE
          )
        },
        error = function(e) {
          lubridate::ymd_hms(dt_str, quiet = TRUE)
        }
      )
      
      # Determine matrix block size dynamically
      matrix_start <- i + 2
      matrix_end <- matrix_start
      while (matrix_end <= length(lines) && nchar(stringr::str_trim(lines[matrix_end])) > 0 && !stringr::str_detect(lines[matrix_end], "(?i)^Plate:")) {
        matrix_end <- matrix_end + 1
      }
      matrix_end <- matrix_end - 1
      
      if (matrix_end >= matrix_start) {
        matrix_lines_clean <- stringr::str_trim(lines[matrix_start:matrix_end])
        n_rows <- length(matrix_lines_clean)
        
        matrix_df <- tryCatch({
            readr::read_tsv(
              I(matrix_lines_clean),
              col_names = FALSE, show_col_types = FALSE
            )
          }, error = function(e) NULL
        )
        
        if (!is.null(matrix_df) && ncol(matrix_df) > 0) {
          matrix_df$Row <- LETTERS[1:n_rows]
          long_df <- matrix_df %>%
            tidyr::pivot_longer(
              cols = starts_with("X"),
              names_to = "Column_Code", values_to = "OD"
            ) %>%
            dplyr::mutate(
              Column   = sprintf("%02d", as.integer(stringr::str_remove(Column_Code, "X"))),
              Plate    = plate_id,
              DateTime = dt_object
            ) %>%
            dplyr::select(Plate, DateTime, Row, Column, OD)
          data_list[[length(data_list) + 1]] <- long_df
        }
      }
      i <- matrix_end
    }
    i <- i + 1
  }
  if (length(data_list) == 0) return(data.frame())
  dplyr::bind_rows(data_list)
}

#' Parse 96 or 384 well Plate Layout
#'
#' @param filepath character.
#' @return A list containing format info and the parsed data.frame.
#' @export
parse_plate_layout <- function(filepath) {
  if (!file.exists(filepath)) stop(paste("File not found:", filepath))
  ext <- tools::file_ext(filepath)
  if (tolower(ext) %in% c("xlsx", "xls")) {
    raw_map <- readxl::read_excel(filepath)
  } else {
    raw_map <- readr::read_tsv(filepath, col_types = readr::cols(.default = "c"))
  }
  
  colnames(raw_map)[1] <- "Row"
  
  # Clean up row names (in case they have spaces)
  raw_map$Row <- stringr::str_trim(as.character(raw_map$Row))
  
  long_map <- raw_map %>%
    tidyr::pivot_longer(cols = -Row, names_to = "Column", values_to = "Gene") %>%
    dplyr::filter(!is.na(Gene), stringr::str_trim(Gene) != "") %>%
    dplyr::mutate(
      Column = sprintf("%02d", as.integer(Column)),
      Row = toupper(Row)
    )
    
  n_wells <- nrow(long_map)
  cols_int <- as.integer(long_map$Column)
  
  # Format detection
  if (n_wells == 96 && all(long_map$Row %in% LETTERS[1:8]) && all(cols_int %in% 1:12)) {
    format <- "96"
    n_rows <- 8L
    n_cols <- 12L
  } else if (n_wells == 384 && all(long_map$Row %in% LETTERS[1:16]) && all(cols_int %in% 1:24)) {
    format <- "384"
    n_rows <- 16L
    n_cols <- 24L
  } else {
    stop(sprintf("Expected 96 or 384 wells; found %d. Ensure layout dimensions are valid.", n_wells))
  }
  
  processed_map <- long_map %>%
    dplyr::arrange(Row, Column) %>%
    dplyr::mutate(Well_ID = sprintf("%s%s", Row, as.integer(Column))) %>%
    dplyr::mutate(
      Type = dplyr::case_when(
        Gene %in% c("Blank", "neg") ~ "Blank",
        Gene %in% c("WT", "pos") ~ "WT_Control",
        TRUE ~ "Mutant"
      ),
      Is_Edge = (Row == LETTERS[1] | Row == LETTERS[n_rows] | 
                 as.integer(Column) == 1 | as.integer(Column) == n_cols),
      Is_NoDrug_Control = FALSE # Can be set later
    ) %>%
    dplyr::group_by(Gene) %>%
    dplyr::mutate(Rep_Suffix = sprintf("Rep_%02d", dplyr::row_number())) %>%
    dplyr::ungroup() %>%
    dplyr::select(Row, Column, Well_ID, Type, Gene, Rep_Suffix, Is_Edge, Is_NoDrug_Control)
    
  list(
    data = processed_map,
    format = format,
    n_wells = as.integer(n_wells),
    n_rows = n_rows,
    n_cols = n_cols
  )
}
