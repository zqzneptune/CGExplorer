#' @include classes.R
NULL

#' Constructor for Plate
#'
#' @param slot_id character. Unique identifier for the plate slot (e.g. "S1L1").
#' @param growth_data data.frame. Raw growth assay readings.
#' @param layout data.frame. Well layout definition.
#' @param timezone character. Timezone string, defaults to "America/Halifax".
#' @return A \code{\link{Plate-class}} object.
#' @export
new_plate <- function(slot_id, growth_data, layout, timezone = "America/Halifax") {
  if (nrow(growth_data) == 0) stop("growth_data cannot be empty")
  
  t0 <- min(growth_data$DateTime, na.rm = TRUE)
  t_end <- max(growth_data$DateTime, na.rm = TRUE)
  
  p <- new("Plate",
           uuid = uuid::UUIDgenerate(),
           slot_id = slot_id,
           layout = layout,
           assays = list(growth = growth_data, final = NULL, staining = NULL),
           t0 = t0,
           t_end = t_end,
           timezone = timezone)
  
  # Basic validity check
  validObject(p)
  p
}

#' Constructor for Batch
#'
#' @param name character. Name of the batch.
#' @param plate_uuids character vector. UUIDs of plates included in this batch.
#' @return A \code{\link{Batch-class}} object.
#' @export
new_batch <- function(name, plate_uuids) {
  b <- new("Batch",
           uuid = uuid::UUIDgenerate(),
           name = name,
           plate_uuids = plate_uuids)
  validObject(b)
  b
}

#' Constructor for PlateRegistry
#'
#' @param project_name character. Project display name.
#' @param project_dir character. Path to project directory storing registry.
#' @return A \code{\link{PlateRegistry-class}} object.
#' @export
new_registry <- function(project_name = "Untitled Project", project_dir = getwd()) {
  r <- new("PlateRegistry",
           project_name = project_name,
           project_dir = project_dir)
  r
}

#' Load a PlateRegistry from RDS
#'
#' @param path character. Path to the RDS file or directory containing cgexplorer_registry.rds.
#' @return A \code{\link{PlateRegistry-class}} object loaded from disk.
#' @export
load_registry <- function(path) {
  file_path <- if (dir.exists(path)) file.path(path, "cgexplorer_registry.rds") else path
  if (!file.exists(file_path)) stop("Registry file not found at ", file_path)
  
  reg <- readRDS(file_path)
  
  # Verify object
  if (!inherits(reg, "PlateRegistry")) {
    stop("File is not a PlateRegistry object.")
  }
  
  # Check schema migration here later if needed
  current_version <- utils::packageVersion("CGExplorer")
  if (reg@schema_version != as.character(current_version)) {
    warning(sprintf("Registry schema version (%s) differs from current package version (%s).", 
                    reg@schema_version, current_version))
  }
  
  reg
}
