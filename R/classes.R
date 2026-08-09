#' @include generics.R
NULL

# -- 1. Plate Class ----------------------------------------------------------

#' S4 Class representing a Chemo-Genomic Plate
#'
#' @slot uuid character. Internal UUID.
#' @slot slot_id character. E.g. "S1L1".
#' @slot label character. User-facing name.
#' @slot treatment character. E.g. "NoDrug", "Drug1".
#' @slot replicate integer. Biological replicate number.
#' @slot media character. Growth media label.
#' @slot layout data.frame. Well map.
#' @slot assays list. Named list of assay data.frames (growth, final, staining).
#' @slot staining_hr numeric. Hours post t0 when staining was performed (NA if unknown).
#' @slot staining_confirmed logical. TRUE after user confirms staining detection.
#' @slot t0 POSIXct. Earliest DateTime in assays$growth.
#' @slot t_end POSIXct. Latest DateTime across all non-NULL assays.
#' @slot timezone character. Timezone, fixed to "America/Halifax".
#' @slot metrics list. Populated by compute_metrics().
#' @slot qc_flags list. QC results.
#' @slot provenance list. Audit trail records.
#' @slot is_merged logical. TRUE if created from merge_plates().
#' @export
setClass(
  "Plate",
  slots = c(
    uuid = "character",
    slot_id = "character",
    label = "character",
    treatment = "character",
    replicate = "integer",
    media = "character",
    layout = "data.frame",
    assays = "list",
    staining_hr = "numeric",
    staining_confirmed = "logical",
    t0 = "POSIXct",
    t_end = "POSIXct",
    timezone = "character",
    metrics = "list",
    qc_flags = "list",
    provenance = "list",
    is_merged = "logical"
  ),
  prototype = list(
    uuid = NA_character_,
    slot_id = NA_character_,
    label = NA_character_,
    treatment = NA_character_,
    replicate = 1L,
    media = NA_character_,
    layout = data.frame(),
    assays = list(growth = data.frame(), final = NULL, staining = NULL),
    staining_hr = NA_real_,
    staining_confirmed = FALSE,
    timezone = "America/Halifax",
    metrics = list(),
    qc_flags = list(),
    provenance = list(),
    is_merged = FALSE
  )
)

setValidity("Plate", function(object) {
  errors <- character()
  if (is.null(object@assays$growth) || nrow(object@assays$growth) == 0) {
    errors <- c(errors, "@assays$growth must be a non-empty data.frame")
  }
  if (length(object@uuid) != 1 || is.na(object@uuid)) {
    errors <- c(errors, "uuid must be a valid length-1 character string")
  }
  if (length(errors) == 0) TRUE else errors
})

# -- 2. Batch Class ----------------------------------------------------------

#' S4 Class representing a Scoring Batch
#'
#' @slot uuid character.
#' @slot name character. User-assigned name.
#' @slot plate_uuids character. Ordered Plate UUIDs included in this batch.
#' @slot scores list. Scoring results per strategy/metric.
#' @slot scoring_params list. Params used in last score_batch().
#' @slot created_at POSIXct.
#' @export
setClass(
  "Batch",
  slots = c(
    uuid = "character",
    name = "character",
    plate_uuids = "character",
    scores = "list",
    scoring_params = "list",
    created_at = "POSIXct"
  ),
  prototype = list(
    uuid = NA_character_,
    name = NA_character_,
    plate_uuids = character(),
    scores = list(),
    scoring_params = list(),
    created_at = Sys.time()
  )
)

setValidity("Batch", function(object) {
  errors <- character()
  if (length(object@plate_uuids) < 2) {
    errors <- c(errors, "A Batch must reference at least 2 plate UUIDs")
  }
  if (length(errors) == 0) TRUE else errors
})

# -- 3. PlateRegistry Class --------------------------------------------------

#' S4 Class representing the Project Plate Registry
#'
#' @slot plates list. Named list: UUID -> Plate S4 object.
#' @slot batches list. Named list: batch_name -> Batch S4 object.
#' @slot layouts list. Named layout registry.
#' @slot project_name character. User-assigned project name.
#' @slot project_dir character. Absolute path to directory holding registry.rds.
#' @slot schema_version character. Package version at save time.
#' @slot created_at POSIXct.
#' @slot modified_at POSIXct.
#' @export
setClass(
  "PlateRegistry",
  slots = c(
    plates = "list",
    batches = "list",
    layouts = "list",
    project_name = "character",
    project_dir = "character",
    schema_version = "character",
    created_at = "POSIXct",
    modified_at = "POSIXct"
  ),
  prototype = list(
    plates = list(),
    batches = list(),
    layouts = list(),
    project_name = "Untitled Project",
    project_dir = ".",
    schema_version = "0.3.0",
    created_at = Sys.time(),
    modified_at = Sys.time()
  )
)
