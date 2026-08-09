#' S4 Generics for CGExplorer
#'
#' @name CGExplorer-generics
#' @docType methods
NULL

#' Add layout to registry
#'
#' @param registry PlateRegistry object
#' @param layout_name character name of layout
#' @param layout_data data.frame of layout
#' @param ... additional arguments
#' @export
setGeneric("add_layout", function(registry, layout_name, layout_data, ...) standardGeneric("add_layout"))

#' Get layout from registry
#'
#' @param registry PlateRegistry object
#' @param layout_name character name of layout
#' @param ... additional arguments
#' @export
setGeneric("get_layout", function(registry, layout_name, ...) standardGeneric("get_layout"))

#' Remove layout from registry
#'
#' @param registry PlateRegistry object
#' @param layout_name character name of layout
#' @param ... additional arguments
#' @export
setGeneric("remove_layout", function(registry, layout_name, ...) standardGeneric("remove_layout"))

#' Add assay data to plate
#'
#' @param plate Plate object
#' @param assay_name character name of assay
#' @param slot_data data.frame of assay data
#' @param ... additional arguments
#' @noRd
setGeneric("add_assay", function(plate, assay_name, slot_data, ...) standardGeneric("add_assay"))

#' Compute metrics for plate
#'
#' @param plate Plate object
#' @param ... additional arguments
#' @noRd
setGeneric("compute_metrics", function(plate, ...) standardGeneric("compute_metrics"))

#' Get metrics from plate
#'
#' @param plate Plate object
#' @param metric_name character metric name
#' @param ... additional arguments
#' @noRd
setGeneric("get_metrics", function(plate, metric_name, ...) standardGeneric("get_metrics"))

#' Merge plates in registry
#'
#' @param registry PlateRegistry object
#' @param uuids character vector of plate UUIDs
#' @param ... additional arguments
#' @export
setGeneric("merge_plates", function(registry, uuids, ...) standardGeneric("merge_plates"))

#' Add plate to registry
#'
#' @param registry PlateRegistry object
#' @param plate Plate object
#' @param ... additional arguments
#' @export
setGeneric("add_plate", function(registry, plate, ...) standardGeneric("add_plate"))

#' Remove plate from registry
#'
#' @param registry PlateRegistry object
#' @param uuid character plate UUID
#' @param ... additional arguments
#' @export
setGeneric("remove_plate", function(registry, uuid, ...) standardGeneric("remove_plate"))

#' Update plate in registry
#'
#' @param registry PlateRegistry object
#' @param plate Plate object
#' @param ... additional arguments
#' @export
setGeneric("update_plate", function(registry, plate, ...) standardGeneric("update_plate"))

#' Save registry to disk
#'
#' @param registry PlateRegistry object
#' @param ... additional arguments
#' @export
setGeneric("save_registry", function(registry, ...) standardGeneric("save_registry"))

#' Add batch to registry
#'
#' @param registry PlateRegistry object
#' @param batch Batch object
#' @param ... additional arguments
#' @export
setGeneric("add_batch", function(registry, batch, ...) standardGeneric("add_batch"))

#' Remove batch from registry
#'
#' @param registry PlateRegistry object
#' @param batch_name character batch name
#' @param ... additional arguments
#' @export
setGeneric("remove_batch", function(registry, batch_name, ...) standardGeneric("remove_batch"))

#' Score batch dataset
#'
#' @param batch Batch object
#' @param registry PlateRegistry object
#' @param ... additional arguments
#' @noRd
setGeneric("score_batch", function(batch, registry, ...) standardGeneric("score_batch"))

#' Get scores from batch
#'
#' @param batch Batch object
#' @param strategy character strategy
#' @param metric character metric
#' @param ... additional arguments
#' @noRd
setGeneric("get_scores", function(batch, strategy, metric, ...) standardGeneric("get_scores"))

#' Detect and confirm staining
#'
#' @param plate Plate object
#' @param ... additional arguments
#' @noRd
setGeneric("detect_and_confirm_staining", function(plate, ...) standardGeneric("detect_and_confirm_staining"))
