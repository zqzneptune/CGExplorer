#' @include classes.R generics.R
NULL

#' @describeIn add_layout Add layout to PlateRegistry
#' @export
setMethod("add_layout", "PlateRegistry", function(registry, layout_name, layout_data, ...) {
  if (layout_name %in% names(registry@layouts)) {
    stop("Layout name already exists in registry.")
  }
  
  layout_data$uploaded_at <- Sys.time()
  
  registry@layouts[[layout_name]] <- layout_data
  registry@modified_at <- Sys.time()
  registry
})

#' @describeIn get_layout Get layout from PlateRegistry
#' @export
setMethod("get_layout", "PlateRegistry", function(registry, layout_name, ...) {
  if (!layout_name %in% names(registry@layouts)) {
    stop("Layout not found.")
  }
  registry@layouts[[layout_name]]
})

#' @describeIn remove_layout Remove layout from PlateRegistry
#' @export
setMethod("remove_layout", "PlateRegistry", function(registry, layout_name, ...) {
  if (!layout_name %in% names(registry@layouts)) {
    stop("Layout not found.")
  }
  
  registry@layouts[[layout_name]] <- NULL
  registry@modified_at <- Sys.time()
  registry
})

#' @describeIn save_registry Save PlateRegistry to RDS file
#' @export
setMethod("save_registry", "PlateRegistry", function(registry, ...) {
  registry@modified_at <- Sys.time()
  file_path <- file.path(registry@project_dir, "cgexplorer_registry.rds")
  saveRDS(registry, file_path)
  registry
})

#' @describeIn add_plate Add Plate object to PlateRegistry
#' @export
setMethod("add_plate", "PlateRegistry", function(registry, plate, ...) {
  if (plate@uuid %in% names(registry@plates)) stop("Plate already exists. Use update_plate.")
  registry@plates[[plate@uuid]] <- plate
  registry@modified_at <- Sys.time()
  registry
})

#' @describeIn update_plate Update Plate object in PlateRegistry
#' @export
setMethod("update_plate", "PlateRegistry", function(registry, plate, ...) {
  if (!plate@uuid %in% names(registry@plates)) stop("Plate not found in registry.")
  registry@plates[[plate@uuid]] <- plate
  registry@modified_at <- Sys.time()
  registry
})

#' @describeIn remove_plate Remove Plate object from PlateRegistry
#' @export
setMethod("remove_plate", "PlateRegistry", function(registry, uuid, ...) {
  if (!uuid %in% names(registry@plates)) {
    stop("Plate not found in registry.")
  }
  registry@plates[[uuid]] <- NULL
  registry@modified_at <- Sys.time()
  registry
})

#' @describeIn merge_plates Merge multiple plates in PlateRegistry
#' @export
setMethod("merge_plates", "PlateRegistry", function(registry, uuids, ...) {
  if (length(uuids) < 2) stop("Need at least 2 plates to merge.")
  
  plates <- registry@plates[uuids]
  if (any(sapply(plates, is.null))) stop("One or more plates not found in registry.")
  
  # Order by t0
  plates <- plates[order(sapply(plates, function(p) as.numeric(p@t0)))]
  
  p1 <- plates[[1]]
  p1_genes <- sort(p1@layout$Gene)
  
  # Validate
  for (i in 2:length(plates)) {
    p_prev <- plates[[i-1]]
    p_curr <- plates[[i]]
    
    if (!identical(sort(p_curr@layout$Gene), p1_genes)) stop("Plate layouts do not match.")
    if (p_prev@t_end >= p_curr@t0) stop("Time windows overlap - cannot merge.")
    if (p1@treatment != p_curr@treatment) stop("Treatments differ.")
    if (p1@media != p_curr@media) stop("Media differ.")
    if (p1@replicate != p_curr@replicate) stop("Replicate numbers differ.")
  }
  
  # Merge assays
  growth_list <- lapply(plates, function(p) p@assays$growth)
  growth_df <- dplyr::bind_rows(growth_list) %>% dplyr::arrange(DateTime)
  
  staining_assays <- lapply(plates, function(p) p@assays$staining)
  staining_assays <- staining_assays[!sapply(staining_assays, is.null)]
  if (length(staining_assays) > 1) stop("More than one plate has staining data.")
  staining_df <- if (length(staining_assays) == 1) staining_assays[[1]] else NULL
  
  final_assays <- lapply(plates, function(p) p@assays$final)
  final_assays <- final_assays[!sapply(final_assays, is.null)]
  final_df <- if (length(final_assays) > 0) final_assays[[length(final_assays)]] else NULL
  
  # Construct new plate
  new_p <- new("Plate",
               uuid = uuid::UUIDgenerate(),
               slot_id = paste(sapply(plates, function(p) p@slot_id), collapse = "\u2192"),
               label = p1@label,
               treatment = p1@treatment,
               replicate = p1@replicate,
               media = p1@media,
               layout = p1@layout,
               assays = list(growth = growth_df, final = final_df, staining = staining_df),
               t0 = p1@t0,
               t_end = plates[[length(plates)]]@t_end,
               timezone = p1@timezone,
               is_merged = TRUE)
               
  # Re-compute QC metrics on fly for the merged plate
  new_p <- compute_metrics(new_p)

  # Remove old plates
  for (u in uuids) {
    registry@plates[[u]] <- NULL
  }
  
  # Add new plate
  registry@plates[[new_p@uuid]] <- new_p
  registry@modified_at <- Sys.time()
  
  registry
})

#' @describeIn add_batch Add Batch object to PlateRegistry
#' @export
setMethod("add_batch", "PlateRegistry", function(registry, batch, ...) {
  registry@batches[[batch@name]] <- batch
  registry@modified_at <- Sys.time()
  registry
})

#' @describeIn remove_batch Remove Batch object from PlateRegistry
#' @export
setMethod("remove_batch", "PlateRegistry", function(registry, batch_name, ...) {
  if (!batch_name %in% names(registry@batches)) stop("Batch not found.")
  registry@batches[[batch_name]] <- NULL
  registry@modified_at <- Sys.time()
  registry
})
