#' @include classes.R
#' @import ggplot2
NULL

# Strain names that are always placed (in this order) before mutants.
.cg_control_strains <- c("WT", "PseudoWT", "WT_NoDrug", "PseudoWT_NoDrug",
                         "Blank_NoDrug", "Blank", "CONTROL", "Empty")

#' Build long-format growth curve data for a set of plates
#'
#' @param plates list of Plate objects
#' @param mode "raw" (assays$growth OD) or "corrected" (blank-subtracted OD, unstained reads)
#' @return data.frame with Panel_Title (factor), Clean_Gene, Well_ID, Time_Hours, OD_Plot
#' @noRd
build_growth_curve_data <- function(plates, mode = c("raw", "corrected")) {
  mode <- match.arg(mode)
  plates <- plates[!vapply(plates, is.null, logical(1))]
  if (length(plates) == 0) stop("No plates supplied.")

  trts <- unique(vapply(plates, function(p) as.character(p@treatment), character(1)))
  is_nodrug <- grepl("^no.?drug$", trts, ignore.case = TRUE)
  trt_order <- c(sort(trts[is_nodrug]), sort(trts[!is_nodrug]))

  df_list <- list()
  for (p in plates) {
    if (mode == "raw") {
      df <- p@assays$growth
      if (is.null(df) || nrow(df) == 0) next
      df <- dplyr::left_join(df, p@layout, by = c("Row", "Column"))
      t0 <- if (length(p@t0) == 1 && !is.na(p@t0)) p@t0 else min(df$DateTime, na.rm = TRUE)
      df$Time_Hours <- round(as.numeric(difftime(df$DateTime, t0, units = "hours")), digits = 2)
      df$OD_Plot <- as.numeric(df$OD)
    } else {
      if (is.null(p@metrics$corrected)) p <- compute_metrics(p)
      df <- p@metrics$corrected
      if (is.null(df) || nrow(df) == 0) next
      if ("Is_Stained" %in% names(df)) df <- df[df$Is_Stained == "N", , drop = FALSE]
      t0 <- if (length(p@t0) == 1 && !is.na(p@t0)) p@t0 else min(df$DateTime, na.rm = TRUE)
      # Plain blank subtraction, no floor/clipping: OD_Raw - plate blank median
      df$OD_Plot <- if (all(c("OD_Raw", "Blank_Median_OD") %in% names(df))) {
        as.numeric(df$OD_Raw) - as.numeric(df$Blank_Median_OD)
      } else {
        as.numeric(df$OD_Corrected)
      }
    }

    df <- df[!is.na(df$Gene) & nzchar(stringr::str_trim(df$Gene)), , drop = FALSE]
    if (nrow(df) == 0) next
    if (!"Well_ID" %in% names(df)) df$Well_ID <- sprintf("%s%d", df$Row, as.integer(df$Column))

    date_str <- format(as.Date(t0), "%Y-%m-%d")
    batch_prefix <- if (!is.na(p@label)) stringr::str_extract(p@label, "Batch\\d+") else NA_character_
    line1 <- paste0(if (!is.na(batch_prefix)) paste0(batch_prefix, " ") else "",
                    p@treatment, " R", p@replicate)

    df$Plate_UUID <- p@uuid
    df$Date_Str <- date_str
    df$Treatment <- factor(as.character(p@treatment), levels = trt_order)
    df$Replicate <- as.integer(p@replicate)
    df$Clean_Gene <- stringr::str_trim(stringr::str_remove(df$Gene, "\\s*\\(.*\\)"))
    df$Panel_Title <- paste0(line1, "\n", date_str, " (", p@slot_id, ")")
    df_list[[p@uuid]] <- df[, c("Plate_UUID", "Date_Str", "Treatment", "Replicate", "Panel_Title",
                                "Clean_Gene", "Well_ID", "Time_Hours", "OD_Plot")]
  }
  if (length(df_list) == 0) stop("Selected plates contain no growth data.")
  combined <- dplyr::bind_rows(df_list)

  # One panel per plate: order by treatment, date, replicate; disambiguate duplicate titles
  panels <- combined %>%
    dplyr::distinct(Plate_UUID, Treatment, Date_Str, Replicate, Panel_Title) %>%
    dplyr::arrange(Treatment, Date_Str, Replicate, Panel_Title)
  panels$Panel_Unique <- make.unique(panels$Panel_Title, sep = " #")
  combined$Panel_Title <- factor(panels$Panel_Unique[match(combined$Plate_UUID, panels$Plate_UUID)],
                                 levels = panels$Panel_Unique)
  combined
}

#' Order strains: known controls first, then mutants alphabetically
#' @noRd
order_growth_curve_strains <- function(strains) {
  c(intersect(.cg_control_strains, strains), sort(setdiff(strains, .cg_control_strains)))
}

#' Pick facet columns and page size from panel count (matches reference report layout)
#' @noRd
growth_curve_page_dims <- function(n_panels) {
  if (n_panels <= 2)       list(ncol = 2, w = 9.5,  h = 5.5)
  else if (n_panels <= 8)  list(ncol = 4, w = 11.5, h = 7.5)
  else if (n_panels <= 12) list(ncol = 4, w = 12,   h = 8.5)
  else if (n_panels <= 16) list(ncol = 4, w = 12.5, h = 10)
  else                     list(ncol = 6, w = 15.5, h = 10.5)
}

#' Export per-strain growth curves for a batch to a multi-page PDF
#'
#' Writes one page per strain (WT/control strains first, then mutants
#' alphabetically). Each page has one facet per plate in the batch and one
#' line per well, coloured by well ID.
#'
#' @param registry PlateRegistry
#' @param batch_name character. Name of a batch in \code{registry@batches}.
#' @param file character. Output PDF path.
#' @param mode "raw" (raw OD values from the growth assay, as read) or
#'   "corrected" (OD_Raw minus the plate blank median at each timepoint, with no
#'   floor or clipping). Values are never normalized or rescaled; all wells,
#'   including blanks, are plotted in both modes. Staining reads are excluded.
#' @param strains optional character vector restricting which strains are drawn.
#' @param progress optional function(i, n, strain) called after each page.
#' @return Invisibly, the number of pages written.
#' @export
export_batch_growth_curves_pdf <- function(registry, batch_name, file,
                                           mode = c("raw", "corrected"),
                                           strains = NULL, progress = NULL) {
  mode <- match.arg(mode)
  batch <- registry@batches[[batch_name]]
  if (is.null(batch)) stop(sprintf("Batch '%s' not found in registry.", batch_name))
  plates <- registry@plates[batch@plate_uuids]
  plates <- plates[!vapply(plates, is.null, logical(1))]
  if (length(plates) == 0) stop(sprintf("Batch '%s' references no plates in the registry.", batch_name))

  df <- build_growth_curve_data(plates, mode = mode)

  gene_splits <- split(df, df$Clean_Gene)
  ordered <- order_growth_curve_strains(names(gene_splits))
  if (!is.null(strains)) ordered <- ordered[ordered %in% strains]
  if (length(ordered) == 0) stop("No strains to plot.")

  media <- paste(unique(stats::na.omit(vapply(plates, function(p) as.character(p@media), character(1)))),
                 collapse = "/")
  if (!nzchar(media)) media <- "NA"
  mode_label <- if (mode == "raw") "Raw_OD" else "Corrected_OD"
  mode_title <- if (mode == "raw") "Raw Optical Density (OD 594nm)" else "Corrected Optical Density (OD 594nm - Blank)"

  all_panels <- levels(df$Panel_Title)
  dims <- growth_curve_page_dims(length(all_panels))
  dummy_df <- data.frame(Panel_Title = factor(all_panels, levels = all_panels),
                         Time_Hours = 0, OD_Plot = NA_real_)
  x_max <- max(df$Time_Hours, na.rm = TRUE)
  x_step <- if (x_max > 72) 12 else 6
  x_breaks <- seq(0, ceiling(x_max / x_step) * x_step, by = x_step)

  theme_cg <- theme_bw(base_size = 8.5) +
    theme(
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
      plot.subtitle = element_text(size = 8.5, hjust = 0.5, color = "grey30"),
      strip.background = element_rect(fill = "#f0f4f8", color = "#cbd5e1"),
      strip.text = element_text(face = "bold", size = 7.5),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#f1f5f9", linewidth = 0.35),
      axis.title = element_text(size = 8, face = "bold"),
      axis.text = element_text(size = 7),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 7.5),
      legend.text = element_text(size = 7),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.margin = margin(t = -2, b = 2)
    )

  grDevices::pdf(file, width = dims$w, height = dims$h, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)

  n <- length(ordered)
  for (idx in seq_along(ordered)) {
    strain <- ordered[idx]
    s_df <- gene_splits[[strain]]
    finite_od <- s_df$OD_Plot[is.finite(s_df$OD_Plot)]
    max_od <- if (length(finite_od)) max(finite_od) else 0
    min_od <- if (length(finite_od)) min(finite_od) else 0

    n_wells <- length(unique(s_df$Well_ID))
    well_note <- if (n_wells > 1) sprintf("Replicate Wells on Plate (%d wells)", n_wells) else "Single Well per Plate"

    # Dynamic per-page y-axis range fitted to strain's own data range
    # Add a comfortable 5% margin or small round-up so data doesn't clip
    if (length(finite_od) > 0) {
      span <- max_od - min_od
      padding <- if (span > 0) span * 0.08 else 0.05
      y_min <- if (mode == "raw") max(0, min_od - padding) else (min_od - padding)
      y_max <- max_od + padding
    } else {
      y_min <- 0
      y_max <- 1.0
    }

    g <- ggplot() +
      geom_blank(data = dummy_df, aes(x = Time_Hours, y = OD_Plot)) +
      geom_line(data = s_df, aes(x = Time_Hours, y = OD_Plot, color = Well_ID,
                                 group = interaction(Panel_Title, Well_ID)),
                linewidth = 0.75, alpha = 0.85, na.rm = TRUE) +
      geom_point(data = s_df, aes(x = Time_Hours, y = OD_Plot, color = Well_ID),
                 size = 0.6, alpha = 0.7, na.rm = TRUE) +
      facet_wrap(~ Panel_Title, ncol = dims$ncol, drop = FALSE) +
      scale_y_continuous(limits = c(y_min, y_max), expand = expansion(mult = c(0.02, 0.05))) +
      scale_x_continuous(breaks = x_breaks, expand = expansion(mult = c(0.02, 0.02))) +
      labs(
        title = sprintf("%s | Media: %s | Strain: %s (%s)", batch_name, media, strain, mode_label),
        subtitle = sprintf("Batch: %s | Mode: %s | %s | Max OD: %.3f", batch_name, mode_title, well_note, max_od),
        x = "Time (Hours post-t0)", y = mode_title, color = "Well ID / Replicate"
      ) +
      theme_cg +
      guides(color = guide_legend(nrow = if (n_wells > 10) 2 else 1, byrow = TRUE,
                                  override.aes = list(linewidth = 1.2, size = 1.8)))
    print(g)
    if (is.function(progress)) progress(idx, n, strain)
  }
  invisible(n)
}
