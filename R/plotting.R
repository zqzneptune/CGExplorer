#' @import ggplot2
NULL

#' @noRd
plot_od_distribution <- function(df) {
  if (is.null(df) || nrow(df) == 0) {
    return(plotly::plot_ly() %>%
      plotly::add_annotations(text = "No data", showarrow = FALSE) %>%
      plotly::config(displayModeBar = FALSE))
  }
  
  df_f <- df %>%
    dplyr::filter(!is.na(Type)) %>%
    dplyr::mutate(
      Type_Display = dplyr::case_when(
        Type == "Blank" ~ "Blank",
        stringr::str_detect(Type, "^WT") ~ "WT",
        TRUE ~ "Mutant"
      ),
      Type_Display = factor(Type_Display, levels = c("Blank", "WT", "Mutant"))
    )
    
  pal_solid <- c(Blank = "#757575", WT = "#1E88E5", Mutant = "#E53935")
  pal_fill  <- c(Blank = "rgba(117,117,117,0.45)", WT = "rgba(30,136,229,0.45)", Mutant = "rgba(229,57,53,0.45)")
  
  p <- plotly::plot_ly()
  for (t in levels(df_f$Type_Display)) {
    sub <- df_f %>% dplyr::filter(Type_Display == t)
    if (nrow(sub) == 0) next
    
    p <- p %>% plotly::add_trace(
      data = sub, x = ~OD_Raw, type = "violin", name = t,
      fillcolor = pal_fill[[t]],
      line = list(color = pal_solid[[t]], width = 1.5),
      marker = list(color = pal_solid[[t]]),
      spanmode = "soft", side = "positive",
      box = list(visible = TRUE, fillcolor = pal_fill[[t]], line = list(color = pal_solid[[t]])), 
      meanline = list(visible = TRUE, color = pal_solid[[t]]),
      hovertemplate = paste0("<b>", t, "</b><br>OD: %{x:.3f}<extra></extra>")
    )
  }
  p %>%
    plotly::layout(
      xaxis = list(title = "Raw OD Value"),
      yaxis = list(title = "", showticklabels = FALSE),
      legend = list(title = list(text = "Well Type"), orientation = "h", x = 0, y = 1.15),
      margin = list(l = 30, r = 10, t = 30, b = 40)
    ) %>%
    plotly::config(displayModeBar = FALSE)
}

#' @noRd
plot_blank_stats <- function(blank_stats, removed_blanks, threshold = 0.05) {
  if (is.null(blank_stats) || nrow(blank_stats) == 0) {
    return(plotly::plot_ly() %>%
      plotly::add_annotations(text = "No blank wells", showarrow = FALSE) %>%
      plotly::config(displayModeBar = FALSE))
  }
  
  # Protect against empty removed_blanks
  if (is.null(removed_blanks) || nrow(removed_blanks) == 0) {
    flagged_wells <- character(0)
  } else {
    flagged_wells <- paste0(removed_blanks$Row, removed_blanks$Column)
  }
  
  df <- blank_stats %>%
    dplyr::mutate(
      Well    = paste0(Row, Column),
      Flagged = Well %in% flagged_wells,
      Color   = ifelse(Flagged, "red", "steelblue")
    )
  plotly::plot_ly(df,
    x = ~Well, y = ~Growth_Diff, type = "scatter", mode = "markers",
    marker = list(color = ~Color, size = 7),
    text = ~ paste0(
      "Well: ", Well, "<br>\u0394OD: ", round(Growth_Diff, 4),
      "<br>", ifelse(Flagged, "FLAGGED", "OK")
    ),
    hoverinfo = "text"
  ) %>%
    plotly::add_segments(
      x = 0, xend = nrow(df) + 1,
      y = threshold, yend = threshold,
      line = list(color = "red", dash = "dot", width = 1),
      showlegend = FALSE, inherit = FALSE
    ) %>%
    plotly::layout(
      xaxis = list(title = "Well", tickangle = -45),
      yaxis = list(title = "\u0394OD"),
      margin = list(l = 40, r = 10, t = 10, b = 60)
    ) %>%
    plotly::config(displayModeBar = FALSE)
}

#' @noRd
plot_wt_stats <- function(wt_stats, threshold = 0.02) {
  if (is.null(wt_stats) || nrow(wt_stats) == 0) {
    return(plotly::plot_ly() %>%
      plotly::add_annotations(text = "No WT wells", showarrow = FALSE) %>%
      plotly::config(displayModeBar = FALSE))
  }
  df <- wt_stats %>%
    dplyr::mutate(
      Well  = paste0(Row, Column),
      Delta = Growth_Max_OD - Growth_Initial_OD,
      Color = ifelse(Is_Not_Growing == "Y", "red", "steelblue")
    )
  plotly::plot_ly(df,
    x = ~Well, y = ~Delta, type = "bar",
    marker = list(color = ~Color),
    text = ~ paste0(
      "Well: ", Well, "<br>\u0394OD: ", round(Delta, 4),
      "<br>", ifelse(Is_Not_Growing == "Y", "NO-GROW", "Growing")
    ),
    hoverinfo = "text"
  ) %>%
    plotly::add_segments(
      x = -0.5, xend = nrow(df) - 0.5,
      y = threshold, yend = threshold,
      line = list(color = "red", dash = "dot", width = 1),
      showlegend = FALSE, inherit = FALSE
    ) %>%
    plotly::layout(
      xaxis = list(title = "WT Well", tickangle = -45),
      yaxis = list(title = "\u0394OD"),
      margin = list(l = 40, r = 10, t = 10, b = 60)
    ) %>%
    plotly::config(displayModeBar = FALSE)
}

#' @noRd
plot_growth_curves <- function(plate_df, stain_hr = 26) {
  if (is.null(plate_df) || nrow(plate_df) == 0) {
    return(plotly::plot_ly() %>%
      plotly::add_annotations(text = "No data", showarrow = FALSE) %>%
      plotly::config(displayModeBar = FALSE))
  }

  pal <- c(Blank = "#888888", WT_Control = "#2196F3", Mutant = "#E53935")

  df <- plate_df %>%
    dplyr::filter(!is.na(Type), !is.na(OD_Raw)) %>%
    dplyr::mutate(
      Type = factor(Type, levels = c("Blank", "WT_Control", "Mutant")),
      Well = paste0(Row, Column)
    )

  med_df <- df %>%
    dplyr::group_by(Type, Time_Hours) %>%
    dplyr::summarise(Med_OD = median(OD_Raw, na.rm = TRUE), .groups = "drop")

  max_od <- max(df$OD_Raw, na.rm = TRUE)
  max_time <- max(df$Time_Hours, na.rm = TRUE)

  p <- plotly::plot_ly()

  for (typ in c("Blank", "WT_Control", "Mutant")) {
    col <- pal[[typ]]
    sub <- df %>% dplyr::filter(Type == typ)
    med_sub <- med_df %>% dplyr::filter(Type == typ)
    if (nrow(sub) == 0) next

    p <- p %>% plotly::add_trace(
      data = sub,
      x = ~Time_Hours, y = ~OD_Raw,
      type = "violin",
      name = typ,
      legendgroup = typ,
      fillcolor = paste0(col, "28"),
      line = list(color = col, width = 0.8),
      meanline = list(visible = FALSE),
      points = FALSE,
      spanmode = "soft",
      scalemode = "width",
      width = 0.9,
      side = "positive",
      showlegend = FALSE,
      hoverinfo = "skip"
    )

    p <- p %>% plotly::add_trace(
      data = sub,
      x = ~Time_Hours, y = ~OD_Raw,
      type = "scatter", mode = "markers",
      name = typ,
      legendgroup = typ,
      marker = list(color = col, size = 3, opacity = 0.2),
      text = ~ paste0(Well, "<br>OD: ", round(OD_Raw, 3), "<br>t=", Time_Hours, "h"),
      hoverinfo = "text",
      showlegend = TRUE
    )

    p <- p %>% plotly::add_trace(
      data = med_sub,
      x = ~Time_Hours, y = ~Med_OD,
      type = "scatter", mode = "lines",
      name = paste0(typ, " median"),
      legendgroup = typ,
      line = list(color = col, width = 2.5),
      hoverinfo = "skip",
      showlegend = FALSE
    )
  }

  if (!is.na(stain_hr)) {
    p <- p %>%
      plotly::add_segments(
        x = stain_hr, xend = stain_hr,
        y = 0, yend = max_od * 1.15,
        line = list(color = "#7B1FA2", dash = "dash", width = 1.5),
        name = paste0("Staining (", stain_hr, "h)"),
        showlegend = TRUE,
        inherit = FALSE
      )
  }

  p %>%
    plotly::layout(
      xaxis = list(
        title = "Time (hours)", zeroline = FALSE,
        range = c(-0.5, max_time + 0.5)
      ),
      yaxis = list(title = "OD (raw)", zeroline = FALSE),
      legend = list(
        orientation = "h", x = 0, y = 1.08,
        font = list(size = 12), traceorder = "normal"
      ),
      margin = list(l = 60, r = 20, t = 45, b = 55)
    ) %>%
    plotly::config(displayModeBar = FALSE)
}

#' Plot Plate QC Layout
#'
#' @param plate Plate object
#' @param metric character. Column name in metrics to plot.
#' @param colorscale character. Color scheme key.
#' @param cmin numeric. Minimum value for color scale.
#' @param cmax numeric. Maximum value for color scale.
#' @return plotly object
#' @export
plot_plate_qc_layout <- function(plate, metric = "Growth_Endpoint_OD", colorscale = "plasma", cmin = NULL, cmax = NULL) {
  # Get data
  if (is.null(plate@metrics$growth)) return(plotly::plot_ly() %>% plotly::add_annotations(text="No metrics computed", showarrow=FALSE))
  
  df <- plate@metrics$growth
  if (metric == "Biofilm_Stain_OD") {
    if (!is.null(plate@metrics$biofilm)) {
      df <- dplyr::left_join(df, plate@metrics$biofilm, by = intersect(names(df), names(plate@metrics$biofilm)))
    } else {
       return(plotly::plot_ly() %>% plotly::add_annotations(text="No Biofilm metrics computed", showarrow=FALSE))
    }
  }
  
  if (!metric %in% names(df)) {
    return(plotly::plot_ly() %>% plotly::add_annotations(text=paste("Metric", metric, "not found"), showarrow=FALSE))
  }
  
  # For metrics other than Raw OD, set Blank values to 0 to avoid misleading metrics
  if (metric != "Growth_Endpoint_OD_Raw") {
    df[[metric]][df$Type == "Blank"] <- 0
  }
  
  # Parse Row and Column to numeric for plotting
  df$Row_Num <- match(df$Row, LETTERS)
  df$Col_Num <- as.numeric(df$Column)
  
  # Dynamic grid bounds for 96-well vs 384-well
  max_r <- max(df$Row_Num, na.rm = TRUE)
  max_c <- max(df$Col_Num, na.rm = TRUE)
  r_ticks <- if (max_r <= 8) 1:8 else 1:16
  c_ticks <- if (max_c <= 12) 1:12 else 1:24
  
  # Legend names & Shapes: square for Blank, triangle for WT, round circle for Mutant
  df$Well_Type <- "Mutant"
  df$Well_Type[df$Type == "Blank"] <- "Blank"
  df$Well_Type[grep("^WT", df$Type)] <- "WT"
  df$Well_Type <- factor(df$Well_Type, levels = c("Blank", "WT", "Mutant"))
  
  # Palette selection
  pal <- switch(colorscale,
    "staining_purple" = c("#E0E0E0", "#9C27B0", "#4A148C"),
    "growth_heatmap"  = c("#313695", "#E0F3F8", "#FEE090", "#A50026"),
    "greens"          = c("#F7FCF5", "#74C476", "#00441B"),
    "blues"           = c("#F7FBFF", "#6BAED6", "#08306B"),
    "plasma"          = "Plasma",
    "magma"           = "Magma",
    "inferno"         = "Inferno",
    "cividis"         = "Cividis",
    "Plasma"
  )
  
  # Tooltip
  df$HoverText <- sprintf("<b>%s</b><br>Well: %s%s<br>Type: %s<br>Value: %.3f", 
                          df$Gene, df$Row, df$Column, df$Well_Type, df[[metric]])
  
  marker_opts <- list(size = 14, line = list(color = 'rgba(0,0,0,0.5)', width = 1))
  symbols_map <- c("Blank" = "square", "WT" = "triangle-up", "Mutant" = "circle")
  
  vals <- df[[metric]]
  vals <- vals[!is.na(vals)]
  val_min <- if (!is.null(cmin) && !is.na(cmin)) cmin else if (length(vals) > 0) min(vals) else 0
  val_max <- if (!is.null(cmax) && !is.na(cmax)) cmax else if (length(vals) > 0) max(vals) else 1
  
  p <- plotly::plot_ly()
  
  for (wt in c("Blank", "WT", "Mutant")) {
    sub_df <- df %>% dplyr::filter(Well_Type == wt)
    if (nrow(sub_df) == 0) next
    
    m_opts <- marker_opts
    m_opts$symbol <- symbols_map[[wt]]
    m_opts$cauto <- FALSE
    m_opts$cmin <- val_min
    m_opts$cmax <- val_max
    
    p <- p %>% plotly::add_trace(
      data = sub_df,
      x = ~Col_Num, y = ~Row_Num,
      color = stats::as.formula(paste0("~", metric)),
      colors = pal,
      marker = m_opts,
      name = wt,
      text = ~HoverText, hoverinfo = "text",
      type = 'scatter', mode = 'markers',
      showlegend = TRUE
    )
  }
  
  p %>% plotly::layout(
    yaxis = list(title = "", tickvals = r_ticks, ticktext = LETTERS[r_ticks], autorange = "reversed", 
                 zeroline = FALSE, showgrid = FALSE),
    xaxis = list(title = "", tickvals = c_ticks, side = "top", 
                 zeroline = FALSE, showgrid = FALSE),
    plot_bgcolor = "#f9f9f9",
    margin = list(t = 40, b = 20, l = 40, r = 20)
  ) %>%
  plotly::colorbar(title = metric)
}

#' Plot Multi-Strain Growth Curves
#'
#' @param df data.frame containing Time_Hours, OD_Corrected, OD_Raw, Gene, Type, Rep_Suffix, Row, Column
#' @param selected_genes character vector of selected genes/strains
#' @param use_raw_od logical. TRUE for OD_Raw, FALSE for OD_Corrected
#' @param color_by character. "strain", "replicate", or "type"
#' @return plotly object
#' @export
plot_multi_strain_growth_curves <- function(df, selected_genes = NULL, use_raw_od = FALSE, color_by = "strain") {
  if (is.null(df) || nrow(df) == 0) {
    return(plotly::plot_ly() %>% plotly::add_annotations(text = "No growth data available", showarrow = FALSE))
  }
  
  y_col <- if (use_raw_od && "OD_Raw" %in% names(df)) "OD_Raw" else if ("OD_Corrected" %in% names(df)) "OD_Corrected" else "OD"
  
  if (is.null(selected_genes) || length(selected_genes) == 0) {
    avail_genes <- sort(unique(df$Gene[!is.na(df$Gene) & nzchar(df$Gene)]))
    if (length(avail_genes) > 0) {
      selected_genes <- avail_genes[1]
    }
  }
  
  if (!is.null(selected_genes) && length(selected_genes) > 0) {
    df <- df %>% dplyr::filter(Gene %in% selected_genes | Type %in% selected_genes)
  }
  
  if (nrow(df) == 0) {
    return(plotly::plot_ly() %>% plotly::add_annotations(text = "No data for selected strains", showarrow = FALSE))
  }
  
  if (!"Well" %in% names(df)) df$Well <- paste0(df$Row, df$Column)
  if (!"Plate_Label" %in% names(df) && "label" %in% names(df)) df$Plate_Label <- df$label
  if (!"Plate_Label" %in% names(df)) df$Plate_Label <- "Plate"
  
  if (color_by == "replicate") {
    df$Group_Var <- paste(df$Plate_Label, ifelse("Rep_Suffix" %in% names(df), df$Rep_Suffix, ""), sep = " ")
  } else if (color_by == "type") {
    df$Group_Var <- df$Type
  } else {
    df$Group_Var <- df$Gene
  }
  
  df$Trace_ID <- paste(df$Plate_Label, df$Gene, df$Well, sep = " - ")
  
  p <- plotly::plot_ly()
  
  groups <- sort(unique(df$Group_Var))
  n_g <- length(groups)
  pal_colors <- if (n_g <= 8) {
    RColorBrewer::brewer.pal(max(3, n_g), "Set1")[1:n_g]
  } else {
    grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(n_g)
  }
  colors_map <- setNames(pal_colors, groups)
  
  for (g in groups) {
    sub_df <- df %>% dplyr::filter(Group_Var == g)
    traces <- unique(sub_df$Trace_ID)
    
    first_trace <- TRUE
    for (tr in traces) {
      tr_df <- sub_df %>% dplyr::filter(Trace_ID == tr) %>% dplyr::arrange(Time_Hours)
      
      p <- p %>% plotly::add_trace(
        data = tr_df,
        x = ~Time_Hours,
        y = stats::as.formula(paste0("~", y_col)),
        type = "scatter",
        mode = "lines+markers",
        name = g,
        legendgroup = g,
        showlegend = first_trace,
        line = list(color = colors_map[[g]], width = 2),
        marker = list(color = colors_map[[g]], size = 4),
        text = ~paste0("<b>Strain: </b>", Gene, 
                      "<br><b>Well: </b>", Well, 
                      "<br><b>Plate: </b>", Plate_Label,
                      "<br><b>Time: </b>", round(Time_Hours, 1), " h",
                      "<br><b>OD: </b>", round(get(y_col), 3)),
        hoverinfo = "text"
      )
      first_trace <- FALSE
    }
  }
  
  p %>% plotly::layout(
    xaxis = list(title = "Time (hours)", zeroline = FALSE),
    yaxis = list(title = if (y_col == "OD_Raw") "Raw OD" else "Corrected OD", zeroline = FALSE),
    legend = list(title = list(text = paste("Grouped by", color_by)), orientation = "h", x = 0, y = 1.15),
    margin = list(l = 50, r = 20, t = 40, b = 50)
  )
}

#' Plot Replicate Growth Curves
#'
#' @param replicate_data data.frame containing Time_Hours, OD_Raw or OD_Corrected, Row, Column, Gene, Rep_Suffix
#' @param stain_hr numeric staining hour threshold
#' @return plotly object
#' @export
plot_replicate_growth_curves <- function(replicate_data, stain_hr = 26) {
  if (is.null(replicate_data) || nrow(replicate_data) == 0) {
    return(plotly::plot_ly() %>%
      plotly::add_annotations(text = "No replicate data available", showarrow = FALSE) %>%
      plotly::config(displayModeBar = FALSE))
  }

  if (!"Well" %in% names(replicate_data)) {
    replicate_data <- replicate_data %>%
      dplyr::mutate(Well = paste0(Row, Column))
  }
  replicate_data <- replicate_data %>% dplyr::arrange(Well, Time_Hours)

  y_col <- if ("OD_Raw" %in% names(replicate_data)) "OD_Raw" else if ("OD_Corrected" %in% names(replicate_data)) "OD_Corrected" else if ("OD" %in% names(replicate_data)) "OD" else NULL
  if (is.null(y_col) || !y_col %in% names(replicate_data)) {
    return(plotly::plot_ly() %>%
      plotly::add_annotations(text = "No valid OD column for growth curves", showarrow = FALSE) %>%
      plotly::config(displayModeBar = FALSE))
  }

  wells <- unique(replicate_data$Well)
  colors_list <- c("#2EC4B6", "#FF9F1C", "#E71D36", "#9C27B0", "#4CAF50", "#00BCD4", "#3F51B5", "#FFEB3B")

  p <- plotly::plot_ly()

  for (idx in seq_along(wells)) {
    w <- wells[idx]
    w_data <- replicate_data %>% dplyr::filter(Well == w)
    color_val <- colors_list[(idx - 1) %% length(colors_list) + 1]
    rep_lbl <- if ("Rep_Suffix" %in% names(w_data)) w_data$Rep_Suffix[1] else if ("Plate_UUID" %in% names(w_data)) w_data$Plate_UUID[1] else ""

    p <- p %>% plotly::add_trace(
      data = w_data,
      x = ~Time_Hours, y = stats::as.formula(paste0("~", y_col)),
      type = "scatter", mode = "lines+markers",
      name = paste0(w, ifelse(nzchar(rep_lbl), paste0(" (", rep_lbl, ")"), "")),
      line = list(color = color_val, width = 2),
      marker = list(color = color_val, size = 6),
      text = ~ paste0("Well: ", Well, "<br>Time: ", Time_Hours, "h<br>OD: ", round(get(y_col), 3)),
      hoverinfo = "text"
    )
  }

  if (!is.null(stain_hr) && !is.na(stain_hr)) {
    max_od <- max(replicate_data[[y_col]], na.rm = TRUE)

    p <- p %>%
      plotly::add_segments(
        x = stain_hr, xend = stain_hr,
        y = 0, yend = max_od * 1.1,
        line = list(color = "#7B1FA2", dash = "dash", width = 1.5),
        name = paste0("Staining (", stain_hr, "h)"),
        showlegend = TRUE,
        inherit = FALSE
      )
  }

  p %>%
    plotly::layout(
      xaxis = list(title = "Time (hours)", zeroline = FALSE),
      yaxis = list(title = if (y_col == "OD_Raw") "OD (raw)" else "OD", zeroline = FALSE),
      legend = list(orientation = "h", x = 0, y = -0.15),
      margin = list(l = 60, r = 20, t = 40, b = 60)
    ) %>%
    plotly::config(displayModeBar = TRUE)
}

