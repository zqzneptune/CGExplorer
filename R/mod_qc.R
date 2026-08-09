mod_plate_explorer_ui <- function(id) {
  ns <- NS(id)
  box(
    title = tagList(
      textOutput(ns("box_title"), inline = TRUE),
      div(
        style = "float: right;",
        actionButton(ns("btn_close_explorer"), "Close", class = "btn btn-default btn-xs", style = "color: #333;", icon = icon("times"))
      )
    ),
    width = 12, status = "info", solidHeader = TRUE, collapsible = TRUE,
    tabsetPanel(
      tabPanel("Virtual Plate Layout",
               br(),
               fluidRow(
                 column(3, selectInput(ns("metric"), "Metric", 
                                       choices = c("Corrected OD (Endpoint)" = "Growth_Endpoint_OD", 
                                                   "Raw Final Growth OD" = "Growth_Endpoint_OD_Raw", 
                                                   "Growth AUC" = "Growth_AUC", 
                                                   "Final Staining OD" = "Biofilm_Stain_OD"))),
                 column(3, selectInput(ns("palette"), "Color Scheme",
                                       choices = c("Staining (Purple)" = "staining_purple",
                                                   "Growth Heatmap (Blue-Red)" = "growth_heatmap",
                                                   "Greens" = "greens",
                                                   "Blues" = "blues"))),
                 column(3, numericInput(ns("cmin"), "Min Value Color", value = NA, step = 0.05)),
                 column(3, numericInput(ns("cmax"), "Max Value Color", value = NA, step = 0.05))
               ),
               plotly::plotlyOutput(ns("plot_layout"), height = "500px")
      ),
      tabPanel("Stats",
               br(),
               tabsetPanel(
                 tabPanel("Raw Data Stats", plotly::plotlyOutput(ns("plot_raw_dist"), height = "450px")),
                 tabPanel("Blanks Flagging", plotly::plotlyOutput(ns("plot_blanks"), height = "450px")),
                 tabPanel("WT Controls Flagging", plotly::plotlyOutput(ns("plot_wt"), height = "450px")),
                 tabPanel("Flagged Summary", DT::dataTableOutput(ns("dt_flagged")))
               )
      ),
      tabPanel("Growth Curve",
               br(),
               uiOutput(ns("growth_curve_controls_ui")),
               plotly::plotlyOutput(ns("plot_growth_curves_multi"), height = "500px")
      )
    )
  )
}

#' @import shiny
#' @noRd
mod_plate_explorer_server <- function(id, plate, rv = NULL, on_close = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    observeEvent(input$btn_close_explorer, {
      if (!is.null(on_close)) {
        on_close()
      }
    })
    
    output$box_title <- renderText({
      if (is.null(plate())) return("No Plate")
      p <- plate()
      m_str <- if (is.null(p@media) || is.na(p@media) || !nzchar(p@media)) "Media" else p@media
      t_str <- if (is.null(p@treatment) || is.na(p@treatment) || !nzchar(p@treatment)) "NoDrug" else p@treatment
      r_str <- if (is.null(p@replicate) || is.na(p@replicate)) "1" else as.character(p@replicate)
      s_str <- if (is.null(p@slot_id) || is.na(p@slot_id) || !nzchar(p@slot_id)) "" else sprintf("(%s)", p@slot_id)
      
      sprintf("Plate QC: %s + %s + Rep%s %s", m_str, t_str, r_str, s_str)
    })
    
    # Auto-populate default min/max when metric changes
    observeEvent(input$metric, {
      p <- plate()
      req(p, p@metrics$growth)
      m <- input$metric
      df <- p@metrics$growth
      if (m == "Biofilm_Stain_OD" && !is.null(p@metrics$biofilm)) {
        df <- dplyr::left_join(df, p@metrics$biofilm, by = intersect(names(df), names(p@metrics$biofilm)))
      }
      if (m %in% names(df)) {
        vals <- df[[m]]
        if (m != "Growth_Endpoint_OD_Raw") vals[df$Type == "Blank"] <- 0
        vals <- vals[!is.na(vals)]
        if (length(vals) > 0) {
          updateNumericInput(session, "cmin", value = round(min(vals, na.rm = TRUE), 3))
          updateNumericInput(session, "cmax", value = round(max(vals, na.rm = TRUE), 3))
        }
      }
    })
    
    output$plot_layout <- plotly::renderPlotly({
      p <- plate()
      req(p)
      plot_plate_qc_layout(p, metric = input$metric, colorscale = input$palette, cmin = input$cmin, cmax = input$cmax)
    })
    
    output$plot_raw_dist <- plotly::renderPlotly({
      p <- plate()
      req(p, p@metrics$corrected)
      plot_od_distribution(p@metrics$corrected)
    })
    
    output$plot_blanks <- plotly::renderPlotly({
      p <- plate()
      req(p, p@qc_flags$blank_stats)
      plot_blank_stats(p@qc_flags$blank_stats, p@qc_flags$removed_blanks)
    })
    
    output$plot_wt <- plotly::renderPlotly({
      p <- plate()
      req(p, p@qc_flags$wt_stats)
      plot_wt_stats(p@qc_flags$wt_stats)
    })
    
    output$dt_flagged <- DT::renderDataTable({
      p <- plate()
      req(p, p@qc_flags$removed_blanks, p@qc_flags$wt_stats)
      
      bad_blanks <- p@qc_flags$removed_blanks
      if (nrow(bad_blanks) > 0) {
        bad_blanks$Reason <- "Contaminated Blank"
        bad_blanks$Well <- paste0(bad_blanks$Row, bad_blanks$Column)
      }
      
      bad_wt <- p@qc_flags$wt_stats %>% dplyr::filter(Is_Not_Growing == "Y")
      if (nrow(bad_wt) > 0) {
        bad_wt$Reason <- "Non-Growing WT"
        bad_wt$Well <- paste0(bad_wt$Row, bad_wt$Column)
      }
      
      cols_to_keep <- c("Well", "Reason")
      
      df_blanks <- if (nrow(bad_blanks) > 0) bad_blanks[, cols_to_keep] else data.frame(Well=character(), Reason=character())
      df_wt <- if (nrow(bad_wt) > 0) bad_wt[, cols_to_keep] else data.frame(Well=character(), Reason=character())
      
      df <- rbind(df_blanks, df_wt)
      if (nrow(df) == 0) return(data.frame(Message = "No flagged wells"))
      
      DT::datatable(df, rownames = FALSE, options = list(pageLength = 10, dom = 'ftip'))
    })
    
    output$growth_curve_controls_ui <- renderUI({
      p <- plate()
      req(p)
      
      df <- p@metrics$corrected
      if (is.null(df)) df <- p@assays$growth
      if (is.null(df)) return(p("No growth data available"))
      
      if (!"Gene" %in% names(df) && nrow(p@layout) > 0) {
        df <- df %>% dplyr::left_join(p@layout, by = c("Row", "Column"))
      }
      
      genes <- sort(unique(df$Gene[!is.na(df$Gene) & nzchar(df$Gene)]))
      first_gene <- if (length(genes) > 0) genes[1] else character(0)
      
      fluidRow(
        column(5, 
               selectInput(ns("growth_strains"), "Select Strains / Mutants to Display",
                           choices = genes, selected = first_gene, multiple = TRUE)
        ),
        column(3,
               selectInput(ns("growth_y_metric"), "OD Metric",
                           choices = c("Corrected OD" = "corrected", "Raw OD" = "raw"))
        ),
        column(4,
               selectInput(ns("growth_color_by"), "Color Lines By",
                           choices = c("Strain / Mutant" = "strain", "Replicate / Plate" = "replicate", "Well Type" = "type"))
        )
      )
    })
    
    output$plot_growth_curves_multi <- plotly::renderPlotly({
      p <- plate()
      req(p)
      
      df <- p@metrics$corrected
      if (is.null(df)) df <- p@assays$growth
      if (is.null(df)) return(plotly::plot_ly() %>% plotly::add_annotations(text = "No growth data", showarrow = FALSE))
      
      if (!"Gene" %in% names(df) && nrow(p@layout) > 0) {
        df <- df %>% dplyr::left_join(p@layout, by = c("Row", "Column"))
      }
      
      df$Plate_Label <- sprintf("%s + %s + Rep%s", p@media, p@treatment, p@replicate)
      
      strains <- input$growth_strains
      if (is.null(strains)) {
        genes <- sort(unique(df$Gene[!is.na(df$Gene) & nzchar(df$Gene)]))
        strains <- if (length(genes) > 0) genes[1] else character(0)
      }
      
      use_raw <- isTRUE(input$growth_y_metric == "raw")
      color_by <- if (!is.null(input$growth_color_by)) input$growth_color_by else "strain"
      
      plot_multi_strain_growth_curves(df, selected_genes = strains, use_raw_od = use_raw, color_by = color_by)
    })
  })
}

#' @import shiny
#' @noRd
mod_qc_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    tags$head(
      tags$style(HTML("
        .qc-top-panel {
          background: #ffffff;
          border-radius: 8px;
          padding: 16px;
          margin-bottom: 20px;
          box-shadow: 0 1px 3px rgba(0,0,0,0.05);
          border: 1px solid #e2e8f0;
        }
        .qc-scroll-container {
          display: flex;
          gap: 16px;
          overflow-x: auto;
          overflow-y: hidden;
          padding: 6px 2px 12px 2px;
          scrollbar-width: thin;
          scrollbar-color: #cbd5e1 #f1f5f9;
        }
        .qc-scroll-container::-webkit-scrollbar {
          height: 6px;
        }
        .qc-scroll-container::-webkit-scrollbar-track {
          background: #f1f5f9;
          border-radius: 3px;
        }
        .qc-scroll-container::-webkit-scrollbar-thumb {
          background: #cbd5e1;
          border-radius: 3px;
        }
        .qc-date-group {
          flex: 0 0 auto;
          background: #f8fafc;
          border-radius: 8px;
          padding: 10px 12px;
          border: 1px solid #e2e8f0;
          display: flex;
          flex-direction: column;
          gap: 8px;
        }
        .qc-date-header {
          font-weight: 600;
          font-size: 13px;
          color: #3b82f6;
          padding-bottom: 4px;
          border-bottom: 1px solid #e2e8f0;
          display: flex;
          align-items: center;
          gap: 6px;
          white-space: nowrap;
        }
        .qc-cards-row {
          display: flex;
          gap: 8px;
          flex-wrap: nowrap;
        }
        .qc-plate-card-btn {
          border-radius: 6px !important;
          padding: 6px 12px !important;
          font-size: 12px !important;
          font-weight: 500 !important;
          cursor: pointer !important;
          display: inline-flex !important;
          align-items: center !important;
          gap: 6px !important;
          border: 1px solid #cbd5e1 !important;
          background-color: #ffffff !important;
          color: #334155 !important;
          white-space: nowrap !important;
          transition: all 0.15s ease-in-out !important;
          box-shadow: 0 1px 2px rgba(0,0,0,0.04) !important;
        }
        .qc-plate-card-btn:hover {
          border-color: #3c8dbc !important;
          color: #3c8dbc !important;
          background-color: #f0f9ff !important;
        }
        .qc-plate-card-btn.active {
          background-color: #3c8dbc !important;
          color: #ffffff !important;
          border-color: #3c8dbc !important;
          box-shadow: 0 2px 4px rgba(60,141,188,0.3) !important;
        }
        .qc-plate-card-btn.active:hover {
          background-color: #357ca5 !important;
          border-color: #357ca5 !important;
          color: #ffffff !important;
        }
        .qc-explorers-container {
          display: flex;
          flex-wrap: wrap;
          gap: 16px;
          align-items: flex-start;
          width: 100%;
          margin-top: 10px;
        }
        .qc-explorer-card {
          box-sizing: border-box;
          background: #ffffff;
          border-radius: 8px;
          box-shadow: 0 1px 4px rgba(0,0,0,0.08);
        }
      "))
    ),
    fluidRow(
      column(12,
        div(class = "qc-top-panel",
          div(style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 8px;",
            h4(strong("Plate Selection"), style = "margin: 0; color: #1e293b; font-size: 16px;"),
            uiOutput(ns("active_count_status"), inline = TRUE)
          ),
          uiOutput(ns("plate_cards_scroll_ui"))
        )
      )
    ),
    fluidRow(
      column(12,
        uiOutput(ns("qc_explorers_ui"))
      )
    )
  )
}

safe_extract_plate_date <- function(t0_val) {
  if (is.null(t0_val) || length(t0_val) == 0 || all(is.na(t0_val))) {
    return(Sys.Date())
  }
  if (inherits(t0_val, "Date")) {
    return(t0_val[1])
  }
  if (inherits(t0_val, "POSIXt")) {
    return(as.Date(t0_val[1]))
  }
  str_val <- as.character(t0_val[1])
  if (!nzchar(str_val) || is.na(str_val)) {
    return(Sys.Date())
  }
  d <- tryCatch(as.Date(str_val, format = "%Y-%m-%d"), error = function(e) NA)
  if (!is.na(d)) return(d)
  d <- tryCatch(as.Date(as.POSIXct(str_val)), error = function(e) NA)
  if (!is.na(d)) return(d)
  return(Sys.Date())
}

#' @import shiny
#' @noRd
mod_qc_server <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Store currently toggled plate UUIDs for QC
    qc_plates_selected <- reactiveVal(character(0))
    
    # Track registered observer UUIDs to avoid duplicate observers
    registered_uuids <- reactiveVal(character(0))
    
    # Track initialized explorer module servers
    initialized_explorers <- reactiveVal(character(0))

    # Keep selected UUIDs valid if plates are deleted/merged
    observe({
      plates <- rv$registry@plates
      sel <- qc_plates_selected()
      valid_sel <- intersect(sel, names(plates))
      if (length(valid_sel) != length(sel)) {
        qc_plates_selected(valid_sel)
      }
    })
    
    # Register dynamic click observers for all plates in registry
    observe({
      plates <- rv$registry@plates
      uuids <- names(plates)
      reg <- registered_uuids()
      new_uuids <- setdiff(uuids, reg)
      
      for (u in new_uuids) {
        local({
          my_uuid <- u
          btn_id <- paste0("btn_toggle_qc_", my_uuid)
          
          observeEvent(input[[btn_id]], {
            curr_selected <- qc_plates_selected()
            if (my_uuid %in% curr_selected) {
              # Toggle OFF
              qc_plates_selected(setdiff(curr_selected, my_uuid))
            } else {
              # Toggle ON: check and compute QC metrics if needed
              p <- rv$registry@plates[[my_uuid]]
              if (!is.null(p)) {
                if (is.null(p@metrics) || length(p@metrics) == 0) {
                  tryCatch({
                    p_qc <- compute_metrics(p)
                    rv$registry@plates[[my_uuid]] <- p_qc
                    save_registry(rv$registry)
                  }, error = function(e) {
                    shinyjs::alert(paste("QC computation failed for plate:", e$message))
                  })
                }
              }
              qc_plates_selected(c(curr_selected, my_uuid))
            }
          }, ignoreInit = TRUE)
        })
      }
      
      if (length(new_uuids) > 0) {
        registered_uuids(c(reg, new_uuids))
      }
    })
    
    # Status text and manual format refresh for active selection count
    output$active_count_status <- renderUI({
      sel <- qc_plates_selected()
      n_total <- length(rv$registry@plates)
      n_sel <- length(sel)
      
      div(style = "display: flex; align-items: center; gap: 12px;",
        span(
          style = "color: #64748b; font-size: 13px;",
          sprintf("%d of %d plates selected", n_sel, n_total)
        ),
        actionButton(
          ns("btn_force_refresh_format"),
          "Refresh Layout",
          icon = icon("sync"),
          class = "btn btn-default btn-xs",
          style = "color: #3b82f6; border-color: #cbd5e1;"
        )
      )
    })
    
    # Observe plate addition/removal to force refresh panel format and plot sizing
    observeEvent(qc_plates_selected(), {
      shinyjs::runjs("
        setTimeout(function() {
          window.dispatchEvent(new Event('resize'));
          $('.js-plotly-plot').each(function() {
            if (window.Plotly && this.layout) {
              Plotly.Plots.resize(this);
            }
          });
          if (typeof makeExplorersInteractive === 'function') makeExplorersInteractive();
        }, 200);
      ")
    }, ignoreInit = TRUE)
    
    observeEvent(input$btn_force_refresh_format, {
      shinyjs::runjs("
        window.dispatchEvent(new Event('resize'));
        $('.js-plotly-plot').each(function() {
          if (window.Plotly && this.layout) {
            Plotly.Plots.resize(this);
          }
        });
        if (typeof makeExplorersInteractive === 'function') makeExplorersInteractive();
      ")
    })
    
    # Top Panel: Horizontally scrollable compact plate cards, grouped & ordered by date (newest to oldest)
    output$plate_cards_scroll_ui <- renderUI({
      plates <- rv$registry@plates
      if (length(plates) == 0) {
        return(div(style = "text-align: center; color: #94a3b8; padding: 20px 0;",
                   h5("No plates available in registry.")))
      }
      
      # Extract plate info
      plate_info <- lapply(names(plates), function(u) {
        p <- plates[[u]]
        d_val <- safe_extract_plate_date(p@t0)
        d_str <- format(d_val, "%Y-%m-%d")
        d_label <- format(d_val, "%B %d, %Y")
        
        m_str <- if (!is.null(p@media) && !is.na(p@media) && nzchar(p@media)) p@media else "Media"
        t_str <- if (!is.null(p@treatment) && !is.na(p@treatment) && nzchar(p@treatment)) p@treatment else "NoDrug"
        r_str <- if (!is.null(p@replicate) && !is.na(p@replicate)) as.character(p@replicate) else "1"
        s_str <- if (!is.null(p@slot_id) && !is.na(p@slot_id) && nzchar(p@slot_id)) p@slot_id else ""
        
        label <- if (nzchar(s_str)) {
          sprintf("%s: %s + %s (Rep%s)", s_str, m_str, t_str, r_str)
        } else {
          sprintf("%s + %s (Rep%s)", m_str, t_str, r_str)
        }
        
        list(uuid = u, day_str = d_str, day_label = d_label, media = m_str, treatment = t_str, replicate = r_str, slot_id = s_str, label = label)
      })
      
      # Unique date strings sorted newest to oldest (descending order)
      all_dates_str <- unique(sapply(plate_info, `[[`, "day_str"))
      unique_days_str <- sort(all_dates_str, decreasing = TRUE)
      
      selected_uuids <- qc_plates_selected()
      
      # Build date groups
      ui_date_groups <- lapply(unique_days_str, function(d_str) {
        day_plates <- plate_info[sapply(plate_info, `[[`, "day_str") == d_str]
        day_label <- day_plates[[1]]$day_label
        
        # Sort plates within the same date by slot_id, media, treatment, replicate
        day_plates <- day_plates[order(sapply(day_plates, `[[`, "slot_id"),
                                       sapply(day_plates, `[[`, "media"),
                                       sapply(day_plates, `[[`, "treatment"),
                                       sapply(day_plates, `[[`, "replicate"))]
        
        div(class = "qc-date-group",
          div(class = "qc-date-header",
            icon("calendar-alt"),
            day_label
          ),
          div(class = "qc-cards-row",
            lapply(day_plates, function(info) {
              is_sel <- info$uuid %in% selected_uuids
              btn_cls <- if (is_sel) "btn qc-plate-card-btn active" else "btn qc-plate-card-btn"
              ic <- if (is_sel) "check-circle" else "circle"
              
              actionButton(
                ns(paste0("btn_toggle_qc_", info$uuid)),
                label = tagList(
                  icon(ic),
                  span(info$label)
                ),
                class = btn_cls
              )
            })
          )
        )
      })
      
      div(class = "qc-scroll-container",
        do.call(tagList, ui_date_groups)
      )
    })
    
    # Main Panel: Show plate explorers for selected plates in draggable & resizable layout (max 2 per row)
    output$qc_explorers_ui <- renderUI({
      selected <- qc_plates_selected()
      if (length(selected) == 0) {
        return(div(
          style = "text-align: center; color: #94a3b8; padding: 60px 20px; background: #ffffff; border-radius: 8px; border: 2px dashed #cbd5e1; margin-top: 10px;",
          icon("hand-pointer", class = "fa-3x", style = "margin-bottom: 15px; color: #cbd5e1;"),
          h4("No Plate Selected", style = "color: #475569; font-weight: 600; margin-bottom: 8px;"),
          p("Click any compact plate card above to toggle its QC Explorer down below.", style = "color: #94a3b8; max-width: 480px; margin: 0 auto;")
        ))
      }
      
      # Enforce maximum 2 plates per row (calc(50% - 10px) width when 2 or more plates selected)
      card_width_style <- if (length(selected) == 1) {
        "flex: 0 0 100%; width: 100%; max-width: 100%;"
      } else {
        "flex: 0 0 calc(50% - 10px); width: calc(50% - 10px); max-width: calc(50% - 10px);"
      }
      
      explorers <- lapply(selected, function(uuid) {
        div(
          class = "qc-explorer-card",
          style = paste(card_width_style, "min-width: 350px;"),
          mod_plate_explorer_ui(ns(paste0("explorer_", uuid)))
        )
      })
      
      div(class = "qc-explorers-container",
        do.call(tagList, explorers)
      )
    })
    
    # Initialize Plate Explorer server modules on demand
    observe({
      selected <- qc_plates_selected()
      init <- initialized_explorers()
      to_init <- setdiff(selected, init)
      
      for (uuid in to_init) {
        local({
          my_uuid <- uuid
          mod_plate_explorer_server(
            paste0("explorer_", my_uuid),
            plate = reactive({ rv$registry@plates[[my_uuid]] }),
            on_close = function() {
              qc_plates_selected(setdiff(qc_plates_selected(), my_uuid))
            }
          )
        })
      }
      
      if (length(to_init) > 0) {
        initialized_explorers(c(init, to_init))
      }
    })
    
  })
}
