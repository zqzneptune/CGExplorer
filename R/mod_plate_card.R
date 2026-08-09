#' @import shiny
#' @noRd
mod_plate_card_ui <- function(id, plate, selected = FALSE) {
  ns <- NS(id)
  
  # Format dates
  t0_str <- format(plate@t0, "%Y-%m-%d %H:%M")
  tend_str <- format(plate@t_end, "%H:%M")
  
  # Status badges
  has_staining <- !is.null(plate@assays$staining)
  staining_badge <- if (has_staining && !is.na(plate@staining_hr)) {
    span(class = "label label-success", sprintf("Stained at %.1f hr", plate@staining_hr))
  } else if (has_staining) {
    span(class = "label label-success", "Stained \u2713")
  } else {
    span(class = "label label-warning", "No Staining")
  }
  
  metrics_badge <- if (length(plate@metrics) > 0) {
    span(class = "label label-primary", "\u2713 Computed")
  } else {
    span(class = "label label-default", "\u23f3 Pending")
  }
  
  qc_flags <- plate@qc_flags
  blanks_badge <- if (!is.null(qc_flags$removed_blanks) && nrow(qc_flags$removed_blanks) > 0) {
    span(class = "label label-danger", sprintf(" %d Flagged", nrow(qc_flags$removed_blanks)))
  } else {
    span(class = "label label-success", "Blanks OK")
  }
  
  merged_tag <- if (plate@is_merged) span(class = "label label-info", "Merged") else NULL
  
  plate_label <- sprintf("%s + %s + Rep%s", plate@media, plate@treatment, plate@replicate)
  
  div(class = "box box-solid box-primary", style = "margin-bottom: 15px; width: 320px; display: inline-block; vertical-align: top; margin-right: 15px;",
      div(class = "box-header with-border", style = "padding: 8px 10px;",
          h3(class = "box-title", 
             style = "display: flex; align-items: center; gap: 8px; font-size: 14px; margin: 0;",
             checkboxInput(ns("chk_merge"), "", value = selected, width = "20px"),
             span(strong(plate_label), span(style = "color:#888; font-size:12px;", sprintf("(%s)", plate@slot_id)), merged_tag))
      ),
      div(class = "box-body", style = "padding: 10px;",
          p(strong("Time: "), t0_str, " \u2192 ", tend_str, style = "margin: 0; font-size: 13px;"),
          p(strong("Staining: "), staining_badge, style = "margin: 0; font-size: 13px;")
      ),
      div(class = "box-footer", style = "padding: 8px 10px;",
          actionButton(ns("btn_edit"), "\u270f\ufe0f Edit", class = "btn-xs btn-default"),
          actionButton(ns("btn_delete"), "\U0001f5d1 Delete", class = "btn-xs btn-danger pull-right")
      )
  )
}

#' @import shiny
#' @noRd
mod_plate_card_server <- function(id, plate, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    observeEvent(input$btn_delete, {
      plate_disp <- if (is.na(plate@label) || plate@label == "") plate@slot_id else plate@label
      showModal(modalDialog(
        title = "Delete Plate",
        sprintf("Are you sure you want to delete plate %s (%s)?", plate_disp, plate@uuid),
        footer = tagList(
          modalButton("Cancel"),
          actionButton(ns("btn_confirm_delete"), "Yes, Delete", class = "btn-danger")
        )
      ))
    })
    
    observeEvent(input$btn_confirm_delete, {
      rv$registry <- remove_plate(rv$registry, plate@uuid)
      save_registry(rv$registry)
      removeModal()
      shinyjs::info("Plate deleted.")
    })
    
    observeEvent(input$btn_edit, {
      current_p <- rv$registry@plates[[plate@uuid]]
      if (is.null(current_p)) return()
      
      known_media <- character(0)
      known_treat <- character(0)
      
      if (length(rv$registry@plates) > 0) {
        known_media <- unique(sapply(rv$registry@plates, function(p) p@media))
        known_treat <- unique(sapply(rv$registry@plates, function(p) p@treatment))
      }
      
      machine_dir <- "machine"
      if (dir.exists(machine_dir)) {
        csv_files <- list.files(machine_dir, pattern = "\\.csv$", full.names = TRUE)
        for (f in csv_files) {
          try({
            df <- utils::read.csv(f, stringsAsFactors = FALSE)
            if ("Media" %in% colnames(df)) known_media <- c(known_media, df$Media)
            if ("Drug" %in% colnames(df)) known_treat <- c(known_treat, df$Drug)
            if ("Treatment" %in% colnames(df)) known_treat <- c(known_treat, df$Treatment)
          }, silent = TRUE)
        }
      }
      
      known_media <- unique(stats::na.omit(known_media))
      known_treat <- unique(stats::na.omit(known_treat))
      
      known_media <- unique(c("SCFM", "AUM", "CFA", known_media))
      known_treat <- unique(c("NoDrug", "Duloxetine", "Mirtazapine", "Proflavine", "Sertraline", known_treat))
      
      plate_disp <- if (is.na(current_p@label) || current_p@label == "") current_p@slot_id else current_p@label
      
      showModal(modalDialog(
        title = paste("Edit Plate:", plate_disp),
        fluidPage(
          fluidRow(
            column(6, selectizeInput(ns("edit_media"), "Media", choices = known_media, selected = current_p@media, options = list(create = TRUE))),
            column(6, selectizeInput(ns("edit_treat"), "Treatment", choices = known_treat, selected = current_p@treatment, options = list(create = TRUE)))
          ),
          fluidRow(
            column(6, numericInput(ns("edit_rep"), "Replicate", value = current_p@replicate, min = 1)),
            column(6, 
                   div(style = "margin-bottom: 0;",
                       radioButtons(ns("edit_stain_status"), "Staining?", 
                                    choices = c("Yes" = "yes", "No Staining" = "no"), 
                                    inline = TRUE, selected = if(!is.na(current_p@staining_hr)) "yes" else "no")
                   ),
                   conditionalPanel(
                     condition = sprintf("input['%s'] == 'yes'", ns("edit_stain_status")),
                     div(style = "display: flex; align-items: center; gap: 8px;",
                       span("Staining hour:"),
                       numericInput(ns("edit_stain_hr"), NULL, 
                                    value = if(!is.na(current_p@staining_hr)) current_p@staining_hr else 24.0, 
                                    min = 0, step = 0.1, width = "80px")
                     )
                   )
            )
          )
        ),
        footer = tagList(
          modalButton("Cancel"),
          actionButton(ns("btn_save_edit"), "Save Changes", class = "btn-primary")
        )
      ))
    })
    
    observeEvent(input$btn_save_edit, {
      current_p <- rv$registry@plates[[plate@uuid]]
      if (is.null(current_p)) return()
      
      current_p@media <- input$edit_media
      current_p@treatment <- input$edit_treat
      current_p@replicate <- as.integer(input$edit_rep)
      
      if (input$edit_stain_status == "yes" && !is.na(input$edit_stain_hr)) {
        current_p@staining_hr <- input$edit_stain_hr
      } else {
        current_p@staining_hr <- as.numeric(NA)
      }
      
      rv$registry <- update_plate(rv$registry, current_p)
      save_registry(rv$registry)
      removeModal()
      shinyjs::info("Plate metadata updated.")
    })
    
    # Manage selection state
    observeEvent(input$chk_merge, {
      if (input$chk_merge) {
        rv$selected_plates <- unique(c(rv$selected_plates, plate@uuid))
      } else {
        rv$selected_plates <- setdiff(rv$selected_plates, plate@uuid)
      }
    }, ignoreInit = TRUE)
  })
}
