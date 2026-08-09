#' The application server-side
#'
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {
  
  rv <- reactiveValues(
    registry = new_registry(project_name = "Default Project"),
    initialized_plates = character(0),
    selected_plates = character(0)
  )
  
  # Try to load existing registry on startup
  observe({
    tryCatch({
      loaded <- load_registry(getwd())
      rv$registry <- loaded
    }, error = function(e) {
      # No existing registry, keep default
    })
  })
  
  output$sidebar_project_info <- renderUI({
    tagList(
      p(paste("Project:", rv$registry@project_name), style = "padding: 10px;"),
      p(paste("Plates:", length(rv$registry@plates), "| Batches:", length(rv$registry@batches), "| Layouts:", length(rv$registry@layouts)), style = "padding: 10px;")
    )
  })
  
  output$plate_count_text <- renderText({
    paste("Showing", length(rv$registry@plates), "plates")
  })
  
  output$plate_cards_ui <- renderUI({
    plates <- rv$registry@plates
    if (length(plates) == 0) {
      div(style = "text-align: center; color: #888; margin: 50px 0;",
        h3("(empty - no plates yet)")
      )
    } else {
      # Group and sort plates
      plate_info <- lapply(names(plates), function(u) {
        p <- plates[[u]]
        list(uuid = u, day = as.Date(p@t0), slot_id = p@slot_id, t0 = p@t0)
      })
      
      # Order by day ascending, then slot_id ascending
      plate_info <- plate_info[order(as.Date(sapply(plate_info, `[[`, "day")), sapply(plate_info, `[[`, "slot_id"))]
      
      # Determine day order based on sort toggle
      unique_days <- unique(sapply(plate_info, function(x) as.character(x$day)))
      if (req(input$plate_sort) == "Newest first") {
        unique_days <- rev(unique_days)
      }
      
      ui_groups <- lapply(unique_days, function(d) {
        # Get uuids for this day
        day_uuids <- sapply(plate_info[sapply(plate_info, function(x) as.character(x$day)) == d], function(x) x$uuid)
        
        # Build cards for this day
        cards <- lapply(day_uuids, function(uuid) {
          if (!uuid %in% rv$initialized_plates) {
            mod_plate_card_server(paste0("plate_", uuid), plate = plates[[uuid]], rv = rv)
            rv$initialized_plates <- c(rv$initialized_plates, uuid)
          }
          is_selected <- uuid %in% rv$selected_plates
          mod_plate_card_ui(paste0("plate_", uuid), plate = plates[[uuid]], selected = is_selected)
        })
        
        # Wrapper for the day group
        div(style = "margin-bottom: 30px;",
            h4(strong(format(as.Date(d), "%B %d, %Y")), style = "border-bottom: 2px solid #3c8dbc; padding-bottom: 5px; color: #3c8dbc;"),
            div(style = "display: flex; flex-wrap: wrap;",
                do.call(tagList, cards)
            )
        )
      })
      
      do.call(tagList, ui_groups)
    }
  })
  
  output$floating_merge_bar <- renderUI({
    req(length(rv$selected_plates) >= 2)
    div(
      style = "position: fixed; bottom: 20px; left: 50%; transform: translateX(-50%); background: #3c8dbc; padding: 15px; border-radius: 8px; box-shadow: 0 4px 8px rgba(0,0,0,0.3); z-index: 9999;",
      actionButton("btn_merge_selected", sprintf("Merge Selected Plates (%d)", length(rv$selected_plates)), 
                   class = "btn-warning btn-lg", icon = icon("object-group"))
    )
  })
  
  observeEvent(input$btn_merge_selected, {
    showModal(modalDialog(
      title = "Merge Plates",
      "Are you sure you want to merge these plates? The original plates will be removed.",
      footer = tagList(
        modalButton("Cancel"),
        actionButton("btn_confirm_merge", "Merge", class = "btn-warning")
      )
    ))
  })
  
  observeEvent(input$btn_confirm_merge, {
    tryCatch({
      rv$registry <- merge_plates(rv$registry, rv$selected_plates)
      rv$selected_plates <- character(0)
      save_registry(rv$registry)
      removeModal()
      shinyjs::info("Plates successfully merged.")
    }, error = function(e) {
      shinyjs::alert(paste("Merge failed:", e$message))
    })
  })
  
  mod_data_upload_server("data_upload", trigger = reactive(input$btn_upload_data), rv = rv)
  
  mod_layout_upload_server("layout_upload", trigger = reactive(input$btn_upload_layout), rv = rv)
  
  mod_batch_builder_server("batch_builder", rv = rv)
  
  mod_batch_overview_server("batch_overview", rv = rv)
  
  mod_scoring_server("scoring", rv = rv)
  
  mod_qc_server("qc", rv = rv)
}
