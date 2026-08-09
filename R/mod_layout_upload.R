#' @import shiny
#' @noRd
mod_layout_upload_ui <- function(id) {
  # The UI for the modal content can be rendered dynamically in the server
  # But we can place a hidden div or just rely on showModal
  ns <- NS(id)
  tagList()
}

#' @import shiny
#' @noRd
mod_layout_upload_server <- function(id, trigger, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Internal state for the modal
    modal_state <- reactiveValues(
      parsed_data = NULL,
      error_msg = NULL,
      success = FALSE
    )
    
    observeEvent(trigger(), {
      req(trigger() > 0)
      modal_state$parsed_data <- NULL
      modal_state$error_msg <- NULL
      modal_state$success <- FALSE
      
      showModal(
        modalDialog(
          title = "Upload Layout",
          uiOutput(ns("modal_content")),
          footer = uiOutput(ns("modal_footer")),
          size = "m",
          easyClose = FALSE
        )
      )
    })
    
    output$modal_content <- renderUI({
      if (modal_state$success && !is.null(modal_state$parsed_data)) {
        # Success state
        p_data <- modal_state$parsed_data
        
        wt_count <- sum(p_data$data$Type == "WT_Control", na.rm = TRUE)
        blank_count <- sum(p_data$data$Type == "Blank", na.rm = TRUE)
        mutant_count <- sum(p_data$data$Type == "Mutant", na.rm = TRUE)
        
        tagList(
          div(style = "color: green; font-weight: bold; margin-bottom: 10px;", "\u2713 PARSE SUCCESSFUL"),
          p(sprintf("Format detected: %s-well (%d rows x %d cols)", 
                    p_data$format, p_data$n_rows, p_data$n_cols)),
          p(sprintf("Wells mapped: %d   Genes: %d mutants, %d WT, %d Blank", 
                    p_data$n_wells, mutant_count, wt_count, blank_count)),
          hr(),
          textInput(ns("layout_name"), "Name this layout: (required, unique)", placeholder = "e.g., Lib384_v1"),
          textInput(ns("layout_notes"), "Notes (optional):", placeholder = "e.g., Standard library")
        )
      } else {
        # Upload / Failure state
        tagList(
          if (!is.null(modal_state$error_msg)) {
            div(style = "color: red; font-weight: bold; margin-bottom: 15px;",
                " PARSE FAILED",
                p(modal_state$error_msg, style = "font-weight: normal; margin-top: 5px;")
            )
          },
          fileInput(ns("layout_file"), "Select layout file (.xlsx or .tsv):", 
                    accept = c(".xlsx", ".tsv", ".xls")),
          actionButton(ns("btn_parse"), "Upload & Parse", class = "btn-primary")
        )
      }
    })
    
    output$modal_footer <- renderUI({
      if (modal_state$success) {
        tagList(
          modalButton("Cancel"),
          actionButton(ns("btn_save"), "Save Layout", class = "btn-success")
        )
      } else {
        modalButton("Cancel")
      }
    })
    
    observeEvent(input$btn_parse, {
      req(input$layout_file)
      
      tryCatch({
        res <- parse_plate_layout(input$layout_file$datapath)
        modal_state$parsed_data <- res
        modal_state$parsed_data$filename <- input$layout_file$name
        modal_state$success <- TRUE
        modal_state$error_msg <- NULL
      }, error = function(e) {
        modal_state$success <- FALSE
        modal_state$parsed_data <- NULL
        modal_state$error_msg <- e$message
      })
    })
    
    observeEvent(input$btn_save, {
      req(input$layout_name)
      
      lname <- trimws(input$layout_name)
      if (lname == "") {
        shinyjs::alert("Please provide a name for the layout.")
        return()
      }
      
      if (lname %in% names(rv$registry@layouts)) {
        shinyjs::alert("A layout with this name already exists. Please choose a different name.")
        return()
      }
      
      layout_entry <- modal_state$parsed_data
      layout_entry$notes <- trimws(input$layout_notes)
      
      rv$registry <- add_layout(rv$registry, lname, layout_entry)
      
      # Auto-save
      save_registry(rv$registry)
      
      removeModal()
      
      # Inform user
      shinyjs::info("Layout saved successfully.")
    })
    
  })
}
