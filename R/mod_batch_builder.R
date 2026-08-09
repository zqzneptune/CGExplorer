#' @import shiny
#' @noRd
mod_batch_builder_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    fluidRow(
      column(4,
        box(title = "Create New Batch", width = 12, status = "success", solidHeader = TRUE,
            textInput(ns("batch_name"), "Batch Name", placeholder = "e.g. Aug_Experiment"),
            uiOutput(ns("selected_plates_ui")),
            hr(),
            actionButton(ns("btn_create"), "\u2713 Create Batch", class = "btn-success", width = "100%")
        ),
        box(title = "Existing Batches", width = 12, status = "info", solidHeader = TRUE,
            uiOutput(ns("existing_batches_ui"))
        )
      ),
      column(8,
        box(title = "Available Plates", width = 12, status = "primary", solidHeader = TRUE,
            DT::dataTableOutput(ns("dt_plates"))
        )
      )
    )
  )
}

#' @import shiny
#' @noRd
mod_batch_builder_server <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Reactive to get ordered UUIDs for consistent indexing
    ordered_uuids <- reactive({
      plates <- rv$registry@plates
      if (length(plates) == 0) return(character(0))
      
      t0 <- sapply(plates, function(p) p@t0)
      media <- sapply(plates, function(p) p@media)
      treatment <- sapply(plates, function(p) p@treatment)
      replicate <- sapply(plates, function(p) p@replicate)
      
      ord <- order(t0, media, treatment, replicate, na.last = TRUE)
      names(plates)[ord]
    })
    
    # 1. Available Plates Table
    output$dt_plates <- DT::renderDataTable({
      uuids <- ordered_uuids()
      if (length(uuids) == 0) return(data.frame(Message = "No plates available."))
      
      plates <- rv$registry@plates[uuids]
      
      df <- data.frame(
        ` ` = rep("", length(uuids)),
        `Start Time` = sapply(plates, function(p) format(p@t0, "%Y-%m-%d %H:%M")),
        Treatment = sapply(plates, function(p) p@treatment),
        Media = sapply(plates, function(p) p@media),
        Replicate = sapply(plates, function(p) p@replicate),
        Label = sapply(plates, function(p) p@label),
        Slot_ID = sapply(plates, function(p) p@slot_id),
        UUID = uuids,
        stringsAsFactors = FALSE, check.names = FALSE
      )
      
      DT::datatable(df, selection = 'none', rownames = FALSE, 
                    extensions = 'Select',
                    options = list(
                      pageLength = 10, dom = 'ftip',
                      columnDefs = list(
                        list(orderable = FALSE, className = 'select-checkbox', targets = 0),
                        list(visible = FALSE, targets = c(6, 7))
                      ),
                      select = list(style = 'multi', selector = 'td:first-child')
                    ))
    }, server = FALSE)
    
    # 2. Selected Plates preview
    output$selected_plates_ui <- renderUI({
      sel <- input$dt_plates_rows_selected
      if (is.null(sel) || length(sel) == 0) {
        p(style = "color: #888;", "Select >= 2 plates from the table to form a batch.")
      } else {
        plates <- rv$registry@plates
        uuids <- ordered_uuids()[sel]
        tagList(
          p(strong(sprintf("Selected %d plates:", length(uuids)))),
          tags$ul(
            lapply(uuids, function(u) {
              tags$li(sprintf("%s + %s + Rep%s", plates[[u]]@media, plates[[u]]@treatment, plates[[u]]@replicate))
            })
          )
        )
      }
    })
    
    # 3. Existing Batches List UI
    output$existing_batches_ui <- renderUI({
      batches <- rv$registry@batches
      if (length(batches) == 0) {
        return(p(style = "color: #888;", "No batches created yet."))
      }
      
      items <- lapply(names(batches), function(b_name) {
        b <- batches[[b_name]]
        n_p <- length(b@plate_uuids)
        div(style = "display: flex; justify-content: space-between; align-items: center; padding: 8px 0; border-bottom: 1px solid #eee;",
          div(
            strong(b_name),
            br(),
            span(style = "color: #666; font-size: 0.9em;", sprintf("%d plates", n_p))
          ),
          actionButton(ns(paste0("btn_del_", b_name)), "Delete", class = "btn-danger btn-xs", icon = icon("trash"))
        )
      })
      do.call(tagList, items)
    })
    
    # Observe dynamic delete button clicks for existing batches
    observe({
      batches <- rv$registry@batches
      for (b_name in names(batches)) {
        local({
          name <- b_name
          btn_id <- paste0("btn_del_", name)
          observeEvent(input[[btn_id]], {
            rv$registry <- remove_batch(rv$registry, name)
            save_registry(rv$registry)
            shinyjs::info(sprintf("Batch '%s' deleted.", name))
          }, ignoreInit = TRUE, once = TRUE)
        })
      }
    })
    
    # 4. Create Batch
    observeEvent(input$btn_create, {
      req(input$batch_name)
      sel <- input$dt_plates_rows_selected
      if (is.null(sel) || length(sel) < 2) {
        shinyjs::alert("Please select at least 2 plates.")
        return()
      }
      
      uuids <- ordered_uuids()[sel]
      
      b <- new_batch(name = input$batch_name, plate_uuids = uuids)
      rv$registry <- add_batch(rv$registry, b)
      save_registry(rv$registry)
      
      shinyjs::info("Batch created successfully.")
      updateTextInput(session, "batch_name", value = "")
      DT::dataTableProxy("dt_plates") %>% DT::selectRows(NULL)
    })
    
  })
}
