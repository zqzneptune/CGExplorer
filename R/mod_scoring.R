#' @import shiny
#' @noRd
mod_scoring_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    fluidRow(
      column(4,
        box(title = "Scoring Parameters", width = 12, status = "primary", solidHeader = TRUE,
            uiOutput(ns("batch_selector")),
            selectInput(ns("strategy"), "Strategy", 
                        choices = c("NoDrug Plates (Inter-plate)" = "NoDrug_Plates", 
                                    "NoDrug Wells (Intra-plate)" = "NoDrug_Wells")),
            selectInput(ns("metric"), "Metric to Score", 
                        choices = c("Growth_AUC", "Growth_AUC_Rate", "Growth_Endpoint_OD", 
                                    "Biofilm_Stain_OD", "Biofilm_Score_RB")),
            textInput(ns("control_drug"), "Control Drug Name", value = "NoDrug"),
            hr(),
            actionButton(ns("btn_score"), "\u25b6 Run Scoring", class = "btn-success btn-lg", width = "100%")
        )
      ),
      column(8,
        box(title = "Scoring Results", width = 12, status = "info", solidHeader = TRUE,
            div(style = "display: flex; justify-content: space-between; align-items: center;",
              uiOutput(ns("results_summary")),
              downloadButton(ns("btn_download_csv"), "Download CSV", class = "btn-primary btn-sm")
            ),
            hr(),
            DT::dataTableOutput(ns("dt_scores"))
        )
      )
    )
  )
}

#' @import shiny
#' @noRd
mod_scoring_server <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    output$batch_selector <- renderUI({
      batches <- rv$registry@batches
      if (length(batches) == 0) {
        selectInput(ns("selected_batch"), "Select Batch", choices = "No batches available")
      } else {
        selectInput(ns("selected_batch"), "Select Batch", choices = names(batches))
      }
    })
    
    observeEvent(input$btn_score, {
      b_name <- input$selected_batch
      if (is.null(b_name) || b_name == "No batches available") {
        shinyjs::alert("Please create and select a batch first.")
        return()
      }
      
      batch <- rv$registry@batches[[b_name]]
      
      shinyjs::disable("btn_score")
      id_msg <- showNotification("Scoring batch...", duration = NULL, type = "message")
      
      tryCatch({
        b_scored <- score_batch(batch, rv$registry, 
                                strategies = c(input$strategy), 
                                metrics = c(input$metric), 
                                control_drug = input$control_drug)
                                
        rv$registry <- add_batch(rv$registry, b_scored)
        save_registry(rv$registry)
        
        showNotification("Scoring completed.", type = "message")
      }, error = function(e) {
        showNotification(paste("Scoring failed:", e$message), type = "error")
      }, finally = {
        removeNotification(id_msg)
        shinyjs::enable("btn_score")
      })
    })
    
    scored_data <- reactive({
      b_name <- input$selected_batch
      if (is.null(b_name) || b_name == "No batches available") return(NULL)
      
      batch <- rv$registry@batches[[b_name]]
      get_scores(batch, strategy = input$strategy, metric = input$metric)
    })
    
    output$results_summary <- renderUI({
      df <- scored_data()
      if (is.null(df) || nrow(df) == 0) return(p("No scoring results for this combination. Click 'Run Scoring'."))
      
      n_total <- nrow(df)
      n_sig <- sum(df$Is_Significant == "Y", na.rm = TRUE)
      
      tagList(
        p(strong("Total Interactions Scored: "), n_total, style = "margin-bottom: 2px;"),
        p(strong("Significant Hits (FDR < 0.05, |log2CG| > 1): "), n_sig, style = "margin-bottom: 0;")
      )
    })
    
    output$btn_download_csv <- downloadHandler(
      filename = function() {
        paste0("CGExplorer_Scores_", input$selected_batch, "_", input$strategy, "_", input$metric, ".csv")
      },
      content = function(file) {
        df <- scored_data()
        if (!is.null(df) && nrow(df) > 0) {
          readr::write_csv(df, file)
        }
      }
    )
    
    output$dt_scores <- DT::renderDataTable({
      df <- scored_data()
      if (is.null(df) || nrow(df) == 0) return(NULL)
      
      display_df <- df %>%
        dplyr::select(Treatment, Gene, Replicate_Pair, cg_score, log_cg_score, adj_p, Is_Significant) %>%
        dplyr::mutate(
          cg_score = round(cg_score, 3),
          log_cg_score = round(log_cg_score, 3),
          adj_p = signif(adj_p, 3)
        )
        
      DT::datatable(display_df, filter = "top", rownames = FALSE,
                    options = list(pageLength = 15, scrollX = TRUE))
    })
    
  })
}
