#' @import shiny
#' @noRd
mod_batch_overview_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    fluidRow(
      column(4,
        box(title = "Select Batch", width = 12, status = "primary", solidHeader = TRUE,
            uiOutput(ns("batch_selector")),
            hr(),
            uiOutput(ns("batch_info"))
        )
      ),
      column(8,
        box(title = "Replicate Growth Curves", width = 12, status = "info", solidHeader = TRUE,
            uiOutput(ns("gene_selector")),
            plotly::plotlyOutput(ns("plot_replicate_curves"), height = "400px")
        )
      )
    )
  )
}

#' @import shiny
#' @noRd
mod_batch_overview_server <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    output$batch_selector <- renderUI({
      batches <- rv$registry@batches
      if (length(batches) == 0) {
        selectInput(ns("sel_batch"), "Batch", choices = "No batches available")
      } else {
        selectInput(ns("sel_batch"), "Batch", choices = names(batches))
      }
    })
    
    active_batch <- reactive({
      b_name <- input$sel_batch
      if (is.null(b_name) || b_name == "No batches available") return(NULL)
      rv$registry@batches[[b_name]]
    })
    
    output$batch_info <- renderUI({
      b <- active_batch()
      if (is.null(b)) return(p(style = "color: #888;", "No batch selected."))
      
      uuids <- b@plate_uuids
      plates <- rv$registry@plates[uuids]
      plates <- plates[!sapply(plates, is.null)]
      
      treatments <- paste(unique(sapply(plates, function(p) p@treatment)), collapse = ", ")
      media <- paste(unique(sapply(plates, function(p) p@media)), collapse = ", ")
      
      tagList(
        h4(strong(b@name)),
        p(strong("Total Plates: "), length(plates)),
        p(strong("Media: "), media),
        p(strong("Treatments: "), treatments),
        tags$ul(
          lapply(plates, function(p) {
            tags$li(sprintf("%s (%s + %s + Rep%s)", p@label, p@media, p@treatment, p@replicate))
          })
        )
      )
    })
    
    # Pre-compute the full phenotype data for the active batch
    batch_data <- reactive({
      b <- active_batch()
      if (is.null(b)) return(NULL)
      
      df_list <- list()
      for (u in b@plate_uuids) {
        p <- rv$registry@plates[[u]]
        if (is.null(p)) next
        
        # Prefer corrected, fallback to growth
        df <- p@metrics$corrected
        if (is.null(df)) df <- p@assays$growth
        if (is.null(df)) next
        
        # Ensure layout columns are present
        if (!"Gene" %in% names(df) && nrow(p@layout) > 0) {
          join_cols <- intersect(names(df), c("Row", "Column"))
          cols_to_add <- setdiff(names(p@layout), names(df))
          if (length(cols_to_add) > 0) {
            df <- df %>% dplyr::left_join(p@layout[, c(join_cols, cols_to_add), drop = FALSE], by = join_cols)
          }
        }
        
        df$Batch <- b@name
        df$Plate_UUID <- u
        df$Rep_Suffix <- paste0("P", which(b@plate_uuids == u))
        
        df_list[[u]] <- df
      }
      if (length(df_list) == 0) return(NULL)
      dplyr::bind_rows(df_list)
    })
    
    output$gene_selector <- renderUI({
      df <- batch_data()
      if (is.null(df)) {
        return(selectInput(ns("sel_gene"), "Select Mutant/Gene", choices = NULL))
      }
      
      genes <- sort(unique(df$Gene[!is.na(df$Gene) & nzchar(df$Gene) & df$Type != "Blank"]))
      selectInput(ns("sel_gene"), "Select Mutant/Gene", choices = genes)
    })
    
    output$plot_replicate_curves <- plotly::renderPlotly({
      df <- batch_data()
      g <- input$sel_gene
      if (is.null(df) || is.null(g) || g == "") {
        return(plotly::plot_ly() %>% plotly::add_annotations(text="No data available", showarrow=FALSE))
      }
      
      sub_df <- df %>% dplyr::filter(Gene == g)
      
      plot_replicate_growth_curves(sub_df, stain_hr = 26)
    })
    
  })
}
