#' The application User-Interface
#'
#' @import shiny
#' @import shinydashboard
#' @noRd
app_ui <- function(request) {
  tagList(
    shinyjs::useShinyjs(),
    dashboardPage(
      dashboardHeader(title = "CGExplorer v0.3"),
      dashboardSidebar(
        sidebarMenu(
          id = "tabs",
          menuItem("Plate Inventory", tabName = "plates", icon = icon("home")),
          menuItem("Batch Builder", tabName = "batch_builder", icon = icon("cogs")),
          # menuItem("Batch Overview", tabName = "batch_overview", icon = icon("microscope")),
          # menuItem("Scoring", tabName = "scoring", icon = icon("chart-bar")),
          menuItem("Plate QC", tabName = "qc", icon = icon("search")),
          menuItem("About", tabName = "about", icon = icon("info-circle"))
        ),
        hr(),
        uiOutput("sidebar_project_info")
      ),
      dashboardBody(
        tabItems(
          tabItem(tabName = "plates",
            fluidRow(
              column(12,
                div(style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 20px;",
                  div(
                    selectInput("plate_sort", NULL, choices = c("Oldest first", "Newest first"), width = "150px")
                  ),
                  div(textOutput("plate_count_text"))
                ),
                uiOutput("plate_cards_ui"),
                hr(),
                div(style = "text-align: center; margin-top: 20px;",
                  actionButton("btn_upload_data", " Upload Data", class = "btn-primary btn-lg", style = "margin-right: 15px;"),
                  actionButton("btn_upload_layout", "\U0001f4d0 Upload Layout", class = "btn-info btn-lg")
                ),
                uiOutput("floating_merge_bar")
              )
            )
          ),
          tabItem(tabName = "batch_builder", mod_batch_builder_ui("batch_builder")),
          tabItem(tabName = "batch_overview", mod_batch_overview_ui("batch_overview")),
          tabItem(tabName = "scoring", mod_scoring_ui("scoring")),
          tabItem(tabName = "qc", mod_qc_ui("qc")),
          tabItem(tabName = "about", 
            div(style = "padding: 20px;",
              h2("About CGExplorer"),
              p("CGExplorer is a Shiny application for preprocessing and visualizing Chemo-Genomic Screening dataset."),
              p("Features include standard pipeline processing, batch analysis, interaction scoring, and interactive visualizations.")
            )
          )
        )
      )
    )
  )
}
