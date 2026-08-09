#' @import shiny
#' @noRd
mod_data_upload_ui <- function(id) {
  ns <- NS(id)
  tagList()
}

#' @import shiny
#' @noRd
mod_data_upload_server <- function(id, trigger, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Steps: 1 = File Select, 2 = Parsing, 3 = Annotation
    wizard_state <- reactiveValues(
      step = 1,
      files = NULL,
      parsed_data = NULL,
      slots_found = character(0),
      is_parsing = FALSE
    )
    
    observeEvent(trigger(), {
      req(trigger() > 0)
      wizard_state$step <- 1
      wizard_state$files <- NULL
      wizard_state$parsed_data <- NULL
      wizard_state$slots_found <- character(0)
      wizard_state$is_parsing <- FALSE
      
      showModal(
        modalDialog(
          title = "Upload Machine Data",
          uiOutput(ns("wizard_content")),
          footer = uiOutput(ns("wizard_footer")),
          size = "l",
          easyClose = FALSE
        )
      )
    })
    
    output$wizard_content <- renderUI({
      if (wizard_state$step == 1) {
        tagList(
          h4("Step 1: Select Files"),
          p("Upload raw Tecan machine .txt files."),
          fileInput(ns("raw_files"), "Choose Files", multiple = TRUE, accept = ".txt")
        )
      } else if (wizard_state$step == 2) {
        tagList(
          h4("Step 2: Parsing Progress"),
          if (wizard_state$is_parsing) {
            div(style = "text-align: center; margin: 20px;",
                shinycssloaders::withSpinner(uiOutput(ns("parsing_spinner"))),
                p("Parsing files...")
            )
          } else {
            div(style = "color: green; font-weight: bold; margin-bottom: 10px;",
                "\u2713 PARSE SUCCESSFUL",
                p(sprintf("Slots found: %d", length(wizard_state$slots_found)))
            )
          }
        )
      } else if (wizard_state$step == 3) {
        tagList(
          h4("Step 3: Slot Annotation"),
          if (length(rv$registry@layouts) == 0) {
            div(style = "color: red; font-weight: bold; padding: 10px; border: 1px solid red;",
                "(No layouts found - please close this and click '\U0001f4d0 Upload Layout' to add one first)")
          } else {
            selectInput(ns("selected_layout"), "Layout to apply:", 
                        choices = names(rv$registry@layouts))
          },
          hr(),
          uiOutput(ns("slot_annotations_ui"))
        )
      }
    })
    
    output$parsing_spinner <- renderUI({ "" })
    
    output$slot_annotations_ui <- renderUI({
      req(wizard_state$slots_found)
      
      # Extract historical choices from registry and machine directory
      known_media <- character(0)
      known_treat <- character(0)
      
      # 1. From registry
      if (length(rv$registry@plates) > 0) {
        known_media <- c(known_media, sapply(rv$registry@plates, function(p) p@media))
        known_treat <- c(known_treat, sapply(rv$registry@plates, function(p) p@treatment))
      }
      
      # 2. From machine directory (legacy)
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
      
      # Combine dynamically learned options with the default requested options
      known_media <- unique(c("SCFM", "AUM", "CFA", known_media))
      known_treat <- unique(c("NoDrug", "Duloxetine", "Mirtazapine", "Proflavine", "Sertraline", known_treat))
      
      rows <- lapply(wizard_state$slots_found, function(slot) {
        slot_data <- wizard_state$parsed_data %>% dplyr::filter(Plate == slot)
        stain_info <- detect_staining_jump(slot_data)
        
        stain_val <- if (stain_info$has_staining) round(stain_info$staining_hr, 1) else 24.0
        stain_label <- if (stain_info$has_staining) {
          sprintf("\u26a1 Detected jump at %s", format(stain_info$jump_time, "%H:%M"))
        } else {
          ""
        }
        
        stain_ui <- div(
          style = "margin-top: 5px; margin-bottom: 15px; padding-left: 15px; border-left: 3px solid #f39c12; display: flex; align-items: center; gap: 15px;",
          div(style = "margin-bottom: 0;",
            radioButtons(ns(paste0("stain_status_", slot)), "Staining?", 
                         choices = c("Yes" = "yes", "No Staining" = "no"), 
                         inline = TRUE, selected = if (stain_info$has_staining) "yes" else "no")
          ),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'yes'", ns(paste0("stain_status_", slot))),
            div(style = "display: flex; align-items: center; gap: 8px; margin-top: 25px;",
              span("Staining hour:"),
              numericInput(ns(paste0("stain_hr_", slot)), NULL, 
                           value = stain_val, min = 0, step = 0.1, width = "80px"),
              span(stain_label, style = "color: #888; font-size: 13px;")
            )
          )
        )
        
        div(
          fluidRow(
            style = "align-items: center;",
            column(2, strong(slot)),
            column(3, selectizeInput(ns(paste0("media_", slot)), "Media", choices = known_media, options = list(create = TRUE))),
            column(3, selectizeInput(ns(paste0("treat_", slot)), "Treatment", choices = known_treat, options = list(create = TRUE))),
            column(2, numericInput(ns(paste0("rep_", slot)), "Replicate", value = 1, min = 1)),
            column(2, textInput(ns(paste0("label_", slot)), "Label", placeholder = "Optional"))
          ),
          stain_ui,
          hr(style = "margin-top: 5px; margin-bottom: 5px;")
        )
      })
      do.call(tagList, rows)
    })
    
    output$wizard_footer <- renderUI({
      if (wizard_state$step == 1) {
        tagList(
          modalButton("Cancel"),
          actionButton(ns("btn_next1"), "Next \u2192", class = "btn-primary")
        )
      } else if (wizard_state$step == 2) {
        if (!wizard_state$is_parsing) {
          tagList(
            modalButton("Cancel"),
            actionButton(ns("btn_next2"), "Next \u2192", class = "btn-primary")
          )
        } else {
          modalButton("Cancel")
        }
      } else if (wizard_state$step == 3) {
        can_create <- length(rv$registry@layouts) > 0
        tagList(
          modalButton("Cancel"),
          actionButton(ns("btn_create"), "\u2713 Create Plates", class = "btn-success", 
                       disabled = if (!can_create) TRUE else NULL)
        )
      }
    })
    
    observeEvent(input$btn_next1, {
      req(input$raw_files)
      wizard_state$files <- input$raw_files
      wizard_state$step <- 2
      wizard_state$is_parsing <- TRUE
      
      # Mock parsing delay to show progress, then actually parse
      shinyjs::delay(500, {
        parsed_list <- list()
        slots <- c()
        for (i in 1:nrow(wizard_state$files)) {
          df <- parse_tecan_txt(wizard_state$files$datapath[i])
          if (nrow(df) > 0) {
            parsed_list[[i]] <- df
            slots <- unique(c(slots, unique(df$Plate)))
          }
        }
        wizard_state$parsed_data <- dplyr::bind_rows(parsed_list)
        # Sort slots by S and L if they follow the 'SxLy' pattern
        sort_slots <- function(s) {
          s_num <- suppressWarnings(as.numeric(gsub("^S(\\d+)L\\d+$", "\\1", s)))
          l_num <- suppressWarnings(as.numeric(gsub("^S\\d+L(\\d+)$", "\\1", s)))
          s[order(s_num, l_num, s, na.last = TRUE)]
        }
        wizard_state$slots_found <- sort_slots(slots)
        wizard_state$is_parsing <- FALSE
      })
    })
    
    observeEvent(input$btn_next2, {
      wizard_state$step <- 3
    })
    
    observeEvent(input$btn_create, {
      req(input$selected_layout)
      
      layout_name <- input$selected_layout
      layout_data <- rv$registry@layouts[[layout_name]]$data
      
      for (slot_id in wizard_state$slots_found) {
        treat_val <- input[[paste0("treat_", slot_id)]]
        rep_val <- input[[paste0("rep_", slot_id)]]
        media_val <- input[[paste0("media_", slot_id)]]
        label_val <- input[[paste0("label_", slot_id)]]
        
        slot_data <- wizard_state$parsed_data %>%
          dplyr::filter(Plate == slot_id)
        
        # Staining logic
        stain_status <- input[[paste0("stain_status_", slot_id)]]
        stain_hr_val <- input[[paste0("stain_hr_", slot_id)]]
        
        growth_df <- slot_data
        stain_df <- NULL
        stain_hr_final <- as.numeric(NA)
        stain_confirmed_val <- FALSE
        
        if (!is.null(stain_status) && stain_status == "yes" && !is.null(stain_hr_val)) {
          t0 <- min(slot_data$DateTime, na.rm = TRUE)
          stain_time <- t0 + as.difftime(stain_hr_val, units = "hours")
          growth_df <- slot_data %>% dplyr::filter(DateTime < stain_time)
          stain_df <- slot_data %>% dplyr::filter(DateTime >= stain_time)
          stain_hr_final <- stain_hr_val
          stain_confirmed_val <- TRUE
        }
        
        p <- new_plate(slot_id = slot_id, growth_data = growth_df, layout = layout_data)
        
        if (!is.null(stain_df) && nrow(stain_df) > 0) {
          p@assays$staining <- stain_df
          p@staining_hr <- stain_hr_final
          p@staining_confirmed <- stain_confirmed_val
        }
        
        if (!is.null(treat_val) && treat_val != "") p@treatment <- treat_val
        if (!is.null(rep_val)) p@replicate <- as.integer(rep_val)
        if (!is.null(media_val) && media_val != "") p@media <- media_val
        if (!is.null(label_val) && label_val != "") p@label <- label_val
        
        rv$registry <- add_plate(rv$registry, p)
      }
      
      save_registry(rv$registry)
      removeModal()
      shinyjs::info("Plates successfully created and saved.")
    })
    
  })
}
