#' @import shiny
#' @importFrom stats sd pnorm p.adjust
#' @noRd
mod_scoring_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    fluidRow(
      column(
        3,
        box(
          title = "Scoring Configuration", width = 12, status = "primary", solidHeader = TRUE,
          uiOutput(ns("batch_selector")),
          selectInput(
            ns("metric"), "Phenotype Metric:",
            choices = c(
              "Growth AUC" = "Growth_AUC",
              "Growth AUC Rate" = "Growth_AUC_Rate",
              "Growth Endpoint OD" = "Growth_Endpoint_OD",
              "Biofilm Stain OD" = "Biofilm_Stain_OD",
              "Relative Biofilm Score" = "Biofilm_Score_RB"
            ),
            selected = "Growth_AUC"
          ),
          textInput(ns("control_drug"), "Control Drug Name:", value = "NoDrug"),
          hr(style = "margin: 10px 0;"),
          radioButtons(
            ns("score_mode"), "Scoring Mode:",
            choices = c(
              "One-vs-One" = "one_vs_one",
              "Merge (Metric-Level)" = "merge"
            ),
            selected = "one_vs_one"
          ),
          uiOutput(ns("mode_options_ui")),
          hr(style = "margin: 10px 0;"),
          numericInput(
            ns("score_alpha"), "FDR Highlight Cutoff (adj. p):",
            value = 0.05, min = 0.001, max = 1, step = 0.01
          ),
          actionButton(
            ns("btn_score"), "\u25b6 Run Scoring",
            class = "btn btn-primary btn-block",
            style = "margin-top: 10px; font-weight: bold;"
          ),
          downloadButton(
            ns("btn_download_csv"), "Download All Results (CSV)",
            class = "btn btn-success btn-block",
            style = "margin-top: 8px;"
          )
        )
      ),
      column(
        9,
        uiOutput(ns("run_status_ui")),
        fluidRow(
          valueBoxOutput(ns("box_total_scored"), width = 4),
          valueBoxOutput(ns("box_sig_pos"), width = 4),
          valueBoxOutput(ns("box_sig_neg"), width = 4)
        ),
        tabBox(
          width = 12, title = "Results (Unfiltered)", id = ns("results_tabbox"),
          tabPanel(
            "Ranked Hits Plot",
            plotly::plotlyOutput(ns("scoring_plot"), height = "520px") %>% shinycssloaders::withSpinner(type = 6)
          ),
          tabPanel(
            "Scoring Table",
            DT::dataTableOutput(ns("dt_scores")) %>% shinycssloaders::withSpinner(type = 6)
          ),
          tabPanel(
            "4-Strategy Comparison",
            uiOutput(ns("strategy_comparison_info")),
            plotly::plotlyOutput(ns("strategy_comparison_plot"), height = "360px") %>% shinycssloaders::withSpinner(type = 6),
            hr(),
            h4(icon("table"), " Significant Hits Overlap Across Merge Strategies"),
            DT::dataTableOutput(ns("strategy_comparison_table")) %>% shinycssloaders::withSpinner(type = 6)
          )
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

    local_rv <- reactiveValues(
      scoring_run      = NULL,
      scoring_run_all4 = NULL,
      scoring_run_meta = NULL
    )

    output$batch_selector <- renderUI({
      batches <- rv$registry@batches
      if (length(batches) == 0L) {
        selectInput(ns("selected_batch"), "Select Batch:", choices = "No batches available")
      } else {
        selectInput(ns("selected_batch"), "Select Batch:", choices = names(batches))
      }
    })

    # Detect whether selected batch has intra-plate _NoDrug control wells
    batch_has_nodrug_wells <- reactive({
      b_name <- input$selected_batch
      if (is.null(b_name) || b_name == "No batches available") return(FALSE)
      batch <- rv$registry@batches[[b_name]]
      if (is.null(batch) || length(batch@plate_uuids) == 0L) return(FALSE)
      for (u in batch@plate_uuids) {
        p <- rv$registry@plates[[u]]
        if (!is.null(p) && !is.null(p@layout) && "Is_NoDrug_Control" %in% names(p@layout)) {
          if (any(p@layout$Is_NoDrug_Control == TRUE | p@layout$Is_NoDrug_Control == "Y", na.rm = TRUE)) {
            return(TRUE)
          }
        }
      }
      FALSE
    })

    # Reset displayed run when batch changes
    observeEvent(input$selected_batch, {
      local_rv$scoring_run <- NULL
      local_rv$scoring_run_all4 <- NULL
      local_rv$scoring_run_meta <- NULL
    })

    # Dynamic mode-specific options (One-vs-One vs Merge)
    output$mode_options_ui <- renderUI({
      req(input$score_mode)

      if (input$score_mode == "merge") {
        return(tagList(
          radioButtons(
            ns("combine_strategy"), "Merge Strategy:",
            choices = c(
              "Average (Mean)" = "Average",
              "Min (Worst-Case)" = "Min",
              "Max (Best-Case)" = "Max",
              "Variance-Aware (Huber + EB)" = "Robust_Variance_Aware"
            ),
            selected = "Average"
          ),
          checkboxInput(
            ns("exclude_low_growth"),
            "Exclude low-growth NoDrug mutants (QC)",
            value = FALSE
          ),
          tags$small(
            style = "display: block; color: #7f8c8d; margin-top: -6px;",
            "Combines plate-level phenotypes across replicates first, then scores all mutants without filtration."
          )
        ))
      }

      # One-vs-One mode
      ovo_choices <- if (isTRUE(batch_has_nodrug_wells())) {
        c("NoDrug Plates (Inter-plate)" = "NoDrug_Plates", "NoDrug Wells (Intra-plate)" = "NoDrug_Wells")
      } else {
        c("NoDrug Plates (Inter-plate)" = "NoDrug_Plates")
      }

      tagList(
        radioButtons(ns("strategy"), "Control Basis:", choices = ovo_choices, selected = "NoDrug_Plates"),
        uiOutput(ns("pair_selector_ui"))
      )
    })

    # Available replicate pairs for One-vs-One comparison
    output$pair_selector_ui <- renderUI({
      req(input$selected_batch, input$strategy)
      b_name <- input$selected_batch
      if (b_name == "No batches available") {
        return(selectInput(ns("replicate_pair"), "Select Pair to Compare:", choices = "All Pairs"))
      }

      batch <- rv$registry@batches[[b_name]]
      ctrl_drug <- if (is.null(input$control_drug) || !nzchar(input$control_drug)) "NoDrug" else input$control_drug
      pairs <- character(0)

      if (!is.null(batch) && length(batch@plate_uuids) > 0L) {
        plates <- rv$registry@plates[batch@plate_uuids]
        trts <- vapply(plates, function(p) if (is.null(p)) NA_character_ else as.character(p@treatment), character(1))
        reps <- vapply(plates, function(p) if (is.null(p)) NA_character_ else as.character(p@replicate), character(1))
        ctrl_reps <- sort(unique(reps[!is.na(trts) & trts == ctrl_drug]))
        drug_reps <- sort(unique(reps[!is.na(trts) & trts != ctrl_drug]))

        if (input$strategy == "NoDrug_Plates" && length(ctrl_reps) > 0L && length(drug_reps) > 0L) {
          pairs <- as.vector(outer(ctrl_reps, drug_reps, function(j, i) paste0("Control_Rep_", j, "_vs_Drug_Rep_", i)))
        } else if (input$strategy == "NoDrug_Wells" && length(drug_reps) > 0L) {
          pairs <- paste0("IntraPlate_Drug_Rep_", drug_reps)
        }
      }

      choices <- c("All Pairs", pairs)
      selectInput(
        ns("replicate_pair"), "Select Pair to Compare:",
        choices = choices,
        selected = if (length(pairs) > 0L) pairs[1] else "All Pairs"
      )
    })

    # Run Scoring button handler -> generates ALL unfiltered results
    observeEvent(input$btn_score, {
      b_name <- input$selected_batch
      if (is.null(b_name) || b_name == "No batches available") {
        shinyjs::alert("Please create and select a batch first.")
        return()
      }

      batch <- rv$registry@batches[[b_name]]
      if (is.null(batch)) return()

      shinyjs::disable("btn_score")
      id_msg <- showNotification("Running scoring...", duration = NULL, type = "message")

      tryCatch({
        mode <- if (is.null(input$score_mode)) "one_vs_one" else input$score_mode
        met  <- input$metric
        ctrl <- if (is.null(input$control_drug) || !nzchar(input$control_drug)) "NoDrug" else input$control_drug

        if (mode == "merge") {
          comb_strat <- if (is.null(input$combine_strategy)) "Average" else input$combine_strategy
          excl_lg    <- isTRUE(input$exclude_low_growth)

          b_scored <- score_batch(
            batch, rv$registry,
            strategies = "Metric_Combined",
            metrics = c(met),
            control_drug = ctrl,
            combine_strategies = c("Min", "Max", "Average", "Robust_Variance_Aware"),
            exclude_low_growth = excl_lg
          )

          rv$registry <- add_batch(rv$registry, b_scored)
          tryCatch(save_registry(rv$registry), error = function(e) NULL)

          df_all4 <- get_scores(b_scored, strategy = "Metric_Combined", metric = met)
          if (is.null(df_all4) || nrow(df_all4) == 0L) {
            stop("No metric-level combined scores produced for this batch and metric.")
          }

          local_rv$scoring_run_all4 <- df_all4
          local_rv$scoring_run <- df_all4 %>% dplyr::filter(Combine_Strategy == comb_strat)

          strat_lbl <- switch(comb_strat,
            "Average" = "Average (Mean)",
            "Min" = "Min (Worst-Case)",
            "Max" = "Max (Best-Case)",
            "Robust_Variance_Aware" = "Variance-Aware (Huber + EB)"
          )
          local_rv$scoring_run_meta <- list(
            mode = "merge",
            label = paste0("Merge \u2014 ", strat_lbl),
            strategy_id = comb_strat,
            metric = met,
            n_rows = nrow(local_rv$scoring_run),
            excluded_genes = attr(df_all4, "excluded_low_growth_genes"),
            timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
          )
        } else {
          ovo_strat <- if (is.null(input$strategy)) "NoDrug_Plates" else input$strategy
          pair_sel  <- if (is.null(input$replicate_pair)) "All Pairs" else input$replicate_pair

          b_scored <- score_batch(
            batch, rv$registry,
            strategies = c(ovo_strat),
            metrics = c(met),
            control_drug = ctrl
          )

          rv$registry <- add_batch(rv$registry, b_scored)
          tryCatch(save_registry(rv$registry), error = function(e) NULL)

          df_raw <- get_scores(b_scored, strategy = ovo_strat, metric = met)
          if (is.null(df_raw) || nrow(df_raw) == 0L) {
            stop("No One-vs-One scores produced for this batch and metric.")
          }

          if (pair_sel != "All Pairs") {
            df_out <- df_raw %>%
              dplyr::filter(Replicate_Pair == pair_sel) %>%
              dplyr::group_by(Media, Treatment) %>%
              dplyr::mutate(
                mean_log_cg = mean(log_cg_score, na.rm = TRUE),
                sd_log_cg = stats::sd(log_cg_score, na.rm = TRUE),
                z_score = dplyr::if_else(is.na(sd_log_cg) | sd_log_cg == 0, NA_real_, (log_cg_score - mean_log_cg) / sd_log_cg),
                p_value = 2 * stats::pnorm(-abs(z_score)),
                adj_p = stats::p.adjust(p_value, method = "BH"),
                Is_Significant = dplyr::if_else(!is.na(adj_p) & adj_p < 0.05, "Y", "N")
              ) %>%
              dplyr::ungroup()
          } else {
            df_out <- df_raw %>%
              dplyr::group_by(Media, Treatment, Replicate_Pair) %>%
              dplyr::mutate(
                mean_log_cg = mean(log_cg_score, na.rm = TRUE),
                sd_log_cg = stats::sd(log_cg_score, na.rm = TRUE),
                z_score = dplyr::if_else(is.na(sd_log_cg) | sd_log_cg == 0, NA_real_, (log_cg_score - mean_log_cg) / sd_log_cg),
                p_value = 2 * stats::pnorm(-abs(z_score)),
                adj_p = stats::p.adjust(p_value, method = "BH"),
                Is_Significant = dplyr::if_else(!is.na(adj_p) & adj_p < 0.05, "Y", "N")
              ) %>%
              dplyr::ungroup()
          }

          local_rv$scoring_run_all4 <- NULL
          local_rv$scoring_run <- df_out
          local_rv$scoring_run_meta <- list(
            mode = "one_vs_one",
            label = paste0("One-vs-One (", ovo_strat, " | ", pair_sel, ")"),
            strategy_id = paste0(ovo_strat, "_", gsub("[^A-Za-z0-9_]", "_", pair_sel)),
            metric = met,
            n_rows = nrow(df_out),
            excluded_genes = character(0),
            timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
          )
        }

        showNotification(
          sprintf("Scoring completed: %d unfiltered results generated.", local_rv$scoring_run_meta$n_rows),
          type = "message"
        )
      }, error = function(e) {
        showNotification(paste("Scoring failed:", e$message), type = "error")
      }, finally = {
        removeNotification(id_msg)
        shinyjs::enable("btn_score")
      })
    })

    # Status banner
    output$run_status_ui <- renderUI({
      m <- local_rv$scoring_run_meta
      if (is.null(m)) {
        return(wellPanel(
          style = "background: #f8f9fa; border-left: 4px solid #3c8dbc; padding: 12px 15px; margin-bottom: 15px;",
          icon("info-circle", style = "color: #3c8dbc;"),
          tags$b(" Ready to score. "),
          "Select ", tags$b("One-vs-One"), " or ", tags$b("Merge (Metric-Level)"),
          " on the left and click ", tags$b("Run Scoring"), " to generate, plot, and download all unfiltered results."
        ))
      }

      excl_txt <- if (!is.null(m$excluded_genes) && length(m$excluded_genes) > 0L) {
        sprintf(" | Low-growth QC excluded: %d (%s)", length(m$excluded_genes), paste(m$excluded_genes, collapse = ", "))
      } else {
        " | No mutant filtration"
      }

      tags$div(
        style = "background: #eaf6ee; border-left: 4px solid #00a65a; padding: 10px 15px; margin-bottom: 15px; border-radius: 4px; font-size: 13px;",
        icon("check-circle", style = "color: #00a65a;"),
        tags$b(" Active Run: "), m$label,
        " | ", tags$b("Metric: "), m$metric,
        " | ", tags$b("Rows: "), m$n_rows,
        excl_txt,
        tags$span(style = "float: right; color: #7f8c8d;", m$timestamp)
      )
    })

    # Summary value boxes
    output$box_total_scored <- renderValueBox({
      val <- if (!is.null(local_rv$scoring_run)) nrow(local_rv$scoring_run) else 0L
      valueBox(val, "Total Scored (Unfiltered)", icon = icon("list"), color = "aqua")
    })

    output$box_sig_pos <- renderValueBox({
      alpha <- if (is.null(input$score_alpha)) 0.05 else input$score_alpha
      val <- if (!is.null(local_rv$scoring_run)) {
        sum(local_rv$scoring_run$log_cg_score > 0 & !is.na(local_rv$scoring_run$adj_p) & local_rv$scoring_run$adj_p < alpha, na.rm = TRUE)
      } else {
        0L
      }
      valueBox(val, "Significant Positive Hits", icon = icon("plus-circle"), color = "green")
    })

    output$box_sig_neg <- renderValueBox({
      alpha <- if (is.null(input$score_alpha)) 0.05 else input$score_alpha
      val <- if (!is.null(local_rv$scoring_run)) {
        sum(local_rv$scoring_run$log_cg_score < 0 & !is.na(local_rv$scoring_run$adj_p) & local_rv$scoring_run$adj_p < alpha, na.rm = TRUE)
      } else {
        0L
      }
      valueBox(val, "Significant Negative Hits", icon = icon("minus-circle"), color = "red")
    })

    # Direct download of ALL unfiltered results
    output$btn_download_csv <- downloadHandler(
      filename = function() {
        m <- local_rv$scoring_run_meta
        if (is.null(m)) {
          paste0("CGExplorer_Scores_", input$selected_batch, "_", input$metric, ".csv")
        } else {
          paste0("CGExplorer_Scores_", input$selected_batch, "_", m$strategy_id, "_", m$metric, ".csv")
        }
      },
      content = function(file) {
        df <- local_rv$scoring_run
        if (!is.null(df) && nrow(df) > 0L) {
          readr::write_csv(df, file)
        }
      }
    )

    # Ranked Hits Plot (all unfiltered results)
    output$scoring_plot <- plotly::renderPlotly({
      df <- local_rv$scoring_run
      if (is.null(df) || nrow(df) == 0L) {
        return(plotly::plotly_empty(type = "scatter", mode = "markers") %>%
          plotly::add_annotations(text = "Select options on the left and click 'Run Scoring' to generate plots.", showarrow = FALSE) %>%
          plotly::config(displayModeBar = FALSE))
      }

      alpha <- if (is.null(input$score_alpha)) 0.05 else input$score_alpha
      df_plot <- df %>%
        dplyr::arrange(log_cg_score) %>%
        dplyr::mutate(
          Rank = dplyr::row_number(),
          Significance = dplyr::case_when(
            !is.na(adj_p) & adj_p < alpha & log_cg_score > 0 ~ "Significant Positive",
            !is.na(adj_p) & adj_p < alpha & log_cg_score < 0 ~ "Significant Negative",
            TRUE ~ "Not Significant"
          )
        )

      has_robust <- "se_log2_cg" %in% names(df_plot) && any(!is.na(df_plot$se_log2_cg))
      df_plot$hover_txt <- paste0(
        "Gene: ", df_plot$Gene,
        "<br>Treatment: ", df_plot$Treatment,
        if ("Media" %in% names(df_plot)) paste0("<br>Media: ", df_plot$Media) else "",
        if ("Replicate_Pair" %in% names(df_plot)) paste0("<br>Pair/Strategy: ", df_plot$Replicate_Pair) else "",
        "<br>Log2 CG: ", round(df_plot$log_cg_score, 4),
        "<br>CG (Linear): ", round(df_plot$cg_score, 4),
        if ("z_score" %in% names(df_plot)) paste0("<br>Z-score: ", round(df_plot$z_score, 3)) else "",
        "<br>p-value: ", format.pval(df_plot$p_value, digits = 4),
        "<br>FDR adj. p: ", format.pval(df_plot$adj_p, digits = 4),
        if (has_robust) paste0(
          "<br>SE (Log2 CG): ", round(df_plot$se_log2_cg, 4),
          "<br>Shrinkage: ", round(df_plot$shrinkage_factor, 3)
        ) else ""
      )

      color_map <- c(
        "Significant Positive" = "#2E7D32",
        "Significant Negative" = "#C62828",
        "Not Significant"      = "#B0BEC5"
      )

      m <- local_rv$scoring_run_meta
      plot_title <- if (!is.null(m)) {
        paste0("Ranked Interaction Scores \u2014 ", m$label, " (", m$metric, ")")
      } else {
        "Ranked Chemo-Genomic Interaction Scores (Log2)"
      }

      plotly::plot_ly(
        data = df_plot,
        x = ~Rank,
        y = ~log_cg_score,
        type = "scatter",
        mode = "markers",
        color = ~Significance,
        colors = color_map,
        marker = list(size = 7, opacity = 0.75),
        text = ~hover_txt,
        hoverinfo = "text"
      ) %>%
        plotly::layout(
          title = list(text = plot_title, font = list(size = 15)),
          xaxis = list(title = "Rank (All Unfiltered Results)", showgrid = FALSE),
          yaxis = list(title = "Log2 Interaction Score", showgrid = TRUE),
          legend = list(orientation = "h", x = 0, y = -0.15),
          margin = list(l = 60, r = 20, t = 45, b = 60)
        ) %>%
        plotly::config(displayModeBar = TRUE)
    })

    # Unfiltered Scoring Table
    output$dt_scores <- DT::renderDataTable({
      df <- local_rv$scoring_run
      if (is.null(df) || nrow(df) == 0L) return(NULL)

      has_robust <- "se_log2_cg" %in% names(df) && any(!is.na(df$se_log2_cg))
      has_n_reps <- "N_Replicates" %in% names(df)

      if (has_robust) {
        display_df <- df %>%
          dplyr::select(
            Media, Treatment, Gene,
            `N Reps` = N_Replicates,
            `Linear CG` = cg_score,
            `Log2 CG` = log_cg_score,
            `SE (Log2 CG)` = se_log2_cg,
            Shrinkage = shrinkage_factor,
            `Z-score` = z_score,
            `p-value` = p_value,
            `FDR adj. p` = adj_p
          ) %>%
          dplyr::arrange(`FDR adj. p`)
        round_cols <- c("Linear CG", "Log2 CG", "SE (Log2 CG)", "Shrinkage", "Z-score")
      } else if (has_n_reps) {
        display_df <- df %>%
          dplyr::select(
            Media, Treatment, Gene,
            `N Reps` = N_Replicates,
            `Linear CG` = cg_score,
            `Log2 CG` = log_cg_score,
            `Z-score` = z_score,
            `p-value` = p_value,
            `FDR adj. p` = adj_p
          ) %>%
          dplyr::arrange(`FDR adj. p`)
        round_cols <- c("Linear CG", "Log2 CG", "Z-score")
      } else {
        display_df <- df %>%
          dplyr::select(
            Media, Treatment,
            Pair = Replicate_Pair,
            Gene,
            `Linear CG` = cg_score,
            `Log2 CG` = log_cg_score,
            `Z-score` = z_score,
            `p-value` = p_value,
            `FDR adj. p` = adj_p
          ) %>%
          dplyr::arrange(`FDR adj. p`)
        round_cols <- c("Linear CG", "Log2 CG", "Z-score")
      }

      DT::datatable(
        display_df,
        extensions = "Buttons",
        filter = "top",
        rownames = FALSE,
        options = list(
          pageLength = 15,
          scrollX = TRUE,
          dom = "Bfrtip",
          buttons = list(list(extend = "csv", text = "Download Table CSV"))
        )
      ) %>%
        DT::formatRound(columns = round_cols, digits = 4) %>%
        DT::formatSignif(columns = c("p-value", "FDR adj. p"), digits = 4)
    }, server = FALSE)

    # ── 4-Strategy Comparison Tab (Merge Mode) ────────────────────────────────
    output$strategy_comparison_info <- renderUI({
      df_mc <- local_rv$scoring_run_all4
      if (is.null(df_mc) || nrow(df_mc) == 0L) {
        return(wellPanel(
          p(icon("info-circle"), " Select ", tags$b("Merge (Metric-Level)"), " and click ", tags$b("Run Scoring"), " to compare Min, Max, Average, and Variance-Aware.")
        ))
      }
      alpha <- if (is.null(input$score_alpha)) 0.05 else input$score_alpha
      excl <- attr(df_mc, "excluded_low_growth_genes")
      excl_msg <- if (!is.null(excl) && length(excl) > 0L) {
        sprintf("Low-growth QC excluded %d mutant(s): %s.", length(excl), paste(excl, collapse = ", "))
      } else {
        "No low-growth mutants excluded (unfiltered)."
      }
      tags$div(
        style = "margin-bottom: 10px; color: #555; font-size: 13px;",
        icon("layer-group"),
        tags$b(" Comparing Min, Max, Average, and Variance-Aware "),
        sprintf("at FDR adj. p < %g. ", alpha),
        tags$span(style = "color: #7f8c8d;", excl_msg)
      )
    })

    output$strategy_comparison_plot <- plotly::renderPlotly({
      df_mc <- local_rv$scoring_run_all4
      if (is.null(df_mc) || nrow(df_mc) == 0L) {
        return(plotly::plotly_empty(type = "scatter", mode = "markers") %>%
          plotly::add_annotations(text = "Run Scoring in Merge mode to compare all 4 strategies.", showarrow = FALSE) %>%
          plotly::config(displayModeBar = FALSE))
      }

      alpha <- if (is.null(input$score_alpha)) 0.05 else input$score_alpha
      strat_levels <- c("Min", "Max", "Average", "Robust_Variance_Aware")
      strat_labels <- c("Min", "Max", "Average", "Variance-Aware")

      counts_df <- df_mc %>%
        dplyr::mutate(Combine_Strategy = factor(Combine_Strategy, levels = strat_levels, labels = strat_labels)) %>%
        dplyr::group_by(Combine_Strategy) %>%
        dplyr::summarise(
          Sensitive = sum(!is.na(adj_p) & adj_p < alpha & log_cg_score < 0, na.rm = TRUE),
          Resistant = sum(!is.na(adj_p) & adj_p < alpha & log_cg_score > 0, na.rm = TRUE),
          .groups = "drop"
        )

      plotly::plot_ly(counts_df, x = ~Combine_Strategy) %>%
        plotly::add_bars(
          y = ~Sensitive, name = sprintf("Sensitive / Negative (adj. p < %g)", alpha),
          marker = list(color = "#2166ac"),
          hovertemplate = "%{x}<br>Sensitive Hits: %{y}<extra></extra>"
        ) %>%
        plotly::add_bars(
          y = ~Resistant, name = sprintf("Resistant / Positive (adj. p < %g)", alpha),
          marker = list(color = "#b2182b"),
          hovertemplate = "%{x}<br>Resistant Hits: %{y}<extra></extra>"
        ) %>%
        plotly::layout(
          barmode = "stack",
          title = list(text = sprintf("Significant Hit Counts Across 4 Merge Strategies (adj. p < %g)", alpha), font = list(size = 15)),
          xaxis = list(title = "Merge Strategy"),
          yaxis = list(title = "Number of Significant Hits"),
          legend = list(orientation = "h", x = 0, y = -0.18),
          margin = list(l = 60, r = 20, t = 45, b = 70)
        ) %>%
        plotly::config(displayModeBar = TRUE)
    })

    output$strategy_comparison_table <- DT::renderDataTable({
      df_mc <- local_rv$scoring_run_all4
      if (is.null(df_mc) || nrow(df_mc) == 0L) return(NULL)

      alpha <- if (is.null(input$score_alpha)) 0.05 else input$score_alpha

      sig_genes <- df_mc %>%
        dplyr::filter(!is.na(adj_p) & adj_p < alpha) %>%
        dplyr::select(Media, Treatment, Gene) %>%
        dplyr::distinct()

      if (nrow(sig_genes) == 0L) {
        return(DT::datatable(
          data.frame(Message = sprintf("No significant hits found across any of the 4 strategies at FDR adj. p < %g.", alpha)),
          options = list(dom = "t"), rownames = FALSE
        ))
      }

      overlap_df <- df_mc %>%
        dplyr::inner_join(sig_genes, by = c("Media", "Treatment", "Gene")) %>%
        dplyr::group_by(Media, Treatment, Gene) %>%
        dplyr::summarise(
          `N Strategies Hit` = sum(!is.na(adj_p) & adj_p < alpha, na.rm = TRUE),
          `Strategies Hit` = paste(Combine_Strategy[!is.na(adj_p) & adj_p < alpha], collapse = ", "),
          `Min log2(CG)` = dplyr::first(log_cg_score[Combine_Strategy == "Min"]),
          `Max log2(CG)` = dplyr::first(log_cg_score[Combine_Strategy == "Max"]),
          `Average log2(CG)` = dplyr::first(log_cg_score[Combine_Strategy == "Average"]),
          `Variance-Aware log2(CG)` = dplyr::first(log_cg_score[Combine_Strategy == "Robust_Variance_Aware"]),
          `Best adj. p` = min(adj_p, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        dplyr::arrange(dplyr::desc(`N Strategies Hit`), `Best adj. p`)

      DT::datatable(
        overlap_df,
        extensions = "Buttons",
        filter = "top",
        rownames = FALSE,
        options = list(
          pageLength = 15,
          scrollX = TRUE,
          dom = "Bfrtip",
          buttons = list(list(extend = "csv", text = "Download Overlap CSV"))
        )
      ) %>%
        DT::formatRound(columns = c("Min log2(CG)", "Max log2(CG)", "Average log2(CG)", "Variance-Aware log2(CG)"), digits = 4) %>%
        DT::formatSignif(columns = "Best adj. p", digits = 4)
    }, server = FALSE)
  })
}
