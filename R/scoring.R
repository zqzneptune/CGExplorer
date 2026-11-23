#' @import dplyr tidyr
#' @importFrom stats median mad sd pnorm p.adjust
NULL

#' Huber IRLS Robust Phenotypic Location and Weighted Variance
#'
#' Computes an outlier-resistant location estimate and weighted variance across
#' plate replicates at the phenotypic metric level using Huber iteratively
#' reweighted least squares (IRLS).
#'
#' @param p_vals numeric vector of plate-level phenotypic measurements
#' @param c_val numeric tuning constant for Huber weighting (default `1.345`)
#' @param n_iter integer number of IRLS iterations (default `6`)
#' @return A list with elements `mu_p` (robust mean), `var_p` (weighted variance),
#'   and `k_eff` (effective replicate count / sum of weights).
#' @export
huber_pheno_stats <- function(p_vals, c_val = 1.345, n_iter = 6L) {
  x <- p_vals[!is.na(p_vals)]
  n <- length(x)
  if (n == 0L) {
    return(list(mu_p = NA_real_, var_p = NA_real_, k_eff = 0))
  }
  if (n == 1L) {
    return(list(mu_p = x[1], var_p = 0, k_eff = 1))
  }
  mu <- stats::median(x)
  s <- max(1.4826 * stats::mad(x, constant = 1), 0.02 * max(abs(mu), 0.05))
  w <- rep(1, n)
  for (it in seq_len(n_iter)) {
    r <- abs(x - mu) / s
    w <- ifelse(r <= c_val, 1, c_val / r)
    mu <- sum(w * x) / sum(w)
  }
  k_eff <- sum(w)
  var_p <- sum(w * (x - mu)^2) / max(k_eff - 1, 1)
  list(mu_p = mu, var_p = var_p, k_eff = k_eff)
}

#' Filter Non-Growing Mutants in Untreated Control Plates
#'
#' Identifies and removes mutants that fail to grow under untreated control
#' conditions (`control_drug`), preventing division-by-near-zero artifacts when
#' computing relative fitness ratios.
#'
#' @param clean_df data.frame of phenotypic metrics (with `Treatment` or `Drug`)
#' @param control_drug character, name of untreated control treatment (default `"NoDrug"`)
#' @param min_endpoint_od numeric minimum mean untreated endpoint OD (default `0.10`)
#' @param min_delta_od numeric minimum mean untreated delta OD (default `0.05`)
#' @param min_rel_fitness numeric minimum relative untreated fitness vs WT (default `0.15`)
#' @param min_auc numeric minimum mean untreated Growth AUC (default `1.00`)
#' @return A list with `data` (filtered data.frame) and `excluded_genes` (character vector).
#' @export
filter_low_growth_mutants <- function(clean_df,
                                      control_drug = "NoDrug",
                                      min_endpoint_od = 0.10,
                                      min_delta_od = 0.05,
                                      min_rel_fitness = 0.15,
                                      min_auc = 1.00) {
  req_cols <- c("Growth_Endpoint_OD", "Growth_Max_OD", "Growth_Initial_OD", "Growth_AUC")
  if (!all(req_cols %in% names(clean_df))) {
    return(list(data = clean_df, excluded_genes = character(0)))
  }

  trt_col <- if ("Treatment" %in% names(clean_df)) "Treatment" else if ("Drug" %in% names(clean_df)) "Drug" else NULL
  if (is.null(trt_col) || !any(clean_df[[trt_col]] == control_drug, na.rm = TRUE)) {
    return(list(data = clean_df, excluded_genes = character(0)))
  }

  df_work <- clean_df
  if (!"Media" %in% names(df_work)) df_work$Media <- "Default"
  df_work$.__trt__ <- df_work[[trt_col]]

  wt_ref <- df_work %>%
    dplyr::filter(Type %in% c("WT_Control", "WT"), .__trt__ == control_drug) %>%
    dplyr::group_by(Media) %>%
    dplyr::summarise(
      WT_Endpoint_OD = mean(Growth_Endpoint_OD, na.rm = TRUE),
      WT_AUC         = mean(Growth_AUC, na.rm = TRUE),
      .groups = "drop"
    )

  if (nrow(wt_ref) == 0L) {
    return(list(data = clean_df, excluded_genes = character(0)))
  }

  low_growth_df <- df_work %>%
    dplyr::filter(Type == "Mutant", .__trt__ == control_drug) %>%
    dplyr::group_by(Media, Gene) %>%
    dplyr::summarise(
      Mean_NoDrug_Endpoint_OD = mean(Growth_Endpoint_OD, na.rm = TRUE),
      Mean_NoDrug_Delta_OD    = mean(Growth_Max_OD - Growth_Initial_OD, na.rm = TRUE),
      Mean_NoDrug_AUC         = mean(Growth_AUC, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::left_join(wt_ref, by = "Media") %>%
    dplyr::mutate(
      Rel_NoDrug_Endpoint_OD = dplyr::if_else(
        !is.na(WT_Endpoint_OD) & WT_Endpoint_OD > 0,
        Mean_NoDrug_Endpoint_OD / WT_Endpoint_OD,
        NA_real_
      ),
      Is_Low_Growth = (
        (!is.na(Mean_NoDrug_Endpoint_OD) & Mean_NoDrug_Endpoint_OD < min_endpoint_od) |
        (!is.na(Mean_NoDrug_Delta_OD)    & Mean_NoDrug_Delta_OD    < min_delta_od)    |
        (!is.na(Rel_NoDrug_Endpoint_OD)  & Rel_NoDrug_Endpoint_OD  < min_rel_fitness) |
        (!is.na(Mean_NoDrug_AUC)         & Mean_NoDrug_AUC         < min_auc)
      )
    ) %>%
    dplyr::filter(Is_Low_Growth)

  if (nrow(low_growth_df) == 0L) {
    return(list(data = clean_df, excluded_genes = character(0)))
  }

  excl_keys <- paste(low_growth_df$Media, low_growth_df$Gene, sep = "||")
  media_vec <- if ("Media" %in% names(clean_df)) clean_df$Media else rep("Default", nrow(clean_df))
  filtered_df <- clean_df[!(clean_df$Type == "Mutant" & paste(media_vec, clean_df$Gene, sep = "||") %in% excl_keys), , drop = FALSE]

  list(data = filtered_df, excluded_genes = sort(unique(low_growth_df$Gene)))
}

#' Metric-Level Combining Chemo-Genomic Scoring Engine
#'
#' Combines plate-level mutant and WT phenotypic metrics across replicates first
#' using one or more strategies (`"Min"`, `"Max"`, `"Average"`, and
#' `"Robust_Variance_Aware"`), then derives relative fitness, linear `cg_score`,
#' `log_cg_score`, Z-scores, and Benjamini-Hochberg FDR (`adj_p`).
#'
#' @param pheno_df data.frame containing well-level phenotypic metrics across a batch
#' @param metrics character vector of target metrics to score
#' @param control_drug character, name of untreated control treatment (default `"NoDrug"`)
#' @param combine_strategies character vector of combining strategies to compute.
#'   Subset of `c("Min", "Max", "Average", "Robust_Variance_Aware")`.
#' @param exclude_low_growth logical, whether to exclude non-growing mutants in
#'   `control_drug` prior to scoring (default `FALSE` for unfiltered scoring).
#' @return data.frame of scored interactions with `Combine_Strategy` column and
#'   attribute `excluded_low_growth_genes`.
#' @export
score_batch_metric_combined <- function(pheno_df,
                                        metrics = c("Growth_AUC", "Growth_AUC_Rate", "Growth_Endpoint_OD", "Biofilm_Stain_OD", "Biofilm_Score_RB"),
                                        control_drug = "NoDrug",
                                        combine_strategies = c("Min", "Max", "Average", "Robust_Variance_Aware"),
                                        exclude_low_growth = FALSE) {
  if (is.null(pheno_df) || nrow(pheno_df) == 0L) return(data.frame())

  valid_strategies <- c("Min", "Max", "Average", "Robust_Variance_Aware")
  combine_strategies <- intersect(combine_strategies, valid_strategies)
  if (length(combine_strategies) == 0L) return(data.frame())

  clean_df <- pheno_df
  if (!"Treatment" %in% names(clean_df) && "Drug" %in% names(clean_df)) {
    clean_df$Treatment <- clean_df$Drug
  }
  if (!"Treatment" %in% names(clean_df)) {
    stop("pheno_df must contain a 'Treatment' (or 'Drug') column.")
  }
  if ("Exclude_From_Scoring" %in% names(clean_df)) {
    clean_df <- clean_df %>% dplyr::filter(Exclude_From_Scoring == "N")
  }
  if ("Is_NoDrug_Control" %in% names(clean_df)) {
    clean_df <- clean_df %>%
      dplyr::filter(!(Treatment != control_drug & (Is_NoDrug_Control == "Y" | Is_NoDrug_Control == TRUE)))
  }

  clean_df <- clean_df %>%
    dplyr::filter(
      Type %in% c("Mutant", "WT_Control", "WT"),
      !Gene %in% c("Blank_NoDrug", "WT_NoDrug", "PseudoWT_NoDrug")
    )

  if (nrow(clean_df) == 0L) return(data.frame())

  batch_val <- if ("Batch" %in% names(clean_df)) unique(clean_df$Batch)[1] else "Batch_1"
  if (!"Batch" %in% names(clean_df)) clean_df$Batch <- batch_val
  if (!"Media" %in% names(clean_df)) clean_df$Media <- "Default"
  if (!"Plate_UUID" %in% names(clean_df)) {
    clean_df$Plate_UUID <- paste(clean_df$Batch, clean_df$Media, clean_df$Treatment, clean_df$Replicate, sep = "_")
  }

  excluded_genes <- character(0)
  if (isTRUE(exclude_low_growth)) {
    lg_res <- filter_low_growth_mutants(clean_df, control_drug = control_drug)
    clean_df <- lg_res$data
    excluded_genes <- lg_res$excluded_genes
  }

  all_results <- list()

  for (m in metrics) {
    if (!m %in% names(clean_df)) next

    metric_df <- clean_df %>%
      dplyr::filter(!is.na(!!dplyr::sym(m))) %>%
      dplyr::mutate(Phenotype_Val = !!dplyr::sym(m))

    if (nrow(metric_df) == 0L) next
    if (!any(metric_df$Treatment == control_drug)) next

    drugs <- sort(unique(metric_df$Treatment[metric_df$Treatment != control_drug]))
    if (length(drugs) == 0L) next

    # Step A: Aggregate intra-plate technical duplicate wells -> 1 value per plate
    plate_mut <- metric_df %>%
      dplyr::filter(Type == "Mutant") %>%
      dplyr::group_by(Batch, Plate_UUID, Media, Treatment, Replicate, Gene) %>%
      dplyr::summarise(Phenotype = mean(Phenotype_Val, na.rm = TRUE), .groups = "drop")

    plate_wt <- metric_df %>%
      dplyr::filter(Type %in% c("WT_Control", "WT")) %>%
      dplyr::group_by(Batch, Plate_UUID, Media, Treatment, Replicate) %>%
      dplyr::summarise(WT_Plate_P = mean(Phenotype_Val, na.rm = TRUE), .groups = "drop")

    if (nrow(plate_mut) == 0L || nrow(plate_wt) == 0L) next

    wt_summary <- plate_wt %>%
      dplyr::group_by(Media, Treatment) %>%
      dplyr::summarise(
        WT_Mean_P = max(mean(WT_Plate_P, na.rm = TRUE), 0.001),
        WT_SE_P   = ifelse(dplyr::n() > 1L, stats::sd(WT_Plate_P, na.rm = TRUE) / sqrt(dplyr::n()), 0),
        WT_Rob    = list(huber_pheno_stats(WT_Plate_P)),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        WT_Rob_P  = pmax(sapply(WT_Rob, `[[`, "mu_p"), 0.001),
        WT_Rob_SE = sqrt(pmax(sapply(WT_Rob, `[[`, "var_p"), 0) / pmax(sapply(WT_Rob, `[[`, "k_eff"), 1))
      ) %>%
      dplyr::select(-WT_Rob)

    media_vals <- sort(unique(plate_mut$Media))

    for (med in media_vals) {
      wt_med <- wt_summary %>% dplyr::filter(Media == med)
      wt_nd_row <- wt_med %>% dplyr::filter(Treatment == control_drug)
      if (nrow(wt_nd_row) == 0L) next

      wt_nd_mean <- wt_nd_row$WT_Mean_P[1]
      wt_nd_rob  <- wt_nd_row$WT_Rob_P[1]
      wt_nd_rse  <- ifelse(is.na(wt_nd_row$WT_Rob_SE[1]), 0, wt_nd_row$WT_Rob_SE[1])

      pm_med <- plate_mut %>% dplyr::filter(Media == med)
      med_drugs <- intersect(drugs, unique(pm_med$Treatment))
      if (length(med_drugs) == 0L) next

      # Step B: Min, Max, Average
      simple_strats <- intersect(combine_strategies, c("Min", "Max", "Average"))
      for (st in simple_strats) {
        comb_fn <- switch(st,
          "Min"     = function(x) min(x, na.rm = TRUE),
          "Max"     = function(x) max(x, na.rm = TRUE),
          "Average" = function(x) mean(x, na.rm = TRUE)
        )

        mut_comb <- pm_med %>%
          dplyr::group_by(Treatment, Gene) %>%
          dplyr::summarise(
            N_Replicates       = dplyr::n(),
            Phenotype_Combined = comb_fn(Phenotype),
            Phenotype_SD       = ifelse(dplyr::n() > 1L, stats::sd(Phenotype, na.rm = TRUE), NA_real_),
            .groups = "drop"
          )

        mut_nd <- mut_comb %>%
          dplyr::filter(Treatment == control_drug) %>%
          dplyr::select(Gene, Phenotype_NoDrug = Phenotype_Combined, SD_NoDrug = Phenotype_SD)

        if (nrow(mut_nd) == 0L) next

        for (d in med_drugs) {
          wt_d_row <- wt_med %>% dplyr::filter(Treatment == d)
          wt_d_mean <- if (nrow(wt_d_row) > 0L) wt_d_row$WT_Mean_P[1] else wt_nd_mean

          scored_d <- mut_comb %>%
            dplyr::filter(Treatment == d) %>%
            dplyr::inner_join(mut_nd, by = "Gene") %>%
            dplyr::mutate(
              Batch            = batch_val,
              Media            = med,
              Type             = "Mutant",
              Metric           = m,
              WT_NoDrug_P      = wt_nd_mean,
              WT_Drug_P        = wt_d_mean,
              w_drug           = Phenotype_Combined / WT_NoDrug_P,
              w_control        = Phenotype_NoDrug / WT_NoDrug_P,
              w_wt_drug        = WT_Drug_P / WT_NoDrug_P,
              expected_w_drug  = w_control * w_wt_drug,
              cg_score         = w_drug - expected_w_drug,
              log_cg_score     = log2(pmax(w_drug, 0.001)) - log2(pmax(expected_w_drug, 0.001)),
              mean_log_cg      = mean(log_cg_score, na.rm = TRUE),
              sd_log_cg        = stats::sd(log_cg_score, na.rm = TRUE),
              z_score          = dplyr::if_else(is.na(sd_log_cg) | sd_log_cg == 0, NA_real_, (log_cg_score - mean_log_cg) / sd_log_cg),
              p_value          = 2 * stats::pnorm(-abs(z_score)),
              adj_p            = stats::p.adjust(p_value, method = "BH"),
              Is_Significant   = dplyr::if_else(!is.na(adj_p) & adj_p < 0.05, "Y", "N"),
              se_log2_cg       = NA_real_,
              shrinkage_factor = NA_real_,
              log2_cg_raw      = log_cg_score,
              Strategy         = "Metric_Combined",
              Combine_Strategy = st,
              Replicate_Pair   = paste0("Combined_", st)
            )

          all_results[[length(all_results) + 1L]] <- scored_d
        }
      }

      # Step C: Robust Variance-Aware (Huber IRLS + Empirical Bayes + Delta Method + James-Stein)
      if ("Robust_Variance_Aware" %in% combine_strategies) {
        mut_rob <- pm_med %>%
          dplyr::group_by(Treatment, Gene) %>%
          dplyr::summarise(
            N_Replicates = dplyr::n(),
            res          = list(huber_pheno_stats(Phenotype)),
            Phenotype_SD = ifelse(dplyr::n() > 1L, stats::sd(Phenotype, na.rm = TRUE), NA_real_),
            .groups = "drop"
          ) %>%
          dplyr::mutate(
            Phenotype_Rob = sapply(res, `[[`, "mu_p"),
            var_p         = sapply(res, `[[`, "var_p"),
            k_eff         = sapply(res, `[[`, "k_eff")
          ) %>%
          dplyr::select(-res) %>%
          dplyr::group_by(Treatment) %>%
          dplyr::mutate(
            s0_sq_p   = dplyr::coalesce(stats::median(var_p, na.rm = TRUE), 0),
            d0_p      = 3,
            mod_var_p = (d0_p * s0_sq_p + pmax(k_eff - 1, 0) * dplyr::coalesce(var_p, 0)) / (d0_p + pmax(k_eff - 1, 0)),
            se_rob_p  = sqrt(pmax(mod_var_p, 0) / pmax(k_eff, 1))
          ) %>%
          dplyr::ungroup()

        mut_nd_rob <- mut_rob %>%
          dplyr::filter(Treatment == control_drug) %>%
          dplyr::select(Gene, Phenotype_NoDrug = Phenotype_Rob, SE_NoDrug_P = se_rob_p, SD_NoDrug = Phenotype_SD)

        if (nrow(mut_nd_rob) > 0L) {
          for (d in med_drugs) {
            wt_d_row <- wt_med %>% dplyr::filter(Treatment == d)
            wt_d_rob <- if (nrow(wt_d_row) > 0L) wt_d_row$WT_Rob_P[1] else wt_nd_rob
            wt_d_rse <- if (nrow(wt_d_row) > 0L && !is.na(wt_d_row$WT_Rob_SE[1])) wt_d_row$WT_Rob_SE[1] else 0

            scored_rob <- mut_rob %>%
              dplyr::filter(Treatment == d) %>%
              dplyr::inner_join(mut_nd_rob, by = "Gene") %>%
              dplyr::mutate(
                Batch              = batch_val,
                Media              = med,
                Type               = "Mutant",
                Metric             = m,
                Phenotype_Combined = Phenotype_Rob,
                WT_NoDrug_P        = wt_nd_rob,
                WT_Drug_P          = wt_d_rob,
                w_drug             = Phenotype_Combined / WT_NoDrug_P,
                w_control          = Phenotype_NoDrug / WT_NoDrug_P,
                w_wt_drug          = WT_Drug_P / WT_NoDrug_P,
                expected_w_drug    = w_control * w_wt_drug,
                cg_score_raw       = w_drug - expected_w_drug,
                log2_cg_raw        = log2(pmax(w_drug, 0.001)) - log2(pmax(expected_w_drug, 0.001)),
                se_log2_cg = (1 / log(2)) * sqrt(
                  (se_rob_p    / pmax(Phenotype_Combined, 0.001))^2 +
                  (SE_NoDrug_P / pmax(Phenotype_NoDrug, 0.001))^2 +
                  (wt_d_rse    / pmax(WT_Drug_P, 0.001))^2 +
                  (wt_nd_rse   / pmax(WT_NoDrug_P, 0.001))^2
                ),
                sigma_screen     = 1.4826 * stats::mad(log2_cg_raw, constant = 1, na.rm = TRUE),
                sigma_sq         = sigma_screen^2,
                se_sq            = se_log2_cg^2,
                tot_var          = sigma_sq + se_sq,
                shrinkage_factor = dplyr::if_else(is.na(tot_var) | tot_var <= 0, 1, sigma_sq / tot_var),
                log_cg_score     = log2_cg_raw * shrinkage_factor,
                cg_score         = cg_score_raw * shrinkage_factor,
                mean_log_cg      = stats::median(log2_cg_raw, na.rm = TRUE),
                sd_log_cg        = sqrt(tot_var),
                z_score          = dplyr::if_else(is.na(tot_var) | tot_var <= 0, NA_real_, (log2_cg_raw - stats::median(log2_cg_raw, na.rm = TRUE)) / sqrt(tot_var)),
                p_value          = 2 * stats::pnorm(-abs(z_score)),
                adj_p            = stats::p.adjust(p_value, method = "BH"),
                Is_Significant   = dplyr::if_else(!is.na(adj_p) & adj_p < 0.05, "Y", "N"),
                Strategy         = "Metric_Combined",
                Combine_Strategy = "Robust_Variance_Aware",
                Replicate_Pair   = "Combined_Robust_Variance_Aware"
              ) %>%
              dplyr::select(
                Batch, Media, Treatment, Gene, Type, Metric,
                N_Replicates, Phenotype_Combined, Phenotype_SD, Phenotype_NoDrug, SD_NoDrug,
                WT_NoDrug_P, WT_Drug_P, w_drug, w_control, w_wt_drug, expected_w_drug,
                cg_score, log_cg_score, mean_log_cg, sd_log_cg, z_score, p_value, adj_p,
                Is_Significant, se_log2_cg, shrinkage_factor, log2_cg_raw,
                Strategy, Combine_Strategy, Replicate_Pair
              )

            all_results[[length(all_results) + 1L]] <- scored_rob
          }
        }
      }
    }
  }

  if (length(all_results) == 0L) return(data.frame())

  res_df <- dplyr::bind_rows(all_results)
  attr(res_df, "excluded_low_growth_genes") <- excluded_genes
  return(res_df)
}

#' Standardized Scoring Engine
#'
#' @param pheno_df data.frame containing aggregated phenotypic metrics across a batch
#' @param strategy_type character, one of `"NoDrug_Plates"`, `"NoDrug_Wells"`, or `"Metric_Combined"`
#' @param metrics character vector of target metrics
#' @param control_drug character, name of control treatment
#' @param ... additional arguments passed to `score_batch_metric_combined()` when
#'   `strategy_type == "Metric_Combined"` (e.g., `combine_strategies`, `exclude_low_growth`).
#' @return data.frame of scores
#' @export
score_batch_dataset <- function(pheno_df, strategy_type = c("NoDrug_Plates", "NoDrug_Wells", "Metric_Combined"), metrics, control_drug = "NoDrug", ...) {
  strategy_type <- match.arg(strategy_type)
  if (strategy_type == "Metric_Combined") {
    return(score_batch_metric_combined(pheno_df, metrics = metrics, control_drug = control_drug, ...))
  }
  scored_records <- list()
  
  clean_df <- pheno_df %>% dplyr::filter(Exclude_From_Scoring == "N")
  
  for (metric in metrics) {
    if (!metric %in% names(clean_df)) next
    
    metric_df <- clean_df %>%
      dplyr::filter(!is.na(!!dplyr::sym(metric))) %>%
      dplyr::rename(Phenotype = !!dplyr::sym(metric))
    
    if (nrow(metric_df) == 0) next
    
    drugs <- unique(metric_df$Treatment[metric_df$Treatment != control_drug])
    
    # -- STRATEGY 1: USING NODRUG PLATES --------------
    if (strategy_type == "NoDrug_Plates") {
      valid_metric_df <- if ("Is_NoDrug_Control" %in% names(metric_df)) {
        metric_df %>%
          dplyr::filter(!(Treatment != control_drug & (Is_NoDrug_Control == "Y" | Is_NoDrug_Control == TRUE)))
      } else {
        metric_df
      }
      
      nodrug_plates_df <- valid_metric_df %>% dplyr::filter(Treatment == control_drug)
      if (nrow(nodrug_plates_df) == 0) next
      
      well_agg <- valid_metric_df %>%
        dplyr::group_by(Batch, Plate_UUID, Media, Treatment, Replicate, Gene, Type) %>%
        dplyr::summarise(Phenotype = mean(Phenotype, na.rm = TRUE), .groups = "drop")
      
      nodrug_reps <- sort(unique(nodrug_plates_df$Replicate))
      
      for (d in drugs) {
        drug_df <- well_agg %>% dplyr::filter(Treatment == d)
        drug_reps <- sort(unique(drug_df$Replicate))
        
        for (i_d in drug_reps) {
          d_rep_data <- drug_df %>% dplyr::filter(Replicate == i_d)
          
          for (j_c in nodrug_reps) {
            ctrl_rep_data <- well_agg %>% dplyr::filter(Treatment == control_drug & Replicate == j_c)
            if (nrow(ctrl_rep_data) == 0) next
            
            wt_ctrl_val <- mean(ctrl_rep_data$Phenotype[ctrl_rep_data$Type == "WT_Control"], na.rm = TRUE)
            wt_ctrl_val <- pmax(ifelse(is.na(wt_ctrl_val), 1, wt_ctrl_val), 0.001)
            
            mut_ctrl_avg <- ctrl_rep_data %>%
              dplyr::filter(Type == "Mutant") %>%
              dplyr::group_by(Gene) %>%
              dplyr::summarise(Control_P = mean(Phenotype, na.rm = TRUE), .groups = "drop")
            
            wt_drug_val <- median(d_rep_data$Phenotype[d_rep_data$Type == "WT_Control"], na.rm = TRUE)
            wt_drug_val <- ifelse(is.na(wt_drug_val), wt_ctrl_val, wt_drug_val)
            
            scored_pair <- d_rep_data %>%
              dplyr::filter(Type == "Mutant") %>%
              dplyr::left_join(mut_ctrl_avg, by = "Gene") %>%
              dplyr::mutate(
                Mean_WT_Control_P = wt_ctrl_val,
                Mean_Control_P = dplyr::coalesce(Control_P, Mean_WT_Control_P),
                WT_Drug_P = wt_drug_val,
                
                w_drug = Phenotype / Mean_WT_Control_P,
                w_control = Mean_Control_P / Mean_WT_Control_P,
                w_wt_drug = WT_Drug_P / Mean_WT_Control_P,
                
                cg_score = w_drug - (w_control * w_wt_drug),
                log_cg_score = log2(pmax(w_drug, 0.001)) - (log2(pmax(w_control, 0.001)) + log2(pmax(w_wt_drug, 0.001))),
                
                Strategy = "NoDrug_Plates",
                Metric = metric,
                Replicate_Pair = paste0("Control_Rep_", j_c, "_vs_Drug_Rep_", i_d)
              )
            
            scored_records[[length(scored_records) + 1]] <- scored_pair
          }
        }
      }
    }
    
    # -- STRATEGY 2: USING NODRUG WELLS ---------------
    if (strategy_type == "NoDrug_Wells") {
      has_nodrug_wells <- "Is_NoDrug_Control" %in% names(metric_df) && 
        any(metric_df$Is_NoDrug_Control == "Y" | metric_df$Is_NoDrug_Control == TRUE, na.rm = TRUE)
      
      if (!has_nodrug_wells) {
        stop("Intra-plate calculation (NoDrug_Wells) is unavailable because no '_NoDrug' control wells are present on the plate layout.")
      }
      
      drug_only_df <- metric_df %>% dplyr::filter(Treatment != control_drug)
      
      well_agg <- drug_only_df %>%
        dplyr::group_by(Batch, Plate_UUID, Media, Treatment, Replicate, Gene, Type, Is_NoDrug_Control) %>%
        dplyr::summarise(Phenotype = mean(Phenotype, na.rm = TRUE), .groups = "drop")
      
      for (d in drugs) {
        d_df <- well_agg %>% dplyr::filter(Treatment == d)
        plates <- unique(d_df$Plate_UUID)
        
        for (p_id in plates) {
          plate_data <- d_df %>% dplyr::filter(Plate_UUID == p_id)
          rep_val <- plate_data$Replicate[1]
          
          intra_wt_nodrug <- mean(plate_data$Phenotype[plate_data$Is_NoDrug_Control == "Y"], na.rm = TRUE)
          if (is.na(intra_wt_nodrug) || intra_wt_nodrug <= 0) {
            intra_wt_nodrug <- mean(plate_data$Phenotype[plate_data$Type == "WT_Control"], na.rm = TRUE)
          }
          intra_wt_nodrug <- pmax(dplyr::coalesce(intra_wt_nodrug, 1), 0.001)
          
          wt_drug_val <- median(plate_data$Phenotype[plate_data$Type == "WT_Control"], na.rm = TRUE)
          wt_drug_val <- ifelse(is.na(wt_drug_val), intra_wt_nodrug, wt_drug_val)
          
          mut_ctrl_avg <- well_agg %>%
            dplyr::filter(Type == "Mutant" & Is_NoDrug_Control == "Y") %>%
            dplyr::group_by(Gene) %>%
            dplyr::summarise(Control_P = mean(Phenotype, na.rm = TRUE), .groups = "drop")
          
          scored_plate <- plate_data %>%
            dplyr::filter(Type == "Mutant" & Is_NoDrug_Control == "N") %>%
            dplyr::left_join(mut_ctrl_avg, by = "Gene") %>%
            dplyr::mutate(
              Mean_WT_Control_P = intra_wt_nodrug,
              Mean_Control_P = dplyr::coalesce(Control_P, Mean_WT_Control_P),
              WT_Drug_P = wt_drug_val,
              
              w_drug = Phenotype / Mean_WT_Control_P,
              w_control = Mean_Control_P / Mean_WT_Control_P,
              w_wt_drug = WT_Drug_P / Mean_WT_Control_P,
              
              cg_score = w_drug - (w_control * w_wt_drug),
              log_cg_score = log2(pmax(w_drug, 0.001)) - (log2(pmax(w_control, 0.001)) + log2(pmax(w_wt_drug, 0.001))),
              
              Strategy = "NoDrug_Wells",
              Metric = metric,
              Replicate_Pair = paste0("IntraPlate_Drug_Rep_", rep_val)
            )
          
          scored_records[[length(scored_records) + 1]] <- scored_plate
        }
      }
    }
  }
  
  if (length(scored_records) == 0) return(data.frame())
  
  combined_scores <- dplyr::bind_rows(scored_records)
  
  # FDR
  combined_sig <- combined_scores %>%
    dplyr::group_by(Treatment, Metric) %>%
    dplyr::mutate(
      mean_log_cg = mean(log_cg_score, na.rm = TRUE),
      sd_log_cg = sd(log_cg_score, na.rm = TRUE),
      z_score = dplyr::if_else(is.na(sd_log_cg) | sd_log_cg == 0, NA_real_, (log_cg_score - mean_log_cg) / sd_log_cg),
      p_value = 2 * pnorm(-abs(z_score))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      adj_p = p.adjust(p_value, method = "BH"),
      Is_Significant = ifelse(adj_p < 0.05 & abs(log_cg_score) > 1.0, "Y", "N")
    )
  
  return(combined_sig)
}
