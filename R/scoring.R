#' @import dplyr tidyr
NULL

#' Standardized Scoring Engine
#'
#' @param pheno_df data.frame containing aggregated phenotypic metrics across a batch
#' @param strategy_type character, either "NoDrug_Plates" or "NoDrug_Wells"
#' @param metrics character vector of target metrics
#' @param control_drug character, name of control treatment
#' @return data.frame of scores
#' @export
score_batch_dataset <- function(pheno_df, strategy_type = c("NoDrug_Plates", "NoDrug_Wells"), metrics, control_drug = "NoDrug") {
  strategy_type <- match.arg(strategy_type)
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
      valid_metric_df <- metric_df %>%
        dplyr::filter(!(Treatment != control_drug & Is_NoDrug_Control == "Y"))
      
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
