if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "Batch", "Batch_Mean_Control", "Batch_Mean_Target", "Batch_SD_Control",
    "Biofilm_CV", "Biofilm_Mean_OD", "Biofilm_OD", "Biofilm_SD_OD",
    "Biofilm_Score_RB", "Biofilm_Stain_OD", "Blank_Median_OD", "Column",
    "Column_Code", "Control_P", "Control_Value", "Cycle", "DateTime",
    "Drug", "Elapsed_hr", "Exclude_From_Scoring", "Excluded", "FDR",
    "Final_OD", "Flagged", "Gene", "Genotype", "Genotype_Base",
    "Genotype_Suffix", "Group_Var", "Growth_AUC", "Growth_AUC_Rate",
    "Growth_Delta", "Growth_Diff", "Growth_Duration", "Growth_Endpoint_OD",
    "Growth_Endpoint_OD_Raw", "Growth_Failed", "Growth_Initial_OD",
    "Growth_Max_OD", "Growth_Max_Rate", "Growth_OD", "Growth_Time_Duration",
    "Initial_OD", "Interaction_Log2", "Interaction_Score",
    "Is_Contaminated_Blank", "Is_Edge", "Is_NoDrug_Control",
    "Is_NonGrowing_Control", "Is_Not_Growing", "Is_Significant", "Is_Stained",
    "Max_Jump_Index", "Mean_Control_P", "Mean_WT_Control_P", "Media",
    "Median_Delta", "Median_OD", "Metric", "OD", "OD_Corrected", "OD_Delta",
    "OD_Diff", "OD_Raw", "Pass_Blank", "Pass_WT", "Phenotype", "Plate",
    "PlateID", "Plate_Key", "Plate_UUID", "RawFile", "Relative_Biofilm",
    "Rep_Suffix", "Replicate", "ReplicateID", "Replicate_Pair", "Row",
    "Staining_Hour", "Strategy", "Target_Value", "Time_Hours", "Trace_ID",
    "Treatment", "Type", "Type_Display", "WT_Drug_P", "Well", "Well_ID",
    "WellType", "Well_Label", "Well_Type", "Z_Score", "adj_p", "cg_score",
    "fill_val", "label_text", "log_cg_score", "mean_log_cg", "n_bad_blanks",
    "n_mutants", "n_total_blanks", "p_value", "pass_blank", "pass_wt",
    "sd_log_cg", "tooltip_text", "w_control", "w_drug", "w_wt_drug",
    "x_coord", "y_coord", "z_score"
  ))
}
