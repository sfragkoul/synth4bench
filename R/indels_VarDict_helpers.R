call_indels_VarDict <- function(path, merged_file, pick_gt_stdz){
    
    VarDict_somatic <- load_VarDict_vcf(path, merged_file)
    VarDict_indels <-select_indels(VarDict_somatic)
    VarDict_indels <- VarDict_indels[,c("POS", "REF", "ALT", "DP", "VD", "AF" ,"mut" )]
    colnames(VarDict_indels) <- c("POS", "REF",  "ALT",  "DP", "AD", "AF","mut" )
    VarDict_indels$AF = as.numeric(VarDict_indels$AF)
    VarDict_indels$POS = as.numeric(VarDict_indels$POS)
    
    #TP ------------------
    tp_var = define_tp(VarDict_indels, pick_gt_stdz)
    tp_var$category = NaN
    
    #FN ------------------
    fn_var = define_fn(VarDict_indels, pick_gt_stdz)
    #categorize INDELs
    ##Same POS
    same_POS <- merge(fn_var, VarDict_indels, by = c("POS"))
    fn_var[, category := ifelse(POS %in% same_POS$POS, "diff REF", "diff POS")]
    ##Same POS & REF
    same_POS_REF <- merge(fn_var, VarDict_indels, by = c("POS", "REF"))
    ##Update only rows where POS and REF match
    fn_var[POS %in% same_POS_REF$POS & REF %in% same_POS_REF$REF, 
           category := "diff ALT"]
    
    #FP ------------------
    fp_var = define_fp(VarDict_indels, pick_gt_stdz)
    #categorize INDELs
    ##Same POS
    same_POS <- merge(fp_var, pick_gt_stdz, by = c("POS"))
    fp_var[, category := ifelse(POS %in% same_POS$POS, "diff REF", "diff POS")]
    
    ##Same POS & REF
    same_POS_REF <- merge(fp_var, pick_gt_stdz, by = c("POS", "REF"))
    ##Update only rows where POS and REF match
    fp_var[POS %in% same_POS_REF$POS & REF %in% same_POS_REF$REF, 
           category := "diff ALT"]
    
    
    recall = nrow(tp_var)/(nrow(tp_var) + nrow(fn_var))
    precision = nrow(tp_var)/(nrow(tp_var) + nrow(fp_var))
    
    return(list(
        "fp" = fp_var,
        "fn" = fn_var,
        "tp" = tp_var,
        "indel_recall" = recall,
        "indel_precision" = precision)
    )
}