call_indels_gatk <- function(path, merged_file, pick_gt_stdz){
    
    Mutect2_somatic <- load_gatk_vcf(path, merged_file)
    Mutect2_indels <-select_indels(Mutect2_somatic)
    Mutect2_indels[, AD := as.numeric(sapply(strsplit(gt_AD, ","), function(x) x[2]))]#######
    Mutect2_indels <- Mutect2_indels[,c("POS", "REF", "ALT", "gt_DP", "AD", "gt_AF" ,"mut" )]
    Mutect2_indels <- Mutect2_indels[,c("POS", "REF", "ALT", "gt_DP", "AD", "gt_AF" ,"mut" )]
    colnames(Mutect2_indels) <- c("POS", "REF",  "ALT",  "DP", "AD", "AF","mut" )
    Mutect2_indels$AF = as.numeric(Mutect2_indels$AF)
    Mutect2_indels$POS = as.numeric(Mutect2_indels$POS)
    
    #TP ------------------
    tp_var = define_tp(Mutect2_indels, pick_gt_stdz)
    tp_var$category = NaN
    
    #FN ------------------
    fn_var = define_fn(Mutect2_indels, pick_gt_stdz)
    #categorize INDELs
    ##Same POS
    same_POS <- merge(fn_var, Mutect2_indels, by = c("POS"))
    fn_var[, category := ifelse(POS %in% same_POS$POS, "diff REF", "not exist")]
    ##Same POS & REF
    same_POS_REF <- merge(fn_var, Mutect2_indels, by = c("POS", "REF"))
    ##Update only rows where POS and REF match
    fn_var[POS %in% same_POS_REF$POS & REF %in% same_POS_REF$REF, 
           category := "diff ALT"]
    
    #FP ------------------
    fp_var = define_fp(Mutect2_indels, pick_gt_stdz)
    #categorize INDELs
    ##Same POS
    same_POS <- merge(fp_var, pick_gt_stdz, by = c("POS"))
    fp_var[, category := ifelse(POS %in% same_POS$POS, "diff REF", "not exist")]
    
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