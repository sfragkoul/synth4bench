noise_snvs_VarDict <- function(path, merged_file, gt_load, gt_tv){
    
    VarDict_somatic <- load_VarDict_vcf(path, merged_file)
    VarDict_somatic_snvs <-select_snvs(VarDict_somatic)
    VarDict_somatic_snvs <- VarDict_somatic_snvs[,c("POS", "REF", "ALT", "DP", "VD", "AF" ,"mut" )]
    colnames(VarDict_somatic_snvs) <- c("POS", "REF",  "ALT",  "DP", "AD", "AF","mut" )
    VarDict_somatic_snvs$AF = as.numeric(VarDict_somatic_snvs$AF)######
    VarDict_somatic_snvs <- VarDict_somatic_snvs[!mut %in% gt_tv$mut]
    
    fp_var = define_fp(VarDict_somatic_snvs, gt_load)
    fn_var = define_fn(VarDict_somatic_snvs, gt_load)
    tp_var = define_tp(VarDict_somatic_snvs, gt_load)
    
    recall = nrow(tp_var)/(nrow(tp_var) + nrow(fn_var))
    precision = nrow(tp_var)/(nrow(tp_var) + nrow(fp_var))
    
    return(list(
        "fp" = fp_var,
        "fn" = fn_var,
        "tp" = tp_var,
        "noise_recall" = recall,
        "noise_precision" = precision)
    )
}