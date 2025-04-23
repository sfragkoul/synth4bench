noise_snvs_Freebayes <- function(path, merged_file, gt_load, gt_tv){
    
    Freebayes_somatic <- load_Freebayes_vcf(path, merged_file)
    Freebayes_somatic_snvs <-select_snvs(Freebayes_somatic)
    Freebayes_somatic_snvs <- Freebayes_somatic_snvs[,c("POS", "REF", "ALT", "DP", "AO", "AF" ,"mut" )]
    colnames(Freebayes_somatic_snvs) <- c("POS", "REF",  "ALT",  "DP", "AD", "AF","mut" )
    Freebayes_somatic_snvs$AF = as.numeric(Freebayes_somatic_snvs$AF)
    Freebayes_somatic_snvs <- Freebayes_somatic_snvs[!mut %in% gt_tv$mut]
    
    fp_var = define_fp(Freebayes_somatic_snvs, gt_load)
    fn_var = define_fn(Freebayes_somatic_snvs, gt_load)
    tp_var = define_tp(Freebayes_somatic_snvs, gt_load)
    
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