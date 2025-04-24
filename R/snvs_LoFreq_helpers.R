noise_snvs_LoFreq <- function(path, merged_file, gt_load, gt_tv){
    
    LoFreq_somatic <- load_LoFreq_vcf(path, merged_file)
    LoFreq_somatic_snvs <-select_snvs(LoFreq_somatic)
    LoFreq_somatic_snvs <- LoFreq_somatic_snvs[,c("POS", "REF", "ALT", "DP", "AD", "AF" ,"mut" )]
    colnames(LoFreq_somatic_snvs) <- c("POS", "REF",  "ALT",  "DP", "AD", "AF","mut" )
    LoFreq_somatic_snvs$AF = as.numeric(LoFreq_somatic_snvs$AF)######
    LoFreq_somatic_snvs <- LoFreq_somatic_snvs[!mut %in% gt_tv$mut]
    
    fp_var = define_fp(LoFreq_somatic_snvs, gt_load)
    fn_var = define_fn(LoFreq_somatic_snvs, gt_load)
    tp_var = define_tp(LoFreq_somatic_snvs, gt_load)
    
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