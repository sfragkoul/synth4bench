noise_snvs_VarScan <- function(path, merged_file, gt_load, gt_tv){
    
    VarScan_somatic <- load_VarScan_vcf(path, merged_file)
    VarScan_somatic_snvs <-select_snvs(VarScan_somatic)
    VarScan_somatic_snvs <- VarScan_somatic_snvs[,c("POS", "REF", "ALT", "DP", "Strands2", "AF" ,"mut" )]
    colnames(VarScan_somatic_snvs) <- c("POS", "REF",  "ALT",  "DP", "AD", "AF","mut" )
    VarScan_somatic_snvs$AF = as.numeric(VarScan_somatic_snvs$AF)######
    VarScan_somatic_snvs <- VarScan_somatic_snvs[!mut %in% gt_tv$mut]
    
    fp_var = define_fp(VarScan_somatic_snvs, gt_load)
    fn_var = define_fn(VarScan_somatic_snvs, gt_load)
    tp_var = define_tp(VarScan_somatic_snvs, gt_load)
    
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