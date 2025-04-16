
noise_snvs_gatk <- function(path, merged_file, gt_load, gt_tv){
    
    Mutect2_somatic <- load_gatk_vcf(path, merged_file)
    Mutect2_somatic_snvs <-select_snvs(Mutect2_somatic)
    Mutect2_somatic_snvs[, AD := as.numeric(sapply(strsplit(gt_AD, ","), function(x) x[2]))]#######
    Mutect2_somatic_snvs <- Mutect2_somatic_snvs[,c("POS", "REF", "ALT", "gt_DP", "AD", "gt_AF" ,"mut" )]
    colnames(Mutect2_somatic_snvs) <- c("POS", "REF",  "ALT",  "DP", "AD", "AF","mut" )
    Mutect2_somatic_snvs$AF = as.numeric(Mutect2_somatic_snvs$AF)######
    Mutect2_somatic_snvs <- Mutect2_somatic_snvs[!mut %in% gt_tv$mut]
    
    fp_var = define_fp(Mutect2_somatic_snvs, gt_load)
    fn_var = define_fn(Mutect2_somatic_snvs, gt_load)
    tp_var = define_tp(Mutect2_somatic_snvs, gt_load)
    
    return(list(
        "fp" = fp_var,
        "fn" = fn_var,
        "tp" = tp_var)
    )
}