
fp_snvs_gatk <- function(Mutect2_somatic_snvs, pick_gt, gt_all){#term snvs is redundant
    #find MUtect2 FP variants
    fp_var = define_fp(Mutect2_somatic_snvs, pick_gt)
    fp_var$gt_AF = as.numeric(fp_var$gt_AF)
    colnames(fp_var) = c("CHROM", "POS","ID", "Mutect2 REF",	
                         "Mutect2 ALT", "Mutect2 QUAL",	"Mutect2 FILTER",
                         "key", "Indiv", "Mutect2 AD", "Mutect2 AF",
                         "Mutect2 DP", "gt_F1R2", "gt_F2R1", "gt_FAD",	
                         "gt_GQ", "gt_GT",	"gt_PGT",	"gt_PID",	"gt_PL",
                         "gt_PS",	"gt_SB",	"gt_GT_alleles", "mut")
    
    #find DP of FP variants'  location in GT
    tmp = gt_all[which(POS %in% unique(fp_var$POS))]
    a = unique(tmp, by = "POS")
    #to include the presence multiple variants in a POS
    index = match(fp_var$POS, a$POS)
    fp_var$`Ground Truth DP` = a[index]$DP
    fp_var$`DP Percentage` = fp_var$`Mutect2 DP`/fp_var$`Ground Truth DP`
    fp_var$type = "FP"
    return(fp_var)
}

final_fp_snvs_gatk <- function(path, merged_file, pick_gt, gt_all){
    
    Mutect2_somatic <- load_gatk_vcf(path, merged_file)
    Mutect2_somatic_snvs <-select_snvs(Mutect2_somatic)
    fp_var <- fp_snvs_gatk(Mutect2_somatic_snvs, pick_gt, gt_all)
    
    return(fp_var)
}

final_fn_snvs_gatk <- function(path, merged_file, pick_gt){
    
    Mutect2_somatic <- load_gatk_vcf(path, merged_file)
    Mutect2_somatic_snvs <- select_snvs(Mutect2_somatic)
    fn_var <- define_fn(Mutect2_somatic_snvs, pick_gt)
    colnames(fn_var) = c("POS", "Ground Truth REF", "Ground Truth DP",
                         "Ground Truth ALT", "Count", "Ground Truth AF",
                         "mut", "type")
    return(fn_var)
}


