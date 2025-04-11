
fp_snvs_LoFreq <- function(LoFreq_somatic_snvs, pick_gt, gt_all){
    #find LoFreq FP variants
    fp_var = define_fp(LoFreq_somatic_snvs, pick_gt)
    fp_var$AF = as.numeric(fp_var$AF)
    colnames(fp_var) = c("CHROM", "POS","ID", "LoFreq REF",	
                         "LoFreq ALT", "LoFreq QUAL",	"LoFreq FILTER",
                         "LoFreq DP", "LoFreq AF", "mut")
    
    #find DP of FP variants'  location in GT
    tmp = gt_all[which(POS %in% unique(fp_var$POS))]
    a = unique(tmp, by = "POS")
    #to include the presence multiple variants in a POS
    index = match(fp_var$POS, a$POS)
    fp_var$`Ground Truth DP` = a[index]$DP
    fp_var$`DP Percentage` = fp_var$`LoFreq DP`/fp_var$`Ground Truth DP`
    fp_var$type = "FP"
    return(fp_var)
}

final_fp_snvs_LoFreq <- function(path, merged_file, pick_gt, gt_all){
    
    LoFreq_somatic <- load_LoFreq_vcf(path, merged_file)
    LoFreq_somatic_snvs <-select_snvs(LoFreq_somatic)
    fp_var = fp_snvs_LoFreq(LoFreq_somatic_snvs, pick_gt, gt_all)
    
    return(fp_var)
}

final_fn_snvs_LoFreq <- function(path, merged_file, pick_gt){
    
    LoFreq_somatic <- load_LoFreq_vcf(path, merged_file)
    LoFreq_somatic_snvs <-select_snvs(LoFreq_somatic)
    fn_var = define_fn(LoFreq_somatic_snvs, pick_gt)
    colnames(fn_var) = c("POS", "Ground Truth REF", "Ground Truth DP", 
                         "Ground Truth ALT", "Count", "Ground Truth AF", "mut", "type")
    
    return(fn_var)
}
