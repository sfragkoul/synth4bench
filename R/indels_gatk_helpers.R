
categorize_fns_gatk <- function(caller, fn_var) {
    #function to identify FN categories
    
    caller$POS = as.numeric(caller$POS)
    fn_var$POS = as.numeric(fn_var$POS)
    colnames(fn_var) = c("POS","REF", "Ground Truth DP",  "ALT",
                         "Count", "Ground Truth AF","mut","type")
    #Same POS
    same_POS <- merge(fn_var, caller, by = c("POS"))
    fn_var[, category := ifelse(POS %in% same_POS$POS, "diff REF", "not exist")]
    
    #Same POS & REF
    same_POS_REF <- merge(fn_var, caller, by = c("POS", "REF"))
    # Update only rows where POS and REF match
    fn_var[POS %in% same_POS_REF$POS & REF %in% same_POS_REF$REF, 
           category := "diff ALT"]
    
    return(fn_var)
}

categorize_fps_gatk <- function(pick_gt_stdz, fp_indels_gatk) {
    #function to identify FP categories
    pick_gt_stdz$POS = as.numeric(pick_gt_stdz$POS)
    fp_indels_gatk$POS = as.numeric(fp_indels_gatk$POS)
    
    colnames(fp_indels_gatk) = c("CHROM", "POS", "ID", "REF", "ALT", "Mutect2 QUAL",
                                 "Mutect2 FILTER", "key", "Indiv", "Mutect2 AD", 
                                 "Mutect2 AF", "Mutect2 DP", "gt_F1R2", "gt_F2R1", 
                                 "gt_FAD", "gt_GQ", "gt_GT", "gt_PGT", "gt_PID", 
                                 "gt_PL" , "gt_PS", "gt_SB", "gt_GT_alleles", 
                                 "mut", "Ground Truth DP","DP Percentage", "type")
    #Same POS
    same_POS <- merge(fp_indels_gatk, pick_gt_stdz, by = c("POS"))
    fp_indels_gatk[, category := ifelse(POS %in% same_POS$POS, "diff REF", "not exist")]
    
    #Same POS & REF
    same_POS_REF <- merge(fp_indels_gatk, pick_gt_stdz, by = c("POS", "REF"))
    # Update only rows where POS and REF match
    fp_indels_gatk[POS %in% same_POS_REF$POS & REF %in% same_POS_REF$REF, 
                   category := "diff ALT"]
    
    return(fp_indels_gatk)
}

final_fp_indels_gatk <- function(path, merged_file, pick_gt, gt_all){
    #function to identify FP indels
    Mutect2_somatic <- load_gatk_vcf(path, merged_file)
    Mutect2_somatic_indels <-select_indels(Mutect2_somatic)
    fp_var = fp_snvs_gatk(Mutect2_somatic_indels, pick_gt, gt_all)
    return(fp_var)
}

final_fn_indels_gatk <- function(path, merged_file, pick_gt){
    #function to identify FN indels
    Mutect2_somatic <- load_gatk_vcf(path, merged_file)
    Mutect2_somatic_indels <-select_indels(Mutect2_somatic)
    fn_var = define_fn(Mutect2_somatic_indels, pick_gt)
    colnames(fn_var) = c("POS", "Ground Truth REF", "Ground Truth DP", 
                         "Ground Truth ALT", "Count", "Ground Truth AF", "mut", "type")
    return(fn_var)
}

final_tp_indels_gatk <- function(path, merged_file, pick_gt){
    #function to identify TP indels
    Mutect2_somatic <- load_gatk_vcf(path, merged_file)
    Mutect2_somatic_indels <-select_indels(Mutect2_somatic)
    tp_var = define_tp(Mutect2_somatic_indels, pick_gt)
    return(tp_var)
}

call_fn_indels_gatk <- function(path, merged_file, pick_gt_stdz){
  #function to output categorized FN indels
  fn_indels_gatk = final_fn_indels_gatk(path, merged_file, pick_gt_stdz)
  Mutect2_somatic = load_gatk_vcf(path, merged_file)
  Mutect2_indels = select_indels(Mutect2_somatic)
  fn_indels_gatk_categories = categorize_fns_gatk(Mutect2_indels, fn_indels_gatk)
  
  return(fn_indels_gatk_categories)
}

call_fp_indels_gatk <- function(path, merged_file, pick_gt_stdz){
  #function to output categorized FP indels
  gt_all = load_gt_report_indels(path, merged_file)$all |> standardize_indels()
  fp_indels_gatk = final_fp_indels_gatk(path, merged_file, pick_gt_stdz, gt_all)
  fp_indels_gatk_categories = categorize_fps_gatk(pick_gt_stdz, fp_indels_gatk)
  
  return(fp_indels_gatk_categories)
}
