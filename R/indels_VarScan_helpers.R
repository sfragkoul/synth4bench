
#INDELS------------------------------------------------------------------------

categorize_fns_VarScan <- function(caller, fn_var) {
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

categorize_fps_VarScan <- function(pick_gt_stdz, fp_indels_VarScan) {
    #function to identify FP categories
    pick_gt_stdz$POS = as.numeric(pick_gt_stdz$POS)
    fp_indels_VarScan$POS = as.numeric(fp_indels_VarScan$POS)
    
    colnames(fp_indels_VarScan) = c("CHROM", "POS", "ID", "REF", 
                                    "ALT", "VarScan QUAL", "VarScan FILTER", "VarScan DP", 
                                    "VarScan AF", "mut", "Ground Truth DP","DP Percentage", "type")
    #Same POS
    same_POS <- merge(fp_indels_VarScan, pick_gt_stdz, by = c("POS"))
    fp_indels_VarScan[, category := ifelse(POS %in% same_POS$POS, "diff REF", "not exist")]
    
    #Same POS & REF
    same_POS_REF <- merge(fp_indels_VarScan, pick_gt_stdz, by = c("POS", "REF"))
    # Update only rows where POS and REF match
    fp_indels_VarScan[POS %in% same_POS_REF$POS & REF %in% same_POS_REF$REF, 
                      category := "diff ALT"]
    
    return(fp_indels_VarScan)
}


final_tp_indels_VarScan <- function(path, merged_file, pick_gt_stdz){
    #function to identify TP indels
    VarScan_somatic_indels <- load_VarScan_vcf(path, merged_file) |> select_indels()
    tp_var = define_tp(VarScan_somatic_indels, pick_gt_stdz)
    return(tp_var)
}

final_fn_indels_VarScan <- function(path, merged_file, pick_gt_stdz){
    #function to identify FN indels
    VarScan_somatic_indels <- load_VarScan_vcf(path, merged_file) |> select_indels()
    fn_var = define_fn(VarScan_somatic_indels, pick_gt_stdz)
    colnames(fn_var) = c("POS", "Ground Truth REF", "Ground Truth DP", 
                         "Ground Truth ALT", "Count", "Ground Truth AF", "mut", "type")
    return(fn_var)
}

final_fp_indels_VarScan <- function(path, merged_file, pick_gt, gt_all){
    #function to identify FP indels
    VarScan_somatic_indels <- load_VarScan_vcf(path, merged_file) |> select_indels()
    fp_var = fp_snvs_VarScan(VarScan_somatic_indels, pick_gt, gt_all)
    return(fp_var)
}


call_fn_indels_VarScan <- function(path, merged_file, pick_gt_stdz){
    #function to output categorized FN indels
    fn_indels_VarScan = final_fn_indels_VarScan(path, merged_file, pick_gt_stdz)
    VarScan_indels = load_VarScan_vcf(path, merged_file) |> select_indels()
    fn_indels_VarScan_categories = categorize_fns_VarScan(VarScan_indels, fn_indels_VarScan)
    
    return(fn_indels_VarScan_categories)
}

call_fp_indels_VarScan <- function(path, merged_file, pick_gt_stdz){
    #function to output categorized FP indels
    gt_all = load_gt_report_indels(path, merged_file)$all |> standardize_indels()
    fp_indels_VarScan = final_fp_indels_VarScan(path, merged_file, pick_gt_stdz, gt_all)
    fp_indels_VarScan_categories = categorize_fps_VarScan(pick_gt_stdz, fp_indels_VarScan)
    
    return(fp_indels_VarScan_categories)
}
