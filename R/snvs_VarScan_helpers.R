

fp_snvs_VarScan <- function(VarScan_somatic_snvs, pick_gt, gt_all){
    #find VarScan FP variants
    fp_var = define_fp(VarScan_somatic_snvs, pick_gt)
    fp_var$AF = as.numeric(fp_var$AF)
    colnames(fp_var) = c("CHROM", "POS","ID", "VarScan REF",	
                         "VarScan ALT", "VarScan QUAL",	"VarScan FILTER",
                         "VarScan DP", "VarScan AF", "mut")
    
    #find DP of FP variants'  location in GT
    tmp = gt_all[which(POS %in% unique(fp_var$POS))]
    tmp = tmp[nchar(tmp$REF) == nchar(tmp$ALT)]
    a = unique(tmp, by = "POS")
    #to include the presence multiple variants in a POS
    index = match(fp_var$POS, a$POS)
    fp_var$`Ground Truth DP` = a[index]$DP
    fp_var$`DP Percentage` = fp_var$`VarScan DP`/fp_var$`Ground Truth DP`
    fp_var$type = "FP"
    return(fp_var)
}

final_fp_snvs_VarScan <- function(path, merged_file, pick_gt, gt_all){
    
    VarScan_somatic <- load_VarScan_vcf(path, merged_file)
    VarScan_somatic_snvs <-select_snvs(VarScan_somatic)
    fp_var = fp_snvs_VarScan(VarScan_somatic_snvs, pick_gt, gt_all)
    
    return(fp_var)
}

final_fn_snvs_VarScan <- function(path, merged_file, pick_gt){
    
    VarScan_somatic <- load_VarScan_vcf(path, merged_file)
    VarScan_somatic_snvs <-select_snvs(VarScan_somatic)
    fn_var = define_fn(VarScan_somatic_snvs, pick_gt)
    colnames(fn_var) = c("POS", "Ground Truth REF", "Ground Truth DP", 
                         "Ground Truth ALT", "Count", "Ground Truth AF", "mut", "type")
    
    return(fn_var)
}

