source("R/libraries.R")
source("R/common_helpers.R")
#source("R/helpers_gatk.R")
source("R/helpers_VarDict.R")

pick_gt_stdz = gt_stdz_indels("results", "Merged")

#Caller------------------------------------------------------------------------
#TP
tp_indels_gatk = final_tp_indels_gatk("results", "Merged", pick_gt_stdz)

#FN
call_fn_indels_gatk <- function(path, merged_file, pick_gt_stdz){
    #function to output categorized FN indels
    fn_indels_gatk = final_fn_indels_gatk(path, merged_file, pick_gt_stdz)
    Mutect2_somatic = load_gatk_vcf(path, merged_file)
    Mutect2_indels = select_indels(Mutect2_somatic)
    fn_indels_gatk_categories = categorize_fns_gatk(Mutect2_indels, fn_indels_gatk)
    
    return(fn_indels_gatk_categories)
}
new_fn = call_fn_indels_gatk("results", "Merged")

call_fp_indels_gatk <- function(path, merged_file){
    #function to output categorized FP indels
    gt_all = load_gt_report_indels(path, merged_file)$all |> standardize_indels()
    fp_indels_gatk = final_fp_indels_gatk(path, merged_file, pick_gt_stdz, gt_all)
    fp_indels_gatk_categories = categorize_fps_gatk(pick_gt_stdz, fp_indels_gatk)
    
    return(fp_indels_gatk_categories)
}
new_fp = call_fp_indels_gatk("results", "Merged")
#------------------------------------------------------------------------------
#TP
final_tp_indels_VarDict <- function(path, merged_file, pick_gt_stdz){
    #function to identify TP indels
    VarDict_somatic_indels <- load_VarDict_vcf(path, merged_file) |> select_indels()
    tp_var = define_tp(VarDict_somatic_indels, pick_gt_stdz)
    return(tp_var)
}

new_tp = final_tp_indels_VarDict("results", "Merged", pick_gt_stdz)

#FN
final_fn_indels_VarDict <- function(path, merged_file, pick_gt_stdz){
    #function to identify FN indels
    VarDict_somatic_indels <- load_VarDict_vcf(path, merged_file) |> select_indels()
    fn_var = define_fn(VarDict_somatic_indels, pick_gt_stdz)
    colnames(fn_var) = c("POS", "Ground Truth REF", "Ground Truth DP", 
                         "Ground Truth ALT", "Count", "Ground Truth AF", "mut", "type")
    return(fn_var)
}

categorize_fns_VarDict <- function(caller, fn_var) {
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


call_fn_indels_VarDict <- function(path, merged_file, pick_gt_stdz){
    #function to output categorized FN indels
    fn_indels_VarDict = final_fn_indels_VarDict(path, merged_file, pick_gt_stdz)
    VarDict_indels = load_VarDict_vcf(path, merged_file) |> select_indels()
    fn_indels_VarDict_categories = categorize_fns_VarDict(VarDict_indels, fn_indels_VarDict)
    
    return(fn_indels_VarDict_categories)
}

new_fn = call_fn_indels_VarDict("results", "Merged", pick_gt_stdz)



#FP
final_fp_indels_VarDict <- function(path, merged_file, pick_gt, gt_all){
    #function to identify FP indels
    VarDict_somatic_indels <- load_VarDict_vcf(path, merged_file) |> select_indels()
    fp_var = fp_snvs_VarDict(VarDict_somatic_indels, pick_gt, gt_all)
    return(fp_var)
}


categorize_fps_VarDict <- function(pick_gt_stdz, fp_indels_VarDict) {
    #function to identify FP categories
    pick_gt_stdz$POS = as.numeric(pick_gt_stdz$POS)
    fp_indels_VarDict$POS = as.numeric(fp_indels_VarDict$POS)
    
    colnames(fp_indels_VarDict) = c("CHROM", "POS", "ID", "REF", 
                                      "ALT", "VarDict QUAL", "VarDict FILTER", "VarDict DP", 
                                      "VarDict AF", "mut", "Ground Truth DP","DP Percentage", "type")
    #Same POS
    same_POS <- merge(fp_indels_VarDict, pick_gt_stdz, by = c("POS"))
    fp_indels_VarDict[, category := ifelse(POS %in% same_POS$POS, "diff REF", "not exist")]
    
    #Same POS & REF
    same_POS_REF <- merge(fp_indels_VarDict, pick_gt_stdz, by = c("POS", "REF"))
    # Update only rows where POS and REF match
    fp_indels_VarDict[POS %in% same_POS_REF$POS & REF %in% same_POS_REF$REF, 
                        category := "diff ALT"]
    
    return(fp_indels_VarDict)
}

call_fp_indels_VarDict <- function(path, merged_file, pick_gt_stdz){
    #function to output categorized FP indels
    gt_all = load_gt_report_indels(path, merged_file)$all |> standardize_indels()
    fp_indels_VarDict = final_fp_indels_VarDict(path, merged_file, pick_gt_stdz, gt_all)
    fp_indels_VarDict_categories = categorize_fps_VarDict(pick_gt_stdz, fp_indels_VarDict)
    
    return(fp_indels_VarDict_categories)
}

new_fp = call_fp_indels_VarDict("results", "Merged", pick_gt_stdz)


path = "results"
merged_file = "Merged"
caller =  "VarDict"

p = circular_plot_VarDict(path, merged_file, caller)

fwrite(
    new_tp, paste0("Merged_VarDict_indels_TP.tsv"),
    row.names = FALSE, quote = FALSE, sep = "\t"
)

fwrite(
    new_fp, paste0("Merged_VarDict_indels_FP.tsv"),
    row.names = FALSE, quote = FALSE, sep = "\t"
)

fwrite(
    new_fn, paste0("Merged_VarDict_indels_FN.tsv"),
    row.names = FALSE, quote = FALSE, sep = "\t"
)



