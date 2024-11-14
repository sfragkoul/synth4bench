source("R/libraries.R")
source("R/common_helpers.R")
source("R/helpers_gatk.R")

pick_gt_stdz = gt_stdz_indels("results", "Merged")

#Mutect2-----------------------------------------------------------------------
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



# fwrite(
#     new_fn, paste0("fn_var_new.tsv"),
#     row.names = FALSE, quote = FALSE, sep = "\t"
# )
# 
# fwrite(
#     new_fp, paste0("fp_var_new.tsv"),
#     row.names = FALSE, quote = FALSE, sep = "\t"
# )



