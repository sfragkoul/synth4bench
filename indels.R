source("R/libraries.R")
source("R/common_helpers.R")
#source("R/helpers_gatk.R")
source("R/helpers_freebayes.R")

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

final_tp_indels_freebayes <- function(path, merged_file, pick_gt_stdz){
    #function to identify TP indels
    Freebayes_somatic_indels <- load_Freebayes_vcf(path, merged_file) |> select_indels()
    tp_var = define_tp(Freebayes_somatic_indels, pick_gt_stdz)
    return(tp_var)
}

tp_indels_freebayes = final_tp_indels_freebayes("results", "Merged", pick_gt_stdz)



#FN
call_fn_indels_Freebayes <- function(path, merged_file, pick_gt_stdz){
    #function to output categorized FN indels
    fn_indels_Freebayes = final_fn_indels_Freebayes(path, merged_file, pick_gt_stdz)
    Freebayes_somatic = load_Freebayes_vcf(path, merged_file)
    Freebayes_indels = select_indels(Freebayes_somatic)
    fn_indels_Freebayes_categories = categorize_fns_Freebayes(Freebayes_indels, fn_indels_Freebayes)
    
    return(fn_indels_Freebayes_categories)
}
new_fn = call_fn_indels_Freebayes("results", "Merged")


















