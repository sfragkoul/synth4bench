
categorize_fps_gatk <- function(pick_gt_stdz, fp_indels_gatk) { #!!!! NEW FUNCTION
    #Function to identify FP categories
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


#-------
source("R/libraries.R")

caller -> pick_gt_stdz
fn_var -> fp_indels_gatk



fwrite(
    fp_indels_gatk, "fp_var.tsv",
     
                 row.names = FALSE, quote = FALSE, sep = "\t"
         )




