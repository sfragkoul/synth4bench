# Function to identify FN variants with specific categories
categorize_fns <- function(caller, fn_var) {
    
    caller = Mutect2_indels 
    fn_var = fn_indels_gatk

    return()
}
 


caller$POS = as.numeric(caller$POS)  
fn_var$POS = as.numeric(fn_var$POS) 

colnames(fn_var) = c("POS","REF", "Ground Truth DP",  "ALT",
                                           "Count", "Ground Truth AF","mut","type") 


#Same POS
same_POS <- merge(fn_var, caller, by = c("POS"))
fn_var[, category := ifelse(POS %in% same_POS$POS, "same POS", "diff POS")]

#Same POS & REF
same_POS_REF <- merge(fn_var, caller, by = c("POS", "REF"))
fn_var[, category := ifelse(POS %in% same_POS_REF$POS &
                            REF %in% same_POS_REF$REF, "same POS & REF",,
                            category)]  # Keep existing category if it's "same POS"

#Same POS & REF & ALT
same_POS_REF_ALT <- merge(fn_var, caller, by = c("POS", "REF", "ALT"))
