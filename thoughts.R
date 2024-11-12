# Function to identify FN variants with specific categories
categorize_fns <- function(caller, fn_var) {
    
    caller = Mutect2_indels 
    fn_var = fn_indels_gatk
#     setDT(caller)
#     setDT(fn_var)
#     
#     # For each FN, check the categories
#     categorized_fns <- fn_var[, .(
#         category = if (!POS %in% caller$POS) {
#             "Different_POS"
#         } else {
#             # Subset based on current POS
#             subset_caller <- caller[POS == .SD$POS]
#             # Check if subset is empty for each condition
#             if (nrow(subset_caller) == 0) {
#                 "Different_POS"
#             } else if (!any(`Ground Truth REF` %in% subset_caller$REF) && any(`Ground Truth ALT` %in% subset_caller$ALT)) {
#                 "Different_REF"
#             } else if (any(`Ground Truth REF` %in% subset_caller$REF) && !any(`Ground Truth ALT` %in% subset_caller[REF == .SD$`Ground Truth REF`]$ALT)) {
#                 "Different_ALT"
#             } else {
#                 "Uncategorized"
#             }
#         }
#     ), by = .(POS, `Ground Truth REF`, `Ground Truth ALT`)]
#     
#     return(categorized_fns)
# }
}
 
  

caller$POS = as.numeric(caller$POS)  
fn_var$POS = as.numeric(fn_var$POS) 


# Same POS?
#same_POS <- merge(fn_indels_gatk, caller, by = "POS", all = FALSE)
#same_POS = fn_indels_gatk[which(fn_indels_gatk$POS %in% caller$POS)]
            
#Same POS?
fn_var_same_POS <- fn_indels_gatk[POS %in% caller$POS]
caller_same_POS <- caller[POS %in% fn_indels_gatk$POS]


colnames(fn_var) = c("POS","REF", "Ground Truth DP",  "ALT",
                                           "Count", "Ground Truth AF","mut","type") 


#Same REF?
fn_var_same_POS_REF <- fn_var_same_POS[POS %in% caller_same_POS$POS & 
                                       REF %in% caller_same_POS$REF]

caller_same_POS_REF <- caller_same_POS[POS %in% fn_var_same_POS$POS & 
                                       REF %in% fn_var_same_POS$REF]

same_POS <- merge(fn_var, caller, by = c("POS"))
same_POS_REF <- merge(fn_var, caller, by = c("POS", "REF"))
same_POS_REF_ALT <- merge(fn_var, caller, by = c("POS", "REF", "ALT"))




    
    