# Function to identify FN variants with specific categories
categorize_fns <- function(caller, fn_var) {

    caller$POS = as.numeric(caller$POS)
    fn_var$POS = as.numeric(fn_var$POS)
    colnames(fn_var) = c("POS","REF", "Ground Truth DP",  "ALT",
                         "Count", "Ground Truth AF","mut","type")
    #Same POS
    same_POS <- merge(fn_var, caller, by = c("POS"))
    fn_var[, category := ifelse(POS %in% same_POS$POS, "diff REF", "not exist")]
    
    #Same POS & REF
    same_POS_REF <- merge(fn_var, caller, by = c("POS", "REF"))
    fn_var[POS %in% same_POS_REF$POS & REF %in% same_POS_REF$REF, 
           category := "diff ALT"]
    
    return(fn_var)
}


#-------
source("R/libraries.R")

caller <- fread("caller.tsv", header = TRUE, sep = "\t") 
fn_var <- fread("fn_var.tsv", header = TRUE, sep = "\t") 

new = categorize_fns(caller, fn_var)


# fwrite(
#     fn_var, "fn_var.tsv",
# 
#         row.names = FALSE, quote = FALSE, sep = "\t"
#     )
