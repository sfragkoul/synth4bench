

#!/usr/bin/env Rscript
source("R/libraries.R")
source("R/common_helpers.R")

#Parse arguments from command line
options <- list(
    make_option(c("-f", "--reference"), action = "store", type = "character", help="Choose caller name (freebayes, gatk, LoFreq, VarDict, VarScan)"),
    make_option(c("-r", "--runs"), action = "store", type = "integer", help="Choose caller name (freebayes, gatk, LoFreq, VarDict, VarScan)"),
    make_option(c("-w", "--working_directory"), action = "store", type = "character", help="Choose caller name (freebayes, gatk, LoFreq, VarDict, VarScan)"),
    make_option(c("-c", "--caller"), action = "store", type = "character", help="Choose caller name (freebayes, gatk, LoFreq, VarDict, VarScan)")
)

arguments <- parse_args(OptionParser(option_list = options))

df_list <- report_varbp(arguments$runs, arguments$working_directory, arguments$reference)


for(i in length(df_list)) {
    
    fwrite(
        df_list[[i]][["report"]], 
        paste0(arguments$working_directory, "/", i, "/", i, "_position_report.csv"),
        row.names = FALSE, quote = TRUE, sep = ","
    )
    
    fwrite(
        df_list[[i]][["report2"]], 
        paste0(arguments$working_directory, "/", i, "/", i, "_position_report2.csv"),
        row.names = FALSE, quote = TRUE, sep = ","
    )
    
}


df_list2 <- explore_mut_pos(arguments$runs, arguments$working_directory, arguments$caller)

