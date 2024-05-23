

#!/usr/bin/env Rscript
library(optparse)

#Parse arguments from command line
options <- list(
  make_option(c("-p", "--vcf_path"), action = "store", type = "character", help="VCF files directory path."),
  make_option(c("-c", "--caller"), action = "store", type = "character", help="Choose caller name (freebayes, gatk, LoFreq, VarDict, VarScan)"),
  make_option(c("-r", "--runs"), action = "store", type = "integer", help="Choose caller name (freebayes, gatk, LoFreq, VarDict, VarScan)"),
  make_option(c("-w", "--working_directory"), action = "store", type = "character", help="Choose caller name (freebayes, gatk, LoFreq, VarDict, VarScan)")
  
)


arguments <- parse_args(OptionParser(option_list = options))

print(arguments)

# PART 1 ----------------

source("R/common_helpers.R")

gt <- gt_analysis(arguments$runs, arguments$working_directory)

# PART 2 ----------------------

source("R/2_downstream_analysis.R")

out_df <- read_vcf(arguments$vcf_path, arguments$caller, gt)

fwrite(
  out_df, paste0(arguments$working_directory, "/Ground_truth_vs_", arguments$caller, ".clean_norm.tsv"),
  row.names = FALSE, quote = FALSE, sep = "\t"
)

