

#!/usr/bin/env Rscript
source("libraries.R")
source("common_helpers.R")
source("helpers_freebayes.R")
source("helpers_gatk.R")
source("helpers_LoFreq.R")
source("helpers_VarDict.R")
source("helpers_VarScan.R")

#Parse arguments from command line
options <- list(
  make_option(c("-p", "--vcf_path"), 
              action = "store", 
              type = "character", 
              help="Directory path where VCF files are located."),
  
  make_option(c("-c", "--caller"), 
              action = "store", 
              type = "character", 
              help="Choose caller name (Freebayes, Mutect2, LoFreq, VarDict, VarScan)"),
  
  make_option(c("-r", "--runs"), 
              action = "store", 
              type = "integer", 
              help="Number of individual runs to produce synthetic data which will then be combined to form the final Merged ground truth file."),
  
  make_option(c("-w", "--working_directory"), 
              action = "store", 
              type = "character", 
              help="Path of working directory."),
  
  make_option(c("-m", "--merged_file"), 
              action = "store", 
              type = "character", 
              help="Indicate the name given to the final merged ground truth file.")
  
)

arguments <- parse_args(OptionParser(option_list = options))

print(arguments)

# PART 1 ----------------

gt <- gt_analysis(arguments$runs, arguments$working_directory, arguments$merged_file)

# PART 2 ---------------------con-

out_df <- read_vcf(arguments$vcf_path, arguments$caller, gt, arguments$merged_file)

fwrite(
  out_df, paste0(arguments$working_directory, "/", arguments$merged_file, "_Ground_truth_vs_", arguments$caller, ".clean_norm.tsv"),
  row.names = FALSE, quote = FALSE, sep = "\t"
)

