

#!/usr/bin/env Rscript
source("R/libraries.R")
source("R/common_helpers.R")
source("R/helpers_freebayes.R")
source("R/helpers_gatk.R")
source("R/helpers_LoFreq.R")
source("R/helpers_VarDict.R")
source("R/helpers_VarScan.R")

#Parse arguments from command line
options <- list(
  make_option(c("-t", "--gt_comparison"), 
              action = "store", 
              type = "character", 
              help="Directory path where Ground Truth vs Caller file tsv file is located."),
  
  make_option(c("-v", "--vcf_path"), 
              action = "store", 
              type = "character", 
              help="Directory path where VCF files are located."),
  
  make_option(c("-g", "--gt_path"), 
              action = "store", 
              type = "character", 
              help="Directory path where ground truth vcf file is located."),
  
  make_option(c("-c", "--caller"), 
              action = "store", 
              type = "character", 
              help="Choose caller name (Freebayes, Mutect2, LoFreq, VarDict, VarScan)"),
  
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

#print(arguments)

plots <- plot_synth4bench(arguments$gt_comparison, arguments$vcf_path, arguments$gt_path, arguments$caller, arguments$merged_file)

dir.create(paste0(arguments$working_directory, "/Plots"))

ggsave(
  plot = plots[[1]], filename = paste0(arguments$working_directory, "/Plots/Poster_", arguments$caller, ".png"),
  width = 16, height = 12, units = "in", dpi = 600
)

ggsave(
  plot = plots[[2]], filename = paste0(arguments$working_directory, "/Plots/Venn_", arguments$caller, ".png"),
  width = 8, height = 8, units = "in", dpi = 600
)