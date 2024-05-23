

#!/usr/bin/env Rscript
library(optparse)

#Parse arguments from command line
options <- list(
  make_option(c("-m", "--gt_comparison"), action = "store", type = "character", help="VCF files directory path."),
  make_option(c("-v", "--vcf_path"), action = "store", type = "character", help="Choose caller name (freebayes, gatk, LoFreq, VarDict, VarScan)"),
  make_option(c("-g", "--gt_path"), action = "store", type = "character", help="Choose caller name (freebayes, gatk, LoFreq, VarDict, VarScan)"),
  make_option(c("-c", "--caller"), action = "store", type = "character", help="Choose caller name (freebayes, gatk, LoFreq, VarDict, VarScan)"),
  make_option(c("-w", "--working_directory"), action = "store", type = "character", help="Choose caller name (freebayes, gatk, LoFreq, VarDict, VarScan)")
  
)


arguments <- parse_args(OptionParser(option_list = options))

print(arguments)

source("R/2_downstream_analysis.R")

plots <- plot_synth4bench(arguments$gt_comparison, arguments$vcf_path, arguments$gt_path, arguments$caller)

dir.create(paste0(arguments$working_directory, "/Plots"))

ggsave(
  plot = multi, filename = paste0(arguments$working_directory, "/Plots/Poster_", arguments$caller, ".png"),
  width = 16, height = 12, units = "in", dpi = 600
)

ggsave(
  plot = out4, filename = paste0(arguments$working_directory, "/Plots/Venn_", arguments$caller, ".png"),
  width = 8, height = 8, units = "in", dpi = 600
)