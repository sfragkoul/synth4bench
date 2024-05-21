


#!/usr/bin/env Rscript
library(optparse)

#Parse arguments from command line
options <- list(
  make_option(c("-p", "--vcf_path"), action = "store", type = "character", help="Working directory path."),
  make_option(c("-c", "--caller"), action = "store", type = "character", help="File(s) extension.")
)


arguments <- parse_args(OptionParser(option_list = options))

# TODO ----------------

# source("helpers.R")

# construct gt object

# ---------------------------

source("read_vcf.R")

source("helpers_freebayes.R")
source("helpers_gatk.R")


read_vcf(arguments$vcf_path, arguments$caller, gt)