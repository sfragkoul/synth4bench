#'A script, written in R, that calls the appropriate functions 
#'to perform the comparison between the ground truth and the caller.
#'
#' Input files: ground truth files, caller vcf file
#'
#' Output files:  tsv files with the comparison between the ground truth 
#' and the caller.
#'
#' Authors: Nikos Pechlivanis(github:npechl), Stella Fragkouli(sfragkoul)
#' 


#!/usr/bin/env Rscript
message("Loading Libraries...")
source("R/libraries.R")

source("R/common_helpers.R")
source("R/indels_common_helpers.R")
source("R/snvs_common_helpers.R")

source("R/common_helpers_gatk.R")
source("R/indels_gatk_helpers.R")
source("R/snvs_gatk_helpers.R")

source("R/common_helpers_Freebayes.R")
source("R/indels_Freebayes_helpers.R")
source("R/snvs_Freebayes_helpers.R")

source("R/common_helpers_VarDict.R")
source("R/indels_VarDict_helpers.R")
source("R/snvs_VarDict_helpers.R")

source("R/common_helpers_VarScan.R")
source("R/indels_VarScan_helpers.R")
source("R/snvs_VarScan_helpers.R")

source("R/common_helpers_LoFreq.R")
source("R/indels_LoFreq_helpers.R")
source("R/snvs_LoFreq_helpers.R")

#Parse arguments from command line
options <- list(
  # make_option(c("-v", "--vcf_path"), 
  #             action = "store", 
  #             type = "character", 
  #             help="Directory path where VCF files are located."),
  
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
              help="Path of working directory were all files are located and the results will be generated."),
  
  make_option(c("-m", "--merged_file"), 
              action = "store", 
              type = "character", 
              help="Indicate the name given to the final merged ground truth file.")
  
)

arguments <- parse_args(OptionParser(option_list = options))
arguments$vcf_path <- arguments$working_directory

#SNVS True Variants------------------------------------------------------------
# print("Begin SNVs True Variant Analysis")
# 
# # check if Merged_snvs_GT.tsv exists
# output_file <- file.path(arguments$working_directory,
#                          paste0(arguments$merged_file, "_snvs_TV.tsv"))
# 
# if (!file.exists(output_file)) {
#     gt <- gt_analysis(seq_len(arguments$runs),
#                       arguments$working_directory,
#                       arguments$merged_file)
#     
#     fwrite(gt, output_file, row.names = FALSE, quote = FALSE, sep = "\t")
# } else {
#    message("File already exists: ", output_file)
#    gt <- fread(output_file)
# }
# 
# 
# 
# snvs <- read_vcf_snvs_TP(arguments$vcf_path,
#                                    arguments$caller,
#                                    gt,
#                                    arguments$merged_file)
# 
# out_snvs = snvs$vcf_snvs_cleaned
# recall = snvs$recall
# 
# print(paste("Recall score for True Variants detection:", round(recall, 2)))
# 
# 
# fwrite(
#   out_snvs, paste0(arguments$working_directory,
#                          "/",
#                          arguments$merged_file,
#                          "_",
#                          arguments$caller, "_snvs_TV.tsv"),
# 
#   row.names = FALSE, quote = FALSE, sep = "\t"
# )

#SNVS Noise variants-----------------------------------------------------------
print("Begin SNVs FP Variant Analysis")

gt_load <- load_gt_report(arguments$vcf_path,
                         arguments$merged_file)$snvs



noise_snvs <- noise_variants(arguments$vcf_path,
                                   arguments$caller,
                                   arguments$merged_file,
                                   gt_load)

out_noise_snvs = rbind(noise_snvs$tp, noise_snvs$fp, noise_snvs$fn)
out_noise_snvs$POS = as.numeric(out_noise_snvs$POS)
out_noise_snvs = out_noise_snvs[order(POS)]

fwrite(
    out_noise_snvs, paste0(arguments$working_directory, "/",
                           arguments$merged_file, "_",
                           arguments$caller, "_snvs_Noise.tsv"),
    row.names = FALSE, quote = FALSE, sep = "\t"
)



#INDELs------------------------------------------------------------------
# print("Begin INDELs Variant Analysis")
# 
# pick_gt_stdz <- gt_stdz_indels(arguments$vcf_path,
#                               arguments$merged_file)
# 
# print("Begin TP INDELs Variant Analysis")
# tp_indels <- call_tp_indels(arguments$vcf_path,
#                            arguments$caller,
#                            arguments$merged_file,
#                            pick_gt_stdz)
# fwrite(
#         tp_indels, paste0(arguments$working_directory, "/",
#                          arguments$merged_file, "_",
#                          arguments$caller, "_indels_TP.tsv"),
#   row.names = FALSE, quote = FALSE, sep = "\t"
# )
# 
# print("Begin FN INDELs Variant Analysis")
# fn_indels <- call_fn_indels(arguments$vcf_path,
#                            arguments$caller,
#                            arguments$merged_file,
#                            pick_gt_stdz)
# fwrite(
#   fn_indels, paste0(arguments$working_directory, "/",
#                     arguments$merged_file, "_",
#                     arguments$caller, "_indels_FN.tsv"),
#   row.names = FALSE, quote = FALSE, sep = "\t"
# )
# 
# print("Begin FP INDELs Variant Analysis")
# fp_indels <- call_fp_indels(arguments$vcf_path,
#                            arguments$caller,
#                            arguments$merged_file,
#                            pick_gt_stdz)
# fwrite(
#   fp_indels, paste0(arguments$working_directory, "/",
#                     arguments$merged_file, "_",
#                     arguments$caller, "_indels_FP.tsv"),
#   row.names = FALSE, quote = FALSE, sep = "\t"
# )
# 
# 
# print("Analysis Completed")








