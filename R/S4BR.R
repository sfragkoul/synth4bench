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

#SNVS -------------------------------------------------------------------------
print("Begin SNVs True Variant Analysis")

# check if Merged_snvs_GT.tsv exists
output_file <- file.path(arguments$working_directory,
                         paste0(arguments$merged_file, "_snvs_TV.tsv"))


if (!file.exists(output_file)) {
    gt_tv <- gt_analysis(seq_len(arguments$runs),
                      arguments$working_directory,
                      arguments$merged_file)

    fwrite(gt_tv, output_file, row.names = FALSE, quote = FALSE, sep = "\t")
} else {
   message("File already exists: ", output_file)
    gt_tv <- fread(output_file)
}



snvs <- read_vcf_snvs_TP(arguments$vcf_path,
                                   arguments$caller,
                                   gt_tv,
                                   arguments$merged_file)


fwrite(
    snvs$vcf_snvs_cleaned, paste0(arguments$working_directory,
                         "/",
                         arguments$merged_file,
                         "_",
                         arguments$caller, "_snvs_TV.tsv"),

  row.names = FALSE, quote = FALSE, sep = "\t"
)

print("Begin SNVs Noise Analysis")

gt_load <- load_gt_report(arguments$vcf_path,
                         arguments$merged_file)$snvs
# All-TV=noise
gt_load <- gt_load[!mut %in% gt_tv$mut]

noise_snvs <- noise_variants(arguments$vcf_path,
                                   arguments$caller,
                                   arguments$merged_file,
                                   gt_load,
                                   gt_tv)

out_noise_snvs = rbind(noise_snvs$tp, noise_snvs$fp, noise_snvs$fn)
out_noise_snvs$POS = as.numeric(out_noise_snvs$POS)
out_noise_snvs = out_noise_snvs[order(POS)]


fwrite(
    out_noise_snvs, paste0(arguments$working_directory, "/",
                           arguments$merged_file, "_",
                           arguments$caller, "_snvs_Noise.tsv"),
    row.names = FALSE, quote = FALSE, sep = "\t"
)

#write snvs stats in a txt file
stats <- data.frame(
    true_variants_recall = snvs$recall,
    noise_recall      = noise_snvs$noise_recall,
    noise_precision   = noise_snvs$noise_precision
)


write.table(stats,
            file = paste0(arguments$working_directory, "/",
                               arguments$merged_file, "_",
                               arguments$caller, "_snvs_stats.txt"),
            sep = "\t",
            row.names = FALSE,
            quote     = FALSE)


#INDELs------------------------------------------------------------------
print("Begin INDELs Variant Analysis")

pick_gt_stdz <- gt_stdz_indels(arguments$vcf_path,
                              arguments$merged_file)

indels <- call_indels(arguments$vcf_path,
                             arguments$caller,
                             arguments$merged_file,
                             pick_gt_stdz)



indels_out = rbind(indels$tp, indels$fp, indels$fn)
indels_out$POS = as.numeric(indels_out$POS)
indels_out = indels_out[order(POS)]


fwrite(
    indels_out, paste0(arguments$working_directory, "/",
                           arguments$merged_file, "_",
                           arguments$caller, "_indels.tsv"),
    row.names = FALSE, quote = FALSE, sep = "\t"
)

#write indel stats in a txt file
stats <- data.frame(
    indel_recall      = indels$indel_recall,
    indel_precision   = indels$indel_precision
)


write.table(stats,
            file = paste0(arguments$working_directory, "/",
                          arguments$merged_file, "_",
                          arguments$caller, "_indel_stats.txt"),
            sep = "\t",
            row.names = FALSE,
            quote     = FALSE)



print("Analysis Completed")
