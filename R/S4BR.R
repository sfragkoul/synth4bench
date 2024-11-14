#'A script, written in R, that calls the appropriate functions 
#'to perform the comparison between the ground truth and the caller.
#'
#' Input files: ground truth files, caller vcf file
#'
#' Output files:  tsv files with the comparison between the ground truth 
#' and the caller.
#'
#' Authors: Nikos Pechlivanis(github:npechl), Stella Fragkouli(github:sfragkoul)
#' 


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
  make_option(c("-v", "--vcf_path"), 
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

#SNVS TP-----------------------------------------------------------------------
print("Begin SNVs TP Variant Analysis")

gt <- gt_analysis(seq_len(arguments$runs),
                  arguments$working_directory,
                  arguments$merged_file)

out_df_snvs_tp <- read_vcf_snvs_TP(arguments$vcf_path,
                                   arguments$caller,
                                   gt,
                                   arguments$merged_file)

fwrite(
  out_df_snvs_tp, paste0(arguments$working_directory,
                         "/",
                         arguments$merged_file,
                         "_",
                         arguments$caller, "_snvs_TP.tsv"),

  row.names = FALSE, quote = FALSE, sep = "\t"
)

#SNVS FP & FN------------------------------------------------------------------
print("Begin SNVs FP Variant Analysis")

gt_all = load_gt_report(arguments$vcf_path,
                        arguments$merged_file)$all

gt_snvs = load_gt_report(arguments$vcf_path,
                         arguments$merged_file)$snvs

pick_gt = load_gt_vcf(arguments$vcf_path,
                      arguments$merged_file)

out_df_snvs_fp <- read_vcf_snvs_FP(arguments$vcf_path,
                                   arguments$caller,
                                   arguments$merged_file,
                                   pick_gt,
                                   gt_all)

fwrite(
    out_df_snvs_fp, paste0(arguments$working_directory, "/",
                           arguments$merged_file, "_",
                           arguments$caller, "_snvs_FP.tsv"),
    row.names = FALSE, quote = FALSE, sep = "\t"
)

print("Begin SNVs FN Variant Analysis")
out_df_snvs_fn = read_vcf_snvs_FN(arguments$vcf_path,
                                  arguments$caller,
                                  arguments$merged_file, 
                                  pick_gt)

fwrite(
    out_df_snvs_fn, paste0(arguments$working_directory, "/", 
                               arguments$merged_file, "_", 
                               arguments$caller, "_snvs_FN.tsv"),
    row.names = FALSE, quote = FALSE, sep = "\t"
)

#INDELs------------------------------------------------------------------
print("Begin INDELs Variant Analysis")

pick_gt_stdz = gt_stdz_indels(arguments$vcf_path,
                              arguments$merged_file)

print("Begin TP INDELs Variant Analysis")
tp_indels = call_tp_indels(arguments$vcf_path,
                           arguments$caller,
                           arguments$merged_file,
                           pick_gt_stdz)
fwrite(
        tp_indels, paste0(arguments$working_directory, "/",
                         arguments$merged_file, "_",
                         arguments$caller, "_indels_TP.tsv"),
  row.names = FALSE, quote = FALSE, sep = "\t"
)

print("Begin FN INDELs Variant Analysis")
fn_indels = call_fn_indels(arguments$vcf_path,
                           arguments$caller,
                           arguments$merged_file,
                           pick_gt_stdz)
fwrite(
  fn_indels, paste0(arguments$working_directory, "/",
                    arguments$merged_file, "_",
                    arguments$caller, "_indels_FN.tsv"),
  row.names = FALSE, quote = FALSE, sep = "\t"
)

print("Begin FP INDELs Variant Analysis")
fp_indels = call_fp_indels(arguments$vcf_path,
                           arguments$caller,
                           arguments$merged_file,
                           pick_gt_stdz)
fwrite(
  fp_indels, paste0(arguments$working_directory, "/",
                    arguments$merged_file, "_",
                    arguments$caller, "_indels_FP.tsv"),
  row.names = FALSE, quote = FALSE, sep = "\t"
)











