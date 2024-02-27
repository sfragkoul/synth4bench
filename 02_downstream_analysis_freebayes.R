#'
#'This R script produces the final tsv file with the comparison for each caller
#' against the ground truth.
#'
#'Input: ground truth vcf, caller's vcf (from function)
#'
#'Output: tsv file with the comparison
#'

rm(list = ls())
gc()



source("libraries.R")
source("helpers_freebayes.R")


folder = 'read_length/1000_300'

freebayes_somatic_vcf <- read.vcfR( paste0(folder, "/Merged_freebayes_norm.vcf"), verbose = FALSE )


v0 = gt_analysis(seq(1, 10), folder)

v1 = freebayes_somatic_vcf |>
    merge_freebayes(v0) |>
    clean_freebayes()

fwrite(
    v1, paste0(folder, "/Ground_truth_vs_freebayes.clean_norm.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
)
