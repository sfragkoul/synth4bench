

rm(list = ls())
gc()



source("libraries.R")
source("analysis_helpers_gatk.R")


folder = 'read_length/1000_200'

gatk_somatic_vcf <- read.vcfR( paste0(folder, "/Merged_GATK.vcf"), verbose = FALSE )


v0 = gt_analysis(seq(1, 10), folder)

v1 = gatk_somatic_vcf |>
    merge_gatk(v0) |>
    clean_gatk()
fwrite(
    v1, paste0(folder, "/Ground_truth_vs_Mutect2.clean.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
)
