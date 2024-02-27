

rm(list = ls())
gc()



source("libraries.R")
source("helpers_VarScan.R")


folder = 'read_length/1000_100'

VarScan_somatic_vcf <- read.vcfR( paste0(folder, "/Merged_VarScan_norm.vcf"), verbose = FALSE )


v0 = gt_analysis(seq(1, 10), folder)

v1 = VarScan_somatic_vcf |>
    merge_VarScan(v0) |>
    clean_VarScan()

fwrite(
    v1, paste0(folder, "/Ground_truth_vs_VarScan.clean_norm.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
)
