

rm(list = ls())
gc()



source("libraries.R")
source("helpers_VarDict.R")


folder = 'read_length/1000_200'

VarDIct_somatic_vcf <- read.vcfR( paste0(folder, "/Merged_VarDIct.vcf"), verbose = FALSE )


v0 = gt_analysis(seq(1, 10), folder)

v1 = VarDIct_somatic_vcf |>
    merge_VarDict(v0) |>
    clean_VarDict()

fwrite(
    v1, paste0(folder, "/Ground_truth_vs_VarDict.clean.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
)
