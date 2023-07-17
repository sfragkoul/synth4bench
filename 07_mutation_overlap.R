

# rm(list = ls())
gc()

library(data.table)
library(stringr)
library(vcfR)
library(dplyr)

vcf_read_GT <- read.vcfR("Merged_Ground_Truth.vcf", 
                      verbose = FALSE )

vcf_GT = vcfR::getFIX(vcf_read_GT) |> as.data.frame() |> setDT()
vcf_GT$scenario = "GT"
vcf_GT[which(vcf_GT$POS == unique(vcf_GT$POS))]

vcf_GT = distinct(vcf_GT)

vcf_read_gatk <- read.vcfR("Merged2_GATK.vcf", 
                               verbose = FALSE )

vcf_gatk = vcfR::getFIX(vcf_read_gatk) |> as.data.frame() |> setDT()
vcf_gatk$scenario = "GATK"


rm(vcf_read_GT, vcf_read_gatk)
x = rbind(vcf_GT, vcf_gatk)

rm(vcf_GT, vcf_gatk)


y = x[, c("CHROM", "POS", "REF", "ALT", "scenario"), with = FALSE]



y$mut = paste(y$CHROM, y$POS, y$REF, y$ALT, sep = ":")


z = y[, by = mut, .(
    scenarios = paste(scenario, collapse = "+")
)]

z = z[order(z$scenarios), ]

fwrite(
    z, "All_mutation-overlaps.csv",
    row.names = FALSE, quote = TRUE, sep = ","
)


library(ggvenn)
library(ggplot2)

y = split(y, y$scenario)


y = list(
    'Ground Truth' = y$GT$mut,
    'GATK' = y$GATK$mut
)

gr = ggvenn(
    y,
    fill_color = c("#43ae8d", "#ae4364")
) +
    
    coord_equal(clip = "off")

gr6 = gr

ggsave(
    plot = gr, filename = "Plots/All-overlap-plot.pdf", device = cairo_pdf,
    width = 8, height = 8, units = "in"
)

ggsave(
    plot = gr, filename = "Plots/All-overlap-plot.svg",
    width = 8, height = 8, units = "in"
)

ggsave(
    plot = gr, filename = "Plots/All-overlap-plot.jpeg",
    width = 8, height = 8, units = "in", dpi = 600
)
