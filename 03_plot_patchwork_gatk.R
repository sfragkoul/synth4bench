#'
#'This This R script produces the final Figure of the Benchmarking of the 
#'poster.
#'Input: all produced plots
#'
#'Output: final Figure for the poster
#'

rm(list = ls())
gc()



source("libraries.R")
source("helpers_gatk.R")


folder = 'read_length/1000_200'

df = fread(paste0(folder, "/Ground_truth_vs_Mutect2.clean.tsv"))

vcf_read_GT <- read.vcfR(paste0(folder, "/Merged_ground_truth.vcf"), verbose = FALSE )

vcf_read_gatk <- read.vcfR(paste0(folder, "/Merged_GATK.vcf"), verbose = FALSE )

out1 = bar_plots_gatk(df)
out2 = density_plot_gatk(df)
out3 = bubble_plots_gatk(df)
out4 = venn_plot_gatk(vcf_read_GT, vcf_read_gatk)

library(patchwork)

multi2 = out2$groundtruth / out2$mutect2 &
    
    theme(
        plot.margin = margin(10, 10, 10, 10)
    )




ann1 = (out1$coverage + theme(plot.margin = margin(r = 50))) + 
    (out1$allele + theme(plot.margin = margin(r = 50))) + 
    multi2 +
    
    plot_layout(
        widths = c(1, 1, 3)
    )



ann2 = out3 + out4 +
    
    plot_layout(
        widths = c(2, 1)
    )


multi = ann1 / ann2 +
    
    plot_layout(heights = c(1.5, 1)) + 
    plot_annotation(title = folder)


#ggsave(
#    plot = multi, filename = paste0(folder, "/Plots/Poster_gatk.svg"),
#    width = 16, height = 12, units = "in", dpi = 600
#)

ggsave(
    plot = multi, filename = paste0(folder, "/Plots/Poster_gatk.png"),
    width = 16, height = 12, units = "in", dpi = 600
)

ggsave(
    plot = out4, filename = paste0(folder, "/Plots/Venn_gatk.png"),
    width = 8, height = 8, units = "in", dpi = 600
)




    
    