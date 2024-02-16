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
source("helpers_freebayes.R")


folder = 'read_length/1000_300'

df = fread(paste0(folder, "/Ground_truth_vs_freebayes.clean.tsv"))

vcf_read_GT <- read.vcfR(paste0(folder, "/Merged_ground_truth.vcf"), verbose = FALSE )

vcf_read_freebayes <- read.vcfR(paste0(folder, "/Merged.freebayes.vcf"), verbose = FALSE )

out1 = bar_plots_freebayes(df)
out2 = density_plot_freebayes(df)
out3 = bubble_plots_freebayes(df)
out4 = venn_plot_freebayes(vcf_read_GT, vcf_read_freebayes)

library(patchwork)

multi2 = out2$groundtruth / out2$freebayes &
    
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
#    plot = multi, filename = paste0(folder, "/Plots/Poster_freebayes.svg"),
#    width = 16, height = 12, units = "in", dpi = 600
#)

ggsave(
    plot = multi, filename = paste0(folder, "/Plots/Poster_freebayes.png"),
    width = 16, height = 12, units = "in", dpi = 600
)

ggsave(
    plot = out4, filename = paste0(folder, "/Plots/Venn_freebayes.png"),
    width = 8, height = 8, units = "in", dpi = 600
)




    
    