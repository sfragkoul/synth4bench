#'
#'This R script produces the final Multipanel Figure of the Benchmarking 
#'process.
#'
#'Input: comparison tsv file, ground truth vcf, caller's vcf
#'
#'Output: final Multipanel Figure
#'
#'#'Authors: Nikos Pechlivanis(github:npechl), Stella Fragkouli(github:sfragkoul)
#'

rm(list = ls())
gc()



source("libraries.R")
source("helpers_VarDict.R")


folder = 'read_length/1000_200'

df = fread(paste0(folder, "/Ground_truth_vs_VarDict.clean_norm.tsv"))

vcf_read_GT <- read.vcfR(paste0(folder, "/Merged_ground_truth_norm.vcf"), verbose = FALSE )

vcf_read_VarDict <- read.vcfR(paste0(folder, "/Merged_VarDict_norm.vcf"), verbose = FALSE )

out1 = bar_plots_VarDict(df)
out2 = density_plot_VarDict(df)
out3 = bubble_plots_VarDict(df)
out4 = venn_plot_VarDict(vcf_read_GT, vcf_read_VarDict)

library(patchwork)

multi2 = out2$groundtruth / out2$VarDict &
    
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
#    plot = multi, filename = paste0(folder, "/Plots/Poster_VarDict.svg"),
#    width = 16, height = 12, units = "in", dpi = 600
#)

ggsave(
    plot = multi, filename = paste0(folder, "/Plots/Poster_VarDict_norm.png"),
    width = 16, height = 12, units = "in", dpi = 600
)

ggsave(
    plot = out4, filename = paste0(folder, "/Plots/Venn_VarDict_norm.png"),
    width = 8, height = 8, units = "in", dpi = 600
)




    
    