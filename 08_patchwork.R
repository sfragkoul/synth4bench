#'
#'This This R script produces the final Figure of the Benchmarking of the 
#'poster.
#'Input: all produced plots
#'
#'Output: final Figure for the poster
#'



library(patchwork)


source("03_bar_plots.R")
source("03_density_plot.R")
source("03_bubble_plots.R")
source("04_mutation_overlap.R")

multi2 = gr1 / gr2 &
    
    theme(
        plot.margin = margin(10, 10, 10, 10)
    )




ann1 = (gr3 + theme(plot.margin = margin(r = 50))) + 
    (gr4 + theme(plot.margin = margin(r = 50))) + 
    multi2 +
    
    plot_layout(
        widths = c(1, 1, 3)
    )



ann2 = gr5 + gr6 +
    
    plot_layout(
        widths = c(2, 1)
    )


multi = ann1 / ann2 +
    
    plot_layout(heights = c(1.5, 1))


ggsave(
    plot = multi, filename = "Plots/Poster.svg",
    width = 16, height = 12, units = "in", dpi = 600
)

ggsave(
    plot = multi, filename = "Plots/Poster.jpeg",
    width = 16, height = 12, units = "in", dpi = 600
)






    
    