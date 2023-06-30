


# rm(list = ls())
gc()

library(data.table)
library(stringr)

library(ggsci)

df = fread("Ground_truth_vs_Mutect2.clean.annotated.tsv")

df[which(df$`Mutect2 ALT` == "")]$`Mutect2 ALT` = NA


df = df[, c(
    "POS", 
    "Ground Truth AF",
    "Ground Truth ALT",
    "Mutect2 ALT",
    "Mutect2 AF"
), with = FALSE] |>
    unique()

gr1 = ggplot(data = df[, 1:3], aes(x = `Ground Truth AF`)) +
    
    geom_density(aes(color = `Ground Truth ALT`, fill = `Ground Truth ALT`),
                 alpha = .5) +
    
    scale_x_continuous(expand = c(0, 0), breaks = c(.25, .5, .75, 1), limits = c(0, 1), labels = scales::percent) +
    scale_y_continuous(expand = c(0, 0)) +
    
    scale_fill_npg() +
    scale_color_npg() +
    
    facet_wrap(vars(`Ground Truth ALT`), nrow = 1) +
    
    theme_minimal() +
    
    theme(
        legend.position = "none",
        
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 13),
        strip.text = element_text(face = "bold", size = 13),
        axis.text.y = element_text(face = "bold", size = 13),
        axis.text.x = element_text(face = "bold", size = 13),
        
        panel.spacing = unit(1, "lines"),
        
        panel.grid = element_line(linetype = "dashed")
    ) +
    
    labs(y = "Ground Truth (density)")


gr2 = ggplot(data = df[which(!is.na(`Mutect2 ALT`)), c(1, 4, 5)], aes(x = `Mutect2 AF`)) +
    
    geom_density(aes(color = `Mutect2 ALT`, fill = `Mutect2 ALT`),
                 alpha = .5) +
    
    scale_x_continuous(expand = c(0, 0), breaks = c(.25, .5, .75, 1), limits = c(0, 1), labels = scales::percent) +
    scale_y_continuous(expand = c(0, 0)) +
    
    scale_fill_npg() +
    scale_color_npg() +
    
    facet_wrap(vars(`Mutect2 ALT`), nrow = 1) +
    
    theme_minimal() +
    
    theme(
        legend.position = "none",
        axis.title.x = element_text(face = "bold", size = 13),
        axis.title.y = element_text(face = "bold", size = 13),
        strip.text = element_text(face = "bold", size = 13),
        axis.text.y = element_text(face = "bold", size = 13),
        axis.text.x = element_text(face = "bold", size = 13),
        
        panel.spacing = unit(1, "lines"),
        
        panel.grid = element_line(linetype = "dashed")
    ) +
    
    labs(
        x = "Allele Frequency",
        y = "Mutect2 (density)"
    )


library(patchwork)

multi = gr1 / gr2

ggsave(
    plot = multi, filename = "Plots/density-plot.pdf", device = cairo_pdf,
    width = 12, height = 10, units = "in"
)

ggsave(
    plot = multi, filename = "Plots/density-plot.svg",
    width = 12, height = 10, units = "in"
)

ggsave(
    plot = multi, filename = "Plots/density-plot.jpeg",
    width = 12, height = 10, units = "in", dpi = 600
)





