
# coverage plot ----------------------------------------

# rm(list = ls())
gc()

library(data.table)
library(stringr)

df = fread("Ground_truth_vs_Mutect2.clean.annotated.tsv")

df[which(df$`Mutect2 ALT` == "")]$`Mutect2 ALT` = NA


df = df[, c(
    "POS", 
    "Ground Truth DP",
    "Mutect2 DP"
), with = FALSE] |>
    unique() |>
    
    melt(id.vars = "POS", variable.factor = FALSE, value.factor = FALSE)

library(ggplot2)
library(ggforce)

gr = ggplot(data = df) +
    
    geom_point(aes(x = variable, y = value, fill = variable),
               position = position_jitternormal(sd_x = .01, sd_y = 0),
               shape = 21, stroke = .1, size = 2.5) +
    
    geom_boxplot(aes(x = variable, y = value, fill = variable),
                 width = .25, alpha = .5, outlier.shape = NA) +
    
    scale_fill_manual(
        values = c(
            "Ground Truth DP" = "#43ae8d",
            "Mutect2 DP"      = "#ae4364"
        )
    ) +
    
    scale_x_discrete(
        breaks = c("Ground Truth DP", "Mutect2 DP"),
        labels = c("Ground Truth", "Mutect2")
    ) +
    
    scale_y_continuous(labels = scales::comma) +
    
    theme_minimal() +
    
    theme(
        legend.position = "none",
        
        axis.title.x = element_blank(),
        axis.text.x = element_text(face = "bold", size = 13),
        axis.text.y = element_text(face = "bold", size = 13),
        
        axis.line = element_line(),
        axis.ticks = element_line(),
        
        panel.grid = element_blank(),
        
        plot.margin = margin(20, 20, 20, 20)
    ) +
    
    labs(
        y = "Coverage (No. of reads)"
    )


ggsave(
    plot = gr, filename = "Plots/coverage-plot.pdf", device = cairo_pdf,
    width = 8, height = 8, units = "in"
)

ggsave(
    plot = gr, filename = "Plots/coverage-plot.svg",
    width = 8, height = 8, units = "in"
)

ggsave(
    plot = gr, filename = "Plots/coverage-plot.jpeg",
    width = 8, height = 8, units = "in", dpi = 600
)    

gr3 = gr

# allele frequency plot ----------------------------------------

# rm(list = ls())
gc()

library(data.table)
library(stringr)

df = fread("Ground_truth_vs_Mutect2.clean.annotated.tsv")

df[which(df$`Mutect2 ALT` == "")]$`Mutect2 ALT` = NA


df = df[, c(
    "POS", 
    "Ground Truth AF",
    "Mutect2 AF"
), with = FALSE] |>
    unique() |>
    
    melt(id.vars = "POS", variable.factor = FALSE, value.factor = FALSE)

library(ggplot2)
library(ggforce)

gr = ggplot(data = df[which(!is.na(value) & value != 0)]) +
    
    geom_point(aes(x = variable, y = value, fill = variable),
               position = position_jitternormal(sd_x = .01, sd_y = 0),
               shape = 21, stroke = .1, size = 2.5) +
    
    geom_boxplot(aes(x = variable, y = value, fill = variable),
                 width = .25, alpha = .5, outlier.shape = NA) +
    
    scale_fill_manual(
        values = c(
            "Ground Truth AF" = "#43ae8d",
            "Mutect2 AF"      = "#ae4364"
        )
    ) +
    
    scale_x_discrete(
        breaks = c("Ground Truth AF", "Mutect2 AF"),
        labels = c("Ground Truth", "Mutect2")
    ) +
    
    scale_y_continuous(labels = scales::percent, trans = "log10") +
    
    theme_minimal() +
    
    theme(
        legend.position = "none",
        
        axis.title.x = element_blank(),
        axis.text.x = element_text(face = "bold", size = 13),
        axis.text.y = element_text(face = "bold", size = 13),
        
        axis.line = element_line(),
        axis.ticks = element_line(),
        
        panel.grid = element_blank(),
        
        plot.margin = margin(20, 20, 20, 20)
    ) +
    
    labs(
        y = "Allele Frequency"
    )


ggsave(
    plot = gr, filename = "Plots/AF-plot.pdf", device = cairo_pdf,
    width = 8, height = 8, units = "in"
)

ggsave(
    plot = gr, filename = "Plots/AF-plot.svg",
    width = 8, height = 8, units = "in"
)

ggsave(
    plot = gr, filename = "Plots/AF-plot.jpeg",
    width = 8, height = 8, units = "in", dpi = 600
)

gr4 = gr
