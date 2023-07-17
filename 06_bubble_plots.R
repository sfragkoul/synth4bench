


# rm(list = ls())
gc()

library(data.table)
library(stringr)

df = fread("Ground_truth_vs_Mutect2.clean.annotated.tsv")

df[which(df$`Mutect2 ALT` == "")]$`Mutect2 ALT` = NA

library(ggplot2)
library(ggforce)

gr = ggplot(data = df) +
    
    geom_point(
        aes(x = POS, y = `Ground Truth ALT`, size = `Ground Truth AF`, fill = "Ground Truth"),
        position = position_nudge(y = .15), shape = 21, stroke = .25
    ) +
    
    geom_point(
        aes(x = POS, y = `Mutect2 ALT`, size = `Mutect2 AF`, fill = "Mutect2"),
        position = position_nudge(y = -.15), shape = 21, stroke = .25
    ) +
    
    scale_size_continuous(
        range = c(2, 10), 
        limits = c(0, 1),
        breaks = c(.05, .1, .2, .5, .8),
        labels = scales::percent,
        guide = guide_legend(
            title = "Allele Frequency",
            title.position = "top"
        )
    ) +
    
    scale_fill_manual(
        values = c(
            "Ground Truth" = "#43ae8d",
            "Mutect2"      = "#ae4364"
        ),
        
        guide = guide_legend(
            title = "Category",
            title.position = "top",
            
            override.aes = list(size = 3.5)
        )
    ) +
    
    scale_x_continuous(labels = scales::comma_format(suffix = " bp"), 
                       expand = c(0, 0),
                       limits = c(0, 20000)) +
    
    theme_minimal() +
    
    theme(
        legend.position = "bottom",
        
        axis.line.x = element_line(),
        axis.ticks.x = element_line(),
        
        axis.text.y = element_text(face = "bold", size = 13),
        axis.text.x = element_text(face = "bold", size = 11),
        axis.title.y = element_text(face = "bold", size = 13),
        axis.title.x = element_text(face = "bold", size = 13),
        
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        
        plot.margin = margin(20, 30, 20, 20)
    ) +
    
    labs(
        x = "Chromosomal Position",
        y = "Alterations"
    )



ggsave(
    plot = gr, filename = "Plots/Bubble-plot.pdf", device = cairo_pdf,
    width = 12, height = 6, units = "in"
)

ggsave(
    plot = gr, filename = "Plots/Bubble-plot.svg",
    width = 12, height = 6, units = "in"
)

ggsave(
    plot = gr, filename = "Plots/Bubble-plot.jpeg",
    width = 12, height = 6, units = "in", dpi = 600
)


gr5 = gr
