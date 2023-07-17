#'
#' This R script produces the barplots of the genomic content of the reference.
#'
#'
#'Input: fasta reference file
#'
#'Output: barplots of the genomic content of the reference
#'
#'
#'

rm(list = ls())
gc()


library(data.table)
library(stringr)
library(seqinr)
library(ggplot2)
library(stringr)



read_fasta = read.fasta(file="../reference/TP53.fasta", as.string = FALSE)

test = read_fasta$hg38_knownGene_ENST00000610292.4[c(1:19080)]

test = test |> str_to_upper()


test = test |> table()
test = test |> as.data.frame()

colnames(test) = c("Base", "Frequency")

test$Perc = test$Frequency / sum(test$Frequency)



gr = ggplot(data=test, aes(x = Base, y = Perc, fill = Base)) +
    geom_bar(stat="identity", position="dodge", alpha = .5) +
    scale_fill_npg() +
    scale_color_npg() +
    theme_minimal() +
    theme(
        legend.position = "none",
        axis.title.x = element_text(face = "bold", size = 13),
        strip.text = element_text(face = "bold", size = 13),
        axis.title.y = element_text(face = "bold", size = 13),
        panel.spacing = unit(1, "lines"),
        panel.grid = element_line(linetype = "dashed")
    ) +
    labs(
        x = "Base",
        y = "Percentage"
    )

ggsave(
    plot = gr, filename = "Plots/reference-plot.pdf", device = cairo_pdf,
    width = 12, height = 10, units = "in"
)
ggsave(
    plot = gr, filename = "Plots/reference-plot.svg",
    width = 12, height = 10, units = "in"
)
ggsave(
    plot = gr, filename = "Plots/reference-plot.jpeg",
    width = 12, height = 10, units = "in", dpi = 600
)