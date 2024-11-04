

rm(list = ls())
gc()



# Libraries
library(data.table)
library(stringr)


# load data

tp = fread("GATK_indels_TP.tsv", sep = "\t")

fp = fread("GATK_indels_FP.tsv", sep = "\t")

fn = fread("GATK_indels_FN.tsv", sep = "\t")



tp = tp[, .(POS, REF, ALT, type)]

tp$REF_len <- str_length(tp$REF)

tp$ALT_len <- str_length(tp$ALT)

tp$len_dif <- tp$ALT_len - tp$REF_len




fp = fp[, .(POS, `Mutect2 REF`, `Mutect2 ALT`, type)]

fp$REF_len <- str_length(fp$`Mutect2 REF`)

fp$ALT_len <- str_length(fp$`Mutect2 ALT`)

fp$len_dif <- fp$ALT_len - fp$REF_len


colnames(fp) <- c("POS", "REF", "ALT", "type", "REF_len", "ALT_len", "len_dif")




fn = fn[, .(POS, `Ground Truth REF`, `Ground Truth ALT`, type)]

fn$REF_len <- str_length(fn$`Ground Truth REF`)

fn$ALT_len <- str_length(fn$`Ground Truth ALT`)


fn$len_dif <- fn$ALT_len - fn$REF_len


colnames(fn) <- c("POS", "REF", "ALT", "type", "REF_len", "ALT_len", "len_dif")


# Combine the datasets 

data = rbind(tp, fp)

df = rbind(data, fn)




# plot ---------

library(ggplot2)
library(dplyr)
library(colorspace)


# Adjust data so that each type has its own y-offset
df <- df |>
    mutate(y_cycle = case_when(
        type == "FN" ~ len_dif + 50,   # Shift FN cycle outward
        type == "FP" ~ len_dif + 25,   # Shift FP cycle to middle
        type == "TP" ~ len_dif         # Keep TP at the center
    ))




p = ggplot(df, aes(x = POS, y = y_cycle)) +
    
    # Lollipop segments: start each from the respective baseline to the point
    # geom_segment(
    #     aes(x = POS, xend = POS,
    #         y = ifelse(type == "FN", 50, ifelse(type == "FP", 25, 0)),
    #         yend = y_cycle),
    #     color = "grey75", linewidth = 0.25, lineend = "round"
    # ) +

    # add connecting line
    geom_line(aes(color = type), linewidth = 0.25) + 

    # Dashed lines for separation of each cycle level
    geom_hline(yintercept = 50, color = "grey40") +
    geom_hline(yintercept = 25, color = "grey40") +
    geom_hline(yintercept = 0,  color = "grey40") +
    
    # Add points at the end of each segment for the lollipop head
    geom_point(aes(fill = type, color = type), size = 1.5, shape = 21, stroke = .15) +
    
    
    
    # Define custom colors for each type
    scale_fill_manual(values = c("TP" = "#a78d95", "FP" = "#ae4364", "FN" = "#43ae8d")) +
    scale_color_manual(values = c("TP" = "#a78d95", "FP" = "#ae4364", "FN" = "#43ae8d") |> darken(.25), guide = "none") +
    
    scale_x_continuous(breaks = c(0, 4751, 9503, 14255, 19007), limits = c(0, 19007)) +
    
    
    coord_radial(start = pi / 2.5, inner.radius = .25, end = 2.6 * pi) +
    
    theme_minimal() +
    
    theme(
        axis.text.y = element_blank(),
        
        panel.grid.major = element_line(linewidth = 0.35),
        panel.grid.minor = element_blank(),
        
        plot.margin = margin(20, 20, 20, 20),
        
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)
        
    ) +
    
    labs(
        title = "Ground Truth vs Mutect2 INDELS",
        y = "REF vs ALT Length Difference",
        x = "Chromosomal Position",
        fill = "Type"
    )

p

ggsave(
    plot = p, filename = "Rplot_2.jpeg",
    width = 12, height = 12, units = "in", dpi = 600
)
