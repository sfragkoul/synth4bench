

rm(list = ls())
gc()



# Libraries
library(data.table)
library(stringr)
library(ggplot2)


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



df$condition <- ifelse(df$len_dif > 0, "insertion", "deletion")


df$type <- factor(df$type, levels = c("FN", "TP", "FP"))




# plot ---------
p = ggplot(df, aes(x = POS, y = len_dif)) +
    
    # geom_segment(aes(x = POS, xend = POS, y = 0, yend = len_dif), color = "grey85") +
    
    # geom_point(color= "#b9b8e7", size = 3.5) +
    
    geom_hline(yintercept = 0) +
    geom_point(aes(color =type), size = 2) +
    
    scale_color_manual(values = c("TP" = "#b9b8e7", "FP" = "#0072B5", "FN" = "#BC3C29")) +  # Customize colors here
  
    # geom_density_2d(aes(linetype = condition)) +
    
    # scale_y_continuous(limits = c(-10.5, 10.5)) +
    
    # facet_grid(rows = vars(condition), scales = "free_y") +
    
    facet_wrap(vars(type)) +
    
    theme_minimal() +
    
    theme(
        panel.grid.major = element_line(linewidth = .25),
        panel.grid.minor = element_blank(),
        
        plot.margin = margin(20, 20, 20, 20),
        panel.border = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)
    ) +
    
    labs(
        title = "Ground truth vs Mutect2 INDELS",
        y = "REF vs ALT length"
    )


# ggsave(
#     plot = p, filename = "Rplot.png",
#     width = 14, height = 12, units = "in", dpi = 600
# )  
