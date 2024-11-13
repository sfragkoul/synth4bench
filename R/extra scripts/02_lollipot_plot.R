
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


# plot ---------

df$type <- factor(df$type, levels = c("FN", "TP", "FP"))


fn_count <- nrow(fn)
tp_count <- nrow(tp)
fp_count <- nrow(fp)


facet_labels <- c(
    FN = paste0("FN (n=", fn_count, ")"),
    TP = paste0("TP (n=", tp_count, ")"),
    FP = paste0("FP (n=", fp_count, ")")
)


p=ggplot(df, aes(x = POS, y = len_dif)) +
    
    geom_segment(aes(x = POS, xend = POS, y = 0, yend = len_dif), color = "grey92") +
    
    geom_hline(yintercept = 0) +
    
    geom_point(aes(color = type), size = 1.5) +
    
    scale_color_manual(
        values = c("TP" = "#a78d95", "FP" = "#ae4364", "FN" = "#43ae8d"),
        labels = c("FN", "TP", "FP")
    ) +
    
    
    facet_wrap(vars(type), labeller = as_labeller(facet_labels)) +  # Set facet titles with custom order
    
    scale_x_continuous(breaks = c(0, 10000, 19010), limits = c(0, 19010)) +
    
    theme_minimal() +
    
    theme(
        panel.grid.major = element_line(linewidth = 0.25),
        panel.grid.minor = element_blank(),
        plot.margin = margin(20, 20, 20, 20),
        panel.border = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)
    ) +
    
    labs(
        title = "Ground Truth vs Mutect2 INDELS",
        y = "REF vs ALT length",
        x = "Chromosomal Position", 
        color = "Type"
    )




ggsave(
    plot = p, filename = "Mutect2_indels.jpeg",
    width = 16, height = 12, units = "in", dpi = 600
)


