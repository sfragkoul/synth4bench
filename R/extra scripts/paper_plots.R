#'
#'This R script produces the final Multipanel Figure for the paper.
#'
#'Input: comparison tsv file, ground truth vcf, caller's vcf
#'
#'Output: Paper Multipanel Figure
#'
#'Authors: Nikos Pechlivanis(github:npechl), Stella Fragkouli(github:sfragkoul)
#'

source("R/libraries.R")

#SNVs TP-----------------------------------------------------------------------
folder = 'read_length/1000_300'
name = "1000_300"

df = list()

df[["VarDict"]] = fread(paste0(folder, "/Ground_truth_vs_VarDict.clean_norm.tsv"))
df[["VarScan"]] = fread(paste0(folder, "/Ground_truth_vs_VarScan.clean_norm.tsv"))
df[["Freebayes"]] = fread(paste0(folder, "/Ground_truth_vs_freebayes.clean_norm.tsv"))
df[["Mutect2"]] = fread(paste0(folder, "/Ground_truth_vs_Mutect2.clean_norm.tsv"))
df[["LoFreq"]] = fread(paste0(folder, "/Ground_truth_vs_LoFreq.clean_norm.tsv"))

df$Freebayes$`Freebayes AO` = NULL

# df |> lapply(colnames)

df = df |> lapply(function(x) {
    
    colnames(x) = c(
        "POS", 
        "Ground Truth REF", "Ground Truth ALT", "Ground Truth DP", "Ground Truth AF",
        "Caller REF", "Caller ALT", "Caller DP", "Caller AF"
    )
    
    return(x)
    
}) |> rbindlist(idcol = "caller")


df[which(df$`Caller ALT` == "")]$`Caller ALT` = NA



# DP plot
s1 = df[, c(
    "caller",
    "POS", 
    "Caller DP"
), with = FALSE] 

s2 = data.table(
    "caller" = "Ground\nTruth",
    "POS" = df$POS,
    "Caller DP" = df$`Ground Truth DP`
) 

q = rbind(s1, s2) |> 
    unique() |>
    melt(id.vars = c("caller", "POS"), variable.factor = FALSE, value.factor = FALSE)

q$caller = q$caller |> factor(levels = c("Ground\nTruth", "Freebayes", "LoFreq", "Mutect2", "VarDict", "VarScan"))

gr1 = ggplot(data = q) +
    
    geom_point(aes(x = caller, y = value, fill = caller),
               position = position_jitternormal(sd_x = .025, sd_y = 0),
               shape = 21, stroke = .1, size = 1) +
    
    geom_boxplot(aes(x = caller, y = value, fill = caller),
                 width = .3, alpha = .5, outlier.shape = NA) +
    
    geom_hline(yintercept = mean(q[which(caller == "Ground\nTruth")]$value), linewidth = .5, color = "yellow2") +
    
    scale_fill_manual(
        values = c(
            "VarDict" = "#8d43ae",
            "VarScan" = "#439aae",
            "Freebayes" = "#ae8d43",
            "Mutect2" = "#ae4364",
            "LoFreq" = "#c974ba",
            "Ground\nTruth" = "#43ae8d"
        )
    ) +
    
    # scale_x_discrete(
    #     breaks = c("Ground Truth DP", "Mutect2 DP"),
    #     labels = c("Ground Truth", "Mutect2")
    # ) +
    
    scale_y_continuous(labels = scales::comma) +
    
    theme_minimal() +
    
    theme(
        legend.position = "none",
        
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 13),
        axis.text.x = element_text(face = "bold", size = 12),
        axis.text.y = element_text(face = "bold", size = 12),
        
        axis.line = element_line(),
        axis.ticks = element_line(),
        
        panel.grid = element_blank()
    ) +
    
    labs(
        y = "Coverage (No. of reads)"
    )



# AF plot

s1 = df[, c(
    "caller",
    "POS", 
    "Caller AF"
), with = FALSE] 

s2 = data.table(
    "caller" = "Ground\nTruth",
    "POS" = df$POS,
    "Caller AF" = df$`Ground Truth AF`
) 

q = rbind(s1, s2) |> 
    unique() |>
    melt(id.vars = c("caller", "POS"), variable.factor = FALSE, value.factor = FALSE)

q$caller = q$caller |> factor(levels = c("Ground\nTruth", 
                                         "Freebayes", 
                                         "LoFreq", 
                                         "Mutect2", 
                                         "VarDict", 
                                         "VarScan"))

gr2 = ggplot(data = q) +
    
    geom_point(aes(x = caller, y = value, fill = caller),
               position = position_jitternormal(sd_x = .025, sd_y = 0),
               shape = 21, stroke = .1, size = 1) +
    
    geom_boxplot(aes(x = caller, y = value, fill = caller),
                 width = .3, alpha = .5, outlier.shape = NA) +
    
    geom_hline(yintercept = mean(q[which(caller == "Ground\nTruth")]$value), linewidth = .5, color = "yellow2") +
    
    scale_fill_manual(
        values = c(
            "VarDict" = "#8d43ae",
            "VarScan" = "#439aae",
            "Freebayes" = "#ae8d43",
            "Mutect2" = "#ae4364",
            "LoFreq" = "#c974ba",
            "Ground\nTruth" = "#43ae8d"
        )
    ) +
    
    
    scale_y_continuous(labels = scales::percent, trans = "log10", limits = c(.00001, 1)) +

    
    theme_minimal() +
    
    theme(
        legend.position = "none",
        
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 13),
        axis.text.x = element_text(face = "bold", size = 12),
        axis.text.y = element_text(face = "bold", size = 12),
        
        axis.line = element_line(),
        axis.ticks = element_line(),
        
        panel.grid = element_blank()
    ) +
    
    labs(
        y = "Allele Frequency"
    )

# plot 3 

#vcfs = list()
#
#GroundTruth_vcf = read.vcfR(paste0(folder, "/Merged_ground_truth_norm.vcf"), verbose = FALSE)
#
#vcfs[["LoFreq"]] <- read.vcfR( paste0(folder, "/Merged_LoFreq_norm.vcf"), verbose = FALSE )
#vcfs[["Mutect2"]] <- read.vcfR( paste0(folder, "/Merged_GATK_norm.vcf"), verbose = FALSE )
#vcfs[["Freebayes"]] <- read.vcfR( paste0(folder, "/Merged_freebayes_norm.vcf"), verbose = FALSE )
#vcfs[["VarDict"]] <- read.vcfR( paste0(folder, "/Merged_VarDict_norm.vcf"), verbose = FALSE )
##vcfs[["VarScan"]] <- read.vcfR( paste0(folder, "/Merged_VarScan_norm.vcf"), verbose = FALSE )
#vcfs[["VarScan"]] <- read.vcfR( paste0(folder, "/VarScan_norm.vcf"), verbose = FALSE )
#
#
#q = vcfs |> 
#    lapply(function(p, q) {
#        
#        vcf_GT = vcfR::getFIX(q) |> as.data.frame() |> setDT()
#        vcf_GT$scenario = "GT"
#        
#        vcf_caller = vcfR::getFIX(p) |> as.data.frame() |> setDT()
#        vcf_caller$scenario = "caller"
#        
#        x = rbind(vcf_GT, vcf_caller)
#        y = x[, c("CHROM", "POS", "REF", "ALT", "scenario"), with = FALSE]
#        
#        y$mut = paste(y$CHROM, y$POS, y$REF, y$ALT, sep = ":")
#        
#        y = split(y, y$scenario)
#        
#        y = data.table(
#            "N" = which(y$GT$mut %in% y$caller$mut) |> length(),
#            "Freq" = (which(y$GT$mut %in% y$caller$mut) |> length()) / (c(y$GT$mut, y$caller$mut) |> unique() |> length())
#        )
#        
#        return(y)
#        
#    }, GroundTruth_vcf) |>
#    
#    rbindlist(idcol = "caller")
#
#
#q$gt = "Ground\nTruth"
#
#q$label = paste0(
#    q$N |> scales::comma(), "\n",
#    round(100 * q$Freq, digits = 2), "%"
#)
#
#q$caller = q$caller |> factor(levels = c("Freebayes", "LoFreq", "Mutect2", "VarDict", "VarScan"))
#
#library(shadowtext)
#
#gr3 = q |>
#    ggplot() +
#    geom_raster(aes(caller, gt, fill = Freq), hjust = .5, vjust = .5) +
#    # geom_point(aes(caller, gt, fill = Freq, size = Freq), shape = 21, color = "grey10") +
#    geom_shadowtext(aes(caller, gt, label = label), bg.color = "white", color = "grey10", size = 4, fontface = "bold") +
#    
#    scale_fill_gradient(low = "#EAC0BD", high = "#C11317") +
#   
#    scale_size_continuous(range = c(10, 50)) +
#    
#    
#    theme_minimal() +
#    theme(
#        legend.position = "none",
#        
#        panel.grid = element_blank(),
#        axis.title = element_blank(),
#        axis.text.x = element_text(face = "bold", size = 13),
#        axis.text.y = element_text(face = "bold", size = 13)
#    )


# DAF density

q = df[, c("caller", "POS", "Ground Truth AF", "Caller AF"), with = FALSE] |>
    unique()

q$delta = q$`Caller AF` - q$`Ground Truth AF`

q = q[, by = caller, cut := mean(delta, na.rm = TRUE)]


library(ggridges)
library(ggplot2)
library(ggtext)
library(colorspace)

q$caller = q$caller |>
    factor(levels = c("VarScan", "VarDict", "Mutect2", "LoFreq", "Freebayes"))

gr4 = q[which(!is.na(delta))] |>
    
    ggplot(aes(delta, caller)) +
    
    geom_density_ridges(aes(fill = cut), color = "grey96", linewidth = .15) +
    
    geom_vline(xintercept = 0, linewidth = .5, color = "yellow2") +
    
    scale_x_continuous(limits = c(-.5, .5), expand = c(0, 0)) +
    # xlim(-0.6, 0.8) +
    
    scale_fill_gradient2(
        low = "#2E5A87" |> darken(.5) |> alpha(.85), 
        mid = "#C3C8D2" |> alpha(.85), 
        high = "#A90C38" |> alpha(.85), 
        midpoint = -.0016
    ) +
    
    theme_minimal() +
    
    theme(
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x = element_markdown(),
        axis.text.y = element_text(face = "bold", size = 13),
        axis.text.x = element_text(size = 11)
    ) +
    
    labs(x = "**Δ***AF*")

# patchwork 

#gr = ( (gr1 | gr2) / gr4 / gr3 ) + 
gr = ( (gr1 | gr2) / gr4  ) + 
    plot_layout(heights = c(3, 2)) +
    plot_annotation(tag_levels = "A") &
    theme(
        plot.tag = element_text(face = "bold"),
        plot.margin = margin(5, 10, 10, 5)
    )

ggsave(
    plot = gr, filename = paste0(folder,"/Plots/Final_", name,".jpeg"),
    width = 14, height = 12, units = "in", dpi = 600
)


#SNVs FP & FN -----------------------------------------------------------------
#FP
folder = 'results/'
name = "Merged"

mutect = fread(paste0(folder, "Merged_Mutect2_snvs_FP.tsv"))
mutect = mutect[,c("CHROM",	"POS",	"ID",	"Mutect2 REF",	"Mutect2 ALT", 
                   "Mutect2 QUAL",	"Mutect2 FILTER",	"Mutect2 DP",	"Mutect2 AF",
                   "mut",	"Ground Truth DP",	"DP Percentage",	"type")]

df = list()

df[["VarDict"]] = fread(paste0(folder, "Merged_VarDict_snvs_FP.tsv"))
df[["VarScan"]] = fread(paste0(folder, "Merged_Varscan_snvs_FP.tsv"))
df[["Freebayes"]] = fread(paste0(folder, "Merged_Freebayes_snvs_FP.tsv"))
df[["Mutect2"]] = mutect
df[["LoFreq"]] = fread(paste0(folder, "Merged_LoFreq_snvs_FP.tsv"))


df = df |> lapply(function(x) {
    
    colnames(x) = c(
        "CHROM",	"POS",	"ID",	"Caller REF",	"Caller ALT", 
        "Caller QUAL",	"Caller FILTER",	"Caller DP",	"Caller AF",
        "mut",	"Ground Truth DP",	"DP Percentage",	"type"
    )
    
    return(x)
    
}) |> rbindlist(idcol = "caller")


df[which(df$`Caller ALT` == "0")]$`Caller ALT` = NA

#DP
s1 = df[, c(
    "caller",
    "POS", 
    "Caller DP"
), with = FALSE] 


q = rbind(s1) |> 
    unique() |>
    melt(id.vars = c("caller", "POS"), variable.factor = FALSE, value.factor = FALSE)

q$caller = q$caller |> factor(levels = c("Freebayes", "LoFreq", "Mutect2", "VarDict", "VarScan"))

fp1 = ggplot(data = q) +
    
    geom_point(aes(x = caller, y = value, fill = caller),
               position = position_jitternormal(sd_x = .025, sd_y = 0),
               shape = 21, stroke = .1, size = 1) +
    
    geom_boxplot(aes(x = caller, y = value, fill = caller),
                 width = .3, alpha = .5, outlier.shape = NA) +
    
    scale_fill_manual(
        values = c(
            "VarDict" = "#8d43ae",
            "VarScan" = "#439aae",
            "Freebayes" = "#ae8d43",
            "Mutect2" = "#ae4364",
            "LoFreq" = "#c974ba"
        )
    ) +
    
    
    scale_y_continuous(labels = scales::comma) +
    
    theme_minimal() +
    
    theme(
        legend.position = "none",
        
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 13),
        axis.text.x = element_text(face = "bold", size = 12),
        axis.text.y = element_text(face = "bold", size = 12),
        
        axis.line = element_line(),
        axis.ticks = element_line(),
        
        panel.grid = element_blank()
    ) +
    
    labs(
        y = "Coverage (No. of reads)"
    )

#AF
s2 = df[, c(
    "caller",
    "POS", 
    "Caller AF"
), with = FALSE] 


p = rbind(s2) |> 
    unique() |>
    melt(id.vars = c("caller", "POS"), variable.factor = FALSE, value.factor = FALSE)

p$caller = p$caller |> factor(levels = c("Freebayes", 
                                         "LoFreq", 
                                         "Mutect2", 
                                         "VarDict", 
                                         "VarScan"))

fp2 = ggplot(data = p) +
    
    geom_point(aes(x = caller, y = value, fill = caller),
               position = position_jitternormal(sd_x = .025, sd_y = 0),
               shape = 21, stroke = .1, size = 1) +
    
    geom_boxplot(aes(x = caller, y = value, fill = caller),
                 width = .3, alpha = .5, outlier.shape = NA) +

    scale_fill_manual(
        values = c(
            "VarDict" = "#8d43ae",
            "VarScan" = "#439aae",
            "Freebayes" = "#ae8d43",
            "Mutect2" = "#ae4364",
            "LoFreq" = "#c974ba"
        )
    ) +
    

    scale_y_continuous(labels = scales::percent, 
                       trans = "log10", 
                       #limits = c(.00001, 1),
                       breaks = c(0.01, 0.10, 0.25, 0.50,  1)) +
    
    
    theme_minimal() +
    
    theme(
        legend.position = "none",
        
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 13),
        axis.text.x = element_text(face = "bold", size = 12),
        axis.text.y = element_text(face = "bold", size = 12),
        
        axis.line = element_line(),
        axis.ticks = element_line(),
        
        panel.grid = element_blank()
    ) +
    
    labs(
        y = "Allele Frequency"
    )

#patchwork

fp_paper = (fp1 | fp2)  + 
    plot_layout(heights = c(3, 2)) +
    plot_annotation(tag_levels = "A"
                    # title = "FP Variants",
                    # theme = theme(plot.title = element_text(size = 20, 
                    #                                         hjust = 0.5, 
                    #                                         vjust = -1))
                    ) &
    theme(
        plot.tag = element_text(face = "bold"),
        plot.margin = margin(5, 10, 10, 5)
    )


ggsave(
    plot = fp_paper, filename = paste0(folder,"/Plots/Final_", name,"_SNVs_FP.jpeg"),
    width = 14, height = 12, units = "in", dpi = 600
)

ggsave(
    plot = fp1, filename = paste0(folder,"/Plots/Final_", name,"_SNVs_FP_DP.jpeg"),
    width = 14, height = 12, units = "in", dpi = 600
)

#FN----------------------------------------------------------------------------
#FN
folder = 'results/'
name = "Merged"

df1 = list()

df1[["VarDict"]] = fread(paste0(folder, "Merged_VarDict_snvs_FN.tsv"))
df1[["VarScan"]] = fread(paste0(folder, "Merged_Varscan_snvs_FN.tsv"))
df1[["Freebayes"]] = fread(paste0(folder, "Merged_Freebayes_snvs_FN.tsv"))
df1[["Mutect2"]] = fread(paste0(folder, "Merged_Mutect2_snvs_FN.tsv"))
df1[["LoFreq"]] = fread(paste0(folder, "Merged_LoFreq_snvs_FN.tsv"))


df1 = df1 |> lapply(function(x) {
    
    colnames(x) = c(
        "POS",	"Ground Truth REF",	"Ground Truth DP",	
        "Ground Truth ALT",	"Count",	"Ground Truth AF",
        "mut",	"type"
    )
    
    return(x)
    
}) |> rbindlist(idcol = "caller")


df1[which(df1$`Caller ALT` == "0")]$`Caller ALT` = NA

#DP
s3 = df1[, c(
    "caller",
    "POS", 
    "Ground Truth DP"
), with = FALSE] 


q = rbind(s3) |> 
    unique() |>
    melt(id.vars = c("caller", "POS"), variable.factor = FALSE, value.factor = FALSE)

q$caller = q$caller |> factor(levels = c("Freebayes", "LoFreq", "Mutect2", "VarDict", "VarScan"))

FN1 = ggplot(data = q) +
    
    geom_point(aes(x = caller, y = value, fill = caller),
               position = position_jitternormal(sd_x = .025, sd_y = 0),
               shape = 21, stroke = .1, size = 1) +
    
    geom_boxplot(aes(x = caller, y = value, fill = caller),
                 width = .3, alpha = .5, outlier.shape = NA) +
    
    scale_fill_manual(
        values = c(
            "VarDict" = "#8d43ae",
            "VarScan" = "#439aae",
            "Freebayes" = "#ae8d43",
            "Mutect2" = "#ae4364",
            "LoFreq" = "#c974ba"
        )
    ) +
    
    
    scale_y_continuous(labels = scales::comma) +
    
    theme_minimal() +
    
    theme(
        legend.position = "none",
        
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 13),
        axis.text.x = element_text(face = "bold", size = 12),
        axis.text.y = element_text(face = "bold", size = 12),
        
        axis.line = element_line(),
        axis.ticks = element_line(),
        
        panel.grid = element_blank()
    ) +
    
    labs(
        y = "Coverage (No. of reads)"
    )


#AF
s4 = df1[, c(
    "caller",
    "POS", 
    "Ground Truth AF"
), with = FALSE] 


p = rbind(s4) |> 
    unique() |>
    melt(id.vars = c("caller", "POS"), variable.factor = FALSE, value.factor = FALSE)

p$caller = p$caller |> factor(levels = c("Freebayes", 
                                         "LoFreq", 
                                         "Mutect2", 
                                         "VarDict", 
                                         "VarScan"))

FN2 = ggplot(data = p) +
    
    geom_point(aes(x = caller, y = value, fill = caller),
               position = position_jitternormal(sd_x = .025, sd_y = 0),
               shape = 21, stroke = .1, size = 1) +
    
    geom_boxplot(aes(x = caller, y = value, fill = caller),
                 width = .3, alpha = .5, outlier.shape = NA) +
    
    scale_fill_manual(
        values = c(
            "VarDict" = "#8d43ae",
            "VarScan" = "#439aae",
            "Freebayes" = "#ae8d43",
            "Mutect2" = "#ae4364",
            "LoFreq" = "#c974ba"
        )
    ) +
    
    
    scale_y_continuous(labels = scales::percent, 
                       trans = "log10", 
                       #limits = c(.00001, 1),
                       breaks = c(0.01, 0.10, 0.25, 0.50,  1)) +
    
    
    theme_minimal() +
    
    theme(
        legend.position = "none",
        
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 13),
        axis.text.x = element_text(face = "bold", size = 12),
        axis.text.y = element_text(face = "bold", size = 12),
        
        axis.line = element_line(),
        axis.ticks = element_line(),
        
        panel.grid = element_blank()
    ) +
    
    labs(
        y = "Allele Frequency"
    )

#patchwork
FN_paper = (FN1 | FN2)  + 
    plot_layout(heights = c(3, 2)) +
    plot_annotation(tag_levels = "A"
                    # title = "FN Variants",
                    # theme = theme(plot.title = element_text(size = 20, 
                    #                                         hjust = 0.5, 
                    #                                         vjust = -1))
                    ) &
    theme(
        plot.tag = element_text(face = "bold"),
        plot.margin = margin(5, 10, 10, 5)
    )


ggsave(
    plot = FN_paper, filename = paste0(folder,"/Plots/Final_", name,"_SNVs_FN.jpeg"),
    width = 14, height = 12, units = "in", dpi = 600
)







ggsave(
    plot = FN1, filename = paste0(folder,"/Plots/Final_", name,"_SNVs_FN_DP.jpeg"),
    width = 14, height = 12, units = "in", dpi = 600
)


#BOTH
FP_FN_all_paper = (fp_paper / FN_paper)  + 
    plot_layout(heights = c(1, 1)) +
    plot_annotation(tag_levels = "A",
                    title = folder,
                    
    ) &
    theme(
        plot.tag = element_text(face = "bold"),
        plot.margin = margin(5, 5, 5, 5)
    )


ggsave(
    plot = FP_FN_all_paper, filename = paste0(folder,"/Plots/Final_", name,"_SNVs_FP_FN_all_new.jpeg"),
    width = 14, height = 12, units = "in", dpi = 600
)

