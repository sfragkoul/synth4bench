library(data.table)
library(stringr)
library(vcfR)
library(ggplot2)
library(ggvenn)
library(ggforce)
library(patchwork)

#function "not in" def
`%ni%` <- Negate(`%in%`) 

load_gt_report <- function(path, merged_file) {
    #function to load Ground Truth bam-report 
    a <- paste0(path, merged_file, "_report.tsv") |>
      readLines() |>
      str_split(pattern = "\t", simplify = TRUE) |>
      as.data.frame() |> 
      setDT()
    
    a$V1 = NULL
    a$V5 = NULL
    
    colnames(a) = c("POS", "REF", "DP", paste0("ALT_", 1:(ncol(a) - 3)))
    
    a = melt(
      a, id.vars = c("POS", "REF", "DP"),
      variable.factor = FALSE, value.factor = FALSE,
      variable.name = "ALT", value.name = "Count"
    )
    
    a = a[which(Count != "")]
    
    a$POS = as.numeric(a$POS)
    a$DP = as.numeric(a$DP)
    
    a$ALT = str_split_i(a$Count, "\\:", 1)
    
    a$Count = str_split_i(a$Count, "\\:", 2) |>
      as.numeric()
    
    a$Freq = round(a$Count / a$DP, digits = 6)
    
    a = a[order(POS, -Count)]
    
    a = a[which(REF != a$ALT & Count != 0)]

    return(a)
}

load_gt_vcf <- function(path, merged_file){
    #function to load Ground Truth vcf
    ground_truth_vcf <- read.vcfR( paste0(path, merged_file, 
                                          "_ground_truth_norm.vcf"),
                                   verbose = FALSE )
    
    ground_truth_vcf  = ground_truth_vcf |> vcfR::getFIX() |> as.data.frame() |> setDT()
    
    pick_gt = nt_runs[which(nt_runs$POS %in% ground_truth_vcf$POS)]
    pick_gt$mut = paste(pick_gt$POS, 
                        pick_gt$REF, 
                        pick_gt$ALT, sep = ":")
    return(pick_gt)
}

select_snvs <- function(df){
    # select SNVs from caller based on length of REF and ALT
    snvs = df[nchar(df$REF) == nchar(df$ALT)]
    snvs = snvs[which(nchar(snvs$REF) <2), ]
    snvs$mut = paste(snvs$POS, snvs$REF, snvs$ALT, sep = ":")
    
    return(snvs)
}

define_fp <- function(caller, gt){
    #FP Variants
    fp_var = caller[which(caller$mut %ni% gt$mut)]
    
    return(fp_var)
}

define_fn <- function(caller, gt){
    #FN Variants
    fn_var = gt[which(gt$mut %ni% caller$mut)]
    
    return(fn_var)
}

define_tp <- function(caller, gt){
    #FN Variants
    tp_var = caller[which(caller$mut %in% gt$mut)]
    
    return(tp_var)
}

gt_runs = load_gt_report("results/", "Merged_auto")
# select SNVs
nt_runs = gt_runs[which(ALT %in% c("A", "C", "G", "T")), ]
#filter DEPTH>2
nt_runs = nt_runs[which(nt_runs$Count >2), ]
pick_gt = load_gt_vcf("results/", "Merged_auto")


fn_dp_barplot <- function(q, caller){
  #FP DP plot
  df = q[, c(
    "POS", 
    "Ground Truth DP"
  ), with = FALSE] |>
    unique() |>
    melt(id.vars = "POS", variable.factor = FALSE, value.factor = FALSE)
  
  o3=ggplot(data = df) +
    
    geom_point(aes(x = variable, y = value, fill = variable),
               position = position_jitternormal(sd_x = .01, sd_y = 0),
               shape = 21, stroke = .1, size = 2.5) +
    
    geom_boxplot(aes(x = variable, y = value, fill = variable),
                 width = .25, alpha = .5, outlier.shape = NA) +
    
    scale_fill_manual(
      values = c(
        "Ground Truth DP" = "#43ae8d"
      )
    ) +
    
    scale_x_discrete(
      labels = c(paste0(caller, " FN Variants"))
    ) +
    
    scale_y_continuous(labels = scales::comma) +
    
    theme_minimal() +
    
    theme(
      legend.position = "none",
      
      axis.title.x = element_blank(),
      axis.title.y = element_text(face = "bold", size = 13),
      
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
  return(o3)
  
}
fn_af_barplot <- function(q, caller){
  #FP AF plot
  df = q[, c(
    "POS",
    "Ground Truth AF"
  ), with = FALSE] |>
    unique() |>
    
    melt(id.vars = "POS", variable.factor = FALSE, value.factor = FALSE)
  
  o4 = ggplot(data = df[which(!is.na(value) & value != 0)]) +
    
    geom_point(aes(x = variable, y = value, fill = variable),
               position = position_jitternormal(sd_x = .01, sd_y = 0),
               shape = 21, stroke = .1, size = 2.5) +
    
    geom_boxplot(aes(x = variable, y = value, fill = variable),
                 width = .25, alpha = .5, outlier.shape = NA) +
    
    scale_fill_manual(
      values = c(
        "Ground Truth AF" = "#43ae8d"
      )
    ) +
    
    scale_x_discrete(
      labels = c(paste0(caller, " FN Variants"))
    ) +
    
    scale_y_continuous(labels = scales::percent, trans = "log10") +
    
    theme_minimal() +
    
    theme(
      legend.position = "none",
      
      axis.title.x = element_blank(),
      axis.title.y = element_text(face = "bold", size = 13),
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
  return(o4)
  
}

#GATK--------------------------------------------------------------------------

load_gatk_vcf <- function(path, merged_file){
  #function to load caller vcf
  Mutect2_somatic_vcf <- read.vcfR( paste0(path, merged_file, 
                                           "_Mutect2_norm.vcf"), verbose = FALSE )
  
  Mutect2_s0  = Mutect2_somatic_vcf |> vcfR::getFIX() |> as.data.frame() |> setDT()
  Mutect2_s1  = Mutect2_somatic_vcf |> extract_gt_tidy() |> setDT()
  Mutect2gatk_s21 = Mutect2_somatic_vcf |> extract_info_tidy() |> setDT()
  Mutect2_somatic = cbind(Mutect2_s0[Mutect2_s1$Key, ], Mutect2_s1)
  return(Mutect2_somatic)
}


Mutect2_somatic <- load_gatk_vcf("results/", "Merged_auto")
Mutect2_somatic_snvs <-select_snvs(Mutect2_somatic)

#FP
fp_var = define_fp(Mutect2_somatic_snvs, pick_gt)
fp_var$gt_AF = as.numeric(fp_var$gt_AF)
colnames(fp_var) = c("CHROM", "POS","ID", "Mutect2 REF",	
                                   "Mutect2 ALT", "Mutect2 QUAL",	"Mutect2 FILTER",
                                   "key", "Indiv", "Mutect2 AD", "Mutect2 AF",
                                   "Mutect2 DP", "gt_F1R2", "gt_F2R1", "gt_FAD",	
                                   "gt_GQ", "gt_GT",	"gt_PGT",	"gt_PID",	"gt_PL",
                                   "gt_PS",	"gt_SB",	"gt_GT_alleles", "mut")

teeeeest = gt_runs[which(POS %in% unique(fp_var$POS))]
a = unique(teeeeest, by = "POS")
fp_var$`Ground Truth DP` = a$DP
fp_var$`DP Percentage` = fp_var$`Mutect2 DP`/fp_var$`Ground Truth DP`



# Plot the alluvial diagram---------------------------------------------------
ggplot(data = fp_var,
       aes(axis1 = `Mutect2 DP`, axis2 = `Ground Truth DP`)) +
    geom_alluvium(aes(fill = `Mutect2 DP`), width = 0.1) +
    geom_stratum(width = 0.1, fill = "grey", color = "black") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_x_discrete(limits = c("Mutect2", "Ground Truth"), expand = c(0.15, 0.05)) +
    theme_minimal() +
    labs(title = "Alluvial Diagram",
         x = "Dimensions")

#------------------------------------------
#FN
fn_var = define_fn(Mutect2_somatic_snvs, pick_gt)
colnames(fn_var) = c("POS", "Ground Truth REF", "Ground Truth DP", 
                "Ground Truth ALT", "Count", "Ground Truth AF", "mut")


tp_var = define_tp(Mutect2_somatic_snvs, pick_gt)

# fp_dp_barplot <- function(q ){
#   #FP DP plot
#   df = q[, c(
#     "POS", 
#     "Mutect2 DP"
#   ), with = FALSE] |>
#     unique() |>
#     melt(id.vars = "POS", variable.factor = FALSE, value.factor = FALSE)
#   
#   o1=ggplot(data = df) +
#     
#     geom_point(aes(x = variable, y = value, fill = variable),
#                position = position_jitternormal(sd_x = .01, sd_y = 0),
#                shape = 21, stroke = .1, size = 2.5) +
#     
#     geom_boxplot(aes(x = variable, y = value, fill = variable),
#                  width = .25, alpha = .5, outlier.shape = NA) +
#     
#     scale_fill_manual(
#       values = c(
#        # "Ground Truth DP" = "#43ae8d",
#         "Mutect2 DP"      = "#ae4364"
#       )
#     ) +
#     
#     scale_x_discrete(
#       labels = c("Mutect2 FP Variants")
#     ) +
#     
#     scale_y_continuous(labels = scales::comma) +
#     
#     theme_minimal() +
#     
#     theme(
#       legend.position = "none",
#       
#       axis.title.x = element_blank(),
#       axis.title.y = element_text(face = "bold", size = 13),
#       
#       axis.text.x = element_text(face = "bold", size = 13),
#       axis.text.y = element_text(face = "bold", size = 13),
#       
#       axis.line = element_line(),
#       axis.ticks = element_line(),
#       
#       panel.grid = element_blank(),
#       
#       plot.margin = margin(20, 20, 20, 20)
#     ) +
#     
#     labs(
#       y = "Coverage (No. of reads)"
#     )
#   return(o1)
# 
# }


fp_violin_plots_gatk <- function(q) {
    #function to produce variants' barplots for coverage and AF
    #q[which(q$`Mutect2 ALT` == "")]$`Mutect2 ALT` = NA
    q$POS = as.numeric(q$POS)
    q$`Ground Truth DP` = as.numeric(q$`Ground Truth DP`)
    q$`Mutect2 DP` = as.numeric(q$`Mutect2 DP`)
    
    #DP plot
    df = q[, c(
        "POS", 
        "Ground Truth DP",
        "Mutect2 DP"
    ), with = FALSE] |>
        unique() |>
        
        melt(id.vars = "POS", variable.factor = FALSE, value.factor = FALSE)
    
    o1 = ggplot(data = df) +
        
        geom_point(aes(x = variable, y = value, fill = variable),
                   position = position_jitternormal(sd_x = .01, sd_y = 0),
                   shape = 21, stroke = .1, size = 2.5) +
        
        geom_violin(aes(x = variable, y = value, fill = variable),
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
            axis.title.y = element_text(face = "bold", size = 13),
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
}

fp_af_barplot <- function(q){
  #FP AF plot
  df = q[, c(
    "POS",
    "Mutect2 AF"
  ), with = FALSE] |>
    unique() |>
    
    melt(id.vars = "POS", variable.factor = FALSE, value.factor = FALSE)
  
  o2 = ggplot(data = df[which(!is.na(value) & value != 0)]) +
    
    geom_point(aes(x = variable, y = value, fill = variable),
               position = position_jitternormal(sd_x = .01, sd_y = 0),
               shape = 21, stroke = .1, size = 2.5) +
    
    geom_boxplot(aes(x = variable, y = value, fill = variable),
                 width = .25, alpha = .5, outlier.shape = NA) +
    
    scale_fill_manual(
      values = c(
        #"Ground Truth AF" = "#43ae8d",
        "Mutect2 AF"      = "#ae4364"
      )
    ) +
    
    scale_x_discrete(
      labels = c("Mutect2 FP Variants")
    ) +
    
    scale_y_continuous(labels = scales::percent, trans = "log10") +
    
    theme_minimal() +
    
    theme(
      legend.position = "none",
      
      axis.title.x = element_blank(),
      axis.title.y = element_text(face = "bold", size = 13),
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
  return(o2)
  
}

gatk_fp_plot1 <- fp_violin_plots_gatk(fp_var)
gatk_fp_plot2 <- fp_af_barplot(fp_var)
gatk_fn_plot1 <- fn_dp_barplot(fn_var, caller = "Mutect2")
gatk_fn_plot2 <- fn_af_barplot(fn_var, caller = "Mutect2")

gatk_multi1 = gatk_fp_plot1 + gatk_fp_plot2 +
  
  plot_layout(
    widths = c(1, 1)
  )

gatk_multi2 = gatk_fn_plot1 + gatk_fn_plot2 +
  
  plot_layout(
    widths = c(1, 1)
  )

gatk_multi3 = gatk_multi1 / gatk_multi2 &

theme(
  plot.margin = margin(10, 10, 10, 10)
)

ggsave(
  plot = gatk_multi3, filename = paste0("results/Plots/", "Merged_auto_", "Mutect2_FN_FP.png"),
  width = 16, height = 12, units = "in", dpi = 600
)

#LoFreq------------------------------------------------------------------------

load_LoFreq_vcf <- function(path, merged_file){
  #function to load caller vcf
  LoFreq_somatic_vcf <- read.vcfR( paste0(path, merged_file, 
                                           "_LoFreq_norm.vcf"), verbose = FALSE )
  LoFreq_s0  = LoFreq_somatic_vcf |> vcfR::getFIX() |> as.data.frame() |> setDT()
  #LoFreq_s1  = LoFreq_somatic_vcf |> extract_gt_tidy() |> setDT()
  LoFreq_s2 = LoFreq_somatic_vcf |> extract_info_tidy() |> setDT()
  LoFreq_s2 = LoFreq_s2[,c( "DP", "AF" )]
  LoFreq_somatic = cbind(LoFreq_s0, LoFreq_s2)
  return(LoFreq_somatic)
}

LoFreq_somatic <- load_LoFreq_vcf("results/", "Merged_auto")
LoFreq_somatic_snvs <-select_snvs(LoFreq_somatic)
LoFreq_fp_var = define_fp(LoFreq_somatic_snvs, pick_gt)
#LoFreq_fp_var$gt_AF = as.numeric(fp_var$gt_AF)
LoFreq_fn_var = define_fn(LoFreq_somatic_snvs, pick_gt)
LoFreq_tp_var = define_tp(LoFreq_somatic_snvs, pick_gt)


lofreq_fp_af_barplot <- function(q){
  #FP AF plot
  df = q[, c(
    "POS",
    "AF"
  ), with = FALSE] |>
    unique() |>
    
    melt(id.vars = "POS", variable.factor = FALSE, value.factor = FALSE)
  
  o2 = ggplot(data = df[which(!is.na(value) & value != 0)]) +
    
    geom_point(aes(x = variable, y = value, fill = variable),
               position = position_jitternormal(sd_x = .01, sd_y = 0),
               shape = 21, stroke = .1, size = 2.5) +
    
    geom_boxplot(aes(x = variable, y = value, fill = variable),
                 width = .25, alpha = .5, outlier.shape = NA) +
    
    scale_fill_manual(
      values = c(
        #"Ground Truth AF" = "#43ae8d",
        "AF"      = "#c974ba"
      )
    ) +
    
    scale_x_discrete(
      labels = c("LoFreq FP Variants")
    ) +
    
    scale_y_continuous(labels = scales::percent, trans = "log10") +
    
    theme_minimal() +
    
    theme(
      legend.position = "none",
      
      axis.title.x = element_blank(),
      axis.title.y = element_text(face = "bold", size = 13),
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
  return(o2)
  
}

lofreq_fp_plot1 <- lofreq_fp_dp_barplot(LoFreq_fp_var)
lofreq_fp_plot2 <- lofreq_fp_af_barplot(LoFreq_fp_var)


lofreq_fn_plot1 <- fn_dp_barplot(LoFreq_fn_var, caller = "LoFreq")
lofreq_fn_plot2 <- fn_af_barplot(LoFreq_fn_var, caller = "LoFreq")

lofreq_multi1 = lofreq_fp_plot1 + lofreq_fp_plot2 +
  
  plot_layout(
    widths = c(1, 1)
  )

lofreq_multi2 = lofreq_fn_plot1 + lofreq_fn_plot2 +
  
  plot_layout(
    widths = c(1, 1)
  )

lofreq_multi3 = lofreq_multi1 / lofreq_multi2 &
  
  theme(
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(
  plot = lofreq_multi3, filename = paste0("results/Plots/", "Merged_auto_", "LoFreq_FN_FP.png"),
  width = 16, height = 12, units = "in", dpi = 600
)





#VarDict-----------------------------------------------------------------------

#VarScan-----------------------------------------------------------------------

#Freebayes---------------------------------------------------------------------