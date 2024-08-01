source("R/libraries.R")
source("R/common_helpers.R")

# gt_all = load_gt_report("results/", "Merged_auto")$all
# gt_snvs = load_gt_report("results/", "Merged_auto")$snvs
# pick_gt = load_gt_vcf("results/", "Merged_auto")

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
fp_snvs_gatk <- function(Mutect2_somatic_snvs, pick_gt, gt_all){
    #find MUtect2 FP variants
    fp_var = define_fp(Mutect2_somatic_snvs, pick_gt)
    fp_var$gt_AF = as.numeric(fp_var$gt_AF)
    colnames(fp_var) = c("CHROM", "POS","ID", "Mutect2 REF",	
                         "Mutect2 ALT", "Mutect2 QUAL",	"Mutect2 FILTER",
                         "key", "Indiv", "Mutect2 AD", "Mutect2 AF",
                         "Mutect2 DP", "gt_F1R2", "gt_F2R1", "gt_FAD",	
                         "gt_GQ", "gt_GT",	"gt_PGT",	"gt_PID",	"gt_PL",
                         "gt_PS",	"gt_SB",	"gt_GT_alleles", "mut")
    
    #find DP of FP variants'  location in GT
    tmp = gt_all[which(POS %in% unique(fp_var$POS))]
    a = unique(tmp, by = "POS")
    #to include the presence multiple variants in a POS
    index = match(fp_var$POS, a$POS)
    fp_var$`Ground Truth DP` = a[index]$DP
    fp_var$`DP Percentage` = fp_var$`Mutect2 DP`/fp_var$`Ground Truth DP`
    fp_var$type = "FP"
    return(fp_var)
}

final_fp_snvs_gatk <- function(path, merged_file, pick_gt, gt_all){
    
    Mutect2_somatic <- load_gatk_vcf(path, merged_file)
    Mutect2_somatic_snvs <-select_snvs(Mutect2_somatic)
    fp_var = fp_snvs_gatk(Mutect2_somatic_snvs, pick_gt, gt_all)
    
    return(fp_var)
}

fp_var = final_fp_snvs_gatk("results/", "Merged_auto", pick_gt, gt_all)

# fwrite(
#     fp_var, paste0("results/", "Merged_auto_", "Mutect2_", "snvs_FP.tsv"),
#     row.names = FALSE, quote = FALSE, sep = "\t"
# )


#FN
final_fn_snvs_gatk <- function(path, merged_file, pick_gt){
    
    Mutect2_somatic <- load_gatk_vcf(path, merged_file)
    Mutect2_somatic_snvs <-select_snvs(Mutect2_somatic)
    fn_var = define_fn(Mutect2_somatic_snvs, pick_gt)
    colnames(fn_var) = c("POS", "Ground Truth REF", "Ground Truth DP", 
                         "Ground Truth ALT", "Count", "Ground Truth AF", "mut", "type")
    
    return(fp_var)
}

final_fn_snvs_gatk("results/", "Merged_auto", pick_gt)



# fwrite(
#     fn_var, paste0("results/", "Merged_auto_", "Mutect2_", "snvs_FN.tsv"),
#     row.names = FALSE, quote = FALSE, sep = "\t"
# )

#tp_var = define_tp(Mutect2_somatic_snvs, pick_gt)

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
  plot = gatk_multi3, filename = paste0("results/Plots/", "Merged_auto_", "Mutect2_snvs_FN_FP.png"),
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

fp_snvs_LoFreq <- function(LoFreq_somatic_snvs, pick_gt, gt_all){
    #find LoFreq FP variants
    fp_var = define_fp(LoFreq_somatic_snvs, pick_gt)
    fp_var$AF = as.numeric(fp_var$AF)
    colnames(fp_var) = c("CHROM", "POS","ID", "LoFreq REF",	
                         "LoFreq ALT", "LoFreq QUAL",	"LoFreq FILTER",
                         "LoFreq DP", "LoFreq AF", "mut")
    
    #find DP of FP variants'  location in GT
    tmp = gt_all[which(POS %in% unique(fp_var$POS))]
    a = unique(tmp, by = "POS")
    #to include the presence multiple variants in a POS
    index = match(fp_var$POS, a$POS)
    fp_var$`Ground Truth DP` = a[index]$DP
    fp_var$`DP Percentage` = fp_var$`LoFreq DP`/fp_var$`Ground Truth DP`
    fp_var$type = "FP"
    return(fp_var)
}

fp_violin_plots_LoFreq <- function(q) {
    #function to produce variants' barplots for coverage and AF
    #q[which(q$`LoFreq ALT` == "")]$`LoFreq ALT` = NA
    q$POS = as.numeric(q$POS)
    q$`Ground Truth DP` = as.numeric(q$`Ground Truth DP`)
    q$`LoFreq DP` = as.numeric(q$`LoFreq DP`)
    
    #DP plot
    df = q[, c(
        "POS", 
        "Ground Truth DP",
        "LoFreq DP"
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
                "LoFreq DP"      = "#c974ba"
            )
        ) +
        
        scale_x_discrete(
            breaks = c("Ground Truth DP", "LoFreq DP"),
            labels = c("Ground Truth", "LoFreq")
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
        "LoFreq AF"
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
                "LoFreq AF"      = "#c974ba"
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


LoFreq_somatic <- load_LoFreq_vcf("results/", "Merged_auto")
LoFreq_somatic_snvs <-select_snvs(LoFreq_somatic)
#FP
LoFreq_fp_var = fp_snvs_LoFreq(LoFreq_somatic_snvs, pick_gt, gt_all)

fwrite(
    fp_var, paste0("results/", "Merged_auto_", "LoFreq_", "snvs_FP.tsv"),
    row.names = FALSE, quote = FALSE, sep = "\t"
)

#FN
LoFreq_fn_var = define_fn(LoFreq_somatic_snvs, pick_gt)
colnames(LoFreq_fn_var) = c("POS", "Ground Truth REF", "Ground Truth DP", 
                     "Ground Truth ALT", "Count", "Ground Truth AF", "mut", "type")

fwrite(
    LoFreq_fn_var, paste0("results/", "Merged_auto_", "LoFreq_", "snvs_FN.tsv"),
    row.names = FALSE, quote = FALSE, sep = "\t"
)

#Plots

LoFreq_fp_plot1 <- fp_violin_plots_LoFreq(LoFreq_fp_var)
LoFreq_fp_plot2 <- fp_af_barplot(LoFreq_fp_var)
LoFreq_fn_plot1 <- fn_dp_barplot(LoFreq_fn_var, caller = "LoFreq")
LoFreq_fn_plot2 <- fn_af_barplot(LoFreq_fn_var, caller = "LoFreq")

LoFreq_multi1 = LoFreq_fp_plot1 + LoFreq_fp_plot2 +
    
    plot_layout(
        widths = c(1, 1)
    )

LoFreq_multi2 = LoFreq_fn_plot1 + LoFreq_fn_plot2 +
    
    plot_layout(
        widths = c(1, 1)
    )

LoFreq_multi3 = LoFreq_multi1 / LoFreq_multi2 &
    
    theme(
        plot.margin = margin(10, 10, 10, 10)
    )

ggsave(
    plot = LoFreq_multi3, filename = paste0("results/Plots/", "Merged_auto_", "LoFreq_snvs_FN_FP.png"),
    width = 16, height = 12, units = "in", dpi = 600
)



#VarDict-----------------------------------------------------------------------
load_VarDict_vcf <- function(path, merged_file){
    #function to load caller vcf
    VarDict_somatic_vcf <- read.vcfR( paste0(path, merged_file, 
                                             "_VarDict_norm.vcf"), verbose = FALSE )
    VarDict_s0  = VarDict_somatic_vcf |> vcfR::getFIX() |> as.data.frame() |> setDT()
    #VarDict_s1  = VarDict_somatic_vcf |> extract_gt_tidy() |> setDT()
    VarDict_s2 = VarDict_somatic_vcf |> extract_info_tidy() |> setDT()
    VarDict_s2 = VarDict_s2[,c( "DP", "AF" )]
    VarDict_somatic = cbind(VarDict_s0, VarDict_s2)
    return(VarDict_somatic)
}

fp_snvs_VarDict <- function(VarDict_somatic_snvs, pick_gt, gt_all){
    #find VarDict FP variants
    fp_var = define_fp(VarDict_somatic_snvs, pick_gt)
    fp_var$AF = as.numeric(fp_var$AF)
    colnames(fp_var) = c("CHROM", "POS","ID", "VarDict REF",	
                         "VarDict ALT", "VarDict QUAL",	"VarDict FILTER",
                         "VarDict DP", "VarDict AF", "mut")
    
    #find DP of FP variants'  location in GT
    tmp = gt_all[which(POS %in% unique(fp_var$POS))]
    tmp = tmp[nchar(tmp$REF) == nchar(tmp$ALT)]
    a = unique(tmp, by = "POS")
    #to include the presence multiple variants in a POS
    index = match(fp_var$POS, a$POS)
    fp_var$`Ground Truth DP` = a[index]$DP
    fp_var$`DP Percentage` = fp_var$`VarDict DP`/fp_var$`Ground Truth DP`
    fp_var$type = "FP"
    return(fp_var)
}

fp_violin_plots_VarDict <- function(q) {
    #function to produce variants' barplots for coverage and AF
    #q[which(q$`VarDict ALT` == "")]$`VarDict ALT` = NA
    q$POS = as.numeric(q$POS)
    q$`Ground Truth DP` = as.numeric(q$`Ground Truth DP`)
    q$`VarDict DP` = as.numeric(q$`VarDict DP`)
    
    #DP plot
    df = q[, c(
        "POS", 
        "Ground Truth DP",
        "VarDict DP"
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
                "VarDict DP"      = "#8d43ae"
            )
        ) +
        
        scale_x_discrete(
            breaks = c("Ground Truth DP", "VarDict DP"),
            labels = c("Ground Truth", "VarDict")
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
        "VarDict AF"
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
                "VarDict AF"      = "#8d43ae"
            )
        ) +
        
        scale_x_discrete(
            labels = c("VarDict FP Variants")
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


VarDict_somatic <- load_VarDict_vcf("results/", "Merged_auto")
VarDict_somatic_snvs <-select_snvs(VarDict_somatic)


#FP
VarDict_fp_var = fp_snvs_VarDict(VarDict_somatic_snvs, pick_gt, gt_all)

fwrite(
    VarDict_fp_var, paste0("results/", "Merged_auto_", "VarDict_", "snvs_FP.tsv"),
    row.names = FALSE, quote = FALSE, sep = "\t"
)

#FN
VarDict_fn_var = define_fn(VarDict_somatic_snvs, pick_gt)
colnames(VarDict_fn_var) = c("POS", "Ground Truth REF", "Ground Truth DP", 
                             "Ground Truth ALT", "Count", "Ground Truth AF", "mut", "type")

fwrite(
    VarDict_fn_var, paste0("results/", "Merged_auto_", "VarDict_", "snvs_FN.tsv"),
    row.names = FALSE, quote = FALSE, sep = "\t"
)

#Plots

VarDict_fp_plot1 <- fp_violin_plots_VarDict(VarDict_fp_var)
VarDict_fp_plot2 <- fp_af_barplot(VarDict_fp_var)
VarDict_fn_plot1 <- fn_dp_barplot(VarDict_fn_var, caller = "VarDict")
VarDict_fn_plot2 <- fn_af_barplot(VarDict_fn_var, caller = "VarDict")

VarDict_multi1 = VarDict_fp_plot1 + VarDict_fp_plot2 +
    
    plot_layout(
        widths = c(1, 1)
    )

VarDict_multi2 = VarDict_fn_plot1 + VarDict_fn_plot2 +
    
    plot_layout(
        widths = c(1, 1)
    )

VarDict_multi3 = VarDict_multi1 / VarDict_multi2 &
    
    theme(
        plot.margin = margin(10, 10, 10, 10)
    )

ggsave(
    plot = VarDict_multi3, filename = paste0("results/Plots/", "Merged_auto_", "VarDict_snvs_FN_FP.png"),
    width = 16, height = 12, units = "in", dpi = 600
)
#VarScan-----------------------------------------------------------------------
load_VarScan_vcf <- function(path, merged_file){
    #function to load caller vcf
    VarScan_somatic_vcf <- read.vcfR( paste0(path, merged_file, 
                                             "_VarScan_norm.vcf"), verbose = FALSE )
    VarScan_s0  = VarScan_somatic_vcf |> vcfR::getFIX() |> as.data.frame() |> setDT()
    #VarScan_s1  = VarScan_somatic_vcf |> extract_gt_tidy() |> setDT()
    VarScan_s2 = VarScan_somatic_vcf |> extract_info_tidy() |> setDT()
    VarScan_s2 = VarScan_s2[,c( "DP", "AF" )]
    VarScan_somatic = cbind(VarScan_s0, VarScan_s2)
    return(VarScan_somatic)
}

fp_snvs_VarScan <- function(VarScan_somatic_snvs, pick_gt, gt_all){
    #find VarScan FP variants
    fp_var = define_fp(VarScan_somatic_snvs, pick_gt)
    fp_var$AF = as.numeric(fp_var$AF)
    colnames(fp_var) = c("CHROM", "POS","ID", "VarScan REF",	
                         "VarScan ALT", "VarScan QUAL",	"VarScan FILTER",
                         "VarScan DP", "VarScan AF", "mut")
    
    #find DP of FP variants'  location in GT
    tmp = gt_all[which(POS %in% unique(fp_var$POS))]
    tmp = tmp[nchar(tmp$REF) == nchar(tmp$ALT)]
    a = unique(tmp, by = "POS")
    #to include the presence multiple variants in a POS
    index = match(fp_var$POS, a$POS)
    fp_var$`Ground Truth DP` = a[index]$DP
    fp_var$`DP Percentage` = fp_var$`VarScan DP`/fp_var$`Ground Truth DP`
    fp_var$type = "FP"
    return(fp_var)
}

fp_violin_plots_VarScan <- function(q) {
    #function to produce variants' barplots for coverage and AF
    #q[which(q$`VarScan ALT` == "")]$`VarScan ALT` = NA
    q$POS = as.numeric(q$POS)
    q$`Ground Truth DP` = as.numeric(q$`Ground Truth DP`)
    q$`VarScan DP` = as.numeric(q$`VarScan DP`)
    
    #DP plot
    df = q[, c(
        "POS", 
        "Ground Truth DP",
        "VarScan DP"
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
                "VarScan DP"      = "#439aae"
            )
        ) +
        
        scale_x_discrete(
            breaks = c("Ground Truth DP", "VarScan DP"),
            labels = c("Ground Truth", "VarScan")
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
        "VarScan AF"
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
                "VarScan AF"      = "#439aae"
            )
        ) +
        
        scale_x_discrete(
            labels = c("VarScan FP Variants")
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


VarScan_somatic <- load_VarScan_vcf("results/", "Merged_auto")
VarScan_somatic_snvs <-select_snvs(VarScan_somatic)


#FP
VarScan_fp_var = fp_snvs_VarScan(VarScan_somatic_snvs, pick_gt, gt_all)

fwrite(
    VarScan_fp_var, paste0("results/", "Merged_auto_", "VarScan_", "snvs_FP.tsv"),
    row.names = FALSE, quote = FALSE, sep = "\t"
)

#FN
VarScan_fn_var = define_fn(VarScan_somatic_snvs, pick_gt)
colnames(VarScan_fn_var) = c("POS", "Ground Truth REF", "Ground Truth DP", 
                             "Ground Truth ALT", "Count", "Ground Truth AF", "mut", "type")

fwrite(
    VarScan_fn_var, paste0("results/", "Merged_auto_", "VarScan_", "snvs_FN.tsv"),
    row.names = FALSE, quote = FALSE, sep = "\t"
)

#Plots

VarScan_fp_plot1 <- fp_violin_plots_VarScan(VarScan_fp_var)
VarScan_fp_plot2 <- fp_af_barplot(VarScan_fp_var)
VarScan_fn_plot1 <- fn_dp_barplot(VarScan_fn_var, caller = "VarScan")
VarScan_fn_plot2 <- fn_af_barplot(VarScan_fn_var, caller = "VarScan")

VarScan_multi1 = VarScan_fp_plot1 + VarScan_fp_plot2 +
    
    plot_layout(
        widths = c(1, 1)
    )

VarScan_multi2 = VarScan_fn_plot1 + VarScan_fn_plot2 +
    
    plot_layout(
        widths = c(1, 1)
    )

VarScan_multi3 = VarScan_multi1 / VarScan_multi2 &
    
    theme(
        plot.margin = margin(10, 10, 10, 10)
    )

ggsave(
    plot = VarScan_multi3, filename = paste0("results/Plots/", "Merged_auto_", "VarScan_snvs_FN_FP.png"),
    width = 16, height = 12, units = "in", dpi = 600
)
#Freebayes---------------------------------------------------------------------
load_Freebayes_vcf <- function(path, merged_file){
    #function to load caller vcf
    Freebayes_somatic_vcf <- read.vcfR( paste0(path, merged_file, 
                                               "_Freebayes_norm.vcf"), verbose = FALSE )
    Freebayes_s0  = Freebayes_somatic_vcf |> vcfR::getFIX() |> as.data.frame() |> setDT()
    #Freebayes_s1  = Freebayes_somatic_vcf |> extract_gt_tidy() |> setDT()
    Freebayes_s2 = Freebayes_somatic_vcf |> extract_info_tidy() |> setDT()
    Freebayes_s2 = Freebayes_s2[,c( "DP", "AF" )]
    Freebayes_somatic = cbind(Freebayes_s0, Freebayes_s2)
    return(Freebayes_somatic)
}

fp_snvs_Freebayes <- function(Freebayes_somatic_snvs, pick_gt, gt_all){
    #find Freebayes FP variants
    fp_var = define_fp(Freebayes_somatic_snvs, pick_gt)
    fp_var$AF = as.numeric(fp_var$AF)
    colnames(fp_var) = c("CHROM", "POS","ID", "Freebayes REF",	
                         "Freebayes ALT", "Freebayes QUAL",	"Freebayes FILTER",
                         "Freebayes DP", "Freebayes AF", "mut")
    
    #find DP of FP variants'  location in GT
    tmp = gt_all[which(POS %in% unique(fp_var$POS))]
    tmp = tmp[nchar(tmp$REF) == nchar(tmp$ALT)]
    a = unique(tmp, by = "POS")
    #to include the presence multiple variants in a POS
    index = match(fp_var$POS, a$POS)
    fp_var$`Ground Truth DP` = a[index]$DP
    fp_var$`DP Percentage` = fp_var$`Freebayes DP`/fp_var$`Ground Truth DP`
    fp_var$type = "FP"
    return(fp_var)
}

fp_violin_plots_Freebayes <- function(q) {
    #function to produce variants' barplots for coverage and AF
    #q[which(q$`Freebayes ALT` == "")]$`Freebayes ALT` = NA
    q$POS = as.numeric(q$POS)
    q$`Ground Truth DP` = as.numeric(q$`Ground Truth DP`)
    q$`Freebayes DP` = as.numeric(q$`Freebayes DP`)
    
    #DP plot
    df = q[, c(
        "POS", 
        "Ground Truth DP",
        "Freebayes DP"
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
                "Freebayes DP"      = "#ae8d43"
            )
        ) +
        
        scale_x_discrete(
            breaks = c("Ground Truth DP", "Freebayes DP"),
            labels = c("Ground Truth", "Freebayes")
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
        "Freebayes AF"
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
                "Freebayes AF"      = "#ae8d43"
            )
        ) +
        
        scale_x_discrete(
            labels = c("Freebayes FP Variants")
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


Freebayes_somatic <- load_Freebayes_vcf("results/", "Merged_auto")
Freebayes_somatic_snvs <-select_snvs(Freebayes_somatic)


#FP
Freebayes_fp_var = fp_snvs_Freebayes(Freebayes_somatic_snvs, pick_gt, gt_all)

fwrite(
    Freebayes_fp_var, paste0("results/", "Merged_auto_", "Freebayes_", "snvs_FP.tsv"),
    row.names = FALSE, quote = FALSE, sep = "\t"
)

#FN
Freebayes_fn_var = define_fn(Freebayes_somatic_snvs, pick_gt)
colnames(Freebayes_fn_var) = c("POS", "Ground Truth REF", "Ground Truth DP", 
                               "Ground Truth ALT", "Count", "Ground Truth AF", "mut", "type")

fwrite(
    Freebayes_fn_var, paste0("results/", "Merged_auto_", "Freebayes_", "snvs_FN.tsv"),
    row.names = FALSE, quote = FALSE, sep = "\t"
)

#Plots

Freebayes_fp_plot1 <- fp_violin_plots_Freebayes(Freebayes_fp_var)
Freebayes_fp_plot2 <- fp_af_barplot(Freebayes_fp_var)
Freebayes_fn_plot1 <- fn_dp_barplot(Freebayes_fn_var, caller = "Freebayes")
Freebayes_fn_plot2 <- fn_af_barplot(Freebayes_fn_var, caller = "Freebayes")

Freebayes_multi1 = Freebayes_fp_plot1 + Freebayes_fp_plot2 +
    
    plot_layout(
        widths = c(1, 1)
    )

Freebayes_multi2 = Freebayes_fn_plot1 + Freebayes_fn_plot2 +
    
    plot_layout(
        widths = c(1, 1)
    )

Freebayes_multi3 = Freebayes_multi1 / Freebayes_multi2 &
    
    theme(
        plot.margin = margin(10, 10, 10, 10)
    )

ggsave(
    plot = Freebayes_multi3, filename = paste0("results/Plots/", "Merged_auto_", "Freebayes_snvs_FN_FP.png"),
    width = 16, height = 12, units = "in", dpi = 600
)
