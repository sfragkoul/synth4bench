
fp_snvs_gatk <- function(Mutect2_somatic_snvs, pick_gt, gt_all){#term snvs is redundant
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
    fp_var <- fp_snvs_gatk(Mutect2_somatic_snvs, pick_gt, gt_all)
    
    return(fp_var)
}

final_fn_snvs_gatk <- function(path, merged_file, pick_gt){
    
    Mutect2_somatic <- load_gatk_vcf(path, merged_file)
    Mutect2_somatic_snvs <- select_snvs(Mutect2_somatic)
    fn_var <- define_fn(Mutect2_somatic_snvs, pick_gt)
    colnames(fn_var) = c("POS", "Ground Truth REF", "Ground Truth DP",
                         "Ground Truth ALT", "Count", "Ground Truth AF",
                         "mut", "type")
    return(fn_var)
}

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
    
    return(o1)
}

fp_af_barplot_gatk <- function(q){
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

plot_snvs_FP_gatk <- function(df, merged_file) {
    #plotting function
    out1 = fp_violin_plots_gatk(df)
    out2 = fp_af_barplot_gatk(df)
    
    multi = out1 + out2 +
        
        plot_layout(
            widths = c(1, 1)
        )
    return(multi)
}

