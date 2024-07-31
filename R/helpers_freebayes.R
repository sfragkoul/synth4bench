#'A script, written in R, where all the appropriate functions for 
#'the analysis of Freebayes are located.
#'
#'
#'Authors: Nikos Pechlivanis(github:npechl), Stella Fragkouli(github:sfragkoul)
#'
#'

read_vcf_freebayes <- function(path, gt, merged_file) {
  #takes two files and produce a caller vcf file in a certain format   
  vcf <- read.vcfR( paste0(path, "/", merged_file, "_Freebayes_norm.vcf"), verbose = FALSE )
  
  vcf_df <- vcf |>
    merge_freebayes(gt) |>
    clean_freebayes()
  
  return(vcf_df)
  
}

plot_synth4bench_freebayes <- function(df, vcf_GT, vcf_caller, merged_file) {
    #plotting function
    out1 = bar_plots_freebayes(df)
    out2 = density_plot_freebayes(df)
    out3 = bubble_plots_freebayes(df)
    out4 = venn_plot_freebayes(vcf_GT, vcf_caller)
    
    multi2 = out2$groundtruth / out2$Freebayes &
        
        theme(
            plot.margin = margin(10, 10, 10, 10)
        )
    
    
    ann1 = (out1$coverage + theme(plot.margin = margin(r = 50))) + 
        (out1$allele + theme(plot.margin = margin(r = 50))) + 
        multi2 +
        
        plot_layout(
            widths = c(1, 1, 3)
        )
    
    
    ann2 = out3 + out4 +
        
        plot_layout(
            widths = c(2, 1)
        )
    
    
    multi = ann1 / ann2 +
        
        plot_layout(heights = c(1.5, 1)) + 
        plot_annotation(title = merged_file)
    
    
    return(list(multi, out4))
}


merge_freebayes <- function(freebayes_somatic_vcf, merged_gt) {
    #return cleaned vcf
    freebayes_s0  = freebayes_somatic_vcf |> vcfR::getFIX() |> as.data.frame() |> setDT()
    freebayes_s1  = freebayes_somatic_vcf |> extract_gt_tidy() |> setDT()
    freebayesgatk_s21 = freebayes_somatic_vcf |> extract_info_tidy() |> setDT()
    
    freebayes_somatic = cbind(freebayes_s0[freebayes_s1$Key, ], freebayes_s1)
    
    
    #Merge everything into a common file
    merged_gt$POS = as.character(merged_gt$POS)
    
    merged_bnch = merge(merged_gt, freebayes_somatic,  by = "POS", all.x = TRUE)
    
    merged_bnch$POS = as.numeric(merged_bnch$POS)
    
    merged_bnch = merged_bnch[order(POS)]
    
    colnames(merged_bnch) = c(
        "POS",	"Ground Truth REF",	"Ground Truth DP",
        "Ground Truth ALT", "Ground Truth AD", 
        "Ground Truth AF", "Freebayes CHROM", "Freebayes ID",
        "Freebayes REF", "Freebayes ALT", "Freebayes QUAL", "Freebayes FILTER",
        "Freebayes key", "Freebayes Indiv", "Freebayes GT", "Freebayes GQ",
        "Freebayes GL", "Freebayes DP", "Freebayes RO", "Freebayes QR", 
        "Freebayes AO", "Freebayes QA", "Freebayes alleles"
    )
    
    return(merged_bnch)
    
}


clean_freebayes <- function(df) {
    #function to produce the caller's reported variants in the desired format 
    df2 = df[, c(
        "POS",
        
        "Ground Truth REF",
        "Ground Truth ALT",
        "Ground Truth DP",
        "Ground Truth AF",
        
        "Freebayes REF", 
        "Freebayes ALT", 
        "Freebayes DP",
        "Freebayes AO"
    ), with = FALSE]
    
    
    
    df2 = df2[, by = c(
        "POS",
        "Ground Truth REF",
        "Ground Truth DP",
        "Freebayes REF", 
        "Freebayes ALT", 
        "Freebayes DP",
        "Freebayes AO"
        
    ), .(
        "Ground Truth ALT" = `Ground Truth ALT` |> tstrsplit(",") |> unlist(),
        "Ground Truth AF"  = `Ground Truth AF` |> tstrsplit(",") |> unlist()
    )]



    freebayes_alt = str_split(df2$`Freebayes ALT`, ",")
    freebayes_ao = str_split(df2$`Freebayes AO`, ",")


    cln = mapply(
        function(x, y, z) {

            index = which(y == x)

            return(
                c(y[index], z[index])
            )

        },

        df2$`Ground Truth ALT`, freebayes_alt, freebayes_ao
    )


    df2$`Freebayes ALT` = cln |> lapply(function(x) { return(x [1]) }) |> unlist()
    df2$`Freebayes AO`  = cln |> lapply(function(x) { return(x [2]) }) |> unlist()
    
    df2[which(is.na(`Freebayes AO`))]$`Freebayes DP` = NA
    df2[which(is.na(`Freebayes AO`))]$`Freebayes REF` = NA
    
    df2 = df2[, c(
        "POS", 
        "Ground Truth REF",
        "Ground Truth ALT",
        "Ground Truth DP",
        "Ground Truth AF",
        "Freebayes REF", 
        "Freebayes ALT", 
        "Freebayes DP",
        "Freebayes AO"
    ), with = FALSE]
    
    #AF = AO/DP
    df2$"Freebayes AF" =  as.numeric(format(round(as.numeric(df2$"Freebayes AO")/ as.numeric(df2$"Freebayes DP"), 3)))
    
    return(df2)
    
}


bar_plots_freebayes <- function(q) {
    #function to produce variants' barplots for coverage and AF
    q[which(q$`Freebayes ALT` == "")]$`Freebayes ALT` = NA
    
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
        
        geom_boxplot(aes(x = variable, y = value, fill = variable),
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
    
    
    #AF plot
    df = q[, c(
        "POS",
        "Ground Truth AF",
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
                "Ground Truth AF" = "#43ae8d",
                "Freebayes AF"      = "#ae8d43"
            )
        ) +
        
        scale_x_discrete(
            breaks = c("Ground Truth AF", "Freebayes AF"),
            labels = c("Ground Truth", "Freebayes")
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
    
    return(
        list(
            "coverage" = o1,
            "allele" = o2
        )
    )
    
}


density_plot_freebayes <- function(q) {
    #function to produce AF density plots
    q[which(q$`Freebayes ALT` == "")]$`Freebayes ALT` = NA
    
    df = q[, c(
        "POS", 
        "Ground Truth AF",
        "Ground Truth ALT",
        "Freebayes ALT",
        "Freebayes AF"
    ), with = FALSE] |>
        unique()
    
    #Ground Truth AF density plot
    
    o1 = ggplot(data = df[, 1:3], aes(x = `Ground Truth AF`)) +
        
        geom_density(aes(color = `Ground Truth ALT`, fill = `Ground Truth ALT`),
                     alpha = .5) +
        
        scale_x_continuous(expand = c(0, 0), breaks = c(.25, .5, .75, 1), limits = c(0, 1), labels = scales::percent) +
        scale_y_continuous(expand = c(0, 0), limits = c(0, 15), breaks = c(5, 10, 15)) +
        
        scale_fill_npg() +
        scale_color_npg() +
        
        facet_wrap(vars(`Ground Truth ALT`), nrow = 1) +
        
        theme_minimal() +
        
        theme(
            legend.position = "none",
            
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold", size = 13),
            strip.text = element_text(face = "bold", size = 13),
            axis.text.y = element_text(face = "bold", size = 13),
            axis.text.x = element_blank(),
            
            panel.spacing = unit(1, "lines"),
            
            panel.grid = element_line(linetype = "dashed")
        ) +
        
        labs(y = "Ground Truth (density)")
    
    #Caller AF density plot
    
    o2 = ggplot(data = df[which(!is.na(`Freebayes ALT`)), c(1, 4, 5)], aes(x = `Freebayes AF`)) +
        
        geom_density(aes(color = `Freebayes ALT`, fill = `Freebayes ALT`),
                     alpha = .5) +
        
        scale_x_continuous(expand = c(0, 0), breaks = c(.25, .5, .75, 1), limits = c(0, 1), labels = scales::percent) +
        scale_y_continuous(expand = c(0, 0), limits = c(0, 15), breaks = c(5, 10, 15))+
        
        scale_fill_npg() +
        scale_color_npg() +
        
        facet_wrap(vars(`Freebayes ALT`), nrow = 1) +
        
        theme_minimal() +
        
        theme(
            legend.position = "none",
            axis.title.x = element_text(face = "bold", size = 13),
            axis.title.y = element_text(face = "bold", size = 13),
            strip.text = element_blank(), # element_text(face = "bold", size = 13),
            axis.text.y = element_text(face = "bold", size = 13),
            axis.text.x = element_text(face = "bold", size = 11),
            
            panel.spacing = unit(1, "lines"),
            
            panel.grid = element_line(linetype = "dashed")
        ) +
        
        labs(
            x = "Allele Frequency",
            y = "Freebayes (density)"
        )
    
    
    return(
        list(
            "groundtruth" = o1,
            "Freebayes" = o2
        )
    )
    
}


bubble_plots_freebayes <- function(q) {
    #function to produce SNVs bubble plot
    # q[which(q$`Freebayes ALT` == "")]$`Freebayes ALT` = NA
    
    
    q1 = q[which(q$`Freebayes ALT` != "")]
    
    
    o = ggplot() +
        
        geom_point(
            data = q,
            aes(x = POS, y = `Ground Truth ALT`, size = `Ground Truth AF`, fill = "Ground Truth"),
            position = position_nudge(y = .15), shape = 21, stroke = .25
        ) +
        
        geom_point(
            data = q1,
            aes(x = POS, y = `Freebayes ALT`, size = `Freebayes AF`, fill = "Freebayes"),
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
                "Freebayes"      = "#ae8d43"
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
        
        scale_y_discrete(breaks = c("A", "C", "G", "T")) +
        
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
    
    return(o)
    
}


venn_plot_freebayes <- function(q, p) {
    #function to produce Venn plot for each caller
    vcf_GT = vcfR::getFIX(q) |> as.data.frame() |> setDT()
    vcf_GT$scenario = "GT"
    
    vcf_freebayes = vcfR::getFIX(p) |> as.data.frame() |> setDT()
    vcf_freebayes$scenario = "Freebayes"
    
    x = rbind(vcf_GT, vcf_freebayes)
    y = x[, c("CHROM", "POS", "REF", "ALT", "scenario"), with = FALSE]
    
    y$mut = paste(y$CHROM, y$POS, y$REF, y$ALT, sep = ":")
    
    y = split(y, y$scenario)
    
    y = list(
        'Ground Truth' = y$GT$mut,
        'Freebayes'         = y$Freebayes$mut
    )
    
    gr = ggvenn(y, fill_color = c("#43ae8d", "#ae8d43")) +
        
        coord_equal(clip = "off")
    
    return(gr)
}
