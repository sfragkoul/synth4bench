#'
#'
#'function to locate variants of 100% AF in the individual files and
#'search their POS of interest in the Merged bam file
#'
#'Authors: Nikos Pechlivanis(github:npechl), Stella Fragkouli(github:sfragkoul)
#'

#function to search the POS of interest from the caller's vcf file

read_vcf_VarScan <- function(path, gt) {
  
  vcf <- read.vcfR( path, verbose = FALSE )
  
  vcf_df = vcf |>
    merge_VarScan(gt) |>
    clean_VarScan()
  
  return(vcf_df)
  
}


plot_synth4bench_VarScan(df, vcf_GT, vcf_caller){
    
    out1 = bar_plots_VarScan(df)
    out2 = density_plot_VarScan(df)
    out3 = bubble_plots_VarScan(df)
    out4 = venn_plot_VarScan(vcf_read_GT, vcf_read_VarScan)
    
    multi2 = out2$groundtruth / out2$VarScan &
        
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
        plot_annotation(title = folder)
    
    return(list(multi, out4))
}


merge_VarScan <- function(VarScan_somatic_vcf, merged_gt) {
    
    VarScan_s0  = VarScan_somatic_vcf |> vcfR::getFIX() |> as.data.frame() |> setDT()
    #VarScan_s1  = VarScan_somatic_vcf |> extract_gt_tidy() |> setDT()
    VarScan_s2 = VarScan_somatic_vcf |> extract_info_tidy() |> setDT()
    VarScan_s2 = VarScan_s2[,c( "DP", "Pvalue", "AF" )]
    
    VarScan_s0 = VarScan_s0[which(VarScan_s2$AF>0.0)]
    VarScan_s2 = VarScan_s2[which(VarScan_s2$AF>0.0)]
    
    
    VarScan_somatic = cbind(VarScan_s0, VarScan_s2)
    
    
    #Merge everything into a common file-------------------------------------------
    merged_gt$POS = as.character(merged_gt$POS)
    
    merged_bnch = merge(merged_gt, VarScan_somatic,  by = "POS", all.x = TRUE)
    
    merged_bnch$POS = as.numeric(merged_bnch$POS)
    
    merged_bnch = merged_bnch[order(POS)]
    
    colnames(merged_bnch) = c(
        "POS",	"Ground Truth REF",	"Ground Truth DP",
        "Ground Truth ALT", "Ground Truth AD", 
        "Ground Truth AF", "CHROM", "ID", "VarScan REF",	
        "VarScan ALT", "VarScan QUAL",	"VarScan FILTER", "VarScan DP", "Pvalue","VarScan AF"
    )
    
    return(merged_bnch)
    
}

#function to produce the caller's reported variants in the desired format 
clean_VarScan <- function(df) {
    
    df2 = df[, c(
        "POS",
        
        "Ground Truth REF",
        "Ground Truth ALT",
        "Ground Truth DP",
        "Ground Truth AF",
        
        "VarScan REF", 
        "VarScan ALT", 
        "VarScan DP",
        "VarScan AF"
    ), with = FALSE]
    
    
    
    df2 = df2[, by = c(
        "POS",
        "Ground Truth REF",
        "Ground Truth DP",
        "VarScan REF", 
        "VarScan ALT", 
        "VarScan DP",
        "VarScan AF"
        
    ), .(
        "Ground Truth ALT" = `Ground Truth ALT` |> tstrsplit(",") |> unlist(),
        "Ground Truth AF"  = `Ground Truth AF` |> tstrsplit(",") |> unlist()
    )]



    VarScan_alt = str_split(df2$`VarScan ALT`, ",")
    VarScan_af = str_split(df2$`VarScan AF`, ",")


    cln = mapply(
        function(x, y, z) {

            index = which(y == x)

            return(
                c(y[index], z[index])
            )

        },

        df2$`Ground Truth ALT`, VarScan_alt, VarScan_af
    )


    df2$`VarScan ALT` = cln |> lapply(function(x) { return(x [1]) }) |> unlist()
    df2$`VarScan AF`  = cln |> lapply(function(x) { return(x [2]) }) |> unlist()
    
    df2[which(is.na(`VarScan AF`))]$`VarScan DP` = NA
    df2[which(is.na(`VarScan AF`))]$`VarScan REF` = NA
    
    df2 = df2[, c(
        "POS", 
        "Ground Truth REF",
        "Ground Truth ALT",
        "Ground Truth DP",
        "Ground Truth AF",
        "VarScan REF", 
        "VarScan ALT", 
        "VarScan DP",
        "VarScan AF"
    ), with = FALSE]

    return(df2)
    
}

#function to produce variants' barplots for coverage and AF
bar_plots_VarScan <- function(q) {
    
    q[which(q$`VarScan ALT` == "")]$`VarScan ALT` = NA
    
    # plot 1 ------------------------
    
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
        
        geom_boxplot(aes(x = variable, y = value, fill = variable),
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
    
    
    # plot 2 ---------------------
    
    df = q[, c(
        "POS",
        "Ground Truth AF",
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
                "Ground Truth AF" = "#43ae8d",
                "VarScan AF"      = "#439aae"
            )
        ) +
        
        scale_x_discrete(
            breaks = c("Ground Truth AF", "VarScan AF"),
            labels = c("Ground Truth", "VarScan")
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
    
    # return -------------
    
    return(
        list(
            "coverage" = o1,
            "allele" = o2
        )
    )
    
}

#function to produce AF density plots
density_plot_VarScan <- function(q) {
    
    q[which(df$`VarScan ALT` == "")]$`VarScan ALT` = NA
    
    df = q[, c(
        "POS", 
        "Ground Truth AF",
        "Ground Truth ALT",
        "VarScan ALT",
        "VarScan AF"
    ), with = FALSE] |>
        unique()
    
    # plot 1 ---------------------------
    
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
    
    # plot 2 ----------------------
    
    o2 = ggplot(data = df[which(!is.na(`VarScan ALT`)), c(1, 4, 5)], aes(x = `VarScan AF`)) +
        
        geom_density(aes(color = `VarScan ALT`, fill = `VarScan ALT`),
                     alpha = .5) +
        
        scale_x_continuous(expand = c(0, 0), breaks = c(.25, .5, .75, 1), limits = c(0, 1), labels = scales::percent) +
        scale_y_continuous(expand = c(0, 0), limits = c(0, 15), breaks = c(5, 10, 15))+
        
        scale_fill_npg() +
        scale_color_npg() +
        
        facet_wrap(vars(`VarScan ALT`), nrow = 1) +
        
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
            y = "VarScan (density)"
        )
    
    
    # return -----------------
    
    return(
        list(
            "groundtruth" = o1,
            "VarScan" = o2
        )
    )
    
}

#function to produce SNVs bubble plot
bubble_plots_VarScan <- function(q) {
    
    # q[which(q$`VarScan ALT` == "")]$`VarScan ALT` = NA
    
    
    q1 = q[which(q$`VarScan ALT` != "")]
    
    
    o = ggplot() +
        
        geom_point(
            data = q,
            aes(x = POS, y = `Ground Truth ALT`, size = `Ground Truth AF`, fill = "Ground Truth"),
            position = position_nudge(y = .15), shape = 21, stroke = .25
        ) +
        
        geom_point(
            data = q1,
            aes(x = POS, y = `VarScan ALT`, size = `VarScan AF`, fill = "VarScan"),
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
                "VarScan"      = "#439aae"
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

#function to produce Venn plot for each caller
venn_plot_VarScan <- function(q, p) {
    
    vcf_GT = vcfR::getFIX(q) |> as.data.frame() |> setDT()
    vcf_GT$scenario = "GT"
    
    vcf_VarScan = vcfR::getFIX(p) |> as.data.frame() |> setDT()
    vcf_VarScan$scenario = "VarScan"
    
    x = rbind(vcf_GT, vcf_VarScan)
    y = x[, c("CHROM", "POS", "REF", "ALT", "scenario"), with = FALSE]
    
    y$mut = paste(y$CHROM, y$POS, y$REF, y$ALT, sep = ":")
    
    y = split(y, y$scenario)
    
    y = list(
        'Ground Truth' = y$GT$mut,
        'VarScan'         = y$VarScan$mut
    )
    
    gr = ggvenn(y, fill_color = c("#43ae8d", "#439aae")) +
        
        coord_equal(clip = "off")
    
    return(gr)
}
