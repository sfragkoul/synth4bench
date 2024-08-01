#'A script, written in R, where all the appropriate functions for 
#'the analysis of VarDict are located.
#'
#'Authors: Nikos Pechlivanis(github:npechl), Stella Fragkouli(github:sfragkoul)
#'
#'

#TP SNVS-----------------------------------------------------------------------
read_vcf_VarDict <- function(path, gt, merged_file) {
  #takes two files and produce a caller vcf file in a certain format 
  vcf <- read.vcfR(paste0(path, "/", merged_file, "_VarDict_norm.vcf"), verbose = FALSE )
  
  vcf_df = vcf |>
    merge_VarDict(gt) |>
    clean_VarDict()
  
  return(vcf_df)
  
}


plot_synth4bench_VarDict <- function(df, vcf_GT, vcf_caller, merged_file){
    #plotting function
    out1 = bar_plots_VarDict(df)
    out2 = density_plot_VarDict(df)
    out3 = bubble_plots_VarDict(df)
    out4 = venn_plot_VarDict(vcf_GT, vcf_caller)
    
    multi2 = out2$groundtruth / out2$VarDict &
        
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

merge_VarDict <- function(VarDict_somatic_vcf, merged_gt) {
    #return cleaned vcf
    VarDict_s0  = VarDict_somatic_vcf |> vcfR::getFIX() |> as.data.frame() |> setDT()
    VarDict_s1  = VarDict_somatic_vcf |> extract_gt_tidy() |> setDT()
    VarDictgatk_s21 = VarDict_somatic_vcf |> extract_info_tidy() |> setDT()
    
    VarDict_somatic = cbind(VarDict_s0[VarDict_s1$Key, ], VarDict_s1)
    
    
    #Merge everything into a common file
    merged_gt$POS = as.character(merged_gt$POS)
    
    merged_bnch = merge(merged_gt, VarDict_somatic,  by = "POS", all.x = TRUE)
    
    merged_bnch$POS = as.numeric(merged_bnch$POS)
    
    merged_bnch = merged_bnch[order(POS)]
    
    colnames(merged_bnch) = c(
        "POS",	"Ground Truth REF",	"Ground Truth DP",
        "Ground Truth ALT", "Ground Truth AD", 
        "Ground Truth AF", "CHROM", "ID", "VarDict REF",	
        "VarDict ALT", "VarDict QUAL",	"VarDict FILTER",
        "key", "Indiv", "gt_GT", "VarDict DP", "gt_VD", "VarDict AD", 
        "VarDict AF", "gt_RD", "gt_ALD", "gt_GT_alleles"
    )
    
    return(merged_bnch)
    
}


clean_VarDict <- function(df) {
    #function to produce the caller's reported variants in the desired format 
    df2 = df[, c(
        "POS",
        
        "Ground Truth REF",
        "Ground Truth ALT",
        "Ground Truth DP",
        "Ground Truth AF",
        
        "VarDict REF", 
        "VarDict ALT", 
        "VarDict DP",
        "VarDict AF"
    ), with = FALSE]
    
    
    
    df2 = df2[, by = c(
        "POS",
        "Ground Truth REF",
        "Ground Truth DP",
        "VarDict REF", 
        "VarDict ALT", 
        "VarDict DP",
        "VarDict AF"
        
    ), .(
        "Ground Truth ALT" = `Ground Truth ALT` |> tstrsplit(",") |> unlist(),
        "Ground Truth AF"  = `Ground Truth AF` |> tstrsplit(",") |> unlist()
    )]



    VarDict_alt = str_split(df2$`VarDict ALT`, ",")
    VarDict_af = str_split(df2$`VarDict AF`, ",")


    cln = mapply(
        function(x, y, z) {

            index = which(y == x)

            return(
                c(y[index], z[index])
            )

        },

        df2$`Ground Truth ALT`, VarDict_alt, VarDict_af
    )


    df2$`VarDict ALT` = cln |> lapply(function(x) { return(x [1]) }) |> unlist()
    df2$`VarDict AF`  = cln |> lapply(function(x) { return(x [2]) }) |> unlist()
    
    df2[which(is.na(`VarDict AF`))]$`VarDict DP` = NA
    df2[which(is.na(`VarDict AF`))]$`VarDict REF` = NA
    
    df2 = df2[, c(
        "POS", 
        "Ground Truth REF",
        "Ground Truth ALT",
        "Ground Truth DP",
        "Ground Truth AF",
        "VarDict REF", 
        "VarDict ALT", 
        "VarDict DP",
        "VarDict AF"
    ), with = FALSE]

    return(df2)
    
}


bar_plots_VarDict <- function(q) {
    #function to produce variants' barplots for coverage and AF
    q[which(q$`VarDict ALT` == "")]$`VarDict ALT` = NA
    
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
        
        geom_boxplot(aes(x = variable, y = value, fill = variable),
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
    
    
    #AF plot
    df = q[, c(
        "POS",
        "Ground Truth AF",
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
                "Ground Truth AF" = "#43ae8d",
                "VarDict AF"      = "#8d43ae"
            )
        ) +
        
        scale_x_discrete(
            breaks = c("Ground Truth AF", "VarDict AF"),
            labels = c("Ground Truth", "VarDict")
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


density_plot_VarDict <- function(q) {
    #function to produce AF density plots
    q[which(q$`VarDict ALT` == "")]$`VarDict ALT` = NA
    
    df = q[, c(
        "POS", 
        "Ground Truth AF",
        "Ground Truth ALT",
        "VarDict ALT",
        "VarDict AF"
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
    
    o2 = ggplot(data = df[which(!is.na(`VarDict ALT`)), c(1, 4, 5)], aes(x = `VarDict AF`)) +
        
        geom_density(aes(color = `VarDict ALT`, fill = `VarDict ALT`),
                     alpha = .5) +
        
        scale_x_continuous(expand = c(0, 0), breaks = c(.25, .5, .75, 1), limits = c(0, 1), labels = scales::percent) +
        scale_y_continuous(expand = c(0, 0), limits = c(0, 15), breaks = c(5, 10, 15))+
        
        scale_fill_npg() +
        scale_color_npg() +
        
        facet_wrap(vars(`VarDict ALT`), nrow = 1) +
        
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
            y = "VarDict (density)"
        )
    
    return(
        list(
            "groundtruth" = o1,
            "VarDict" = o2
        )
    )
    
}


bubble_plots_VarDict <- function(q) {
    #function to produce SNVs bubble plot
    # q[which(q$`VarDict ALT` == "")]$`VarDict ALT` = NA
    
    
    q1 = q[which(q$`VarDict ALT` != "")]
    
    
    o = ggplot() +
        
        geom_point(
            data = q,
            aes(x = POS, y = `Ground Truth ALT`, size = `Ground Truth AF`, fill = "Ground Truth"),
            position = position_nudge(y = .15), shape = 21, stroke = .25
        ) +
        
        geom_point(
            data = q1,
            aes(x = POS, y = `VarDict ALT`, size = `VarDict AF`, fill = "VarDict"),
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
                "VarDict"      = "#8d43ae"
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


venn_plot_VarDict <- function(q, p) {
    #function to produce Venn plot for each caller
    vcf_GT = vcfR::getFIX(q) |> as.data.frame() |> setDT()
    vcf_GT$scenario = "GT"
    
    vcf_VarDict = vcfR::getFIX(p) |> as.data.frame() |> setDT()
    vcf_VarDict$scenario = "VarDict"
    
    x = rbind(vcf_GT, vcf_VarDict)
    y = x[, c("CHROM", "POS", "REF", "ALT", "scenario"), with = FALSE]
    
    y$mut = paste(y$CHROM, y$POS, y$REF, y$ALT, sep = ":")
    
    y = split(y, y$scenario)
    
    y = list(
        'Ground Truth' = y$GT$mut,
        'VarDict'         = y$VarDict$mut
    )
    
    gr = ggvenn(y, fill_color = c("#43ae8d", "#8d43ae")) +
        
        coord_equal(clip = "off")
    
    return(gr)
}



#FP & FN SNVS------------------------------------------------------------------

load_VarDict_vcf <- function(path, merged_file){
    #function to load caller vcf
    VarDict_somatic_vcf <- read.vcfR( paste0(path, "/",merged_file, 
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

final_fp_snvs_VarDict<- function(path, merged_file, pick_gt, gt_all){
    
    VarDict_somatic <- load_VarDict_vcf(path, merged_file)
    VarDict_somatic_snvs <-select_snvs(VarDict_somatic)
    fp_var = fp_snvs_VarDict(VarDict_somatic_snvs, pick_gt, gt_all)
    
    return(fp_var)
}

final_fn_snvs_VarDict<- function(path, merged_file, pick_gt){
    
    VarDict_somatic <- load_VarDict_vcf(path, merged_file)
    VarDict_somatic_snvs <-select_snvs(VarDict_somatic)
    fn_var = define_fn(VarDict_somatic_snvs, pick_gt)
    colnames(fn_var) = c("POS", "Ground Truth REF", "Ground Truth DP", 
                         "Ground Truth ALT", "Count", "Ground Truth AF", "mut", "type")
    
    return(fn_var)
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



