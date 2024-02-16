


gt_analysis <- function(runs, folder) {
    
    nt_runs = list()
    
    for(r in runs) {
        
        a <- paste0(folder, "/", r, "/", r, "_report.tsv") |>
            readLines() |>
            str_split(pattern = "\t", simplify = TRUE) |>
            as.data.frame() |> 
            setDT()
        
        a$V1 = NULL
        # a$V3 = NULL
        a$V5 = NULL
        
        colnames(a) = c("POS", "REF", "DP", paste0("Nt_", 1:(ncol(a) - 3)))
        
        a = melt(
            a, id.vars = c("POS", "REF", "DP"),
            variable.factor = FALSE, value.factor = FALSE,
            variable.name = "Nt", value.name = "Count"
        )
        
        a = a[which(Count != "")]
        
        a$POS = as.numeric(a$POS)
        a$DP = as.numeric(a$DP)
        
        a$Nt = str_split_i(a$Count, "\\:", 1)
        
        a$Count = str_split_i(a$Count, "\\:", 2) |>
            as.numeric()
        
        a$Freq = round(100 * a$Count / a$DP, digits = 6)
        
        a = a[order(POS, -Count)]
        
        a = a[which(REF != a$Nt & Count != 0)]
        
        b = a[which(Nt %in% c("A", "C", "G", "T")), ]
        
        nt_runs[[ as.character(r) ]] = b
    }
    
    nt_runs = rbindlist(nt_runs, idcol = "Run")
    
    pos_of_interest = nt_runs[which(Freq == 100)]$POS |> unique()
    
    gt_runs = nt_runs[which(POS %in% pos_of_interest)]
    
    a <- paste0(folder, "/Merged_report.tsv") |> 
        readLines() |>
        str_split(pattern = "\t", simplify = TRUE) |> 
        as.data.frame() |> 
        setDT()
    
    a$V1 = NULL
    # a$V3 = NULL
    a$V5 = NULL
    
    colnames(a) = c("POS", "REF", "DP", paste0("Nt_", 1:(ncol(a) - 3)))
    
    a = melt(
        a, id.vars = c("POS", "REF", "DP"),
        variable.factor = FALSE, value.factor = FALSE,
        variable.name = "Nt", value.name = "Count"
    )
    
    a = a[which(Count != "")]
    
    a$POS = as.numeric(a$POS)
    a$DP = as.numeric(a$DP)
    
    a$Nt = str_split_i(a$Count, "\\:", 1)
    
    a$Count = str_split_i(a$Count, "\\:", 2) |>
        as.numeric()
    
    a$Freq = round(100 * a$Count / a$DP, digits = 6)
    
    a = a[order(POS, -Count)]
    
    a = a[which(REF != a$Nt & Count != 0)]
    
    b = a[which(Nt %in% c("A", "C", "G", "T")), ]
    
    
    merged_gt = b[which(POS %in% gt_runs$POS)]
    merged_gt = merged_gt[order(POS)]
    
    merged_gt$Freq = merged_gt$Freq / 100
    
    merged_gt = merged_gt[, by = .(POS, REF, DP), .(
        Nt = paste(Nt, collapse = ","),
        Count = paste(Count, collapse = ","),
        Freq = paste(round(Freq, digits = 3), collapse = ",")
    )]
    
    
    return(merged_gt)
    
}

merge_VarDict <- function(VarDict_somatic_vcf, merged_gt) {
    
    VarDict_s0  = VarDict_somatic_vcf |> vcfR::getFIX() |> as.data.frame() |> setDT()
    VarDict_s1  = VarDict_somatic_vcf |> extract_gt_tidy() |> setDT()
    VarDictgatk_s21 = VarDict_somatic_vcf |> extract_info_tidy() |> setDT()
    
    VarDict_somatic = cbind(VarDict_s0[VarDict_s1$Key, ], VarDict_s1)
    
    
    #Merge everything into a common file-------------------------------------------
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
    
    q[which(q$`VarDict ALT` == "")]$`VarDict ALT` = NA
    
    # plot 1 ------------------------
    
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
    
    
    # plot 2 ---------------------
    
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
    
    # return -------------
    
    return(
        list(
            "coverage" = o1,
            "allele" = o2
        )
    )
    
}

density_plot_VarDict <- function(q) {
    
    q[which(df$`VarDict ALT` == "")]$`VarDict ALT` = NA
    
    df = q[, c(
        "POS", 
        "Ground Truth AF",
        "Ground Truth ALT",
        "VarDict ALT",
        "VarDict AF"
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
    
    
    # return -----------------
    
    return(
        list(
            "groundtruth" = o1,
            "VarDict" = o2
        )
    )
    
}

bubble_plots_VarDict <- function(q) {
    
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

