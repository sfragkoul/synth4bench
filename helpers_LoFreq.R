


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

merge_LoFreq <- function(LoFreq_somatic_vcf, merged_gt) {
    LoFreq_s0  = LoFreq_somatic_vcf |> vcfR::getFIX() |> as.data.frame() |> setDT()
    #LoFreq_s1  = LoFreq_somatic_vcf |> extract_gt_tidy() |> setDT()
    LoFreq_s2 = LoFreq_somatic_vcf |> extract_info_tidy() |> setDT()
    LoFreq_s2 = LoFreq_s2[,c( "DP", "AF" )]
    
    LoFreq_somatic = cbind(LoFreq_s0, LoFreq_s2)
    
    #Merge everything into a common file-------------------------------------------
    merged_gt$POS = as.character(merged_gt$POS)
    
    merged_bnch = merge(merged_gt, LoFreq_somatic,  by = "POS", all.x = TRUE)
    
    merged_bnch$POS = as.numeric(merged_bnch$POS)
    
    merged_bnch = merged_bnch[order(POS)]
    
    colnames(merged_bnch) = c(
        "POS",	"Ground Truth REF",	"Ground Truth DP",
        "Ground Truth ALT", "Ground Truth AD", 
        "Ground Truth AF", "CHROM", "ID",	"LoFreq REF",	
        "LoFreq ALT", "LoFreq QUAL", "LoFreq FILTER", 
        "LoFreq DP", "LoFreq AF"
    )
    
    return(merged_bnch)
    
}

clean_LoFreq <- function(df) {
    
    df2 = df[, c(
        "POS",
        
        "Ground Truth REF",
        "Ground Truth ALT",
        "Ground Truth DP",
        "Ground Truth AF",
        
        "LoFreq REF",
        "LoFreq ALT",
        "LoFreq DP",
        "LoFreq AF"
    ), with = FALSE]
    
    
    
    df2 = df2[, by = c(
        "POS",
        "Ground Truth REF",
        "Ground Truth DP",
        "LoFreq REF",
        "LoFreq ALT",
        "LoFreq DP",
        "LoFreq AF"

    ), .(
        "Ground Truth ALT" = `Ground Truth ALT` |> tstrsplit(",") |> unlist(),
        "Ground Truth AF"  = `Ground Truth AF` |> tstrsplit(",") |> unlist()
        # "LoFreq REF" = `LoFreq REF` |> tstrsplit(",") |> unlist(),
        # "LoFreq ALT" = `LoFreq ALT` |> tstrsplit(",") |> unlist(),
        # "LoFreq DP"  = `LoFreq DP` |> tstrsplit(",") |> unlist() |> as.integer(),
        # "LoFreq AF"  = `LoFreq AF` |> tstrsplit(",") |> unlist() |> as.numeric()
    )]



    LoFreq_alt = str_split(df2$`LoFreq ALT`, ",")
    LoFreq_af = str_split(df2$`LoFreq AF`, ",")


    cln = mapply(
        function(x, y, z) {

            index = which(y == x)

            return(
                c(y[index], z[index])
            )

        },

        df2$`Ground Truth ALT`, LoFreq_alt, LoFreq_af
    )


    df2$`LoFreq ALT` = cln |> lapply(function(x) { return(x [1]) }) |> unlist()
    df2$`LoFreq AF` = cln |> lapply(function(x) { return(x [2]) }) |> unlist()

    df2[which(is.na(`LoFreq AF`))]$`LoFreq DP` = NA
    df2[which(is.na(`LoFreq AF`))]$`LoFreq REF` = NA

    df2 = df2[, c(
        "POS",
        "Ground Truth REF",
        "Ground Truth ALT",
        "Ground Truth DP",
        "Ground Truth AF",
        "LoFreq REF",
        "LoFreq ALT",
        "LoFreq DP",
        "LoFreq AF"
    ), with = FALSE]
    
    return(df2)
    
}



bar_plots_LoFreq <- function(q) {
    
    q[which(q$`LoFreq ALT` == "")]$`LoFreq ALT` = NA
    
    # plot 1 ------------------------
    
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
        
        geom_boxplot(aes(x = variable, y = value, fill = variable),
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
    
    
    # plot 2 ---------------------
    
    df = q[, c(
        "POS",
        "Ground Truth AF",
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
                "Ground Truth AF" = "#43ae8d",
                "LoFreq AF"      = "#c974ba"
            )
        ) +
        
        scale_x_discrete(
            breaks = c("Ground Truth AF", "LoFreq AF"),
            labels = c("Ground Truth", "LoFreq")
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

density_plot_LoFreq <- function(q) {
    
    q[which(df$`LoFreq ALT` == "")]$`LoFreq ALT` = NA
    
    df = q[, c(
        "POS", 
        "Ground Truth AF",
        "Ground Truth ALT",
        "LoFreq ALT",
        "LoFreq AF"
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
    
    o2 = ggplot(data = df[which(!is.na(`LoFreq ALT`)), c(1, 4, 5)], aes(x = `LoFreq AF`)) +
        
        geom_density(aes(color = `LoFreq ALT`, fill = `LoFreq ALT`),
                     alpha = .5) +
        
        scale_x_continuous(expand = c(0, 0), breaks = c(.25, .5, .75, 1), limits = c(0, 1), labels = scales::percent) +
        scale_y_continuous(expand = c(0, 0), limits = c(0, 15), breaks = c(5, 10, 15))+
        
        scale_fill_npg() +
        scale_color_npg() +
        
        facet_wrap(vars(`LoFreq ALT`), nrow = 1) +
        
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
            y = "LoFreq (density)"
        )
    
    
    # return -----------------
    
    return(
        list(
            "groundtruth" = o1,
            "LoFreq" = o2
        )
    )
    
}

bubble_plots_LoFreq <- function(q) {
    
    # q[which(q$`LoFreq ALT` == "")]$`LoFreq ALT` = NA
    
    
    q1 = q[which(q$`LoFreq ALT` != "")]
    
    
    o = ggplot() +
        
        geom_point(
            data = q,
            aes(x = POS, y = `Ground Truth ALT`, size = `Ground Truth AF`, fill = "Ground Truth"),
            position = position_nudge(y = .15), shape = 21, stroke = .25
        ) +
        
        geom_point(
            data = q1,
            aes(x = POS, y = `LoFreq ALT`, size = `LoFreq AF`, fill = "LoFreq"),
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
                "LoFreq"      = "#c974ba"
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

venn_plot_LoFreq <- function(q, p) {
    
    vcf_GT = vcfR::getFIX(q) |> as.data.frame() |> setDT()
    vcf_GT$scenario = "GT"
    
    vcf_LoFreq = vcfR::getFIX(p) |> as.data.frame() |> setDT()
    vcf_LoFreq$scenario = "LoFreq"
    
    x = rbind(vcf_GT, vcf_LoFreq)
    y = x[, c("CHROM", "POS", "REF", "ALT", "scenario"), with = FALSE]
    
    y$mut = paste(y$CHROM, y$POS, y$REF, y$ALT, sep = ":")
    
    y = split(y, y$scenario)
    
    y = list(
        'Ground Truth' = y$GT$mut,
        'LoFreq'         = y$LoFreq$mut
    )
    
    gr = ggvenn(y, fill_color = c("#43ae8d", "#c974ba")) +
        
        coord_equal(clip = "off")
    
    return(gr)
}
