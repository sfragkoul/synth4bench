
#TP SNVS-----------------------------------------------------------------------

plot_snvs_TP <- function(gt_snv_tp_comparison, vcf_path, gt_path, caller, merged_file) {
    
    df = fread(paste0(gt_snv_tp_comparison, "/", merged_file, "_", caller, "_snvs_TV.tsv"))
    
    #vcf_GT <- read.vcfR(paste0(vcf_path, "/", merged_file, "_ground_truth_norm.vcf"), verbose = FALSE )
    #vcf_caller <- read.vcfR(paste0(vcf_path, "/", merged_file, "_", caller, "_norm.vcf"), verbose = FALSE )
    
    if(caller == "Freebayes") {

        plots <- plot_snvs_TP_freebayes(df, merged_file)

    } else if (caller == "Mutect2") {

        plots <- plot_snvs_TP_gatk(df, merged_file)

    } else if (caller == "LoFreq") {

        plots <- plot_snvs_TP_LoFreq(df, merged_file)

    } else if (caller == "VarDict") {

        plots <- plot_snvs_TP_VarDict(df, merged_file)

    } else if (caller == "VarScan") {

        plots <- plot_snvs_TP_VarScan(df, merged_file)

    }
    
    return(plots)
    
}

#Noise-------------------------------------------------------------------------
#FN
fn_dp_barplot <- function(q, caller){
    #FP DP plot
    df = q[, c(
        "POS", 
        "AD"
    ), with = FALSE] |>
        unique() |>
        melt(id.vars = "POS", variable.factor = FALSE, value.factor = FALSE)
    
    #set color
    if(caller == "Freebayes") {
        
        color <- "#ae8d43"
        
    } else if (caller == "Mutect2") {
        
        color <- "#ae4364"
        
    } else if (caller == "LoFreq") {
        
        color <- "#c974ba"
        
    } else if (caller == "VarDict") {
        
        color <- "#8d43ae"
        
    } else if (caller == "VarScan") {
        
        color <- "#439aae"
    }
    
    
    o3=ggplot(data = df) +
        
        geom_point(aes(x = variable, y = value, fill = variable),
                   position = position_jitternormal(sd_x = .01, sd_y = 0),
                   shape = 21, stroke = .1, size = 2.5) +
        
        geom_boxplot(aes(x = variable, y = value, fill = variable),
                     width = .25, alpha = .5, outlier.shape = NA) +
        
        scale_fill_manual(
            values = c(
                "AD" = color
            )
        ) +
        
        scale_x_discrete(
            labels = c(paste0(caller, " Noise FN Variants"))
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
            y = "Allele Depth (No. of reads)"
        )
    return(o3)
    
}
fn_af_barplot <- function(q, caller){
    #FP AF plot
    df = q[, c(
        "POS",
        "AF"
    ), with = FALSE] |>
        unique() |>
        
        melt(id.vars = "POS", variable.factor = FALSE, value.factor = FALSE)
    
    #set color
    if(caller == "Freebayes") {
        
        color <- "#ae8d43"
        
    } else if (caller == "Mutect2") {
        
        color <- "#ae4364"
        
    } else if (caller == "LoFreq") {
        
        color <- "#c974ba"
        
    } else if (caller == "VarDict") {
        
        color <- "#8d43ae"
        
    } else if (caller == "VarScan") {
        
        color <- "#439aae"
    }
    
    o4 = ggplot(data = df[which(!is.na(value) & value != 0)]) +
        
        geom_point(aes(x = variable, y = value, fill = variable),
                   position = position_jitternormal(sd_x = .01, sd_y = 0),
                   shape = 21, stroke = .1, size = 2.5) +
        
        geom_boxplot(aes(x = variable, y = value, fill = variable),
                     width = .25, alpha = .5, outlier.shape = NA) +
        
        scale_fill_manual(
            values = c(
                "AF" = color
            )
        ) +
        
        scale_x_discrete(
            labels = c(paste0(caller, " Noise FN Variants"))
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
plot_snvs_FN <- function(gt_comparison, caller, merged_file) {
    
    # Construct file path
    file_path <- paste0(gt_comparison, "/", merged_file, "_", caller, "_snvs_Noise.tsv")
    
    # Check if file exists
    if (!file.exists(file_path)) {
        stop(paste("File does not exist:", file_path))
    }
    
    # Read the file
    df <- fread(file_path)
    
    
    df_fn <- df[df$type == "FN", ]
    
    # Check if any false positives are present
    if (nrow(df_fn) == 0) {
        warning(paste("No FN SNVs data for", caller))
        # Return a placeholder plot indicating no FP SNVs
        return(
            ggplot() +
                labs(
                    title = paste("No FN SNVs for", caller),
                    x = NULL,
                    y = NULL
                )
        )
    }
    
    
    # Generate subplots if the file is not empty
    fn_plot1 <- fn_dp_barplot(df_fn, caller)
    fn_plot2 <- fn_af_barplot(df_fn, caller)
    
    # Combine the subplots
    fn_plot <- fn_plot1 + fn_plot2 +
        plot_layout(
            widths = c(1, 1)
        )
    
    return(fn_plot)
}

#FP
fp_dp_barplot <- function(q, caller){
    #FP DP plot
    df = q[, c(
        "POS", 
        "AD"
    ), with = FALSE] |>
        unique() |>
        melt(id.vars = "POS", variable.factor = FALSE, value.factor = FALSE)
    
    #set color
    if(caller == "Freebayes") {
        
        color <- "#ae8d43"
        
    } else if (caller == "Mutect2") {
        
        color <- "#ae4364"
        
    } else if (caller == "LoFreq") {
        
        color <- "#c974ba"
        
    } else if (caller == "VarDict") {
        
        color <- "#8d43ae"
        
    } else if (caller == "VarScan") {
        
        color <- "#439aae"
    }
    
    
    o3=ggplot(data = df) +
        
        geom_point(aes(x = variable, y = value, fill = variable),
                   position = position_jitternormal(sd_x = .01, sd_y = 0),
                   shape = 21, stroke = .1, size = 2.5) +
        
        geom_boxplot(aes(x = variable, y = value, fill = variable),
                     width = .25, alpha = .5, outlier.shape = NA) +
        
        scale_fill_manual(
            values = c(
                "AD" = color
            )
        ) +
        
        scale_x_discrete(
            labels = c(paste0(caller, " Noise FP Variants"))
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
            y = "Allele Depth (No. of reads)"
        )
    return(o3)
    
}
fp_af_barplot <- function(q, caller){
    #FP AF plot
    df = q[, c(
        "POS",
        "AF"
    ), with = FALSE] |>
        unique() |>
        
        melt(id.vars = "POS", variable.factor = FALSE, value.factor = FALSE)
    
    #set color
    if(caller == "Freebayes") {
        
        color <- "#ae8d43"
        
    } else if (caller == "Mutect2") {
        
        color <- "#ae4364"
        
    } else if (caller == "LoFreq") {
        
        color <- "#c974ba"
        
    } else if (caller == "VarDict") {
        
        color <- "#8d43ae"
        
    } else if (caller == "VarScan") {
        
        color <- "#439aae"
    }
    
    o4 = ggplot(data = df[which(!is.na(value) & value != 0)]) +
        
        geom_point(aes(x = variable, y = value, fill = variable),
                   position = position_jitternormal(sd_x = .01, sd_y = 0),
                   shape = 21, stroke = .1, size = 2.5) +
        
        geom_boxplot(aes(x = variable, y = value, fill = variable),
                     width = .25, alpha = .5, outlier.shape = NA) +
        
        scale_fill_manual(
            values = c(
                "AF" = color
            )
        ) +
        
        scale_x_discrete(
            labels = c(paste0(caller, " Noise FP Variants"))
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
plot_snvs_FP <- function(gt_comparison, caller, merged_file) {
    
    # Construct file path
    file_path <- paste0(gt_comparison, "/", merged_file, "_", caller, "_snvs_Noise.tsv")
    
    # Check if file exists
    if (!file.exists(file_path)) {
        stop(paste("File does not exist:", file_path))
    }
    
    # Read the file
    df <- fread(file_path)
    
    
    df_fp <- df[df$type == "FP", ]
    
    # Check if any false positives are present
    if (nrow(df_fp) == 0) {
        warning(paste("No FP SNVs data for", caller))
        # Return a placeholder plot indicating no FP SNVs
        return(
            ggplot() +
                labs(
                    title = paste("No FP SNVs for", caller),
                    x = NULL,
                    y = NULL
                )
        )
    }
    
    
    # Generate subplots if the file is not empty
    fp_plot1 <- fp_dp_barplot(df_fp, caller)
    fp_plot2 <- fp_af_barplot(df_fp, caller)
    
    # Combine the subplots
    fp_plot <- fp_plot1 + fp_plot2 +
        plot_layout(
            widths = c(1, 1)
        )
    
    return(fp_plot)
}

#TP
plot_snvs_TP <- function(gt_comparison, caller, merged_file) {
    
    # Construct file path
    file_path <- paste0(gt_comparison, "/", merged_file, "_", caller, "_snvs_Noise.tsv")
    
    # Check if file exists
    if (!file.exists(file_path)) {
        stop(paste("File does not exist:", file_path))
    }
    
    # Read the file
    df <- fread(file_path)
    df_tp <- df[df$type == "TP", ]
    
    # Check if any false positives are present
    if (nrow(df_tp) == 0) {
        warning(paste("No TP SNVs data for", caller))
        # Return a placeholder plot indicating no FP SNVs
        return(
            ggplot() +
                labs(
                    title = paste("No TP SNVs for", caller),
                    x = NULL,
                    y = NULL
                )
        )
    }
    
    # Generate subplots if the file is not empty
    tp_plot1 <- tp_dp_barplot(df_tp, caller)
    tp_plot2 <- tp_af_barplot(df_tp, caller)
    
    # Combine the subplots
    tp_plot <- tp_plot1 + tp_plot2 +
        plot_layout(
            widths = c(1, 1)
        )
    
    return(tp_plot)
}
tp_dp_barplot <- function(q, caller){
    #FP DP plot
    df = q[, c(
        "POS", 
        "AD"
    ), with = FALSE] |>
        unique() |>
        melt(id.vars = "POS", variable.factor = FALSE, value.factor = FALSE)
    
    #set color
    if(caller == "Freebayes") {
        
        color <- "#ae8d43"
        
    } else if (caller == "Mutect2") {
        
        color <- "#ae4364"
        
    } else if (caller == "LoFreq") {
        
        color <- "#c974ba"
        
    } else if (caller == "VarDict") {
        
        color <- "#8d43ae"
        
    } else if (caller == "VarScan") {
        
        color <- "#439aae"
    }
    
    
    o3=ggplot(data = df) +
        
        geom_point(aes(x = variable, y = value, fill = variable),
                   position = position_jitternormal(sd_x = .01, sd_y = 0),
                   shape = 21, stroke = .1, size = 2.5) +
        
        geom_boxplot(aes(x = variable, y = value, fill = variable),
                     width = .25, alpha = .5, outlier.shape = NA) +
        
        scale_fill_manual(
            values = c(
                "AD" = color
            )
        ) +
        
        scale_x_discrete(
            labels = c(paste0(caller, " Noise TP Variants"))
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
            y = "Allele Depth (No. of reads)"
        )
    return(o3)
    
}
tp_af_barplot <- function(q, caller){
    #FP AF plot
    df = q[, c(
        "POS",
        "AF"
    ), with = FALSE] |>
        unique() |>
        
        melt(id.vars = "POS", variable.factor = FALSE, value.factor = FALSE)
    
    #set color
    if(caller == "Freebayes") {
        
        color <- "#ae8d43"
        
    } else if (caller == "Mutect2") {
        
        color <- "#ae4364"
        
    } else if (caller == "LoFreq") {
        
        color <- "#c974ba"
        
    } else if (caller == "VarDict") {
        
        color <- "#8d43ae"
        
    } else if (caller == "VarScan") {
        
        color <- "#439aae"
    }
    
    o4 = ggplot(data = df[which(!is.na(value) & value != 0)]) +
        
        geom_point(aes(x = variable, y = value, fill = variable),
                   position = position_jitternormal(sd_x = .01, sd_y = 0),
                   shape = 21, stroke = .1, size = 2.5) +
        
        geom_boxplot(aes(x = variable, y = value, fill = variable),
                     width = .25, alpha = .5, outlier.shape = NA) +
        
        scale_fill_manual(
            values = c(
                "AF" = color
            )
        ) +
        
        scale_x_discrete(
            labels = c(paste0(caller, " Noise TP Variants"))
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

#INDELs------------------------------------------------------------------------
plot_indels <- function(path, merged_file, caller) {
    
    if(caller == "Freebayes") {
        
        plots <- circular_plot_Freebayes(path, merged_file, caller)
        
    } else if (caller == "Mutect2") {
        
        plots <- circular_plot_gatk(path, merged_file, caller)
        
    } else if (caller == "LoFreq") {
        
        plots <- circular_plot_LoFreq(path, merged_file, caller)
        
    } else if (caller == "VarDict") {
        
        plots <- circular_plot_VarDict(path, merged_file, caller)
        
    } else if (caller == "VarScan") {
        
        plots <- circular_plot_VarScan(path, merged_file, caller)
        
    }
    
    return(plots)
    
}