
#TP SNVS-----------------------------------------------------------------------

plot_snvs_TP <- function(gt_snv_tp_comparison, vcf_path, gt_path, caller, merged_file) {
    
    df = fread(paste0(gt_snv_tp_comparison, "/", merged_file, "_", caller, "_snvs_TV.tsv"))
    
    vcf_GT <- read.vcfR(paste0(vcf_path, "/", merged_file, "_ground_truth_norm.vcf"), verbose = FALSE )
    
    vcf_caller <- read.vcfR(paste0(vcf_path, "/", merged_file, "_", caller, "_norm.vcf"), verbose = FALSE )
    
    if(caller == "Freebayes") {

        plots <- plot_snvs_TP_freebayes(df, vcf_GT, vcf_caller, merged_file)

    } else if (caller == "Mutect2") {

        plots <- plot_snvs_TP_gatk(df, vcf_GT, vcf_caller, merged_file)

    } else if (caller == "LoFreq") {

        plots <- plot_snvs_TP_LoFreq(df, vcf_GT, vcf_caller, merged_file)

    } else if (caller == "VarDict") {

        plots <- plot_snvs_TP_VarDict(df, vcf_GT, vcf_caller, merged_file)

    } else if (caller == "VarScan") {

        plots <- plot_snvs_TP_VarScan(df, vcf_GT, vcf_caller, merged_file)

    }
    
    return(plots)
    
}

#FP & FN SNVS------------------------------------------------------------------

fn_dp_barplot <- function(q, caller){
    #FP DP plot
    df = q[, c(
        "POS", 
        "Ground Truth DP"
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
                "Ground Truth DP" = color
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
                "Ground Truth AF" = color
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

plot_snvs_FP <- function(gt_comparison, caller, merged_file) {

    # Construct file path
    file_path <- paste0(gt_comparison, "/", merged_file, "_", caller, "_snvs_FP.tsv")

    # Check if file exists
    if (!file.exists(file_path)) {
        stop(paste("File does not exist:", file_path))
    }

    # Read the file
    df <- tryCatch(
        fread(file_path),
        error = function(e) {
            stop(paste("Error reading file:", file_path, "\n", e$message))
        }
    )

    # Check if the file is empty
    if (nrow(df) == 0) {
        warning(paste("File is empty:", file_path))
        # Return a placeholder plot or NULL
        return(ggplot() + labs(title = paste("No FP snvs data for", caller), x = NULL, y = NULL))
    }

    # Call specific plotting function based on the caller
    if (caller == "Freebayes") {
        fp_plot <- plot_snvs_FP_Freebayes(df, merged_file)
    } else if (caller == "Mutect2") {
        fp_plot <- plot_snvs_FP_gatk(df, merged_file)
    } else if (caller == "LoFreq") {
        fp_plot <- plot_snvs_FP_LoFreq(df, merged_file)
    } else if (caller == "VarDict") {
        fp_plot <- plot_snvs_FP_VarDict(df, merged_file)
    } else if (caller == "VarScan") {
        fp_plot <- plot_snvs_FP_VarScan(df, merged_file)
    } else {
        stop(paste("Unknown caller:", caller))
    }

    return(fp_plot)
}

plot_snvs_FN <- function(gt_comparison, caller, merged_file) {
    
    # Construct file path
    file_path <- paste0(gt_comparison, "/", merged_file, "_", caller, "_snvs_FN.tsv")
    
    # Check if file exists
    if (!file.exists(file_path)) {
        stop(paste("File does not exist:", file_path))
    }
    
    # Read the file
    df <- fread(file_path)
    
    # Check if the file is empty
    if (nrow(df) == 0) {
        warning(paste("File is empty:", file_path))
        # Return a placeholder plot or NULL
        return(ggplot() + labs(title = paste("No FN snvs data for", caller), x = NULL, y = NULL))
    }
    
    # Generate subplots if the file is not empty
    fn_plot1 <- fn_dp_barplot(df, caller)
    fn_plot2 <- fn_af_barplot(df, caller)
    
    # Combine the subplots
    fn_plot <- fn_plot1 + fn_plot2 +
        plot_layout(
            widths = c(1, 1)
        )
    
    return(fn_plot)
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