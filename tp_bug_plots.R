source("R/libraries.R")

#df = fread("C:/Users/sfragkoul/Desktop/synth_data/coverage_test/300_30_10/Merged_Mutect2_snvs_TV.tsv", verbose = FALSE) 

gt_comparison <- "C:/Users/sfragkoul/Desktop/synth_data/coverage_test/300_30_10"
caller <- "Mutect2"
merged_file <- "Merged"

plot_snvs_TP <- function(gt_comparison, caller, merged_file) {
    
    # Construct file path
    file_path <- paste0(gt_comparison, "/", merged_file, "_", caller, "_snvs_Noise.tsv")
    
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
        return(ggplot() + labs(title = paste("No FP snvs data for", caller), x = NULL, y = NULL))
    }
    
    
    df_tp <- df[df$type == "TP", ]
    
    
    # Generate subplots if the file is not empty
    tp_plot1 <- tp_dp_barplot(df_tp, caller)
    tp_plot2 <- tp_af_barplot(df_tp, caller)
    
    # Combine the subplots
    tp_plot <- tp_plot1 + tp_plot2 +
        plot_layout(
            widths = c(1, 1)
        )
    
    return(fp_plot)
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
            labels = c(paste0(caller, " TP Variants"))
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
            labels = c(paste0(caller, " TP Variants"))
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

