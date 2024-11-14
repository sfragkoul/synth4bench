
# Libraries
source("R/libraries.R")

circular_plot_gatk <- function(path, merged_file, caller){
    #Load data
    tp = fread(paste0(path, "/", merged_file, "_", caller, "_indels_TP.tsv"), sep = "\t")
    fp = fread(paste0(path, "/", merged_file, "_", caller, "_indels_FP.tsv"), sep = "\t")
    fn = fread(paste0(path, "/", merged_file, "_", caller, "_indels_FN.tsv"), sep = "\t")
    
    tp = tp[, .(POS, REF, ALT, type)]
    tp$REF_len <- str_length(tp$REF)
    tp$ALT_len <- str_length(tp$ALT)
    tp$len_dif <- tp$ALT_len - tp$REF_len
    tp$category <- "not exist"
    
    fp = fp[, .(POS, REF, ALT, type, category)]
    fp$REF_len <- str_length(fp$REF)
    fp$ALT_len <- str_length(fp$ALT)
    fp$len_dif <- fp$ALT_len - fp$REF_len
    
    fn = fn[, .(POS, REF, ALT, type, category)]
    fn$REF_len <- str_length(fn$REF)
    fn$ALT_len <- str_length(fn$ALT)
    fn$len_dif <- fn$ALT_len - fn$REF_len
    
    #Combine the datasets 
    data = rbind(tp, fp)
    df = rbind(data, fn)
    colnames(df) <- c("POS", "REF", "ALT", "Type",  "REF_len", "ALT_len", "len_dif", "Category")
    
    #plot ------------------------------------------------------------------------
    #Adjust data so that each type has its own y-offset
    df <- df |>
        mutate(y_cycle = case_when(
            Type == "FN" ~ len_dif + 50,   # Shift FN cycle outward
            Type == "FP" ~ len_dif + 25,   # Shift FP cycle to middle
            Type == "TP" ~ len_dif         # Keep TP at the center
        ))
    
    #Ensure 'category' is a factor
    df$Category <- factor(df$Category, levels = c("not exist", "diff REF", "diff ALT"))
    
    p <- ggplot(df, aes(x = POS, y = y_cycle)) +
        
        #Lollipop segments: start each from the respective baseline to the point
        geom_segment(
            aes(x = POS, xend = POS,
                y = ifelse(Type == "FN", 50, ifelse(Type == "FP", 25, 0)),
                yend = y_cycle),
            color = "grey75", linewidth = 0.25, lineend = "round"
        ) +
        
        #Dashed lines for separation of each cycle level
        geom_hline(yintercept = 50, color = "grey40") +
        geom_hline(yintercept = 25, color = "grey40") +
        geom_hline(yintercept = 0,  color = "grey40") +
    
    
        # Add points at the end of each segment for the lollipop head
        geom_point(aes(fill = Type, color = Type, shape = Category, 
                       size = ifelse(Category == "not exist", 1.5, 3)), # Increase size for specific categories
                   stroke = .15) +
        
        scale_size_identity() +
        
        #Define specific shapes for each category level
        scale_shape_manual(values = c("diff REF" = 23, "diff ALT" = 24, "not exist" = 21)) +
        
        #Define custom colors for each type
        scale_fill_manual(values = c("TP" = "#a78d95", "FP" = "#ae4364", "FN" = "#43ae8d")) +
        scale_color_manual(values = c("TP" = "#a78d95", "FP" = "#ae4364", "FN" = "#43ae8d") |> darken(.25)) +
    
        #Customize the x-axis and radial coordinates
        scale_x_continuous(breaks = c(0, 4751, 9503, 14255, 19007), limits = c(0, 19007)) +
        
        coord_radial(start = pi / 2.5, inner.radius = .25, end = 2.6 * pi) +
    
        #Remove legend for size if unnecessary
        guides(size = "none") +
        
        #Define minimal theme and other plot aesthetics
        theme_minimal() +
        theme(
            axis.text.y = element_blank(),
            panel.grid.major = element_line(linewidth = 0.35),
            panel.grid.minor = element_blank(),
            plot.margin = margin(20, 20, 20, 20),
            plot.title = element_text(size = 15, face = "bold", hjust = 0.5)
        ) +
        
        #Add labels for the plot
        labs(
            title = "Ground Truth vs Mutect2 INDELS",
            y = "REF vs ALT Length Difference",
            x = "Chromosomal Position"
            # color = "Type"
        )

    return(p)
}

p = circular_plot_gatk("results", "Merged", "Mutect2")


ggsave(
    plot = p, filename = "Rplot_2.jpeg",
    width = 14, height = 12, units = "in", dpi = 600
)
