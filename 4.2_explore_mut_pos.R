#'
#'A script to report variants that were either  detected or not detected
#'by each caller. The read positions were divided in bins to study possible
#'correlated trends.
#'
#'Input: tsv file with reported variants from caller, reports in tsv format 
#'with the variants in each chromosomal position for each individual 
#'"golden" bam file 
#'
#'Output: reports with variants that were either detected or not detected
#'by each caller. The read positions were divided in bins.
#'
#'Author: Nikos Pechlivanis(github:npechl) 
#'
#'Please replace the name of the caller for each use case. 

rm(list = ls())
gc()

# load libraries --------------------------------------------------------------

library(data.table)
library(stringr)
library(ggplot2)
library(ggforce)


files = as.character(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
folder = "coverage_test/3000_300_10/"
name = "3000_300_10"

for (file in files){
    #file= "1"
    message(paste0("working on file ", file))
    # input parameters --------------------------------------------------------
    qv = paste0(folder, "Ground_truth_vs_VarScan.clean_norm.tsv")
    gb = paste0(folder, file, "/", file, "_position_report2.csv")
    
    # golden alignments -------------------------------------------------------
    gb = gb |> fread()
    gb$key = paste(gb$pos_chromosomal, gb$REF, gb$ALT, sep = ":")
    
    # VarScan missing variants ----------------------
    
    qv = qv |> fread()
    qv$key = paste(qv$POS, qv$`Ground Truth REF`, qv$`Ground Truth ALT`, sep = ":")
    
    
    # filtering common mutations ----------------------------
    
    gb_filtered = gb[which(key %in% qv$key)]
    qv_filtered = qv[which(key %in% gb_filtered$key)]
    
    q1 = qv_filtered[which(is.na(`VarScan AF`))]
    q2 = gb_filtered[which(key %in% q1$key)]
    
    q2$bin = floor(q2$pos_read / 10) + 1
    q2 = q2[order(q2$bin), ]
    
    fwrite(
    q2, paste0(folder, file, "/", file, "_VarScan_read_pos_bins_not_found.tsv"),
    row.names = FALSE, quote = FALSE, sep = "\t"
    )
    
    #boxplot-------------------------------------------------------------------
    
    #q2$bin = paste0("bin", q2$bin)
    
    #q2$bin = q2$bin |> factor(levels = unique(q2$bin))
    
    #gr1 = ggplot(data = q2) +
        
        # geom_point(aes(x = `bin`, y = `No. of sequences`, color = `bin`))
        
     #   geom_violin(aes(x = bin, y = `No. of sequences`, fill = bin)) + 
        
     #   ggtitle("Not Found")
    
    
    #ggsave(
    #    plot = gr1, filename = paste0(folder, file, "/", file, "_Box_plots_not_found.jpeg"),
    #    width = 8, height = 8, units = "in", dpi = 600
    #)    
    
    #FOUND---------------------------------------------------------------------
    
    q1 = qv_filtered[which(!is.na(`VarScan AF`))]
    q2 = gb_filtered[which(key %in% q1$key)] 
    q2$bin = floor(q2$pos_read / 10) + 1
    
    q2 = q2[order(q2$bin), ]
    
    fwrite(
        q2, paste0(folder, file, "/", file, "_VarScan_read_pos_bins_found.tsv"),
        row.names = FALSE, quote = FALSE, sep = "\t"
    )
    
    
    
    #q2$bin = paste0("bin", q2$bin)
    
    #q2$bin = q2$bin |> factor(levels = unique(q2$bin))
    
    #gr1 = ggplot(data = q2) +
        
     #   geom_violin(aes(x = bin, y = `No. of sequences`, fill = bin)) + 
        
      #  ggtitle("Found")
    
    
   # ggsave(
    #    plot = gr1, filename = paste0(folder, file, "/", file, "_Box_plots_found.jpeg"),
   #     width = 8, height = 8, units = "in", dpi = 600
   # )  
}


# bind all reports ------------------------------------------------------------

nt_runs = list()

for(r in files) {
  #r=1
  gb = paste0(folder, r, "/", r, "_VarScan_read_pos_bins_found.tsv")
  gb = gb |> fread()
  nt_runs[[ as.character(r) ]] = gb
  
}

nt_runs = rbindlist(nt_runs, idcol = "Run")
nt_runs$Found = "Yes"

#-----------------------------------------------------------------------------

nt_runs2 = list()

for(g in files) {
  
  gb2 = paste0(folder, g, "/", g, "_VarScan_read_pos_bins_not_found.tsv")
  gb2 = gb2 |> fread()
  nt_runs2[[ as.character(g) ]] = gb2
  
}

nt_runs2 = rbindlist(nt_runs2, idcol = "Run")
nt_runs2$Found = "No"


final = rbind(nt_runs, nt_runs2)
final$names = NULL

fwrite(
  final, paste0(folder, name,"_All_VarScan_read_pos_bins_report", ".csv"),
  row.names = FALSE, quote = TRUE, sep = ","
)


# fwrite(
#     nt_runs, paste0(folder, "All_found_report2_", name, ".csv"),
#     row.names = FALSE, quote = TRUE, sep = ","
# )
# 
# fwrite(
#     nt_runs2, paste0(folder, "All_not_found_report2_", name, ".csv"),
#     row.names = FALSE, quote = TRUE, sep = ","
# )

