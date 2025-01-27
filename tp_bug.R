source("R/libraries.R")

folder = "D:/sfragkoul/Synth_Data/Synthesizers/NEAT/testing/TP53/coverage_test/300_30_10"
runs = c(1,2)
#runs = c(1,2,3,4,5,6,7,8,9,10)
merged_file = "Merged"

gt_analysis <- function(runs, folder, merged_file) {
    
    nt_runs = list()
    #ground truth   variants from individual files
    for(r in runs) {
        #folder = "."
        #merged_file = "Merged"
        #process reports.tsv files for individual files
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
    
    gt_runs = nt_runs[POS %in% pos_of_interest & Freq == "100"] #!!!NEW
    
    #same process reports.tsv files for Merged file
    a <- paste0(folder, "/", merged_file , "_report.tsv") |> 
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
    
    
    #merged_gt = b[which(POS %in% gt_runs$POS)]
    merged_gt <- merge(b, gt_runs, by = c("POS", "REF", "Nt")) #ΝEW!!!!!!!
    colnames(merged_gt) = c("POS", "REF", "DP", "Nt", "Count", "Freq",
                             "Run", "DP Indiv", "Count Indiv", "Freq Indiv")
    
    merged_gt = merged_gt[order(POS)]
    
    merged_gt$Freq = merged_gt$Freq / 100
    
    # merged_gt1 = merged_gt[, by = .(POS, REF, DP), .(
    #     Nt = paste(Nt, collapse = ","),
    #     Count = paste(Count, collapse = ","),
    #     Freq = paste(round(Freq, digits = 3), collapse = ",")
    # )]
    
    return(merged_gt)
    
}