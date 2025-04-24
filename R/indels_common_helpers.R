
load_gt_report_indels <- function(path, merged_file) {
    #function to load Ground Truth bam-report 
    a <- paste0(path, "/", merged_file, "_report.tsv") |>
        readLines() |>
        str_split(pattern = "\t", simplify = TRUE) |>
        as.data.frame() |> 
        setDT()
    
    a$V1 = NULL
    a$V5 = NULL
    
    colnames(a) = c("POS", "REF", "DP", paste0("ALT_", 1:(ncol(a) - 3)))
    
    a = melt(
        a, id.vars = c("POS", "REF", "DP"),
        variable.factor = FALSE, value.factor = FALSE,
        variable.name = "ALT", value.name = "Count"
    )
    
    a = a[which(Count != "")]
    a$POS = as.numeric(a$POS)
    a$DP = as.numeric(a$DP)
    a$ALT = str_split_i(a$Count, "\\:", 1)
    a$Count = str_split_i(a$Count, "\\:", 2) |>
        as.numeric()
    a$Freq = round(a$Count / a$DP, digits = 6)
    a = a[order(POS, -Count)]
    a = a[which(REF != a$ALT & Count != 0)]
    
    # select indels
    a_indels = a[which(ALT %ni% c("A", "C", "G", "T")), ]
    #filter DEPTH>2
    #a_indels = a_indels[which(a_indels$Count >2), ]
    
    gt = list(
        all = a,
        indels = a_indels
        
    )
    return(gt)
}

select_indels <- function(df){
    #function to select indels from caller based on length of REF and ALT
    
    #identify indels based on length
    indels = df[nchar(df$REF) != nchar(df$ALT)]
    indels$mut = paste(indels$POS, indels$REF, indels$ALT, sep = ":")
    
    return(indels)
}

standardize_indels <- function(dt) {
    #function to standardize indels
    setDT(dt)
    
    #deletions
    dt[grepl("^-", ALT), `:=` (
        ALT = substring(REF, 1, 1), 
        REF = paste0(REF, substring(ALT, 2)),
        POS = POS - 1  #Adjust POS for deletions
    )]
    
    #insertions
    dt[grepl("^\\+", ALT), ALT := paste0(REF, substring(ALT, 2))]
    
    dt$mut = paste(dt$POS, 
                   dt$REF, 
                   dt$ALT, sep = ":")
    return(dt)
}

gt_stdz_indels <- function(path, merged_file){
    gt_indels = load_gt_report_indels(path, merged_file)$indels
    pick_gt_stdz = standardize_indels(gt_indels)
    colnames(pick_gt_stdz) = c("POS", "REF", "DP", "ALT", "AD", "Freq", "mut" )
    return(pick_gt_stdz)
} 

call_indels <- function(path, caller, merged_file, pick_gt_stdz) {
    
    if(caller == "Freebayes") {
        
        indels_var <- call_indels_Freebayes(path, merged_file, pick_gt_stdz)
        
    } else if (caller == "Mutect2") {
        
        indels_var <- call_indels_gatk(path, merged_file, pick_gt_stdz)
        
    } else if (caller == "LoFreq") {
        
        indels_var <- call_indels_LoFreq(path, merged_file, pick_gt_stdz)
        
    } else if (caller == "VarDict") {
        
        indels_var <- call_indels_VarDict(path, merged_file, pick_gt_stdz)
        
    } else if (caller == "VarScan") {
        
        indels_var <- call_indels_VarScan(path, merged_file, pick_gt_stdz)
    }
    
    return(indels_var)
}