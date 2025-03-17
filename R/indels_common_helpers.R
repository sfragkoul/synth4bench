
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
    a_indels = a_indels[which(a_indels$Count >2), ]
    
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

load_gt_vcf_indels <- function(path, merged_file, gt_indels){
    #function to load Ground Truth vcf
    ground_truth_vcf <- read.vcfR( paste0(path, "/",merged_file, 
                                          "_ground_truth_norm.vcf"),
                                   verbose = FALSE )
    
    ground_truth_vcf  = ground_truth_vcf |> vcfR::getFIX() |> as.data.frame() |> setDT()
    
    pick_gt = gt_indels[which(gt_indels$POS %in% ground_truth_vcf$POS)]
    pick_gt$mut = paste(pick_gt$POS, 
                        pick_gt$REF, 
                        pick_gt$ALT, sep = ":")
    return(pick_gt)
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
    gt_all = load_gt_report_indels(path, merged_file)$all
    gt_indels = load_gt_report_indels(path, merged_file)$indels
    pick_gt = load_gt_vcf_indels(path, merged_file, gt_indels)
    pick_gt_stdz = standardize_indels(pick_gt)
    return(pick_gt_stdz)
} 


call_tp_indels <- function(path, caller, merged_file, pick_gt_stdz) {
  
  if(caller == "Freebayes") {
    
    tp_indels <- final_tp_indels_Freebayes(path, merged_file, pick_gt_stdz)
    
  } else if (caller == "Mutect2") {
    
    tp_indels <- final_tp_indels_gatk(path, merged_file, pick_gt_stdz)
    
  } else if (caller == "LoFreq") {
    
    tp_indels <- final_tp_indels_LoFreq(path, merged_file, pick_gt_stdz)
    
  } else if (caller == "VarDict") {
    
    tp_indels <- final_tp_indels_VarDict(path, merged_file, pick_gt_stdz)
    
  } else if (caller == "VarScan") {
    
    tp_indels <- final_tp_indels_VarScan(path, merged_file, pick_gt_stdz)
    
  }
  
  return(tp_indels)
}



call_fn_indels <- function(path, caller, merged_file, pick_gt_stdz) {
  
  if(caller == "Freebayes") {
    
    fn_indels <- call_fn_indels_Freebayes(path, merged_file, pick_gt_stdz)
    
  } else if (caller == "Mutect2") {
    
    fn_indels <- call_fn_indels_gatk(path, merged_file, pick_gt_stdz)
    
  } else if (caller == "LoFreq") {
    
    fn_indels <- call_fn_indels_LoFreq(path, merged_file, pick_gt_stdz)
    
  } else if (caller == "VarDict") {
    
    fn_indels <- call_fn_indels_VarDict(path, merged_file, pick_gt_stdz)
    
  } else if (caller == "VarScan") {
    
    fn_indels <- call_fn_indels_VarScan(path, merged_file, pick_gt_stdz)
    
  }
  
  return(fn_indels)
}


call_fp_indels <- function(path, caller, merged_file, pick_gt_stdz) {
  
  if(caller == "Freebayes") {
    
    fp_indels <- call_fp_indels_Freebayes(path, merged_file, pick_gt_stdz)
    
  } else if (caller == "Mutect2") {
    
    fp_indels <- call_fp_indels_gatk(path, merged_file, pick_gt_stdz)
    
  } else if (caller == "LoFreq") {
    
    fp_indels <- call_fp_indels_LoFreq(path, merged_file, pick_gt_stdz)
    
  } else if (caller == "VarDict") {
    
    fp_indels <- call_fp_indels_VarDict(path, merged_file, pick_gt_stdz)
    
  } else if (caller == "VarScan") {
    
    fp_indels <- call_fp_indels_VarScan(path, merged_file, pick_gt_stdz)
    
  }
  
  return(fp_indels)
}




