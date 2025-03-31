
#TP ---------------------------------------------------------------------------
gt_analysis <- function(runs, folder, merged_file) {
  
  nt_runs = list()
  
  for(r in runs) {
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
  
  gt_runs = nt_runs[POS %in% pos_of_interest & Freq == "100"]
  
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
  merged_gt <- merge(b, gt_runs, by = c("POS", "REF", "Nt"))
  colnames(merged_gt) = c("POS", "REF", "ALT", "DP", "Count", "Freq",
                          "Run", "DP Indiv", "Count Indiv", "Freq Indiv")
  merged_gt = merged_gt[order(POS)]
  
  merged_gt$Freq = merged_gt$Freq / 100
  
  # merged_gt = merged_gt[, by = .(POS, REF, DP), .(
  #   Nt = paste(Nt, collapse = ","),
  #   Count = paste(Count, collapse = ","),
  #   Freq = paste(round(Freq, digits = 3), collapse = ",")
  # )]
  
  
  return(merged_gt)
  
}

read_vcf_snvs_TP <- function(path, caller, gt, merged_file) {
    
    if(caller == "Freebayes") {
        
        vcf_df <- read_vcf_freebayes(path, gt, merged_file)
        
    } else if (caller == "Mutect2") {
        
        vcf_df <- read_vcf_mutect2(path, gt, merged_file)
        
    } else if (caller == "LoFreq") {
        
        vcf_df <- read_vcf_LoFreq(path, gt, merged_file)
        
    } else if (caller == "VarDict") {
        
        vcf_df <- read_vcf_VarDict(path, gt, merged_file)
        
    } else if (caller == "VarScan") {
        
        vcf_df <- read_vcf_VarScan(path, gt, merged_file)
        
    }
    
    return(vcf_df)
}


#FP & FN ----------------------------------------------------------------------

load_gt_report <- function(path, merged_file) {
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
    
    # select SNVs
    a_snvs = a[which(ALT %in% c("A", "C", "G", "T")), ]
    #filter DEPTH>2
    a_snvs = a_snvs[which(a_snvs$Count >2), ]
    
    
    gt = list(
        all = a,
        snvs = a_snvs
        
    )
    return(gt)
}

load_gt_vcf <- function(path, merged_file, gt_snvs){
    #function to load Ground Truth vcf
    ground_truth_vcf <- read.vcfR( paste0(path, "/",merged_file, 
                                          "_ground_truth_norm.vcf"),
                                   verbose = FALSE )
    
    ground_truth_vcf  = ground_truth_vcf |> vcfR::getFIX() |> as.data.frame() |> setDT()
    
    pick_gt = gt_snvs[which(gt_snvs$POS %in% ground_truth_vcf$POS)]
    pick_gt$mut = paste(pick_gt$POS, 
                        pick_gt$REF, 
                        pick_gt$ALT, sep = ":")
    return(pick_gt)
}

select_snvs <- function(df){
    # select SNVs from caller based on length of REF and ALT
    snvs = df[nchar(df$REF) == nchar(df$ALT)]
    snvs = snvs[which(nchar(snvs$REF) <2), ]
    snvs = snvs[which(nchar(snvs$ALT) <2), ]
    snvs$mut = paste(snvs$POS, snvs$REF, snvs$ALT, sep = ":")
    
    return(snvs)
}

read_vcf_snvs_FP <- function(path, caller, merged_file, pick_gt, gt_all) {
    
    if(caller == "Freebayes") {
        
        fp_var <- final_fp_snvs_Freebayes(path, merged_file, pick_gt, gt_all)
        
    } else if (caller == "Mutect2") { 
        
        fp_var <- final_fp_snvs_gatk(path, merged_file, pick_gt, gt_all)
        
    } else if (caller == "LoFreq") {
        
        fp_var <- final_fp_snvs_LoFreq(path, merged_file, pick_gt, gt_all)
        
    } else if (caller == "VarDict") {
        
        fp_var <- final_fp_snvs_VarDict(path, merged_file, pick_gt, gt_all)
        
    } else if (caller == "VarScan") {
        
        fp_var <- final_fp_snvs_VarScan(path, merged_file, pick_gt, gt_all)
        
    }
    
    return(fp_var)
}

read_vcf_snvs_FN <- function(path, caller, merged_file, pick_gt) {
    
    if(caller == "Freebayes") {
        
        fn_var <- final_fn_snvs_Freebayes(path, merged_file, pick_gt)
        
    } else if (caller == "Mutect2") {
        
        fn_var <- final_fn_snvs_gatk(path, merged_file, pick_gt)
        
    } else if (caller == "LoFreq") {
        
        fn_var <- final_fn_snvs_LoFreq(path, merged_file, pick_gt)
        
    } else if (caller == "VarDict") {
        
        fn_var <- final_fn_snvs_VarDict(path, merged_file, pick_gt)
        
    } else if (caller == "VarScan") {
        
        fn_var <- final_fn_snvs_VarScan(path, merged_file, pick_gt)
        
    }
    
    return(fn_var)
}
