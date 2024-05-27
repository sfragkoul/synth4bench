

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

read_vcf <- function(path, caller, gt) {
    
    if(caller == "freebayes") {
        
        vcf_df <- read_vcf_freebays(path, gt)
        
    } else if (caller == "mutect2") {
        
        vcf_df <- read_vcf_mutect2(path, gt)
        
    } else if (caller == "LoFreq") {
        
        vcf_df <- read_vcf_LoFreq(path, gt)
        
    } else if (caller == "VarDict") {
        
        vcf_df <- read_vcf_VarDict(path, gt)
        
    } else if (caller == "VarScan") {
        
        vcf_df <- read_vcf_VarScan(path, gt)
        
    }
    
    return(vcf_df)
}

plot_synth4bench <- function(gt_comparison, vcf_path, gt_path, caller) {
    
    
    df = fread( gt_comparison )
    
    vcf_GT <- read.vcfR(gt_path, verbose = FALSE )
    
    vcf_caller <- read.vcfR(vcf_path, verbose = FALSE )
    
    if(caller == "freebayes") {
        
        plots <- plot_synth4bench_freebayes(df, vcf_GT, vcf_caller)
        
    } else if (caller == "mutect2") {
        
        plots <- plot_synth4bench_gatk(df, vcf_GT, vcf_caller)
        
    } else if (caller == "LoFreq") {
        
        plots <- plot_synth4bench_LoFreq(df, vcf_GT, vcf_caller)
        
    } else if (caller == "VarDict") {
        
        plots <- plot_synth4bench_VarDict(df, vcf_GT, vcf_caller)
        
    } else if (caller == "VarScan") {
        
        plots <- plot_synth4bench_VarScan(df, vcf_GT, vcf_caller)
        
    }
    
    
    return(plots)
    
}


