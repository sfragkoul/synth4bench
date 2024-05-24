

source("R/libraries.R")

source("R/helpers_freebayes.R")
source("R/helpers_gatk.R")
source("R/helpers_LoFreq.R")
source("R/helpers_VarDict.R")
source("R/helpers_VarScan.R")

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












