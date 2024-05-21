source("libraries.R")

read_vcf <- function(path, caller, gt) {
  
  if(caller == "freebayes") {
    
    vcf_df <- read_vcf_freebays(path, gt)
    
  } else if (caller == "mutect2") {
    
    vcf_df <- read_vcf_mutect2(path, gt)
    
  }else if (caller == "LoFreq") {
      
      vcf_df <- read_vcf_LoFreq(path, gt)
      
  }else if (caller == "VarDict") {
      
      vcf_df <- read_vcf_VarDict(path, gt)
      
  }else if (caller == "VarScan") {
      
      vcf_df <- read_vcf_VarScan(path, gt)
      
  }
  
  
}

source("helpers_freebayes.R")
read_vcf_freebays <- function(path, gt) {
  
  vcf <- read.vcfR( path, verbose = FALSE )
  
  vcf_df <- gt |>
    merge_freebayes(gt) |>
    clean_freebayes()
  
  return(vcf_df)
    
}

source("helpers_gatk.R")
read_vcf_mutect2 <- function(path, gt) {
  
  vcf <- read.vcfR( path, verbose = FALSE )
  
  vcf_df = vcf |>
    merge_gatk(gt) |>
    clean_gatk()
  
  return(vcf_df)
  
}

source("helpers_LoFreq.R")
read_vcf_LoFreq <- function(path, gt) {
    
    vcf <- read.vcfR( path, verbose = FALSE )
    
    vcf_df = vcf |>
        merge_LoFreq(gt) |>
        clean_LoFreq()
    
    return(vcf_df)
    
}

source("helpers_VarDict.R")
read_vcf_VarDict <- function(path, gt) {
    
    vcf <- read.vcfR( path, verbose = FALSE )
    
    vcf_df = vcf |>
        merge_VarDict(gt) |>
        clean_VarDict()
    
    return(vcf_df)
    
}

source("helpers_VarScan.R")
read_vcf_VarScan <- function(path, gt) {
    
    vcf <- read.vcfR( path, verbose = FALSE )
    
    vcf_df = vcf |>
        merge_VarScan(gt) |>
        clean_VarScan()
    
    return(vcf_df)
    
}
