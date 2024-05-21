



read_vcf <- function(path, caller, gt) {
  
  if(caller == "freebayes") {
    
    vcf_df <- read_vcf_freebays(path, gt)
    
  } else if (caller == "mutect2") {
    
    vcf_df <- read_vcf_mutect2(path, gt)
    
  }
  
  
  
}


read_vcf_freebays <- function(path, gt) {
  
  vcf <- read.vcfR( path, verbose = FALSE )
  
  vcf_df <- gt |>
    merge_freebayes(gt) |>
    clean_freebayes()
  
  return(vcf_df)
    
}

read_vcf_mutect2 <- function(path, gt) {
  
  vcf <- read.vcfR( path, verbose = FALSE )
  
  vcf_df = vcf |>
    merge_gatk(gt) |>
    clean_gatk()
  
  
  return(vcf_df)
  
}




