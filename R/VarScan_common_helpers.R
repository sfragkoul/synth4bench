
#TP SNVS-----------------------------------------------------------------------
read_vcf_VarScan <- function(path, gt, merged_file) {
  #takes two files and produce a caller vcf file in a certain format 
  vcf <- read.vcfR(paste0(path, "/", merged_file, "_VarScan_norm.vcf"), verbose = FALSE )
  
  df = (vcf |>
    merge_VarScan(gt) |>
    clean_VarScan())
  
  return(df)
  
}

merge_VarScan <- function(VarScan_somatic_vcf, merged_gt) {
    #return cleaned vcf
    VarScan_s0  = VarScan_somatic_vcf |> vcfR::getFIX() |> as.data.frame() |> setDT()
    #VarScan_s1  = VarScan_somatic_vcf |> extract_gt_tidy() |> setDT()
    VarScan_s2 = VarScan_somatic_vcf |> extract_info_tidy() |> setDT()
    VarScan_s2 = VarScan_s2[,c( "DP", "Pvalue", "AF" )]
    
    VarScan_s0 = VarScan_s0[which(VarScan_s2$AF>0.0)]
    VarScan_s2 = VarScan_s2[which(VarScan_s2$AF>0.0)]
    
    
    VarScan_somatic = cbind(VarScan_s0, VarScan_s2)
    
    
    #Merge everything into a common file
    merged_gt$POS = as.character(merged_gt$POS)
    
    merged_bnch = merge(merged_gt, VarScan_somatic,  by = "POS", all.x = TRUE)
    
    merged_bnch$POS = as.numeric(merged_bnch$POS)
    
    merged_bnch = merged_bnch[order(POS)]
    
    colnames(merged_bnch) = c(
        "POS",	"Ground Truth REF",	"Ground Truth DP",
        "Ground Truth ALT", "Ground Truth AD", 
        "Ground Truth AF", "CHROM", "ID", "VarScan REF",	
        "VarScan ALT", "VarScan QUAL",	"VarScan FILTER", "VarScan DP", "Pvalue","VarScan AF"
    )
    
    return(merged_bnch)
    
}


clean_VarScan <- function(df) {
    #function to produce the caller's reported variants in the desired format 
    df2 = df[, c(
        "POS",
        
        "Ground Truth REF",
        "Ground Truth ALT",
        "Ground Truth DP",
        "Ground Truth AF",
        
        "VarScan REF", 
        "VarScan ALT", 
        "VarScan DP",
        "VarScan AF"
    ), with = FALSE]
    
    
    
    df2 = df2[, by = c(
        "POS",
        "Ground Truth REF",
        "Ground Truth DP",
        "VarScan REF", 
        "VarScan ALT", 
        "VarScan DP",
        "VarScan AF"
        
    ), .(
        "Ground Truth ALT" = `Ground Truth ALT` |> tstrsplit(",") |> unlist(),
        "Ground Truth AF"  = `Ground Truth AF` |> tstrsplit(",") |> unlist()
    )]



    VarScan_alt = str_split(df2$`VarScan ALT`, ",")
    VarScan_af = str_split(df2$`VarScan AF`, ",")


    cln = mapply(
        function(x, y, z) {

            index = which(y == x)

            return(
                c(y[index], z[index])
            )

        },

        df2$`Ground Truth ALT`, VarScan_alt, VarScan_af
    )


    df2$`VarScan ALT` = cln |> lapply(function(x) { return(x [1]) }) |> unlist()
    df2$`VarScan AF`  = cln |> lapply(function(x) { return(x [2]) }) |> unlist()
    
    df2[which(is.na(`VarScan AF`))]$`VarScan DP` = NA
    df2[which(is.na(`VarScan AF`))]$`VarScan REF` = NA
    
    df2 = df2[, c(
        "POS", 
        "Ground Truth REF",
        "Ground Truth ALT",
        "Ground Truth DP",
        "Ground Truth AF",
        "VarScan REF", 
        "VarScan ALT", 
        "VarScan DP",
        "VarScan AF"
    ), with = FALSE]

    return(df2)
    
}

load_VarScan_vcf <- function(path, merged_file){
    #function to load caller vcf
    VarScan_somatic_vcf <- read.vcfR( paste0(path, "/", merged_file, 
                                             "_VarScan_norm.vcf"), verbose = FALSE )
    VarScan_s0  = VarScan_somatic_vcf |> vcfR::getFIX() |> as.data.frame() |> setDT()
    #VarScan_s1  = VarScan_somatic_vcf |> extract_gt_tidy() |> setDT()
    VarScan_s2 = VarScan_somatic_vcf |> extract_info_tidy() |> setDT()
    VarScan_s2 = VarScan_s2[,c( "DP", "AF" )]
    VarScan_somatic = cbind(VarScan_s0, VarScan_s2)
    return(VarScan_somatic)
}
