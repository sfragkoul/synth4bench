
#True Variants SNVS------------------------------------------------------------
read_vcf_VarDict <- function(path, gt, merged_file) {
  #takes two files and produce a caller vcf file in a certain format 
  vcf <- read.vcfR(paste0(path, "/", merged_file, "_VarDict_norm.vcf"), verbose = FALSE )
  
  vcf_df = (vcf |>
    merge_VarDict(gt) |>
    clean_VarDict())
  
  return(vcf_df)
  
}

merge_VarDict <- function(VarDict_somatic_vcf, merged_gt) {
    #return cleaned vcf
    VarDict_s0  = VarDict_somatic_vcf |> vcfR::getFIX() |> as.data.frame() |> setDT()
    VarDict_s1  = VarDict_somatic_vcf |> extract_gt_tidy() |> setDT()
    VarDictgatk_s21 = VarDict_somatic_vcf |> extract_info_tidy() |> setDT()
    
    VarDict_somatic = cbind(VarDict_s0[VarDict_s1$Key, ], VarDict_s1)
    
    
    #Merge everything into a common file
    merged_gt$POS = as.character(merged_gt$POS)
    
    merged_bnch = merge(merged_gt, VarDict_somatic,  by = "POS", all.x = TRUE)
    
    merged_bnch$POS = as.numeric(merged_bnch$POS)
    
    merged_bnch = merged_bnch[order(POS)]
    
    colnames(merged_bnch) = c(
        "POS",	"Ground Truth REF",	"Ground Truth DP",
        "Ground Truth ALT", "Ground Truth AD", 
        "Ground Truth AF", "CHROM", "ID", "VarDict REF",	
        "VarDict ALT", "VarDict QUAL",	"VarDict FILTER",
        "key", "Indiv", "gt_GT", "VarDict DP", "gt_VD", "VarDict AD", 
        "VarDict AF", "gt_RD", "gt_ALD", "gt_GT_alleles"
    )
    
    return(merged_bnch)
    
}

clean_VarDict <- function(df) {
    #function to produce the caller's reported variants in the desired format 
    df2 = df[, c(
        "POS",
        
        "Ground Truth REF",
        "Ground Truth ALT",
        "Ground Truth DP",
        "Ground Truth AF",
        
        "VarDict REF", 
        "VarDict ALT", 
        "VarDict DP",
        "VarDict AF"
    ), with = FALSE]
    
    
    
    df2 = df2[, by = c(
        "POS",
        "Ground Truth REF",
        "Ground Truth DP",
        "VarDict REF", 
        "VarDict ALT", 
        "VarDict DP",
        "VarDict AF"
        
    ), .(
        "Ground Truth ALT" = `Ground Truth ALT` |> tstrsplit(",") |> unlist(),
        "Ground Truth AF"  = `Ground Truth AF` |> tstrsplit(",") |> unlist()
    )]



    VarDict_alt = str_split(df2$`VarDict ALT`, ",")
    VarDict_af = str_split(df2$`VarDict AF`, ",")


    cln = mapply(
        function(x, y, z) {

            index = which(y == x)

            return(
                c(y[index], z[index])
            )

        },

        df2$`Ground Truth ALT`, VarDict_alt, VarDict_af
    )


    df2$`VarDict ALT` = cln |> lapply(function(x) { return(x [1]) }) |> unlist()
    df2$`VarDict AF`  = cln |> lapply(function(x) { return(x [2]) }) |> unlist()
    
    df2[which(is.na(`VarDict AF`))]$`VarDict DP` = NA
    df2[which(is.na(`VarDict AF`))]$`VarDict REF` = NA
    
    df2 = df2[, c(
        "POS", 
        "Ground Truth REF",
        "Ground Truth ALT",
        "Ground Truth DP",
        "Ground Truth AF",
        "VarDict REF", 
        "VarDict ALT", 
        "VarDict DP",
        "VarDict AF"
    ), with = FALSE]

    return(df2)
    
}

load_VarDict_vcf <- function(path, merged_file){
    #function to load caller vcf
    VarDict_somatic_vcf <- read.vcfR( paste0(path, "/",merged_file, 
                                             "_VarDict_norm.vcf"), verbose = FALSE )
    VarDict_s0  = VarDict_somatic_vcf |> vcfR::getFIX() |> as.data.frame() |> setDT()
    #VarDict_s1  = VarDict_somatic_vcf |> extract_gt_tidy() |> setDT()
    VarDict_s2 = VarDict_somatic_vcf |> extract_info_tidy() |> setDT()
    VarDict_s2 = VarDict_s2[,c( "DP", "AF" )]
    VarDict_somatic = cbind(VarDict_s0, VarDict_s2)
    return(VarDict_somatic)
}
