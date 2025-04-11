#True Variants SNVS------------------------------------------------------------
read_vcf_LoFreq <- function(path, gt, merged_file) {
  #takes two files and produce a caller vcf file in a certain format 
  vcf <- read.vcfR(paste0(path, "/", merged_file, "_LoFreq_norm.vcf"), verbose = FALSE )
  
  vcf_df = vcf |>
    merge_LoFreq(gt) |>
    clean_LoFreq()
  
  return(vcf_df)
  
}

merge_LoFreq <- function(LoFreq_somatic_vcf, merged_gt) {
    #return cleaned vcf
    LoFreq_s0  = LoFreq_somatic_vcf |> vcfR::getFIX() |> as.data.frame() |> setDT()
    #LoFreq_s1  = LoFreq_somatic_vcf |> extract_gt_tidy() |> setDT()
    LoFreq_s2 = LoFreq_somatic_vcf |> extract_info_tidy() |> setDT()
    LoFreq_s2 = LoFreq_s2[,c( "DP", "AF" )]
    
    LoFreq_somatic = cbind(LoFreq_s0, LoFreq_s2)
    
    #Merge everything into a common file
    merged_gt$POS = as.character(merged_gt$POS)
    
    merged_bnch = merge(merged_gt, LoFreq_somatic,  by = "POS", all.x = TRUE)
    
    merged_bnch$POS = as.numeric(merged_bnch$POS)
    
    merged_bnch = merged_bnch[order(POS)]
    
    colnames(merged_bnch) = c(
        "POS",	"Ground Truth REF",	"Ground Truth DP",
        "Ground Truth ALT", "Ground Truth AD", 
        "Ground Truth AF", "CHROM", "ID",	"LoFreq REF",	
        "LoFreq ALT", "LoFreq QUAL", "LoFreq FILTER", 
        "LoFreq DP", "LoFreq AF"
    )
    
    return(merged_bnch)
    
}

clean_LoFreq <- function(df) {
    #function to produce the caller's reported variants in the desired format 
    df2 = df[, c(
        "POS",
        
        "Ground Truth REF",
        "Ground Truth ALT",
        "Ground Truth DP",
        "Ground Truth AF",
        
        "LoFreq REF",
        "LoFreq ALT",
        "LoFreq DP",
        "LoFreq AF"
    ), with = FALSE]
    
    
    
    df2 = df2[, by = c(
        "POS",
        "Ground Truth REF",
        "Ground Truth DP",
        "LoFreq REF",
        "LoFreq ALT",
        "LoFreq DP",
        "LoFreq AF"

    ), .(
        "Ground Truth ALT" = `Ground Truth ALT` |> tstrsplit(",") |> unlist(),
        "Ground Truth AF"  = `Ground Truth AF` |> tstrsplit(",") |> unlist()
        # "LoFreq REF" = `LoFreq REF` |> tstrsplit(",") |> unlist(),
        # "LoFreq ALT" = `LoFreq ALT` |> tstrsplit(",") |> unlist(),
        # "LoFreq DP"  = `LoFreq DP` |> tstrsplit(",") |> unlist() |> as.integer(),
        # "LoFreq AF"  = `LoFreq AF` |> tstrsplit(",") |> unlist() |> as.numeric()
    )]



    LoFreq_alt = str_split(df2$`LoFreq ALT`, ",")
    LoFreq_af = str_split(df2$`LoFreq AF`, ",")


    cln = mapply(
        function(x, y, z) {

            index = which(y == x)

            return(
                c(y[index], z[index])
            )

        },

        df2$`Ground Truth ALT`, LoFreq_alt, LoFreq_af
    )


    df2$`LoFreq ALT` = cln |> lapply(function(x) { return(x [1]) }) |> unlist()
    df2$`LoFreq AF` = cln |> lapply(function(x) { return(x [2]) }) |> unlist()

    df2[which(is.na(`LoFreq AF`))]$`LoFreq DP` = NA
    df2[which(is.na(`LoFreq AF`))]$`LoFreq REF` = NA

    df2 = df2[, c(
        "POS",
        "Ground Truth REF",
        "Ground Truth ALT",
        "Ground Truth DP",
        "Ground Truth AF",
        "LoFreq REF",
        "LoFreq ALT",
        "LoFreq DP",
        "LoFreq AF"
    ), with = FALSE]
    
    return(df2)
    
}

load_LoFreq_vcf <- function(path, merged_file){
    #function to load caller vcf
    LoFreq_somatic_vcf <- read.vcfR( paste0(path, "/", merged_file, 
                                            "_LoFreq_norm.vcf"), verbose = FALSE )
    LoFreq_s0  = LoFreq_somatic_vcf |> vcfR::getFIX() |> as.data.frame() |> setDT()
    #LoFreq_s1  = LoFreq_somatic_vcf |> extract_gt_tidy() |> setDT()
    LoFreq_s2 = LoFreq_somatic_vcf |> extract_info_tidy() |> setDT()
    LoFreq_s2 = LoFreq_s2[,c( "DP", "AF" )]
    LoFreq_somatic = cbind(LoFreq_s0, LoFreq_s2)
    return(LoFreq_somatic)
}
