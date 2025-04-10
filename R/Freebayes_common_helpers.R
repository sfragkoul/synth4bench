
#True Variants SNVS------------------------------------------------------------
read_vcf_freebayes <- function(path, gt, merged_file) {
  #takes two files and produce a caller vcf file in a certain format   
  vcf <- read.vcfR( paste0(path, "/", merged_file, "_freebayes_norm.vcf"), verbose = FALSE )
  
  vcf_df <- vcf |>
    merge_freebayes(gt) |>
    clean_freebayes()
  
  return(vcf_df)
  
}

merge_freebayes <- function(freebayes_somatic_vcf, merged_gt) {
    #return cleaned vcf
    freebayes_s0  = freebayes_somatic_vcf |> vcfR::getFIX() |> as.data.frame() |> setDT()
    freebayes_s1  = freebayes_somatic_vcf |> extract_gt_tidy() |> setDT()
    freebayesgatk_s21 = freebayes_somatic_vcf |> extract_info_tidy() |> setDT()
    
    freebayes_somatic = cbind(freebayes_s0[freebayes_s1$Key, ], freebayes_s1)
    
    
    #Merge everything into a common file
    merged_gt$POS = as.character(merged_gt$POS)
    merged_bnch = merge(merged_gt, freebayes_somatic,  by = "POS", all.x = TRUE)
    merged_bnch$POS = as.numeric(merged_bnch$POS)
    merged_bnch = merged_bnch[order(POS)]
    
    colnames(merged_bnch) = c(
        "POS",	"Ground Truth REF",	"Ground Truth DP",
        "Ground Truth ALT", "Ground Truth AD", "Ground Truth AF", 
        
        "Run", "DP Indiv", "Count Indiv", "Freq Indiv", 
        
        "Freebayes CHROM", "Freebayes ID", "Freebayes REF", "Freebayes ALT", 
        "Freebayes QUAL", "Freebayes FILTER", "Freebayes key", 
        "Freebayes Indiv", "Freebayes GT", "Freebayes GQ", "Freebayes GL", 
        "Freebayes DP", "Freebayes RO", "Freebayes QR", "Freebayes AO", 
        "Freebayes QA", "Freebayes alleles"
    )
    
    #after unlisting multiple variants in the same position, we must
    # keep only unique FN POS
    merged_bnch <- merged_bnch[, .SD[1], by = POS]
    
    return(
        list(
            "merged_bnch" = merged_bnch,
            "freebayes_somatic" = freebayes_somatic)
    )
    
}


clean_freebayes <- function(df) {
    #function to produce the caller's reported variants in the desired format 
    df2 = df[, c(
        "POS",
        
        "Ground Truth REF",
        "Ground Truth ALT",
        "Ground Truth DP",
        "Ground Truth AF",
        
        "Freebayes REF", 
        "Freebayes ALT", 
        "Freebayes DP",
        "Freebayes AO"
    ), with = FALSE]
    
    
    
    df2 = df2[, by = c(
        "POS",
        "Ground Truth REF",
        "Ground Truth DP",
        "Freebayes REF", 
        "Freebayes ALT", 
        "Freebayes DP",
        "Freebayes AO"
        
    ), .(
        "Ground Truth ALT" = `Ground Truth ALT` |> tstrsplit(",") |> unlist(),
        "Ground Truth AF"  = `Ground Truth AF` |> tstrsplit(",") |> unlist()
    )]



    freebayes_alt = str_split(df2$`Freebayes ALT`, ",")
    freebayes_ao = str_split(df2$`Freebayes AO`, ",")


    cln = mapply(
        function(x, y, z) {

            index = which(y == x)

            return(
                c(y[index], z[index])
            )

        },

        df2$`Ground Truth ALT`, freebayes_alt, freebayes_ao
    )


    df2$`Freebayes ALT` = cln |> lapply(function(x) { return(x [1]) }) |> unlist()
    df2$`Freebayes AO`  = cln |> lapply(function(x) { return(x [2]) }) |> unlist()
    
    df2[which(is.na(`Freebayes AO`))]$`Freebayes DP` = NA
    df2[which(is.na(`Freebayes AO`))]$`Freebayes REF` = NA
    
    df2 = df2[, c(
        "POS", 
        "Ground Truth REF",
        "Ground Truth ALT",
        "Ground Truth DP",
        "Ground Truth AF",
        "Freebayes REF", 
        "Freebayes ALT", 
        "Freebayes DP",
        "Freebayes AO"
    ), with = FALSE]
    
    #AF = AO/DP
    df2$"Freebayes AF" =  as.numeric(format(round(as.numeric(df2$"Freebayes AO")/ as.numeric(df2$"Freebayes DP"), 3)))
    
    return(df2)
    
}

load_Freebayes_vcf <- function(path, merged_file){
    #function to load caller vcf
    Freebayes_somatic_vcf <- read.vcfR( paste0(path, "/", merged_file, 
                                               "_freebayes_norm.vcf"), verbose = FALSE )
    Freebayes_s0  = Freebayes_somatic_vcf |> vcfR::getFIX() |> as.data.frame() |> setDT()
    #Freebayes_s1  = Freebayes_somatic_vcf |> extract_gt_tidy() |> setDT()
    Freebayes_s2 = Freebayes_somatic_vcf |> extract_info_tidy() |> setDT()
    Freebayes_s2 = Freebayes_s2[,c( "DP", "AF" )]
    Freebayes_somatic = cbind(Freebayes_s0, Freebayes_s2)
    return(Freebayes_somatic)
}
