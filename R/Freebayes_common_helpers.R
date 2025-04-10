
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
    ## Extract relevant columns
    df2 <- df[, c(
        "POS",
        "Ground Truth REF", "Ground Truth ALT", "Ground Truth DP", "Ground Truth AF",
        "Freebayes REF", "Freebayes ALT", "Freebayes DP", "Freebayes AO"
    ), with = FALSE]
    
    ## Expand multiallelic ground-truth columns into separate rows
    df2 <- df2[, by = .(POS, `Ground Truth REF`, `Ground Truth DP`), .(
        "Ground Truth ALT" = tstrsplit(`Ground Truth ALT`, ",") |> unlist(),
        "Ground Truth AF"  = tstrsplit(`Ground Truth AF`, ",")  |> unlist(),
        "Freebayes REF" = `Freebayes REF`[1],
        "Freebayes ALT" = `Freebayes ALT`[1],
        "Freebayes DP"  = `Freebayes DP`[1],
        "Freebayes AO"  = `Freebayes AO`[1]
    )]
    
    ## Match ALT alleles between Ground Truth and Freebayes caller
    df2[, `:=` (
        # Check if the Ground Truth ALT allele is present in the comma‐separated Freebayes ALT field
        `ALT Match` = mapply(function(gt_alt, fb_alt) {
            if (is.na(fb_alt)) return(FALSE)
            return(gt_alt %in% unlist(str_split(fb_alt, ",")))
        }, `Ground Truth ALT`, `Freebayes ALT`),
        
        # For matching rows, get the AO value from the allele that matches the ground truth
        `AO Match` = mapply(function(gt_alt, fb_alt, fb_ao) {
            if (is.na(fb_alt) || is.na(fb_ao)) return(NA)
            alt_list <- str_split(fb_alt, ",")[[1]]
            ao_list  <- str_split(fb_ao, ",")[[1]]
            idx <- which(alt_list == gt_alt)
            if (length(idx) == 0) return(NA)
            return(ao_list[[idx]])
        }, `Ground Truth ALT`, `Freebayes ALT`, `Freebayes AO`)
    )]
    
    ## Retain only matching alleles and adjust caller fields
    df2[, `Freebayes ALT` := ifelse(`ALT Match`, `Ground Truth ALT`, NA)]
    df2[, `Freebayes AO`  := `AO Match`]
    df2[, `Freebayes REF` := ifelse(`ALT Match`, `Freebayes REF`, NA)]
    df2[, `Freebayes DP`  := ifelse(`ALT Match`, `Freebayes DP`, NA)]
    
    ## Compute Freebayes AF as AO/DP
    # Note: this will be NA if either AO or DP is NA.
    df2[, `Freebayes AF` := as.numeric(`Freebayes AO`) / as.numeric(`Freebayes DP`)]
    
    ## Classify as TP or FN based on matching ALT allele
    df2[, type := ifelse(is.na(`Freebayes ALT`), "FN", "TP")]
    
    ## Compute the allele frequency deviation for true positives.
    df2[, `AF Deviation` := NA_real_]
    df2[type == "TP", `AF Deviation` := as.numeric(`Freebayes AF`) - as.numeric(`Ground Truth AF`)]
    
    ## Final output table (you can rearrange columns as needed)
    final_table <- df2[, .(
        POS,
        `Ground Truth REF`, `Ground Truth ALT`, `Ground Truth DP`, `Ground Truth AF`,
        `Freebayes REF`, `Freebayes ALT`, `Freebayes DP`, `Freebayes AO`, `Freebayes AF`,
        type, `AF Deviation`
    )]
    
    ## Calculate recall (e.g. proportion of calls with a non-NA caller REF)
    recall <- sum(!is.na(final_table$`Freebayes REF`)) / nrow(final_table)
    
    return(list(
        "vcf_snvs_cleaned" = final_table,
        "recall" = recall
    ))
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
