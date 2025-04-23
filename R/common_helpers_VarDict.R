
#True Variants SNVS------------------------------------------------------------
read_vcf_VarDict <- function(path, gt, merged_file) {
  #takes two files and produce a caller vcf file in a certain format 
  vcf <- read.vcfR(paste0(path, "/", merged_file, "_VarDict_norm.vcf"), verbose = FALSE )
  
  df = (vcf |>
    merge_VarDict(gt) |>
    clean_VarDict())
  
  return(df)
  
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
        "POS",	"Ground Truth REF",	"Ground Truth ALT",
        "Ground Truth DP", "Ground Truth AD", "Ground Truth AF", 
        
        "Run", "DP Indiv", "Count Indiv", "Freq Indiv", "mut",
        
        "CHROM", "ID", "VarDict REF",	
        "VarDict ALT", "VarDict QUAL",	"VarDict FILTER",
        "key", "Indiv", "gt_GT", "VarDict DP", "gt_VD", "VarDict AD", 
        "VarDict AF", "gt_RD", "gt_ALD", "gt_GT_alleles"
    )
    
    #after unlisting multiple variants in the same position, we must
    # keep only unique FN POS
    merged_bnch <- merged_bnch[, .SD[1], by = POS]
    
    return(
        list(
            "merged_bnch" = merged_bnch,
            "VarDict_somatic" = VarDict_somatic)
    )
    
}

clean_VarDict <- function(df) {
    # Extract relevant columns
    df2 <- df$merged_bnch[, c(
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
    
    # Expand multiallelic GT sites into separate rows
    df2 <- df2[, by = .(POS, `Ground Truth REF`, `Ground Truth DP`), .(
        "Ground Truth ALT" = tstrsplit(`Ground Truth ALT`, ",") |> unlist(),
        "Ground Truth AF"  = tstrsplit(`Ground Truth AF`, ",") |> unlist(),
        "VarDict REF" = `VarDict REF`[1],
        "VarDict ALT" = `VarDict ALT`[1],
        "VarDict DP"  = `VarDict DP`[1],
        "VarDict AF"  = `VarDict AF`[1]
    )]
    
    # Match ALT alleles between GT and VarDict
    df2[, `:=` (
        `ALT Match` = mapply(function(gt_alt, gatk_alt) {
            if (is.na(gatk_alt)) return(FALSE)
            return(gt_alt %in% unlist(str_split(gatk_alt, ",")))
        }, `Ground Truth ALT`, `VarDict ALT`),
        
        `AF Match` = mapply(function(gt_alt, gatk_alt, gatk_af) {
            if (is.na(gatk_alt) | is.na(gatk_af)) return(NA)
            alt_list <- str_split(gatk_alt, ",")[[1]]
            af_list  <- str_split(gatk_af, ",")[[1]]
            idx <- which(alt_list == gt_alt)
            if (length(idx) == 0) return(NA)
            return(af_list[[idx]])
        }, `Ground Truth ALT`, `VarDict ALT`, `VarDict AF`)
    )]
    
    # Keep only matched alleles
    df2[, `VarDict ALT` := ifelse(`ALT Match`, `Ground Truth ALT`, NA)]
    df2[, `VarDict AF`  := `AF Match`]
    df2[, `VarDict REF` := ifelse(`ALT Match`, `VarDict REF`, NA)]
    df2[, `VarDict DP`  := ifelse(`ALT Match`, `VarDict DP`, NA)]
    
    # Classify as TP or FN
    df2[, type := ifelse(is.na(`VarDict ALT`), "FN", "TP")]
    
    # ΔAF: Caller AF - Ground Truth AF (numeric)
    df2[, `AF Deviation ` := NA_real_]
    df2[type == "TP", `AF Deviation` := as.numeric(`VarDict AF`) - as.numeric(`Ground Truth AF`)]
    
    # Final output
    df2 <- df2[, .(
        POS,
        `Ground Truth REF`, `Ground Truth ALT`, `Ground Truth DP`, `Ground Truth AF`,
        `VarDict REF`, `VarDict ALT`, `VarDict DP`, `VarDict AF`,
        type, `AF Deviation`
    )]
    
    
    recall = sum(!is.na(df2$`VarDict REF`)) / nrow(df2)
    
    return(list(
        "vcf_snvs_cleaned" = df2,
        "recall" = recall))
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
