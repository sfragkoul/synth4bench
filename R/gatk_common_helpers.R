
#TP SNVS-----------------------------------------------------------------------
read_vcf_mutect2 <- function(path, gt, merged_file) {
  #takes two files and produce a caller vcf file in a certain format 
  vcf <- read.vcfR(paste0(path, "/", merged_file, "_Mutect2_norm.vcf"), verbose = FALSE )
  
   df = (vcf |>
    merge_gatk(gt) |>
    clean_gatk())
  
  return(df)
  
}

merge_gatk <- function(gatk_somatic_vcf, merged_gt) {
    #return cleaned vcf
    gatk_s0  = gatk_somatic_vcf |> vcfR::getFIX() |> as.data.frame() |> setDT()
    gatk_s1  = gatk_somatic_vcf |> extract_gt_tidy() |> setDT()
    gatk_s21 = gatk_somatic_vcf |> extract_info_tidy() |> setDT()
    gatk_somatic = cbind(gatk_s0[gatk_s1$Key, ], gatk_s1)
    
    #Merge everything into a common file
    merged_gt$POS = as.character(merged_gt$POS)
    merged_bnch = merge(merged_gt, gatk_somatic,  by = "POS", all.x = TRUE)
    merged_bnch$POS = as.numeric(merged_bnch$POS)
    merged_bnch = merged_bnch[order(POS)]

    colnames(merged_bnch) = c(
        "POS",	"Ground Truth REF",	"Ground Truth ALT",
        "Ground Truth DP", "Ground Truth AD", 
        "Ground Truth AF", "Run", "DP Indiv", "Count Indiv", 
        "Freq Indiv", "CHROM", "ID",	"Mutect2 REF",	
        "Mutect2 ALT", "Mutect2 QUAL",	"Mutect2 FILTER",
        "key", "Indiv", "Mutect2 AD", "Mutect2 AF",
        "Mutect2 DP", "gt_F1R2", "gt_F2R1", "gt_FAD",	
        "gt_GQ", "gt_GT",	"gt_PGT",	"gt_PID",	"gt_PL",
        "gt_PS",	"gt_SB",	"gt_GT_alleles"
    )
    
    #after unlisting multiple variants in the same position, we must
    # keep only unique FN POS
    merged_bnch <- merged_bnch[, .SD[1], by = POS]
    
    return(
        list(
            "merged_bnch" = merged_bnch,
            "gatk_somatic" = gatk_somatic)
    )
}


clean_gatk <- function(df) {
    # Extract relevant columns
    df2 <- df$merged_bnch[, c(
        "POS", 
        
        "Ground Truth REF", 
        "Ground Truth ALT", 
        "Ground Truth DP", 
        "Ground Truth AF",
        
        "Mutect2 REF", 
        "Mutect2 ALT", 
        "Mutect2 DP", 
        "Mutect2 AF"
    ), with = FALSE]
    
    # Expand multiallelic GT sites into separate rows
    df2 <- df2[, by = .(POS, `Ground Truth REF`, `Ground Truth DP`), .(
        "Ground Truth ALT" = tstrsplit(`Ground Truth ALT`, ",") |> unlist(),
        "Ground Truth AF"  = tstrsplit(`Ground Truth AF`, ",") |> unlist(),
        "Mutect2 REF" = `Mutect2 REF`[1],
        "Mutect2 ALT" = `Mutect2 ALT`[1],
        "Mutect2 DP"  = `Mutect2 DP`[1],
        "Mutect2 AF"  = `Mutect2 AF`[1]
    )]
    
    # Match ALT alleles between GT and Mutect2
    df2[, `:=` (
        `ALT Match` = mapply(function(gt_alt, gatk_alt) {
            if (is.na(gatk_alt)) return(FALSE)
            return(gt_alt %in% unlist(str_split(gatk_alt, ",")))
        }, `Ground Truth ALT`, `Mutect2 ALT`),
        
        `AF Match` = mapply(function(gt_alt, gatk_alt, gatk_af) {
            if (is.na(gatk_alt) | is.na(gatk_af)) return(NA)
            alt_list <- str_split(gatk_alt, ",")[[1]]
            af_list  <- str_split(gatk_af, ",")[[1]]
            idx <- which(alt_list == gt_alt)
            if (length(idx) == 0) return(NA)
            return(af_list[[idx]])
        }, `Ground Truth ALT`, `Mutect2 ALT`, `Mutect2 AF`)
    )]
    
    # Keep only matched alleles
    df2[, `Mutect2 ALT` := ifelse(`ALT Match`, `Ground Truth ALT`, NA)]
    df2[, `Mutect2 AF`  := `AF Match`]
    df2[, `Mutect2 REF` := ifelse(`ALT Match`, `Mutect2 REF`, NA)]
    df2[, `Mutect2 DP`  := ifelse(`ALT Match`, `Mutect2 DP`, NA)]
    
    # Classify as TP or FN
    df2[, type := ifelse(is.na(`Mutect2 ALT`), "FN", "TP")]
    
    # ΔAF: Caller AF - Ground Truth AF (numeric)
    df2[, `AF Deviation ` := NA_real_]
    df2[type == "TP", `AF Deviation` := as.numeric(`Mutect2 AF`) - as.numeric(`Ground Truth AF`)]
    
    # Final output
    df2 <- df2[, .(
        POS,
        `Ground Truth REF`, `Ground Truth ALT`, `Ground Truth DP`, `Ground Truth AF`,
        `Mutect2 REF`, `Mutect2 ALT`, `Mutect2 DP`, `Mutect2 AF`,
        type, `AF Deviation`
    )]
    
    
    recall = sum(!is.na(df2$`Mutect2 REF`)) / dim(df2)[1]
    
    return(list(
        "df2" = df2,
        "recall" = recall))
}


load_gatk_vcf <- function(path, merged_file){
    #function to load caller vcf
    Mutect2_somatic_vcf <- read.vcfR( paste0(path, "/", merged_file, 
                                             "_Mutect2_norm.vcf"), verbose = FALSE )
    
    Mutect2_s0  = Mutect2_somatic_vcf |> vcfR::getFIX() |> as.data.frame() |> setDT()
    Mutect2_s1  = Mutect2_somatic_vcf |> extract_gt_tidy() |> setDT()
    Mutect2gatk_s21 = Mutect2_somatic_vcf |> extract_info_tidy() |> setDT()
    Mutect2_somatic = cbind(Mutect2_s0[Mutect2_s1$Key, ], Mutect2_s1)
    return(Mutect2_somatic)
}
