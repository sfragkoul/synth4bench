
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
        "POS",	"Ground Truth REF",	"Ground Truth ALT",
        "Ground Truth DP", "Ground Truth AD", "Ground Truth AF", 
        
        "Run", "DP Indiv", "Count Indiv", "Freq Indiv", "mut",
        
        "CHROM", "ID", "VarScan REF",	
        "VarScan ALT", "VarScan QUAL",	"VarScan FILTER", "VarScan DP", "Pvalue","VarScan AF"
    )
    
    #after unlisting multiple variants in the same position, we must
    # keep only unique FN POS
    merged_bnch <- merged_bnch[, .SD[1], by = POS]
    
    return(
        list(
            "merged_bnch" = merged_bnch,
            "VarScan_somatic" = VarScan_somatic)
    )
    
}


clean_VarScan <- function(df) {
    # Extract relevant columns
    df2 <- df$merged_bnch[, c(
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
    
    # Expand multiallelic GT sites into separate rows
    df2 <- df2[, by = .(POS, `Ground Truth REF`, `Ground Truth DP`), .(
        "Ground Truth ALT" = tstrsplit(`Ground Truth ALT`, ",") |> unlist(),
        "Ground Truth AF"  = tstrsplit(`Ground Truth AF`, ",") |> unlist(),
        "VarScan REF" = `VarScan REF`[1],
        "VarScan ALT" = `VarScan ALT`[1],
        "VarScan DP"  = `VarScan DP`[1],
        "VarScan AF"  = `VarScan AF`[1]
    )]
    
    # Match ALT alleles between GT and VarScan
    df2[, `:=` (
        `ALT Match` = mapply(function(gt_alt, gatk_alt) {
            if (is.na(gatk_alt)) return(FALSE)
            return(gt_alt %in% unlist(str_split(gatk_alt, ",")))
        }, `Ground Truth ALT`, `VarScan ALT`),
        
        `AF Match` = mapply(function(gt_alt, gatk_alt, gatk_af) {
            if (is.na(gatk_alt) | is.na(gatk_af)) return(NA)
            alt_list <- str_split(gatk_alt, ",")[[1]]
            af_list  <- str_split(gatk_af, ",")[[1]]
            idx <- which(alt_list == gt_alt)
            if (length(idx) == 0) return(NA)
            return(af_list[[idx]])
        }, `Ground Truth ALT`, `VarScan ALT`, `VarScan AF`)
    )]
    
    # Keep only matched alleles
    df2[, `VarScan ALT` := ifelse(`ALT Match`, `Ground Truth ALT`, NA)]
    df2[, `VarScan AF`  := `AF Match`]
    df2[, `VarScan REF` := ifelse(`ALT Match`, `VarScan REF`, NA)]
    df2[, `VarScan DP`  := ifelse(`ALT Match`, `VarScan DP`, NA)]
    
    # Classify as TP or FN
    df2[, type := ifelse(is.na(`VarScan ALT`), "FN", "TP")]
    
    # ΔAF: Caller AF - Ground Truth AF (numeric)
    df2[, `AF Deviation ` := NA_real_]
    df2[type == "TP", `AF Deviation` := as.numeric(`VarScan AF`) - as.numeric(`Ground Truth AF`)]
    
    # Final output
    df2 <- df2[, .(
        POS,
        `Ground Truth REF`, `Ground Truth ALT`, `Ground Truth DP`, `Ground Truth AF`,
        `VarScan REF`, `VarScan ALT`, `VarScan DP`, `VarScan AF`,
        type, `AF Deviation`
    )]
    
    
    recall = sum(!is.na(df2$`VarScan REF`)) / nrow(df2)
    
    return(list(
        "vcf_snvs_cleaned" = df2,
        "recall" = recall))
}

load_VarScan_vcf <- function(path, merged_file){
    #function to load caller vcf
    VarScan_somatic_vcf <- read.vcfR( paste0(path, "/", merged_file, 
                                             "_VarScan_norm.vcf"), verbose = FALSE )
    VarScan_s0  = VarScan_somatic_vcf |> vcfR::getFIX() |> as.data.frame() |> setDT()
    #VarScan_s1  = VarScan_somatic_vcf |> extract_gt_tidy() |> setDT()
    VarScan_s2 = VarScan_somatic_vcf |> extract_info_tidy() |> setDT()
    VarScan_s2 = VarScan_s2[,c( "DP", "Strands2", "AF" )]
    VarScan_somatic = cbind(VarScan_s0, VarScan_s2)
    return(VarScan_somatic)
}
