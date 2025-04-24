
read_vcf_LoFreq <- function(path, gt, merged_file) {
  #takes two files and produce a caller vcf file in a certain format 
  vcf <- read.vcfR(paste0(path, "/", merged_file, "_LoFreq_norm.vcf"), verbose = FALSE )
  
  df = (vcf |>
    merge_LoFreq(gt) |>
    clean_LoFreq())
  
  return(df)
  
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
        "POS",	"Ground Truth REF",	"Ground Truth ALT",
        "Ground Truth DP", "Ground Truth AD", "Ground Truth AF", 
        
        "Run", "DP Indiv", "Count Indiv", "Freq Indiv",  "mut",
        
        "CHROM", "ID",	"LoFreq REF",	
        "LoFreq ALT", "LoFreq QUAL", "LoFreq FILTER", 
        "LoFreq DP", "LoFreq AF"
    )
    
    #after unlisting multiple variants in the same position, we must
    # keep only unique FN POS
    merged_bnch <- merged_bnch[, .SD[1], by = POS]
    
    return(
        list(
            "merged_bnch" = merged_bnch,
            "LoFreq_somatic" = LoFreq_somatic)
    )
    
}

clean_LoFreq <- function(df) {
    # Extract relevant columns
    df2 <- df$merged_bnch[, c(
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
    
    # Expand multiallelic GT sites into separate rows
    df2 <- df2[, by = .(POS, `Ground Truth REF`, `Ground Truth DP`), .(
        "Ground Truth ALT" = tstrsplit(`Ground Truth ALT`, ",") |> unlist(),
        "Ground Truth AF"  = tstrsplit(`Ground Truth AF`, ",") |> unlist(),
        "LoFreq REF" = `LoFreq REF`[1],
        "LoFreq ALT" = `LoFreq ALT`[1],
        "LoFreq DP"  = `LoFreq DP`[1],
        "LoFreq AF"  = `LoFreq AF`[1]
    )]
    
    # Match ALT alleles between GT and LoFreq
    df2[, `:=` (
        `ALT Match` = mapply(function(gt_alt, LoFreq_alt) {
            if (is.na(LoFreq_alt)) return(FALSE)
            return(gt_alt %in% unlist(str_split(LoFreq_alt, ",")))
        }, `Ground Truth ALT`, `LoFreq ALT`),
        
        `AF Match` = mapply(function(gt_alt, LoFreq_alt, LoFreq_af) {
            if (is.na(LoFreq_alt) | is.na(LoFreq_af)) return(NA)
            alt_list <- str_split(LoFreq_alt, ",")[[1]]
            af_list  <- str_split(LoFreq_af, ",")[[1]]
            idx <- which(alt_list == gt_alt)
            if (length(idx) == 0) return(NA)
            return(af_list[[idx]])
        }, `Ground Truth ALT`, `LoFreq ALT`, `LoFreq AF`)
    )]
    
    # Keep only matched alleles
    df2[, `LoFreq ALT` := ifelse(`ALT Match`, `Ground Truth ALT`, NA)]
    df2[, `LoFreq AF`  := `AF Match`]
    df2[, `LoFreq REF` := ifelse(`ALT Match`, `LoFreq REF`, NA)]
    df2[, `LoFreq DP`  := ifelse(`ALT Match`, `LoFreq DP`, NA)]
    
    # Classify as TP or FN
    df2[, type := ifelse(is.na(`LoFreq ALT`), "FN", "TP")]
    
    # ΔAF: Caller AF - Ground Truth AF (numeric)
    df2[, `AF Deviation ` := NA_real_]
    df2[type == "TP", `AF Deviation` := as.numeric(`LoFreq AF`) - as.numeric(`Ground Truth AF`)]
    
    # Final output
    df2 <- df2[, .(
        POS,
        `Ground Truth REF`, `Ground Truth ALT`, `Ground Truth DP`, `Ground Truth AF`,
        `LoFreq REF`, `LoFreq ALT`, `LoFreq DP`, `LoFreq AF`,
        type, `AF Deviation`
    )]
    
    
    recall = sum(!is.na(df2$`LoFreq REF`)) / nrow(df2)
    
    return(list(
        "vcf_snvs_cleaned" = df2,
        "recall" = recall))
}

load_LoFreq_vcf <- function(path, merged_file){
    #function to load caller vcf
    LoFreq_somatic_vcf <- read.vcfR( paste0(path, "/", merged_file, 
                                            "_LoFreq_norm.vcf"), verbose = FALSE )
    LoFreq_s0  = LoFreq_somatic_vcf |> vcfR::getFIX() |> as.data.frame() |> setDT()
    #LoFreq_s1  = LoFreq_somatic_vcf |> extract_gt_tidy() |> setDT()
    LoFreq_s2 = LoFreq_somatic_vcf |> extract_info_tidy() |> setDT()
    LoFreq_s2$AD <- sapply( strsplit(LoFreq_s2$DP4, ","), function(x) {
        # ensure we have at least four pieces
        if(length(x) >= 4) {
            sum(as.integer(x[3]), as.integer(x[4]))
        } else {
            NA_integer_
        }
    }
    )
    
    LoFreq_s2 = LoFreq_s2[,c( "DP", "AD","AF" )]
    LoFreq_somatic = cbind(LoFreq_s0, LoFreq_s2)
    return(LoFreq_somatic)
}
