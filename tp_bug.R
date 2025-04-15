source("R/libraries.R")

#-----------------------------------------------------------------------------
gt_analysis <- function(runs, folder, merged_file) {
    
    nt_runs = list()
    #ground truth   variants from individual files
    for(r in runs) {
        #folder = "."
        #merged_file = "Merged"
        #process reports.tsv files for individual files
        a <- paste0(folder, "/", r, "/", r, "_report.tsv") |>
            readLines() |>
            str_split(pattern = "\t", simplify = TRUE) |>
            as.data.frame() |> 
            setDT()
        
        a$V1 = NULL
        # a$V3 = NULL
        a$V5 = NULL
        
        colnames(a) = c("POS", "REF", "DP", paste0("Nt_", 1:(ncol(a) - 3)))
        
        a = melt(
            a, id.vars = c("POS", "REF", "DP"),
            variable.factor = FALSE, value.factor = FALSE,
            variable.name = "Nt", value.name = "Count"
        )
        
        a = a[which(Count != "")]
        
        a$POS = as.numeric(a$POS)
        a$DP = as.numeric(a$DP)
        
        a$Nt = str_split_i(a$Count, "\\:", 1)
        
        a$Count = str_split_i(a$Count, "\\:", 2) |>
            as.numeric()
        
        a$Freq = round(100 * a$Count / a$DP, digits = 6)
        
        a = a[order(POS, -Count)]
        
        a = a[which(REF != a$Nt & Count != 0)]
        
        b = a[which(Nt %in% c("A", "C", "G", "T")), ]
        
        nt_runs[[ as.character(r) ]] = b
    }
    
    nt_runs = rbindlist(nt_runs, idcol = "Run")
    
    pos_of_interest = nt_runs[which(Freq == 100)]$POS |> unique()
    
    gt_runs = nt_runs[POS %in% pos_of_interest & Freq == "100"] #!!!NEW
    
    #same process reports.tsv files for Merged file
    a <- paste0(folder, "/", merged_file , "_report.tsv") |> 
        readLines() |>
        str_split(pattern = "\t", simplify = TRUE) |> 
        as.data.frame() |> 
        setDT()
    
    a$V1 = NULL
    # a$V3 = NULL
    a$V5 = NULL
    
    colnames(a) = c("POS", "REF", "DP", paste0("Nt_", 1:(ncol(a) - 3)))
    
    a = melt(
        a, id.vars = c("POS", "REF", "DP"),
        variable.factor = FALSE, value.factor = FALSE,
        variable.name = "Nt", value.name = "Count"
    )
    
    a = a[which(Count != "")]
    
    a$POS = as.numeric(a$POS)
    a$DP = as.numeric(a$DP)
    
    a$Nt = str_split_i(a$Count, "\\:", 1)
    
    a$Count = str_split_i(a$Count, "\\:", 2) |>
        as.numeric()
    
    a$Freq = round(100 * a$Count / a$DP, digits = 6)
    
    a = a[order(POS, -Count)]
    
    a = a[which(REF != a$Nt & Count != 0)]
    
    b = a[which(Nt %in% c("A", "C", "G", "T")), ]
    
    
    #merged_gt = b[which(POS %in% gt_runs$POS)]
    merged_gt <- merge(b, gt_runs, by = c("POS", "REF", "Nt")) #ΝEW!!!!!!!
    colnames(merged_gt) = c("POS", "REF", "ALT", "DP", "Count", "Freq",
                             "Run", "DP Indiv", "Count Indiv", "Freq Indiv")
    
    merged_gt = merged_gt[order(POS)]
    
    merged_gt$Freq = merged_gt$Freq / 100
    
    # merged_gt1 = merged_gt[, by = .(POS, REF, DP), .(
    #     Nt = paste(Nt, collapse = ","),
    #     Count = paste(Count, collapse = ","),
    #     Freq = paste(round(Freq, digits = 3), collapse = ",")
    # )]
    
    return(merged_gt)
    
}
merged_gt = gt_analysis(c(1,2,3,4,5,6,7,8,9,10),
                        "C:/Users/sfragkoul/Desktop/synth_data/coverage_test/300_30_10",
                        "Merged")


#-----------------------------------------------------------------------------

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
        
        "Run", "DP Indiv", "Count Indiv", "Freq Indiv", 
        
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

df = merge_VarDict(read.vcfR("C:/Users/sfragkoul/Desktop/synth_data/coverage_test/300_30_10/Merged_VarDict_norm.vcf", verbose = FALSE ), 
                merged_gt)

#df = merged_bnch$merged_bnch
#freebayes_somatic = df$freebayes_somatic

#-----------------------------------------------------------------------------
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

clean_out = clean_VarDict(df)

df_cleaned = clean_out$vcf_snvs_cleaned
recall = clean_out$recall




#------------------------------------------------
# load_vcf <- function(vcf){
#     freebayes_s0  = freebayes_somatic_vcf |> vcfR::getFIX() |> as.data.frame() |> setDT()
#     freebayes_s1  = freebayes_somatic_vcf |> extract_gt_tidy() |> setDT()
#     freebayes_s21 = freebayes_somatic_vcf |> extract_info_tidy() |> setDT()
#     freebayes_somatic = cbind(freebayes_s0[freebayes_s1$Key, ], freebayes_s1)
# 
#     return(freebayes_somatic)
# }
# 
# freebayes_somatic_vcf = read.vcfR("C:/Users/sfragkoul/Desktop/synth_data/coverage_test/5000_500_10/Merged_Freebayes_norm.vcf")
# freebayes_somatic_vcf = load_vcf(freebayes_somatic_vcf)