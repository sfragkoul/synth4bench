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
        "POS",	"Ground Truth REF",	"Ground Truth ALT",
        "Ground Truth DP", "Ground Truth AD", "Ground Truth AF", 
        
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
    
    return(merged_bnch)
    
}

df = merge_freebayes(read.vcfR("C:/Users/sfragkoul/Desktop/synth_data/coverage_test/300_30_10/Merged_Freebayes_norm.vcf", verbose = FALSE ), 
                merged_gt)

#df = merged_bnch$merged_bnch
#freebayes_somatic = df$freebayes_somatic

#-----------------------------------------------------------------------------
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

clean_out = clean_freebayes(df)

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