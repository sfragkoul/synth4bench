


gt_analysis <- function(runs, folder) {
    
    nt_runs = list()
    
    for(r in runs) {
        
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
    
    gt_runs = nt_runs[which(POS %in% pos_of_interest)]
    
    a <- paste0(folder, "/Merged_report.tsv") |> 
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
    
    
    merged_gt = b[which(POS %in% gt_runs$POS)]
    merged_gt = merged_gt[order(POS)]
    
    merged_gt$Freq = merged_gt$Freq / 100
    
    merged_gt = merged_gt[, by = .(POS, REF, DP), .(
        Nt = paste(Nt, collapse = ","),
        Count = paste(Count, collapse = ","),
        Freq = paste(round(Freq, digits = 3), collapse = ",")
    )]
    
    
    return(merged_gt)
    
}

merge_gatk <- function(gatk_somatic_vcf, merged_gt) {
    
    gatk_s0  = gatk_somatic_vcf |> vcfR::getFIX() |> as.data.frame() |> setDT()
    gatk_s1  = gatk_somatic_vcf |> extract_gt_tidy() |> setDT()
    gatk_s21 = gatk_somatic_vcf |> extract_info_tidy() |> setDT()
    
    gatk_somatic = cbind(gatk_s0[gatk_s1$Key, ], gatk_s1)
    
    
    #Merge everything into a common file-------------------------------------------
    merged_gt$POS = as.character(merged_gt$POS)
    
    merged_bnch = merge(merged_gt, gatk_somatic,  by = "POS", all.x = TRUE)
    
    merged_bnch$POS = as.numeric(merged_bnch$POS)
    
    merged_bnch = merged_bnch[order(POS)]
    
    colnames(merged_bnch) = c(
        "POS",	"Ground Truth REF",	"Ground Truth DP",
        "Ground Truth ALT", "Ground Truth AD", 
        "Ground Truth AF", "CHROM", "ID",	"Mutect2 REF",	
        "Mutect2 ALT", "Mutect2 QUAL",	"Mutect2 FILTER",
        "key", "Indiv", "Mutect2 AD", "Mutect2 AF",
        "Mutect2 DP", "gt_F1R2", "gt_F2R1", "gt_FAD",	
        "gt_GQ", "gt_GT",	"gt_PGT",	"gt_PID",	"gt_PL",
        "gt_PS",	"gt_SB",	"gt_GT_alleles"
    )
    
    return(merged_bnch)
    
}

clean_gatk <- function(df) {
    
    df2 = df[, c(
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
    
    
    
    df2 = df2[, by = c(
        "POS",
        "Ground Truth REF",
        "Ground Truth DP",
        "Mutect2 REF",
        "Mutect2 ALT",
        "Mutect2 DP",
        "Mutect2 AF"

    ), .(
        "Ground Truth ALT" = `Ground Truth ALT` |> tstrsplit(",") |> unlist(),
        "Ground Truth AF"  = `Ground Truth AF` |> tstrsplit(",") |> unlist()
        # "Mutect2 REF" = `Mutect2 REF` |> tstrsplit(",") |> unlist(),
        # "Mutect2 ALT" = `Mutect2 ALT` |> tstrsplit(",") |> unlist(),
        # "Mutect2 DP"  = `Mutect2 DP` |> tstrsplit(",") |> unlist() |> as.integer(),
        # "Mutect2 AF"  = `Mutect2 AF` |> tstrsplit(",") |> unlist() |> as.numeric()
    )]



    mutect2_alt = str_split(df2$`Mutect2 ALT`, ",")
    mutect2_af = str_split(df2$`Mutect2 AF`, ",")


    cln = mapply(
        function(x, y, z) {

            index = which(y == x)

            return(
                c(y[index], z[index])
            )

        },

        df2$`Ground Truth ALT`, mutect2_alt, mutect2_af
    )


    df2$`Mutect2 ALT` = cln |> lapply(function(x) { return(x [1]) }) |> unlist()
    df2$`Mutect2 AF` = cln |> lapply(function(x) { return(x [2]) }) |> unlist()

    df2[which(is.na(`Mutect2 AF`))]$`Mutect2 DP` = NA
    df2[which(is.na(`Mutect2 AF`))]$`Mutect2 REF` = NA

    df2 = df2[, c(
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
    
    return(df2)
    
}
