source("R/libraries.R")
#load GT-----------------------------------------------------------------------
path <- "C:/Users/sfragkoul/Desktop/synth_data/coverage_test/300_30_10"
merged_file <- "Merged"

output_file <- file.path(path,
                         paste0(merged_file, "_snvs_TV.tsv"))

gt_tv <- fread(output_file)

load_gt_report <- function(path, merged_file) {
    #function to load Ground Truth bam-report 
    a <- paste0(path, "/", merged_file, "_report.tsv") |>
        readLines() |>
        str_split(pattern = "\t", simplify = TRUE) |>
        as.data.frame() |> 
        setDT()
    
    a$V1 = NULL
    a$V5 = NULL
    
    colnames(a) = c("POS", "REF", "DP", paste0("ALT_", 1:(ncol(a) - 3)))
    
    a = melt(
        a, id.vars = c("POS", "REF", "DP"),
        variable.factor = FALSE, value.factor = FALSE,
        variable.name = "ALT", value.name = "Count"
    )
    
    a = a[which(Count != "")]
    
    a$POS = as.numeric(a$POS)
    a$DP = as.numeric(a$DP)
    
    a$ALT = str_split_i(a$Count, "\\:", 1)
    
    a$Count = str_split_i(a$Count, "\\:", 2) |>
        as.numeric()
    
    a$Freq = round(a$Count / a$DP, digits = 6)
    
    a = a[order(POS, -Count)]
    
    a = a[which(REF != a$ALT & Count != 0)]
    a$mut = paste(a$POS, #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  a$REF, 
                  a$ALT, sep = ":")
    
    colnames(a) = c("POS", "REF", "DP","ALT", "AD", "Freq", "mut")################
    
    
    # select SNVs
    a_snvs = a[which(ALT %in% c("A", "C", "G", "T")), ]
    #filter DEPTH>2
    #a_snvs = a_snvs[which(a_snvs$Count >2), ]
    
    
    gt = list(
        all = a,
        snvs = a_snvs
        
    )
    return(gt)
}

gt_load <- load_gt_report(path,
                          merged_file)$snvs

# All-TV=noise
gt_load <- gt_load[!mut %in% gt_tv$mut]
#load common functions---------------------------------------------------------
select_snvs <- function(df){
    # select SNVs from caller based on length of REF and ALT
    snvs = df[nchar(df$REF) == nchar(df$ALT)]
    snvs = snvs[which(nchar(snvs$REF) <2), ]
    snvs = snvs[which(nchar(snvs$ALT) <2), ]
    snvs$mut = paste(snvs$POS, snvs$REF, snvs$ALT, sep = ":")
    
    return(snvs)
}
define_fp <- function(caller, gt){
    #FP Variants
    fp_var = caller[which(caller$mut %ni% gt$mut)]
    fp_var$type = "FP"
    
    return(fp_var)
}
define_fn <- function(caller, gt){
    #FN Variants
    fn_var = gt[which(gt$mut %ni% caller$mut)]
    fn_var = fn_var[,c("POS", "REF",  "ALT",  "DP", "AD", "Freq","mut" )]
    colnames(fn_var) = c("POS", "REF",  "ALT",  "DP", "AD", "AF","mut" )
    fn_var$type = "FN"
    
    return(fn_var)
}
define_tp <- function(caller, gt){
    #FN Variants
    tp_var = caller[which(caller$mut %in% gt$mut)]
    tp_var$type = "TP"
    return(tp_var)
}
`%ni%` <- Negate(`%in%`) 


#load caller functions---------------------------------------------------------
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


#build caller noise function---------------------------------------------------
noise_snvs_LoFreq <- function(path, merged_file, gt_load, gt_tv){
    
    LoFreq_somatic <- load_LoFreq_vcf(path, merged_file)
    LoFreq_somatic_snvs <-select_snvs(LoFreq_somatic)
    LoFreq_somatic_snvs <- LoFreq_somatic_snvs[,c("POS", "REF", "ALT", "DP", "AD", "AF" ,"mut" )]
    colnames(LoFreq_somatic_snvs) <- c("POS", "REF",  "ALT",  "DP", "AD", "AF","mut" )
    LoFreq_somatic_snvs$AF = as.numeric(LoFreq_somatic_snvs$AF)######
    LoFreq_somatic_snvs <- LoFreq_somatic_snvs[!mut %in% gt_tv$mut]
    
    fp_var = define_fp(LoFreq_somatic_snvs, gt_load)
    fn_var = define_fn(LoFreq_somatic_snvs, gt_load)
    tp_var = define_tp(LoFreq_somatic_snvs, gt_load)
    
    recall = nrow(tp_var)/(nrow(tp_var) + nrow(fn_var))
    precision = nrow(tp_var)/(nrow(tp_var) + nrow(fp_var))
    
    return(list(
        "fp" = fp_var,
        "fn" = fn_var,
        "tp" = tp_var,
        "noise_recall" = recall,
        "noise_precision" = precision)
    )
}



# test function ---------------------------------------------------------------

noise = noise_snvs_LoFreq(path, merged_file, gt_load, gt_tv)

print(noise$noise_recall)
print(noise$noise_precision)



#venn method-------------------------------------------------------------------
# 
#     vcf_GT = gt_snvs[, c("POS", "REF", "ALT")]
#     vcf_GT$scenario = "GT"
# 
#     vcf_Freebayes = noise_Freebayes[,c("POS", "REF", "ALT")]
#     vcf_Freebayes$scenario = "Freebayes"
# 
#     x = rbind(vcf_GT, vcf_Freebayes)
#     y = x[, c("POS", "REF", "ALT", "scenario"), with = FALSE]
# 
#     y$mut = paste( y$POS, y$REF, y$ALT, sep = ":")
# 
#     y = split(y, y$scenario)
# 
#     y = list(
#         'Ground Truth' = y$GT$mut,
#         'Freebayes'         = y$Freebayes$mut
#     )
# 
#     gr = ggvenn(y, fill_color = c("#43ae8d", "#ae4364")) +
# 
#         coord_equal(clip = "off")
# 
#     gr


















