source("R/libraries.R")

path <- "C:/Users/sfragkoul/Desktop/synth_data/coverage_test/5000_500_10"
merged_file <- "Merged"
#------------------------------------------------------------------------------
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
                          merged_file)
gt_all<- gt_load$all
gt_snvs<- gt_load$snvs
#------------------------------------------------------------------------------
# load_gt_vcf <- function(path, merged_file, gt_snvs){
#     #function to load Ground Truth vcf
#     ground_truth_vcf <- read.vcfR( paste0(path, "/",merged_file, 
#                                           "_ground_truth_norm.vcf"),
#                                    verbose = FALSE )
#     
#     ground_truth_vcf  = ground_truth_vcf |> vcfR::getFIX() |> as.data.frame() |> setDT()
#     
#     pick_gt = gt_snvs[which(gt_snvs$POS %in% ground_truth_vcf$POS)]
#     pick_gt$mut = paste(pick_gt$POS, 
#                         pick_gt$REF, 
#                         pick_gt$ALT, sep = ":")
#     return(pick_gt)
# }
# 
# pick_gt <- load_gt_vcf(path,
#                        merged_file,
#                        gt_load$snvs)

#------------------------------------------------------------------------------
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

select_snvs <- function(df){
    # select SNVs from caller based on length of REF and ALT
    snvs = df[nchar(df$REF) == nchar(df$ALT)]
    snvs = snvs[which(nchar(snvs$REF) <2), ]
    snvs = snvs[which(nchar(snvs$ALT) <2), ]
    snvs$mut = paste(snvs$POS, snvs$REF, snvs$ALT, sep = ":")
    
    return(snvs)
}

fp_snvs_gatk <- function(Mutect2_somatic_snvs, gt_all){#term snvs is redundant
    #find MUtect2 FP variants
    fp_var = define_fp(Mutect2_somatic_snvs, pick_gt)
    fp_var$gt_AF = as.numeric(fp_var$gt_AF)
    colnames(fp_var) = c("CHROM", "POS","ID", "Mutect2 REF",	
                         "Mutect2 ALT", "Mutect2 QUAL",	"Mutect2 FILTER",
                         "key", "Indiv", "Mutect2 AD", "Mutect2 AF",
                         "Mutect2 DP", "gt_F1R2", "gt_F2R1", "gt_FAD",	
                         "gt_GQ", "gt_GT",	"gt_PGT",	"gt_PID",	"gt_PL",
                         "gt_PS",	"gt_SB",	"gt_GT_alleles", "mut")
    
    #find DP of FP variants'  location in GT
    tmp = gt_all[which(POS %in% unique(fp_var$POS))]
    a = unique(tmp, by = "POS")
    #to include the presence multiple variants in a POS
    index = match(fp_var$POS, a$POS)
    fp_var$`Ground Truth DP` = a[index]$DP
    fp_var$`DP Percentage` = fp_var$`Mutect2 DP`/fp_var$`Ground Truth DP`
    fp_var$type = "FP"
    return(fp_var)
}

define_fp <- function(caller, gt){
    #FP Variants
    fp_var = caller[which(caller$mut %ni% gt$mut)]
    
    return(fp_var)
}


`%ni%` <- Negate(`%in%`) 

final_fp_snvs_gatk <- function(path, merged_file, pick_gt, gt_all){
    
    Mutect2_somatic <- load_gatk_vcf(path, merged_file)
    Mutect2_somatic_snvs <-select_snvs(Mutect2_somatic)
    fp_var <- fp_snvs_gatk(Mutect2_somatic_snvs, pick_gt, gt_all)
    
    return(fp_var)
}

Mutect2_somatic <- load_gatk_vcf(path, merged_file)
Mutect2_somatic_snvs <-select_snvs(Mutect2_somatic)
Mutect2_somatic_snvs <- Mutect2_somatic_snvs[,c("POS", "REF", "ALT", "gt_DP" , "mut" )]

fp_var = define_fp(Mutect2_somatic_snvs, gt_all)



fp_var$gt_AF = as.numeric(fp_var$gt_AF)
colnames(fp_var) = c("CHROM", "POS","ID", "Mutect2 REF",	
                     "Mutect2 ALT", "Mutect2 QUAL",	"Mutect2 FILTER",
                     "key", "Indiv", "Mutect2 AD", "Mutect2 AF",
                     "Mutect2 DP", "gt_F1R2", "gt_F2R1", "gt_FAD",	
                     "gt_GQ", "gt_GT",	"gt_PGT",	"gt_PID",	"gt_PL",
                     "gt_PS",	"gt_SB",	"gt_GT_alleles", "mut")

#find DP of FP variants'  location in GT
tmp = gt_all[which(POS %in% unique(fp_var$POS))]
a = unique(tmp, by = "POS")
#to include the presence multiple variants in a POS
index = match(fp_var$POS, a$POS)
fp_var$`Ground Truth DP` = a[index]$DP
fp_var$`DP Percentage` = fp_var$`Mutect2 DP`/fp_var$`Ground Truth DP`
fp_var$type = "FP"

