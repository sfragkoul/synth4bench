library(data.table)
library(stringr)
library(vcfR)
library(ggplot2)
library(ggvenn)

#function "not in" def
`%ni%` <- Negate(`%in%`) 

load_gt_report <- function(path, merged_file) {
    #function to load Ground Truth bam-report 
    a <- paste0(path, merged_file, "_report.tsv") |>
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
    
    a$Freq = round(100 * a$Count / a$DP, digits = 6)
    
    a = a[order(POS, -Count)]
    
    a = a[which(REF != a$ALT & Count != 0)]

    return(a)
}

load_gt_vcf <- function(path, merged_file){
    #function to load Ground Truth vcf
    ground_truth_vcf <- read.vcfR( paste0(path, merged_file, 
                                          "_ground_truth_norm.vcf"),
                                   verbose = FALSE )
    
    ground_truth_vcf  = ground_truth_vcf |> vcfR::getFIX() |> as.data.frame() |> setDT()
    
    pick_gt = nt_runs[which(nt_runs$POS %in% ground_truth_vcf$POS)]
    pick_gt$mut = paste(pick_gt$POS, 
                        pick_gt$REF, 
                        pick_gt$ALT, sep = ":")
    return(pick_gt)
}

load_gatk_vcf <- function(path, merged_file){
    #function to load caller vcf
    Mutect2_somatic_vcf <- read.vcfR( paste0(path, merged_file, 
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
    snvs$mut = paste(snvs$POS, snvs$REF, snvs$ALT, sep = ":")
    
    return(snvs)
}

define_fp <- function(caller, gt){
    #FP Variants
    fp_var = caller[which(caller$mut %ni% gt$mut)]
    
    return(fp_var)
}

define_fn <- function(caller, gt){
    #FN Variants
    fn_var = gt[which(gt$mut %ni% caller$mut)]
    
    return(fn_var)
}

define_tp <- function(caller, gt){
    #FN Variants
    tp_var = caller[which(caller$mut %in% gt$mut)]
    
    return(tp_var)
}

nt_runs = load_gt_report("results/", "Merged_auto")
# select SNVs
nt_runs = nt_runs[which(ALT %in% c("A", "C", "G", "T")), ]
#filter DEPTH>2
nt_runs = nt_runs[which(nt_runs$Count >2), ]

pick_gt = load_gt_vcf("results/", "Merged_auto")
Mutect2_somatic <- load_gatk_vcf("results/", "Merged_auto")
Mutect2_somatic_snvs <-select_snvs(Mutect2_somatic)

fp_var = define_fp(Mutect2_somatic_snvs, pick_gt)
fn_var = define_fn(Mutect2_somatic_snvs, pick_gt)
tp_var = define_tp(Mutect2_somatic_snvs, pick_gt)



# Venn Plot--------------------------------------------------------------------------
# vcf_GT = pick_gt
# vcf_GT$scenario = "GT"
# vcf_GT = vcf_GT[, c("POS", "REF", "ALT", "scenario")]
# 
# 
# vcf_Mutect2 = Mutect2_somatic_snvs
# vcf_Mutect2$scenario = "Mutect2"
# vcf_Mutect2 = vcf_Mutect2[, c("POS", "REF", "ALT", "scenario")]
# 
# 
# x = rbind(vcf_GT, vcf_Mutect2)
# y = x[, c("POS", "REF", "ALT", "scenario"), with = FALSE]
# 
# y$mut = paste(y$POS, y$REF, y$ALT, sep = ":")
# 
# z = y[, by = mut, .(
#     scenarios = paste(scenario, collapse = "+")
# )]
# 
# z = z[order(z$scenarios), ]
# 
# y = split(y, y$scenario)
# 
# 
# y = list(
#     'Ground Truth' = y$GT$mut,
#     'Mutect2'         = y$Mutect2$mut
# )
# 
# ggvenn(y, fill_color = c("#43ae8d", "#ae4364")) +
# 
# coord_equal(clip = "off") +
# 
# ggtitle("SNVs Venn Diagram \n")  +
#     
# theme(
#     plot.title = element_text(size = 18, hjust = 0.5) # Adjust size and center title
# )      

