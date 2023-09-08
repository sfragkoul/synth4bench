#'
#'This R script compares the variants that Mutect2 reported against the ground 
#'truth. Firstly, it identifies the variants with 100% Allele Frequency (AF) in 
#'the individual bam files and then calculates their AF in the final Merged 
#'bam file.
#'
#'
#'Input:  bam-readcount tsv reports, vcf file from Mutect2
#'
#'Output: tsv file containing infrtomation regarding the ground truth variants
#'
#'

rm(list = ls())
gc()

library(stringr)
library(data.table)
library(vcfR)
library(ggplot2)

runs = as.character(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))

nt_runs = list()

for(r in runs) {
    
    a <- readLines(
        paste0(r, "/", r, "_bam_report.tsv")
    )
    
    a = str_split(a, pattern = "\t", simplify = TRUE)
    a = a |> as.data.frame() |> setDT()
    
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
    
    a$Nt = str_split(a$Count, "\\:") |>
        lapply(function(x) { return(x[1]) }) |>
        unlist()
    
    a$Count = str_split(a$Count, "\\:") |>
        lapply(function(x) { return(x[2]) }) |>
        unlist() |>
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

rm(a, b, r, runs)



a <- readLines(
    "Merged2_report.tsv"
)

a = str_split(a, pattern = "\t", simplify = TRUE)
a = a |> as.data.frame() |> setDT()

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

a$Nt = str_split(a$Count, "\\:") |>
    lapply(function(x) { return(x[1]) }) |>
    unlist()

a$Count = str_split(a$Count, "\\:") |>
    lapply(function(x) { return(x[2]) }) |>
    unlist() |>
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


rm(a, b)


somatic_vcf <- read.vcfR( paste0("Merged2_GATK.vcf"), verbose = FALSE )

s0  = vcfR::getFIX(somatic_vcf) |> as.data.frame() |> setDT()
s1  = extract_gt_tidy(somatic_vcf) |> setDT()
s21 = extract_info_tidy(somatic_vcf) |> setDT()

somatic = cbind(s0[s1$Key, ], s1)
rm(somatic_vcf, s0, s1, s21)

merged_gt$POS = as.character(merged_gt$POS)

merged_bnch = merge(merged_gt, somatic, by = "POS", all.x = TRUE)

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




fwrite(
    merged_bnch, "Ground_truth_vs_Mutect2.tsv",
    row.names = FALSE, quote = FALSE, sep = "\t"
)

