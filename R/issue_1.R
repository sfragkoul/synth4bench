library(data.table)
library(stringr)
library(vcfR)
library(ggplot2)
library(ggforce)
library(ggsci)
library(ggvenn)
library(patchwork)

#Find ALT in Ground Truth------------------------------------------------------    
a <- paste0("results/", "Merged_auto_report.tsv") |>
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

# select SNVs
nt_runs = a[which(Nt %in% c("A", "C", "G", "T")), ]


#Load Caller vcf---------------------------------------------------------------   
#read_vcf_freebayes() function

freebayes_somatic_vcf <- read.vcfR( paste0("results/", 
                                           "Merged_auto_freebayes_norm.vcf"), 
                                            verbose = FALSE )

freebayes_s0  = freebayes_somatic_vcf |> vcfR::getFIX() |> as.data.frame() |> setDT()
freebayes_s1  = freebayes_somatic_vcf |> extract_gt_tidy() |> setDT()
freebayesgatk_s21 = freebayes_somatic_vcf |> extract_info_tidy() |> setDT()
freebayes_somatic = cbind(freebayes_s0[freebayes_s1$Key, ], freebayes_s1)
remove(freebayes_somatic_vcf,freebayes_s0, freebayes_s1, freebayesgatk_s21)

new = freebayes_somatic[which(REF %in% c("A", "C", "G", "T")), ]
new = new[which(ALT %in% c("A", "C", "G", "T")), ]

#FP Variants -----------------------------------------------------------------
`%ni%` <- Negate(`%in%`)

fp_var = freebayes_somatic[which(freebayes_somatic$POS %ni% nt_runs$POS)]

fn_var = freebayes_somatic[which(nt_runs$POS %ni% freebayes_somatic$POS)]


#Compare-----------------------------------------------------------------------

