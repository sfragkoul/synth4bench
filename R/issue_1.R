library(data.table)
library(stringr)
library(vcfR)
library(ggplot2)
library(ggvenn)

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
#filter DEPTH>2
nt_runs = nt_runs[which(nt_runs$Count >2), ]
remove(a)
#Load GT vcf---------------------------------------------------------------   
ground_truth_vcf <- read.vcfR( paste0("results/",
                                      "Merged_auto_ground_truth_norm.vcf"),
                               verbose = FALSE )

ground_truth_vcf  = ground_truth_vcf |> vcfR::getFIX() |> as.data.frame() |> setDT()

pick_gt = nt_runs[which(nt_runs$POS %in% ground_truth_vcf$POS)]

#Load Caller vcf---------------------------------------------------------------   
#read_vcf_Mutect2() function

Mutect2_somatic_vcf <- read.vcfR( paste0("results/", 
                                           "Merged_auto_Mutect2_norm.vcf"), 
                                            verbose = FALSE )

Mutect2_s0  = Mutect2_somatic_vcf |> vcfR::getFIX() |> as.data.frame() |> setDT()
Mutect2_s1  = Mutect2_somatic_vcf |> extract_gt_tidy() |> setDT()
Mutect2gatk_s21 = Mutect2_somatic_vcf |> extract_info_tidy() |> setDT()
Mutect2_somatic = cbind(Mutect2_s0[Mutect2_s1$Key, ], Mutect2_s1)
remove(Mutect2_somatic_vcf,Mutect2_s0, Mutect2_s1, Mutect2gatk_s21)

# select SNVs from caller based on length of RED and ALT
Mutect2_somatic_snvs = Mutect2_somatic[nchar(Mutect2_somatic$REF) == nchar(Mutect2_somatic$ALT)]
Mutect2_somatic_snvs = Mutect2_somatic_snvs[which(nchar(Mutect2_somatic_snvs$REF) <2), ]

#function "not in" def --------------------------------------------------------
`%ni%` <- Negate(`%in%`)
#Mini example----------------------------------------------------------------
# ground <- c(1,2,3,4)
# call <- c(3,4,5,6)
# 
# fn1 <- ground[which(ground %ni% call)]
# fp2 <- call[which(call %ni% ground)]

#FP and FN Variants -----------------------------------------------------------
fp_var = Mutect2_somatic_snvs[which(Mutect2_somatic_snvs$POS %ni% pick_gt$POS)]

fn_var = pick_gt[which(pick_gt$POS %ni% Mutect2_somatic_snvs$POS)]


#tp_var = pick_gt[which(pick_gt$POS %in% Mutect2_somatic_snvs$POS)]



#Plot--------------------------------------------------------------------------
vcf_GT = pick_gt
vcf_GT$scenario = "GT"
vcf_GT = vcf_GT[, c("POS", "REF", "Nt", "scenario")]
colnames(vcf_GT) = c("POS", "REF", "ALT", "scenario")



vcf_Mutect2 = Mutect2_somatic_snvs
vcf_Mutect2$scenario = "Mutect2"
vcf_Mutect2 = vcf_Mutect2[, c("POS", "REF", "ALT", "scenario")]


x = rbind(vcf_GT, vcf_Mutect2)
y = x[, c("POS", "REF", "scenario"), with = FALSE]

y$mut = paste(y$POS, y$REF, sep = ":")
# 
# z = y[, by = mut, .(
#     scenarios = paste(scenario, collapse = "+")
# )]
# 
# z = z[order(z$scenarios), ]

y = split(y, y$scenario)


y = list(
    'Ground Truth' = y$GT$mut,
    'Mutect2'         = y$Mutect2$mut
)

ggvenn(y, fill_color = c("#43ae8d", "#ae4364")) +

coord_equal(clip = "off")

