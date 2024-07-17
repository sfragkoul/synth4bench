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
#------------------------------------------------------------------------------

venn_plot_freebayes <- function(q, p) {
    
    vcf_GT = vcfR::getFIX(q) |> as.data.frame() |> setDT()
    vcf_GT$scenario = "GT"
    
    vcf_freebayes = vcfR::getFIX(p) |> as.data.frame() |> setDT()
    vcf_freebayes$scenario = "Freebayes"
    
    x = rbind(vcf_GT, vcf_freebayes)
    y = x[, c("CHROM", "POS", "REF", "ALT", "scenario"), with = FALSE]
    
    y$mut = paste(y$CHROM, y$POS, y$REF, y$ALT, sep = ":")
    
    y = split(y, y$scenario)
    
    y = list(
        'Ground Truth' = y$GT$mut,
        'Freebayes'         = y$Freebayes$mut
    )
    
    gr = ggvenn(y, fill_color = c("#43ae8d", "#ae8d43")) +
        
        coord_equal(clip = "off")
    
    return(gr)
}

q = read.vcfR( paste0("results/","Merged_auto_ground_truth_norm.vcf"),verbose = FALSE )
p = read.vcfR( paste0("results/","Merged_auto_freebayes_norm.vcf"),verbose = FALSE )


venn = venn_plot_freebayes(q,p)
venn
