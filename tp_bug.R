source("R/libraries.R")

path <- "C:/Users/sfragkoul/Desktop/synth_data/coverage_test/300_30_10"
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

gt_snvs <- load_gt_report(path,
                          merged_file)$snvs

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
#load functions
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

#################################################
# fp_snvs_gatk <- function(Mutect2_somatic_snvs, gt_snvs){#term snvs is redundant
#     #find MUtect2 FP variants
#     fp_var = define_fp(Mutect2_somatic_snvs, gt_snvs)
#     # fp_var$AF = as.numeric(fp_var$AF)######################################
#     fp_var$type = "FP"
#     return(fp_var)
# }

define_fp <- function(caller, gt){
    #FP Variants
    fp_var = caller[which(caller$mut %ni% gt$mut)]
    fp_var$type = "FP"#################################
    
    return(fp_var)
}

define_fn <- function(caller, gt){
    #FN Variants
    fn_var = gt[which(gt$mut %ni% caller$mut)]
    fn_var = fn_var[,c("POS", "REF",  "ALT",  "DP", "AD", "Freq","mut" )]#####
    colnames(fn_var) = c("POS", "REF",  "ALT",  "DP", "AD", "AF","mut" )#####
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

final_fp_snvs_gatk <- function(path, merged_file, gt_snvs){
    
    Mutect2_somatic <- load_gatk_vcf(path, merged_file)
    Mutect2_somatic_snvs <-select_snvs(Mutect2_somatic)
    Mutect2_somatic_snvs[, AD := as.numeric(sapply(strsplit(gt_AD, ","), function(x) x[2]))]#######
    Mutect2_somatic_snvs <- Mutect2_somatic_snvs[,c("POS", "REF", "ALT", "gt_DP", "AD", "gt_AF" ,"mut" )]
    colnames(Mutect2_somatic_snvs) <- c("POS", "REF",  "ALT",  "DP", "AD", "AF","mut" )
    Mutect2_somatic_snvs$AF = as.numeric(Mutect2_somatic_snvs$AF)######
    fp_var <- define_fp(Mutect2_somatic_snvs, gt_snvs)###################
    
    return(fp_var)
}

fp_gatk <- final_fp_snvs_gatk(path, merged_file, gt_snvs)



#see how vars are treated
Mutect2_somatic <- load_gatk_vcf(path, merged_file)
Mutect2_somatic_snvs <-select_snvs(Mutect2_somatic)
Mutect2_somatic_snvs[, AD := as.numeric(sapply(strsplit(gt_AD, ","), function(x) x[2]))]
Mutect2_somatic_snvs <- Mutect2_somatic_snvs[,c("POS", "REF", "ALT", "gt_DP", "AD", "gt_AF" ,"mut" )]
colnames(Mutect2_somatic_snvs) <- c("POS", "REF", "ALT", "DP" , "AD", "AF","mut" )
Mutect2_somatic_snvs$AF = as.numeric(Mutect2_somatic_snvs$AF)
rm(Mutect2_somatic)



fp_var = define_fp(Mutect2_somatic_snvs, gt_snvs)

fn_var = define_fn(Mutect2_somatic_snvs, gt_snvs)

tp_var = define_tp(Mutect2_somatic_snvs, gt_snvs)





















#venn method-------------------------------------------------------------------

    vcf_GT = gt_snvs[, c("POS", "REF", "ALT")]
    vcf_GT$scenario = "GT"
    
    vcf_gatk = Mutect2_somatic_snvs[,c("POS", "REF", "ALT")]
    vcf_gatk$scenario = "GATK"
    
    x = rbind(vcf_GT, vcf_gatk)
    y = x[, c("POS", "REF", "ALT", "scenario"), with = FALSE]
    
    y$mut = paste( y$POS, y$REF, y$ALT, sep = ":")
    
    y = split(y, y$scenario)
    
    y = list(
        'Ground Truth' = y$GT$mut,
        'GATK'         = y$GATK$mut
    )
    
    gr = ggvenn(y, fill_color = c("#43ae8d", "#ae4364")) +
        
        coord_equal(clip = "off")
    

    
    
    















