source("R/libraries.R")
#load GT-----------------------------------------------------------------------
`%ni%` <- Negate(`%in%`) 

path <- "C:/Users/sfragkoul/Desktop/synth_data/coverage_test/700_70_10"
merged_file <- "Merged"

load_gt_report_indels <- function(path, merged_file) {
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
    
    # select indels
    a_indels = a[which(ALT %ni% c("A", "C", "G", "T")), ]
    #filter DEPTH>2
    #a_indels = a_indels[which(a_indels$Count >2), ]
    
    gt = list(
        all = a,
        indels = a_indels
        
    )
    return(gt)
}

standardize_indels <- function(dt) {
    #function to standardize indels
    setDT(dt)
    
    #deletions
    dt[grepl("^-", ALT), `:=` (
        ALT = substring(REF, 1, 1), 
        REF = paste0(REF, substring(ALT, 2)),
        POS = POS - 1  #Adjust POS for deletions
    )]
    
    #insertions
    dt[grepl("^\\+", ALT), ALT := paste0(REF, substring(ALT, 2))]
    
    dt$mut = paste(dt$POS, 
                   dt$REF, 
                   dt$ALT, sep = ":")
    return(dt)
}

gt_stdz_indels <- function(path, merged_file){
    #gt_all = load_gt_report_indels(path, merged_file)$all
    gt_indels = load_gt_report_indels(path, merged_file)$indels
    #pick_gt = load_gt_vcf_indels(path, merged_file, gt_indels)
    pick_gt_stdz = standardize_indels(gt_indels)
    colnames(pick_gt_stdz) = c("POS", "REF", "DP", "ALT", "AD", "Freq", "mut" )
    return(pick_gt_stdz)
} 


pick_gt_stdz <- gt_stdz_indels(path, merged_file)

#load common functions---------------------------------------------------------
select_indels <- function(df){
    #function to select indels from caller based on length of REF and ALT
    
    #identify indels based on length
    indels = df[nchar(df$REF) != nchar(df$ALT)]
    indels$mut = paste(indels$POS, indels$REF, indels$ALT, sep = ":")
    
    return(indels)
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

#load caller functions---------------------------------------------------------
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

#build caller function---------------------------------------------------------

call_indels_gatk <- function(path, merged_file, pick_gt_stdz){
    
    Mutect2_somatic <- load_gatk_vcf(path, merged_file)
    Mutect2_indels <-select_indels(Mutect2_somatic)
    Mutect2_indels[, AD := as.numeric(sapply(strsplit(gt_AD, ","), function(x) x[2]))]#######
    Mutect2_indels <- Mutect2_indels[,c("POS", "REF", "ALT", "gt_DP", "AD", "gt_AF" ,"mut" )]
    Mutect2_indels <- Mutect2_indels[,c("POS", "REF", "ALT", "gt_DP", "AD", "gt_AF" ,"mut" )]
    colnames(Mutect2_indels) <- c("POS", "REF",  "ALT",  "DP", "AD", "AF","mut" )
    Mutect2_indels$AF = as.numeric(Mutect2_indels$AF)
    Mutect2_indels$POS = as.numeric(Mutect2_indels$POS)
    
    #TP
    tp_var = define_tp(Mutect2_indels, pick_gt_stdz)
    tp_var$category = NaN
    
    #FN
    fn_var = define_fn(Mutect2_indels, pick_gt_stdz)
    #categorize INDELs
    ##Same POS
    same_POS <- merge(fn_var, Mutect2_indels, by = c("POS"))
    fn_var[, category := ifelse(POS %in% same_POS$POS, "diff REF", "not exist")]
    ##Same POS & REF
    same_POS_REF <- merge(fn_var, Mutect2_indels, by = c("POS", "REF"))
    ##Update only rows where POS and REF match
    fn_var[POS %in% same_POS_REF$POS & REF %in% same_POS_REF$REF, 
           category := "diff ALT"]
    
    
    
    #FP
    fp_var = define_fp(Mutect2_indels, pick_gt_stdz)
    #categorize INDELs
    #Same POS
    same_POS <- merge(fp_var, pick_gt_stdz, by = c("POS"))
    fp_var[, category := ifelse(POS %in% same_POS$POS, "diff REF", "not exist")]
    
    #Same POS & REF
    same_POS_REF <- merge(fp_var, pick_gt_stdz, by = c("POS", "REF"))
    # Update only rows where POS and REF match
    fp_var[POS %in% same_POS_REF$POS & REF %in% same_POS_REF$REF, 
                   category := "diff ALT"]
   
    
    recall = nrow(tp_var)/(nrow(tp_var) + nrow(fn_var))
    precision = nrow(tp_var)/(nrow(tp_var) + nrow(fp_var))
    
    return(list(
        "fp" = fp_var,
        "fn" = fn_var,
        "tp" = tp_var,
        "indel_recall" = recall,
        "indel_precision" = precision)
    )
}


# test function ---------------------------------------------------------------

indels <- call_indels_gatk(path, merged_file, pick_gt_stdz)


indels_out = rbind(indels$tp, indels$fp, indels$fn)

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


















