source("R/libraries.R")

#coverage
folder1 = "D:/sfragkoul/Synth_Data/Synthesizers/NEAT/testing/TP53/coverage_test/"

#SNVs
tp_300 = fread(paste0(folder1, "700_70_10","/Merged_", "LoFreq", "_snvs_TP.tsv"))
tp_300= tp_300[!is.na(tp_300$`LoFreq DP`),]
fp_300 = fread(paste0(folder1, "700_70_10","/Merged_", "LoFreq", "_snvs_FP.tsv"))
fn_300 = fread(paste0(folder1, "700_70_10","/Merged_", "LoFreq", "_snvs_FN.tsv"))

snv_sum = as.numeric(nrow(fn_300)+nrow(tp_300)+nrow(fp_300))

#Indels
indels_tp_300 = fread(paste0(folder1, "700_70_10","/Merged_", "LoFreq", "_indels_TP.tsv"))
indels_fp_300 = fread(paste0(folder1, "700_70_10","/Merged_", "LoFreq", "_indels_FP.tsv"))
indels_fn_300 = fread(paste0(folder1, "700_70_10","/Merged_", "LoFreq", "_indels_FN.tsv"))

indels_sum = as.numeric(nrow(indels_fn_300)+nrow(indels_tp_300)+nrow(indels_fp_300))

print("SNVS ")
print(paste0("Number of FN: ", nrow(fn_300)))
print(paste0("Percentage of FN: ", percent(nrow(fn_300)/snv_sum, accuracy = 0.1)))
print(paste0("Number of TP: ", nrow(tp_300)))
print(paste0("Percentage of TP: ", percent(nrow(tp_300)/snv_sum, accuracy = 0.1)))
print(paste0("Number of FP: ", nrow(fp_300)))
print(paste0("Percentage of FP: ", percent(nrow(fp_300)/snv_sum, accuracy = 0.1)))
print("Indels ")
print(paste0("Number of FN: ", nrow(indels_fn_300)))
print(paste0("Percentage of FN: ", percent(nrow(indels_fn_300)/indels_sum, accuracy = 0.1)))
print(paste0("Number of TP: ", nrow(indels_tp_300)))
print(paste0("Percentage of TP: ", percent(nrow(indels_tp_300)/indels_sum, accuracy = 0.1)))
print(paste0("Number of FP: ", nrow(indels_fp_300)))
print(paste0("Percentage of FP: ", percent(nrow(indels_fp_300)/indels_sum, accuracy = 0.1)))






      