#'Description: This script calculates the Kandall's tau coefficient of the 
#'Coverage and Read length with the modified Accuracy (mAC), 
#'False negative (FN) and False positive (FP) rates.
#' 
#'Input: an xlsx file with a column containing the callers (named Caller),
#'a column specifying the AC, FN, FP rate (named Variants), and a column for 
#'each coverage (in a Sheet named Coverage)
#'or read length (in a Sheet named Read_length) with the respective rates.
#'
#'Output: A xlsx file with the results in two sheets: Coverage and Read_length
#' 
#'Author: Georgios Karakatsoulis
#' 
library(dplyr)
library(rstatix)
library(ggplot2)

myfile_wd = gsub('Rscripts', '', getwd())

df_coverage = readxl::read_xlsx(paste0(myfile_wd, 'Coverage_Read_length.xlsx'),
                                sheet = 'Coverage') %>%
  
  gather(., key = 'Coverage', value = 'prop', -Caller, -Variants) %>%
  
  mutate(Coverage = gsub('x', '', Coverage) |> as.numeric(),
         ID = paste0(Caller, '_', Variants))

tmp_coverage = split(df_coverage, f = df_coverage$ID)

df_Read_length = readxl::read_xlsx(paste0(myfile_wd, 'Coverage_Read_length.xlsx'),
                                   sheet = 'Read_length') %>%
  
  gather(., key = 'Read_length', value = 'prop', -Caller, -Variants) %>%
  
  mutate(Read_length = Read_length |> as.numeric(),
         ID = paste0(Caller, '_', Variants))

tmp_Read_length = split(df_Read_length, f = df_Read_length$ID)



results = list()

results[['Coverage']] = lapply(tmp_coverage, function(x){
  
  z = cor.test(x$prop, x$Coverage, use = 'pairwise.complete', method = 'kendall')
  
  results = cbind('tau' = z$estimate, 'p' = z$p.value) |> as.data.frame()
  
  return(results)
  
}) %>%
  
  data.table::rbindlist(idcol = 'ID') %>%
  
  mutate(Caller = gsub('_.*', '', ID),
         Variant = gsub('.*_', '', ID)) %>%
  
  select(Caller, Variant, tau, p)


results[['Read_length']] = lapply(tmp_Read_length, function(x){
  
  z = cor.test(x$prop, x$Read_length, use = 'pairwise.complete', method = 'kendall')
  
  results = cbind('tau' = z$estimate, 'p' = z$p.value) |> as.data.frame()
  
  return(results)
  
}) %>%
  
  data.table::rbindlist(idcol = 'ID') %>%
  
  mutate(Caller = gsub('_.*', '', ID),
         Variant = gsub('.*_', '', ID)) %>%
  
  select(Caller, Variant, tau, p)


myfile_save = gsub('Rscripts', 'Results', getwd())

openxlsx::write.xlsx(results, paste0(myfile_save, '/Correlation_results.xlsx'))
