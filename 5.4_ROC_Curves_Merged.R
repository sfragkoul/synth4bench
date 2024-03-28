#'Description: This script produces a merged figure with the ROC curves for 
#'a given caller, each curve corresponding to a specific Coverage and 
#'Read length.
#'
#'Input: A folder with the files containing all the datasets and the 
#'"5.3_ROC_Curves.R" file.
#'
#'Output: A figure with the ROC curves for each caller.
#'
#'Author: Georgios Karakatsoulis
#'

library(dplyr)
library(rstatix)
library(ggplot2)
library(pROC)
library(cutpointr)
library(glmnet)
library(caret)
library(patchwork)

myfile_wd = gsub('Rscripts', 'Desktop', getwd())


myfiles = list.files(myfile_wd)

mycallers = gsub('.*All_', '', myfiles)
mycallers = gsub('_.*', '', mycallers)
mycallers = unique(mycallers) |> toupper()

for (caller in mycallers){
  
  myfiles_tmp = myfiles[grepl(caller, toupper(myfiles))]
  
  for(i in myfiles_tmp){
    
    myname = i
    
    source('Report_AUC_Curves_Clean.R')
    
    assign(
      
      x = paste0('p_',gsub('_All.*', '', i)),
      value = p
    )
    
  }
  
  
  p_cov_300 = p_300_30_10 + annotate('text', x = 0.15, y = 0.9, label = 'Coverage: 300x', color = 'black')
  p_cov_1000 = p_1000_100_10 + annotate('text', x = 0.15, y = 0.9, label = 'Coverage: 1000x', color = 'black')
  p_cov_3000 = p_3000_300_10 + annotate('text', x = 0.15, y = 0.9, label = 'Coverage: 3000x', color = 'black')
  
  p_read_50 = p_1000_50 + annotate('text', x = 0.15, y = 0.9, label = 'Read length: 50', color = 'black')
  p_read_100 = p_1000_100 + annotate('text', x = 0.15, y = 0.9, label = 'Read length: 100', color = 'black')
  p_read_300 = p_1000_300 + annotate('text', x = 0.15, y = 0.9, label = 'Read length: 300', color = 'black')
  
  
  global_p = (p_cov_300 | p_cov_1000 | p_cov_3000) /
    (p_read_50 | p_read_100 | p_read_300) +
    plot_layout(guides = 'collect') &
    theme(legend.position = "bottom")
  
  
  myfile_save = gsub('Rscripts', 'Results/ROC_Curves', getwd())
  
  ggsave(paste0(myfile_save, '/', caller, '_ROC_Paper.jpeg'), global_p, width = 16, height = 10, units = "in", dpi = 600)
  
  
  
}
