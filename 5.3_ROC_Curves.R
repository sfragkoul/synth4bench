#'Description: This script produces a ROC curve for a given dataset.
#' 
#'Input: A csv file with variants that were either detected or not detected
#'by each caller.
#' 
#'Output: A ROC curve.
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

df = read.csv(paste0(myfile_wd, '/', myname))

colnames(df)[6] = 'sequences'

df$ref_alt = paste(df$REF, df$ALT)


df$outcome = ifelse(df$Found == 'Yes', 1, 0)


# Interactions

rfs = c('pos_chromosomal', 'pos_read', 'REF', 'ALT', 'sequences')


df1 = df %>%
  
  select(Found, all_of(rfs)) %>%
  
  na.omit() %>%
  
  mutate(
    
    pos_chromosomal = scale(pos_chromosomal),
    pos_read = scale(pos_read)
    
  )


set.seed(123)

training.samples <- df1$Found %>% 
  createDataPartition(p = 0.8, list = FALSE)

train.data  <- df1[training.samples, ]

test.data <- df1[-training.samples, ]

# Dummy code categorical predictor variables
x <- model.matrix(Found~ pos_chromosomal * REF + pos_chromosomal * ALT + 
                    
                    pos_read * REF + pos_read * ALT + 
                    
                    sequences * REF + sequences * ALT, train.data)[,-1]

# Convert the outcome (class) to a numerical variable
y <- train.data$Found

# Find the best lambda using cross-validation
set.seed(123) 

cv.lasso <- cv.glmnet(x, y, alpha = 1, family = "binomial")

# Fit the final model on the training data
model <- glmnet(x, y, alpha = 1, family = "binomial",
                lambda = cv.lasso$lambda.min)


# Make predictions on the test data
x.test <- model.matrix(Found ~ pos_chromosomal * REF + pos_chromosomal * ALT + 
                         
                         pos_read * REF + pos_read * ALT + 
                         
                         sequences * REF + sequences * ALT, test.data)[,-1]

probabilities <- model %>%
  
  predict(newx = x.test) %>%
  
  as.vector()

myroc_int = pROC::roc(
  
  test.data$Found,
  probabilities,
  plot      = F,
  print.auc = T,
  ci        = T
  
)

# No interactions

rfs = c('pos_chromosomal', 'pos_read', 'REF', 'ALT', 'sequences')


df1 = df %>%
  
  select(Found, all_of(rfs)) %>%
  
  na.omit() %>%
  
  mutate(
    
    pos_chromosomal = scale(pos_chromosomal),
    pos_read = scale(pos_read)
    
  )


set.seed(123)

training.samples <- df1$Found %>% 
  createDataPartition(p = 0.8, list = FALSE)

train.data  <- df1[training.samples, ]

test.data <- df1[-training.samples, ]

# Dummy code categorical predictor variables
x <- model.matrix(
  
  Found~ pos_chromosomal + REF + ALT + pos_read + sequences,
  train.data)[,-1]

# Convert the outcome (class) to a numerical variable
y <- train.data$Found

# Find the best lambda using cross-validation
set.seed(123) 

cv.lasso <- cv.glmnet(x, y, alpha = 1, family = "binomial")

# Fit the final model on the training data
model <- glmnet(x, y, alpha = 1, family = "binomial",
                lambda = cv.lasso$lambda.min)

# Make predictions on the test data
x.test <- model.matrix(Found ~ pos_chromosomal + REF + ALT + 
                         
                         pos_read + sequences, test.data)[,-1]

probabilities <- model %>%
  
  predict(newx = x.test) %>%
  
  as.vector()

myroc_no_int = pROC::roc(
  
  test.data$Found,
  probabilities,
  plot      = F,
  print.auc = T,
  ci        = T
  
)



dd1 = data.frame(
  
  myroc_int$sensitivities,
  myroc_int$specificities
  
)

colnames(dd1) = c('sensitivities', 'specificities')

dd1$Group = 'Yes'

auc1 = paste0(
  
  'AUC: ',
  as.numeric(myroc_int$auc) |> round(3),
  ' (',
  round(as.numeric(myroc_int$ci), 3)[1],
  
  ', ',
  round(as.numeric(myroc_int$ci), 3)[3],
  ')'
  
)

dd2 = data.frame(
  
  myroc_no_int$sensitivities,
  myroc_no_int$specificities
  
)

colnames(dd2) = c('sensitivities', 'specificities')

dd2$Group = 'No'

auc2 = paste0(
  
  'AUC: ',
  as.numeric(myroc_no_int$auc) |> round(3),
  ' (',
  round(as.numeric(myroc_no_int$ci), 3)[1],
  
  ', ',
  round(as.numeric(myroc_no_int$ci), 3)[3],
  ')'
  
)


ddd = rbind(dd1, dd2)


p = ggplot(data = ddd, aes(x = 1- specificities, y = sensitivities)) + 
  
  geom_point(aes(colour = Group), size = 1) +
  
  geom_segment(
    aes(x = 0, y = 0, xend = 1, yend = 1),
    linetype = 'dashed', colour = 'darkgrey'
  ) +
  
  annotate("text", x = 0.75, y = 0.45, label = auc1, color = "red") +
  annotate("text", x = 0.75, y = .4, label = auc2, color = "steelblue") +
  
  theme_minimal() +
  
  
  scale_color_manual(values = c('steelblue', 'red')) +
  
  labs(color = 'Interaction terms') +
  
  theme(
    legend.position = 'bottom',
    
    panel.border = element_rect(fill = NA, color = "grey10"),
    
    plot.margin = margin(20, 20, 20, 20),
    
    axis.title.x = element_text(face = "bold", size = 13),
    
    axis.title.y = element_text(face = "bold", size = 13),
    
    legend.title = element_text(face = "bold", size = 13),
    
    legend.text = element_text(size = 11)
    
  )


# ggsave(file = paste0('ROC_', myname, '.svg'),
#        plot = p,
#        width  = 5,
#        height = 5,
#        units = "in"
# )
