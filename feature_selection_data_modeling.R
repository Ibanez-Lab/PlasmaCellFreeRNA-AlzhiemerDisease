# Packages
library(dplyr)
library (tibble)
library (ggplot2)
library (pROC)
library (reshape)
library(readr)
library(readxl)
library(doParallel)
library(stringr)
library(plyr)
library(tidyverse)
library(data.table)
library(glmnet)
library(InformationValue)
library(entropy)
library(ggcorrplot)
library(lares)
library(corrr)
library(ggpubr)
library(kableExtra)
library(scales)
library(beeswarm)
library(ggbeeswarm)
library(enrichR)
library(biomaRt)

# Read files 
discovery = fread("k99presimnotion.tsv") # Training
replication = fread( "ExpansionALLPSvsCOallgenes.tsv" ) # Testing



# Get common genes
inner_dplyr <- inner_join(discovery, replication, by = "V1")
inner_dplyr.menores <- subset(inner_dplyr , inner_dplyr$log2FoldChange.x < 0 & inner_dplyr$log2FoldChange.y < 0)
inner_dplyr.mayores <- subset(inner_dplyr , inner_dplyr$log2FoldChange.x > 0 & inner_dplyr$log2FoldChange.y > 0)
inner_dplyr = rbind(inner_dplyr.menores, inner_dplyr.mayores)

# As a strategy to select common genes from independent datasets to feed Machine-Learning algorithms we computed the Kullback–Leibler divergence (KL)
KLgenedf = function(df1, df2, numBins = 5){
  
  # Order an obtain the same genes  
  df1names =NULL
  df1names = colnames(df1)
  
  df1 <- df1[, df1names]
  df2 <- df2[, df1names]
  
  # Decide range and compute
  df1range = NULL
  df2range = NULL
  
  
  dfkldif = data.frame(matrix(0, nrow = 1, ncol = length(df1names),
                              dimnames = list(NULL, df1names)) )
  
  for( gene in df1names) {
    df1range = NULL
    df1range = range(df1[,gene])
    df1rangemin = df1range[1]
    df1rangemax = df1range[2]
    
    
    df2range = range(df2[,gene])
    df2rangemin = df2range[1]
    df2rangemax = df2range[2]
    
    # Obtain range
    minrange = ifelse(df1rangemin < df2rangemin, df1rangemin, df2rangemin)
    maxrange = ifelse(df1rangemax > df2rangemax, df1rangemax, df2rangemax)
    ourrange = c(minrange,maxrange)
    
    # Discretize gene distribution
    disdf1 = discretize(df1[,gene], numBins = numBins, r = ourrange)
    disdf2 = discretize(df2[,gene], numBins = numBins, r = ourrange)
    
    
    pdisdf1 = disdf1/sum(disdf1)
    pdisdf2 = disdf2/sum(disdf2)
    
    # Sum kl
    kldf1 = KL.plugin(pdisdf1,pdisdf2)
    kldf2 = KL.plugin(pdisdf2, pdisdf1)
    kldif = kldf1+kldf2
    
    dfkldif[,gene] = kldif
    
  }
  
  return(dfkldif)
  
}


# Rlog matrix discovery
matrix_phase = readRDS("matrix_rlogk99qcADPreclinicalandControls_sex_agep_notion")

# Rlog matrix replication
test = readRDS("rlogmatrix_allpsold38qcp_notion")



# Select common gene names
naboonames = colnames(matrix_phase) %in% inner_dplyr$V1
matrix_phase = as.data.frame(matrix_phase)
naboo = matrix_phase[,naboonames]


oldnames = colnames(test) %in% inner_dplyr$V1
old = test[,oldnames]



# Z score and KL
naboo <- as.data.frame(lapply(naboo, function(x) (x - mean(x))/sd(x) ))
old <- as.data.frame(lapply(old, function(x) (x - mean(x))/sd(x) ))


kl_naboo_old = KLgenedf(old,naboo)



# KL to select genes
klsd = seq(0.01, 0.50, by= 0.01)

# klsdr = seq(1, 10, by= 1)

bestacu = as.data.frame(klsd)
bestacu$accuracy = 0
bestacu$n_genes = 0
bestacu$klsd = as.numeric(bestacu$klsd)
# bestacu$nat =  klsdr

# kltr = 0.06


for(i in  1:50){
  
  
  kltr =  bestacu[i,1]
  
  difnaboold = kl_naboo_old[,kl_naboo_old < kltr]
  
  
  genestodoml = colnames(difnaboold)
  
  
  bestacu[i,3] = length(genestodoml)
  
  
  # Load Matrix
  matrix_phase = readRDS("matrix_rlogk99qcADPreclinicalandControls_sex_agep_notion")
  train = matrix_phase
  
  
  # Select features
  genestodoml = c(genestodoml, "Status")
  commontrain = colnames(train) %in% genestodoml
  train = train[,..commontrain]
  
  
  
  # Select only genes
  x_train <- train[,-1]
  
  # Compute Z score
  x_train <- as.data.frame(lapply(x_train, function(x) (x - mean(x))/sd(x) ))
  x_train <- as.matrix(x_train )
  
  # Code status
  numeros = as.numeric(train$Status)
  y_train <-replace(numeros, numeros == 2,0)
  
  
  # Set seed
  set.seed(123)
  
  # Train ridge
  cv_error <- cv.glmnet(
    x      = x_train,
    y      = factor(y_train),
    alpha  = 0,
    nfolds = 5, #
    type.measure = "mse",
    standardize  = F,
    family = "binomial",
    intercept = F
  )
  
  
  
  # Best model lambda
  modelo <- glmnet(
    x           = x_train,
    y           = y_train,
    alpha       = 0,
    lambda      = cv_error$lambda.min,
    standardize = F,
    family = "binomial",
    intercept = F
  )
  
  
  
  # Training predictions
  # ==============================================================================
  predicciones_train <- predict(modelo, newx = x_train, type = "response", standardize = F)
  
  
  # MSE (Training)
  # ==============================================================================
  training_mse <- mean((predicciones_train - y_train)^2)
  
  
  
  # Load Matrix
  test = readRDS("rlogmatrix_allpsold38qcp_notion")
  
  
  # Select features
  commontest = colnames(test) %in% genestodoml
  test = test[,commontest]
  
  
  # Select only genes
  x_test <- test[,-1]
  
  # Compute Z score
  x_test <- as.data.frame(lapply(x_test, function(x) (x - mean(x))/sd(x) ))
  x_test <- as.matrix(x_test )
  
  # Code status
  numeros = as.numeric(test$Status)
  y_test <-replace(numeros, numeros == 1,0)
  y_test <-replace(y_test, y_test == 2,1)
  # y_test <-replace(numeros, numeros == 1,0)
  # y_test <-replace(y_test, y_test == 2,1)
  
  
  
  # Testing preditions
  # ==============================================================================
  predicciones_test <- predict(modelo, newx = x_test, type = "response" , standarize = F, intercept = F)
  
  
  # MSE (testing)
  # ==============================================================================
  test_mse_lasso <- mean((predicciones_test - y_test)^2)
  
  
  
  
  # Confusion Matrix
  testconfu = caret::confusionMatrix(data = factor(as.numeric(predicciones_test>0.5)), reference = factor(y_test))
  
  accuracytest = testconfu$overall[1]
  
  bestacu[i,2] = accuracytest[[1]]
  
  print(bestacu)
  
}

write.csv(bestacu, "bestacuprs.csv", row.names = F)

# Order by best accuracy
bestacuprs = fread("bestacuprs.csv")

bestacuprs <- bestacuprs[order(bestacuprs$accuracy, decreasing = T ),]


bestacuprs$n_genes = as.factor(bestacuprs$n_genes)
bestacuprs$KL = bestacuprs$klsd

klplot = ggplot(bestacuprs, aes(x = n_genes, y = accuracy)) + 
  geom_point(aes(color = KL), size = 3) + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  + xlab("Number of genes in the predictive model") + ylab("Accuracy") + scale_color_gradientn(colours = rainbow(5)) + labs(title="Different KL tresholds effect on testing accuracy and the number of genes in the predictive model")



ggsave(filename = "klplotfigure.png",
       plot = klplot,
       device = "png",
       width =12,
       height = 7,
       dpi = 300)

klplot



# Training with the KL limit do you want
difnaboold = kl_naboo_old[,kl_naboo_old < 0.27] # 1281 genes


# 1281 genes 
genestodoml = colnames(difnaboold)

# Load Matrix
matrix_phase = readRDS("matrix_rlogk99qcADPreclinicalandControls_sex_agep_notion")
train = matrix_phase


# Select features
genestodoml = c(genestodoml, "Status")
commontrain = colnames(train) %in% genestodoml
train = train[,..commontrain]



# Select only genes
x_train <- train[,-1]

# Compute Z score
x_train <- as.data.frame(lapply(x_train, function(x) (x - mean(x))/sd(x) ))
x_train <- as.matrix(x_train )

# Code status  
numeros = as.numeric(train$Status)
y_train <-replace(numeros, numeros == 2,0)


# Set seed
set.seed(123)

# Train ridge
cv_error <- cv.glmnet(
  x      = x_train,
  y      = factor(y_train),
  alpha  = 0,
  nfolds = 5, # 
  type.measure = "mse", 
  standardize  = F,
  family = "binomial",
  intercept = F
)



# Best model lambda
modelo <- glmnet(
  x           = x_train,
  y           = y_train,
  alpha       = 0,
  lambda      = cv_error$lambda.min, 
  standardize = F,
  family = "binomial",
  intercept = F
)


# Compute coeficients
df_coeficientes <- coef(modelo) %>%
  as.matrix() %>%
  as_tibble(rownames = "predictor")


# Training predictions
# ==============================================================================
predicciones_train <- predict(modelo, newx = x_train, type = "response", standardize = F)


# misClassError(y_train, predicciones_train)

# MSE (Training)
# ==============================================================================
training_mse <- mean((predicciones_train - y_train)^2)
# paste("Train error (mse):", training_mse)



# plot(predicciones_train, y_train, xlab = "Predicted", ylab = "Actual",
#      main = "Classification")
# 
# plotROC(y_train, predicciones_train)


# Confusion Matrix
caret::confusionMatrix(data = factor(as.numeric(predicciones_train>0.5)), reference = factor(y_train)) 


df_coeficientestraning_rlog = subset(df_coeficientes, df_coeficientes$s0 != 0)



# testing on ps

# Load Matrix
test = readRDS("rlogmatrix_allpsold38qcp_notion")


# Select features
commontest = colnames(test) %in% genestodoml
test = test[,commontest]


# Select only genes
x_test <- test[,-1]

# Compute Z score
x_test <- as.data.frame(lapply(x_test, function(x) (x - mean(x))/sd(x) ))
x_test <- as.matrix(x_test )

# Code status
numeros = as.numeric(test$Status)
y_test <-replace(numeros, numeros == 1,0)
y_test <-replace(y_test, y_test == 2,1)



# Testing preditions
# ==============================================================================
predicciones_test <- predict(modelo, newx = x_test, type = "response" , standarize = F, intercept = F) 


# MSE (testing)
# ==============================================================================
test_mse_lasso <- mean((predicciones_test - y_test)^2)

# paste("Test error (mse)", test_mse_lasso)
# plot(predicciones_test, y_test, xlab = "Predicted", ylab = "Actual",
#      main = "Probabilities of each class")
# misClassError(y_test, predicciones_test)
# plotROC(y_test, predicciones_test)


# Confusion Matrix
caret::confusionMatrix(data = factor(as.numeric(predicciones_test>0.5)), reference = factor(y_test))


# testing on CDR 0.5


difnaboold = kl_naboo_old[,kl_naboo_old < 0.27]


# 1281 genes 
genestodoml = colnames(difnaboold)
genestodoml = c(genestodoml, "Status")
# Load Matrix
test = readRDS("rlogmatrix_allCDR05old38_notion")


genesgenestodoml = as.data.frame(genestodoml)
colnames(genesgenestodoml) = "gene"

# Select features
commontest = colnames(test) %in% genestodoml
test = test[,commontest]

nombrestest = as.data.frame(colnames(test))
colnames(nombrestest) = "gene"


genesnoaparecen = anti_join(genesgenestodoml,nombrestest) # 
genesnoaparecen = as.character(genesnoaparecen$gene)
length(genesnoaparecen)


test = readRDS("rlogmatrix_allCDR05old38_notion")



commontest = colnames(test) %in% genestodoml
test = test[,commontest]


# Select only genes
x_test <- test[,-1]



# Compute Z score
x_test <- as.data.frame(lapply(x_test, function(x) (x - mean(x))/sd(x) ))

# add gene name from cdr05
x_test = cbind(x_test, setNames( lapply(genesnoaparecen, function(x) x=0), genesnoaparecen) )
x_test <- as.matrix(x_test )



# Code status
numeros = as.numeric(test$Status)
y_test <-replace(numeros, numeros == 2,0)



# Testing preditions
# ==============================================================================
predicciones_test <- predict(modelo, newx = x_test, type = "response" , standarize = F, intercept = F) 


# MSE (testing)
# ==============================================================================
test_mse_lasso <- mean((predicciones_test - y_test)^2)

# paste("Test error (mse)", test_mse_lasso)
# plot(predicciones_test, y_test, xlab = "Predicted", ylab = "Actual",
#      main = "Probabilities of each class")
# misClassError(y_test, predicciones_test)
# plotROC(y_test, predicciones_test)


# Confusion Matrix
caret::confusionMatrix(data = factor(as.numeric(predicciones_test>0.5)), reference = factor(y_test))









kltr =  0.31

difnaboold = kl_naboo_old[,kl_naboo_old < kltr]


genestodoml = colnames(difnaboold)



# Load Matrix
matrix_phase = readRDS("matrix_rlogk99qcADPreclinicalandControls_sex_agep_notion")
train = matrix_phase


# Select features
genestodoml = c(genestodoml, "Status")
commontrain = colnames(train) %in% genestodoml
train = train[,..commontrain]



# Select only genes
x_train <- train[,-1]

# Compute Z score
x_train <- as.data.frame(lapply(x_train, function(x) (x - mean(x))/sd(x) ))
x_train <- as.matrix(x_train )

# Code status
numeros = as.numeric(train$Status)
y_train <-replace(numeros, numeros == 2,0)




# Set seed
set.seed(123)

# Train ridge
cv_error <- cv.glmnet(
  x      = x_train,
  y      = factor(y_train),
  alpha  = 0,
  nfolds = 5, #
  type.measure = "mse",
  standardize  = F,
  family = "binomial",
  intercept = F
)



# Best model lambda
modelo <- glmnet(
  x           = x_train,
  y           = y_train,
  alpha       = 0,
  lambda      = cv_error$lambda.min,
  standardize = F,
  family = "binomial",
  intercept = F
)

df_coeficientes <- coef(modelo) %>%
  as.matrix() %>%
  as_tibble(rownames = "predictor")


# Training predictions
# ==============================================================================
predicciones_train <- predict(modelo, newx = x_train, type = "response", standardize = F)


# MSE (Training)
# ==============================================================================
training_mse <- mean((predicciones_train - y_train)^2)





df_coeficientestraning_rlog = subset(df_coeficientes, df_coeficientes$s0 != 0)


df_coeficientestraning_rlog$s0 = abs(df_coeficientestraning_rlog$s0)



df_coeficientestraning_rlog <- df_coeficientestraning_rlog[ order(df_coeficientestraning_rlog$s0, decreasing = T),]



# 90 genes
df_coeficientestraning_rlogp = df_coeficientestraning_rlog[1:90,]
genestodoml = df_coeficientestraning_rlogp$predictor



# Load Matrix
matrix_phase = readRDS("matrix_rlogk99qcADPreclinicalandControls_sex_agep_notion")
train = matrix_phase


# Select features
genestodoml = c(genestodoml, "Status")
commontrain = colnames(train) %in% genestodoml
train = train[,..commontrain]



# Select only genes
x_train <- train[,-1]

# Compute Z score
x_train <- as.data.frame(lapply(x_train, function(x) (x - mean(x))/sd(x) ))
x_train <- as.matrix(x_train )

# Code status
numeros = as.numeric(train$Status)
y_train <-replace(numeros, numeros == 2,0)




# Set seed
set.seed(123)

# Train ridge
cv_error <- cv.glmnet(
  x      = x_train,
  y      = factor(y_train),
  alpha  = 0,
  nfolds = 5, #
  type.measure = "mse",
  standardize  = F,
  family = "binomial",
  intercept = F
)



# Best model lambda
modelo <- glmnet(
  x           = x_train,
  y           = y_train,
  alpha       = 0,
  lambda      = cv_error$lambda.min,
  standardize = F,
  family = "binomial",
  intercept = F
)


# Training predictions
# ==============================================================================
predicciones_train <- predict(modelo, newx = x_train, type = "response", standardize = F)


# MSE (Training)
# ==============================================================================
training_mse <- mean((predicciones_train - y_train)^2)



caret::confusionMatrix(data = factor(as.numeric(predicciones_train>0.5)), reference = factor(y_train)) 



#Testing:  
  
  
  
  
  
  # Load Matrix
  test = readRDS("rlogmatrix_allpsold38qcp_notion")


# Select features
commontest = colnames(test) %in% genestodoml
test = test[,commontest]


# Select only genes
x_test <- test[,-1]

# Compute Z score
x_test <- as.data.frame(lapply(x_test, function(x) (x - mean(x))/sd(x) ))
x_test <- as.matrix(x_test )

# Code status
numeros = as.numeric(test$Status)
y_test <-replace(numeros, numeros == 1,0)
y_test <-replace(y_test, y_test == 2,1)
# y_test <-replace(numeros, numeros == 1,0)
# y_test <-replace(y_test, y_test == 2,1)



# Testing preditions
# ==============================================================================
predicciones_test <- predict(modelo, newx = x_test, type = "response" , standarize = F, intercept = F)


# MSE (testing)
# ==============================================================================
test_mse_lasso <- mean((predicciones_test - y_test)^2)




# Confusion Matrix
caret::confusionMatrix(data = factor(as.numeric(predicciones_test>0.5)), reference = factor(y_test))



df_coeficientes_final <- coef(modelo) %>%
  as.matrix() %>%
  as_tibble(rownames = "predictor")


df_coe_heatmap = subset(df_coeficientes_final, df_coeficientes_final$s0 != 0)

df_coe_heatmap_genes = df_coe_heatmap$predictor




#This is the model with 40 genes:




kltr =  0.08

difnaboold = kl_naboo_old[,kl_naboo_old < kltr]


genestodoml = colnames(difnaboold)



# Load Matrix
matrix_phase = readRDS("matrix_rlogk99qcADPreclinicalandControls_sex_agep_notion")
train = matrix_phase


# Select features
genestodoml = c(genestodoml, "Status")
commontrain = colnames(train) %in% genestodoml
train = train[,..commontrain]



# Select only genes
x_train <- train[,-1]

# Compute Z score
x_train <- as.data.frame(lapply(x_train, function(x) (x - mean(x))/sd(x) ))
x_train <- as.matrix(x_train )

# Code status
numeros = as.numeric(train$Status)
y_train <-replace(numeros, numeros == 2,0)




# Set seed
set.seed(123)

# Train ridge
cv_error <- cv.glmnet(
  x      = x_train,
  y      = factor(y_train),
  alpha  = 0,
  nfolds = 5, #
  type.measure = "mse",
  standardize  = F,
  family = "binomial",
  intercept = F
)



# Best model lambda
modelo <- glmnet(
  x           = x_train,
  y           = y_train,
  alpha       = 0,
  lambda      = cv_error$lambda.min,
  standardize = F,
  family = "binomial",
  intercept = F
)

df_coeficientes <- coef(modelo) %>%
  as.matrix() %>%
  as_tibble(rownames = "predictor")


# Training predictions
# ==============================================================================
predicciones_train <- predict(modelo, newx = x_train, type = "response", standardize = F)


# MSE (Training)
# ==============================================================================
training_mse <- mean((predicciones_train - y_train)^2)






df_coeficientestraning_rlog = subset(df_coeficientes, df_coeficientes$s0 != 0)


df_coeficientestraning_rlog$s0 = abs(df_coeficientestraning_rlog$s0)



df_coeficientestraning_rlog <- df_coeficientestraning_rlog[ order(df_coeficientestraning_rlog$s0, decreasing = T),]




df_coeficientestraning_rlogp = df_coeficientestraning_rlog[1:40,]
genestodoml = df_coeficientestraning_rlogp$predictor



# Load Matrix
matrix_phase = readRDS("matrix_rlogk99qcADPreclinicalandControls_sex_agep_notion")
train = matrix_phase


# Select features
genestodoml = c(genestodoml, "Status")
commontrain = colnames(train) %in% genestodoml
train = train[,..commontrain]



# Select only genes
x_train <- train[,-1]

# Compute Z score
x_train <- as.data.frame(lapply(x_train, function(x) (x - mean(x))/sd(x) ))
x_train <- as.matrix(x_train )

# Code status
numeros = as.numeric(train$Status)
y_train <-replace(numeros, numeros == 2,0)




# Set seed
set.seed(123)

# Train ridge
cv_error <- cv.glmnet(
  x      = x_train,
  y      = factor(y_train),
  alpha  = 0,
  nfolds = 5, #
  type.measure = "mse",
  standardize  = F,
  family = "binomial",
  intercept = F
)



# Best model lambda
modelo <- glmnet(
  x           = x_train,
  y           = y_train,
  alpha       = 0,
  lambda      = cv_error$lambda.min,
  standardize = F,
  family = "binomial",
  intercept = F
)


# Training predictions
# ==============================================================================
predicciones_train <- predict(modelo, newx = x_train, type = "response", standardize = F)


# MSE (Training)
# ==============================================================================
training_mse <- mean((predicciones_train - y_train)^2)



# Load Matrix
test = readRDS("rlogmatrix_allpsold38qcp_notion")


# Select features
commontest = colnames(test) %in% genestodoml
test = test[,commontest]


# Select only genes
x_test <- test[,-1]

# Compute Z score
x_test <- as.data.frame(lapply(x_test, function(x) (x - mean(x))/sd(x) ))
x_test <- as.matrix(x_test )

# Code status
numeros = as.numeric(test$Status)
y_test <-replace(numeros, numeros == 1,0)
y_test <-replace(y_test, y_test == 2,1)
# y_test <-replace(numeros, numeros == 1,0)
# y_test <-replace(y_test, y_test == 2,1)



# Testing preditions
# ==============================================================================
predicciones_test <- predict(modelo, newx = x_test, type = "response" , standarize = F, intercept = F)


# MSE (testing)
# ==============================================================================
test_mse_lasso <- mean((predicciones_test - y_test)^2)




# Confusion Matrix
caret::confusionMatrix(data = factor(as.numeric(predicciones_test>0.5)), reference = factor(y_test))



df_coeficientes_final <- coef(modelo) %>%
  as.matrix() %>%
  as_tibble(rownames = "predictor")


df_coe_heatmap = subset(df_coeficientes_final, df_coeficientes_final$s0 != 0)

df_coe_heatmap_genes = df_coe_heatmap$predictor







# Models further analyses

## 40 genes model



kltr =  0.08

difnaboold = kl_naboo_old[,kl_naboo_old < kltr]


genestodoml = colnames(difnaboold)



# Load Matrix
matrix_phase = readRDS("matrix_rlogk99qcADPreclinicalandControls_sex_agep_notion")
train = matrix_phase


# Select features
genestodoml = c(genestodoml, "Status")
commontrain = colnames(train) %in% genestodoml
train = train[,..commontrain]



# Select only genes
x_train <- train[,-1]

# Compute Z score
x_train <- as.data.frame(lapply(x_train, function(x) (x - mean(x))/sd(x) ))
x_train <- as.matrix(x_train )

# Code status
numeros = as.numeric(train$Status)
y_train <-replace(numeros, numeros == 2,0)




# Set seed
set.seed(123)

# Train ridge
cv_error <- cv.glmnet(
  x      = x_train,
  y      = factor(y_train),
  alpha  = 0,
  nfolds = 5, #
  type.measure = "mse",
  standardize  = F,
  family = "binomial",
  intercept = F
)



# Best model lambda
modelo <- glmnet(
  x           = x_train,
  y           = y_train,
  alpha       = 0,
  lambda      = cv_error$lambda.min,
  standardize = F,
  family = "binomial",
  intercept = F
)

df_coeficientes <- coef(modelo) %>%
  as.matrix() %>%
  as_tibble(rownames = "predictor")


# Training predictions
# ==============================================================================
predicciones_train <- predict(modelo, newx = x_train, type = "response", standardize = F)


# MSE (Training)
# ==============================================================================
training_mse <- mean((predicciones_train - y_train)^2)





df_coeficientestraning_rlog = subset(df_coeficientes, df_coeficientes$s0 != 0)


df_coeficientestraning_rlog$s0 = abs(df_coeficientestraning_rlog$s0)



df_coeficientestraning_rlog <- df_coeficientestraning_rlog[ order(df_coeficientestraning_rlog$s0, decreasing = T),]




df_coeficientestraning_rlogp = df_coeficientestraning_rlog[1:40,]
genestodoml = df_coeficientestraning_rlogp$predictor

genestodoml40 = genestodoml

# Load Matrix
matrix_phase = readRDS("matrix_rlogk99qcADPreclinicalandControls_sex_agep_notion")
train = matrix_phase


# Select features
genestodoml = c(genestodoml, "Status")
commontrain = colnames(train) %in% genestodoml
train = train[,..commontrain]



# Select only genes
x_train <- train[,-1]

# Compute Z score
x_train <- as.data.frame(lapply(x_train, function(x) (x - mean(x))/sd(x) ))
x_train <- as.matrix(x_train )

# Code status
numeros = as.numeric(train$Status)
y_train <-replace(numeros, numeros == 2,0)




# Set seed
set.seed(123)

# Train ridge
cv_error <- cv.glmnet(
  x      = x_train,
  y      = factor(y_train),
  alpha  = 0,
  nfolds = 5, #
  type.measure = "mse",
  standardize  = F,
  family = "binomial",
  intercept = F
)



# Best model lambda
modelo <- glmnet(
  x           = x_train,
  y           = y_train,
  alpha       = 0,
  lambda      = cv_error$lambda.min,
  standardize = F,
  family = "binomial",
  intercept = F
)


# Training predictions
# ==============================================================================
predicciones_train <- predict(modelo, newx = x_train, type = "response", standardize = F)


# MSE (Training)
# ==============================================================================
training_mse <- mean((predicciones_train - y_train)^2)



caret::confusionMatrix(data = factor(as.numeric(predicciones_train>0.5)), reference = factor(y_train)) 





genestodoml = genestodoml40
genestodoml = c(genestodoml, "Status")
# Load Matrix
test = readRDS("rlogmatrix_allpsold38qcp_notion")


# Select features
commontest = colnames(test) %in% genestodoml
test = test[,commontest]


# Select only genes
x_test <- test[,-1]

# Compute Z score
x_test <- as.data.frame(lapply(x_test, function(x) (x - mean(x))/sd(x) ))
x_test <- as.matrix(x_test )

# Code status
numeros = as.numeric(test$Status)
y_test <-replace(numeros, numeros == 1,0)
y_test <-replace(y_test, y_test == 2,1)
# y_test <-replace(numeros, numeros == 1,0)
# y_test <-replace(y_test, y_test == 2,1)



# Testing preditions
# ==============================================================================
predicciones_test <- predict(modelo, newx = x_test, type = "response" , standarize = F, intercept = F)


# MSE (testing)
# ==============================================================================
test_mse_lasso <- mean((predicciones_test - y_test)^2)




# Confusion Matrix
caret::confusionMatrix(data = factor(as.numeric(predicciones_test>0.5)), reference = factor(y_test))


# 
# df_coeficientes_final <- coef(modelo) %>%
#   as.matrix() %>%
#   as_tibble(rownames = "predictor")
# 
# 
# df_coe_heatmap = subset(df_coeficientes_final, df_coeficientes_final$s0 != 0)

# df_coe_heatmap_genes = df_coe_heatmap$predictor







## 90 genes model






kltr =  0.31

difnaboold = kl_naboo_old[,kl_naboo_old < kltr]


genestodoml = colnames(difnaboold)



# Load Matrix
matrix_phase = readRDS("matrix_rlogk99qcADPreclinicalandControls_sex_agep_notion")
train = matrix_phase


# Select features
genestodoml = c(genestodoml, "Status")
commontrain = colnames(train) %in% genestodoml
train = train[,..commontrain]



# Select only genes
x_train <- train[,-1]

# Compute Z score
x_train <- as.data.frame(lapply(x_train, function(x) (x - mean(x))/sd(x) ))
x_train <- as.matrix(x_train )

# Code status
numeros = as.numeric(train$Status)
y_train <-replace(numeros, numeros == 2,0)




# Set seed
set.seed(123)

# Train ridge
cv_error <- cv.glmnet(
  x      = x_train,
  y      = factor(y_train),
  alpha  = 0,
  nfolds = 5, #
  type.measure = "mse",
  standardize  = F,
  family = "binomial",
  intercept = F
)



# Best model lambda
modelo <- glmnet(
  x           = x_train,
  y           = y_train,
  alpha       = 0,
  lambda      = cv_error$lambda.min,
  standardize = F,
  family = "binomial",
  intercept = F
)

df_coeficientes <- coef(modelo) %>%
  as.matrix() %>%
  as_tibble(rownames = "predictor")


# Training predictions
# ==============================================================================
predicciones_train <- predict(modelo, newx = x_train, type = "response", standardize = F)


# MSE (Training)
# ==============================================================================
training_mse <- mean((predicciones_train - y_train)^2)





df_coeficientestraning_rlog = subset(df_coeficientes, df_coeficientes$s0 != 0)


df_coeficientestraning_rlog$s0 = abs(df_coeficientestraning_rlog$s0)



df_coeficientestraning_rlog <- df_coeficientestraning_rlog[ order(df_coeficientestraning_rlog$s0, decreasing = T),]




df_coeficientestraning_rlogp = df_coeficientestraning_rlog[1:90,]
genestodoml = df_coeficientestraning_rlogp$predictor

genestodoml90 = genestodoml

# Load Matrix
matrix_phase = readRDS("matrix_rlogk99qcADPreclinicalandControls_sex_agep_notion")
train = matrix_phase


# Select features
genestodoml = c(genestodoml, "Status")
commontrain = colnames(train) %in% genestodoml
train = train[,..commontrain]



# Select only genes
x_train <- train[,-1]

# Compute Z score
x_train <- as.data.frame(lapply(x_train, function(x) (x - mean(x))/sd(x) ))
x_train <- as.matrix(x_train )

# Code status
numeros = as.numeric(train$Status)
y_train <-replace(numeros, numeros == 2,0)




# Set seed
set.seed(123)

# Train ridge
cv_error <- cv.glmnet(
  x      = x_train,
  y      = factor(y_train),
  alpha  = 0,
  nfolds = 5, #
  type.measure = "mse",
  standardize  = F,
  family = "binomial",
  intercept = F
)



# Best model lambda
modelo <- glmnet(
  x           = x_train,
  y           = y_train,
  alpha       = 0,
  lambda      = cv_error$lambda.min,
  standardize = F,
  family = "binomial",
  intercept = F
)


# Training predictions
# ==============================================================================
predicciones_train <- predict(modelo, newx = x_train, type = "response", standardize = F)


# MSE (Training)
# ==============================================================================
training_mse <- mean((predicciones_train - y_train)^2)



caret::confusionMatrix(data = factor(as.numeric(predicciones_train>0.5)), reference = factor(y_train)) 



#Testing:  




genestodoml = genestodoml90
genestodoml = c(genestodoml, "Status")
# Load Matrix
test = readRDS("rlogmatrix_allpsold38qcp_notion")


# Select features
commontest = colnames(test) %in% genestodoml
test = test[,commontest]


# Select only genes
x_test <- test[,-1]

# Compute Z score
x_test <- as.data.frame(lapply(x_test, function(x) (x - mean(x))/sd(x) ))
x_test <- as.matrix(x_test )

# Code status
numeros = as.numeric(test$Status)
y_test <-replace(numeros, numeros == 1,0)
y_test <-replace(y_test, y_test == 2,1)
# y_test <-replace(numeros, numeros == 1,0)
# y_test <-replace(y_test, y_test == 2,1)



# Testing preditions
# ==============================================================================
predicciones_test <- predict(modelo, newx = x_test, type = "response" , standarize = F, intercept = F)


# MSE (testing)
# ==============================================================================
test_mse_lasso <- mean((predicciones_test - y_test)^2)




# Confusion Matrix
caret::confusionMatrix(data = factor(as.numeric(predicciones_test>0.5)), reference = factor(y_test))


# 
# df_coeficientes_final <- coef(modelo) %>%
#   as.matrix() %>%
#   as_tibble(rownames = "predictor")
# 
# 
# df_coe_heatmap = subset(df_coeficientes_final, df_coeficientes_final$s0 != 0)

# df_coe_heatmap_genes = df_coe_heatmap$predictor








## 220 genes model


#Training:  


kltr =  0.49

difnaboold = kl_naboo_old[,kl_naboo_old < kltr]


genestodoml = colnames(difnaboold)



# Load Matrix
matrix_phase = readRDS("matrix_rlogk99qcADPreclinicalandControls_sex_agep_notion")
train = matrix_phase


# Select features
genestodoml = c(genestodoml, "Status")
commontrain = colnames(train) %in% genestodoml
train = train[,..commontrain]



# Select only genes
x_train <- train[,-1]

# Compute Z score
x_train <- as.data.frame(lapply(x_train, function(x) (x - mean(x))/sd(x) ))
x_train <- as.matrix(x_train )

# Code status
numeros = as.numeric(train$Status)
y_train <-replace(numeros, numeros == 2,0)




# Set seed
set.seed(123)

# Train ridge
cv_error <- cv.glmnet(
  x      = x_train,
  y      = factor(y_train),
  alpha  = 0,
  nfolds = 5, #
  type.measure = "mse",
  standardize  = F,
  family = "binomial",
  intercept = F
)



# Best model lambda
modelo <- glmnet(
  x           = x_train,
  y           = y_train,
  alpha       = 0,
  lambda      = cv_error$lambda.min,
  standardize = F,
  family = "binomial",
  intercept = F
)

df_coeficientes <- coef(modelo) %>%
  as.matrix() %>%
  as_tibble(rownames = "predictor")


# Training predictions
# ==============================================================================
predicciones_train <- predict(modelo, newx = x_train, type = "response", standardize = F)


# MSE (Training)
# ==============================================================================
training_mse <- mean((predicciones_train - y_train)^2)





df_coeficientestraning_rlog = subset(df_coeficientes, df_coeficientes$s0 != 0)


df_coeficientestraning_rlog$s0 = abs(df_coeficientestraning_rlog$s0)



df_coeficientestraning_rlog <- df_coeficientestraning_rlog[ order(df_coeficientestraning_rlog$s0, decreasing = T),]




df_coeficientestraning_rlogp = df_coeficientestraning_rlog[1:220,]
genestodoml = df_coeficientestraning_rlogp$predictor

genestodoml220 = genestodoml

# Load Matrix
matrix_phase = readRDS("matrix_rlogk99qcADPreclinicalandControls_sex_agep_notion")
train = matrix_phase


# Select features
genestodoml = c(genestodoml, "Status")
commontrain = colnames(train) %in% genestodoml
train = train[,..commontrain]



# Select only genes
x_train <- train[,-1]

# Compute Z score
x_train <- as.data.frame(lapply(x_train, function(x) (x - mean(x))/sd(x) ))
x_train <- as.matrix(x_train )

# Code status
numeros = as.numeric(train$Status)
y_train <-replace(numeros, numeros == 2,0)




# Set seed
set.seed(123)

# Train ridge
cv_error <- cv.glmnet(
  x      = x_train,
  y      = factor(y_train),
  alpha  = 0,
  nfolds = 5, #
  type.measure = "mse",
  standardize  = F,
  family = "binomial",
  intercept = F
)



# Best model lambda
modelo <- glmnet(
  x           = x_train,
  y           = y_train,
  alpha       = 0,
  lambda      = cv_error$lambda.min,
  standardize = F,
  family = "binomial",
  intercept = F
)


# Training predictions
# ==============================================================================
predicciones_train <- predict(modelo, newx = x_train, type = "response", standardize = F)


# MSE (Training)
# ==============================================================================
training_mse <- mean((predicciones_train - y_train)^2)



caret::confusionMatrix(data = factor(as.numeric(predicciones_train>0.5)), reference = factor(y_train)) 




#Testing:  



genestodoml = genestodoml220
genestodoml = c(genestodoml, "Status")
# Load Matrix
test = readRDS("rlogmatrix_allpsold38qcp_notion")


# Select features
commontest = colnames(test) %in% genestodoml
test = test[,commontest]


# Select only genes
x_test <- test[,-1]

# Compute Z score
x_test <- as.data.frame(lapply(x_test, function(x) (x - mean(x))/sd(x) ))
x_test <- as.matrix(x_test )

# Code status
numeros = as.numeric(test$Status)
y_test <-replace(numeros, numeros == 1,0)
y_test <-replace(y_test, y_test == 2,1)
# y_test <-replace(numeros, numeros == 1,0)
# y_test <-replace(y_test, y_test == 2,1)



# Testing preditions
# ==============================================================================
predicciones_test <- predict(modelo, newx = x_test, type = "response" , standarize = F, intercept = F)


# MSE (testing)
# ==============================================================================
test_mse_lasso <- mean((predicciones_test - y_test)^2)




# Confusion Matrix
caret::confusionMatrix(data = factor(as.numeric(predicciones_test>0.5)), reference = factor(y_test))



plot(predicciones_test, y_test, xlab = "Predicted", ylab = "Actual",
     main = "Probabilities of each class")
misClassError(y_test, predicciones_test)
plotROC(y_test, predicciones_test)


# 
# df_coeficientes_final <- coef(modelo) %>%
#   as.matrix() %>%
#   as_tibble(rownames = "predictor")
# 
# 
# df_coe_heatmap = subset(df_coeficientes_final, df_coeficientes_final$s0 != 0)

# df_coe_heatmap_genes = df_coe_heatmap$predictor







