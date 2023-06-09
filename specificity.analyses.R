#Load the necessary libraries
library(data.table)
library(dplyr)
library(readr)
library(ggplot2)
library(scales)
library(glmnet)
library(pROC)
myPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999", "#000000", "#6C2EB8")
#Define a function that scales two data.frames
i = 1
scaletwosets = function(trainingreference, dftoscale){
  # As data frame
  trainingreference = as.data.frame(trainingreference)
  dftoscale = as.data.frame(dftoscale)
  columnscommon = which(colnames(trainingreference) %in% colnames(dftoscale))
  trainingreference = trainingreference[,columnscommon]
  # Make data frame with the same number of rows and columns
  dfscalated = data.frame(matrix(0, nrow = nrow(dftoscale), ncol = ncol(dftoscale),
                                 dimnames = list(NULL, colnames(dftoscale))) )
  for(i in 1:length(colnames(trainingreference))){
    # Compute max and min for each feature in the training reference
    df1range = NULL
    df1range = range(trainingreference[,i])
    df1rangemin = df1range[1]
    df1rangemax = df1range[2]
    # Scale to this min and max the values in the dftoscale
    vec_range <- rescale(dftoscale[,i], to = c(df1rangemin, df1rangemax))
    dfscalated[,i] = vec_range
  }
  return(dfscalated)
}
#####################################################################
### Compare AD and other diseases (PD, DLB and FTD) 
### with and without APOE 
#####################################################################
##### 220-genes model #####
genestodoml220 = fread("220genes.csv")
genestodoml220 = genestodoml220$ENSG
genestodoml = genestodoml220
genestodoml = c(genestodoml, "Status")
# Load Matrix
matrix_phase = readRDS("discovery_data")
train = matrix_phase
APOEtrain<-train$apoe
trainconpred = train
# Select features
commontrain = colnames(train) %in% genestodoml
train<-as.data.frame(train)[,colnames(train) %in% genestodoml]
# Select only genes
x_train <- train[,-1]
# Compute Z score
modelo_means <- as.data.frame(lapply(x_train, mean))
modelo_sd <- as.data.frame(lapply(x_train, sd))
x_train <- as.data.frame(lapply(x_train, function(x) (x - mean(x))/sd(x) ))
x_train <- as.matrix(x_train )
# Code status
numeros = as.numeric(train$Status)
y_train <-replace(numeros, numeros == 2,0)
y_train_220 = y_train 
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
predicciones_train_220 <- predict(modelo, newx = x_train, type = "response", standardize = F)
# Load Matrix
test = readRDS("replication_data")
APOEtest<-test$APOE
tocorrelationtest = test[,c(1,2,3,4,6)]
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
y_test_220 = y_test
# Testing preditions
# ==============================================================================
predicciones_test <- predict(modelo, newx = x_test, type = "response" , standarize = F, intercept = F)
predicciones_test_220 = predicciones_test
# Model220 roc comes from above code (220 genes model)
model220rocte = cbind(as.data.frame(y_test_220),APOEtest)
model220rocte$predictor = predicciones_test_220
names(model220rocte) = c("status", "APOE","model220")
model220roctr = cbind(as.data.frame(y_train_220),APOEtrain)
model220roctr$predictor = predicciones_train_220
names(model220roctr) = c("status", "APOE","model220")
model220ad = rbind(model220rocte,model220roctr)
model220ad <- model220ad[model220ad$status == 0,]

### Specificity PD individuals
# Load Matrix
matrix_phase = readRDS("discovery_data")
train = matrix_phase
train<-as.data.frame(train)[,colnames(train) %in% genestodoml]
# Select only genes
x_train <- train[,-1]
#Load testing data
test = readRDS("specificity_data_PD")
test$SampleID<-sub("MARS_", "MARS", test$SampleID)
PDapoe<-read.csv("PD_MovementCenter_APOE.csv")
test<-merge(test, PDapoe,by="SampleID")
APOEtest <- test$APOE
test<-as.data.frame(test)[,colnames(test) %in% genestodoml]
x_test <- test[,-1]
# Select only genes
pdscalated = scaletwosets(x_train, x_test)
PDmatrix_means <- as.data.frame(lapply(pdscalated, mean))
x_test <- pdscalated
#dfmodelomenas = as.data.frame(modelo_means[col(x_test[])])
x_testprueba = x_test[]-modelo_means[col(x_test[])]
x_testpruebafinal = x_testprueba[]/modelo_sd[col(x_testprueba[])]
x_test = x_testpruebafinal
x_test = as.data.frame(x_test)
# add gene name from cdr05
#x_test = cbind(x_test, setNames( lapply(genesnoaparecen, function(x) x=0), genesnoaparecen) )
x_test <- as.matrix(x_test)
# Code status
numeros = as.numeric(test$Status)
y_test <-replace(numeros, numeros == 1,0)
y_test <-replace(y_test, y_test == 2,1)
y_test_pd = y_test
pdtoeatmap = as.data.frame(x_test)
# Testing preditions
# ==============================================================================
predicciones_test <- predict(modelo, newx = x_test, type = "response" , standarize = F, intercept = F) 
predicciones_test_pd = predicciones_test
model220pd = cbind(as.data.frame(y_test_pd),APOEtest)
model220pd$predictor = predicciones_test_pd
names(model220pd) = c("status", "APOE","model220")
model220pd$status<-1
model220pd<-rbind(model220pd,model220ad)

###DLB Specificity
test = readRDS("specificity_data")
genestodoml = c(genestodoml, "Status")
test = subset(test, test$Status == "DLB")
APOEtest <- test$apoe
commontest = colnames(test) %in% genestodoml
test = test[,..commontest]
x_test <- test[,-1]
dlbscalated = scaletwosets(x_train, x_test)
genesgenestodoml = as.data.frame(genestodoml)
colnames(genesgenestodoml) = "gene"
nombrestest = as.data.frame(colnames(test))
colnames(nombrestest) = "gene"
genesnoaparecen = anti_join(genesgenestodoml,nombrestest)
genesnoaparecen = as.character(genesnoaparecen$gene)
# Select only genes and the scalated df
x_test <- dlbscalated
modelo_meansftdp = modelo_means
which(colnames(modelo_meansftdp) %in% genesnoaparecen)  
modelo_meansftd = modelo_meansftdp[,-c(170,178)]
modelo_sdftdp = modelo_sd
modelo_sdftdpd = modelo_sdftdp[,-c(170,178)]
x_testprueba = x_test[]-modelo_meansftd[col(x_test[])]
x_testpruebafinal = x_testprueba[]/modelo_sdftdpd[col(x_testprueba[])]
x_test = x_testpruebafinal
x_test = as.data.frame(x_test)
# add gene name from cdr05
x_test = cbind(x_test, setNames( lapply(genesnoaparecen, function(x) x=0), genesnoaparecen))
x_test <- as.matrix(x_test )
dlbtoeatmap = as.data.frame(x_test)
# Code status
numeros = as.numeric(test$Status)
y_test <-replace(numeros, numeros == 2,1)
y_test_dlb = y_test
# Testing preditions
# ==============================================================================
predicciones_test <- predict(modelo, newx = x_test, type = "response" , standarize = F, intercept = F) 
predicciones_test_dlb = predicciones_test
model220dlb = cbind(as.data.frame(y_test_dlb), APOEtest)
model220dlb$predictor = predicciones_test_dlb
names(model220dlb) = c("status", "APOE","model220")
model220dlb$status<-1
model220dlb <-rbind(model220dlb, model220ad)

###FTD Specificity
test = readRDS("specificity_data")
genestodoml = c(genestodoml, "Status")
test = subset(test, test$Status == "FTD")
APOEtest <- test$apoe
commontest = colnames(test) %in% genestodoml
test = test[,..commontest]
x_test <- test[,-1]
ftdscalated = scaletwosets(x_train, x_test)
nombrestest = as.data.frame(colnames(test))
colnames(nombrestest) = "gene"
genesnoaparecen = anti_join(genesgenestodoml,nombrestest)
genesnoaparecen = as.character(genesnoaparecen$gene)
# Select only genes and the scalated df
x_test <- ftdscalated
modelo_meansftdp = modelo_means
modelo_meansftd = modelo_meansftdp[,-c(170,178)]
modelo_sdftdp = modelo_sd
modelo_sdftdpd = modelo_sdftdp[,-c(170,178)]
x_testprueba = x_test[]-modelo_meansftd[col(x_test[])]
x_testpruebafinal = x_testprueba[]/modelo_sdftdpd[col(x_testprueba[])]
x_test = x_testpruebafinal
x_test = as.data.frame(x_test)
# add gene name from cdr05
x_test = cbind(x_test, setNames( lapply(genesnoaparecen, function(x) x=0), genesnoaparecen))
x_test <- as.matrix(x_test )
ftdtoeatmap = as.data.frame(x_test)
# Code status
numeros = as.numeric(test$Status)
y_test <-replace(numeros, numeros == 2,1)
y_test_ftd = y_test
# Testing preditions
# ==============================================================================
predicciones_test <- predict(modelo, newx = x_test, type = "response" , standarize = F, intercept = F) 
predicciones_test_ftd = predicciones_test
model220ftd = cbind(as.data.frame(y_test_ftd),APOEtest)
model220ftd$predictor = predicciones_test_ftd
names(model220ftd) = c("status", "APOE","model220")
model220ftd$status<-1
model220ftd<-rbind(model220ftd, model220ad)

###90 genes###
genestodoml90 = fread("90genes.csv")
genestodoml90 = genestodoml90[-1,1]
genestodoml90 = genestodoml90$V1
genestodoml = genestodoml90
genestodoml = c(genestodoml, "Status")
# Load Matrix
matrix_phase = readRDS("discovery_data")
train = matrix_phase
trainconpred = train
# Select features
commontrain = colnames(train) %in% genestodoml
train<-as.data.frame(train)[,colnames(train) %in% genestodoml]
# Select only genes
x_train <- train[,-1]
# Compute Z score
modelo_means <- as.data.frame(lapply(x_train, mean))
modelo_sd <- as.data.frame(lapply(x_train, sd))
x_train <- as.data.frame(lapply(x_train, function(x) (x - mean(x))/sd(x) ))
x_train <- as.matrix(x_train )
# Code status
numeros = as.numeric(train$Status)
y_train <-replace(numeros, numeros == 2,0)
y_train_90 = y_train 
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
predicciones_train_90 <- predict(modelo, newx = x_train, type = "response", standardize = F)
# Load Matrix
test = readRDS("replication_data")
tocorrelationtest = test[,c(1,2,3,4,6)]
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
y_test_90 = y_test
# Testing preditions
# ==============================================================================
predicciones_test <- predict(modelo, newx = x_test, type = "response" , standarize = F, intercept = F)
predicciones_test_90 = predicciones_test
# Model90 roc comes from above code (90 genes model)
model90rocte = as.data.frame(y_test_90)
model90rocte$predictor = predicciones_test_90
names(model90rocte) = c("status", "model90")
model90roctr = as.data.frame(y_train_90)
model90roctr$predictor = predicciones_train_90
names(model90roctr) = c("status", "model90")
model90ad = rbind(model90rocte,model90roctr)
model90ad<-model90ad[model90ad$status == 0,]

### Specificity PD individuals
# Load Matrix
matrix_phase = readRDS("discovery_data")
train = matrix_phase
train<-as.data.frame(train)[,colnames(train) %in% genestodoml]
# Select only genes
x_train <- train[,-1]
#Load testing data
test = readRDS("specificity_data_PD")
test$SampleID<-sub("MARS_", "MARS", test$SampleID)
PDapoe<-read.csv("PD_MovementCenter_APOE.csv")
test<-merge(test, PDapoe,by="SampleID")
APOEtest <- test$APOE
test<-as.data.frame(test)[,colnames(test) %in% genestodoml]
# Select only genes
x_test <- test[,-1]
# Esto es una prueba
pdscalated = scaletwosets(x_train, x_test)
PDmatrix_means <- as.data.frame(lapply(pdscalated, mean))
x_test <- pdscalated
#dfmodelomenas = as.data.frame(modelo_means[col(x_test[])])
x_testprueba = x_test[]-modelo_means[col(x_test[])]
x_testpruebafinal = x_testprueba[]/modelo_sd[col(x_testprueba[])]
x_test = x_testpruebafinal
x_test = as.data.frame(x_test)
# add gene name from cdr05
#x_test = cbind(x_test, setNames( lapply(genesnoaparecen, function(x) x=0), genesnoaparecen) )
x_test <- as.matrix(x_test )
# Code status
numeros = as.numeric(test$Status)
y_test <-replace(numeros, numeros == 1,0)
y_test <-replace(y_test, y_test == 2,1)
y_test_pd = y_test
pdtoeatmap = as.data.frame(x_test)
# Testing preditions
# ==============================================================================
predicciones_test <- predict(modelo, newx = x_test, type = "response" , standarize = F, intercept = F) 
predicciones_test_pd = predicciones_test
model90pd = as.data.frame(y_test_pd)
model90pd$predictor = predicciones_test_pd
names(model90pd) = c("status", "model90")
model90pd$status<-1
model90pd<-rbind(model90pd, model90ad)

### Specificity DLB individuals
test = readRDS("specificity_data")
test = subset(test, test$Status == "DLB")
commontest = colnames(test) %in% genestodoml
test = test[,..commontest]
# Select only genes
x_test <- test[,-1]
dlbscalated = scaletwosets(x_train, x_test)
x_test <- dlbscalated
genesgenestodoml = as.data.frame(genestodoml)
colnames(genesgenestodoml) = "gene"
nombrestest = as.data.frame(colnames(test))
colnames(nombrestest) = "gene"
genesnoaparecen = anti_join(genesgenestodoml,nombrestest)
genesnoaparecen = as.character(genesnoaparecen$gene)
modelo_meansftdp = modelo_means
which(colnames(modelo_meansftdp) %in% genesnoaparecen)  
modelo_meansftd = modelo_meansftdp[,-66]
modelo_sdftdp = modelo_sd
modelo_sdftdpd = modelo_sdftdp[,-66]
x_testprueba = x_test[]-modelo_meansftd[col(x_test[])]
x_testpruebafinal = x_testprueba[]/modelo_sdftdpd[col(x_testprueba[])]
x_test = x_testpruebafinal
x_test = as.data.frame(x_test)
# add gene name from cdr05
x_test = cbind(x_test, setNames( lapply(genesnoaparecen, function(x) x=0), genesnoaparecen))
x_test <- as.matrix(x_test )
dlbtoeatmap = as.data.frame(x_test)
# Code status
numeros = as.numeric(test$Status)
y_test <-replace(numeros, numeros == 2,1)
y_test_dlb = y_test
# Testing preditions
# ==============================================================================
predicciones_test <- predict(modelo, newx = x_test, type = "response" , standarize = F, intercept = F) 
predicciones_test_dlb = predicciones_test
model90dlb = as.data.frame(y_test_dlb)
model90dlb$predictor = predicciones_test_dlb
names(model90dlb) = c("status", "model90")
model90dlb$status<-1
model90dlb<-rbind(model90dlb, model90ad)

### Specificity FTD individuals
test = readRDS("specificity_data")
test = subset(test, test$Status == "FTD")
test = as.data.frame(test)[,colnames(test) %in% genestodoml]
# Select only genes
x_test <- test[,-1]
ftdscalated = scaletwosets(x_train, x_test)
x_test <- ftdscalated
nombrestest = as.data.frame(colnames(test))
colnames(nombrestest) = "gene"
genesnoaparecen = anti_join(genesgenestodoml,nombrestest)
genesnoaparecen = as.character(genesnoaparecen$gene)
modelo_meansftdp = modelo_means
which(colnames(modelo_meansftdp) %in% genesnoaparecen)  
modelo_meansftd = modelo_meansftdp[,-66]
modelo_sdftdp = modelo_sd
modelo_sdftdpd = modelo_sdftdp[,-66]
x_testprueba = x_test[]-modelo_meansftd[col(x_test[])]
x_testpruebafinal = x_testprueba[]/modelo_sdftdpd[col(x_testprueba[])]
x_test = x_testpruebafinal
x_test = as.data.frame(x_test)
# add gene name from cdr05
x_test = cbind(x_test, setNames( lapply(genesnoaparecen, function(x) x=0), genesnoaparecen))
x_test <- as.matrix(x_test )
ftdtoeatmap = as.data.frame(x_test)
# Code status
numeros = as.numeric(test$Status)
y_test <-replace(numeros, numeros == 2,1)
y_test_ftd = y_test
# Testing preditions
# ==============================================================================
predicciones_test <- predict(modelo, newx = x_test, type = "response" , standarize = F, intercept = F) 
predicciones_test_ftd = predicciones_test
model90ftd = as.data.frame(y_test_ftd)
model90ftd$predictor = predicciones_test_ftd
names(model90ftd) = c("status", "model90")
model90ftd$status<-1
model90ftd<-rbind(model90ftd,model90ad)


###40-genes model####
genestodoml40 = fread("40genes.csv")
genestodoml40 = genestodoml40[-1,1]
genestodoml40 = genestodoml40$V1
genestodoml = genestodoml40
genestodoml = c(genestodoml, "Status")
# Load Matrix
matrix_phase = readRDS("discovery_data")
train = matrix_phase
trainconpred = train
# Select features
commontrain = colnames(train) %in% genestodoml
train<-as.data.frame(train)[,colnames(train) %in% genestodoml]
# Select only genes
x_train <- train[,-1]
# Compute Z score
modelo_means <- as.data.frame(lapply(x_train, mean))
modelo_sd <- as.data.frame(lapply(x_train, sd))
x_train <- as.data.frame(lapply(x_train, function(x) (x - mean(x))/sd(x) ))
x_train <- as.matrix(x_train )
# Code status
numeros = as.numeric(train$Status)
y_train <-replace(numeros, numeros == 2,0)
y_train_40 = y_train 
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
predicciones_train_40 <- predict(modelo, newx = x_train, type = "response", standardize = F)
# Load Matrix
test = readRDS("replication_data")
tocorrelationtest = test[,c(1,2,3,4,6)]
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
y_test_40 = y_test
# Testing preditions
# ==============================================================================
predicciones_test <- predict(modelo, newx = x_test, type = "response" , standarize = F, intercept = F)
predicciones_test_40 = predicciones_test
# Model40 roc comes from above code (40 genes model)
model40rocte = as.data.frame(y_test_40)
model40rocte$predictor = predicciones_test_40
names(model40rocte) = c("status", "model40")
model40roctr = as.data.frame(y_train_40)
model40roctr$predictor = predicciones_train_40
names(model40roctr) = c("status", "model40")
model40ad = rbind(model40rocte,model40roctr)
model40ad <- model40ad[model40ad$status == 0,]

### Specificity PD individuals
# Load Matrix
matrix_phase = readRDS("discovery_data")
train = matrix_phase
train<-as.data.frame(train)[,colnames(train) %in% genestodoml]
# Select only genes
x_train <- train[,-1]
#Load testing data
test = readRDS("specificity_data_PD")
test$SampleID<-sub("MARS_", "MARS", test$SampleID)
PDapoe<-read.csv("PD_MovementCenter_APOE.csv")
test<-merge(test, PDapoe,by="SampleID")
APOEtest <- test$APOE
test<-as.data.frame(test)[,colnames(test) %in% genestodoml]
x_test <- test[,-1]
pdscalated = scaletwosets(x_train, x_test)
PDmatrix_means <- as.data.frame(lapply(pdscalated, mean))
# Select only genes
x_test <- pdscalated
#dfmodelomenas = as.data.frame(modelo_means[col(x_test[])])
x_testprueba = x_test[]-modelo_means[col(x_test[])]
x_testpruebafinal = x_testprueba[]/modelo_sd[col(x_testprueba[])]
x_test = x_testpruebafinal
x_test = as.data.frame(x_test)
x_test <- as.matrix(x_test )
# Code status
numeros = as.numeric(test$Status)
y_test <-replace(numeros, numeros == 1,0)
y_test <-replace(y_test, y_test == 2,1)
y_test_pd = y_test
pdtoeatmap = as.data.frame(x_test)
# Testing preditions
# ==============================================================================
predicciones_test <- predict(modelo, newx = x_test, type = "response" , standarize = F, intercept = F) 
predicciones_test_pd = predicciones_test
model40pd = as.data.frame(y_test_pd)
model40pd$predictor = predicciones_test_pd
names(model40pd) = c("status", "model40")
model40pd$status<-1
model40pd<-rbind(model40pd, model40ad)

### Specificity DLB individuals
test = readRDS("specificity_data")
test = subset(test, test$Status == "DLB")
commontest = colnames(test) %in% genestodoml
test = test[,..commontest]
# Select only genes
x_test <- test[,-1]
dlbscalated = scaletwosets(x_train, x_test)
x_test <- dlbscalated
genesgenestodoml = as.data.frame(genestodoml)
colnames(genesgenestodoml) = "gene"
nombrestest = as.data.frame(colnames(test))
colnames(nombrestest) = "gene"
genesnoaparecen = anti_join(genesgenestodoml,nombrestest)
genesnoaparecen = as.character(genesnoaparecen$gene)
modelo_meansftdp = modelo_means
which(colnames(modelo_meansftdp) %in% genesnoaparecen)  
modelo_meansftd = modelo_meansftdp[,-c(170,178)]
modelo_sdftdp = modelo_sd
modelo_sdftdpd = modelo_sdftdp[,-c(170,178)]
x_testprueba = x_test[]-modelo_meansftd[col(x_test[])]
x_testpruebafinal = x_testprueba[]/modelo_sdftdpd[col(x_testprueba[])]
x_test = x_testpruebafinal
x_test = as.data.frame(x_test)
# add gene name from cdr05
# x_test = cbind(x_test, setNames( lapply(genesnoaparecen, function(x) x=0), genesnoaparecen))
x_test <- as.matrix(x_test )
dlbtoeatmap = as.data.frame(x_test)
# Code status
numeros = as.numeric(test$Status)
y_test <-replace(numeros, numeros == 2,1)
y_test_dlb = y_test
# Testing preditions
# ==============================================================================
predicciones_test <- predict(modelo, newx = x_test, type = "response" , standarize = F, intercept = F) 
predicciones_test_dlb = predicciones_test
model40dlb = as.data.frame(y_test_dlb)
model40dlb$predictor = predicciones_test_dlb
names(model40dlb) = c("status", "model40")
model40dlb$status<-1
model40dlb<-rbind(model40dlb,model40ad)

### Specificity FTD individuals
test = readRDS("specificity_data")
test = subset(test, test$Status == "FTD")
commontest = colnames(test) %in% genestodoml
test = test[,..commontest]
# Select only genes
x_test <- test[,-1]
ftdscalated = scaletwosets(x_train, x_test)
x_test <- ftdscalated
nombrestest = as.data.frame(colnames(test))
colnames(nombrestest) = "gene"
genesnoaparecen = anti_join(genesgenestodoml,nombrestest)
genesnoaparecen = as.character(genesnoaparecen$gene)
modelo_meansftdp = modelo_means
modelo_meansftd = modelo_meansftdp[,-c(170,178)]
modelo_sdftdp = modelo_sd
modelo_sdftdpd = modelo_sdftdp[,-c(170,178)]
x_testprueba = x_test[]-modelo_meansftd[col(x_test[])]
x_testpruebafinal = x_testprueba[]/modelo_sdftdpd[col(x_testprueba[])]
x_test = x_testpruebafinal
x_test = as.data.frame(x_test)
# add gene name from cdr05
x_test <- as.matrix(x_test )
ftdtoeatmap = as.data.frame(x_test)
# Code status
numeros = as.numeric(test$Status)
y_test <-replace(numeros, numeros == 2,1)
y_test_ftd = y_test
# Testing preditions
# ==============================================================================
predicciones_test <- predict(modelo, newx = x_test, type = "response" , standarize = F, intercept = F) 
predicciones_test_ftd = predicciones_test
model40ftd = as.data.frame(y_test_ftd)
model40ftd$predictor = predicciones_test_ftd
names(model40ftd) = c("status", "model40")
model40ftd$status<-1
model40ftd<-rbind(model40ftd, model40ad)


# Combine all data for PD, DLB, and FTD into 
# a single data.frame per disease
modelpd<-cbind(model220pd, as.data.frame(model90pd$model90))
modelpd<-cbind(modelpd, as.data.frame(model40pd$model40))
names(modelpd)<-c("Status", "APOE", "model220genes", "model90genes", "model40genes")
modeldlb<-cbind(model220dlb, as.data.frame(model90dlb$model90))
modeldlb<-cbind(modeldlb, as.data.frame(model40dlb$model40))
names(modeldlb)<-c("Status", "APOE", "model220genes", "model90genes", "model40genes")
modelftd<-cbind(model220ftd, as.data.frame(model90ftd$model90))
modelftd<-cbind(modelftd, as.data.frame(model40ftd$model40))
names(modelftd)<-c("Status", "APOE", "model220genes", "model90genes", "model40genes")
modelpd$model220genes<-as.numeric(modelpd$model220genes)
modeldlb$model220genes<-as.numeric(modeldlb$model220genes)
modelftd$model220genes<-as.numeric(modelftd$model220genes)
# Add columns to code for number of apoe2 and apoe4 allelles
# for each subject in each data frame
#PD
modelpd$apoe4 = 0
modelpd$apoe4 [modelpd$APOE == 44] <- 2
modelpd$apoe4 [modelpd$APOE== 43] <- 1
modelpd$apoe4 [modelpd$APOE== 34 ] <- 1
modelpd$apoe4 [modelpd$APOE== 24] <- 1
modelpd$apoe4 [modelpd$APOE== 42] <- 1
modelpd$apoe2 = 0
modelpd$apoe2 [modelpd$APOE == 22] <- 2
modelpd$apoe2 [modelpd$APOE== 23] <- 1
modelpd$apoe2 [modelpd$APOE== 32 ] <- 1
modelpd$apoe2 [modelpd$APOE== 24] <- 1
#DLB
modeldlb$apoe4 = 0
modeldlb$apoe4 [modeldlb$APOE == 44] <- 2
modeldlb$apoe4 [modeldlb$APOE== 43] <- 1
modeldlb$apoe4 [modeldlb$APOE== 34 ] <- 1
modeldlb$apoe4 [modeldlb$APOE== 24] <- 1
modeldlb$apoe4 [modeldlb$APOE== 42] <- 1
modeldlb$apoe2 = 0
modeldlb$apoe2 [modeldlb$APOE == 22] <- 2
modeldlb$apoe2 [modeldlb$APOE== 23] <- 1
modeldlb$apoe2 [modeldlb$APOE== 32 ] <- 1
modeldlb$apoe2 [modeldlb$APOE== 24] <- 1
#FTD
modelftd$apoe4 = 0
modelftd$apoe4 [modelftd$APOE == 44] <- 2
modelftd$apoe4 [modelftd$APOE== 43] <- 1
modelftd$apoe4 [modelftd$APOE== 34 ] <- 1
modelftd$apoe4 [modelftd$APOE== 24] <- 1
modelftd$apoe4 [modelftd$APOE== 42] <- 1
modelftd$apoe2 = 0
modelftd$apoe2 [modelftd$APOE == 22] <- 2
modelftd$apoe2 [modelftd$APOE== 23] <- 1
modelftd$apoe2 [modelftd$APOE== 32 ] <- 1
modelftd$apoe2 [modelftd$APOE== 24] <- 1
# Fit the new datasets to the models, without APOE
#PD
fit220pd <- glm( Status ~  model220genes, data = modelpd,family=binomial())
preds = predict(fit220pd, type="response")
rocpd220 = as.data.frame(modelpd$Status)
rocpd220$predictor = preds
rocpd220$Classifier = "220 genes AD vs PD, without APOE (AUC = 0.69)"
names(rocpd220) <- c("Status", "Predictor", "Classifier")
fit90pd <- glm( Status ~  model90genes, data = modelpd,family=binomial())
preds = predict(fit90pd, type="response")
rocpd90 = as.data.frame(modelpd$Status)
rocpd90$predictor = preds
rocpd90$Classifier = "90 genes AD vs PD, without APOE (AUC = 0.65)"
names(rocpd90) <- c("Status", "Predictor", "Classifier")
fit40pd <- glm( Status ~  model40genes, data = modelpd,family=binomial())
preds = predict(fit40pd, type="response")
rocpd40 = as.data.frame(modelpd$Status)
rocpd40$predictor = preds
rocpd40$Classifier = "40 genes AD vs PD, without APOE (AUC = 0.55)"
names(rocpd40) <- c("Status", "Predictor", "Classifier")
#DLB
fit220dlb <- glm( Status ~  model220genes, data = modeldlb,family=binomial())
preds = predict(fit220dlb, type="response")
rocdlb220 = as.data.frame(modeldlb$Status)
rocdlb220$predictor = preds
rocdlb220$Classifier = "220 genes AD vs DLB, without APOE (AUC = 0.69)"
names(rocdlb220) <- c("Status", "Predictor", "Classifier")
fit90dlb <- glm( Status ~  model90genes, data = modeldlb,family=binomial())
preds = predict(fit90dlb, type="response")
rocdlb90 = as.data.frame(modeldlb$Status)
rocdlb90$predictor = preds
rocdlb90$Classifier = "90 genes AD vs DLB, without APOE (AUC = 0.65)"
names(rocdlb90) <- c("Status", "Predictor", "Classifier")
fit40dlb <- glm( Status ~  model40genes, data = modeldlb,family=binomial())
preds = predict(fit40dlb, type="response")
rocdlb40 = as.data.frame(modeldlb$Status)
rocdlb40$predictor = preds
rocdlb40$Classifier = "40 genes AD vs DLB, without APOE (AUC = 0.55)"
names(rocdlb40) <- c("Status", "Predictor", "Classifier")
#FTD
fit220ftd <- glm( Status ~  model220genes, data = modelftd,family=binomial())
preds = predict(fit220ftd, type="response")
rocftd220 = as.data.frame(modelftd$Status)
rocftd220$predictor = preds
rocftd220$Classifier = "220 genes AD vs FTD, without APOE (AUC = 0.76)"
names(rocftd220) <- c("Status", "Predictor", "Classifier")
fit90ftd <- glm( Status ~  model90genes, data = modelftd,family=binomial())
preds = predict(fit90ftd, type="response")
rocftd90 = as.data.frame(modelftd$Status)
rocftd90$predictor = preds
rocftd90$Classifier = "90 genes AD vs FTD, without APOE (AUC = 0.75)"
names(rocftd90) <- c("Status", "Predictor", "Classifier")
fit40ftd <- glm( Status ~  model40genes, data = modelftd,family=binomial())
preds = predict(fit40ftd, type="response")
rocftd40 = as.data.frame(modelftd$Status)
rocftd40$predictor = preds
rocftd40$Classifier = "40 genes AD vs FTD, without APOE (AUC = 0.86)"
names(rocftd40) <- c("Status", "Predictor", "Classifier")
# Fit the new datasets to the models, with APOE
#PD
fit220pd <- glm( Status ~  model220genes+apoe2+apoe4, data = modelpd,family=binomial())
preds = predict(fit220pd, type="response")
rocpd220apoe = as.data.frame(modelpd$Status)
rocpd220apoe$predictor = preds
rocpd220apoe$Classifier = "220 genes AD vs PD, with APOE (AUC = 0.72)"
names(rocpd220apoe) <- c("Status", "Predictor", "Classifier")
fit90pd <- glm( Status ~  model90genes+apoe2+apoe4, data = modelpd,family=binomial())
preds = predict(fit90pd, type="response")
rocpd90apoe = as.data.frame(modelpd$Status)
rocpd90apoe$predictor = preds
rocpd90apoe$Classifier = "90 genes AD vs PD, with APOE (AUC = 0.69)"
names(rocpd90apoe) <- c("Status", "Predictor", "Classifier")
fit40pd <- glm( Status ~  model40genes+apoe2+apoe4, data = modelpd,family=binomial())
preds = predict(fit40pd, type="response")
rocpd40apoe = as.data.frame(modelpd$Status)
rocpd40apoe$predictor = preds
rocpd40apoe$Classifier = "40 genes AD vs PD, with APOE (AUC = 0.65)"
names(rocpd40apoe) <- c("Status", "Predictor", "Classifier")
#DLB
fit220dlb <- glm( Status ~  model220genes+apoe2+apoe4, data = modeldlb,family=binomial())
preds = predict(fit220dlb, type="response")
rocdlb220apoe = as.data.frame(modeldlb$Status)
rocdlb220apoe$predictor = preds
rocdlb220apoe$Classifier = "220 genes AD vs DLB, with APOE (AUC = 0.72)"
names(rocdlb220apoe) <- c("Status", "Predictor", "Classifier")
fit90dlb <- glm( Status ~  model90genes+apoe2+apoe4, data = modeldlb,family=binomial())
preds = predict(fit90dlb, type="response")
rocdlb90apoe = as.data.frame(modeldlb$Status)
rocdlb90apoe$predictor = preds
rocdlb90apoe$Classifier = "90 genes AD vs DLB, with APOE (AUC = 0.69)"
names(rocdlb90apoe) <- c("Status", "Predictor", "Classifier")
fit40dlb <- glm( Status ~  model40genes+apoe2+apoe4, data = modeldlb,family=binomial())
preds = predict(fit40dlb, type="response")
rocdlb40apoe = as.data.frame(modeldlb$Status)
rocdlb40apoe$predictor = preds
rocdlb40apoe$Classifier = "40 genes AD vs DLB, with APOE (AUC = 0.65)"
names(rocdlb40apoe) <- c("Status", "Predictor", "Classifier")
#FTD
fit220ftd <- glm( Status ~  model220genes+apoe2+apoe4, data = modelftd,family=binomial())
preds = predict(fit220ftd, type="response")
rocftd220apoe = as.data.frame(modelftd$Status)
rocftd220apoe$predictor = preds
rocftd220apoe$Classifier = "220 genes AD vs FTD, with APOE (AUC = 0.80)"
names(rocftd220apoe) <- c("Status", "Predictor", "Classifier")
fit90ftd <- glm( Status ~  model90genes+apoe2+apoe4, data = modelftd,family=binomial())
preds = predict(fit90ftd, type="response")
rocftd90apoe = as.data.frame(modelftd$Status)
rocftd90apoe$predictor = preds
rocftd90apoe$Classifier = "90 genes AD vs FTD,  with APOE (AUC = 0.78)"
names(rocftd90apoe) <- c("Status", "Predictor", "Classifier")
fit40ftd <- glm( Status ~  model40genes+apoe2+apoe4, data = modelftd,family=binomial())
preds = predict(fit40ftd, type="response")
rocftd40apoe = as.data.frame(modelftd$Status)
rocftd40apoe$predictor = preds
rocftd40apoe$Classifier = "40 genes AD vs FTD, with APOE (AUC = 0.94)"
names(rocftd40apoe) <- c("Status", "Predictor", "Classifier")

#Get accuracy, specificity, positive and negative prediction value
#information foreach model and dataset
acuspec<-as.data.frame(t(caret::confusionMatrix(data = factor(as.numeric(modelpd$ model220genes>0.5)), reference = factor(modelpd$Status))$byClass))
acuspec<-rbind(acuspec,as.data.frame(t(caret::confusionMatrix(data = factor(as.numeric(modelpd$ model90genes>0.5)), reference = factor(modelpd$Status))$byClass)))
acuspec<-rbind(acuspec,as.data.frame(t(caret::confusionMatrix(data = factor(as.numeric(modelpd$ model40genes>0.5)), reference = factor(modelpd$Status))$byClass)))
acuspec<-rbind(acuspec,as.data.frame(t(caret::confusionMatrix(data = factor(as.numeric(modeldlb$ model220genes>0.5)), reference = factor(modeldlb$Status))$byClass)))
acuspec<-rbind(acuspec,as.data.frame(t(caret::confusionMatrix(data = factor(as.numeric(modeldlb$ model90genes>0.5)), reference = factor(modeldlb$Status))$byClass)))
acuspec<-rbind(acuspec,as.data.frame(t(caret::confusionMatrix(data = factor(as.numeric(modeldlb$ model40genes>0.5)), reference = factor(modeldlb$Status))$byClass)))
acuspec<-rbind(acuspec,as.data.frame(t(caret::confusionMatrix(data = factor(as.numeric(modelftd$ model220genes>0.5)), reference = factor(modelftd$Status))$byClass)))
acuspec<-rbind(acuspec,as.data.frame(t(caret::confusionMatrix(data = factor(as.numeric(modelftd$ model90genes>0.5)), reference = factor(modelftd$Status))$byClass)))
acuspec<-rbind(acuspec,as.data.frame(t(caret::confusionMatrix(data = factor(as.numeric(modelftd$ model40genes>0.5)), reference = factor(modelftd$Status))$byClass)))
acuspec<-rbind(acuspec,as.data.frame(t(caret::confusionMatrix(data = factor(as.numeric(rocpd220apoe$Predictor>0.5)), reference = factor(rocpd220apoe$Status))$byClass)))
acuspec<-rbind(acuspec,as.data.frame(t(caret::confusionMatrix(data = factor(as.numeric(rocpd90apoe$Predictor>0.5)), reference = factor(rocpd90apoe$Status))$byClass)))
acuspec<-rbind(acuspec,as.data.frame(t(caret::confusionMatrix(data = factor(as.numeric(rocpd40apoe$Predictor>0.5)), reference = factor(rocpd40apoe$Status))$byClass)))
acuspec<-rbind(acuspec,as.data.frame(t(caret::confusionMatrix(data = factor(as.numeric(rocdlb220apoe$Predictor>0.5)), reference = factor(rocdlb220apoe$Status))$byClass)))
acuspec<-rbind(acuspec,as.data.frame(t(caret::confusionMatrix(data = factor(as.numeric(rocdlb90apoe$Predictor>0.5)), reference = factor(rocdlb90apoe$Status))$byClass)))
acuspec<-rbind(acuspec,as.data.frame(t(caret::confusionMatrix(data = factor(as.numeric(rocdlb40apoe$Predictor>0.5)), reference = factor(rocdlb40apoe$Status))$byClass)))
acuspec<-rbind(acuspec,as.data.frame(t(caret::confusionMatrix(data = factor(as.numeric(rocftd220apoe$Predictor>0.5)), reference = factor(rocftd220apoe$Status))$byClass)))
acuspec<-rbind(acuspec,as.data.frame(t(caret::confusionMatrix(data = factor(as.numeric(rocftd90apoe$Predictor>0.5)), reference = factor(rocftd90apoe$Status))$byClass)))
acuspec<-rbind(acuspec,as.data.frame(t(caret::confusionMatrix(data = factor(as.numeric(rocftd40apoe$Predictor>0.5)), reference = factor(rocftd40apoe$Status))$byClass)))
acuspec$Dataset<-rep(c(rep("PDvAD",3),rep("DLBvAD",3),rep("FTDvAD",3)),2)
acuspec$APOE<-c(rep(0,9),rep(1,9))
acuspec$Gene_model<-rep(c(220,90,40),6)

##AUC whisker plot
t<-ci.auc(as.numeric(rocpd220$Status), as.numeric(rocpd220$Predictor),conf.level = 0.9)
aucdf<-data.frame(low = round(t[1],3), AUC = round(t[2],3), high = round(t[3],3), Classifier = "PD-220-noAPOE", Genes = "model with 220 genes", APOE = "without APOE")
t<-ci.auc(as.numeric(rocpd90$Status), as.numeric(rocpd90$Predictor),conf.level = 0.9)
t<-data.frame(low = round(t[1],3), AUC = round(t[2],3), high = round(t[3],3), Classifier = "PD-90-noAPOE",Genes = "model with 90 genes", APOE = "without APOE")
aucdf<-rbind(aucdf, t)
t<-ci.auc(as.numeric(rocpd40$Status), as.numeric(rocpd40$Predictor),conf.level = 0.9)
t<-data.frame(low = round(t[1],3), AUC = round(t[2],3), high = round(t[3],3), Classifier = "PD-40-noAPOE",Genes = "model with 40 genes", APOE = "without APOE")
aucdf<-rbind(aucdf, t)
t<-ci.auc(as.numeric(rocpd220apoe$Status), as.numeric(rocpd220apoe$Predictor),conf.level = 0.9)
t<-data.frame(low = round(t[1],3), AUC = round(t[2],3), high = round(t[3],3), Classifier = "PD-220-APOE",Genes = "model with 220 genes", APOE = "with APOE")
aucdf<-rbind(aucdf, t)
t<-ci.auc(as.numeric(rocpd90apoe$Status), as.numeric(rocpd90apoe$Predictor),conf.level = 0.9)
t<-data.frame(low = round(t[1],3), AUC = round(t[2],3), high = round(t[3],3), Classifier = "PD-90-APOE",Genes = "model with 90 genes", APOE = "with APOE")
aucdf<-rbind(aucdf, t)
t<-ci.auc(as.numeric(rocpd40apoe$Status), as.numeric(rocpd40apoe$Predictor),conf.level = 0.9)
t<-data.frame(low = round(t[1],3), AUC = round(t[2],3), high = round(t[3],3), Classifier = "PD-40-APOE",Genes = "model with 40 genes", APOE = "with APOE")
aucdf<-rbind(aucdf, t)
t<-ci.auc(as.numeric(rocdlb220$Status), as.numeric(rocdlb220$Predictor),conf.level = 0.9)
t<-data.frame(low = round(t[1],3), AUC = round(t[2],3), high = round(t[3],3), Classifier = "DLB-220-noAPOE", Genes = "model with 220 genes", APOE = "without APOE")
aucdf<-rbind(aucdf, t)
t<-ci.auc(as.numeric(rocdlb90$Status), as.numeric(rocdlb90$Predictor),conf.level = 0.9)
t<-data.frame(low = round(t[1],3), AUC = round(t[2],3), high = round(t[3],3), Classifier = "DLB-90-noAPOE",Genes = "model with 90 genes", APOE = "without APOE")
aucdf<-rbind(aucdf, t)
t<-ci.auc(as.numeric(rocdlb40$Status), as.numeric(rocdlb40$Predictor),conf.level = 0.9)
t<-data.frame(low = round(t[1],3), AUC = round(t[2],3), high = round(t[3],3), Classifier = "DLB-40-noAPOE",Genes = "model with 40 genes", APOE = "without APOE")
aucdf<-rbind(aucdf, t)
t<-ci.auc(as.numeric(rocdlb220apoe$Status), as.numeric(rocdlb220apoe$Predictor),conf.level = 0.9)
t<-data.frame(low = round(t[1],3), AUC = round(t[2],3), high = round(t[3],3), Classifier = "DLB-220-APOE",Genes = "model with 220 genes", APOE = "with APOE")
aucdf<-rbind(aucdf, t)
t<-ci.auc(as.numeric(rocdlb90apoe$Status), as.numeric(rocdlb90apoe$Predictor),conf.level = 0.9)
t<-data.frame(low = round(t[1],3), AUC = round(t[2],3), high = round(t[3],3), Classifier = "DLB-90-APOE",Genes = "model with 90 genes", APOE = "with APOE")
aucdf<-rbind(aucdf, t)
t<-ci.auc(as.numeric(rocdlb40apoe$Status), as.numeric(rocdlb40apoe$Predictor),conf.level = 0.9)
t<-data.frame(low = round(t[1],3), AUC = round(t[2],3), high = round(t[3],3), Classifier = "DLB-40-APOE",Genes = "model with 40 genes", APOE = "with APOE")
aucdf<-rbind(aucdf, t)
t<-ci.auc(as.numeric(rocftd220$Status), as.numeric(rocftd220$Predictor),conf.level = 0.9)
t<-data.frame(low = round(t[1],3), AUC = round(t[2],3), high = round(t[3],3), Classifier = "FTD-220-noAPOE", Genes = "model with 220 genes", APOE = "without APOE")
aucdf<-rbind(aucdf, t)
t<-ci.auc(as.numeric(rocftd90$Status), as.numeric(rocftd90$Predictor),conf.level = 0.9)
t<-data.frame(low = round(t[1],3), AUC = round(t[2],3), high = round(t[3],3), Classifier = "FTD-90-noAPOE", Genes = "model with 90 genes", APOE = "without APOE")
aucdf<-rbind(aucdf, t)
t<-ci.auc(as.numeric(rocftd40$Status), as.numeric(rocftd40$Predictor),conf.level = 0.9)
t<-data.frame(low = round(t[1],3), AUC = round(t[2],3), high = round(t[3],3), Classifier = "FTD-40-noAPOE", Genes = "model with 40 genes", APOE = "without APOE")
aucdf<-rbind(aucdf, t)
t<-ci.auc(as.numeric(rocftd220apoe$Status), as.numeric(rocftd220apoe$Predictor),conf.level = 0.9)
t<-data.frame(low = round(t[1],3), AUC = round(t[2],3), high = round(t[3],3), Classifier = "FTD-220-APOE", Genes = "model with 220 genes", APOE = "with APOE")
aucdf<-rbind(aucdf, t)
t<-ci.auc(as.numeric(rocftd90apoe$Status), as.numeric(rocftd90apoe$Predictor),conf.level = 0.9)
t<-data.frame(low = round(t[1],3), AUC = round(t[2],3), high = round(t[3],3), Classifier = "FTD-90-APOE", Genes = "model with 90 genes", APOE = "with APOE")
aucdf<-rbind(aucdf, t)
t<-ci.auc(as.numeric(rocftd40apoe$Status), as.numeric(rocftd40apoe$Predictor),conf.level = 0.9)
t<-data.frame(low = round(t[1],3), AUC = round(t[2],3), high = round(t[3],3), Classifier = "FTD-40-APOE", Genes = "model with 40 genes", APOE = "with APOE")
aucdf<-rbind(aucdf, t)
aucdf$Classifier<-factor(aucdf$Classifier, levels=c("PD-40-noAPOE", "PD-40-APOE","PD-90-noAPOE", "PD-90-APOE", 
                                                    "PD-220-noAPOE", "PD-220-APOE","DLB-40-noAPOE","DLB-40-APOE",
                                                    "DLB-90-noAPOE","DLB-90-APOE","DLB-220-noAPOE", "DLB-220-APOE", 
                                                    "FTD-40-noAPOE","FTD-40-APOE","FTD-90-noAPOE","FTD-90-APOE",
                                                    "FTD-220-noAPOE", "FTD-220-APOE"))  

whiskerPlot <- ggplot(aucdf, aes(x=AUC, y=Classifier,color = Genes))+
  geom_point(aes(shape=APOE), size =4)+
  scale_colour_manual(values=myPalette)+
  theme_bw()+
  theme(legend.position = "bottom")+
  scale_y_discrete(limits=rev)+
  scale_x_continuous("ROC AUC", limits = c(0.3,1), n.breaks = 9, expand = c(0,0)) +
  geom_errorbarh(aes(xmin=low,xmax=high),height=0)+
  geom_vline(xintercept = 0.5, linetype = "dashed")+
  guides(color = guide_legend(reverse = TRUE, override.aes = list(shape = 15)), shape = guide_legend(reverse = TRUE))
