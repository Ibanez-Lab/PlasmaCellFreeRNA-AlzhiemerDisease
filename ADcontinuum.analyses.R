#Load the necessary libraries
library(data.table)
library(dplyr)
library(readr)
library(ggplot2)
library(scales)
library(plotROC)
library(glmnet)
library(pROC)
myPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999", "#000000", "#6C2EB8")
# AD continuum
## Scale two sets function
#First, we need to rescale each column in the testing datasets for specifity analyses to the same range (Min, Max) in the training population.
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

## 220 genes
genestodoml220 = fread("220genes.csv")
genestodoml220 = genestodoml220$ENSG
genestodoml = genestodoml220
# Load Matrix
matrix_phase = readRDS("discovery_data")
train = matrix_phase
apoetrain<-train$apoe
# Select features
genestodoml = c(genestodoml, "Status")
commontrain = colnames(train) %in% genestodoml
train = train[,..commontrain]
# Select only genes
x_train <- train[,-1]
# Compute Z score
modelo_means <- as.data.frame(lapply(x_train, mean))
modelo_sd <- as.data.frame(lapply(x_train, sd))
### Preclinicals all
genestodoml = genestodoml220
genestodoml = c(genestodoml, "Status")
# Load Matrix
matrix_phase = readRDS("discovery_data")
train = matrix_phase
trainconpred = train
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
predicciones_train_220 <- predict(modelo, newx = x_train, type = "response", standardize = F)
genestodoml = genestodoml220
genestodoml = c(genestodoml, "Status")
# Load Matrix
test = readRDS("replication_data")
apoetest<-test$APOE
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
predicciones_test <- predict(modelo, newx = x_test, type = "response" , standarize = F, intercept = F)
predicciones_test_220 = predicciones_test
# Model220 roc comes from above code (220 genes model)
model220rocte = as.data.frame(y_test_220)
model220rocte$predictor = predicciones_test_220
model220rocte$APOE <-apoetest
names(model220rocte) = c("status", "predictor", "APOE")
model220roctr = as.data.frame(y_train_220)
model220roctr$predictor = predicciones_train_220
model220roctr$APOE <- apoetrain
names(model220roctr) = c("status", "predictor", "APOE")
model220rocadcontinuum = rbind(model220rocte,model220roctr)
# Extract control samples
controlsforspe = model220rocadcontinuum[model220rocadcontinuum$status == 0,]
### cdr 0.5 k99
genestodoml = genestodoml220
genestodoml = c(genestodoml, "Status")
# Load Matrix
test = readRDS("continuum_cdr0.5_data1")
test = as.data.frame(test)
apoetest<-test$apoe
# Select features
commontest = colnames(test) %in% genestodoml
test = test[,commontest]
nombrestest = as.data.frame(colnames(test))
colnames(nombrestest) = "gene"
genestodomldf<-as.data.frame(genestodoml)
names(genestodomldf)<-"gene"
genesnoaparecen = anti_join(genestodomldf,nombrestest)
genesnoaparecen = as.character(genesnoaparecen$gene)
test = readRDS("continuum_cdr0.5_data1")
test = as.data.frame(test)
commontest = colnames(test) %in% genestodoml
test = test[,commontest]
# Select only genes
x_test <- test[,-1]
# Compute Z score
x_test <- as.data.frame(lapply(x_test, function(x) (x - mean(x))/sd(x) ))
# add gene name from cdr05
x_test = cbind(x_test, setNames( lapply(genesnoaparecen, function(x) x=0), genesnoaparecen))
x_test <- as.matrix(x_test )
# Code status
numeros = as.numeric(test$Status)
y_test <-replace(numeros, numeros == 2,0)
y_test_220_cdr05k99 = y_test
# Testing preditions
predicciones_test <- predict(modelo, newx = x_test, type = "response" , standarize = F, intercept = F)
predicciones_test_220_cdr05k99 = predicciones_test
# Model220 roc comes from above code (220 genes model)
model220allcdr05k99 = as.data.frame(y_test_220_cdr05k99)
model220allcdr05k99$predictor = predicciones_test_220_cdr05k99
model220allcdr05k99$APOE <- apoetest
names(model220allcdr05k99) = c("status", "predictor", "APOE")
model220allcdr05k99 = model220allcdr05k99[model220allcdr05k99$status == 1,]
model220allcdr05k99 = rbind(model220allcdr05k99,controlsforspe)
### AD CDR 0.5 EXPANSION values
#Scale CDR 0.5 individuals
genestodoml = genestodoml220
# Load Matrix
matrix_phase = readRDS("discovery_data")
train = matrix_phase
apoetrain<-train$apoe
# Select features
genestodoml = c(genestodoml, "Status")
commontrain = colnames(train) %in% genestodoml
train = train[,..commontrain]
# Select only genes
x_train <- train[,-1]
test = readRDS("continuum_cdr0.5_data2")
genestodoml = c(genestodoml, "Status")
test = subset(test, test$Status == "CA")
apoetest<-test$APOE
commontest = colnames(test) %in% genestodoml
test = test[,commontest]
# Select only genes
x_test <- test[,-1]
cdr05scalated = scaletwosets(x_train, x_test)
genestodoml = c(genestodoml, "Status")
# Load Matrix
test = readRDS("continuum_cdr0.5_data2")
apoetest<-test$APOE
genesgenestodoml = as.data.frame(genestodoml)
colnames(genesgenestodoml) = "gene"
# Select features
commontest = colnames(test) %in% genestodoml
test = test[,commontest]
nombrestest = as.data.frame(colnames(test))
colnames(nombrestest) = "gene"
genesnoaparecen = anti_join(genesgenestodoml,nombrestest)
genesnoaparecen = as.character(genesnoaparecen$gene)
test = readRDS("continuum_cdr0.5_data2")
test = subset(test, test$Status == "CA")
apoetest<-test$APOE
commontest = colnames(test) %in% genestodoml
test = test[,commontest]
# Select only genes
x_test <-cdr05scalated
modelo_meansftdp = modelo_means
which(colnames(modelo_meansftdp) %in% genesnoaparecen)  
modelo_meansftd = modelo_meansftdp[,-192]
modelo_sdftdp = modelo_sd
modelo_sdftdpd = modelo_sdftdp[,-192]
x_testprueba = x_test[]-modelo_meansftd[col(x_test[])]
x_testpruebafinal = x_testprueba[]/modelo_sdftdpd[col(x_testprueba[])]
x_test = x_testpruebafinal
x_test = as.data.frame(x_test)
# add gene name from cdr05
x_test = cbind(x_test, setNames( lapply(genesnoaparecen, function(x) x=0), genesnoaparecen))
x_test <- as.matrix(x_test )
protoeatmap = as.data.frame(x_test)
# Code status
numeros = as.numeric(test$Status)
y_test <-replace(numeros, numeros == 2,0)
y_test_cdr5exp =  y_test
# Testing preditions
predicciones_test <- predict(modelo, newx = x_test, type = "response" , standarize = F, intercept = F) 
predicciones_test_cdr5exp = predicciones_test
model220cdr5exp = as.data.frame(y_test_cdr5exp)
model220cdr5exp$predictor = predicciones_test_cdr5exp
model220cdr5exp$APOE <-apoetest
names(model220cdr5exp) = c("status", "predictor", "APOE")
model220allcdr05roc = rbind(model220allcdr05k99,model220cdr5exp)
### cdr 1 k99 (Esquema cdr0.5)
genestodoml = genestodoml220
genestodoml = c(genestodoml, "Status")
# Load Matrix
test = readRDS("continuum_cdr1_data1")
test = as.data.frame(test)
apoetest<-test$apoe
# Select features
commontest = colnames(test) %in% genestodoml
test = test[,commontest]
nombrestest = as.data.frame(colnames(test))
colnames(nombrestest) = "gene"
genestodomldf<-as.data.frame(genestodoml)
names(genestodomldf)<-"gene"
genesnoaparecen = anti_join(genestodomldf,nombrestest)
genesnoaparecen = as.character(genesnoaparecen$gene)
test = readRDS("continuum_cdr1_data1")
test = as.data.frame(test)
apoetest<-test$apoe
commontest = colnames(test) %in% genestodoml
test = test[,commontest]
# Select only genes
x_test <- test[,-1]
# Compute Z score
x_test <- as.data.frame(lapply(x_test, function(x) (x - mean(x))/sd(x) ))
# add gene name from cdr05
x_test = cbind(x_test, setNames( lapply(genesnoaparecen, function(x) x=0), genesnoaparecen))
x_test <- as.matrix(x_test )
# Code status
numeros = as.numeric(test$Status)
y_test <-replace(numeros, numeros == 2,0)
y_test_220_cdr1k99 = y_test
# Testing preditions
predicciones_test <- predict(modelo, newx = x_test, type = "response" , standarize = F, intercept = F)
predicciones_test_220_cdr1k99 = predicciones_test
# Model220 roc comes from above code (220 genes model)
model220allcdr1k99 = as.data.frame(y_test_220_cdr1k99)
model220allcdr1k99$predictor = predicciones_test_220_cdr1k99
model220allcdr1k99$APOE <- apoetest
names(model220allcdr1k99) = c("status", "predictor", "APOE")
model220allcdr1k99 = model220allcdr1k99[model220allcdr1k99$status == 1,]
model220allcdr1k99 = rbind(model220allcdr1k99,controlsforspe)
### AD CDR 1 EXPANSION values
#cale CDR 1 individuals
genestodoml = genestodoml220
# Load Matrix
matrix_phase = readRDS("discovery_data")
train = matrix_phase
apoetrain<-train$apoe
# Select features
genestodoml = c(genestodoml, "Status")
commontrain = colnames(train) %in% genestodoml
train = train[,..commontrain]
# Select only genes
x_train <- train[,-1]
test = readRDS("continuum_cdr1_data2")
genestodoml = c(genestodoml, "Status")
test = subset(test, test$Status == "CA")
apoetest<-test$APOE
commontest = colnames(test) %in% genestodoml
test = test[,commontest]
# Select only genes
x_test <- test[,-1]
cdr1scalated = scaletwosets(x_train, x_test)
test = readRDS("continuum_cdr1_data2")
test = subset(test, test$Status == "CA")
apoetest<-test$APOE
commontest = colnames(test) %in% genestodoml
test = test[,commontest]
# Select only genes
x_test <-cdr1scalated
x_testprueba = x_test[]-modelo_means[col(x_test[])]
x_testpruebafinal = x_testprueba[]/modelo_sd[col(x_testprueba[])]
x_test = x_testpruebafinal
x_test = as.data.frame(x_test)
# add gene name from cdr05
#x_test = cbind(x_test, setNames( lapply(genesnoaparecen, function(x) x=0), genesnoaparecen) )
x_test <- as.matrix(x_test )
clinictoeatmap = as.data.frame(x_test)
# Code status
numeros = as.numeric(test$Status)
y_test <-replace(numeros, numeros == 2,0)
y_test_cdr1exp = y_test
# Testing preditions
predicciones_test <- predict(modelo, newx = x_test, type = "response" , standarize = F, intercept = F) 
predicciones_test_cdr1exp = predicciones_test
model220cdr1exp = as.data.frame(y_test_cdr1exp)
model220cdr1exp$predictor = predicciones_test_cdr1exp
model220cdr1exp$APOE <-apoetest
names(model220cdr1exp) = c("status", "predictor", "APOE")
model220allcdr1roc = rbind(model220allcdr1k99,model220cdr1exp)
## 40 genes
genestodoml40 = fread("40genes.csv")
genestodoml40 = genestodoml40[-1,1]
genestodoml40 = genestodoml40$V1
genestodoml = genestodoml40
# Load Matrix
matrix_phase = readRDS("discovery_data")
train = matrix_phase
apoetrain<-train$apoe
# Select features
genestodoml = c(genestodoml, "Status")
commontrain = colnames(train) %in% genestodoml
train = train[,..commontrain]
# Select only genes
x_train <- train[,-1]
# Compute Z score
modelo_means <- as.data.frame(lapply(x_train, mean))
modelo_sd <- as.data.frame(lapply(x_train, sd))
genestodoml = genestodoml40
genestodoml = c(genestodoml40, "Status")
# Load Matrix
matrix_phase = readRDS("discovery_data")
train = matrix_phase
apoetrain<-train$apoe
trainconpred = train
# Select features
commontrain = colnames(train) %in% genestodoml
train = train[,..commontrain]
# Select only genes
x_train <- train[,-1]
# Compute Z score
x_train <- as.data.frame(lapply(x_train, function(x) (x - mean(x))/sd(x) ))
x_train <- as.matrix(x_train)
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
predicciones_train_40 <- predict(modelo, newx = x_train, type = "response", standardize = F)
genestodoml = genestodoml40
genestodoml = c(genestodoml, "Status")
# Load Matrix
test = readRDS("replication_data")
apoetest<-test$APOE
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
y_test_40 = y_test
# Testing preditions
predicciones_test <- predict(modelo, newx = x_test, type = "response" , standarize = F, intercept = F)
predicciones_test_40 = predicciones_test
# Model220 roc comes from above code (220 genes model)
model40rocte = as.data.frame(y_test_40)
model40rocte$predictor = predicciones_test_40
model40rocte$APOE <- apoetest
names(model40rocte) = c("status", "predictor", "APOE")
model40roctr = as.data.frame(y_train_40)
model40roctr$predictor = predicciones_train_40
model40roctr$APOE <- apoetrain
names(model40roctr) = c("status", "predictor", "APOE")
model40rocadcontinuum = rbind(model40rocte,model40roctr)
# Extract control samples
controlsforspe = model40rocadcontinuum[model40rocadcontinuum$status == 0,]
# Add this dataset to the final
### cdr 0.5 k99
genestodoml = genestodoml40
genesgenestodoml = as.data.frame(genestodoml)
colnames(genesgenestodoml) = "gene"
genesgenestodoml = as.data.frame(genestodoml)
colnames(genesgenestodoml) = "gene"
genestodoml = c(genestodoml, "Status")
# Load Matrix
test = readRDS("continuum_cdr0.5_data1")
test = as.data.frame(test)
apoetest<-test$apoe
# Select features
commontest = colnames(test) %in% genestodoml
test = test[,commontest]
nombrestest = as.data.frame(colnames(test))
colnames(nombrestest) = "gene"
genesnoaparecen = anti_join(genesgenestodoml,nombrestest)
genesnoaparecen = as.character(genesnoaparecen$gene)
test = readRDS("continuum_cdr0.5_data1")
test = as.data.frame(test)
apoetest<-test$apoe
commontest = colnames(test) %in% genestodoml
test = test[,commontest]
# Select only genes
x_test <- test[,-1]
# Compute Z score
x_test <- as.data.frame(lapply(x_test, function(x) (x - mean(x))/sd(x) ))
x_test <- as.matrix(x_test )
# Code status
numeros = as.numeric(test$Status)
y_test <-replace(numeros, numeros == 2,0)
y_test_40_cdr05k99 = y_test
# Testing preditions
predicciones_test <- predict(modelo, newx = x_test, type = "response" , standarize = F, intercept = F)
predicciones_test_40_cdr05k99 = predicciones_test
# Model40 roc comes from above code (40 genes model)
model40allcdr05k99 = as.data.frame(y_test_40_cdr05k99)
model40allcdr05k99$predictor = predicciones_test_40_cdr05k99
model40allcdr05k99$APOE <- apoetest
names(model40allcdr05k99) = c("status", "predictor", "APOE")
model40allcdr05k99 = model40allcdr05k99[model40allcdr05k99$status == 1,]
model40allcdr05k99 = rbind(model40allcdr05k99,controlsforspe)
### AD CDR 0.5 EXPANSION values
#Scale CDR 0.5 individuals
genestodoml = genestodoml40
# Load Matrix
matrix_phase = readRDS("discovery_data")
train = matrix_phase
apoetrain<-train$apoe
# Select features
genestodoml = c(genestodoml, "Status")
commontrain = colnames(train) %in% genestodoml
train = train[,..commontrain]
# Select only genes
x_train <- train[,-1]
test = readRDS("continuum_cdr0.5_data2")
genestodoml = c(genestodoml, "Status")
test = subset(test, test$Status == "CA")
apoetest<-test$APOE
commontest = colnames(test) %in% genestodoml
test = test[,commontest]
# Select only genes
x_test <- test[,-1]
cdr05scalated = scaletwosets(x_train, x_test)
genestodoml = c(genestodoml, "Status")
# Load Matrix
test = readRDS("continuum_cdr0.5_data2")
apotest<-test$APOE
genesgenestodoml = as.data.frame(genestodoml)
colnames(genesgenestodoml) = "gene"
# Select features
commontest = colnames(test) %in% genestodoml
test = test[,commontest]
nombrestest = as.data.frame(colnames(test))
colnames(nombrestest) = "gene"
genestodomldf<-as.data.frame(genestodoml)
names(genestodomldf)<-"gene"
genesnoaparecen = anti_join(genestodomldf,nombrestest)
genesnoaparecen = as.character(genesnoaparecen$gene)
test = readRDS("continuum_cdr0.5_data2")
test = subset(test, test$Status == "CA")
apoetest<-test$APOE
commontest = colnames(test) %in% genestodoml
test = test[,commontest]
# Select only genes
x_test <-cdr05scalated
modelo_meansftdp = modelo_means
which(colnames(modelo_meansftdp) %in% genesnoaparecen)  
# Aqui - el que no es gen
modelo_meansftd = modelo_meansftdp[,]
modelo_sdftdp = modelo_sd
modelo_sdftdpd = modelo_sdftdp[,]
x_testprueba = x_test[]-modelo_meansftd[col(x_test[])]
x_testpruebafinal = x_testprueba[]/modelo_sdftdpd[col(x_testprueba[])]
x_test = x_testpruebafinal
x_test = as.data.frame(x_test)
# add gene name from cdr05
x_test <- as.matrix(x_test )
# Code status
numeros = as.numeric(test$Status)
y_test <-replace(numeros, numeros == 2,0)
y_test_cdr5exp =  y_test
# Testing preditions
predicciones_test <- predict(modelo, newx = x_test, type = "response" , standarize = F, intercept = F) 
predicciones_test_cdr5exp = predicciones_test
model40cdr5exp = as.data.frame(y_test_cdr5exp)
model40cdr5exp$predictor = predicciones_test_cdr5exp
model40cdr5exp$APOE <- apoetest
names(model40cdr5exp) = c("status", "predictor", "APOE")
model40allcdr05roc = rbind(model40allcdr05k99,model40cdr5exp)
### cdr 1 k99 (Esquema cdr0.5)
genestodoml = genestodoml40
genestodoml = c(genestodoml, "Status")
# Load Matrix
test = readRDS("continuum_cdr1_data1")
test = as.data.frame(test)
apoetest<-test$apoe
# Select features
commontest = colnames(test) %in% genestodoml
test = test[,commontest]
nombrestest = as.data.frame(colnames(test))
colnames(nombrestest) = "gene"
genesnoaparecen = anti_join(genesgenestodoml,nombrestest)
genesnoaparecen = as.character(genesnoaparecen$gene)
test = readRDS("continuum_cdr1_data1")
test = as.data.frame(test)
apoetest<-test$apoe
commontest = colnames(test) %in% genestodoml
test = test[,commontest]
# Select only genes
x_test <- test[,-1]
# Compute Z score
x_test <- as.data.frame(lapply(x_test, function(x) (x - mean(x))/sd(x) ))
# add gene name from cdr05
x_test <- as.matrix(x_test )
# Code status
numeros = as.numeric(test$Status)
y_test <-replace(numeros, numeros == 2,0)
y_test_40_cdr1k99 = y_test
# Testing preditions
predicciones_test <- predict(modelo, newx = x_test, type = "response" , standarize = F, intercept = F)
predicciones_test_40_cdr1k99 = predicciones_test
# Model40 roc comes from above code (40 genes model)
model40allcdr1k99 = as.data.frame(y_test_40_cdr1k99)
model40allcdr1k99$predictor = predicciones_test_40_cdr1k99
model40allcdr1k99$APOE <-apoetest
names(model40allcdr1k99) = c("status", "predictor", "APOE")
model40allcdr1k99 = model40allcdr1k99[model40allcdr1k99$status == 1,]
model40allcdr1k99 = rbind(model40allcdr1k99,controlsforspe)
### AD CDR 1 EXPANSION values
#Scale CDR 1 individuals
genestodoml = genestodoml40
# Load Matrix
matrix_phase = readRDS("discovery_data")
train = matrix_phase
apoetrain<-train$apoe
# Select features
genestodoml = c(genestodoml, "Status")
commontrain = colnames(train) %in% genestodoml
train = train[,..commontrain]
# Select only genes
x_train <- train[,-1]
test = readRDS("continuum_cdr1_data2")
genestodoml = c(genestodoml, "Status")
test = subset(test, test$Status == "CA")
apotest<-test$APOE
commontest = colnames(test) %in% genestodoml
test = test[,commontest]
# Select only genes
x_test <- test[,-1]
cdr1scalated = scaletwosets(x_train, x_test)
test = readRDS("continuum_cdr1_data2")
test = subset(test, test$Status == "CA")
apoetest<-test$APOE
commontest = colnames(test) %in% genestodoml
test = test[,commontest]
# Select only genes
x_test <-cdr1scalated
x_testprueba = x_test[]-modelo_means[col(x_test[])]
x_testpruebafinal = x_testprueba[]/modelo_sd[col(x_testprueba[])]
x_test = x_testpruebafinal
x_test = as.data.frame(x_test)
# add gene name from cdr05
x_test <- as.matrix(x_test )
clinictoeatmap = as.data.frame(x_test)
# Code status
numeros = as.numeric(test$Status)
y_test <-replace(numeros, numeros == 2,0)
y_test_cdr1exp = y_test
# Testing preditions
predicciones_test <- predict(modelo, newx = x_test, type = "response" , standarize = F, intercept = F) 
predicciones_test_cdr1exp = predicciones_test
model40cdr1exp = as.data.frame(y_test_cdr1exp)
model40cdr1exp$predictor = predicciones_test_cdr1exp
model40cdr1exp$APOE = apoetest
names(model40cdr1exp) = c("status", "predictor", "APOE")
model40allcdr1roc = rbind(model40allcdr1k99,model40cdr1exp)
## 90 genes
genestodoml90 = fread("90genes.csv")
genestodoml90 = genestodoml90[-1,1]
genestodoml90 = genestodoml90$V1
genestodoml = genestodoml90
# Load Matrix
matrix_phase = readRDS("discovery_data")
train = matrix_phase
apoetrain<-train$apoe
# Select features
genestodoml = c(genestodoml, "Status")
commontrain = colnames(train) %in% genestodoml
train = train[,..commontrain]
# Select only genes
x_train <- train[,-1]
# Compute Z score
modelo_means <- as.data.frame(lapply(x_train, mean))
modelo_sd <- as.data.frame(lapply(x_train, sd))
genestodoml = genestodoml90
genestodoml = c(genestodoml90, "Status")
# Load Matrix
matrix_phase = readRDS("discovery_data")
train = matrix_phase
apoetrain<-train$apoe
trainconpred = train
# Select features
commontrain = colnames(train) %in% genestodoml
train = train[,..commontrain]
# Select only genes
x_train <- train[,-1]
# Compute Z score
x_train <- as.data.frame(lapply(x_train, function(x) (x - mean(x))/sd(x) ))
x_train <- as.matrix(x_train)
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
predicciones_train_90 <- predict(modelo, newx = x_train, type = "response", standardize = F)
genestodoml = genestodoml90
genestodoml = c(genestodoml, "Status")
# Load Matrix
test = readRDS("replication_data")
apoetest<-test$APOE
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
y_test_90 = y_test
# Testing preditions
predicciones_test <- predict(modelo, newx = x_test, type = "response" , standarize = F, intercept = F)
predicciones_test_90 = predicciones_test
# Model220 roc comes from above code (220 genes model)
model90rocte = as.data.frame(y_test_90)
model90rocte$predictor = predicciones_test_90
model90rocte$APOE = apoetest
names(model90rocte) = c("status", "predictor", "APOE")
model90roctr = as.data.frame(y_train_90)
model90roctr$predictor = predicciones_train_90
model90roctr$APOE = apoetrain
names(model90roctr) = c("status", "predictor", "APOE")
model90rocadcontinuum = rbind(model90rocte,model90roctr)
# Extract control samples
controlsforspe = model90rocadcontinuum[model90rocadcontinuum$status == 0,]
# Add this dataset to the final
### cdr 0.5 k99
genestodoml = genestodoml90
genesgenestodoml = as.data.frame(genestodoml)
colnames(genesgenestodoml) = "gene"
genesgenestodoml = as.data.frame(genestodoml)
colnames(genesgenestodoml) = "gene"
genestodoml = c(genestodoml, "Status")
# Load Matrix
test = readRDS("continuum_cdr0.5_data1")
test = as.data.frame(test)
apoetest<-test$apoe
# Select features
commontest = colnames(test) %in% genestodoml
test = test[,commontest]
nombrestest = as.data.frame(colnames(test))
colnames(nombrestest) = "gene"
genesnoaparecen = anti_join(genesgenestodoml,nombrestest)
genesnoaparecen = as.character(genesnoaparecen$gene)
test = readRDS("continuum_cdr0.5_data1")
test = as.data.frame(test)
apoetest<-test$apoe
commontest = colnames(test) %in% genestodoml
test = test[,commontest]
# Select only genes
x_test <- test[,-1]
# Compute Z score
x_test <- as.data.frame(lapply(x_test, function(x) (x - mean(x))/sd(x) ))
# add gene name from cdr05
x_test = cbind(x_test, setNames( lapply(genesnoaparecen, function(x) x=0), genesnoaparecen))
x_test <- as.matrix(x_test )
# Code status
numeros = as.numeric(test$Status)
y_test <-replace(numeros, numeros == 2,0)
y_test_90_cdr05k99 = y_test
# Testing preditions
predicciones_test <- predict(modelo, newx = x_test, type = "response" , standarize = F, intercept = F)
predicciones_test_90_cdr05k99 = predicciones_test
# Model90 roc comes from above code (90 genes model)
model90allcdr05k99 = as.data.frame(y_test_90_cdr05k99)
model90allcdr05k99$predictor = predicciones_test_90_cdr05k99
model90allcdr05k99$APOE = apoetest
names(model90allcdr05k99) = c("status", "predictor", "APOE")
model90allcdr05k99 = model90allcdr05k99[model90allcdr05k99$status == 1,]
model90allcdr05k99 = rbind(model90allcdr05k99,controlsforspe)
### AD CDR 0.5 EXPANSION values
#Scale CDR 0.5 individuals
genestodoml = genestodoml90
# Load Matrix
matrix_phase = readRDS("discovery_data")
train = matrix_phase
apoetrain<-train$apoe
# Select features
genestodoml = c(genestodoml, "Status")
commontrain = colnames(train) %in% genestodoml
train = train[,..commontrain]
# Select only genes
x_train <- train[,-1]
test = readRDS("continuum_cdr0.5_data2")
genestodoml = c(genestodoml, "Status")
test = subset(test, test$Status == "CA")
apotest<-test$APOE
commontest = colnames(test) %in% genestodoml
test = test[,commontest]
# Select only genes
x_test <- test[,-1]
cdr05scalated = scaletwosets(x_train, x_test)
genestodoml = c(genestodoml, "Status")
# Load Matrix
test = readRDS("continuum_cdr0.5_data2")
apoetest<-test$APOE
genesgenestodoml = as.data.frame(genestodoml)
colnames(genesgenestodoml) = "gene"
# Select features
commontest = colnames(test) %in% genestodoml
test = test[,commontest]
nombrestest = as.data.frame(colnames(test))
colnames(nombrestest) = "gene"
genestodomldf<-as.data.frame(genestodoml)
names(genestodomldf)<-"gene"
genesnoaparecen = anti_join(genestodomldf,nombrestest)
genesnoaparecen = as.character(genesnoaparecen$gene)
test = readRDS("continuum_cdr0.5_data2")
test = subset(test, test$Status == "CA")
apoetest<-test$APOE
commontest = colnames(test) %in% genestodoml
test = test[,commontest]
# Select only genes
x_test <-cdr05scalated
modelo_meansftdp = modelo_means
which(colnames(modelo_meansftdp) %in% genesnoaparecen)  
# Aqui - el que no es gen
modelo_meansftd = modelo_meansftdp[,]
modelo_sdftdp = modelo_sd
modelo_sdftdpd = modelo_sdftdp[,]
x_testprueba = x_test[]-modelo_meansftd[col(x_test[])]
x_testpruebafinal = x_testprueba[]/modelo_sdftdpd[col(x_testprueba[])]
x_test = x_testpruebafinal
x_test = as.data.frame(x_test)
# add gene name from cdr05
x_test <- as.matrix(x_test )
# Code status
numeros = as.numeric(test$Status)
y_test <-replace(numeros, numeros == 2,0)
y_test_cdr5exp =  y_test
# Testing preditions
predicciones_test <- predict(modelo, newx = x_test, type = "response" , standarize = F, intercept = F) 
predicciones_test_cdr5exp = predicciones_test
model90cdr5exp = as.data.frame(y_test_cdr5exp)
model90cdr5exp$predictor = predicciones_test_cdr5exp
model90cdr5exp$APOE <-apoetest
names(model90cdr5exp) = c("status", "predictor", "APOE")
model90allcdr05roc = rbind(model90allcdr05k99,model90cdr5exp)
### cdr 1 k99 (Esquema cdr0.5)
genestodoml = genestodoml90
genestodoml = c(genestodoml, "Status")
# Load Matrix
test = readRDS("continuum_cdr1_data1")
test = as.data.frame(test)
apoetest<-test$apoe
# Select features
commontest = colnames(test) %in% genestodoml
test = test[,commontest]
nombrestest = as.data.frame(colnames(test))
colnames(nombrestest) = "gene"
genesnoaparecen = anti_join(genesgenestodoml,nombrestest)
genesnoaparecen = as.character(genesnoaparecen$gene)
test = readRDS("continuum_cdr1_data1")
test = as.data.frame(test)
apoetest<-test$apoe
commontest = colnames(test) %in% genestodoml
test = test[,commontest]
# Select only genes
x_test <- test[,-1]
# Compute Z score
x_test <- as.data.frame(lapply(x_test, function(x) (x - mean(x))/sd(x) ))
# add gene name from cdr05
x_test = cbind(x_test, setNames( lapply(genesnoaparecen, function(x) x=0), genesnoaparecen))
x_test <- as.matrix(x_test )
# Code status
numeros = as.numeric(test$Status)
y_test <-replace(numeros, numeros == 2,0)
y_test_90_cdr1k99 = y_test
# Testing preditions
predicciones_test <- predict(modelo, newx = x_test, type = "response" , standarize = F, intercept = F)
predicciones_test_90_cdr1k99 = predicciones_test
# Model90 roc comes from above code (90 genes model)
model90allcdr1k99 = as.data.frame(y_test_90_cdr1k99)
model90allcdr1k99$predictor = predicciones_test_90_cdr1k99
model90allcdr1k99$APOE <-apoetest
names(model90allcdr1k99) = c("status", "predictor", "APOE")
model90allcdr1k99 = model90allcdr1k99[model90allcdr1k99$status == 1,]
model90allcdr1k99 = rbind(model90allcdr1k99,controlsforspe)
### AD CDR 1 EXPANSION values
#Scale CDR 1 individuals
genestodoml = genestodoml90
# Load Matrix
matrix_phase = readRDS("discovery_data")
train = matrix_phase
apotr<-train$apoe
# Select features
genestodoml = c(genestodoml, "Status")
commontrain = colnames(train) %in% genestodoml
train = train[,..commontrain]
# Select only genes
x_train <- train[,-1]
test = readRDS("continuum_cdr1_data2")
genestodoml = c(genestodoml, "Status")
test = subset(test, test$Status == "CA")
apoetest<-test$APOE
commontest = colnames(test) %in% genestodoml
test = test[,commontest]
# Select only genes
x_test <- test[,-1]
cdr1scalated = scaletwosets(x_train, x_test)
test = readRDS("continuum_cdr1_data2")
test = subset(test, test$Status == "CA")
apotest<-test$APOE
commontest = colnames(test) %in% genestodoml
test = test[,commontest]
# Select only genes
x_test <-cdr1scalated
x_testprueba = x_test[]-modelo_means[col(x_test[])]
x_testpruebafinal = x_testprueba[]/modelo_sd[col(x_testprueba[])]
x_test = x_testpruebafinal
x_test = as.data.frame(x_test)
# add gene name from cdr05
x_test <- as.matrix(x_test )
clinictoeatmap = as.data.frame(x_test)
# Code status
numeros = as.numeric(test$Status)
y_test <-replace(numeros, numeros == 2,0)
y_test_cdr1exp = y_test
# Testing preditions
predicciones_test <- predict(modelo, newx = x_test, type = "response" , standarize = F, intercept = F) 
predicciones_test_cdr1exp = predicciones_test
model90cdr1exp = as.data.frame(y_test_cdr1exp)
model90cdr1exp$predictor = predicciones_test_cdr1exp
model90cdr1exp$APOE <- apoetest
names(model90cdr1exp) = c("status", "predictor", "APOE")
model90allcdr1roc = rbind(model90allcdr1k99,model90cdr1exp)
#combine data into 1 data.frame per disease status
modelcontr <- cbind(model220rocadcontinuum, as.data.frame(model90rocadcontinuum$predictor))
modelcontr <-cbind(modelcontr, as.data.frame(model40rocadcontinuum$predictor))
names(modelcontr)<-c("Status","model220genes","APOE", "model90genes", "model40genes")
modelcdr05 <- cbind(model220allcdr05roc, as.data.frame(model90allcdr05roc$predictor))
modelcdr05 <- cbind(modelcdr05, as.data.frame(model40allcdr05roc$predictor))
names(modelcdr05)<-c("Status","model220genes","APOE", "model90genes", "model40genes")
modelcdr1 <- cbind(model220allcdr1roc, as.data.frame(model90allcdr1roc$predictor))
modelcdr1 <-cbind(modelcdr1, as.data.frame(model40allcdr1roc$predictor))
names(modelcdr1)<-c("Status","model220genes","APOE", "model90genes", "model40genes")
# Add columns to code for number of apoe2 and apoe4 allelles
# for each subject in each data frame
#Presymptomatic
modelcontr$apoe4 = 0
modelcontr$apoe4 [modelcontr$APOE == 44] <- 2
modelcontr$apoe4 [modelcontr$APOE== 43] <- 1
modelcontr$apoe4 [modelcontr$APOE== 34 ] <- 1
modelcontr$apoe4 [modelcontr$APOE== 24] <- 1
modelcontr$apoe4 [modelcontr$APOE== 42] <- 1
modelcontr$apoe2 = 0
modelcontr$apoe2 [modelcontr$APOE == 22] <- 2
modelcontr$apoe2 [modelcontr$APOE== 23] <- 1
modelcontr$apoe2 [modelcontr$APOE== 32 ] <- 1
modelcontr$apoe2 [modelcontr$APOE== 24] <- 1
#Early Symptomatic
modelcdr05$apoe4 = 0
modelcdr05$apoe4 [modelcdr05$APOE == 44] <- 2
modelcdr05$apoe4 [modelcdr05$APOE== 43] <- 1
modelcdr05$apoe4 [modelcdr05$APOE== 34 ] <- 1
modelcdr05$apoe4 [modelcdr05$APOE== 24] <- 1
modelcdr05$apoe4 [modelcdr05$APOE== 42] <- 1
modelcdr05$apoe2 = 0
modelcdr05$apoe2 [modelcdr05$APOE == 22] <- 2
modelcdr05$apoe2 [modelcdr05$APOE== 23] <- 1
modelcdr05$apoe2 [modelcdr05$APOE== 32 ] <- 1
modelcdr05$apoe2 [modelcdr05$APOE== 24] <- 1
#Symptomatic
modelcdr1$apoe4 = 0
modelcdr1$apoe4 [modelcdr1$APOE == 44] <- 2
modelcdr1$apoe4 [modelcdr1$APOE== 43] <- 1
modelcdr1$apoe4 [modelcdr1$APOE== 34 ] <- 1
modelcdr1$apoe4 [modelcdr1$APOE== 24] <- 1
modelcdr1$apoe4 [modelcdr1$APOE== 42] <- 1
modelcdr1$apoe2 = 0
modelcdr1$apoe2 [modelcdr1$APOE == 22] <- 2
modelcdr1$apoe2 [modelcdr1$APOE== 23] <- 1
modelcdr1$apoe2 [modelcdr1$APOE== 32 ] <- 1
modelcdr1$apoe2 [modelcdr1$APOE== 24] <- 1

# Fit the new datasets to the models, without APOE
#Presymptomatic
fit220contr <- glm( Status ~  model220genes, data = modelcontr,family=binomial())
preds = predict(fit220contr, type="response")
roccontr220 = as.data.frame(modelcontr$Status)
roccontr220$predictor = preds
roccontr220$Classifier = "220 genes Presymptomatic, without APOE"
names(roccontr220) <- c("Status", "Predictor", "Classifier")
fit90contr <- glm( Status ~  model90genes, data = modelcontr,family=binomial())
preds = predict(fit90contr, type="response")
roccontr90 = as.data.frame(modelcontr$Status)
roccontr90$predictor = preds
roccontr90$Classifier = "90 genes Presymptomatic, without APOE"
names(roccontr90) <- c("Status", "Predictor", "Classifier")
fit40contr <- glm( Status ~  model40genes, data = modelcontr,family=binomial())
preds = predict(fit40contr, type="response")
roccontr40 = as.data.frame(modelcontr$Status)
roccontr40$predictor = preds
roccontr40$Classifier = "40 genes Presymptomatic, without APOE"
names(roccontr40) <- c("Status", "Predictor", "Classifier")
#Early Symptomatic
fit220cdr05 <- glm( Status ~  model220genes, data = modelcdr05,family=binomial())
preds = predict(fit220cdr05, type="response")
roccdr05220 = as.data.frame(modelcdr05$Status)
roccdr05220$predictor = preds
roccdr05220$Classifier = "220 genes Early Symptomatic, without APOE"
names(roccdr05220) <- c("Status", "Predictor", "Classifier")
fit90cdr05 <- glm( Status ~  model90genes, data = modelcdr05,family=binomial())
preds = predict(fit90cdr05, type="response")
roccdr0590 = as.data.frame(modelcdr05$Status)
roccdr0590$predictor = preds
roccdr0590$Classifier = "90 genes Early Symptomatic, without APOE"
names(roccdr0590) <- c("Status", "Predictor", "Classifier")
fit40cdr05 <- glm( Status ~  model40genes, data = modelcdr05,family=binomial())
preds = predict(fit40cdr05, type="response")
roccdr0540 = as.data.frame(modelcdr05$Status)
roccdr0540$predictor = preds
roccdr0540$Classifier = "40 genes Early Symptomatic, without APOE"
names(roccdr0540) <- c("Status", "Predictor", "Classifier")
#Symptomatic
fit220cdr1 <- glm( Status ~  model220genes, data = modelcdr1,family=binomial())
preds = predict(fit220cdr1, type="response")
roccdr1220 = as.data.frame(modelcdr1$Status)
roccdr1220$predictor = preds
roccdr1220$Classifier = "220 genes Symptomatic, without APOE"
names(roccdr1220) <- c("Status", "Predictor", "Classifier")
fit90cdr1 <- glm( Status ~  model90genes, data = modelcdr1,family=binomial())
preds = predict(fit90cdr1, type="response")
roccdr190 = as.data.frame(modelcdr1$Status)
roccdr190$predictor = preds
roccdr190$Classifier = "90 genes Symptomatic, without APOE"
names(roccdr190) <- c("Status", "Predictor", "Classifier")
fit40cdr1 <- glm( Status ~  model40genes, data = modelcdr1,family=binomial())
preds = predict(fit40cdr1, type="response")
roccdr140 = as.data.frame(modelcdr1$Status)
roccdr140$predictor = preds
roccdr140$Classifier = "40 genes Symptomatic, without APOE"
names(roccdr140) <- c("Status", "Predictor", "Classifier")

# Fit the new datasets to the models, with APOE
#Presymptomatic
fit220contr <- glm( Status ~  model220genes+apoe2+apoe4, data = modelcontr,family=binomial())
preds = predict(fit220contr, type="response")
roccontr220apoe = as.data.frame(modelcontr$Status)
roccontr220apoe$predictor = preds
roccontr220apoe$Classifier = "220 genes Presymptomatic, with APOE"
names(roccontr220apoe) <- c("Status", "Predictor", "Classifier")
fit90contr <- glm( Status ~  model90genes+apoe2+apoe4, data = modelcontr,family=binomial())
preds = predict(fit90contr, type="response")
roccontr90apoe = as.data.frame(modelcontr$Status)
roccontr90apoe$predictor = preds
roccontr90apoe$Classifier = "90 genes Presymptomatic, with APOE"
names(roccontr90apoe) <- c("Status", "Predictor", "Classifier")
fit40contr <- glm( Status ~  model40genes+apoe2+apoe4, data = modelcontr,family=binomial())
preds = predict(fit40contr, type="response")
roccontr40apoe = as.data.frame(modelcontr$Status)
roccontr40apoe$predictor = preds
roccontr40apoe$Classifier = "40 genes Presymptomatic, with APOE"
names(roccontr40apoe) <- c("Status", "Predictor", "Classifier")
#Early Symptomatic
fit220cdr05 <- glm( Status ~  model220genes+apoe2+apoe4, data = modelcdr05,family=binomial())
preds = predict(fit220cdr05, type="response")
roccdr05220apoe = as.data.frame(modelcdr05$Status)
roccdr05220apoe$predictor = preds
roccdr05220apoe$Classifier = "220 genes Early Symptomatic, with APOE"
names(roccdr05220apoe) <- c("Status", "Predictor", "Classifier")
fit90cdr05 <- glm( Status ~  model90genes+apoe2+apoe4, data = modelcdr05,family=binomial())
preds = predict(fit90cdr05, type="response")
roccdr0590apoe = as.data.frame(modelcdr05$Status)
roccdr0590apoe$predictor = preds
roccdr0590apoe$Classifier = "90 genes Early Symptomatic, with APOE"
names(roccdr0590apoe) <- c("Status", "Predictor", "Classifier")
fit40cdr05 <- glm( Status ~  model40genes+apoe2+apoe4, data = modelcdr05,family=binomial())
preds = predict(fit40cdr05, type="response")
roccdr0540apoe = as.data.frame(modelcdr05$Status)
roccdr0540apoe$predictor = preds
roccdr0540apoe$Classifier = "40 genes Early Symptomatic, with APOE"
names(roccdr0540apoe) <- c("Status", "Predictor", "Classifier")
#Symptomatic
fit220cdr1 <- glm( Status ~  model220genes+apoe2+apoe4, data = modelcdr1,family=binomial())
preds = predict(fit220cdr1, type="response")
roccdr1220apoe = as.data.frame(modelcdr1$Status)
roccdr1220apoe$predictor = preds
roccdr1220apoe$Classifier = "220 genes Symptomatic, with APOE"
names(roccdr1220apoe) <- c("Status", "Predictor", "Classifier")
fit90cdr1 <- glm( Status ~  model90genes+apoe2+apoe4, data = modelcdr1,family=binomial())
preds = predict(fit90cdr1, type="response")
roccdr190apoe = as.data.frame(modelcdr1$Status)
roccdr190apoe$predictor = preds
roccdr190apoe$Classifier = "90 genes Symptomatic, with APOE"
names(roccdr190apoe) <- c("Status", "Predictor", "Classifier")
fit40cdr1 <- glm( Status ~  model40genes+apoe2+apoe4, data = modelcdr1,family=binomial())
preds = predict(fit40cdr1, type="response")
roccdr140apoe = as.data.frame(modelcdr1$Status)
roccdr140apoe$predictor = preds
roccdr140apoe$Classifier = "40 genes Symptomatic, with APOE"
names(roccdr140apoe) <- c("Status", "Predictor", "Classifier")

#Get accuracy, specificity, positive and negative prediction value
#information foreach model and dataset
acuspec<-as.data.frame(t(caret::confusionMatrix(data = factor(as.numeric(modelcontr$ model220genes>0.5)), reference = factor(modelcontr$Status))$byClass))
acuspec<-rbind(acuspec,as.data.frame(t(caret::confusionMatrix(data = factor(as.numeric(modelcontr$ model90genes>0.5)), reference = factor(modelcontr$Status))$byClass)))
acuspec<-rbind(acuspec,as.data.frame(t(caret::confusionMatrix(data = factor(as.numeric(modelcontr$ model40genes>0.5)), reference = factor(modelcontr$Status))$byClass)))
acuspec<-rbind(acuspec,as.data.frame(t(caret::confusionMatrix(data = factor(as.numeric(modelcdr05$ model220genes>0.5)), reference = factor(modelcdr05$Status))$byClass)))
acuspec<-rbind(acuspec,as.data.frame(t(caret::confusionMatrix(data = factor(as.numeric(modelcdr05$ model90genes>0.5)), reference = factor(modelcdr05$Status))$byClass)))
acuspec<-rbind(acuspec,as.data.frame(t(caret::confusionMatrix(data = factor(as.numeric(modelcdr05$ model40genes>0.5)), reference = factor(modelcdr05$Status))$byClass)))
acuspec<-rbind(acuspec,as.data.frame(t(caret::confusionMatrix(data = factor(as.numeric(modelcdr1$ model220genes>0.5)), reference = factor(modelcdr1$Status))$byClass)))
acuspec<-rbind(acuspec,as.data.frame(t(caret::confusionMatrix(data = factor(as.numeric(modelcdr1$ model90genes>0.5)), reference = factor(modelcdr1$Status))$byClass)))
acuspec<-rbind(acuspec,as.data.frame(t(caret::confusionMatrix(data = factor(as.numeric(modelcdr1$ model40genes>0.5)), reference = factor(modelcdr1$Status))$byClass)))
acuspec<-rbind(acuspec,as.data.frame(t(caret::confusionMatrix(data = factor(as.numeric(roccontr220apoe$Predictor>0.5)), reference = factor(roccontr220apoe$Status))$byClass)))
acuspec<-rbind(acuspec,as.data.frame(t(caret::confusionMatrix(data = factor(as.numeric(roccontr90apoe$Predictor>0.5)), reference = factor(roccontr90apoe$Status))$byClass)))
acuspec<-rbind(acuspec,as.data.frame(t(caret::confusionMatrix(data = factor(as.numeric(roccontr40apoe$Predictor>0.5)), reference = factor(roccontr40apoe$Status))$byClass)))
acuspec<-rbind(acuspec,as.data.frame(t(caret::confusionMatrix(data = factor(as.numeric(roccdr05220apoe$Predictor>0.5)), reference = factor(roccdr05220apoe$Status))$byClass)))
acuspec<-rbind(acuspec,as.data.frame(t(caret::confusionMatrix(data = factor(as.numeric(roccdr0590apoe$Predictor>0.5)), reference = factor(roccdr0590apoe$Status))$byClass)))
acuspec<-rbind(acuspec,as.data.frame(t(caret::confusionMatrix(data = factor(as.numeric(roccdr0540apoe$Predictor>0.5)), reference = factor(roccdr0540apoe$Status))$byClass)))
acuspec<-rbind(acuspec,as.data.frame(t(caret::confusionMatrix(data = factor(as.numeric(roccdr1220apoe$Predictor>0.5)), reference = factor(roccdr1220apoe$Status))$byClass)))
acuspec<-rbind(acuspec,as.data.frame(t(caret::confusionMatrix(data = factor(as.numeric(roccdr190apoe$Predictor>0.5)), reference = factor(roccdr190apoe$Status))$byClass)))
acuspec<-rbind(acuspec,as.data.frame(t(caret::confusionMatrix(data = factor(as.numeric(roccdr140apoe$Predictor>0.5)), reference = factor(roccdr140apoe$Status))$byClass)))
acuspec$Dataset<-rep(c(rep("Presymptomatic",3),rep("Early Symptomatic",3),rep("Symptomatic",3)),2)
acuspec$APOE<-c(rep(0,9),rep(1,9))
acuspec$Gene_model<-rep(c(220,90,40),6)
##AUC whisker plot
t<-ci.auc(as.numeric(roccdr05220$Status), as.numeric(roccdr05220$Predictor),conf.level = 0.9)
aucdf<-data.frame(low = round(t[1],2), AUC = round(t[2],2), high = round(t[3],2), Classifier = "Early Symptomatic-220-noAPOE",Genes = "model with 220 transcripts", APOE = "without APOE")
t<-ci.auc(as.numeric(roccdr0590$Status), as.numeric(roccdr0590$Predictor),conf.level = 0.9)
t<-data.frame(low = round(t[1],2), AUC = round(t[2],2), high = round(t[3],2), Classifier = "Early Symptomatic-90-noAPOE",Genes = "model with 90 transcripts", APOE = "without APOE")
aucdf<-rbind(aucdf, t)
t<-ci.auc(as.numeric(roccdr0540$Status), as.numeric(roccdr0540$Predictor),conf.level = 0.9)
t<-data.frame(low = round(t[1],2), AUC = round(t[2],2), high = round(t[3],2), Classifier = "Early Symptomatic-40-noAPOE",Genes = "model with 40 transcripts", APOE = "without APOE")
aucdf<-rbind(aucdf, t)
t<-ci.auc(as.numeric(roccdr05220apoe$Status), as.numeric(roccdr05220apoe$Predictor),conf.level = 0.9)
t<-data.frame(low = round(t[1],2), AUC = round(t[2],2), high = round(t[3],2), Classifier = "Early Symptomatic-220-APOE",Genes = "model with 220 transcripts", APOE = "with APOE")
aucdf<-rbind(aucdf, t)
t<-ci.auc(as.numeric(roccdr0590apoe$Status), as.numeric(roccdr0590apoe$Predictor),conf.level = 0.9)
t<-data.frame(low = round(t[1],2), AUC = round(t[2],2), high = round(t[3],2), Classifier = "Early Symptomatic-90-APOE",Genes = "model with 90 transcripts", APOE = "with APOE")
aucdf<-rbind(aucdf, t)
t<-ci.auc(as.numeric(roccdr0540apoe$Status), as.numeric(roccdr0540apoe$Predictor),conf.level = 0.9)
t<-data.frame(low = round(t[1],2), AUC = round(t[2],2), high = round(t[3],2), Classifier = "Early Symptomatic-40-APOE",Genes = "model with 40 transcripts", APOE = "with APOE")
aucdf<-rbind(aucdf, t)
t<-ci.auc(as.numeric(roccdr1220$Status), as.numeric(roccdr1220$Predictor),conf.level = 0.9)
t<-data.frame(low = round(t[1],2), AUC = round(t[2],2), high = round(t[3],2), Classifier = "Symptomatic-220-noAPOE", Genes = "model with 220 transcripts", APOE = "without APOE")
aucdf<-rbind(aucdf, t)
t<-ci.auc(as.numeric(roccdr190$Status), as.numeric(roccdr190$Predictor),conf.level = 0.9)
t<-data.frame(low = round(t[1],2), AUC = round(t[2],2), high = round(t[3],2), Classifier = "Symptomatic-90-noAPOE",Genes = "model with 90 transcripts", APOE = "without APOE")
aucdf<-rbind(aucdf, t)
t<-ci.auc(as.numeric(roccdr140$Status), as.numeric(roccdr140$Predictor),conf.level = 0.9)
t<-data.frame(low = round(t[1],2), AUC = round(t[2],2), high = round(t[3],2), Classifier = "Symptomatic-40-noAPOE",Genes = "model with 40 transcripts", APOE = "without APOE")
aucdf<-rbind(aucdf, t)
t<-ci.auc(as.numeric(roccdr1220apoe$Status), as.numeric(roccdr1220apoe$Predictor),conf.level = 0.9)
t<-data.frame(low = round(t[1],2), AUC = round(t[2],2), high = round(t[3],2), Classifier = "Symptomatic-220-APOE",Genes = "model with 220 transcripts", APOE = "with APOE")
aucdf<-rbind(aucdf, t)
t<-ci.auc(as.numeric(roccdr190apoe$Status), as.numeric(roccdr190apoe$Predictor),conf.level = 0.9)
t<-data.frame(low = round(t[1],2), AUC = round(t[2],2), high = round(t[3],2), Classifier = "Symptomatic-90-APOE",Genes = "model with 90 transcripts", APOE = "with APOE")
aucdf<-rbind(aucdf, t)
t<-ci.auc(as.numeric(roccdr140apoe$Status), as.numeric(roccdr140apoe$Predictor),conf.level = 0.9)
t<-data.frame(low = round(t[1],2), AUC = round(t[2],2), high = round(t[3],2), Classifier = "Symptomatic-40-APOE",Genes = "model with 40 transcripts", APOE = "with APOE")
aucdf<-rbind(aucdf, t)

aucdf$Classifier<-factor(aucdf$Classifier, levels=c("Early Symptomatic-40-noAPOE","Early Symptomatic-40-APOE",
                                                   "Early Symptomatic-90-noAPOE","Early Symptomatic-90-APOE","Early Symptomatic-220-noAPOE","Early Symptomatic-220-APOE",
                                                   "Symptomatic-40-noAPOE","Symptomatic-40-APOE","Symptomatic-90-noAPOE","Symptomatic-90-APOE","Symptomatic-220-noAPOE","Symptomatic-220-APOE"))

whiskerPlot <- ggplot(aucdf, aes(x=AUC, y=Classifier,color = Genes))+
  geom_point(aes(shape=APOE), size =4)+
  scale_colour_manual(values=myPalette)+
  theme_bw()+
  scale_x_continuous("ROC AUC", limits = c(0.4,1), n.breaks = 8, expand = c(0,0.005)) +
  geom_errorbarh(aes(xmin=low,xmax=high),height=0)+
  scale_y_discrete(limits=rev)+
  geom_vline(xintercept = 0.5, linetype = "dashed")+
  theme(legend.position = "bottom")+
  guides(color = guide_legend(reverse = TRUE, override.aes = list(shape = 15)), shape = guide_legend(reverse = TRUE))
