#Load the necessary libraries
library(dplyr)
library(readr)
library(ggplot2)
library(scales)
library(plotROC)
library(glmnet)
library(pROC)
library(data.table)
library(RColorBrewer)
myPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999", "#000000", "#6C2EB8")

## load phenotype file containing AT_Status infor
library("readxl")
pheno <- as.data.frame(read_excel("CSF-Values_Paper_ATN.xlsx"))
dim(pheno) # 72 21
table(pheno$ATStatus...17)
###########################
##### 220-genes model #####
###########################
genestodoml220 = fread("220genes.csv")
genestodoml220 = genestodoml220$ENSG
genestodoml = genestodoml220
# Load Matrix
matrix_phase = readRDS("discovery_data")
# Change training df colnames to make them consistent with testing df colnames
colnames(matrix_phase)[67] <- "ID" # from "SampleID" to "ID"
colnames(matrix_phase)[70] <- "Age" # from "age_at_draw" to "Age"
colnames(matrix_phase)[74] <- "APOE" # from "apoe" to "APOE"
matrix_phase$Gender <- ifelse(matrix_phase$sex == "Male", 1, 2)

train = matrix_phase
APOEtrain<-train$APOE
head(APOEtrain)
# Select features
genestodoml = c(genestodoml, "Status")
commontrain = colnames(train) %in% genestodoml
train <- as.data.frame(train)
train <- train[,colnames(train) %in% genestodoml]
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
predicciones_train_220 <- predict(modelo, newx = x_train, type = "response", standardize = F)
genestodoml = genestodoml220
genestodoml = c(genestodoml, "Status")
# Load test Matrix
test = readRDS("replication_data")
#Capture only the columns shared between the testing and training data 
testdf <- test[,colnames(test) %in% colnames(matrix_phase)]
matrix_phase <- as.data.frame(matrix_phase)
traindf <- matrix_phase[,colnames(matrix_phase) %in% colnames(test)]
# Make sure columns in both dataset are in the same order
traindf <- traindf %>% select(colnames(testdf))
identical(colnames(testdf), colnames(traindf)) # TRUE
#Combine the two datasets
test <- rbind(testdf, traindf)
table(test$ID %in% pheno$ID)

pheno_test <- pheno[pheno$ID %in% test$ID,]
dim(pheno_test) # 72 21
table(pheno_test$ATStatus...17)

pheno_test_AT <- pheno_test[,c("ID","ATStatus...17")]
colnames(pheno_test_AT)[2] <- "ATStatus"
pheno_test_AT <- pheno_test_AT[pheno_test_AT$ATStatus %in% c("AT-Neg","AT-Pos"),]
pheno_test_AT$AT_Status <- ifelse(pheno_test_AT$ATStatus == "AT-Neg", 0, 1)
dim(pheno_test_AT) # 26  3
table(pheno_test_AT$ATStatus)

table(pheno_test_AT$AT_Status)
pheno_test_AT$ATStatus <- NULL

test_AT <- test[test$ID %in% pheno_test_AT$ID,]
test_AT <- inner_join(test_AT, pheno_test_AT, by="ID")
# get model covariates
APOEtest <- test_AT$APOE
tocorrelationtest <- test_AT[,c("ID", "AT_Status", "Gender", "Age", "APOE")]
# Select features
commontest = colnames(test_AT) %in% genestodoml
test_AT_subset = test_AT[,commontest]
# Select only genes
x_test <- test_AT_subset[,-1]
# Compute Z score
x_test <- as.data.frame(lapply(x_test, function(x) (x - mean(x))/sd(x) ))
x_test <- as.matrix(x_test)
# Code status
y_test = as.numeric(test_AT$AT_Status) # binarized AT_Status (0,1)
y_test_220 = y_test
# Testing preditions
predicciones_test <- predict(modelo, newx = x_test, type = "response" , standarize = F, intercept = F)
predicciones_test_220 = predicciones_test
# Model220 roc comes from above code (220 genes model)
model220rocte = cbind(as.data.frame(y_test_220),APOEtest)
model220rocte$predictor = predicciones_test_220
names(model220rocte) = c("status", "APOE","model220")

model220ad = model220rocte

###########################
##### 90-genes model #####
###########################
genestodoml90 = fread("90genes.csv")
genestodoml90 = genestodoml90[-1,1]
genestodoml90 = genestodoml90$V1
genestodoml = genestodoml90
# Load Matrix
matrix_phase = readRDS("discovery_data")
# Change training df colnames to make them consistent with testing df colnames
colnames(matrix_phase)[67] <- "ID" # from "SampleID" to "ID"
colnames(matrix_phase)[70] <- "Age" # from "age_at_draw" to "Age"
colnames(matrix_phase)[74] <- "APOE" # from "apoe" to "APOE"
matrix_phase$Gender <- ifelse(matrix_phase$sex == "Male", 1, 2)

train = matrix_phase
APOEtrain<-train$APOE
head(APOEtrain)
# Select features
genestodoml = c(genestodoml, "Status")
commontrain = colnames(train) %in% genestodoml
train <- as.data.frame(train)
train <- train[,colnames(train) %in% genestodoml]
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
predicciones_train_90 <- predict(modelo, newx = x_train, type = "response", standardize = F)
genestodoml = genestodoml90
genestodoml = c(genestodoml, "Status")
# Load test Matrix
test = readRDS("replication_data")
#Capture only the columns shared between the testing and training data 
testdf <- test[,colnames(test) %in% colnames(matrix_phase)]
matrix_phase <- as.data.frame(matrix_phase)
traindf <- matrix_phase[,colnames(matrix_phase) %in% colnames(test)]
# Make sure columns in both dataset are in the same order
traindf <- traindf %>% select(colnames(testdf))
identical(colnames(testdf), colnames(traindf)) # TRUE
#Combine the two datasets
test <- rbind(testdf, traindf)
table(test$ID %in% pheno$ID)

pheno_test <- pheno[pheno$ID %in% test$ID,]
table(pheno_test$ATStatus...17)
pheno_test_AT <- pheno_test[,c("ID","ATStatus...17")]
colnames(pheno_test_AT)[2] <- "ATStatus"
pheno_test_AT <- pheno_test_AT[pheno_test_AT$ATStatus %in% c("AT-Neg","AT-Pos"),]
pheno_test_AT$AT_Status <- ifelse(pheno_test_AT$ATStatus == "AT-Neg", 0, 1)
table(pheno_test_AT$ATStatus)
table(pheno_test_AT$AT_Status)
pheno_test_AT$ATStatus <- NULL

test_AT <- test[test$ID %in% pheno_test_AT$ID,]
test_AT <- inner_join(test_AT, pheno_test_AT, by="ID")
# get model covariates
APOEtest <- test_AT$APOE
tocorrelationtest <- test_AT[,c("ID", "AT_Status", "Gender", "Age", "APOE")]
# Select features
commontest = colnames(test_AT) %in% genestodoml
test_AT_subset = test_AT[,commontest]
# Select only genes
x_test <- test_AT_subset[,-1]
# Compute Z score
x_test <- as.data.frame(lapply(x_test, function(x) (x - mean(x))/sd(x) ))
x_test <- as.matrix(x_test)
# Code status
y_test = as.numeric(test_AT$AT_Status) # binarized AT_Status (0,1)
y_test_90 = y_test
# Testing preditions
predicciones_test <- predict(modelo, newx = x_test, type = "response" , standarize = F, intercept = F)
predicciones_test_90 = predicciones_test
# Model90 roc comes from above code (90 genes model)
model90rocte = cbind(as.data.frame(y_test_90),APOEtest)
model90rocte$predictor = predicciones_test_90
names(model90rocte) = c("status", "APOE","model90")
model90ad = model90rocte
###########################
##### 40-genes model #####
###########################
genestodoml40 = fread("40genes.csv")
genestodoml40 = genestodoml40[-1,1]
genestodoml40 = genestodoml40$V1
genestodoml = genestodoml40
# Load Matrix
matrix_phase = readRDS("discovery_data")
# Change training df colnames to make them consistent with testing df colnames
colnames(matrix_phase)[67] <- "ID" # from "SampleID" to "ID"
colnames(matrix_phase)[70] <- "Age" # from "age_at_draw" to "Age"
colnames(matrix_phase)[74] <- "APOE" # from "apoe" to "APOE"
matrix_phase$Gender <- ifelse(matrix_phase$sex == "Male", 1, 2)

train = matrix_phase
APOEtrain<-train$APOE
head(APOEtrain)
# Select features
genestodoml = c(genestodoml, "Status")
commontrain = colnames(train) %in% genestodoml
train <- as.data.frame(train)
train <- train[,colnames(train) %in% genestodoml]
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
predicciones_train_40 <- predict(modelo, newx = x_train, type = "response", standardize = F)
genestodoml = genestodoml40
genestodoml = c(genestodoml, "Status")
# Load test Matrix
test = readRDS("replication_data")
#Capture only the columns shared between the testing and training data 
testdf <- test[,colnames(test) %in% colnames(matrix_phase)]
matrix_phase <- as.data.frame(matrix_phase)
traindf <- matrix_phase[,colnames(matrix_phase) %in% colnames(test)]
# Make sure columns in both dataset are in the same order
traindf <- traindf %>% select(colnames(testdf))
identical(colnames(testdf), colnames(traindf)) # TRUE
#Combine the two datasets
test <- rbind(testdf, traindf)
table(test$ID %in% pheno$ID)

pheno_test <- pheno[pheno$ID %in% test$ID,]
table(pheno_test$ATStatus...17)

pheno_test_AT <- pheno_test[,c("ID","ATStatus...17")]
colnames(pheno_test_AT)[2] <- "ATStatus"
pheno_test_AT <- pheno_test_AT[pheno_test_AT$ATStatus %in% c("AT-Neg","AT-Pos"),]
pheno_test_AT$AT_Status <- ifelse(pheno_test_AT$ATStatus == "AT-Neg", 0, 1)
table(pheno_test_AT$ATStatus)

table(pheno_test_AT$AT_Status)

pheno_test_AT$ATStatus <- NULL

test_AT <- test[test$ID %in% pheno_test_AT$ID,]
test_AT <- inner_join(test_AT, pheno_test_AT, by="ID")

# get model covariates
APOEtest <- test_AT$APOE
tocorrelationtest <- test_AT[,c("ID", "AT_Status", "Gender", "Age", "APOE")]

# Select features
commontest = colnames(test_AT) %in% genestodoml
test_AT_subset = test_AT[,commontest]
# Select only genes
x_test <- test_AT_subset[,-1]
# Compute Z score
x_test <- as.data.frame(lapply(x_test, function(x) (x - mean(x))/sd(x) ))
x_test <- as.matrix(x_test)
# Code status
y_test = as.numeric(test_AT$AT_Status) # binarized AT_Status (0,1)
y_test_40 = y_test
# Testing preditions
predicciones_test <- predict(modelo, newx = x_test, type = "response" , standarize = F, intercept = F)
predicciones_test_40 = predicciones_test

# Model40 roc comes from above code (40 genes model)
model40rocte = cbind(as.data.frame(y_test_40),APOEtest)
model40rocte$predictor = predicciones_test_40
names(model40rocte) = c("status", "APOE","model40")
model40ad = model40rocte

## Combine data for all 3 models into a single data.frame
modelad <- cbind(model220ad, as.data.frame(model90ad$model90))
modelad <- cbind(modelad, as.data.frame(model40ad$model40))
modelad$s0 <- as.numeric(unlist(modelad$s0))
names(modelad)<-c("Status", "APOE", "model220genes", "model90genes", "model40genes")
# Add columns to code for number of apoe2 and apoe4 allelles
# for each subject in each data frame
modelad$apoe4 = 0
modelad$apoe4 [modelad$APOE == 44] <- 2
modelad$apoe4 [modelad$APOE== 43] <- 1
modelad$apoe4 [modelad$APOE== 34 ] <- 1
modelad$apoe4 [modelad$APOE== 24] <- 1
modelad$apoe4 [modelad$APOE== 42] <- 1
modelad$apoe2 = 0
modelad$apoe2 [modelad$APOE == 22] <- 2
modelad$apoe2 [modelad$APOE== 23] <- 1
modelad$apoe2 [modelad$APOE== 32 ] <- 1
modelad$apoe2 [modelad$APOE== 24] <- 1
# Fit the dataset to the models, without APOE
fit220ad <- glm( Status ~  model220genes, data = modelad,family=binomial())
preds = predict(fit220ad, type="response")
rocad220 = as.data.frame(modelad$Status)
rocad220$predictor = preds
rocad220$Classifier = "220 genes AD, without APOE (AUC = 0.69)"
names(rocad220) <- c("Status", "Predictor", "Classifier")
fit90ad <- glm( Status ~  model90genes, data = modelad,family=binomial())
preds = predict(fit90ad, type="response")
rocad90 = as.data.frame(modelad$Status)
rocad90$predictor = preds
rocad90$Classifier = "90 genes AD, without APOE (AUC = 0.65)"
names(rocad90) <- c("Status", "Predictor", "Classifier")
fit40ad <- glm( Status ~  model40genes, data = modelad,family=binomial())
preds = predict(fit40ad, type="response")
rocad40 = as.data.frame(modelad$Status)
rocad40$predictor = preds
rocad40$Classifier = "40 genes AD, without APOE (AUC = 0.55)"
names(rocad40) <- c("Status", "Predictor", "Classifier")
# Fit the dataset to the models, with APOE
fit220ad <- glm( Status ~  model220genes+apoe2+apoe4, data = modelad,family=binomial())
preds = predict(fit220ad, type="response")
rocad220apoe = as.data.frame(modelad$Status)
rocad220apoe$predictor = preds
rocad220apoe$Classifier = "220 genes AD, with APOE (AUC = 0.72)"
names(rocad220apoe) <- c("Status", "Predictor", "Classifier")
fit90ad <- glm( Status ~  model90genes+apoe2+apoe4, data = modelad,family=binomial())
preds = predict(fit90ad, type="response")
rocad90apoe = as.data.frame(modelad$Status)
rocad90apoe$predictor = preds
rocad90apoe$Classifier = "90 genes AD, with APOE (AUC = 0.69)"
names(rocad90apoe) <- c("Status", "Predictor", "Classifier")
fit40ad <- glm( Status ~  model40genes+apoe2+apoe4, data = modelad,family=binomial())
preds = predict(fit40ad, type="response")
rocad40apoe = as.data.frame(modelad$Status)
rocad40apoe$predictor = preds
rocad40apoe$Classifier = "40 genes AD, with APOE (AUC = 0.65)"
names(rocad40apoe) <- c("Status", "Predictor", "Classifier")

##AUC whisker plot
t<-ci.auc(as.numeric(rocad220$Status), as.numeric(rocad220$Predictor),conf.level = 0.9)
aucdf<-data.frame(low = round(t[1],2), AUC = round(t[2],2), high = round(t[3],2), Classifier = "AT-220-noAPOE", Genes = "model with 220 genes", APOE = "without APOE")
t<-ci.auc(as.numeric(rocad90$Status), as.numeric(rocad90$Predictor),conf.level = 0.9)
t<-data.frame(low = round(t[1],2), AUC = round(t[2],2), high = round(t[3],2), Classifier = "AT-90-noAPOE",Genes = "model with 90 genes", APOE = "without APOE")
aucdf<-rbind(aucdf, t)
t<-ci.auc(as.numeric(rocad40$Status), as.numeric(rocad40$Predictor),conf.level = 0.9)
t<-data.frame(low = round(t[1],2), AUC = round(t[2],2), high = round(t[3],2), Classifier = "AT-40-noAPOE",Genes = "model with 40 genes", APOE = "without APOE")
aucdf<-rbind(aucdf, t)
t<-ci.auc(as.numeric(rocad220apoe$Status), as.numeric(rocad220apoe$Predictor),conf.level = 0.9)
t<-data.frame(low = round(t[1],2), AUC = round(t[2],2), high = round(t[3],2), Classifier = "AT-220-APOE",Genes = "model with 220 genes", APOE = "with APOE")
aucdf<-rbind(aucdf, t)
t<-ci.auc(as.numeric(rocad90apoe$Status), as.numeric(rocad90apoe$Predictor),conf.level = 0.9)
t<-data.frame(low = round(t[1],2), AUC = round(t[2],2), high = round(t[3],2), Classifier = "AT-90-APOE",Genes = "model with 90 genes", APOE = "with APOE")
aucdf<-rbind(aucdf, t)
t<-ci.auc(as.numeric(rocad40apoe$Status), as.numeric(rocad40apoe$Predictor),conf.level = 0.9)
t<-data.frame(low = round(t[1],2), AUC = round(t[2],2), high = round(t[3],2), Classifier = "AT-40-APOE",Genes = "model with 40 genes", APOE = "with APOE")
aucdf<-rbind(aucdf, t)
aucdf$Classifier<-factor(aucdf$Classifier, levels=c("AT-40-noAPOE", "AT-40-APOE","AT-90-noAPOE", "AT-90-APOE", 
                                                    "AT-220-noAPOE", "AT-220-APOE")) 

whiskerPlot <- ggplot(aucdf, aes(x=AUC, y=Classifier,color = Genes))+
  geom_point(aes(shape=APOE), size =4)+
  scale_colour_manual(values=myPalette)+
  theme_bw()+
  theme(legend.position = "bottom",axis.text.x = element_text(angle = 90))+
  scale_y_discrete(limits=rev, breaks=NULL)+
  scale_x_continuous("ROC AUC", limits = c(0.4,1), breaks = c(0.4, 0.5, aucdf$AUC,1), expand = c(0,0)) +
  geom_errorbarh(aes(xmin=low,xmax=high),height=0)+
  geom_vline(xintercept = 0.5, linetype = "dashed")+
  guides(color = guide_legend(reverse = TRUE, override.aes = list(shape = 15)), shape = guide_legend(reverse = TRUE))

## load phenotype file containing AT_Status infor
library("readxl")
pheno <- as.data.frame(read_excel("CSF-Values_Paper_ATN.xlsx"))
table(pheno$Astatus)

###########################
##### 220-genes model #####
###########################

genestodoml220 = fread("220genes.csv")
genestodoml220 = genestodoml220$ENSG
genestodoml = genestodoml220
# Load Matrix
matrix_phase = readRDS("discovery_data")
# Change training df colnames to make them consistent with testing df colnames
colnames(matrix_phase)[67] <- "ID" # from "SampleID" to "ID"
colnames(matrix_phase)[70] <- "Age" # from "age_at_draw" to "Age"
colnames(matrix_phase)[74] <- "APOE" # from "apoe" to "APOE"
matrix_phase$Gender <- ifelse(matrix_phase$sex == "Male", 1, 2)

train = matrix_phase
APOEtrain<-train$APOE
head(APOEtrain) 

# Select features
genestodoml = c(genestodoml, "Status")
commontrain = colnames(train) %in% genestodoml
train <- as.data.frame(train)
train <- train[,colnames(train) %in% genestodoml]

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
predicciones_train_220 <- predict(modelo, newx = x_train, type = "response", standardize = F)
genestodoml = genestodoml220
genestodoml = c(genestodoml, "Status")

# Load test Matrix
test = readRDS("replication_data")

#Capture only the columns shared between the testing and training data 
testdf <- test[,colnames(test) %in% colnames(matrix_phase)]
matrix_phase <- as.data.frame(matrix_phase)
traindf <- matrix_phase[,colnames(matrix_phase) %in% colnames(test)]
# Make sure columns in both dataset are in the same order
traindf <- traindf %>% select(colnames(testdf))
identical(colnames(testdf), colnames(traindf)) # TRUE
#Combine the two datasets
test <- rbind(testdf, traindf)
table(test$ID %in% pheno$ID)

pheno_test <- pheno[pheno$ID %in% test$ID,]
table(pheno_test$Astatus)

pheno_test_AT <- pheno_test[,c("ID","Astatus")]
pheno_test_AT <- pheno_test_AT[pheno_test_AT$Astatus %in% c("Negative","Positive"),]
pheno_test_AT$A_Status <- ifelse(pheno_test_AT$Astatus == "Negative", 0, 1)
table(pheno_test_AT$A_Status)

table(pheno_test_AT$Astatus)
# Negative Positive
pheno_test_AT$Astatus <- NULL
test_AT <- test[test$ID %in% pheno_test_AT$ID,]
test_AT <- inner_join(test_AT, pheno_test_AT, by="ID")
# get model covariates
APOEtest <- test_AT$APOE
tocorrelationtest <- test_AT[,c("ID", "A_Status", "Gender", "Age", "APOE")]
head(tocorrelationtest, 2)

# Select features
commontest = colnames(test_AT) %in% genestodoml
test_AT_subset = test_AT[,commontest]
dim(test_AT_subset) # 72 221
# Select only genes
x_test <- test_AT_subset[,-1]
# Compute Z score
x_test <- as.data.frame(lapply(x_test, function(x) (x - mean(x))/sd(x) ))
x_test <- as.matrix(x_test)
# Code status
y_test = as.numeric(test_AT$A_Status) # binarized A_Status (0,1)
y_test_220 = y_test
# Testing preditions
predicciones_test <- predict(modelo, newx = x_test, type = "response" , standarize = F, intercept = F)
predicciones_test_220 = predicciones_test

# Model220 roc comes from above code (220 genes model)
model220rocte = cbind(as.data.frame(y_test_220),APOEtest)
model220rocte$predictor = predicciones_test_220
names(model220rocte) = c("status", "APOE","model220")
model220ad = model220rocte

###########################
##### 90-genes model #####
###########################
genestodoml90 = fread("90genes.csv")
genestodoml90 = genestodoml90[-1,1]
genestodoml90 = genestodoml90$V1
genestodoml = genestodoml90
# Load Matrix
matrix_phase = readRDS("discovery_data")
# Change training df colnames to make them consistent with testing df colnames
colnames(matrix_phase)[67] <- "ID" # from "SampleID" to "ID"
colnames(matrix_phase)[70] <- "Age" # from "age_at_draw" to "Age"
colnames(matrix_phase)[74] <- "APOE" # from "apoe" to "APOE"
matrix_phase$Gender <- ifelse(matrix_phase$sex == "Male", 1, 2)

train = matrix_phase
APOEtrain<-train$APOE
head(APOEtrain)

# Select features
genestodoml = c(genestodoml, "Status")
commontrain = colnames(train) %in% genestodoml
train <- as.data.frame(train)
train <- train[,colnames(train) %in% genestodoml]

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
predicciones_train_90 <- predict(modelo, newx = x_train, type = "response", standardize = F)
genestodoml = genestodoml90
genestodoml = c(genestodoml, "Status")
# Load test Matrix
test = readRDS("replication_data")
# Load test matrix prepared for A-neg A-pos individuals
load("testData_for_Aposneg.RData")
# get model covariates
APOEtest <- test_AT$APOE
tocorrelationtest <- test_AT[,c("ID", "A_Status", "Gender", "Age", "APOE")]
# Select features
commontest = colnames(test_AT) %in% genestodoml
test_AT_subset = test_AT[,commontest]
# Select only genes
x_test <- test_AT_subset[,-1]
# Compute Z score
x_test <- as.data.frame(lapply(x_test, function(x) (x - mean(x))/sd(x) ))
x_test <- as.matrix(x_test)
# Code status
y_test = as.numeric(test_AT$A_Status) # binarized AT_Status (0,1)
y_test_90 = y_test
# Testing preditions
predicciones_test <- predict(modelo, newx = x_test, type = "response" , standarize = F, intercept = F)
predicciones_test_90 = predicciones_test
# Model90 roc comes from above code (90 genes model)
model90rocte = cbind(as.data.frame(y_test_90),APOEtest)
model90rocte$predictor = predicciones_test_90
names(model90rocte) = c("status", "APOE","model90")
model90ad = model90rocte
###########################
##### 40-genes model #####
###########################
genestodoml40 = fread("40genes.csv")
genestodoml40 = genestodoml40[-1,1]
genestodoml40 = genestodoml40$V1
genestodoml = genestodoml40
length(genestodoml) # 40
# Load Matrix
matrix_phase = readRDS("discovery_data")
# Change training df colnames to make them consistent with testing df colnames
colnames(matrix_phase)[67] <- "ID" # from "SampleID" to "ID"
colnames(matrix_phase)[70] <- "Age" # from "age_at_draw" to "Age"
colnames(matrix_phase)[74] <- "APOE" # from "apoe" to "APOE"
matrix_phase$Gender <- ifelse(matrix_phase$sex == "Male", 1, 2)

train = matrix_phase
APOEtrain<-train$APOE
head(APOEtrain)
# Select features
genestodoml = c(genestodoml, "Status")
commontrain = colnames(train) %in% genestodoml
train <- as.data.frame(train)
train <- train[,colnames(train) %in% genestodoml]
dim(train) #  73 41

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
predicciones_train_40 <- predict(modelo, newx = x_train, type = "response", standardize = F)
genestodoml = genestodoml40
genestodoml = c(genestodoml, "Status")
# Load test Matrix
test = readRDS("replication_data")
# Load test matrix prepared for A-neg A-pos individuals
load("testData_for_Aposneg.RData")
# get model covariates
APOEtest <- test_AT$APOE
tocorrelationtest <- test_AT[,c("ID", "A_Status", "Gender", "Age", "APOE")]
# Select features
commontest = colnames(test_AT) %in% genestodoml
test_AT_subset = test_AT[,commontest]
dim(test_AT_subset) # 72 41
# Select only genes
x_test <- test_AT_subset[,-1]
# Compute Z score
x_test <- as.data.frame(lapply(x_test, function(x) (x - mean(x))/sd(x) ))
x_test <- as.matrix(x_test)
dim(x_test) #  72 40
# Code status
y_test = as.numeric(test_AT$A_Status) # binarized A_Status (0,1)
y_test_40 = y_test
# Testing preditions
predicciones_test <- predict(modelo, newx = x_test, type = "response" , standarize = F, intercept = F)
predicciones_test_40 = predicciones_test
# Model40 roc comes from above code (40 genes model)
model40rocte = cbind(as.data.frame(y_test_40),APOEtest)
model40rocte$predictor = predicciones_test_40
names(model40rocte) = c("status", "APOE","model40")
model40ad = model40rocte
## Combine data for all 3 models into a single data.frame
modelad <- cbind(model220ad, as.data.frame(model90ad$model90))
modelad <- cbind(modelad, as.data.frame(model40ad$model40))
modelad$s0 <- as.numeric(unlist(modelad$s0))
names(modelad)<-c("Status", "APOE", "model220genes", "model90genes", "model40genes")
dim(modelad) # 72   5
# Add columns to code for number of apoe2 and apoe4 allelles
# for each subject in each data frame
modelad$apoe4 = 0
modelad$apoe4 [modelad$APOE == 44] <- 2
modelad$apoe4 [modelad$APOE== 43] <- 1
modelad$apoe4 [modelad$APOE== 34 ] <- 1
modelad$apoe4 [modelad$APOE== 24] <- 1
modelad$apoe4 [modelad$APOE== 42] <- 1
modelad$apoe2 = 0
modelad$apoe2 [modelad$APOE == 22] <- 2
modelad$apoe2 [modelad$APOE== 23] <- 1
modelad$apoe2 [modelad$APOE== 32 ] <- 1
modelad$apoe2 [modelad$APOE== 24] <- 1
# Fit the dataset to the models, without APOE
fit220ad <- glm( Status ~  model220genes, data = modelad,family=binomial())
preds = predict(fit220ad, type="response")
rocad220 = as.data.frame(modelad$Status)
rocad220$predictor = preds
rocad220$Classifier = "220 genes AD, without APOE (AUC = 0.69)"
names(rocad220) <- c("Status", "Predictor", "Classifier")
fit90ad <- glm( Status ~  model90genes, data = modelad,family=binomial())
preds = predict(fit90ad, type="response")
rocad90 = as.data.frame(modelad$Status)
rocad90$predictor = preds
rocad90$Classifier = "90 genes AD, without APOE (AUC = 0.65)"
names(rocad90) <- c("Status", "Predictor", "Classifier")
fit40ad <- glm( Status ~  model40genes, data = modelad,family=binomial())
preds = predict(fit40ad, type="response")
rocad40 = as.data.frame(modelad$Status)
rocad40$predictor = preds
rocad40$Classifier = "40 genes AD, without APOE (AUC = 0.55)"
names(rocad40) <- c("Status", "Predictor", "Classifier")
# Fit the dataset to the models, with APOE
fit220ad <- glm( Status ~  model220genes+apoe2+apoe4, data = modelad,family=binomial())
preds = predict(fit220ad, type="response")
rocad220apoe = as.data.frame(modelad$Status)
rocad220apoe$predictor = preds
rocad220apoe$Classifier = "220 genes AD, with APOE (AUC = 0.72)"
names(rocad220apoe) <- c("Status", "Predictor", "Classifier")
fit90ad <- glm( Status ~  model90genes+apoe2+apoe4, data = modelad,family=binomial())
preds = predict(fit90ad, type="response")
rocad90apoe = as.data.frame(modelad$Status)
rocad90apoe$predictor = preds
rocad90apoe$Classifier = "90 genes AD, with APOE (AUC = 0.69)"
names(rocad90apoe) <- c("Status", "Predictor", "Classifier")
fit40ad <- glm( Status ~  model40genes+apoe2+apoe4, data = modelad,family=binomial())
preds = predict(fit40ad, type="response")
rocad40apoe = as.data.frame(modelad$Status)
rocad40apoe$predictor = preds
rocad40apoe$Classifier = "40 genes AD, with APOE (AUC = 0.65)"
names(rocad40apoe) <- c("Status", "Predictor", "Classifier")
##AUC whisker plot
t<-ci.auc(as.numeric(rocad220$Status), as.numeric(rocad220$Predictor),conf.level = 0.9)
aucdf<-data.frame(low = round(t[1],2), AUC = round(t[2],2), high = round(t[3],2), Classifier = "A-220-noAPOE", Genes = "model with 220 genes", APOE = "without APOE")
t<-ci.auc(as.numeric(rocad90$Status), as.numeric(rocad90$Predictor),conf.level = 0.9)
t<-data.frame(low = round(t[1],2), AUC = round(t[2],2), high = round(t[3],2), Classifier = "A-90-noAPOE",Genes = "model with 90 genes", APOE = "without APOE")
aucdf<-rbind(aucdf, t)
t<-ci.auc(as.numeric(rocad40$Status), as.numeric(rocad40$Predictor),conf.level = 0.9)
t<-data.frame(low = round(t[1],2), AUC = round(t[2],2), high = round(t[3],2), Classifier = "A-40-noAPOE",Genes = "model with 40 genes", APOE = "without APOE")
aucdf<-rbind(aucdf, t)
t<-ci.auc(as.numeric(rocad220apoe$Status), as.numeric(rocad220apoe$Predictor),conf.level = 0.9)
t<-data.frame(low = round(t[1],2), AUC = round(t[2],2), high = round(t[3],2), Classifier = "A-220-APOE",Genes = "model with 220 genes", APOE = "with APOE")
aucdf<-rbind(aucdf, t)
t<-ci.auc(as.numeric(rocad90apoe$Status), as.numeric(rocad90apoe$Predictor),conf.level = 0.9)
t<-data.frame(low = round(t[1],2), AUC = round(t[2],2), high = round(t[3],2), Classifier = "A-90-APOE",Genes = "model with 90 genes", APOE = "with APOE")
aucdf<-rbind(aucdf, t)
t<-ci.auc(as.numeric(rocad40apoe$Status), as.numeric(rocad40apoe$Predictor),conf.level = 0.9)
t<-data.frame(low = round(t[1],2), AUC = round(t[2],2), high = round(t[3],2), Classifier = "A-40-APOE",Genes = "model with 40 genes", APOE = "with APOE")
aucdf<-rbind(aucdf, t)
write.csv(aucdf, "AUC_Astatus.csv", row.names = F)
aucdf$Classifier<-factor(aucdf$Classifier, levels=c("A-40-noAPOE", "A-40-APOE","A-90-noAPOE", "A-90-APOE", 
                                                    "A-220-noAPOE", "A-220-APOE")) 

whiskerPlot <- ggplot(aucdf, aes(x=AUC, y=Classifier,color = Genes))+
  geom_point(aes(shape=APOE), size =4)+
  scale_colour_manual(values=myPalette)+
  theme_bw()+
  theme(legend.position = "bottom",axis.text.x = element_text(angle = 90))+
  scale_y_discrete(limits=rev, breaks=NULL)+
  scale_x_continuous("ROC AUC", limits = c(0.4,1), breaks = c(0.4, 0.5, aucdf$AUC,1), expand = c(0,0)) +
  geom_errorbarh(aes(xmin=low,xmax=high),height=0)+
  geom_vline(xintercept = 0.5, linetype = "dashed")+
  guides(color = guide_legend(reverse = TRUE, override.aes = list(shape = 15)), shape = guide_legend(reverse = TRUE))
