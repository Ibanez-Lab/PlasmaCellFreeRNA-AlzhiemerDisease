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
library(DESeq2)
library("tximport")
library(splitstackshape)
library(readxl)
library(doParallel)
library(caret)
library(Rcpp)
library(stringr)
library(glmnet) # to fit a LASSO model
library(MASS) # to use the Boston dataset

library(plyr)

library(tidyverse)
library(skimr)
library(scales)
library(corrr)

library(pls)
library(InformationValue)
library(hrbrthemes)
library(data.table)


setwd("/40/AD/CellFree_Seq_Data/01-cfRNA/202110_NabooAnalyses/ACG")




#Expansion
pheno.expansion <- read.table ("/40/AD/CellFree_Seq_Data/01-cfRNA/072018_Expansion/Phenotype.txt", header=T)
base_dir.expansion<-"/40/AD/CellFree_Seq_Data/01-cfRNA"

pheno.expansion.A = pheno.expansion

pheno.expansion.A = subset(pheno.expansion.A, !(pheno.expansion.A$ID %in% c("65112",
                                                                            "61931",
                                                                            "66055",
                                                                            "65834",
                                                                            "63995")))

pheno.expansion.A$Status [pheno.expansion.A$Group=="CDR1"] <- "cdr"
pheno.expansion.A$Status [pheno.expansion.A$Group=="Control"] <- "CO"
pheno.expansion.A$Status [pheno.expansion.A$Group=="PreSymptomatic"] <- "PreSymptomatic"
pheno.expansion.A$Status [pheno.expansion.A$Group=="CDR05"] <- "CA"

pheno.expansion.A = subset(pheno.expansion.A, Status == "CA" | Status == "CO")


pheno.expansion.A$Phase <- "Expansion"
pheno.expansion.small <- pheno.expansion.A[c(1,18,11,8,19,12,16)]
colnames (pheno.expansion.small) <- c("ID", "Status", "Gender", "Age", "Phase","APOE", "freezertime")



files.expansion <- file.path(base_dir.expansion, "Expansion_hg38", paste0(pheno.expansion.A$ID,".salmonp"),  "quant.genes.sf")

x.expansion<- as.character(pheno.expansion.A$ID)

pheno <-  pheno.expansion.small



#pheno = read.table ("/40/AD/CellFree_Seq_Data/01-cfRNA/082021_ADGroups_DE/phenotype_sampleage.txt", header=T)

rownames (pheno) <- NULL
x <- as.character (pheno$ID)

files.expansion.2 <- as.data.frame (files.expansion)
colnames (files.expansion.2) <- "Files"

files.A <- files.expansion.2
files <- as.character (files.A$Files)
names(files) <- x
all(file.exists(files))



pheno$Status = as.factor(pheno$Status)

library("BiocParallel")
register(MulticoreParam(4))

txi.tx <- tximport(files, type = "salmon", txOut = TRUE, dropInfReps=TRUE)


ddsTxi_Phase <- DESeqDataSetFromTximport(txi.tx, colData = pheno, design= ~ Age + Gender + Status)


dds_phase <- DESeq(ddsTxi_Phase)



saveRDS(dds_phase, "ddsallCDR05old38_notion")


#
# res <- results(dds_phase, pAdjustMethod = "fdr", alpha = .05)
# summary(res)
#
# resdeg = as.data.frame(res)
#
# DEG_PDvsControl = subset(PDvsControl, PDvsControl$padj <= 0.05)
#
#  DEG_DoD_PDvsControl.txt |
# write.table (DEG_PDvsControl, "DEG_DoD_PD+5vsPD-5.txt", quote=FALSE, col.names=T, row.names=T, sep="\t")
# 
# 
# 
# dds_phase = readRDS("ddsallpsold38")
# 
# 
# 
# 
# # # Rlog

 dds_phase = readRDS("ddsallCDR05old38_notion")

genes <- as.data.frame(rownames(dds_phase))
colnames (genes) <- c("GeneName")
rownames(dds_phase) <- genes$GeneName
matrix <- as.matrix (counts(dds_phase))




fewCounts_matrix <- matrix < 10
temp <- rowCounts(fewCounts_matrix)
keepgene <- temp < 0.9*ncol(matrix)
clean.dds = dds_phase[keepgene,]

freezergenes = fread("freezertimefactorrawpvalues.tsv")
freezergenes = freezergenes$V1


clean.dds2 = clean.dds[!rownames(clean.dds) %in% freezergenes, ]

clean.dds = clean.dds2


countsmatri <- as.data.frame(assay(clean.dds)) # 6395

plotDispEsts(clean.dds)
boxplot(log10(assays(clean.dds)[["cooks"]]), range=0, las=2)

res <- results(clean.dds, pAdjustMethod = "fdr", alpha = .05, contrast = c("Status",  "CA", "CO"))
summary(res)


resdf = as.data.frame(res)

# write.table (res, "ExpansionALLCDR05vsCOallgenes.tsv", quote=FALSE, col.names=T, row.names=T, sep="\t")


# res <- results(clean.dds, pAdjustMethod = "fdr", alpha = .1, contrast = c("Status",  "PreSymptomatic", "CO"))
# summary(res)
# 
# 
# write.table (res, "ExpansionALLPSvsCOallgenes.tsv", quote=FALSE, col.names=T, row.names=T, sep="\t")
# 
# 
# resp = subset(res, res$padj < 0.1)
# 
# write.table (resp, "ExpansionALLPSvsCODEG01.tsv", quote=FALSE, col.names=T, row.names=T, sep="\t")
# 
# 
# resraw = subset(res, res$pvalue < 0.05)
# 
# write.table (resraw, "ExpansionALLPSvsCODEGraw.tsv", quote=FALSE, col.names=T, row.names=T, sep="\t")


# rownames(clean.dds) = str_replace_all(str_extract(rownames(clean.dds) , pattern = "\\|(.*?)\\|"), "\\|", "")
# 
# commongenes = fread("commongenesold&new")
# 
# 
# 
# keepgene = rownames(clean.dds) %in% commongenes$gene
# 
# clean.ddsp = clean.dds[keepgene,]
# 
# 
# matrix <- as.matrix (counts(clean.ddsp))
# 
# 
# duplicatedma = matrix[duplicated(row.names(matrix))]
# 
# 
# row.names(matrix)[25]
# 
# encuentraeldu = subset(matrix, row.names(matrix) == "ENSG00000228383.7")
# 
# 
# 
# matrix2 = unique(matrix)
# 
# nombesmatri = row.names(matrix2)
# 
# nombesmatri2 = unique(row.names(matrix2))

# Find duplicate row dates



#clean.dds = dds_phase[rownames(dds_phase) %in% commongenes$gene,]






clean.dds2 <- estimateSizeFactors(clean.dds)


saveRDS(clean.dds2, "ddsallCDR05old38qc_notion")







#rld <- rlog (clean.dds2, blind=FALSE)
rld <- rlog(clean.dds2, blind=FALSE)


saveRDS(rld, "rlogallCDR05old38_notion")

counts <- as.data.frame(assay(rld))
# Counts have negative values
summary(counts)
tdata <- as.data.frame(t(counts))
data <- rownames_to_column (tdata, var = "ID")
data$ID <- sub ("^X", "", data$ID)
data$ID <- sub (".fq.gz", "", data$ID)
pheno$ID <- sub ("^X", "", pheno$ID)
pheno$ID <- sub (".fq.gz", "", pheno$ID)
data.pheno <- merge (pheno, data, by="ID")



saveRDS(data.pheno, file = "rlogmatrix_allCDR05old38_notion") #~/cfRNA_Biology/rlogLauraqcwhole


