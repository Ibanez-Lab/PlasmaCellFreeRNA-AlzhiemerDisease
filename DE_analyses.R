# Load the necessary libraries
library(dplyr)
library(tibble)
library(ggplot2)
library(pROC)
library(reshape)
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
library("BiocParallel")
# Load phenotype and tecnical data
pheno = fread("pheno.csv")
# Remove low quality samples
pheno = subset(pheno, !(pheno$Sample %in% badSamples))
# Remove samples with less than 10% coding bases
pheno = subset(pheno, pheno$PCT_CODING_BASES > 10 )
# Remove samples with less than 1M reads
pheno = subset(pheno, pheno$PF_READS > 1000000 )
# Subset only controls and AD cases of CDR=0
pheno = subset(pheno, pheno$CDR == 0)
# Load salmon data
samlmonFiles <- file.path("path/", "quant.genes.sf"))
salmonData <- tximport(files, type = "salmon", txOut = TRUE, dropInfReps=TRUE)
# Find differentially expressed genes
ddsSalmon <- DESeqDataSetFromTximport(salmonData, colData = pheno, design= ~ sex + age_at_draw + DiseaseStatus)
dds <- DESeq(ddsSalmon)
# Extract normalized counts from teh DESeq2 object
matrix <- as.matrix (counts(dds))
# Remove low count genes from the data
fewCounts_matrix <- matrix < 10
temp <- rowCounts(fewCounts_matrix)
keepgene <- temp < 0.9*ncol(matrix)
clean.dds = dds[keepgene,]
# Remove genes affected by time in freezer of each sample
freezergenes = fread("freezertimefactorrawpvalues.tsv")
freezergenes = freezergenes[,1]
clean.dds <- clean.dds[!rownames(clean.dds) %in% freezergenes,]
# Extract deifferential expression results, excluding low count 
# and freezer time-associated genes
res <- results(clean.dds, pAdjustMethod = "fdr", alpha = .05, contrast = c("DiseaseStatus", "AD", "CO"))
# Subset genes with FDR corrected p-value <0.05
DEgenes = subset(res, res$padj < 0.05)
# Transform gene counts using rlog
clean.dds2 <- estimateSizeFactors(clean.dds)
rld <- rlog(clean.dds2, blind=FALSE)
# PCA
pcaData <- plotPCA(rld, intgroup = "Status", returnData=TRUE)
ggplot(pcaData, aes(PC1, PC2, color=Status)) +
    geom_point(size=1) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed()