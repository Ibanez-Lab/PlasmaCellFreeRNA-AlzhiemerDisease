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
library(stringr)

setwd("/40/AD/CellFree_Seq_Data/01-cfRNA/202110_NabooAnalyses/ACG")

# Load K99 IDs 269 samples
k99 = fread("../k99phenoandtechnical.csv")

# 255 samples
k99 = subset(k99, !(k99$Sample %in% c("GTAC-106700-2102690514.H3HM2DMXY_ATGGCTGT-CCTGTCAA_L001_R1",
                                      "GTAC-106701-2102690512.H3HM2DMXY_TGTTCGCC-AGGTAGGA_L002_R1",
                                      "GTAC-106701-2102690515.H3HM2DMXY_AGGTAGGA-TGTTCGCC_L002_R1",
                                      "GTAC-106703-2102689974.H3J3NDMXY_ATGTTCCT-TCAGCGCC_L002_R1",
                                      "GTAC-106701-2102690373.H3HM2DMXY_GACGTCAT-CCAGTTGA_L002_R1",
                                      "GTAC-106705-2102690346.H3GCVDMXY_AGCCTATC-TGCGTAAC_L002_R1",
                                      "GTAC-106705-2102690358.H3GCVDMXY_TCATCTCC-CTTGCTAG_L002_R1",
                                      "GTAC-106696-2102689983.H3HTCDMXY_CTCGAAAT-CAGGTAAG_L001_R1",
                                      "GTAC-106696-2102690003.H3HTCDMXY_CACTAGAC-TGAGGACT_L001_R1",
                                      "GTAC-106696-2102690267.H3HTCDMXY_AATTAGAC-GTTCTTAT_L001_R1",
                                      "GTAC-106700-2102690053.H3HM2DMXY_ACCCTGAC-GCGCTAAT_L001_R1",
                                      "GTAC-106702-2102690516.H3J3NDMXY_AAGGAAGG-ACCGGAGT_L001_R1",
                                      "GTAC-106702-2102690517.H3J3NDMXY_TCCCACGA-TTCAATAG_L001_R1",
                                      "GTAC-106702-2102690041.H3J3NDMXY_ATCAGAGC-CTGAGCTC_L001_R1")))

k99 = subset(k99, k99$Analysis == "Ethnic" )

summary(as.factor(k99$Group))

k99$Status = "CO" 
k99$Status[k99$Group == "AfACA" ] <- "CA"
k99$Status[k99$Group == "AmIndCA" ] <- "CA"
k99$Status[k99$Group == "AsCA" ] <- "CA"




#219 samples
k99 = subset(k99, k99$PCT_CODING_BASES > 10 )

# 180 samples
k99 = subset(k99, k99$PF_READS > 1000000 )


k99$sex = as.factor(k99$sex)
k99$Status = as.factor(k99$Status)

summary(as.factor(k99$Status))

base_dir<-"/40/AD/CellFree_Seq_Data/01-cfRNA/202110_NabooSalmonFiles/salmon/"
files <- file.path(paste0(base_dir, k99$Sample,".salmon/", "quant.genes.sf"))

x<- as.character(k99$Sample)
names(files) <- x
all(file.exists(files))


library("BiocParallel")
register(MulticoreParam(4))

txi.tx <- tximport(files, type = "salmon", txOut = TRUE, dropInfReps=TRUE)

ddsTxi <- DESeqDataSetFromTximport(txi.tx, colData = k99, design= ~  sex + age_at_draw + Status) # That could be one only for pca plots


dds <- DESeq(ddsTxi)


saveRDS(dds, "ethnicqc") 

#dds <- readRDS("k99allddsstep6")

genes <- as.data.frame(rownames(dds))
colnames (genes) <- c("GeneName")
rownames(dds) <- genes$GeneName
matrix <- as.matrix (counts(dds))



fewCounts_matrix <- matrix < 10
temp <- rowCounts(fewCounts_matrix)
keepgene <- temp < 0.9*ncol(matrix)
clean.dds = dds[keepgene,]


saveRDS(clean.dds, "ethnicqcgenes") 


clean.dds2 <- estimateSizeFactors(clean.dds)


rld <- rlog(clean.dds2, blind=F) # TRUE for QC analyses


saveRDS(rld, "rlogethnicqc") 

counts <- as.data.frame(assay(rld))
# Counts have negative values
summary(counts)
tdata <- as.data.frame(t(counts))
data <- rownames_to_column (tdata, var = "Sample")
data.pheno <- merge (k99, data, by="Sample")


saveRDS(data.pheno, file = "matrix_ethnicqc_sex_age_notion")



