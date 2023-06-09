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



#################### 
#Load K99 IDs
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
                                      "GTAC-106702-2102690041.H3J3NDMXY_ATCAGAGC-CTGAGCTC_L001_R1",
                                      "GTAC-106705-2102690432.H3GCVDMXY_ATTACCCA-CACTGTAG_L002_R1",
                                      "GTAC-106704-2102690439.H3GCVDMXY_AAGGAAGG-ACCGGAGT_L001_R1",
                                      "GTAC-106703-2102690399.H3J3NDMXY_ATAACGCC-CAGGTTCA_L002_R1")))

#219 samples
k99 = subset(k99, k99$PCT_CODING_BASES > 10 )

# 180 samples
k99 = subset(k99, k99$PF_READS > 1000000 )

# Only AD samples 
k99 = subset(k99, k99$Analysis == "AD")


k99 = subset(k99, k99$Group %in% c("Control", "CDR05"))


k99$Group = as.factor(k99$Group)
k99$sex = as.factor(k99$sex)

summary(k99$Group)

k99$Status = "AD" 
k99$Status[k99$Group == "Control" ] <- "CO"


k99$Status = as.factor(k99$Status)


base_dir<-"/40/AD/CellFree_Seq_Data/01-cfRNA/202110_NabooSalmonFiles/salmon/"
files <- file.path(paste0(base_dir, k99$Sample,".salmon/", "quant.genes.sf"))

x<- as.character(k99$Sample)
names(files) <- x
all(file.exists(files))

# # We dont have the following files
# # GTAC-106696-2102690433.H3HTCDMXY_TTCAATAG-TCCCACGA_L001_R1 this one comes from inner QC_K99rlog
# # GTAC-106700-2102690514.H3HM2DMXY_ATGGCTGT-CCTGTCAA_L001_R1
# # GTAC-106701-2102690512.H3HM2DMXY_TGTTCGCC-AGGTAGGA_L002_R1
# # GTAC-106701-2102690515.H3HM2DMXY_AGGTAGGA-TGTTCGCC_L002_R1
#
#
txi.tx <- tximport(files, type = "salmon", txOut = TRUE, dropInfReps=TRUE)

ddsTxi <- DESeqDataSetFromTximport(txi.tx, colData = k99, design= ~ sex + age_at_draw + Status) # sex + age_at_drwa + Status


dds <- DESeq(ddsTxi)




saveRDS(dds, "K99ddsADALLCDR05vsControls_sex_age_notion")

dds <- readRDS("K99ddsADALLCDR05vsControls_sex_age_notion")

genes <- as.data.frame(rownames(dds))
colnames (genes) <- c("GeneName")
#genes.2 <- cSplit(genes, 'GeneName', sep="_", type.convert=FALSE)
rownames(dds) <- genes$GeneName
matrix <- as.matrix (counts(dds))










fewCounts_matrix <- matrix < 10
temp <- rowCounts(fewCounts_matrix)
keepgene <- temp < 0.9*ncol(matrix)
clean.dds = dds[keepgene,]


freezergenes = fread("freezertimefactorrawpvalues.tsv")
freezergenes = freezergenes$V1


clean.dds2 = clean.dds[!rownames(clean.dds) %in% freezergenes, ]

clean.dds = clean.dds2
countsmatri <- as.data.frame(assay(clean.dds))

saveRDS(clean.dds, "K99ddsADALLCDR05vsControlsqc_sex_age_notion") #

clean.dds = readRDS("K99ddsADALLCDR05vsControlsqc_sex_age_notion")

 
 res <- results(clean.dds, pAdjustMethod = "fdr", alpha = .05, contrast = c("Status",  "AD", "CO"))
 summary(res)
 
 boxplot(log10(assays(clean.dds)[["cooks"]]), range=0, las=2)
 
 # 
 # 
 # plotDispEsts(clean.dds)
 # a = boxplot(log10(assays(clean.dds)[["cooks"]]), range=0, las=2)
 # #  
 # a= boxplot(log10(assays(clean.dds)[["cooks"]]), range=0, las=2)

# clean.dds= readRDS('K99ddsADALLCDR05vsControlsqc_sex_age')




# 
# res <- results(clean.dds, pAdjustMethod = "fdr", alpha = .05, contrast = c("Status",  "AD", "CO"))
# summary(res)
# write.table (res, "k99nabooALLCDR05vsCOallgenes.tsv", quote=FALSE, col.names=T, row.names=T, sep="\t")

# 
# resp = subset(res, res$pvalue < 0.05)
# 
# write.table (resp, "k99adqcallcasesDEGrawp.tsv", quote=FALSE, col.names=T, row.names=T, sep="\t")


clean.dds2 <- estimateSizeFactors(clean.dds)


rld <- rlog(clean.dds2, blind=FALSE)


saveRDS(rld, "rldK99ddsADALLCDR05vsControls_sex_age_notion")


# rld <- readRDS("rldK99ddsADALLCDR05vsControls_sex_age")

counts <- as.data.frame(assay(rld))
# Counts have negative values
summary(counts)
tdata <- as.data.frame(t(counts))
data <- rownames_to_column (tdata, var = "Sample")
data.pheno <- merge (k99, data, by="Sample")


saveRDS(data.pheno, file = "matrix_rldK99ddsADALLCDR05vsControls_sex_age_notion")






# rld = readRDS("rldK99ddsADALLCDR05vsControls_sex_age_notion")
# 
# 
# 
# rld_mat <- assay(rld)
# 
# 
# rld_cor <- cor(rld_mat)    ## cor() is a base R function
# 
# head(rld_cor)
# 
# library(pheatmap)
# pheatmap(rld_cor)
# 
# #GTAC-106699-2102690022.H3GHTDMXY_ACGATATG-GACAATTC_L002_R1
# 
# 
# 
# rld <- readRDS("rldK99ddsADALLCDR05vsControls_sex_age_notion")
# 
# 
# 
# 
# 
# k99$PCT_NOT_ALIGNED = (k99$PF_NOT_ALIGNED_BASES/k99$PF_BASES)*100
# 
# k99$yeardrawdate = substr(k99$DrawDate, 1, 4)
# 
# 
# pcaData <- plotPCA(rld, intgroup = "PCT_CODING_BASES", returnData=TRUE)
# percentVar <- round(100 * attr(pcaData, "percentVar"))
# 
# 
# png(file="pca_coding_k99cdr05notion.png", width = 2000, height = 1000, units = "px", pointsize = 12,
#     res = 100)
# ggplot(pcaData, aes(PC1, PC2, color=PCT_CODING_BASES)) +
#   geom_point(size=2) + ggtitle("PCA by % of coding bases") +
#   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#   ylab(paste0("PC2: ",percentVar[2],"% variance")) +
#   coord_fixed() + scale_color_gradientn(colours = rainbow(5))
# dev.off()
# 
# 
# 
# pcaData <- plotPCA(rld, intgroup = "Status", returnData=TRUE)
# percentVar <- round(100 * attr(pcaData, "percentVar"))

# Molarity?
# GTAC-106705-2102690432.H3GCVDMXY_ATTACCCA-CACTGTAG_L002_R1
# GTAC-106704-2102690439.H3GCVDMXY_AAGGAAGG-ACCGGAGT_L001_R1



# c2 <- c(
#     "dodgerblue2", "#E31A1C"
#   )
# 
# 
#   png(file="pca_status_k99_cdr1notion.png", width = 1000, height = 1000, units = "px", pointsize = 12,
#       res = 100)
#   ggplot(pcaData, aes(PC1, PC2, color=Status)) +
#     geom_point(size=1) +
#     xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#     ylab(paste0("PC2: ",percentVar[2],"% variance")) +
#     coord_fixed() + scale_color_manual(values = c2 )
#   dev.off()





# pctpca = plotPCA(rld, intgroup = "PCT_CODING_BASES")
# 
# 
# pca = pctpca$data
# 
# pcadata = pca[,c(1,2,5)]
# 
# colnames(pcadata)[3] = "Sample"
# 
# 
# pheno_qc_pca = inner_join(pcadata,k99, by = "Sample")
# 
# corr = pheno_qc_pca[,c(1,2,19:25,82,27:29,38,66:68,73,83,71)]
# 
# corr$Group = as.numeric(as.factor(corr$Group))
# 
# corr$sex = as.numeric(as.factor(corr$sex))
# 
# 
# corr$yeardrawdate = as.numeric(corr$yeardrawdate)
# 
# library(ggcorrplot)
# 
# 
# p.mat <- cor_pmat(corr)
# 
# corr <- cor(corr)
# 
# 
# png(file="K99step6_correlation.png", width = 1400, height = 1000, units = "px", pointsize = 12,
#     res = 100)
# ggcorrplot(corr,lab = TRUE,lab_size = 3, tl.cex = 10, p.mat=p.mat)
# dev.off()

