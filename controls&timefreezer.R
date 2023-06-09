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

setwd("/40/AD/CellFree_Seq_Data/01-cfRNA/202110_NabooAnalyses/ACG")



# Load K99 IDs
dod = fread("../DoDphenoandtechnical.csv")


dod = subset(dod, dod$PCT_CODING_BASES > 10 )

dod = subset(dod, dod$TOTAL_READS > 1000000 )

summary(as.factor(dod$Group))


dod$Status = "PD" 
dod$Status[dod$Group == "Control" ] <- "CO"
dod$Status[dod$Group == "CO" ] <- "CO"
dod$Status[dod$Group == "Control-PD" ] <- "CO"
dod$Status[dod$Group == "AfACO" ] <- "AfACO"

dod = subset(dod, Status == "CO")

dod$Status = as.factor(dod$Status)
dod$sex = as.factor(dod$sex)

summary(as.factor(dod$Status))




k99 = fread("../k99phenoandtechnical.csv")

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
                                      "GTAC-106702-2102690041.H3J3NDMXY_ATCAGAGC-CTGAGCTC_L001_R1"
)))



k99 = subset(k99, k99$PCT_CODING_BASES > 10 )

k99 = subset(k99, k99$TOTAL_READS > 1000000 )

summary(as.factor(k99$Group))

k99$Status = "other" 
k99$Status[k99$Group == "Control" ] <- "CO"

k99 = subset(k99, Group == "Control")

summary(as.factor(k99$Status))



k99anddod = rbind(k99,dod)

k99anddod$freezertime = 2021 - as.numeric(substr(k99anddod$DrawDate, 1, 4))
k99anddod$sex = k99anddod(fusion$sex)

base_dir<-"/40/AD/CellFree_Seq_Data/01-cfRNA/202110_NabooSalmonFiles/salmon/"
files <- file.path(paste0(base_dir, k99anddod$Sample,".salmon/", "quant.genes.sf"))

x<- as.character(k99anddod$Sample)
names(files) <- x
all(file.exists(files))





library("BiocParallel")
register(MulticoreParam(4))

txi.tx <- tximport(files, type = "salmon", txOut = TRUE, dropInfReps=TRUE)

ddsTxi <- DESeqDataSetFromTximport(txi.tx, colData = k99anddod, design= ~ sex + age_at_draw + freezertime) # That could be one only for pca plots


dds <- DESeq(ddsTxi)

saveRDS(dds, "ddsmulticontrols_sex_age")



# dds <- readRDS("ddsmulticontrols")

genes <- as.data.frame(rownames(dds))
colnames (genes) <- c("GeneName")
#genes.2 <- cSplit(genes, 'GeneName', sep="_", type.convert=FALSE)
rownames(dds) <- genes$GeneName
matrix <- as.matrix (counts(dds))



fewCounts_matrix <- matrix < 10
temp <- rowCounts(fewCounts_matrix)
keepgene <- temp < 0.9*ncol(matrix)
clean.dds = dds[keepgene,]


# saveRDS(clean.dds, "ddsdodPDqc_sex_age") # 
# 
# 
# res <- results(dds, pAdjustMethod = "bonferroni", alpha = .05)
# summary(res) 

#write.table (res, "DOD.tsv", quote=FALSE, col.names=T, row.names=T, sep="\t")

#clean.dds= readRDS('DoDqcdds')


clean.dds2 <- estimateSizeFactors(clean.dds)


rld <- rlog(clean.dds2, blind=F)



saveRDS(rld, "rlogmulticontrols_sex_age")




counts <- as.data.frame(assay(rld))
# Counts have negative values
summary(counts)
tdata <- as.data.frame(t(counts))
data <- rownames_to_column (tdata, var = "Sample")
data.pheno <- merge (k99anddod, data, by="Sample")


saveRDS(data.pheno, file = "matrix_multicontrols_sex_age")


rld = readRDS("rlogmulticontrols_sex_age")



png(file="pcacontrols3.png", width = 1000, height = 1000, units = "px", pointsize = 12,
    res = 100)
plotPCA(rld, intgroup = "Source")

dev.off()





pcaData <- plotPCA(rld, intgroup = "Source", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

library(RColorBrewer)

# c25 <- c(
#   "dodgerblue2", "#E31A1C", # red
#   "green4",
#   "#6A3D9A", # purple
#   "#FF7F00", # orange
#   "black", "gold1",
#   "skyblue2", "#FB9A99", # lt pink
#   "palegreen2",
#   "#CAB2D6", # lt purple
#   "#FDBF6F", # lt orange
#   "gray70", "khaki2",
#   "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
#   "darkturquoise", "green1", "yellow4", "yellow3",
#   "darkorange4", "brown"
# )
# 
# 
# c14 <- c(
#   "dodgerblue2", "#E31A1C", # red
#   "green4",
#   "#6A3D9A", # purple
#   "#FF7F00", # orange
#   "black", "gold1",
#   "skyblue2", "#FB9A99", # lt pink
#   "palegreen2",
#   "#CAB2D6", # lt purple
#   "#FDBF6F", # lt orange
#   "gray70", "brown"
# )


c3 <- c(
  "#E31A1C", # red
  #"green4",
  "#6A3D9A", # purple
  # "#FF7F00", # orange
  "black" #"gold1",
  # "skyblue2", "#FB9A99", # lt pink
  # "palegreen2",
  # "#CAB2D6", # lt purple
  # "#FDBF6F", # lt orange
  # "gray70", "brown"
)


png(file="pcacontrols3.png", width = 1000, height = 1000, units = "px", pointsize = 12,
    res = 100)
ggplot(pcaData, aes(PC1, PC2, color=Source)) +
  geom_point(size=2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() + scale_color_manual(values = c3 )
dev.off()

