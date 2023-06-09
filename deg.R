# "DEG analyses paper"

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
library(ggrepel)



clean.dds= readRDS('K99ddsADPreclinicalandControlsqc_sex_agep_notion')

res <- results(clean.dds, pAdjustMethod = "fdr", alpha = .05, contrast = c("Status",  "AD", "CO"))
summary(res)

write.table (res, "k99presimnotion.tsv", quote=FALSE, col.names=T, row.names=T, sep="\t")
write.table (res, "k99preclinical.tsv", quote=FALSE, col.names=T, row.names=T, sep="\t")




dds_phase = readRDS("ddsallpsold38_notion")

genes <- as.data.frame(rownames(dds_phase))
colnames (genes) <- c("GeneName")
rownames(dds_phase) <- genes$GeneName
matrix <- as.matrix (counts(dds_phase))


fewCounts_matrix <- matrix < 10
temp <- rowCounts(fewCounts_matrix)
keepgene <- temp < 0.9*ncol(matrix)
clean.dds = dds_phase[keepgene,]

freezergenes = fread("freezertimefactorrawpvalues.tsv")
names(freezergenes) = c("num", "V1",             "baseMean",       "log2FoldChange", "lfcSE", "stat", "pvalue", "padj" )
freezergenes = freezergenes$V1


clean.dds2 = clean.dds[!rownames(clean.dds) %in% freezergenes, ]

clean.dds = clean.dds2



res <- results(clean.dds, pAdjustMethod = "fdr", alpha = .05, contrast = c("Status",  "PreSymptomatic", "CO"))
summary(res)
write.table (res, "ExpansionALLPSvsCOallgenes.tsv", quote=FALSE, col.names=T, row.names=T, sep="\t")




dds <- readRDS("ddsmulticontrols_sex_age_factor")

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
res <- results(clean.dds, pAdjustMethod = "fdr", alpha = .05, contrast = c("freezertimefactor",  "OLD", "NEW"))
summary(res)
write.table (res, "freezertimefactor.tsv", quote=FALSE, col.names=T, row.names=T, sep="\t")


k99dodfrezer = fread("freezertimefactor.tsv") 
k99dodfrezer = subset(k99dodfrezer, k99dodfrezer$pvalue < 0.05)
write.table (k99dodfrezer, "freezertimefactorrawpvalues.tsv", quote=FALSE, col.names=T, row.names=T, sep="\t")








# Load DEG

## Load preclinicals AD vs controls DEG from naboo




genesdegk99 = fread("k99preclinical.tsv")

volcanoplotdata = genesdegk99

genesdegk99 = subset(genesdegk99, genesdegk99$padj < 0.05)

# write.table (genesdegk99, "k99preclinicaldeg.tsv", quote=FALSE, col.names=T, row.names=F, sep="\t")

deggenesensg = c("ENSG00000100034.14", "ENSG00000115904.13", "ENSG00000148843.15", "ENSG00000159314.11", "ENSG00000182809.10", "ENSG00000186350.12", "ENSG00000198844.12")





preclinicalgenes = fread("preclinicaldegsymbols.csv", header = F)

preclinicalgenes = preclinicalgenes[order(preclinicalgenes$V2),]

colnames(preclinicalgenes) = c("ENSG", "GENE_SYMBOL")

# write.table (preclinicalgenes, "preclinicalgenessort.tsv", quote=FALSE, col.names=T, row.names=F, sep="\t")


colnames(genesdegk99)[1] = "ENSG"

finaltablegene = inner_join(genesdegk99,preclinicalgenes, by = ("ENSG"))

finaltablegene = finaltablegene[order(finaltablegene$GENE_SYMBOL),]

finaltablegene = finaltablegene[,c(1,8,2:7)]


# write.table (finaltablegene, "degpreclinicaltableforpaper.tsv", quote=FALSE, col.names=T, row.names=F, sep="\t")

genesexpan = fread("expansionpreclinicals.tsv")

genesexpan = genesexpan[,c(1:6)]

colnames(genesexpan) = c("ENSG", "Replication_baseMean", "Replication_ log2FoldChange", "Replication_lfcSE", "Replication_stat", "Replication_pvalue")

commondeg = inner_join(finaltablegene,genesexpan, by = ("ENSG"))

commondeg = subset(commondeg, commondeg$Replication_pvalue < 0.05)

commondeg = subset(commondeg, commondeg$log2FoldChange > 0 & commondeg$`Replication_ log2FoldChange` > 0 | commondeg$log2FoldChange < 0 & commondeg$`Replication_ log2FoldChange` < 0)

# write.table (commondeg, "commondegpreclinicaltableforpaper.tsv", quote=FALSE, col.names=T, row.names=F, sep="\t")





## Brain specific



gtexdata = fread("gtexdata.gct")



mediacortexgtex = mean(gtexdata$`Brain - Cortex`)


colnames(gtexdata)[1] = "ENSG"


tabletogtex = finaltablegene[,1:2]

write.table (tabletogtex, "checknames.tsv", quote=FALSE, col.names=T, row.names=F, sep="\t")

infogtex = inner_join(tabletogtex,gtexdata, by = "ENSG")

colnames(gtexdata)[2] = "GENE_SYMBOL"




infogtexsymbol = inner_join(tabletogtex,gtexdata, by = "GENE_SYMBOL")





infogtexsymbol = subset(infogtexsymbol, rowSums(infogtexsymbol[,4:57]) > 0)

gtexexpressioncortex = subset(infogtexsymbol, rowSums(infogtexsymbol[,16]) > 0)

gtexexpressioncortex = gtexexpressioncortex[,c(1,2,16)]


write.table (gtexexpressioncortex, "gtexexpressioncortex.tsv", quote=FALSE, col.names=T, row.names=F, sep="\t")


mediacortexgenes = mean(infogtexsymbol$`Brain - Cortex`)





gtexdataplot = gtexdata


gtexdataplot$source = "GTEX"

gtexdataplot$source[gtexdataplot$GENE_SYMBOL %in%  tabletogtex$GENE_SYMBOL] <- "DE"




#gtexplot = 
# ggplot(gtexdataplot, aes(x = source, y = `Brain - Cortex`, fill = source )) +
# geom_boxplot(alpha = 0.5)  +
# theme_minimal()  + xlab("") + ylab("Gene expression in brain cortex (TPM)") + labs(title="Gene expression in cortex TPM for all genes ") + scale_y_log10(breaks = c(0, 5, 10, 20, 40, 100, 500, 1000, 10000, 60000))
# 





gtexnewdata = subset(gtexdataplot, rowSums(gtexdataplot[,3:56]) > 0)

# gtexdataplotb = gtexdataplot[,-c(3:9,23:56)]


gtexdataplotg =  gather(gtexnewdata, "tissiue", "expression", 3:56)

gtexplot =    ggplot(gtexdataplotg, aes(x = tissiue, y = expression, fill = source )) +
  geom_boxplot(alpha = 0.5)  +
  theme_minimal()  +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("") + ylab("Gene expression in brain cortex (TPM)") + labs(title="Gene expression (TPM) in cortex for DE and all genes in GTEX") + scale_y_log10()

gtexplot

ggsave(filename = "gtexdataplot.png",
       plot = gtexplot,
       device = "png",
       width =12,
       height = 9,
       dpi = 300)




symbolcommon = c("ARHGAP27", "ARHGEF15", "CRIP2", "PDCD11", "PPM1F", "RXRA", "SOS1")



gtexdata = subset(gtexdata, gtexdata$GENE_SYMBOL %in% symbolcommon)


phyper(23 - 1, 190, 20000, 190, lower.tail = FALSE)


# Brain fold change


gtexdata = fread("gtexdata.gct")

summary(gtexdata)

prueba = gtexdata[gtexdata$`Brain - Cortex` >= 3*gtexdata$`Whole Blood`]

prueba = gtexdata[gtexdata$`Brain - Cortex` >= 3*gtexdata[,c(3:9,22:56)]]




gtexdatamp = subset(gtexdata, gtexdata$`Brain - Cortex` >= 2*gtexdata$`Adipose - Subcutaneous` & gtexdata$`Brain - Cortex` >= 2*gtexdata$`Adipose - Visceral (Omentum)` & gtexdata$`Brain - Cortex` >= 2*gtexdata$`Adrenal Gland` & gtexdata$`Brain - Cortex` >= 2*gtexdata$`Artery - Aorta` & gtexdata$`Brain - Cortex` >= 2*gtexdata$`Artery - Coronary` & gtexdata$`Brain - Cortex` >= 2*gtexdata$`Artery - Tibial` & gtexdata$`Brain - Cortex` >= 2*gtexdata$Bladder & gtexdata$`Brain - Cortex` >= 2*gtexdata$`Breast - Mammary Tissue` & gtexdata$`Brain - Cortex` >= 2*gtexdata$`Cells - Cultured fibroblasts` & gtexdata$`Brain - Cortex` >= 2*gtexdata$`Cells - EBV-transformed lymphocytes` & gtexdata$`Brain - Cortex` >= 2*gtexdata$`Cervix - Ectocervix` & gtexdata$`Brain - Cortex` >= 2*gtexdata$`Cervix - Endocervix` & gtexdata$`Brain - Cortex` >= 2*gtexdata$`Colon - Sigmoid` & gtexdata$`Brain - Cortex` >= 2*gtexdata$`Colon - Transverse` & gtexdata$`Brain - Cortex` >= 2*gtexdata$`Esophagus - Gastroesophageal Junction`)


gtexdatamp = subset(gtexdatamp, gtexdatamp$Description %in% infogtexsymbol$GENE_SYMBOL)


gtexdatam = scale(gtexdata[,3:56])

gtexdatam = as.data.frame(gtexdatam)

gtexdatam$symbol = gtexdata$Description






## Brain DE results

### Sunshine



sunshine = fread("results_Sunshine.csv")

# sunshinep = sunshine
# 
# sunshinep$ajustado = p.adjust(sunshine$pvalue, method = "fdr", n = length(sunshine$pvalue))


# number of genes below 0.05

hipersunshine = fread("results_all_Sunshine.csv")

hipersunshine = subset(hipersunshine, hipersunshine$pvalue < 0.05)

GENESNAME = fread("19000genes.csv")

colnames(GENESNAME)[2] = "V1"

hipersunshinef = inner_join(hipersunshine, GENESNAME, by = "V1")


hipersunshine = subset(hipersunshine, hipersunshine$padj < 0.05)

hipersunshinef = inner_join(hipersunshine, GENESNAME, by = "V1")

colnames(sunshine) = c("GENE_SYMBOL",             "sunshine_baseMean",       "sunshine_log2FoldChange", "sunshine_lfcSE",          "sunshine_stat", "sunshine_pvalue" ,        "sunshine_padj")



sunshineandde = inner_join(finaltablegene,sunshine, by = "GENE_SYMBOL")

sunshineandde = na.omit(sunshineandde)


sunshineanddec = subset(sunshineandde, sunshineandde$sunshine_log2FoldChange > 0 & sunshineandde$log2FoldChange > 0 |  sunshineandde$sunshine_log2FoldChange < 0 & sunshineandde$log2FoldChange < 0)

sunshineanddecruce = sunshineanddec
sunshineanddec = subset(sunshineanddec, sunshineanddec$sunshine_pvalue < 0.05)



# write.table (sunshineanddec, "sunshineanddec.tsv", quote=FALSE, col.names=T, row.names=F, sep="\t")


gtexdatasunshine = subset(gtexdata, gtexdata$GENE_SYMBOL %in% sunshineanddec$GENE_SYMBOL)

#  | commondeg$delog < 0 & commondeg$log2FoldChange < 0


# MvsS = fread("results_MvsS.csv")
# 
# colnames(MvsS) = c("GENE_SYMBOL",             "MvsS_baseMean",       "MvsS_log2FoldChange", "MvsS_lfcSE",          "MvsS_stat", "MvsS_pvalue" ,        "MvsS_padj")
# 
# 
# 
# MvsSdde = inner_join(finaltablegene,MvsS, by = "GENE_SYMBOL")
# 
# MvsSdde = na.omit(MvsSdde)
# 
# 
# MvsSddec = subset(MvsSdde, MvsSdde$MvsS_log2FoldChange > 0 & MvsSdde$log2FoldChange > 0 |  MvsSdde$MvsS_log2FoldChange < 0 & MvsSdde$log2FoldChange < 0)
# 
# MvsSddecruce = MvsSdde
# 
# MvsSddec = subset(MvsSddec, MvsSddec$MvsS_pvalue < 0.05)
# 
# write.table (MvsSddec, "MvsSddec.tsv", quote=FALSE, col.names=T, row.names=F, sep="\t")
# 
# #  | commondeg$delog < 0 & commondeg$log2FoldChange < 0
# 
# 
# 
# 
# sunshineanddecruce = sunshineanddecruce[,c(2,9:14)]
# 
# crucedebrain = inner_join(MvsSddecruce, sunshineanddecruce, by = "GENE_SYMBOL")
# 
# crucedebrain = subset(crucedebrain, crucedebrain$MvsS_log2FoldChange > 0 & crucedebrain$sunshine_log2FoldChange > 0 |  crucedebrain$MvsS_log2FoldChange < 0 & crucedebrain$sunshine_log2FoldChange < 0)
# 
# write.table (crucedebrain, "crucedebrain.tsv", quote=FALSE, col.names=T, row.names=F, sep="\t")






# Volcano plot



volcanoplotdata = fread("volcanoplot.tsv")

volcanoplotdata = na.omit(volcanoplotdata)

volcanoplotdata$expression = ifelse(volcanoplotdata$padj < 0.05 & abs(volcanoplotdata$log2FoldChange) >= 0.45, 
                                    ifelse(volcanoplotdata$log2FoldChange> 0.45 ,'Up','Down'),
                                    'Stable')


# write.table (volcanoplotdata, "volcanoplot.tsv", quote=FALSE, col.names=T, row.names=F, sep="\t")

# Included the gene symbol from https://www.biotools.fr/human/ensembl_symbol_converter




# genes in common with brain or expansion c("SYNPO", "TRAK2", "SPTBN1", "JCAD", "FHDC1", "RN7SL4P", "MYL6", "ZMYND8", "H3F3A", "MBOAT2")
interestingenens = unique(c(sunshineanddec$GENE_SYMBOL, commondeg$GENE_SYMBOL))


downgenestop5 = c("RN7SL4P", "ZMYND8", "H3F3A" , "RBX1" ,"CENPBD1P1")

interestingenens = c(interestingenens,"SYNPO", "TRAK2", "SPTBN1", "JCAD", "FHDC1", "FP671120.3", "AIF1L", "KIF26A", "ERG" , "FP236383.2",downgenestop5)


interestingenens = unique(interestingenens)

volcanoplotdata$expression = ifelse(volcanoplotdata$padj < 0.05 & abs(volcanoplotdata$log2FoldChange) >= 0.45, 
                                    ifelse(volcanoplotdata$log2FoldChange> 0.45 ,'Up','Down'),
                                    'Stable')

# Count top genes
volcanoplotdata$delabel <- NA

volcanoplotdata$delabel[volcanoplotdata$gene_symbol %in% interestingenens] <-volcanoplotdata$gene_symbol[volcanoplotdata$gene_symbol %in% interestingenens]





# Source of the genes
volcanoplotdata$sourcecol <- "No replicated"

#interestingenens = unique(c(sunshineanddec$GENE_SYMBOL, commondeg$GENE_SYMBOL))

volcanoplotdata$sourcecol[volcanoplotdata$gene_symbol %in% sunshineanddec$GENE_SYMBOL] <- "Brain replicated"

volcanoplotdata$sourcecol[volcanoplotdata$gene_symbol %in% commondeg$GENE_SYMBOL] <- "Plasma replicated"




volcanoplotdata$sourcecol = as.factor(volcanoplotdata$sourcecol)

p <- ggplot(data = volcanoplotdata, 
            aes(x = volcanoplotdata$log2FoldChange, 
                y = -log10(volcanoplotdata$padj), 
                color=volcanoplotdata$sourcecol,
                shape=volcanoplotdata$sourcecol,
                label = volcanoplotdata$delabel)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_shape_manual(values=c(19, 18, 15))+
  scale_color_manual(values=c('blue','grey', 'black'))+
  xlim(c(-2.55, 2.55)) +  
  geom_text_repel(size = 4, max.overlaps = 30 ,min.segment.length = 0 , segment.size = 0.35, segment.alpha	= 0.7, 
                  segment.curvature = -0.1,
                  segment.ncp = 3,
                  segment.angle = 25, force =1, show.legend = FALSE) +
  geom_vline(xintercept=c(-0.45,0.45),lty=1,col="black",lwd=0.2) +
  geom_hline(yintercept = 1.301,lty=1,col="black",lwd=0.2) +
  labs(x="log2(fold change)",
       y="-log10 (adj.p-value)",
       title="Differential expression between preclinical AD participants and controls")  +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="bottom", 
        legend.title = element_blank(),text = element_text(size = 15)) +guides(label = FALSE)

p

ggsave(filename = "volcanoplot.png",
       plot = p,
       device = "png",
       width =14,
       height = 9,
       dpi = 300)




genesdegk99 = fread("k99preclinical.tsv")

genesdegk99 = subset(genesdegk99, genesdegk99$padj < 0.05)

# write.table (genesdegk99, "k99preclinicaldeg.tsv", quote=FALSE, col.names=T, row.names=F, sep="\t")



preclinicalgenes = fread("preclinicaldegsymbols.csv", header = F)

preclinicalgenes = preclinicalgenes[order(preclinicalgenes$V2),]

colnames(preclinicalgenes) = c("ENSG", "GENE_SYMBOL")

# write.table (preclinicalgenes, "preclinicalgenessort.tsv", quote=FALSE, col.names=T, row.names=F, sep="\t")


colnames(genesdegk99)[1] = "ENSG"

finaltablegene = inner_join(genesdegk99,preclinicalgenes, by = ("ENSG"))

finaltablegene = finaltablegene[order(finaltablegene$GENE_SYMBOL),]

finaltablegene = finaltablegene[,c(1,8,2:7)]




dbs <- listEnrichrDbs()


websiteLive <- TRUE

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")




genestoenrich = finaltablegene$GENE_SYMBOL
#genestoenrich = sunshineanddec$GENE_SYMBOL




dbsu <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021", "DisGeNET", "GTEx_Tissue_Sample_Gene_Expression_Profiles_up", "KEGG_2021_Human", "Azimuth_Cell_Types_2021")
if (websiteLive) {
  enriched <- enrichr(genestoenrich, dbsu)
}







gobio = if (websiteLive) enriched[["GO_Biological_Process_2021"]]
gocelular = if (websiteLive) enriched[["GO_Cellular_Component_2021"]]
gomolecular = if (websiteLive) enriched[["GO_Molecular_Function_2021"]]
disgenet = if (websiteLive) enriched[["DisGeNET"]]
gtex = if (websiteLive) enriched[["GTEx_Tissue_Sample_Gene_Expression_Profiles_up"]]
kegg = if (websiteLive) enriched[["KEGG_2021_Human"]]
azimuthcell = if (websiteLive) enriched[["Azimuth_Cell_Types_2021"]]


# kegg




if (websiteLive) plotEnrich(enriched[["KEGG_2021_Human"]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")

# Venn diagram 

preclinicalgenes = fread("preclinicaldegsymbols.csv", header = F)
preclinicalgenes = preclinicalgenes$V2
genes40  = fread("40genes.csv", header = T)
genes40  = genes40$V6
genes90 = fread("90genes.csv", header = T)
genes90 = genes90$V6
genes220 = fread("220genes.csv", header = T)
genes220 = genes220$SYMBOL_GENE


listade <- list(
  model40 = genes40, 
  model90 = genes90,
  model220 = genes220
)


library(ggvenn)

png(file="modelsbetween.png",
    width=2200, height=1700, res = 300)
ggvenn(
  listade, 
  fill_color = c("#0073C2FF", "#EFC000FF",  "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)
dev.off()




listade <- list(
  DEG = preclinicalgenes, 
  model40 = genes40
)

png(file="modelde40.png",
    width=2200, height=1700, res = 300)
ggvenn(
  listade, 
  fill_color = c("#868686FF", "#0073C2FF"),
  stroke_size = 0.5, set_name_size = 4
)
dev.off()


listade <- list(
  DEG = preclinicalgenes, 
  model90 = genes90
)

png(file="modelde90.png",
    width=2200, height=1700, res = 300)
ggvenn(
  listade, 
  fill_color = c("#868686FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 4
)
dev.off()



listade <- list(
  DEG = preclinicalgenes, 
  model220 = genes220
)

png(file="modelde220.png",
    width=2200, height=1700, res = 300)
ggvenn(
  listade, 
  fill_color = c("#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)
dev.off()



finaltablegene

genes220 = fread("220genes.csv", header = T)
genes220 = genes220$SYMBOL_GENE


genesoverlap220yde = subset(finaltablegene, finaltablegene$GENE_SYMBOL %in% genes220)




# Venn 


genesdegk99 = fread("k99preclinical.tsv")

preclinicalgenes = fread("preclinicaldegsymbols.csv", header = F)


preclinicalnumber = inner_join(preclinicalgenes, genesdegk99, by = "V1")

preupregulated = subset(preclinicalnumber, preclinicalnumber$log2FoldChange > 0)

preupregulateddf = preupregulated

preupregulated = preupregulated$V2


predownregulated = subset(preclinicalnumber, preclinicalnumber$log2FoldChange < 0)

predownregulateddf = predownregulated

predownregulated = predownregulated$V2



ibarraup = read_excel("ibarrasupplementarygenes.xlsx", "Upregulated genes", col_names = T)

colnames(ibarraup) = c("ensg", "V2")


commonup = inner_join(preupregulateddf, ibarraup, by = "V2")


ibarradown = read_excel("ibarrasupplementarygenes.xlsx", "Downregulated genes", col_names = T)


colnames(ibarradown) = c("ensg", "V2")


commondown = inner_join(predownregulateddf, ibarradown, by = "V2")


ibarraup = ibarraup[-1,]


ibarradown = ibarradown[-1,]


ibarradown = ibarradown$V2

ibarraup = ibarraup$V2

ibarracommon = c(commondown$V2, commonup$V2)

writeLines(ibarracommon, "ibarracommon.csv")

writeLines(sunshineanddec$GENE_SYMBOL, "braincommon.csv")


checknames = fread("checknames.tsv")

antinamles = fread("antijoint.csv")



colnames(antinamles) = "GENE_SYMBOL"


CSEAana = anti_join(checknames, antinamles, by = "GENE_SYMBOL")


write_tsv(CSEAana, "genescsea.tsv")

listade <- list(
  ibarraup = ibarraup, 
  preupregulated = preupregulated
)

png(file="ibarraupvenn.png",
    width=2200, height=1700, res = 300)
ggvenn(
  listade, 
  fill_color = c("#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)
dev.off()




listade <- list(
  ibarradown = ibarradown, 
  predownregulated = predownregulated
)

png(file="ibarradownvenn.png",
    width=2200, height=1700, res = 300)
ggvenn(
  listade, 
  fill_color = c("#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)
dev.off()



# Volcano paper

volcanoplotdata = fread("volcanoplot.tsv")

volcanoplotdata = na.omit(volcanoplotdata)

volcanoplotdata$expression = ifelse(volcanoplotdata$padj < 0.05 & abs(volcanoplotdata$log2FoldChange) >= 0.45, 
                                    ifelse(volcanoplotdata$log2FoldChange> 0.45 ,'Up','Down'),
                                    'Stable')


# write.table (volcanoplotdata, "volcanoplot.tsv", quote=FALSE, col.names=T, row.names=F, sep="\t")

# Included the gene symbol from https://www.biotools.fr/human/ensembl_symbol_converter




# genes in common with brain or expansion c("SYNPO", "TRAK2", "SPTBN1", "JCAD", "FHDC1", "RN7SL4P", "MYL6", "ZMYND8", "H3F3A", "MBOAT2")
interestingenens = unique(c(sunshineanddec$GENE_SYMBOL, ibarracommon))


downgenestop5 = c("RN7SL4P", "ZMYND8", "H3F3A" , "RBX1" ,"CENPBD1P1")
interestingenens = c(interestingenens,"SYNPO", "TRAK2", "SPTBN1", "JCAD", "FHDC1", "FP671120.3", "AIF1L", "KIF26A", "ERG" , "FP236383.2",downgenestop5)


interestingenens = unique(interestingenens)

volcanoplotdata$expression = ifelse(volcanoplotdata$padj < 0.05 & abs(volcanoplotdata$log2FoldChange) >= 0.45, 
                                    ifelse(volcanoplotdata$log2FoldChange> 0.45 ,'Up','Down'),
                                    'Stable')

# Count top genes
volcanoplotdata$delabel <- NA

volcanoplotdata$delabel[volcanoplotdata$gene_symbol %in% interestingenens] <-volcanoplotdata$gene_symbol[volcanoplotdata$gene_symbol %in% interestingenens]





# Source of the genes
volcanoplotdata$sourcecol <- "No Replicated"

#interestingenens = unique(c(sunshineanddec$GENE_SYMBOL, commondeg$GENE_SYMBOL))

volcanoplotdata$sourcecol[volcanoplotdata$gene_symbol %in% sunshineanddec$GENE_SYMBOL] <- "Replicated in Brain"

volcanoplotdata$sourcecol[volcanoplotdata$gene_symbol %in% ibarracommon] <- "Replicated in Plasma"


volcanoplotdata$sourcecol[volcanoplotdata$gene_symbol %in% ibarracommon &  volcanoplotdata$gene_symbol %in% sunshineanddec$GENE_SYMBOL] <- "Replicated in Plasma and Brain"



volcanoplotdata$sourcecol = as.factor(volcanoplotdata$sourcecol)


volcanoplotdata$sourcecol <- factor(volcanoplotdata$sourcecol, levels = c("No Replicated", "Replicated in Plasma", "Replicated in Brain", "Replicated in Plasma and Brain"))

library(oneclust)

cud(1)

cud(3)




p <- ggplot(data = volcanoplotdata, 
            aes(x = volcanoplotdata$log2FoldChange, 
                y = -log10(volcanoplotdata$padj), 
                color=volcanoplotdata$sourcecol,
                shape=volcanoplotdata$sourcecol,
                label = volcanoplotdata$delabel)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_shape_manual(values=c(18, 18, 18, 19))+
  scale_color_manual(values=c('grey',  '#0072b2', "#009e73", 'black'))+
  xlim(c(-2.55, 2.55)) +  
  geom_text_repel(size = 3.5, max.overlaps = 30 ,min.segment.length = 0 , segment.size = 0.35, segment.alpha	= 0.7, 
                  segment.curvature = -0.1,
                  segment.ncp = 3,
                  segment.angle = 25, force =1, show.legend = FALSE) +
  geom_vline(xintercept=c(-0.45,0.45),lty=1,col="black",lwd=0.2) +
  geom_hline(yintercept = 1.301,lty=1,col="black",lwd=0.2) +
  labs(x="log2(fold change)",
       y="-log10 (adj.p-value)",
       title="Differential expression between preclinical AD participants and controls")  +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="bottom", 
        legend.title = element_blank(),text = element_text(size = 15)) +guides(label = FALSE)  + annotate(geom="text", x=2.25, y=6.7, label="Up-regulated",
                                                                                                          color="black", size = 4) + annotate(geom="text", x=-2.25, y=6.7, label="Down-regulated",
                                                                                                                                              color="black", size = 4)

p

ggsave(filename = "volcanoplot3.png",
       plot = p,
       device = "png",
       width =14,
       height = 9,
       dpi = 600)


