##################################################
############# CopyscAT Setup #####################
##################################################
library(devtools)
library(GenomicRanges)
library(future)
library(harmony)
library(SeuratDisk)
library(SeuratObject)
library(Seurat)
library(Signac)
library(ggplot2)
library(glue)
library(scDblFinder)
library(RColorBrewer)
library(dplyr)
library(ggridges)
library(CopyscAT)
library(BSgenome.Hsapiens.UCSC.hg38)
library(RColorBrewer)
#install_github("spcdot/copyscat")

bs <- BSgenome.Hsapiens.UCSC.hg38
generateReferences(BSgenome.Hsapiens.UCSC.hg38,genomeText = "hg38",tileWidth = 1e6,outputDir = glue("/Volumes/VJ_Data/single cell ERV work/scATAC-seq outputs/QC_2"))
table(combined@meta.data$dataset)

initialiseEnvironment(genomeFile=".../hg38_chrom_sizes.tsv",
                      cytobandFile=".../hg38_1e+06_cytoband_densities_granges.tsv",
                      cpgFile=".../hg38_1e+06_cpg_densities.tsv",
                      binSize=1e6,
                      minFrags=3e4,
                      cellSuffix=c("-1","-2"),
                      lowerTrim=0.5,
                      upperTrim=0.8)

summaryFunction <- cutAverage
setwd("...")

##########################################################################################
############# Run sample by sample ############### 84, 86 and 96 failed,  ################
##########################################################################################

Idents(combined) <- "dataset"
factor(levels(combined))
combined_RCC87 <- subset(combined, idents="RCC87")
scData<-readInputTable(glue(".../RCC87_copyscat2.tsv"))
scData_k_norm <- normalizeMatrixN(scData,logNorm = FALSE,maxZero=2000,imputeZeros = FALSE,blacklistProp = 0.8,blacklistCutoff=125,dividingFactor=1,upperFilterQuantile = 0.95)
scData_collapse<-collapseChrom3N(scData_k_norm,summaryFunction=summaryFunction,binExpand = 1,minimumChromValue = 100,logTrans = FALSE,tssEnrich = 1,logBase=2,minCPG=300,powVal=0.73) 
scData_collapse<-filterCells(scData_collapse,minimumSegments = 5,minDensity = 0.1)
norm_barcodes = row.names(combined_RCC87[[]] %>% filter(seurat_clusters %in% c('1','2','3','4','5'))) #Normal clusters were empirically determined using CovPlot signal
norm_barcodes <- trimws(norm_barcodes, whitespace="^RCC_87_")
normal_barcodes_manual_qc = intersect(norm_barcodes,colnames(scData_collapse))

graphCNVDistribution(scData_collapse,outputSuffix = "unscaled")

median_iqr <- computeCenters(scData_collapse %>% dplyr::select(chrom,normal_barcodes_manual_qc),summaryFunction=summaryFunction)
candidate_cnvs<-identifyCNVClusters(scData_collapse,median_iqr,useDummyCells = TRUE,propDummy=0.25,minMix=0.01,deltaMean = 0.03,deltaBIC2 = 0.25,bicMinimum = 0.1, subsetSize=800,fakeCellSD = 0.09, uncertaintyCutoff = 0.65,summaryFunction=summaryFunction,maxClust = 4,mergeCutoff = 3,IQRCutoff = 0.25,medianQuantileCutoff = -1,normalCells=normal_barcodes_manual_qc)
candidate_cnvs_clean<-clusterCNV(initialResultList = candidate_cnvs,medianIQR = candidate_cnvs[[3]],minDiff=1.0) #= 1.5)
final_cnv_list<-annotateCNV4(candidate_cnvs_clean, saveOutput=TRUE,outputSuffix = "RCC87_new_clean_cnv",sdCNV = 0.6,filterResults=TRUE,filterRange=0.4)

cnvs <- read.table(glue('.../RCC87_new_clean_cnv_cnv_scores.csv'), sep = ',', stringsAsFactors =  F, header = T)
cnvs$rowname <- paste0("RCC_87_",cnvs$rowname)
rownames(cnvs) <- cnvs$rowname
cnvs <- cnvs[,-1]
combined_RCC87 = AddMetaData(combined_RCC87, cnvs)
colnames(combined_RCC87)

##############################################################################################################################
############# Post CNV Calling Data Viz ######################################################################################
##############################################################################################################################

Idents(combined_RCC87) <- "celltype"
DimPlot(object = combined_RCC87, label = TRUE, pt.size = 0.3, label.size = 5, label.color = 'black') + NoLegend() +
  theme(text = element_text(size = 14, family = 'Helvetica'),axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.text.x=element_blank()) +
  ylab('UMAP 2') + xlab('UMAP 1')


options(repr.plot.width=10, repr.plot.height=4)
f = FeaturePlot(
  object = combined_RCC99,
  features = c('chr3p'),
  pt.size = 0.8,
)   + theme(text = element_text(size = 14, family = 'Helvetica'),axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.text.x=element_blank(), plot.title = element_text(face = "plain")) + ylab('UMAP 2') + xlab('UMAP 1')

d = DimPlot(combined_RCC99, label = TRUE) + theme(text = element_text(size = 14, family = 'Helvetica'),axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.text.x=element_blank(), plot.title = element_text(face = "plain")) + ylab('UMAP 2') + xlab('UMAP 1')
cowplot::plot_grid(plotlist = list(d,f))

##############################################################################################################################
############# Combine ########################################################################################################
##############################################################################################################################

Idents(combined) <- "dataset"
RCC81_cnvs <- read.table(glue('.../RCC81_new_clean_cnv_cnv_scores.csv'), sep = ',', stringsAsFactors =  F, header = T)
RCC81_cnvs$rowname <- paste0("RCC_81_",RCC81_cnvs$rowname)
rownames(RCC81_cnvs) <- RCC81_cnvs$rowname
RCC81_cnvs <- RCC81_cnvs[,-1]
combined = AddMetaData(combined, RCC81_cnvs)

RCC94_cnvs <- read.table(glue('.../RCC94_new_clean_cnv_cnv_scores.csv'), sep = ',', stringsAsFactors =  F, header = T)
RCC94_cnvs$rowname <- paste0("RCC_94_",RCC94_cnvs$rowname)
rownames(RCC94_cnvs) <- RCC94_cnvs$rowname
RCC94_cnvs <- RCC94_cnvs[,-1]
combined = AddMetaData(combined, RCC94_cnvs)

RCC99_cnvs <- read.table(glue('.../RCC99_new_clean_cnv_cnv_scores.csv'), sep = ',', stringsAsFactors =  F, header = T)
RCC99_cnvs$rowname <- paste0("RCC_99_",RCC99_cnvs$rowname)
rownames(RCC99_cnvs) <- RCC99_cnvs$rowname
RCC99_cnvs <- RCC99_cnvs[,-1]
combined = AddMetaData(combined, RCC99_cnvs)

RCC101_cnvs <- read.table(glue('.../RCC101_new_clean_cnv_cnv_scores.csv'), sep = ',', stringsAsFactors =  F, header = T)
RCC101_cnvs$rowname <- paste0("RCC_101_",RCC101_cnvs$rowname)
rownames(RCC101_cnvs) <- RCC101_cnvs$rowname
RCC101_cnvs <- RCC101_cnvs[,-1]
combined = AddMetaData(combined, RCC101_cnvs)

RCC104_cnvs <- read.table(glue('.../RCC104_new_clean_cnv_cnv_scores.csv'), sep = ',', stringsAsFactors =  F, header = T)
RCC104_cnvs$rowname <- paste0("RCC_104_",RCC104_cnvs$rowname)
rownames(RCC104_cnvs) <- RCC104_cnvs$rowname
RCC104_cnvs <- RCC104_cnvs[,-1]
combined = AddMetaData(combined, RCC104_cnvs)

RCC106_cnvs <- read.table(glue('.../RCC106_new_clean_cnv_cnv_scores.csv'), sep = ',', stringsAsFactors =  F, header = T)
RCC106_cnvs$rowname <- paste0("RCC_106_",RCC106_cnvs$rowname)
rownames(RCC106_cnvs) <- RCC106_cnvs$rowname
RCC106_cnvs <- RCC106_cnvs[,-1]
combined = AddMetaData(combined, RCC106_cnvs)

RCC112_cnvs <- read.table(glue('.../RCC112_new_clean_cnv_cnv_scores.csv'), sep = ',', stringsAsFactors =  F, header = T)
RCC112_cnvs$rowname <- paste0("RCC_112_",RCC112_cnvs$rowname)
rownames(RCC112_cnvs) <- RCC112_cnvs$rowname
RCC112_cnvs <- RCC112_cnvs[,-1]
combined = AddMetaData(combined, RCC112_cnvs)

RCC113_cnvs <- read.table(glue('.../RCC113_new_clean_cnv_cnv_scores.csv'), sep = ',', stringsAsFactors =  F, header = T)
RCC113_cnvs$rowname <- paste0("RCC_113_",RCC113_cnvs$rowname)
rownames(RCC113_cnvs) <- RCC113_cnvs$rowname
RCC113_cnvs <- RCC113_cnvs[,-1]
combined = AddMetaData(combined, RCC113_cnvs)

RCC115_cnvs <- read.table(glue('.../RCC115_new_clean_cnv_cnv_scores.csv'), sep = ',', stringsAsFactors =  F, header = T)
RCC115_cnvs$rowname <- paste0("RCC_115_",RCC115_cnvs$rowname)
rownames(RCC115_cnvs) <- RCC115_cnvs$rowname
RCC115_cnvs <- RCC115_cnvs[,-1]
combined = AddMetaData(combined, RCC115_cnvs)

RCC119_cnvs <- read.table(glue('.../RCC119_new_clean_cnv_cnv_scores.csv'), sep = ',', stringsAsFactors =  F, header = T)
RCC119_cnvs$rowname <- paste0("RCC_119_",RCC119_cnvs$rowname)
rownames(RCC119_cnvs) <- RCC119_cnvs$rowname
RCC119_cnvs <- RCC119_cnvs[,-1]
combined = AddMetaData(combined, RCC119_cnvs)

RCC120_cnvs <- read.table(glue('.../RCC120_new_clean_cnv_cnv_scores.csv'), sep = ',', stringsAsFactors =  F, header = T)
RCC120_cnvs$rowname <- paste0("RCC_120_",RCC120_cnvs$rowname)
rownames(RCC120_cnvs) <- RCC120_cnvs$rowname
RCC120_cnvs <- RCC120_cnvs[,-1]
combined = AddMetaData(combined, RCC120_cnvs)

RCC114_cnvs <- read.table(glue('.../RCC114_new_clean_cnv_cnv_scores.csv'), sep = ',', stringsAsFactors =  F, header = T)
RCC114_cnvs$rowname <- paste0("RCC_114_",RCC114_cnvs$rowname)
rownames(RCC114_cnvs) <- RCC114_cnvs$rowname
RCC114_cnvs <- RCC114_cnvs[,-1]
combined = AddMetaData(combined, RCC114_cnvs)

RCC103_cnvs <- read.table(glue('.../RCC103_new_clean_cnv_cnv_scores.csv'), sep = ',', stringsAsFactors =  F, header = T)
RCC103_cnvs$rowname <- paste0("RCC_103_",RCC103_cnvs$rowname)
rownames(RCC103_cnvs) <- RCC103_cnvs$rowname
RCC103_cnvs <- RCC103_cnvs[,-1]
combined = AddMetaData(combined, RCC103_cnvs)

RCC94_cnvs <- read.table(glue('.../RCC94_new_clean_cnv_cnv_scores.csv'), sep = ',', stringsAsFactors =  F, header = T)
RCC94_cnvs$rowname <- paste0("RCC_94_",RCC94_cnvs$rowname)
rownames(RCC94_cnvs) <- RCC94_cnvs$rowname
RCC94_cnvs <- RCC94_cnvs[,-1]
combined = AddMetaData(combined, RCC94_cnvs)

##########################################################################################
############# Run sample by sample ############### 84, 86 and 96 failed,  ################
##########################################################################################

VlnPlot(combined, features="chr3p", split.by="celltype")
combined_filt <- subset(combined, idents=c("RCC101","RCC103","RCC104","RCC106","RCC112","RCC113","RCC114","RCC115","RCC119","RCC120","RCC81","RCC94","RCC99"))

FeaturePlot(
  object = combined_filt,
  features = c('chr3q'),
  pt.size = 0.5, min.cutoff=0, max.cutoff = 2, cols=c("white","blue")
)  + theme(text = element_text(size = 14, family = 'Helvetica'),axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.text.x=element_blank(), plot.title = element_text(face = "plain")) + ylab('UMAP 2') + xlab('UMAP 1')

