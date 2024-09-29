setwd("/Volumes/VJ_Data/Cell Manuscript_single cell ERV work/scRNA-seq_scTE")
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(tidyverse)
#remotes::install_github("pmbio/MuDataSeurat")
library(MuDataSeurat)
library(scDblFinder)

######################################################
################_Load_H5AD_########################
######################################################

RCC99 <- ReadH5AD("RCC99.h5ad")
RCC101 <- ReadH5AD("RCC101.h5ad")
RCC103 <- ReadH5AD("RCC103.h5ad")
RCC104 <- ReadH5AD("RCC104.h5ad")
RCC106 <- ReadH5AD("RCC106.h5ad")
RCC112 <- ReadH5AD("RCC112.h5ad")
RCC113 <- ReadH5AD("RCC113.h5ad")
RCC114 <- ReadH5AD("RCC114.h5ad")
RCC81 <- ReadH5AD("RCC81.h5ad")
RCC84 <- ReadH5AD("RCC84.h5ad")
RCC119 <- ReadH5AD("RCC119.h5ad")
RCC120 <- ReadH5AD("RCC120.h5ad")
RCC115 <- ReadH5AD("RCC115.h5ad")
RCC86 <- ReadH5AD("RCC86.h5ad")
RCC94 <- ReadH5AD("RCC94.h5ad")
RCC96 <- ReadH5AD("RCC96.h5ad")

################################################
################_Harmony_#######################
################################################
#Remove RCC100, RCC116, RCC87

combined <- merge(RCC101, y=c(RCC103, RCC104, RCC106, RCC112, RCC113, RCC114, RCC81, RCC84, RCC119, RCC120, RCC115, RCC86, RCC94, RCC96, RCC99), add.cell.ids=c("RCC101", "RCC103", "RCC104", "RCC106", "RCC112", "RCC113", "RCC114","RCC81", "RCC84", "RCC119", "RCC120", "RCC115", "RCC86", "RCC94", "RCC96", "RCC99"), project="Combined ccRCC")  
combined$column <- NA
combined$column[grepl("RCC100",colnames(combined))] <- "RCC100"
combined$column[grepl("RCC101",colnames(combined))] <- "RCC101"
combined$column[grepl("RCC103",colnames(combined))] <- "RCC103"
combined$column[grepl("RCC104",colnames(combined))] <- "RCC104"
combined$column[grepl("RCC106",colnames(combined))] <- "RCC106"
combined$column[grepl("RCC112",colnames(combined))] <- "RCC112"
combined$column[grepl("RCC113",colnames(combined))] <- "RCC113"
combined$column[grepl("RCC114",colnames(combined))] <- "RCC114"
combined$column[grepl("RCC116",colnames(combined))] <- "RCC116"
combined$column[grepl("RCC81",colnames(combined))] <- "RCC81"
combined$column[grepl("RCC84",colnames(combined))] <- "RCC84"
combined$column[grepl("RCC119",colnames(combined))] <- "RCC119"
combined$column[grepl("RCC120",colnames(combined))] <- "RCC120"
combined$column[grepl("RCC115",colnames(combined))] <- "RCC115"
combined$column[grepl("RCC86",colnames(combined))] <- "RCC86"
combined$column[grepl("RCC87",colnames(combined))] <- "RCC87"
combined$column[grepl("RCC94",colnames(combined))] <- "RCC94"
combined$column[grepl("RCC96",colnames(combined))] <- "RCC96"
combined$column[grepl("RCC99",colnames(combined))] <- "RCC99"

#a <- rownames(combined)
#a <- as.data.frame(a)
#write.csv(a, file="genes_and_ervs.csv")

peak_assay = GetAssayData(combined, layer="counts")
sce = scDblFinder(peak_assay)
combined$doublet_class = sce@colData$scDblFinder.class
combined$doublet_class = factor(combined$doublet_class, levels = c('singlet', 'doublet'))
combined$doublet_score = sce@colData$scDblFinder.score
rm(sce)
rm(peak_assay)

g <- ggplot(combined[[]], aes_string(x = 'doublet_score', color = 'doublet_class', fill = 'doublet_class')) + geom_histogram()
g + scale_color_manual(name = 'Final classification', values=c("#0571b0", "#ca0020"))+
  scale_fill_manual(name = 'Final classification',values=c("#0571b0", "#ca0020")) + theme_classic() + ylab('Number of cells') + xlab('Probayounglity of cell being a\ndoublet (score)')+ 
  theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), text = element_text(size = 16, family = 'Helvetica'), legend.title=element_text(size=14 ,family = 'Helvetica'), legend.position = 'right')

combined = subset(
  x = combined,
  subset = doublet_class == 'singlet'
)

combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^MT-")
VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=F)

combined <- subset(combined, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 10)

combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
combined <- ScaleData(combined, assay='RNA', features=rownames(combined))
Idents(combined) <- "column"
levels(combined)
combined <- RunPCA(combined, features = VariableFeatures(object = combined))
DimPlot(combined, reduction = "pca")

library(harmony)
combined <- combined %>% 
  RunHarmony("column", plot_convergence = TRUE)

harmony_embeddings <- Embeddings(combined, 'harmony')
harmony_embeddings[1:5, 1:5]
DimPlot(object = combined, reduction = "harmony", group.by = "column")

combined <- combined %>% 
  RunUMAP(reduction = "harmony", dims = 1:10) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:10) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

Idents(combined) <- "column"
Idents(combined) <- "seurat_clusters"
DimPlot(combined, reduction = "umap",raster=FALSE, label=F) +
  theme(text = element_text(size = 14, family = 'Helvetica'),axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.text.x=element_blank()) +
  ylab('UMAP 2') + xlab('UMAP 1') + NoLegend()
DotPlot(combined, features = c("CA9","KRT18","PECAM1","MCAM","PTPRC","NKG7"), group.by = "seurat_clusters")
VlnPlot(combined, features=c("565"), group.by="celltype", pt.size=0.5)


#clusters: Immune, Cancer, Immune, Immune, Immune, Immune, Endothelial, Endothelial, Interstitial, Cancer, Immune, Immune, Immune, Endothelial, Immune, Endothelial, Interstitial, Immune
Idents(combined) <- "celltype"
DotPlot(combined, features = c("CA9","KRT18","PECAM1","MCAM","PTPRC","NKG7"))
DotPlot(combined, features = c("3254", "673", "876", "838", "5875","3797","504","505"))
FeaturePlot(combined, features = c("3547"),order=F,raster=F, pt.size=0.1) + NoLegend() +
  theme(text = element_text(size = 14, family = 'Helvetica'),axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.text.x=element_blank()) +
  ylab('UMAP 2') + xlab('UMAP 1')

new.cluster.ids <- c("Immune", "Cancer", "Immune", "Immune", "Immune", "Immune", "Endothelial", "Immune", "Interstitial", "Cancer", "Immune", "Unknown", "Immune", "Endothelial", "Immune", "Immune", "Interstitial", "Immune")
names(new.cluster.ids) <- levels(combined)
combined <- RenameIdents(combined, new.cluster.ids)
combined$celltype <- Idents(combined)

DimPlot(combined, reduction = "umap",raster=FALSE, label=T, group.by = "celltype")  +
  theme(text = element_text(size = 14, family = 'Helvetica'),axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.text.x=element_blank()) +
  ylab('UMAP 2') + xlab('UMAP 1')

Idents(combined) <- "celltype"
FeaturePlot(combined, features = c("4444", "6078", "2334", "4784", "5875","3797","504","505"),order=TRUE,raster=F, pt.size=0.1) +
  theme(text = element_text(size = 14, family = 'Helvetica'),axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.text.x=element_blank()) +
  ylab('UMAP 2') + xlab('UMAP 1')

rm(combined)

setwd("/Volumes/VJ_Data/Cell Manuscript_DiffAcc")
atac <- read.csv("Ann_Diff_ERV_ATAC_spike_in.csv")
sum(atac$qj=="Y" & atac$change=="Down")
sum(atac$qj=="Y" & atac$change=="Up")

cancer <- subset(combined, idents="Cancer")
immune <- subset(combined, idents="Immune")
DoHeatmap(object = cancer, features = c("4444", "6078", "PTPRC", "CA9", "5875","3797","2256","505"), slot="scale.data", group.by="column", draw.lines=F)+ scale_fill_gradientn(colors = c("blue", "white", "red"))
DoHeatmap(object = immune, features = c("4444", "6078", "PTPRC", "CA9", "5875","3797","2256","505"), slot="data", group.by="column", draw.lines=F)+ scale_fill_gradientn(colors = c("blue", "white", "red"))

filtered <- subset(combined, idents=c("Cancer","Immune"))
VlnPlot(filtered, features=c("PTPRC"), group.by="celltype", pt.size=0.03) + geom_violin(aes(fill = "gray", alpha = 0.5)) + NoLegend()


combined <- ScaleData(combined, assay='RNA', features=rownames(combined))
a <- rownames(combined) 

