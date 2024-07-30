##################################################################
##################################################################
############. Computing Differential Accessibility  ##############
##################################################################
##################################################################

library(ggrepel)
library(tidyverse)

DefaultAssay(combined_filt) <- 'ATAC'
Idents(combined_filt) <- "celltype"
DimPlot(object = combined_filt, label = TRUE, pt.size = 0.3, label.size = 5, label.color = 'black') + NoLegend() +
  theme(text = element_text(size = 14, family = 'Helvetica'),axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.text.x=element_blank()) +
  ylab('UMAP 2') + xlab('UMAP 1')

filtered <- subset(combined_filt, idents=c("Cancer","Immune"))

da_peaks <- FindMarkers(
  object = filtered,
  ident.1 = "Cancer",
  ident.2 = "Immune",
  test.use = 'LR',
  latent.vars = 'peak_region_fragments', only.pos=F, logfc.threshold=0)

setwd(".../DiffATAC")
da_peaks$logpadj <- -log10(da_peaks$p_val_adj)
da_peaks$significant <- ifelse((da_peaks$p_val_adj <= 0.05 & da_peaks$avg_log2FC > 1 | da_peaks$p_val_adj <= 0.05 & da_peaks$avg_log2FC < -1), "Sig", "Not Sig")
da_peaks$change <- ifelse((da_peaks$significant=="Sig" & da_peaks$avg_log2FC > 1), "Up", ifelse((da_peaks$significant=="Sig" & da_peaks$avg_log2FC < -1), "Down", "Unchanged"))
sum(da_peaks$change=="Up") #This is 35359  
sum(da_peaks$change=="Down") #This is 35585  
write.csv(da_peaks, file="all_da_peaks_07242024.csv")                  

#Rough first plot
ggplot(data = da_peaks,
       aes(x = avg_log2FC,
           y = logpadj)) + geom_point(alpha = 0.2, 
           shape = 16,
           size = 1) + labs(title = "Differential accessibility betnween Tumor and Immune",
           x = "log2(fold change)",
           y = "-log10(adjusted P-value)",
          colour = "Expression \nchange") + theme_bw() + theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank()) + geom_hline(yintercept = -log10(0.05),
          linetype = "dashed") + 
  geom_vline(xintercept = 0,
             linetype = "dashed")

#Used bedtools intersect to get coordinates for differential peaks at ERVs using computing cluster

library(dplyr)
setwd(".../DiffATAC")
reference <- read.csv("0724_all_intersect_erv_edit.csv")
reference=reference$coordinates

da.new <- da_peaks %>% dplyr::filter(row.names(da_peaks) %in% reference)
da.new$significant <- ifelse((da.new$p_val_adj <= 0.05 & da.new$avg_log2FC > 1 | da.new$p_val_adj <= 0.05 & da.new$avg_log2FC < -1), "Sig", "Not Sig")
da.new$change <- ifelse((da.new$significant=="Sig" & da.new$avg_log2FC > 1), "Up", ifelse((da.new$significant=="Sig" & da.new$avg_log2FC < -1), "Down", "Unchanged"))
sum(da.new$change=="Up") #This is 322
sum(da.new$change=="Down") #This is 267

##################################################################
##################################################################
#######. Make GRanges object out of ERV bed & Intersecting Bed  ##
##################################################################
##################################################################

write.csv(da.new, file="07242024_Diff_ERV_ATAC.csv")

#Edit columns to have 'Chr', 'Start', 'Stop'

atac_erv <- read.csv("07242024_Diff_ERV_ATAC.csv")

all_erv <- read.delim("hervquant_hg38_3kb.bed.txt", header=TRUE)
all_erv <- makeGRangesFromDataFrame(all_erv,
                                    keep.extra.columns=TRUE,
                                    ignore.strand=FALSE,
                                    seqinfo=NULL,
                                    seqnames.field=c("seqnames", "seqname",
                                                     "chromosome", "chrom",
                                                     "Chr", "chromosome_name",
                                                     "seqid"),
                                    start.field="Start",
                                    end.field=c("end", "Stop"),
                                    strand.field="strand")
all_erv

atac <- makeGRangesFromDataFrame(atac_erv,
                                 keep.extra.columns=FALSE,
                                 ignore.strand=FALSE,
                                 seqinfo=NULL,
                                 seqnames.field=c("seqnames", "seqname",
                                                  "chromosome", "chrom",
                                                  "Chr", "chromosome_name",
                                                  "seqid"),
                                 start.field="Start",
                                 end.field=c("end", "Stop"),
                                 strand.field="strand")

overlaps <- findOverlaps(query = atac, subject = all_erv)
queryHits(overlaps)
subjectHits(overlaps)
mcols(atac)$ERV <- NA
mcols(atac)$ERV[queryHits(overlaps)] <- mcols(all_erv)$Annotation[subjectHits(overlaps)]
atac
atac <- as.data.frame(atac)
da.new$ERV <- atac$ERV
write.csv(da.new, file="07242024_ATAC_ERV_Master_List.csv")


#######################################
#######################################
############. Volcano Plot  ###########
#######################################
#######################################

da.new <- read.csv("07242024_ATAC_ERV_Master_List.csv")
da.new$logpadj <- -log10(da.new$p_val_adj)
qj <- read.csv("81_ervs.csv")
colnames(qj) <- "X"

da.new$imp <- ifelse(da.new$ERV %in% qj$X, "Y", "N")

ggplot(data = da.new,
       aes(x = avg_log2FC, y = logpadj)) + geom_point(alpha = 0.7, shape = 16, size = 0.5) + 
  labs(title = "Differential accessibility at ERVs +/- 3 kb between Tumor and CD45+ Cells",
       x = "log2(fold change)",
       y = "-log10(adjusted P-value)",
       colour = "HIF2-regulated ERVs") + 
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill = NA, size= 1), 
                     panel.grid.minor = element_blank(), 
                     panel.grid.major = element_blank()) + NoLegend() +
  geom_hline(yintercept = -log10(0.05),                                                                                                                                                  linetype = "dashed") + 
  geom_vline(xintercept = 0, linetype = "dashed") + geom_point(aes(colour=imp), size=0.8, ) + scale_color_manual(values = c("grey", "blue")) +
  xlim(-5,5) + geom_label_repel(data=dplyr::filter(da.new, imp=="Y"), aes(label=ERV),
                                box.padding = 1,
                                nudge_x = 0.1,
                                nudge_y = 1,
                                segment.curvature = -0.05,
                                segment.ncp = 3,
                                segment.angle = 20) 
sum(da.new$change=="Up") #This is 323
sum(da.new$change=="Down") #This is 268

##################################################################
##################################################################
############. Add CA9 and PTPRC  #################################
##################################################################
##################################################################

#BiocManager::install("ChIPseeker")
#BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
#BiocManager::install("org.Hs.eg.db")

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

all <- read.csv("GR_all_da_peaks_07242024.csv")
wrong <- which(all$Start > all$Stop)
new_end   <- all[wrong,]$Start
new_start <- all[wrong,]$Stop
all[wrong,]$Start <- new_start
all[wrong,]$Stop   <- new_end


gr <- makeGRangesFromDataFrame(all,
                               keep.extra.columns=FALSE,
                               ignore.strand=FALSE,
                               seqinfo=NULL,
                               seqnames.field=c("seqnames", "seqname",
                                                "chromosome", "chrom",
                                                "Chr", "chromosome_name",
                                                "seqid"),
                               start.field="Start",
                               end.field=c("end", "Stop"),
                               strand.field="strand")

peakAnno <- annotatePeak(gr, tssRegion=c(-3000,3000), TxDb=txdb, annoDb="org.Hs.eg.db")
Anno <- as.data.frame(peakAnno)
symbol <- Anno$SYMBOL
annotation <- Anno$annotation
all$SYMBOL <- symbol
all$annotation <- annotation

write.csv(all, "All_peaks_annotated_07242024.csv")

#################################################
#################################################
#######GGplot after Spike in of CA9 and CD45#####
#################################################
#################################################

da.new <- read.csv("07192024_ATAC_ERV_Master_List_with_spike.csv")
qj <- read.csv("81_ervs.csv")
colnames(qj) <- "X"
da.new$hif <- ifelse(da.new$ERV %in% qj$X, "Y", "N")
write.csv(da.new, file="07192024_ATAC_ERV_Master_List_with_spike.csv")
ggplot(data = da.new,
       aes(x = avg_log2FC, y = logpadj)) + geom_point(alpha = 0.7, shape = 16, size = 1.75) + 
  labs(title = "Differential accessibilitys",
       x = "log2(fold change)",
       y = "-log10(adjusted P-value)",
       colour = "HIF2-regulated ERVs") + 
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill = NA, size= 1), 
                     panel.grid.minor = element_blank(), 
                     panel.grid.major = element_blank()) + NoLegend() +
  geom_hline(yintercept = -log10(0.05),                                                                                                                                                  linetype = "dashed") + 
  geom_vline(xintercept = c(-1,1), linetype = "dashed") + geom_point(aes(colour=hif), size=0.8, ) + scale_color_manual(values = c("grey", "blue")) +
  xlim(-8,8)  + geom_label_repel(data=dplyr::filter(da.new, imp=="Y"), aes(label=ERV),
                                box.padding = 0.5,
                                nudge_x = 0.1,
                                nudge_y = 0.1,
                                segment.curvature = -0.1,
                                segment.ncp = 1,
                                segment.angle = 5) + geom_label_repel(data=dplyr::filter(da.new, spike=="Y"), aes(label=ERV),
                                                                       box.padding = 1,
                                                                       nudge_x = 0.1,
                                                                       nudge_y = 1,
                                                                       segment.curvature = -0.05,
                                                                       segment.ncp = 3,
                                                                       segment.angle = 1,color = "red",     # text color
                                                                       bg.color = "white", # shadow color
                                                                       bg.r = 0.1  )


da.new_hif <- da.new %>%
  filter(hif=="Y" | spike=="Y")
da.new_hif_unc <- da.new_hif %>%
  filter(change=="Unchanged")
unique(da.new_hif_unc$ERV) #10 Unchanged
da.new_hif_up <- da.new_hif %>%
  filter(change=="Up")
unique(da.new_hif_up$ERV)
da.new_hif_dn <- da.new_hif %>% #26 Up
  filter(change=="Down")
unique(da.new_hif_dn$ERV) #5 Down

#18 undetected






