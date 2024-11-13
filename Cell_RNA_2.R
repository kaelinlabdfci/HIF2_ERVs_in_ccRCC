######################################################
################_scTE Diff_########################
######################################################
library(presto)
setwd(".../...")
de_genes <- FindMarkers(
  object = combined,
  ident.1 = "Cancer",
  ident.2 = "Immune",
  test.use = 'wilcox', only.pos=F, min.pct=0.01, logfc.threshold=0)

de_genes$logpadj <- -log10(de_genes$p_val_adj)
de_genes$significant <- ifelse((de_genes$p_val_adj <= 0.05 & de_genes$avg_log2FC > 1 | de_genes$p_val_adj <= 0.05 & de_genes$avg_log2FC < -1), "Sig", "Not Sig")
de_genes$change <- ifelse((de_genes$significant=="Sig" & de_genes$avg_log2FC > 1), "Up", ifelse((de_genes$significant=="Sig" & de_genes$avg_log2FC < -1), "Down", "Unchanged"))

sum(de_genes$change=="Up") #This is 4552
sum(de_genes$change=="Down") #This is 3431
write.csv(de_genes, file="raw_de_genes_and_erv_07252024.csv")  

######################################################
################_Focus on ERVs + gene controls########
######################################################
library(ggrepel)
#Plot All ERVs and highlight HIF regulated ERVs. Spike is a column for CA9 and PTPRC
de_ervs <- read.csv("raw_de_genes_and_erv_07252024.csv")
de_ervs <- na.omit(de_ervs)
qj <- read.csv("81_ervs.csv")
colnames(qj) <- "X"
sum(de_ervs$change=="Up") #This is 86
sum(de_ervs$change=="Down") #This is 42
de_ervs$hif <- ifelse(de_ervs$ERV %in% qj$X, "Y", "N")

write.csv(de_ervs, file="detected_hif_ervs.csv")
de_ervs <- read.csv("07282024_ERV_scRNA.csv")
ggplot(data = de_ervs, aes(x = avg_log2FC, y = logpadj)) + geom_point(alpha = 0.7, shape = 16, size = 0.5) + 
  labs(title = "Differential Expression of ERVs in Tumor and CD45+ Cells",
       x = "log2(fold change)",
       y = "-log10(adjusted P-value)",
       colour = "HIF2-regulated ERVs") + theme_bw() + 
  theme(panel.border = element_rect(colour = "black", fill = NA, size= 2), 
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank()) + NoLegend() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") + 
  geom_vline(xintercept = c(-1,1), linetype = "dashed") + geom_point(aes(colour=hif), size=0.8) + scale_color_manual(values = c("grey", "blue")) +
  xlim(-10,10) + geom_label_repel(data=dplyr::filter(de_ervs, imp=="Y"), aes(label=ERV,
                                                                                     box.padding = 1,
                                                                                     nudge_x = 2,
                                                                                     segment.ncp = 2,
                                                                                     segment.angle = 1)) + geom_label_repel(data=dplyr::filter(de_ervs, spike=="Y"), aes(label=ERV),
                                                                                                                            box.padding = 1,
                                                                                                                            nudge_x = 0.25,
                                                                                                                            nudge_y = 1,
                                                                                                                            segment.curvature = -0.05,
                                                                                                                            segment.ncp = 3,
                                                                                                                            segment.angle = 1,color = "red",     # text color
                                                                                                                            bg.color = "white", # shadow color
                                                                                                                            bg.r = 0.1 ) 
#Only plot the HIF regulated ERVs
de_ervs_hif <- de_ervs %>%
  filter(hif=="Y" | spike=="Y")
table(de_ervs_hif$hif=="Y" & de_ervs_hif$change=="Up") #This is 16
table(de_ervs_hif$hif=="Y" & de_ervs_hif$change=="Down") #This is 2
table(de_ervs_hif$hif=="Y" & de_ervs_hif$change=="Unchanged") #This 5
write.csv(de_ervs_hif,"detected_hif_ervs.csv")

ggplot(data = de_ervs_hif, aes(x = avg_log2FC, y = logpadj)) + geom_point(alpha = 0.7, shape = 16, size = 0.5) + 
  labs(title = "Differential Expression of ERVs in Tumor and CD45+ Cells",
       x = "log2(fold change)",
       y = "-log10(adjusted P-value)",
       colour = "HIF2-regulated ERVs") + theme_bw() + 
  theme(panel.border = element_rect(colour = "black", fill = NA, size= 2), 
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank()) + NoLegend() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") + 
  geom_vline(xintercept = c(-1,1), linetype = "dashed") + geom_point(aes(colour=hif), size=0.8) + scale_color_manual(values = c("grey","blue")) +
  xlim(-10,10) + geom_label_repel(data=dplyr::filter(de_ervs, imp=="Y"), aes(label=ERV,
                                                                             box.padding = 1,
                                                                             nudge_x = 2,
                                                                             segment.ncp = 2,
                                                                             segment.angle = 1)) + geom_label_repel(data=dplyr::filter(de_ervs, spike=="Y"), aes(label=ERV),
                                                                                                                    box.padding = 1,
                                                                                                                    nudge_x = 0.25,
                                                                                                                    nudge_y = 1,
                                                                                                                    segment.curvature = -0.05,
                                                                                                                    segment.ncp = 3,
                                                                                                                    segment.angle = 1,color = "red",     # text color
                                                                                                                    bg.color = "white", # shadow color
                                                                                                                    bg.r = 0.1 ) 
