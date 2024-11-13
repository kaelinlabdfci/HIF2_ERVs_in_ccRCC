#Plot for scRNA-seq correlated to ATAC up, down, unc | Figures 5E to G
setwd("...")

library(dplyr)
library(tidyverse)
library(ggrepel)

atac <- read.csv("07192024_ATAC_ERV_Master_List_with_spike.csv")
rna <-read.csv("detected_hif_ervs.csv")

atac_down <- filter(atac, change=="Down")
atac_up <- filter(atac, change=="Up")
atac_unchanged <- filter(atac, change=="Unchanged")

atac_down <- atac_down$ERV
atac_up <- atac_up$ERV
atac_unchanged <- atac_unchanged$ERV
atac_all <- atac$ERV

rna_down <- rna[rna$ERV %in% atac_down, ]
rna_up <- rna[rna$ERV %in% atac_up, ]
rna_unchanged <- rna[rna$ERV %in% atac_unchanged, ]
rna <- rna[rna$ERV %in% atac_all, ]

write.csv(rna_down, file="0725_rna_atac_down.csv")
write.csv(rna_up, file="0725_rna_atac_up.csv")
write.csv(rna_unchanged, file="0725_rna_atac_unc.csv")
write.csv(rna, file="0725_rna_atac_all.csv")

#XvY scatter plot, Figure 5E
rna$logrna <- rna$avg_log2FC
atac$logatac <- atac$avg_log2FC

combined <- merge(x = rna, y = atac, by = "ERV", all.y = TRUE)
combined <- na.omit(combined)
combined2 <- filter(combined, significant.x=="Sig" & significant.y=="Sig")

colnames(qj) <- "X"
combined2$hif <- ifelse(combined2$ERV %in% qj$X, "Y", "N")

library(ggrepel)

write.csv(combined2, file="Correlation_Plot.csv")
combined2 <- read.csv("Correlation_Plot.csv")

combined2 <- filter(combined2, combined2$hif=="Y" | combined2$spike.x=="Y")

ggplot(data = combined2,
       aes(x = logrna, y = logatac)) + geom_point(alpha = 0.7, shape = 16, size = 1.75) + 
  labs(title = "scATAC-seq vs scRNA-seq",
       x = "Log2FC scRNA-seq (Malignant/Immune)",
       y = "Log2FC scATAC-seq (Malignant/Immune)",
       colour = "HIF2-regulated ERVs") + theme_bw() + theme(panel.border = element_rect(colour = "black", fill = NA, size= 1.5), 
                     panel.grid.minor = element_blank(), 
                     panel.grid.major = element_blank()) + geom_hline(yintercept = c(-1,1), linetype='dashed',                                                                                                                                                  linetype = "dashed") + 
  geom_vline(xintercept = c(-1,1), linetype='dashed') + geom_point(aes(colour=hif), size=0.8, ) + scale_color_manual(values = c("grey", "blue")) +
  xlim(-7,7) + geom_label_repel(data=dplyr::filter(combined2, imp.x=="Y"), aes(label=ERV),
                                box.padding = 0.1,
                                nudge_x = 0.5,
                                nudge_y = 1,
                                segment.curvature = -0.1,
                                segment.ncp = 2,
                                segment.angle = 20) + geom_label_repel(data=dplyr::filter(combined2, spike.x=="Y"), aes(label=ERV),
                                                                       box.padding = 5,
                                                                       nudge_x = 0.1,
                                                                       nudge_y = 0.5,
                                                                       segment.curvature = -0.05,
                                                                       segment.ncp = 3,
                                                                       segment.angle = 1,color = "red",     # text color
                                                                       bg.color = "white", # shadow color
                                                                       bg.r = 0.1  ) + NoLegend()

colnames(qj) <- "X"
combined2$hif <- ifelse(combined2$ERV %in% qj$X, "Y", "N")

