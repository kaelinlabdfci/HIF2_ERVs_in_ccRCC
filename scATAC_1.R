####scATAC-seq
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)
library(GenomicRanges)
library(future)
library(harmony)
library(devtools)
library(scDblFinder)
library(biovizBase)
library(dplyr)
library(tidyverse)

##################################################################
##################################################################
############. Remove RCC87, RCC116, RCC100 as VHL wild type ######
##################################################################
##################################################################

peaks.RCC81 <- read.table(
  file = ".../scATAC-seq/RCC81/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.RCC84 <- read.table(
  file = ".../scATAC-seq/RCC84/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.RCC86 <- read.table(
  file = ".../scATAC-seq/RCC86/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.RCC94 <- read.table(
  file = ".../scATAC-seq/RCC94/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.RCC96 <- read.table(
  file = ".../scATAC-seq/RCC96/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.RCC99 <- read.table(
  file = ".../scATAC-seq/RCC99/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.RCC101 <- read.table(
  file = ".../scATAC-seq/RCC101/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.RCC103 <- read.table(
  file = ".../scATAC-seq/RCC103/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.RCC104 <- read.table(
  file = ".../scATAC-seq/RCC104/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.RCC106 <- read.table(
  file = ".../scATAC-seq/RCC106/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.RCC112 <- read.table(
  file = ".../scATAC-seq/RCC112/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.RCC113 <- read.table(
  file = ".../scATAC-seq/RCC113/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.RCC114 <- read.table(
  file = ".../scATAC-seq/RCC114/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.RCC115 <- read.table(
  file = ".../scATAC-seq/RCC115/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.RCC119 <- read.table(
  file = ".../scATAC-seq/RCC119/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.RCC120 <- read.table(
  file = ".../scATAC-seq/RCC120/peaks.bed",
  col.names = c("chr", "start", "end")
)

gr.RCC81 <- makeGRangesFromDataFrame(peaks.RCC81)
gr.RCC84 <- makeGRangesFromDataFrame(peaks.RCC84)
gr.RCC86 <- makeGRangesFromDataFrame(peaks.RCC86)
gr.RCC94 <- makeGRangesFromDataFrame(peaks.RCC94)
gr.RCC96 <- makeGRangesFromDataFrame(peaks.RCC96)
gr.RCC99 <- makeGRangesFromDataFrame(peaks.RCC99)
gr.RCC101 <- makeGRangesFromDataFrame(peaks.RCC101)
gr.RCC103 <- makeGRangesFromDataFrame(peaks.RCC103)
gr.RCC104 <- makeGRangesFromDataFrame(peaks.RCC104)
gr.RCC106 <- makeGRangesFromDataFrame(peaks.RCC106)
gr.RCC112 <- makeGRangesFromDataFrame(peaks.RCC112)
gr.RCC113 <- makeGRangesFromDataFrame(peaks.RCC113)
gr.RCC114 <- makeGRangesFromDataFrame(peaks.RCC114)
gr.RCC115 <- makeGRangesFromDataFrame(peaks.RCC115)
gr.RCC119 <- makeGRangesFromDataFrame(peaks.RCC119)
gr.RCC120 <- makeGRangesFromDataFrame(peaks.RCC120)

combined.peaks <- reduce(x = c(gr.RCC81, gr.RCC84, gr.RCC86,gr.RCC94,gr.RCC96,gr.RCC99, gr.RCC101,gr.RCC103,gr.RCC104,gr.RCC106,gr.RCC112,gr.RCC113,gr.RCC114,gr.RCC115,gr.RCC119,gr.RCC120))
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

md.RCC81 <- read.table(
  file = ".../scATAC-seq/RCC81/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]
md.RCC84 <- read.table(
  file = ".../scATAC-seq/RCC84/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]
md.RCC86 <- read.table(
  file = ".../scATAC-seq/RCC86/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] 
md.RCC94 <- read.table(
  file = ".../scATAC-seq/RCC94/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]
md.RCC96 <- read.table(
  file = ".../scATAC-seq/RCC96/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] 
md.RCC99 <- read.table(
  file = ".../scATAC-seq/RCC99/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]
md.RCC101 <- read.table(
  file = ".../scATAC-seq/RCC101/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] 
md.RCC103 <- read.table(
  file = ".../scATAC-seq/RCC103/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]
md.RCC104 <- read.table(
  file = ".../scATAC-seq/RCC104/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]
md.RCC106 <- read.table(
  file = ".../scATAC-seq/RCC106/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] 
md.RCC112 <- read.table(
  file = ".../scATAC-seq/RCC112/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] 
md.RCC113 <- read.table(
  file = ".../scATAC-seq/RCC113/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]
md.RCC114 <- read.table(
  file = ".../scATAC-seq/RCC114/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]
md.RCC115 <- read.table(
  file = ".../scATAC-seq/RCC115/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] 
md.RCC119 <- read.table(
  file = ".../scATAC-seq/RCC119/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] 
md.RCC120 <- read.table(
  file = ".../scATAC-seq/RCC120/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] 

md.RCC81 <- md.RCC81[md.RCC81$is__cell_barcode > 0.5, ]
md.RCC84 <- md.RCC84[md.RCC84$is__cell_barcode > 0.5, ]
md.RCC86 <- md.RCC86[md.RCC86$is__cell_barcode > 0.5, ]
md.RCC94 <- md.RCC94[md.RCC94$is__cell_barcode > 0.5, ]
md.RCC96 <- md.RCC96[md.RCC96$is__cell_barcode > 0.5, ]
md.RCC99 <- md.RCC99[md.RCC99$is__cell_barcode > 0.5, ]
md.RCC101 <- md.RCC101[md.RCC101$is__cell_barcode > 0.5, ]
md.RCC103 <- md.RCC103[md.RCC103$is__cell_barcode > 0.5, ]
md.RCC104 <- md.RCC104[md.RCC104$is__cell_barcode > 0.5, ]
md.RCC106 <- md.RCC106[md.RCC106$is__cell_barcode > 0.5, ]
md.RCC112 <- md.RCC112[md.RCC112$is__cell_barcode > 0.5, ]
md.RCC113 <- md.RCC113[md.RCC113$is__cell_barcode > 0.5, ]
md.RCC114 <- md.RCC114[md.RCC114$is__cell_barcode > 0.5, ]
md.RCC115 <- md.RCC115[md.RCC115$is__cell_barcode > 0.5, ]
md.RCC119 <- md.RCC119[md.RCC119$is__cell_barcode > 0.5, ]
md.RCC120 <- md.RCC120[md.RCC120$is__cell_barcode > 0.5, ]

frags.RCC81 <- CreateFragmentObject(
  path = ".../scATAC-seq/RCC81/fragments.tsv.gz",
  cells = rownames(md.RCC81)
)
frags.RCC84 <- CreateFragmentObject(
  path = ".../scATAC-seq/RCC84/fragments.tsv.gz",
  cells = rownames(md.RCC84)
)
frags.RCC86 <- CreateFragmentObject(
  path = ".../scATAC-seq/RCC86/fragments.tsv.gz",
  cells = rownames(md.RCC86)
)
frags.RCC94 <- CreateFragmentObject(
  path = ".../scATAC-seq/RCC94/fragments.tsv.gz",
  cells = rownames(md.RCC94)
)
frags.RCC96 <- CreateFragmentObject(
  path = ".../scATAC-seq/RCC96/fragments.tsv.gz",
  cells = rownames(md.RCC96)
)
frags.RCC99 <- CreateFragmentObject(
  path = ".../scATAC-seq/RCC99/fragments.tsv.gz",
  cells = rownames(md.RCC99)
)
frags.RCC101 <- CreateFragmentObject(
  path = ".../scATAC-seq/RCC101/fragments.tsv.gz",
  cells = rownames(md.RCC101)
)
frags.RCC103 <- CreateFragmentObject(
  path = ".../scATAC-seq/RCC103/fragments.tsv.gz",
  cells = rownames(md.RCC103)
)
frags.RCC104 <- CreateFragmentObject(
  path = ".../scATAC-seq/RCC104/fragments.tsv.gz",
  cells = rownames(md.RCC104)
)
frags.RCC106 <- CreateFragmentObject(
  path = ".../scATAC-seq/RCC106/fragments.tsv.gz",
  cells = rownames(md.RCC106)
)
frags.RCC112 <- CreateFragmentObject(
  path = ".../scATAC-seq/RCC112/fragments.tsv.gz",
  cells = rownames(md.RCC112)
)
frags.RCC113 <- CreateFragmentObject(
  path = ".../scATAC-seq/RCC113/fragments.tsv.gz",
  cells = rownames(md.RCC113)
)
frags.RCC114 <- CreateFragmentObject(
  path = ".../scATAC-seq/RCC114/fragments.tsv.gz",
  cells = rownames(md.RCC114)
)
frags.RCC115 <- CreateFragmentObject(
  path = ".../scATAC-seq/RCC115/fragments.tsv.gz",
  cells = rownames(md.RCC115)
)
frags.RCC119 <- CreateFragmentObject(
  path = ".../scATAC-seq/RCC119/fragments.tsv.gz",
  cells = rownames(md.RCC119)
)
frags.RCC120 <- CreateFragmentObject(
  path = ".../scATAC-seq/RCC120/fragments.tsv.gz",
  cells = rownames(md.RCC120)
)

##### Next Step

RCC81.counts <- FeatureMatrix(
  fragments = frags.RCC81,
  features = combined.peaks,
  cells = rownames(md.RCC81)
)
RCC84.counts <- FeatureMatrix(
  fragments = frags.RCC84,
  features = combined.peaks,
  cells = rownames(md.RCC84)
)
RCC86.counts <- FeatureMatrix(
  fragments = frags.RCC86,
  features = combined.peaks,
  cells = rownames(md.RCC86)
)
RCC94.counts <- FeatureMatrix(
  fragments = frags.RCC94,
  features = combined.peaks,
  cells = rownames(md.RCC94)
)
RCC96.counts <- FeatureMatrix(
  fragments = frags.RCC96,
  features = combined.peaks,
  cells = rownames(md.RCC96)
)
RCC99.counts <- FeatureMatrix(
  fragments = frags.RCC99,
  features = combined.peaks,
  cells = rownames(md.RCC99)
)
RCC101.counts <- FeatureMatrix(
  fragments = frags.RCC101,
  features = combined.peaks,
  cells = rownames(md.RCC101)
)
RCC103.counts <- FeatureMatrix(
  fragments = frags.RCC103,
  features = combined.peaks,
  cells = rownames(md.RCC103)
)
RCC104.counts <- FeatureMatrix(
  fragments = frags.RCC104,
  features = combined.peaks,
  cells = rownames(md.RCC104)
)
RCC106.counts <- FeatureMatrix(
  fragments = frags.RCC106,
  features = combined.peaks,
  cells = rownames(md.RCC106)
)
RCC112.counts <- FeatureMatrix(
  fragments = frags.RCC112,
  features = combined.peaks,
  cells = rownames(md.RCC112)
)
RCC113.counts <- FeatureMatrix(
  fragments = frags.RCC113,
  features = combined.peaks,
  cells = rownames(md.RCC113)
)
RCC114.counts <- FeatureMatrix(
  fragments = frags.RCC114,
  features = combined.peaks,
  cells = rownames(md.RCC114)
)
RCC115.counts <- FeatureMatrix(
  fragments = frags.RCC115,
  features = combined.peaks,
  cells = rownames(md.RCC115)
)
RCC119.counts <- FeatureMatrix(
  fragments = frags.RCC119,
  features = combined.peaks,
  cells = rownames(md.RCC119)
)
RCC120.counts <- FeatureMatrix(
  fragments = frags.RCC120,
  features = combined.peaks,
  cells = rownames(md.RCC120)
)

metadata_RCC81 <- read.csv(
  file = ".../scATAC-seq/RCC81/singlecell.csv",
  header = TRUE,
  row.names = 1
)
metadata_RCC84 <- read.csv(
  file = ".../scATAC-seq/RCC84/singlecell.csv",
  header = TRUE,
  row.names = 1
)
metadata_RCC86 <- read.csv(
  file = ".../scATAC-seq/RCC86/singlecell.csv",
  header = TRUE,
  row.names = 1
)
metadata_RCC94 <- read.csv(
  file = ".../scATAC-seq/RCC94/singlecell.csv",
  header = TRUE,
  row.names = 1
)
metadata_RCC96 <- read.csv(
  file = ".../scATAC-seq/RCC96/singlecell.csv",
  header = TRUE,
  row.names = 1
)
metadata_RCC99 <- read.csv(
  file = ".../scATAC-seq/RCC99/singlecell.csv",
  header = TRUE,
  row.names = 1
)
metadata_RCC101 <- read.csv(
  file = ".../scATAC-seq/RCC101/singlecell.csv",
  header = TRUE,
  row.names = 1
)
metadata_RCC103 <- read.csv(
  file = ".../scATAC-seq/RCC103/singlecell.csv",
  header = TRUE,
  row.names = 1
)
metadata_RCC104 <- read.csv(
  file = ".../scATAC-seq/RCC104/singlecell.csv",
  header = TRUE,
  row.names = 1
)
metadata_RCC106 <- read.csv(
  file = ".../scATAC-seq/RCC106/singlecell.csv",
  header = TRUE,
  row.names = 1
)
metadata_RCC112 <- read.csv(
  file = ".../scATAC-seq/RCC112/singlecell.csv",
  header = TRUE,
  row.names = 1
)
metadata_RCC113 <- read.csv(
  file = ".../scATAC-seq/RCC113/singlecell.csv",
  header = TRUE,
  row.names = 1
)
metadata_RCC114 <- read.csv(
  file = ".../scATAC-seq/RCC114/singlecell.csv",
  header = TRUE,
  row.names = 1
)
metadata_RCC115 <- read.csv(
  file = ".../scATAC-seq/RCC115/singlecell.csv",
  header = TRUE,
  row.names = 1
)
metadata_RCC119 <- read.csv(
  file = ".../scATAC-seq/RCC119/singlecell.csv",
  header = TRUE,
  row.names = 1
)
metadata_RCC120 <- read.csv(
  file = ".../scATAC-seq/RCC120/singlecell.csv",
  header = TRUE,
  row.names = 1
)

RCC81_assay <- CreateChromatinAssay(RCC81.counts, fragments = frags.RCC81)
RCC81 <- CreateSeuratObject(RCC81_assay, assay = "ATAC", meta.data = metadata_RCC81)
RCC84_assay <- CreateChromatinAssay(RCC84.counts, fragments = frags.RCC84)
RCC84 <- CreateSeuratObject(RCC84_assay, assay = "ATAC", meta.data = metadata_RCC84)
RCC86_assay <- CreateChromatinAssay(RCC86.counts, fragments = frags.RCC86)
RCC86 <- CreateSeuratObject(RCC86_assay, assay = "ATAC", meta.data = metadata_RCC86)
RCC94_assay <- CreateChromatinAssay(RCC94.counts, fragments = frags.RCC94)
RCC94 <- CreateSeuratObject(RCC94_assay, assay = "ATAC", meta.data = metadata_RCC94)
RCC96_assay <- CreateChromatinAssay(RCC96.counts, fragments = frags.RCC96)
RCC96 <- CreateSeuratObject(RCC96_assay, assay = "ATAC", meta.data = metadata_RCC96)
RCC99_assay <- CreateChromatinAssay(RCC99.counts, fragments = frags.RCC99)
RCC99 <- CreateSeuratObject(RCC99_assay, assay = "ATAC", meta.data = metadata_RCC99)
RCC101_assay <- CreateChromatinAssay(RCC101.counts, fragments = frags.RCC101)
RCC101 <- CreateSeuratObject(RCC101_assay, assay = "ATAC", meta.data = metadata_RCC101)
RCC103_assay <- CreateChromatinAssay(RCC103.counts, fragments = frags.RCC103)
RCC103 <- CreateSeuratObject(RCC103_assay, assay = "ATAC", meta.data = metadata_RCC103)
RCC104_assay <- CreateChromatinAssay(RCC104.counts, fragments = frags.RCC104)
RCC104 <- CreateSeuratObject(RCC104_assay, assay = "ATAC", meta.data = metadata_RCC104)
RCC106_assay <- CreateChromatinAssay(RCC106.counts, fragments = frags.RCC106)
RCC106 <- CreateSeuratObject(RCC106_assay, assay = "ATAC", meta.data = metadata_RCC106)
RCC112_assay <- CreateChromatinAssay(RCC112.counts, fragments = frags.RCC112)
RCC112 <- CreateSeuratObject(RCC112_assay, assay = "ATAC", meta.data = metadata_RCC112)
RCC113_assay <- CreateChromatinAssay(RCC113.counts, fragments = frags.RCC113)
RCC113 <- CreateSeuratObject(RCC113_assay, assay = "ATAC", meta.data = metadata_RCC113)
RCC114_assay <- CreateChromatinAssay(RCC114.counts, fragments = frags.RCC114)
RCC114 <- CreateSeuratObject(RCC114_assay, assay = "ATAC", meta.data = metadata_RCC114)
RCC115_assay <- CreateChromatinAssay(RCC115.counts, fragments = frags.RCC115)
RCC115 <- CreateSeuratObject(RCC115_assay, assay = "ATAC", meta.data = metadata_RCC115)
RCC119_assay <- CreateChromatinAssay(RCC119.counts, fragments = frags.RCC119)
RCC119 <- CreateSeuratObject(RCC119_assay, assay = "ATAC", meta.data = metadata_RCC119)
RCC120_assay <- CreateChromatinAssay(RCC120.counts, fragments = frags.RCC120)
RCC120 <- CreateSeuratObject(RCC120_assay, assay = "ATAC", meta.data = metadata_RCC120)

RCC81$dataset <- 'RCC81'
RCC84$dataset <- 'RCC84'
RCC86$dataset <- 'RCC86'
RCC94$dataset <- 'RCC94'
RCC96$dataset <- 'RCC96'
RCC99$dataset <- 'RCC99'
RCC101$dataset <- 'RCC101'
RCC103$dataset <- 'RCC103'
RCC104$dataset <- 'RCC104'
RCC106$dataset <- 'RCC106'
RCC112$dataset <- 'RCC112'
RCC113$dataset <- 'RCC113'
RCC114$dataset <- 'RCC114'
RCC115$dataset <- 'RCC115'
RCC119$dataset <- 'RCC119'
RCC120$dataset <- 'RCC120'

##################################################################
##################################################################
###################. Count Fragments #############################
##################################################################
##################################################################

total.frags.RCC81 <- CountFragments(".../scATAC-seq/RCC81/fragments.tsv.gz")
rownames(total.frags.RCC81) <- total.frags.RCC81$CB
RCC81$fragments <- total.frags.RCC81[colnames(RCC81), "frequency_count"]

total.frags.RCC84 <- CountFragments(".../scATAC-seq/RCC84/fragments.tsv.gz")
rownames(total.frags.RCC84) <- total.frags.RCC84$CB
RCC84$fragments <- total.frags.RCC84[colnames(RCC84), "frequency_count"]

total.frags.RCC94 <- CountFragments(".../scATAC-seq/RCC94/fragments.tsv.gz")
rownames(total.frags.RCC94) <- total.frags.RCC94$CB
RCC94$fragments <- total.frags.RCC94[colnames(RCC94), "frequency_count"]

total.frags.RCC96 <- CountFragments(".../scATAC-seq/RCC96/fragments.tsv.gz")
rownames(total.frags.RCC96) <- total.frags.RCC96$CB
RCC96$fragments <- total.frags.RCC96[colnames(RCC96), "frequency_count"]

total.frags.RCC99 <- CountFragments(".../scATAC-seq/RCC99/fragments.tsv.gz")
rownames(total.frags.RCC99) <- total.frags.RCC99$CB
RCC99$fragments <- total.frags.RCC99[colnames(RCC99), "frequency_count"]

total.frags.RCC101 <- CountFragments(".../scATAC-seq/RCC101/fragments.tsv.gz")
rownames(total.frags.RCC101) <- total.frags.RCC101$CB
RCC101$fragments <- total.frags.RCC101[colnames(RCC101), "frequency_count"]

total.frags.RCC103 <- CountFragments(".../scATAC-seq/RCC103/fragments.tsv.gz")
rownames(total.frags.RCC103) <- total.frags.RCC103$CB
RCC103$fragments <- total.frags.RCC103[colnames(RCC103), "frequency_count"]

total.frags.RCC104 <- CountFragments(".../scATAC-seq/RCC104/fragments.tsv.gz")
rownames(total.frags.RCC104) <- total.frags.RCC104$CB
RCC104$fragments <- total.frags.RCC104[colnames(RCC104), "frequency_count"]

total.frags.RCC106 <- CountFragments(".../scATAC-seq/RCC106/fragments.tsv.gz")
rownames(total.frags.RCC106) <- total.frags.RCC106$CB
RCC106$fragments <- total.frags.RCC106[colnames(RCC106), "frequency_count"]

total.frags.RCC112 <- CountFragments(".../scATAC-seq/RCC112/fragments.tsv.gz")
rownames(total.frags.RCC112) <- total.frags.RCC112$CB
RCC112$fragments <- total.frags.RCC112[colnames(RCC112), "frequency_count"]

total.frags.RCC113 <- CountFragments(".../scATAC-seq/RCC113/fragments.tsv.gz")
rownames(total.frags.RCC113) <- total.frags.RCC113$CB
RCC113$fragments <- total.frags.RCC113[colnames(RCC113), "frequency_count"]

total.frags.RCC114 <- CountFragments(".../scATAC-seq/RCC114/fragments.tsv.gz")
rownames(total.frags.RCC114) <- total.frags.RCC114$CB
RCC114$fragments <- total.frags.RCC114[colnames(RCC114), "frequency_count"]

total.frags.RCC115 <- CountFragments(".../scATAC-seq/RCC115/fragments.tsv.gz")
rownames(total.frags.RCC115) <- total.frags.RCC115$CB
RCC115$fragments <- total.frags.RCC115[colnames(RCC115), "frequency_count"]

total.frags.RCC119 <- CountFragments(".../scATAC-seq/RCC119/fragments.tsv.gz")
rownames(total.frags.RCC119) <- total.frags.RCC119$CB
RCC119$fragments <- total.frags.RCC119[colnames(RCC119), "frequency_count"]

total.frags.RCC120 <- CountFragments(".../scATAC-seq/RCC120/fragments.tsv.gz")
rownames(total.frags.RCC120) <- total.frags.RCC120$CB
RCC120$fragments <- total.frags.RCC120[colnames(RCC120), "frequency_count"]


##################################################################
##################################################################
###################. Combine Samples.#############################
##################################################################
##################################################################

combined <- merge(x = RCC81,y = list(RCC84, RCC86, RCC94, RCC96, RCC99, RCC101, RCC103, RCC104, RCC106, RCC112, RCC113, RCC114, RCC115, RCC119, RCC120),add.cell.ids = c("RCC_81", "RCC_84", "RCC_86", "RCC_94", "RCC_96", "RCC_99", "RCC_101", "RCC_103", "RCC_104", "RCC_106", "RCC_112", "RCC_113", "RCC_114", "RCC_115", "RCC_119", "RCC_120"))

combined[["ATAC"]]

granges(combined)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"
Annotation(combined) <- annotations
combined <- subset(combined, subset = nFeature_ATAC > 200)

###################################################################
###################################################################
#####. Removal of Doublets using scDblFinder. #####################
###################################################################
###################################################################

peak_assay = GetAssayData(combined, layer="counts")
sce = scDblFinder(peak_assay, aggregateFeatures=TRUE, nfeatures=25, processing="normFeatures")
combined$doublet_class = sce@colData$scDblFinder.class
combined$doublet_class = factor(combined$doublet_class, levels = c('singlet', 'doublet'))
combined$doublet_score = sce@colData$scDblFinder.score

g <- ggplot(combined[[]], aes_string(x = 'doublet_score', color = 'doublet_class', fill = 'doublet_class')) + geom_histogram()
g + scale_color_manual(name = 'Final classification', values=c("#0571b0", "#ca0020"))+
  scale_fill_manual(name = 'Final classification',values=c("#0571b0", "#ca0020")) + theme_classic() + ylab('Number of cells') + xlab('Probability of cell being a\ndoublet (score)')+ 
  theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), text = element_text(size = 16, family = 'Helvetica'), legend.title=element_text(size=14 ,family = 'Helvetica'), legend.position = 'right')

combined = subset(
  x = combined,
  subset = doublet_class == 'singlet'
)

rm(sce)

############################################################
############################################################
###################. Remaining QC . ########################
############################################################
############################################################

combined <- FRiP(
  object = combined,
  assay = 'ATAC',
  total.fragments = 'fragments', verbose=TRUE)

combined$blacklist_fraction <- FractionCountsInRegion(
  object = combined,
  assay = 'ATAC',
  regions = blacklist_hg38)

combined <- NucleosomeSignal(object = combined)
combined <- TSSEnrichment(object = combined, fast = FALSE)

#QC Visualization

combined$high.tss <- ifelse(combined$TSS.enrichment > 2, 'TSS enrichment > 2', 'TSS enrichment < 2')
TSSPlot(combined, group.by = 'high.tss') + NoLegend()  & scale_color_manual(values = c('black', 'black')) & 
  theme(plot.title = element_blank(), text = element_text(size = 18)) 

combined$nucleosome_group <- ifelse(combined$nucleosome_signal > 4, 'NS > 4\n(low nucleosome-free fragments)', 'NS < 4')
FragmentHistogram(object = combined, group.by = 'nucleosome_group') + theme(text = element_text(size = 18))

options(repr.plot.width=7, repr.plot.height=6)
DensityScatter(combined, x = 'peak_region_fragments', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

options(repr.plot.width=20, repr.plot.height=6)
titles <- c('Fraction of reads\nin peaks', 'Peak region fragments', 'Transcription start site\nenrichment', 'Blacklist fraction', 'Nucleosome signal')
colnames <- c('FRiP','peak_region_fragments','TSS.enrichment', 'blacklist_fraction', 'nucleosome_signal')
plot_lst <- vector('list', length = 0)
for (i in 1:length(colnames)) {
  g <- ggplot(combined[[]], aes_string(x = factor(0), y = combined[[]][,colnames[i]])) + theme_classic() + geom_violin(color = 'black', fill = 'grey') +
    theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_text(size = 15), text = 
            element_text(size = 16, family = 'Helvetica'), axis.title.y=element_blank(), legend.position = "none") + geom_boxplot(width = 0.1) +
    ggtitle(titles[i]) + theme(plot.title = element_text(hjust = 0.5, family = 'Helvetica'))
  plot_lst[[i]] <- g 
}
cowplot::plot_grid(plotlist = plot_lst, nrow= 1, align = 'h')

VlnPlot(
  object = combined,
  features = c('nucleosome_signal'),
  pt.size = F, ncol = 3)

options(repr.plot.width=20, repr.plot.height=4)
titles <- c('Fraction of reads\nin peaks', 'Peak region fragments', 'Transcription start site\nenrichment', 'Blacklist fraction', 'Nucleosome signal')
colnames <- c('FRiP','peak_region_fragments','TSS.enrichment', 'blacklist_fraction', 'nucleosome_signal')

plot_lst <- vector('list', length = 0)
for (i in 1:length(colnames)) {
  g <- ggplot(combined[[]], aes_string(x = colnames[i])) + geom_density(fill = 'grey', color = '#616161')
  g <- g + geom_vline(aes_string(xintercept = median(combined[[colnames[i]]][,1]) + (3 * mad(combined[[colnames[i]]][,1]))))
  g <- g + geom_vline(aes_string(xintercept = median(combined[[colnames[i]]][,1]) - (3 * mad(combined[[colnames[i]]][,1]))))
  g <- g + ggtitle(paste('mincutoff = ', round(median(combined[[colnames[i]]][,1]) - (3 * mad(combined[[colnames[i]]][,1])), 2), '\nmaxcutoff = ', round(median(combined[[colnames[i]]][,1]) + (3 * mad(combined[[colnames[i]]][,1])), 2)))
  g <- g + theme_classic() + theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_text(size = 15), text = 
                                     element_text(size = 14, family = 'Helvetica'), legend.position = "none", plot.title = element_text(family = 'Helvetica', size = 15)) + xlab(titles[i])
  plot_lst[[i]] <- g
}


cowplot::plot_grid(plotlist = plot_lst, nrow= 1, align = 'h')

############################################################
############################################################
#########. After QC, Before Filtering . ####################
############################################################
############################################################

#Filtering
combined <- subset(
  x = combined,
  subset = peak_region_fragments < 10000 &
    peak_region_fragments > 3000 &
    FRiP > 0.20 &
    blacklist_fraction < 0.1 &
    nucleosome_signal < 2.5 &
    TSS.enrichment > 1.4
)

combined

#Processing
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 'q0')
combined <- RunSVD(combined)
options(repr.plot.width=10, repr.plot.height=6)
DepthCor(combined)

DefaultAssay(combined) <- 'ATAC'

combined <- RunUMAP(object = combined, reduction = 'lsi', dims = 2:30)

#Check for dataset effect
DimPlot(combined, group.by="dataset") 

combined <- RunHarmony(object = combined, group.by.vars = 'dataset', reduction.use = 'lsi', assay.use = 'ATAC', project.dim = FALSE)
combined <- RunUMAP(combined, dims = 2:30, reduction = 'harmony')
DimPlot(combined, group.by = 'dataset', pt.size = 0.1)
combined <- FindNeighbors(object = combined, reduction = 'harmony', dims = 2:30)
combined <- FindClusters(object = combined, verbose = FALSE, algorithm = 3, resolution = 0.05)
DimPlot(combined, group.by = 'seurat_clusters', pt.size = 0.1)

Idents(combined) <- "seurat_clusters"

options(repr.plot.width=6, repr.plot.height=6)
DimPlot(object = combined, label = TRUE, pt.size = 0.3, label.size = 5, label.color = 'black') + NoLegend() +
  theme(text = element_text(size = 14, family = 'Helvetica'),axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.text.x=element_blank()) +
  ylab('UMAP 2') + xlab('UMAP 1')

##################################################
##################################################
######### Cell type annotation ###################
##################################################
##################################################

#Broad Cell Type Identities
Idents(combined_filt) <- "seurat_clusters"
DimPlot(object = combined_filt, label = TRUE, pt.size = 0.3, label.size = 5, label.color = 'black') + NoLegend() +
  theme(text = element_text(size = 14, family = 'Helvetica'),axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.text.x=element_blank()) +
  ylab('UMAP 2') + xlab('UMAP 1')

#1, 2 are immune
#0 seems to be tumor
#3, 5 endo
#4 fibroblast/interstitial

############IMMUNE
#PTPRC
CoveragePlot(
  object = combined,
  group.by = 'seurat_clusters',
  region = "chr1-198636713-198759476", extend.upstream = 3000, extend.downstream = 3000,scale.factor = 1e8, heights=15
) & theme(axis.title.y = element_text(size = 20))

#GZMB
CoveragePlot(
  object = combined,
  group.by = 'seurat_clusters', extend.upstream = 3000, extend.downstream = 3000,scale.factor = 1e8, heights=15,
  region = "chr14-24628954-24636190"
) & theme(axis.title.y = element_text(size = 20))

#THEMIS
CoveragePlot(
  object = combined,
  group.by = 'seurat_clusters', extend.upstream = 3000, extend.downstream = 3000,scale.factor = 1e8, heights=15,
  region = "chr6-127694628-127920595
"
) & theme(axis.title.y = element_text(size = 20))

#CD34
CoveragePlot(
  object = combined,
  group.by = 'seurat_clusters', extend.upstream = 3000, extend.downstream = 3000,scale.factor = 1e8, heights=15,
  region = "chr4-15776328-15855232
"
) & theme(axis.title.y = element_text(size = 20))
############TUMOR

#CA9
CoveragePlot(
  object = combined,
  group.by = 'seurat_clusters',
  region = "chr9-35671928-35683159", extend.upstream = 3000, extend.downstream = 3000,scale.factor = 1e8, heights=15,
) & theme(axis.title.y = element_text(size = 20))

#KRT18
CoveragePlot(
  object = combined,
  group.by = 'seurat_clusters', extend.upstream = 3000, extend.downstream = 3000,scale.factor = 1e8, heights=15,
  region = "chr12-52946855-52954906"
) & theme(axis.title.y = element_text(size = 20))

#NDUFA4L2
CoveragePlot(
  object = combined,
  group.by = 'seurat_clusters', extend.upstream = 3000, extend.downstream = 3000,scale.factor = 1e8, heights=15,
  region = "chr12-57232903-57242762
") & theme(axis.title.y = element_text(size = 20))

############OTHER

#VWF
CoveragePlot(
  object = combined,
  group.by = 'seurat_clusters',
  region = "chr12-5946877-6126670", extend.upstream = 3000, extend.downstream = 3000,scale.factor = 1e8, heights=15,
) & theme(axis.title.y = element_text(size = 20))

#PECAM1
CoveragePlot(
  object = combined,
  group.by = 'seurat_clusters', extend.upstream = 3000, extend.downstream = 3000,scale.factor = 1e8, heights=15,
  region = "chr17-64317415-64392860"
) & theme(axis.title.y = element_text(size = 20))

#MCAM
CoveragePlot(
  object = combined,
  group.by = 'seurat_clusters', extend.upstream = 3000, extend.downstream = 3000,scale.factor = 1e8, heights=15,
  region = "chr11-119306529-119319130
") & theme(axis.title.y = element_text(size = 20))

#COL4A1
CoveragePlot(
  object = combined,
  group.by = 'seurat_clusters', extend.upstream = 3000, extend.downstream = 3000,scale.factor = 1e8, heights=15,
  region = "chr13-110146963-110309157
") & theme(axis.title.y = element_text(size = 20))

#SERPINE1
CoveragePlot(
  object = combined,
  group.by = 'seurat_clusters', extend.upstream = 3000, extend.downstream = 3000,scale.factor = 1e8, heights=15,
  region = "chr7-101125104-101141247
") & theme(axis.title.y = element_text(size = 20))

#COL1A1
CoveragePlot(
  object = combined,
  group.by = 'seurat_clusters', extend.upstream = 3000, extend.downstream = 3000,scale.factor = 1e8, heights=15,
  region = "chr17-50182101-50203631
") & theme(axis.title.y = element_text(size = 20))


#Rename seurat clusters
combined$cluster <- combined@active.ident

new.cluster.ids <- c("Cancer", "Immune", "Immune", "Endothelial", "Interstitial", "Endothelial")
names(new.cluster.ids) <- levels(combined)
combined <- RenameIdents(combined, new.cluster.ids)

combined$celltype <- Idents(combined)

Idents(combined) <- "seurat_clusters"
Idents(combined) <- "celltype"
Idents(combined) <- "dataset"

DimPlot(object = combined, label = TRUE, pt.size = 0.3, label.size = 5, label.color = 'black') + NoLegend() +
  theme(text = element_text(size = 14, family = 'Helvetica'),axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.text.x=element_blank()) +
  ylab('UMAP 2') + xlab('UMAP 1')

############IMMUNE

#PTPRC
CoveragePlot(
  object = combined,
  group.by = 'celltype',
  region = "chr1-198636713-198759476", extend.upstream = 3000, extend.downstream = 3000,scale.factor = 1e8
) & theme(axis.title.y = element_text(size = 20))

#GZMB
CoveragePlot(
  object = combined,
  group.by = 'celltype', extend.upstream = 3000, extend.downstream = 3000,
  region = "chr14-24628954-24636190"
) & theme(axis.title.y = element_text(size = 20))

#THEMIS
CoveragePlot(
  object = combined,
  group.by = 'celltype', extend.upstream = 3000, extend.downstream = 3000,
  region = "chr6-127694628-127920595
"
) & theme(axis.title.y = element_text(size = 20))

############TUMOR

#CA9
CoveragePlot(
  object = combined,
  group.by = 'celltype',
  region = "chr9-35671928-35683159", extend.upstream = 3000, extend.downstream = 3000,scale.factor = 1e8
) & theme(axis.title.y = element_text(size = 20))

#KRT18
CoveragePlot(
  object = combined,
  group.by = 'celltype', extend.upstream = 3000, extend.downstream = 3000,
  region = "chr12-52946855-52954906",scale.factor = 1e8
) & theme(axis.title.y = element_text(size = 20))

#MME
CoveragePlot(
  object = combined,
  group.by = 'celltype', extend.upstream = 3000, extend.downstream = 3000,
  region = "chr3-155022202-155185704
")

############OTHER

#VWF
CoveragePlot(
  object = combined,
  group.by = 'celltype',
  region = "chr12-5946877-6126670", extend.upstream = 3000, extend.downstream = 3000
)

#PECAM1
CoveragePlot(
  object = combined,
  group.by = 'celltype', extend.upstream = 3000, extend.downstream = 3000,scale.factor = 1e8,
  region = "chr17-64317415-64392860"
) & theme(axis.title.y = element_text(size = 20))

#MCAM
CoveragePlot(
  object = combined,
  group.by = 'celltype', extend.upstream = 3000, extend.downstream = 3000,
  region = "chr4-15966228-16085729
") & theme(axis.title.y = element_text(size = 20))

############HOUSEKEEPING
#GAPDH
CoveragePlot(
  object = combined,
  group.by = 'celltype', extend.upstream = 3000, extend.downstream = 3000,scale.factor = 1e8,
  region = "chr12-6532517-6540371
") & theme(axis.title.y = element_text(size = 20))

#TUBA1A
CoveragePlot(
  object = combined,
  group.by = 'celltype', extend.upstream = 3000, extend.downstream = 3000,scale.factor = 1e8,
  region = "chr12-49182795-49191080
") & theme(axis.title.y = element_text(size = 20))

#VCL
CoveragePlot(
  object = combined,
  group.by = 'celltype', extend.upstream = 3000, extend.downstream = 3000,scale.factor = 1e8,
  region = "chr10-73996116-74123363
") & theme(axis.title.y = element_text(size = 20))

#ACTB
CoveragePlot(
  object = combined,
  group.by = 'celltype', extend.upstream = 3000, extend.downstream = 3000,scale.factor = 1e8,
  region = "chr7-5525148-5532601

") & theme(axis.title.y = element_text(size = 20))











