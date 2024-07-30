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

####################################################################################################
#Subset out Cancer and Immune#######################################################################
####################################################################################################
filtered <- subset(combined, idents=c("Cancer","Immune"))

ERVs = list(
  'ERV2256' = 'chr6-88662251-88670881',
  'ERV2637' = 'chr7-135174033-135184466',
  'ERV4818' = 'chr20-15980070-15986191',
  'ERV5875' = 'chr1-46325548-46332308',
  'ERV3797' = 'chr12-52109035-52121408',
  'CA9' = 'chr9-35673928-35681159',
  'PTPRC' = 'chr1-198638457-198757476',
  'ACTB' = 'chr7-5526409-5563902',
  'VCL' = 'chr10-73995193-74121363',
  'GAPDH' = 'chr12-6534512-6538374'
)

granges_list = list()

for (i in 1:length(ERVs)) {
  grange = Signac:::FindRegion(combined,ERVs[[i]])
  grange$ERV_loci = names(ERVs)[i]
  granges_list[[i]] = grange
}

erv_ranges = do.call(c, as(granges_list, "GRangesList"))
erv_ranges[1]

mycolors <- colorRampPalette(brewer.pal(4, "Set2"))(4)

#CoveragePlot at ERVs
#ERV2256
CoveragePlot(
  object = filtered,
  group.by = 'celltype',
  region = "chr6-88662251-88670881",ranges=erv_ranges, peaks=F, extend.upstream = 3000, extend.downstream = 3000, scale.factor = 1e8, annotation=F
)  & theme(axis.title.y = element_text(size = 20))

#ERV2637
CoveragePlot(
  object = filtered,
  group.by = 'celltype',
  region = "chr7-135174033-135184466",ranges=erv_ranges, peaks=F, extend.upstream = 3000, extend.downstream = 3000, scale.factor = 1e8, annotation=F
) & theme(axis.title.y = element_text(size = 20))

#ERV4818
CoveragePlot(
  object = filtered,
  group.by = 'celltype',
  region = "chr20-15980070-15986191", ranges=erv_ranges, peaks=F, extend.upstream = 3000, extend.downstream = 3000, scale.factor = 1e8, annotation=F
) & theme(axis.title.y = element_text(size = 20)) 

#ERV5875
CoveragePlot(
  object = filtered,
  group.by = 'celltype',
  region = "chr1-46325548-46332308",ranges=erv_ranges, peaks=F, extend.upstream = 3000, extend.downstream = 3000, scale.factor = 1e8, annotation=F
) & theme(axis.title.y = element_text(size = 20)) 

#ERV3797
CoveragePlot(
  object = filtered,
  group.by = 'celltype',
  region = "chr12-52109035-52121408",ranges=erv_ranges, peaks=F, extend.upstream = 3000, extend.downstream = 3000, scale.factor = 1e8, annotation=F
) & theme(axis.title.y = element_text(size = 20)) 

#CA9
CoveragePlot(
  object = filtered,
  group.by = 'celltype',
  region = "chr9-35673928-35681159", ranges=erv_ranges, peaks=F, extend.upstream = 3000, extend.downstream = 3000, scale.factor = 1e8, annotation=F
) & theme(axis.title.y = element_text(size = 20))

#PTPRC
CoveragePlot(
  object = filtered,
  group.by = 'celltype',
  region = "chr1-198638457-198757476", ranges=erv_ranges, peaks=F, extend.upstream = 3000, extend.downstream = 3000, scale.factor = 1e8, annotation=F
) & theme(axis.title.y = element_text(size = 20))

#ACTB
CoveragePlot(
  object = filtered,
  group.by = 'celltype',
  region = "chr7-5526409-5563902", ranges=erv_ranges, peaks=F, extend.upstream = 3000, extend.downstream = 3000, scale.factor = 1e8, annotation=F
) & theme(axis.title.y = element_text(size = 20))

#VCL
CoveragePlot(
  object = filtered,
  group.by = 'celltype',
  region = "chr10-73995193-74121363", ranges=erv_ranges, peaks=F, extend.upstream = 3000, extend.downstream = 3000, scale.factor = 1e8, annotation=F
) & theme(axis.title.y = element_text(size = 20))

#GAPDH
CoveragePlot(
  object = filtered,
  group.by = 'celltype',
  region = "chr12-6534512-6538374", ranges=erv_ranges, peaks=F, extend.upstream = 3000, extend.downstream = 3000, scale.factor = 1e8, annotation=F
) & theme(axis.title.y = element_text(size = 20))



