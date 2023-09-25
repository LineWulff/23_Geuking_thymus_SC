#' R script for integrating samples and running clustering, ID'ing etc.
#' Author: Line Wulff
#' Date (created): 23-09-22

###### Libraries ######
library(gplots)
library(dplyr)
library(stringr)
library(scales)
library(ggplot2)
library(viridis)
library(rgl)
library(Seurat)
library(pheatmap)

###### variables used through script ######
rm(list=ls())
#date in format year_month_day
dato <- str_sub(str_replace_all(Sys.Date(),"-","_"), 3, -1)
# colour string for imputation and overlays
mycols_b <- c("blue","#bdbdbd","#d9d9d9","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F")
# project
project <- "s7"
proj_dir <- "/Users/linewulff/Documents/work/projects/23_Geuking_thymus_SC"
data_dir <- "/Volumes/Promise RAID/Line/projects/23_Geuking_thymus_SC"

##### Read in samples and integrate based on rpca #####
s6 <- readRDS(paste(data_dir,"s6","s6.rds", sep="/"))
s7 <- readRDS(paste(data_dir,"s7","s7.rds", sep="/"))

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
samp_list <- c(s6,s7)
features <- SelectIntegrationFeatures(object.list = samp_list)
samp_list <- lapply(X = samp_list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = TRUE, 
                 vars.to.regress = c("percent.mt","nFeature_RNA","nCount_RNA","S.Score","G2M.Score"))
  x <- RunPCA(x, features = features, verbose = TRUE)})

anchors <- FindIntegrationAnchors(object.list = samp_list, anchor.features = features, reduction = "rpca")
thymus <- IntegrateData(anchorset = anchors)

#### Initial analysis on merged data ####
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(thymus) <- "integrated"

# Run the standard workflow for visualization and clustering
thymus <- ScaleData(thymus, verbose = TRUE, 
                    vars.to.regress = c("percent.mt","nFeature_RNA","nCount_RNA","S.Score","G2M.Score"))
thymus <- RunPCA(thymus, npcs = 30, verbose = FALSE)
ElbowPlot(thymus)

thymus <- RunUMAP(thymus, reduction = "pca", dims = 1:17)
thymus <- FindNeighbors(thymus, reduction = "pca", dims = 1:17)
res <- seq(0,5,0.1)
thymus <- FindClusters(thymus, resolution = res)

## Reduction with some QC measures
DimPlot(thymus, group.by = "orig.ident")
DimPlot(thymus, group.by = "Phase")
# Quite a bit of effect of cell cycle - worse without the regression
DimPlot(thymus, group.by = "integrated_snn_res.0.1")
DimPlot(thymus, group.by = "integrated_snn_res.1")
P1 <- FeaturePlot(thymus, features = c("percent.mt"))+scale_color_viridis_c(); P2 <- FeaturePlot(thymus, features = c("nFeature_RNA"))+scale_color_viridis_c()
P1+P2

# Aire expression limited to few cells in cl. 8
FeaturePlot(thymus, features = "Aire", order=T)+scale_color_gradientn(colours=mycols_b)

# Myeloid cells or B cells?
FeaturePlot(thymus, features = "Itgax", order=T)+scale_color_gradientn(colours=mycols_b)
FeaturePlot(thymus, features = "Cd79a", order=T)+scale_color_gradientn(colours=mycols_b)

# Split cl. 6 into res.0.2 and rename cells (2 clusters, 6a and 6b)
# check cluster 6 only splits to two clusters
unique(thymus@meta.data[thymus@meta.data$integrated_snn_res.0.1==6,]$integrated_snn_res.0.2)
# make new res0.1
thymus@meta.data$res.0.1 <- as.character(thymus@meta.data$integrated_snn_res.0.1)
thymus@meta.data[thymus@meta.data$integrated_snn_res.0.2==9,]$res.0.1 <- "6a"
thymus@meta.data[thymus@meta.data$integrated_snn_res.0.2==10,]$res.0.1 <- "6b"
thymus@meta.data$res.0.1 <- factor(thymus@meta.data$res.0.1, levels = c("0","1","2","3","4","5","6a","6b","7","8"))
DimPlot(thymus, group.by = "res.0.1")

#### Identifying DEGs ####
Idents(thymus) <- 'res.0.1'
DefaultAssay(thymus) <- 'RNA'
Tot_DEGs <- FindAllMarkers(thymus, only.pos = T, logfc.threshold = 0.5)
write.csv(Tot_DEGs, file = paste(dato,"R1_thymus_DEGs_res.0.1.csv"))

data("cc.genes") ## Tirosh et al. 2015 HUMAN GENES, identify same in mice
sgenes <- rownames(thymus@assays$RNA@counts)[toupper(rownames(thymus@assays$RNA@counts)) %in% cc.genes$s.genes]
g2mgenes <- rownames(thymus@assays$RNA@counts)[toupper(rownames(thymus@assays$RNA@counts)) %in% cc.genes$g2m.genes]

# almost all in cluster 5 for both
Tot_DEGs[Tot_DEGs$gene %in% sgenes,] 
Tot_DEGs[Tot_DEGs$gene %in% g2mgenes,] 

#### Using data from Kernfeld et al.  to ID clusters ####
# https://doi.org/10.1016/j.immuni.2018.04.015
ST1 <- read.csv("/Volumes/Promise RAID/Line/projects/23_Geuking_thymus_SC/refdata/SupTab1_Kernfeld.csv")
Kernfeld <- list()
for (clus in unique(ST1$cluster)){
  Kernfeld[[clus]] <- rownames(thymus@assays$RNA@data)[rownames(thymus@assays$RNA@data) %in% ST1[ST1$cluster==clus,]$gene]
}

# Adding and plotting modulescores for each gene module from SupTab 1 Kernfeld et al.
DefaultAssay(thymus) <- 'RNA'
for (clus in unique(ST1$cluster)){
  thymus <- AddModuleScore(thymus, features = list(Kernfeld[[clus]]), ctrl = length(unlist(Kernfeld[[clus]])), name = clus, seed = 42)
}
for (clus in unique(ST1$cluster)){
  plot1 <- FeaturePlot(thymus, features = paste(clus,"1", sep=""))+
    scale_colour_gradientn(colours=mycols_b)
  print(plot1)
}

# Averaging thymic sc object and modulescores to plot neatly in heatmap
Av_thym <- AverageExpression(thymus, return.seurat = T, group.by = c("res.0.1","orig.ident"))
for (clus in unique(ST1$cluster)){
  Av_thym <- AddModuleScore(Av_thym, features = list(Kernfeld[[clus]]), ctrl = length(unlist(Kernfeld[[clus]])), name = clus, seed = 42)
}
# modulescore averaged at cluster level (and ID split) to check matches from Kernfeld Fig 1
pheatmap(scale(t(Av_thym@meta.data[,c(6:17)]), center = T),
         treeheight_row = 0)



#### Checking subsets within cl 6a - dendritic cells ####
FeatureScatter(subset(thymus, idents = "6a"), feature1 = "H2-Aa",feature2 = "Itgam", group.by = "integrated_snn_res.0.9")
FeatureScatter(subset(thymus, idents = "6a"), feature1 = "H2-Aa",feature2 = "Bst2", group.by = "integrated_snn_res.0.9")
# cl 20 at res.0.9 are pDCs, Itgam+ are cDC2 and H2-Aa+ are cDC1





