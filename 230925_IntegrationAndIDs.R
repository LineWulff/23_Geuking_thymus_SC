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
# sample information
samp_trans <- function(x) {
  ifelse(x == "s6", "gf", ifelse(x == "s7", "sdmdmm2", "spf"))}

##### Read in samples and integrate based on rpca #####
s6 <- readRDS(paste(data_dir,"s6","s6.rds", sep="/"))
s7 <- readRDS(paste(data_dir,"s7","s7.rds", sep="/"))
s8 <- readRDS(paste(data_dir,"s8","s8.rds", sep="/"))

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
samp_list <- c(s6,s7,s8)
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

thymus <- RunUMAP(thymus, reduction = "pca", dims = 1:15)
thymus <- FindNeighbors(thymus, reduction = "pca", dims = 1:15)
res <- seq(0,5,0.1)
thymus <- FindClusters(thymus, resolution = res)

thymus@meta.data$colonization <- samp_trans(thymus@meta.data$orig.ident)

## Reduction with some QC measures
P1 <- DimPlot(thymus, group.by = "colonization")
P2 <- DimPlot(thymus, group.by = "Phase")
P1+P2
# Quite a bit of effect of cell cycle - worse without the regression
P1 <- DimPlot(thymus, group.by = "integrated_snn_res.0.2", label = T, label.size = 6)+NoLegend()
P1+P2
DimPlot(thymus, group.by = "integrated_snn_res.1")
P1 <- FeaturePlot(thymus, features = c("percent.mt"))+scale_color_viridis_c(); P2 <- FeaturePlot(thymus, features = c("nFeature_RNA"))+scale_color_viridis_c()
P1+P2
# set identity to resolution
Idents(thymus) <- 'integrated_snn_res.0.2'

# Aire expression limited to few cells in cl. 8
FeaturePlot(thymus, features = "Aire", order=T)+scale_color_gradientn(colours=mycols_b)

# Myeloid cells or B cells?
FeaturePlot(thymus, features = "Itgax", order=T)+scale_color_gradientn(colours=mycols_b)
FeaturePlot(thymus, features = "Cd79a", order=T)+scale_color_gradientn(colours=mycols_b)
VlnPlot(thymus, features = c("Cd79a","Itgax","Aire"), pt.size = 0)

#### Identifying DEGs ####
DefaultAssay(thymus) <- 'RNA'
Tot_DEGs <- FindAllMarkers(thymus, only.pos = T, logfc.threshold = 0.5)
write.csv(Tot_DEGs, file = paste(dato,"R1_thymus_DEGs_res.0.2.csv"))

data("cc.genes") ## Tirosh et al. 2015 HUMAN GENES, identify same in mice
sgenes <- rownames(thymus@assays$RNA@counts)[toupper(rownames(thymus@assays$RNA@counts)) %in% cc.genes$s.genes]
g2mgenes <- rownames(thymus@assays$RNA@counts)[toupper(rownames(thymus@assays$RNA@counts)) %in% cc.genes$g2m.genes]

# cluster 3 and 8 both heavily cell cycle 
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

# VlnPlot versions per cluster
VlnPlot(thymus, features = paste(unique(ST1$cluster),"1",sep = ""),
        pt.size = 0, ncol = 4)+NoLegend()

#### Kernfeld - T conv - Fig 2 / Sup tab 2 ####
ST2 <- read.csv("/Volumes/Promise RAID/Line/projects/23_Geuking_thymus_SC/refdata/SupTab2_Kernfeld.csv")
ST2 <- ST2[ST2$avg_diff>0,] #from 3645 to 1851
Kernfeld <- list()
for (clus in unique(ST2$cluster)){
  Kernfeld[[clus]] <- rownames(thymus@assays$RNA@data)[rownames(thymus@assays$RNA@data) %in% ST2[ST2$cluster==clus,]$gene]
}

DefaultAssay(thymus) <- 'RNA'
for (clus in unique(ST2$cluster)){
  if (length(unlist(Kernfeld[[clus]]))>500){nctrl=500} else {nctrl=length(unlist(Kernfeld[[clus]]))}
  thymus <- AddModuleScore(thymus, features = list(Kernfeld[[clus]]), ctrl = nctrl, name = clus, seed = 42)
}
for (clus in unique(ST2$cluster)){
  plot1 <- FeaturePlot(thymus, features = paste(clus,"1", sep=""))+
    scale_colour_gradientn(colours=mycols_b)
  print(plot1)
}

# VlnPlot
VlnPlot(subset(thymus, idents = c("0","1","2","3","4","6","8")), 
        features = paste(unique(ST2$cluster),"1",sep = ""),
        pt.size = 0, ncol = 3)+NoLegend()
VlnPlot(subset(thymus, idents =  c("0","1","2","3","4","6","8")),
        features = c("Cd3e","Cd3d","Cd4","Cd8a"),
        pt.size = 0, ncol = 2)

#### Kernfeld - NCL - Fig 3 / Sup tab 4 ####
ST4 <- read.csv("/Volumes/Promise RAID/Line/projects/23_Geuking_thymus_SC/refdata/SupTab4_Kernfeld.csv")
ST4 <- ST4[!is.na(ST4$avg_diff),]
ST4 <- ST4[ST4$avg_diff>0,] #from 2473 to 1620
Kernfeld <- list()
for (clus in unique(ST4$cluster)){
  Kernfeld[[clus]] <- rownames(thymus@assays$RNA@data)[rownames(thymus@assays$RNA@data) %in% ST4[ST4$cluster==clus,]$gene]
}

DefaultAssay(thymus) <- 'RNA'
for (clus in unique(ST4$cluster)){
  if (length(unlist(Kernfeld[[clus]]))>500){nctrl=500} else {nctrl=length(unlist(Kernfeld[[clus]]))}
  thymus <- AddModuleScore(thymus, features = list(Kernfeld[[clus]]), ctrl = nctrl, name = clus, seed = 42)
}
for (clus in unique(ST4$cluster)){
  plot1 <- FeaturePlot(thymus, features = paste(clus,"1", sep=""))+
    scale_colour_gradientn(colours=mycols_b)
  print(plot1)
}

# VlnPlot
VlnPlot(subset(thymus, idents = c("5","7","6")), 
        features = paste(unique(ST4$cluster),"1",sep = ""),
        pt.size = 0, ncol = 3)+NoLegend()

#TCR g module
Tcrg <- rownames(thymus@assays$RNA@counts)[startsWith(rownames(thymus@assays$RNA@counts),"Tcrg")]
thymus <- AddModuleScore(thymus, features = list(Tcrg), ctrl = length(Tcrg), name = "Tcrg", seed = 42)
VlnPlot(subset(thymus, idents = c("5","7","6")), 
        features = c("Tcrg1","Cd4","Klrk1","Tbx21","Klrb1c"),
        pt.size = 0, ncol = 3)+NoLegend()

#### Kernfeld - TEC - Fig 5 / Sup tab 6 ####
ST6 <- read.csv("/Volumes/Promise RAID/Line/projects/23_Geuking_thymus_SC/refdata/SupTab6_Kernfeld.csv")
ST6 <- ST6[!is.na(ST6$avg_diff),]
ST6 <- ST6[ST6$avg_diff>0,]
Kernfeld <- list()
for (clus in unique(ST6$cluster)){
  Kernfeld[[clus]] <- rownames(thymus@assays$RNA@data)[rownames(thymus@assays$RNA@data) %in% ST6[ST6$cluster==clus,]$gene]
}

DefaultAssay(thymus) <- 'RNA'
for (clus in unique(ST6$cluster)){
  if (length(unlist(Kernfeld[[clus]]))>500){nctrl=500} else {nctrl=length(unlist(Kernfeld[[clus]]))}
  thymus <- AddModuleScore(thymus, features = list(Kernfeld[[clus]]), ctrl = nctrl, name = clus, seed = 42)
}
for (clus in unique(ST6$cluster)){
  plot1 <- FeaturePlot(thymus, features = paste(clus,"1", sep=""))+
    scale_colour_gradientn(colours=mycols_b)
  print(plot1)
}

# Plotting genes and modules from Kernfeld et al. to ID different TEC populations
DefaultAssay(thymus) <- 'RNA'
pheatmap(thymus@meta.data[thymus@meta.data$integrated_snn_res.0.2=="11",paste(unique(ST6$cluster),"1",sep = "")],
         show_rownames = F)
FeatureScatter(subset(thymus, idents = "11"),
               feature1 = "Tnfrsf11a", feature2 = "Fezf2")
DoHeatmap(subset(thymus, idents = "11"),
          slot="data",
          features = unlist(Kernfeld),#c("Tnfrsf11a","Fezf2","Aire",
                      # "Ackr4","Ly75","Cd83","Prss16","Tbata","Psmb11"),
          draw.lines = F)+
  scale_fill_gradientn(colours=mycols_b)

# subset cl 11 - plotting genes from text associated with mTEC vs cTEC
cl11 <- subset(thymus, idents = "11")
TECgenes <- c(#mTEC
              "Fezf2","Aire", 
              "Srgn","S100a8","Ccl5","Cd52","Ubd","S100a9","S100a14",
              "Csn2",'Cxcl9',"Cyba","Nts","Tkt","Pyy","Ifi30",
              "Cd74","H2-Ab1","H2-Eb1",
              # mTEC prec./ immature mTEC
              "H2-Aa","Ccl21a","Krt5","Tnfrsf11a",
              "Ccl21a","Gm8113","Krt14","Krt5","Tagln","Col6a1","Fcgbp","Lifr",
              "Flna","Rgs5",'Itga6','Mgp','Col4a1','Hbb-bs','Tpm2','Ascl1','C3','Acta2','Gas1',
              # cTEC
              "Ackr4","Ly75","Cd83","Prss16","Tbata","Psmb11","Trp63","Foxn1",
              'Tbata','Ctsl','Prss16','Ccl25','Cstb','Psmb11','Cxcl12','Ly75',
              'Tmem131l','Crip3','Isg20','Fxyd2','Ndrg3','Pax1','Prxl2b','Pithd1',
              'Psmb4','Lamp2','Sgpl1',
              # Tuft cells
              "Lrmp","Avil","Dclk1","Gng13","L1cam","Sox9","Spib","Alox5ap",
              "Gnat3","Gnb3","Plac8","Gng13","Lrmp", #"Pou2f3"
              "Rgs13","Cited2","Krt19","Ethe1","Fxyd6","Cystm1","Reep5","Ovol3","Ltc4s",
              "Avil","Chil1","Alox5ap","Gm26797","A630073D07Rik","Cap1","Anxa4",
              )
TECgenes <- unique(TECgenes[TECgenes %in% rownames(thymus@assays$RNA@data)])

pheatmap(cl11@assays$RNA@data[TECgenes,],
         show_colnames = T, show_rownames = F,
         cluster_rows = F)





#### Checking subsets within cl 6a - dendritic cells ####
FeatureScatter(subset(thymus, idents = "9"), feature1 = "H2-Aa",feature2 = "Itgam", group.by = "integrated_snn_res.0.9")
FeatureScatter(subset(thymus, idents = "9"), feature1 = "H2-Aa",feature2 = "Bst2", group.by = "integrated_snn_res.0.9")
# cl 20 at res.0.9 are pDCs, Itgam+ are cDC2 and H2-Aa+ are cDC1

#### Distribution per sample ####
# meta_col, cell_vec, Seu_obj,splitgroup
# run the "Distribution_functions.R" script from useful_scripts on LineWulffs github
dist_df <- perc_function_samp("integrated_snn_res.0.2",colnames(thymus),thymus,"orig.ident")
ggplot(data=dist_df, aes(x=samp, y=percent, fill=cluster)) +
  geom_bar(stat= "identity", colour="black")+
  theme_classic()
ggplot(data=dist_df, aes(x=samp , y=percent, fill=cluster)) +
  geom_bar(stat= "identity", colour="black")+
  facet_grid(.~cluster)+
  theme_classic()


#### save object ####
#saveRDS(thymus, file = "2309_R1_Geuking_ThymusSC.rds")
