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
proj_dir <- "/Users/linewulff/Documents/work/projects/23_Geuking_thymus_SC"
data_dir <- "/Volumes/Promise RAID/Line/projects/23_Geuking_thymus_SC"

##### Read in Seurat object #####
# from initial integration adn ID
thymus <- readRDS("2309_Geuking_ThymusSC_allcond.rds")

#### Comaprison between conditions ####
DefaultAssay(thymus) <- "RNA"
for (clus in unique(thymus@meta.data$res.0.2)){
  #reset idents for next round of loop
  Idents(thymus) <- 'res.0.2'
  # subset per cluster
  cell_sub <- subset(thymus, idents = clus)
  Idents(thymus) <- 'colonization'
  
  ### "sdmdmm2" vs "gf"
  # Run DEGA
  DEGs <- FindMarkers(cell_sub, ident.1 = "sdmdmm2", ident.2 = "gf", 
                      logfc.threshold = 0, # slower but necessary for pretty volcano plots
                      test.use = "DESeq2")
  # significance limits as IDs
  
  # plot volcano plot
  volc <- ggplot(DEGs, aes(x = , y = , colour = ))+
    geom_point()+
    geom_hline()+ #sign.
    geom_vline()+ #sign.
    theme_classic()
  
  # Subset significant and save csv
  
  ### "spf" vs "gf"
  # Run DEGA
  DEGs <- FindMarkers(cell_sub, ident.1 = "spf", ident.2 = "gf", 
                      logfc.threshold = 0, # slower but necessary for pretty volcano plots
                      test.use = "DESeq2")
  # significance limits as IDs
  
  # plot volcano plot
  volc <- ggplot(DEGs, aes(x = , y = , colour = ))+
    geom_point()+
    geom_hline()+ #sign.
    geom_vline()+ #sign.
    theme_classic()
  
  # Subset significant and save csv

  
  ### "sdmdmm2" vs "spf"
  # Run DEGA
  DEGs <- FindMarkers(cell_sub, ident.1 = "sdmdmm2", ident.2 = "spf", 
                      logfc.threshold = 0, # slower but necessary for pretty volcano plots
                      test.use = "DESeq2")
  # significance limits as IDs
  
  # plot volcano plot
  volc <- ggplot(DEGs, aes(x = , y = , colour = ))+
    geom_point()+
    geom_hline()+ #sign.
    geom_vline()+ #sign.
    theme_classic()
  
  # Subset significant and save csv
}




