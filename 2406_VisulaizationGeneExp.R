#' R script for visulaizing gene expression
#' Author: Line Wulff
#' Date (created): 23-09-22

###### Libraries ######
# Read in libraries - install first if not installed
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
volc_cols <- c("gf"= "#F8766D","sdmdmm2"="#00BA38","spf"="#619CFF","NA"="lightgrey")
# project
proj_dir <- "/Users/linewulff/Documents/work/projects/23_Geuking_thymus_SC"
out_dir <- paste(proj_dir,"GeneExp",sep = '/')

##### Read in Seurat object #####
# from initial integration adn ID
thymus <- readRDS("2309_Geuking_ThymusSC_allcond.rds")

genenameinlist("Stat",thymus)

RNA_genes <- c("S1pr2","S1pr3","S1pr4","S1pr5","Stat3")
int_genes <- c("S1pr1")


setwd(out_dir)
#### UMAP overlays ####
for (gene in int_genes){
pdf(paste(dato,"_thymusSC_UMAP_GeneExp_",gene,".pdf",sep = ""), height = 4,width = 4)
plot1 <- FeaturePlot(thymus, features = gene )+scale_color_gradientn(colours=mycols_b)
print(plot1)
dev.off()}

#### Violin plots ####
for (gene in int_genes){
  pdf(paste(dato,"_thymusSC_Violin_GeneExpSplitByColonisation_",gene,".pdf",sep = ""), height = 4,width = 8)
plot2 <- VlnPlot(thymus, group.by = "res.0.2", split.by = "colonization", 
        features = gene, pt.size = 0)+
  scale_fill_manual(values=volc_cols)
print(plot2)
dev.off()}

#### Heatmap ####



