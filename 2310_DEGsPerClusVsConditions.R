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
library(ggrepel)

###### variables used through script ######
rm(list=ls())
#date in format year_month_day
dato <- str_sub(str_replace_all(Sys.Date(),"-","_"), 3, -1)
# colour string for imputation and overlays
mycols_b <- c("blue","#bdbdbd","#d9d9d9","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F")
volc_cols <- c("gf"= "#F8766D","sdmdmm2"="#00BA38","spf"="#619CFF","NA"="lightgrey")
# project
proj_dir <- "/Users/linewulff/Documents/work/projects/23_Geuking_thymus_SC"
out_dir <- paste(proj_dir,"ColonizationDEGs",sep = '/')

setwd(out_dir)
##### Read in Seurat object #####
# from initial integration adn ID
thymus <- readRDS("2309_Geuking_ThymusSC_allcond.rds")

#### Comaprison between conditions ####
DefaultAssay(thymus) <- "RNA"
for (clus in unique(thymus@meta.data$res.0.2)){
  print(clus)
  #reset idents for next round of loop
  Idents(thymus) <- 'res.0.2'
  # subset per cluster
  cell_sub <- subset(thymus, idents = clus)
  Idents(cell_sub) <- 'colonization'
  
  ### "sdmdmm2" vs "gf"
  print("sdmdmm2 vs gf")
  if (nrow(cell_sub@meta.data[cell_sub@meta.data$colonization=="sdmdmm2",]) < 3 | nrow(cell_sub@meta.data[cell_sub@meta.data$colonization=="gf",]) < 3){print("Too few cells in one group")
    } else {
  # Run DEGA
  DEGs <- FindMarkers(cell_sub, ident.1 = "sdmdmm2", ident.2 = "gf", 
                      logfc.threshold = 0, # slower but necessary for pretty volcano plots
                      test.use = "wilcox")
  DEGs$gene <- rownames(DEGs)
  # significance limits as IDs
  DEGs$DEG <- "NA"
  if (!identical(DEGs[DEGs$p_val_adj<0.05 & DEGs$avg_log2FC>0.25,]$DEG,character(0))){DEGs[DEGs$p_val_adj<0.05 & DEGs$avg_log2FC>0.25,]$DEG <- "sdmdmm2"}
  if (!identical(DEGs[DEGs$p_val_adj<0.05 & DEGs$avg_log2FC<(-0.25),]$DEG,character(0))){DEGs[DEGs$p_val_adj<0.05 & DEGs$avg_log2FC<(-0.25),]$DEG <- "gf"}
  
  sign <- DEGs[DEGs$DEG!="NA",]
  # top ten markers by avg logFC - if available
  if (length(unique(DEGs$DEG))>1){
    if (nrow(sign[sign$DEG=="sdmdmm2",])>=10 & nrow(sign[sign$DEG=="gf",])>=10){
    top_10 <- sign[order(sign$avg_log2FC, decreasing = T),][c(1:10,(dim(sign)[1]-9):dim(sign)[1]),]
    } else if (nrow(sign[sign$DEG=="sdmdmm2",])<10 & nrow(sign[sign$DEG=="gf",])>=10){
      top_10 <- sign[sign$DEG=="sdmdmm2",]
      top_10 <- top_10 %>% rbind(sign[order(sign$avg_log2FC, decreasing = T),][c((dim(sign)[1]-9):dim(sign)[1]),])
    } else if (nrow(sign[sign$DEG=="sdmdmm2",])>=10 & nrow(sign[sign$DEG=="gf",])<10){
      top_10 <- sign[sign$DEG=="gf",]
      top_10 <- top_10 %>% rbind(sign[order(sign$avg_log2FC, decreasing = T),][c(1:10),])
    } else {
      top_10 <- sign[sign$DEG=="gf",]
      top_10 <- top_10 %>% rbind(sign[sign$DEG=="sdmdmm2",])
    }} else {top_10 <- DEGs[DEGs$DEG!="NA",]}
  
  # plot volcano plot w. top 10 DEGs
  volc <- ggplot(DEGs, aes(x = avg_log2FC, y = -log10(p_val_adj), colour = DEG))+
    geom_point()+
    geom_hline(yintercept = -log10(0.05), linetype = "dashed")+ #sign.
    geom_vline(xintercept = c(-0.25,0.25), linetype = "dashed")+ #sign.
    geom_text_repel(data = top_10,aes(x = avg_log2FC, y = -log10(p_val_adj)), label = top_10$gene, colour="black")+
    scale_colour_manual(values = volc_cols)+
    theme_classic()
  pdf(paste(dato,"volcano",clus,"sdmdmm2vsGF_logFC0.25padj0.05.pdf", sep = "_"), height = 4, width = 5.5)
  print(volc)
  dev.off()
  
  # Subset significant and save csv
  write.csv(DEGs[DEGs$DEG!="NA",], file = paste(dato,"DEG",clus,"sdmdmm2vsGF_logFC0.25padj0.05.csv", sep = "_"))
  }
  
  ### "spf" vs "gf"
  print("spf vs gf")
  if (nrow(cell_sub@meta.data[cell_sub@meta.data$colonization=="spf",]) < 3 | nrow(cell_sub@meta.data[cell_sub@meta.data$colonization=="gf",]) < 3){print("Too few cells in one group")
  } else {
  # Run DEGA
  DEGs <- FindMarkers(cell_sub, ident.1 = "spf", ident.2 = "gf", 
                      logfc.threshold = 0, # slower but necessary for pretty volcano plots
                      test.use = "wilcox")
  DEGs$gene <- rownames(DEGs)
  # significance limits as IDs
  DEGs$DEG <- "NA"
  if (!identical(DEGs[DEGs$p_val_adj<0.05 & DEGs$avg_log2FC>0.25,]$DEG,character(0))){DEGs[DEGs$p_val_adj<0.05 & DEGs$avg_log2FC>0.25,]$DEG <- "spf"}
  if (!identical(DEGs[DEGs$p_val_adj<0.05 & DEGs$avg_log2FC<(-0.25),]$DEG,character(0))){DEGs[DEGs$p_val_adj<0.05 & DEGs$avg_log2FC<(-0.25),]$DEG <- "gf"}
  
  sign <- DEGs[DEGs$DEG!="NA",]
  # top ten markers by avg logFC - if available
  if (length(unique(DEGs$DEG))>1){
    if (nrow(sign[sign$DEG=="sdmdmm2",])>=10 & nrow(sign[sign$DEG=="gf",])>=10){
      top_10 <- sign[order(sign$avg_log2FC, decreasing = T),][c(1:10,(dim(sign)[1]-9):dim(sign)[1]),]
    } else if (nrow(sign[sign$DEG=="sdmdmm2",])<10 & nrow(sign[sign$DEG=="gf",])>=10){
      top_10 <- sign[sign$DEG=="sdmdmm2",]
      top_10 <- top_10 %>% rbind(sign[order(sign$avg_log2FC, decreasing = T),][c((dim(sign)[1]-9):dim(sign)[1]),])
    } else if (nrow(sign[sign$DEG=="sdmdmm2",])>=10 & nrow(sign[sign$DEG=="gf",])<10){
      top_10 <- sign[sign$DEG=="gf",]
      top_10 <- top_10 %>% rbind(sign[order(sign$avg_log2FC, decreasing = T),][c(1:10),])
    } else {
      top_10 <- sign[sign$DEG=="gf",]
      top_10 <- top_10 %>% rbind(sign[sign$DEG=="sdmdmm2",])
    }} else {top_10 <- DEGs[DEGs$DEG!="NA",]}
  
  # plot volcano plot w. top 10 DEGs
  volc <- ggplot(DEGs, aes(x = avg_log2FC, y = -log10(p_val_adj), colour = DEG))+
    geom_point()+
    geom_hline(yintercept = -log10(0.05), linetype = "dashed")+ #sign.
    geom_vline(xintercept = c(-0.25,0.25), linetype = "dashed")+ #sign.
    geom_text_repel(data = top_10,aes(x = avg_log2FC, y = -log10(p_val_adj)), label = top_10$gene, colour="black")+
    scale_colour_manual(values = volc_cols)+
    theme_classic()
  pdf(paste(dato,"volcano",clus,"spfvsGF_logFC0.25padj0.05.pdf", sep = "_"), height = 4, width = 5.5)
  print(volc)
  dev.off()
  
  # Subset significant and save csv
  write.csv(DEGs[DEGs$DEG!="NA",], file = paste(dato,"DEG",clus,"spfvsGF_logFC0.25padj0.05.csv", sep = "_"))
  }
  
  ### "sdmdmm2" vs "spf"
  print("sdmdmm2 vs spf")
  if (nrow(cell_sub@meta.data[cell_sub@meta.data$colonization=="sdmdmm2",]) < 3 | nrow(cell_sub@meta.data[cell_sub@meta.data$colonization=="spf",]) < 3){print("Too few cells in one group")
  } else {
  # Run DEGA
  DEGs <- FindMarkers(cell_sub, ident.1 = "sdmdmm2", ident.2 = "spf", 
                      logfc.threshold = 0, # slower but necessary for pretty volcano plots
                      test.use = "wilcox")
  DEGs$gene <- rownames(DEGs)
  # significance limits as IDs
  DEGs$DEG <- "NA"
  if (!identical(DEGs[DEGs$p_val_adj<0.05 & DEGs$avg_log2FC>0.25,]$DEG,character(0))){DEGs[DEGs$p_val_adj<0.05 & DEGs$avg_log2FC>0.25,]$DEG <- "sdmdmm2"}
  if (!identical(DEGs[DEGs$p_val_adj<0.05 & DEGs$avg_log2FC<(-0.25),]$DEG,character(0))){DEGs[DEGs$p_val_adj<0.05 & DEGs$avg_log2FC<(-0.25),]$DEG <- "spf"}
  
  sign <- DEGs[DEGs$DEG!="NA",]
  # top ten markers by avg logFC - if available
  if (length(unique(DEGs$DEG))>1){
    if (nrow(sign[sign$DEG=="sdmdmm2",])>=10 & nrow(sign[sign$DEG=="gf",])>=10){
      top_10 <- sign[order(sign$avg_log2FC, decreasing = T),][c(1:10,(dim(sign)[1]-9):dim(sign)[1]),]
    } else if (nrow(sign[sign$DEG=="sdmdmm2",])<10 & nrow(sign[sign$DEG=="gf",])>=10){
      top_10 <- sign[sign$DEG=="sdmdmm2",]
      top_10 <- top_10 %>% rbind(sign[order(sign$avg_log2FC, decreasing = T),][c((dim(sign)[1]-9):dim(sign)[1]),])
    } else if (nrow(sign[sign$DEG=="sdmdmm2",])>=10 & nrow(sign[sign$DEG=="gf",])<10){
      top_10 <- sign[sign$DEG=="gf",]
      top_10 <- top_10 %>% rbind(sign[order(sign$avg_log2FC, decreasing = T),][c(1:10),])
    } else {
      top_10 <- sign[sign$DEG=="gf",]
      top_10 <- top_10 %>% rbind(sign[sign$DEG=="sdmdmm2",])
    }} else {top_10 <- DEGs[DEGs$DEG!="NA",]}
  
  # plot volcano plot w. top 10 DEGs
  volc <- ggplot(DEGs, aes(x = avg_log2FC, y = -log10(p_val_adj), colour = DEG))+
    geom_point()+
    geom_hline(yintercept = -log10(0.05), linetype = "dashed")+ #sign.
    geom_vline(xintercept = c(-0.25,0.25), linetype = "dashed")+ #sign.
    geom_text_repel(data = top_10,aes(x = avg_log2FC, y = -log10(p_val_adj)), label = top_10$gene, colour="black")+
    scale_colour_manual(values = volc_cols)+
    theme_classic()
  pdf(paste(dato,"volcano",clus,"sdmdmm2vsSpf_logFC0.25padj0.05.pdf", sep = "_"), height = 4, width = 5.5)
  print(volc)
  dev.off()
  
  # Subset significant and save csv
  write.csv(DEGs[DEGs$DEG!="NA",], file = paste(dato,"DEG",clus,"sdmdmm2vsSpf_lospfC0.25padj0.05.csv", sep = "_"))
  }
  }




