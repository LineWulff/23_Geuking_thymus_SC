#' R script for 10x preprocessing
#' Author: Line Wulff
#' Date (created): 23-09-22


###### Libraries ######
library(gplots)
library(dplyr)
library(stringr)
library(scales)
#library(Rmagic)
library(ggplot2)
library(viridis)
#library(clustree)
library(ccRemover)
library(rgl)
library(Seurat)

###### variables used through script ######
rm(list=ls())
#date in format year_month_day
dato <- str_sub(str_replace_all(Sys.Date(),"-","_"), 3, -1)
# colour string for imputation and overlays
mycols_b <- c("#bdbdbd","#d9d9d9","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F")
# project
project <- "s7"
proj_dir <- "/Users/linewulff/Documents/work/projects/23_Geuking_thymus_SC"
data_dir <- "/Volumes/Promise\ RAID/Markus/Kirsty_scRNA_May2023"

if (!dir.exists(paste(proj_dir,project,sep="/"))){
  dir.create(paste(proj_dir,project,sep="/"))
}else{
  print("dir exists")
}

setwd(paste(proj_dir,project,sep="/"))

#### read in data ####
pat_data <- Read10X(data.dir = paste(data_dir,project,"filtered_feature_bc_matrix/", sep="/"))

###### As Seurat objects ######
pat_data <- CreateSeuratObject(counts = pat_data, project = project, min.cells = 3, min.features = 100)

###### ----------------- Quality Control part 1 ----------------- ######
# MT genes for mice using Broad institutes Mito carta data base as ref.
# see https://www.broadinstitute.org/mitocarta/mitocarta30-inventory-mammalian-mitochondrial-proteins-and-pathways
# for citations etc.
MTref <- read.csv("/Volumes/Promise\ RAID/Line/reference_data/Mouse.MitoCarta3.0.csv", header = T)$Symbol
MTref <- rownames(pat_data@assays$RNA@counts)[rownames(pat_data@assays$RNA@counts) %in% MTref]
pat_data[["percent.mt"]] <- PercentageFeatureSet(pat_data, features = MTref)

#### cell cycle 
data("cc.genes") ## Tirosh et al. 2015 HUMAN GENES, identify same in mice
sgenes <- rownames(pat_data@assays$RNA@counts)[toupper(rownames(pat_data@assays$RNA@counts)) %in% cc.genes$s.genes]
g2mgenes <- rownames(pat_data@assays$RNA@counts)[toupper(rownames(pat_data@assays$RNA@counts)) %in% cc.genes$g2m.genes]

pat_data <- CellCycleScoring(pat_data, s.features = sgenes, g2m.features = g2mgenes, set.ident = FALSE)


###### ----------------- QC plots ----------------- ######
dir.create(paste(proj_dir,project,"QC",sep="/")) 

setwd(paste(proj_dir,project,"QC",sep="/"))
# nFeature = number of genes per barcode
# nCount = number of reads per barcode = count depth

#png(paste(dato,project,"QC_measures1.png",sep="_"),width = 1200, height = 1000, res = 150)
VlnPlot(pat_data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1)
#dev.off()

pat_data_QC <- as.data.frame(pat_data@meta.data)
pat_data_QC <- pat_data_QC[order(pat_data_QC$nCount_RNA),]
pat_data_QC <- cbind(pat_data_QC, rank = c(1:length(pat_data_QC$nCount_RNA)))

## Count depth
#red line = zoom for next plot
png(paste(dato,project,"QC_CountDepthHisto_1.png",sep = "_"),height = 800, width = 1000, res =150)
ggplot(pat_data_QC, aes(x=nCount_RNA))+geom_histogram(fill="grey",bins = 100)+geom_vline(xintercept = 6000, colour = "red")+
  xlab("Count Depth")+ggtitle(paste(project,"_ILF",sep=""))
dev.off()

# change count depth cut off suggestion
png(paste(dato,project,"QC_CountDepthHisto_2.png",sep = "_"),height = 800, width = 1000, res =150)
ggplot(pat_data_QC[pat_data_QC$nCount_RNA<6000,], aes(x=nCount_RNA))+geom_histogram(fill="grey", bins = 30)+
  geom_vline(xintercept = 2500, colour = "red")+xlab("Count Depth")+ggtitle(paste(project,"_ILF",sep=""))
dev.off()

## Number of genes
png(paste(dato,project,"QC_NumberGeneshHisto_1.png",sep = "_"),height = 800, width = 1000, res =150)
ggplot(pat_data_QC, aes(x=nFeature_RNA))+geom_histogram(fill="grey",bins = 100)+
  xlab("Number of genes")+ggtitle(paste(project,"_ILF",sep=""))+geom_vline(xintercept = c(1200,5500), colour = "red")
dev.off()

png(paste(dato,project,"QC_NumberGeneshHisto_2.png",sep = "_"),height = 800, width = 1000, res =150)
ggplot(pat_data_QC[pat_data_QC$nFeature_RNA<2000,], aes(x=nFeature_RNA))+geom_histogram(fill="grey",bins = 100)+
  xlab("Number of genes")+ggtitle(paste(project,"_ILF",sep=""))+geom_vline(xintercept = 1200, colour = "red")
dev.off()

png(paste(dato,project,"QC_CountDepthRanked.png",sep = "_"),height = 800, width = 1000, res =150)
ggplot(pat_data_QC, aes(x=-rank, y=nCount_RNA))+geom_point()+
  xlab("Barcode rank")+ylab("Count Depth")+ggtitle(paste(project,sep=""))+
  geom_hline(yintercept = 2500, colour = "red")
dev.off()

png(paste(dato,project,"QC_CountDepth_GeneCount_1.png",sep = "_"),height = 800, width = 1000, res =150)
ggplot(pat_data_QC, aes(x=nCount_RNA, y=nFeature_RNA, colour=percent.mt))+geom_point()+
  scale_colour_viridis()+xlab("Count Depth")+ylab("Number of genes")+ggtitle(paste(project,sep=""))+
  geom_vline(xintercept = c(2500,25000), colour = "red")+geom_hline(yintercept = c(1200,5500),colour="red")
dev.off()

# crazy high MT? plot this one too
png(paste(dato,project,"QC_CountDepth_GeneCount_1.png",sep = "_"),height = 800, width = 1000, res =150)
ggplot(pat_data_QC, aes(x=percent.mt, y=nFeature_RNA, colour=percent.mt))+geom_point()+
     scale_colour_viridis()+ylab("Number of genes")+ggtitle(paste(project,sep=""))+
  geom_hline(yintercept = c(1200,5500),colour="red")+
  geom_vline(xintercept = c(20),colour="red")
dev.off()

png(paste(dato,project,"QC_CountDepth_GeneCount_2.png",sep = "_"),height = 800, width = 1000, res =150)
ggplot(pat_data_QC[pat_data_QC$nCount_RNA<6000,], aes(x=nCount_RNA, y=nFeature_RNA, colour=percent.mt))+geom_point(size=0.2)+
  scale_colour_viridis()+xlab("Count Depth")+ylab("Number of genes")+ggtitle(paste(project,sep=""))+
  geom_hline(yintercept = 600, colour = "red")
dev.off()

#Set cut offs and subset slightly lower than the lines (can be adjusted later)
MTthres <- 20
CountLow <- NA #include in subsetting if necessary
CountHigh <- NA #include in subsetting if necessary
FeatLow <- 1200
FeatHigh <- 5500

# save some numbers on orig. sample
Pre_ncells <- length(colnames(pat_data))
### Set thresholds and subset data
pat_data <- subset(pat_data, subset = nFeature_RNA > FeatLow & nFeature_RNA < FeatHigh
                   & percent.mt < MTthres)
Post_ncells <- length(colnames(pat_data))
rem_nperc <- round((Pre_ncells-Post_ncells)/Pre_ncells*100)

ThresFile <- file(paste(dato,"SubsetThresholds",project,".txt",sep="_"))
writeLines(c(paste("Pre subsetting there were",Pre_ncells,"cells."),
             paste("Removing ~",rem_nperc,"% of cells during single cell QC."),
             paste("Post subsetting there are",Post_ncells,"cells"),
             "Thresholds were set to:",
             paste("Gene/feature level:",FeatLow,"to",FeatHigh),
             paste("Count/read level:",CountLow,"to",CountHigh),
             paste("MT threshold was set to: <",MTthres,"%")), 
           ThresFile)
close(ThresFile)

### Normalize data
pat_data <- NormalizeData(pat_data)

### save objects separately
pat_data <- FindVariableFeatures(pat_data, selection.method = "vst", nfeatures = 3000)
top10 <- head(VariableFeatures(pat_data), 10)

pat_data <- RenameCells(pat_data, add.cell.id = project)

setwd(paste(proj_dir,project,sep="/"))
saveRDS(object=pat_data, file=paste(project,".rds",sep=""))

plot1 <- VariableFeaturePlot(pat_data)
#png(paste(dato,project,"TopVariableeGenes.png",sep = "_"),height = 1000,width = 1000, res =150)
LabelPoints(plot = plot1, points = top10, repel = TRUE)
#dev.off()
