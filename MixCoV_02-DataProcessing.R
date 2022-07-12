#!/usr/bin/env Rscript
library(stringr)
library(plyr)
library(dplyr)
library(Seurat)
library(SeuratWrappers)
library(glmGamPoi)
library(edgeR)
library(limma)
library(ggrepel)
library(ggrastr)
library(ggplot2)
library(ggridges)
library(patchwork)
library(pheatmap)
library(gprofiler2)
library(reshape2)
source('MixCoV_00-PlotFunctions.R')


## Load the PBMC dataset
data <- Read10X(data.dir = "PanCoV_scRNA_Agg/outs/count/raw_feature_bc_matrix")


## Initialize the Seurat object
so <- CreateSeuratObject(counts = data, project = "PanCoV", min.cells = 3, min.features = 200)


## Annotate dataset
sample <- c('pbs_6','pbs_7','pbs_8','delta_6','delta_7','delta_8','panLow_6','panLow_7','panLow_8','panLow_9','panHigh_6','panHigh_7','panHigh_8','panHigh_9')
treatment <- c(rep('pbs',3), rep('delta',3), rep('panLow', 4), rep('panHigh', 4)) 
lane <- c(rep('L1',10), rep('L2',4))
so[['sample']] <- sample[stringr::str_split(colnames(so), '-', simplify = T)[,2] %>% as.integer()]
so[['treatment']] <- treatment[stringr::str_split(colnames(so), '-', simplify = T)[,2] %>% as.integer()]
so[['lane']] <- lane[stringr::str_split(colnames(so), '-', simplify = T)[,2] %>% as.integer()]
so[['cell.detection.rate']] <- scale(so@meta.data$nFeature_RNA)
so[['percent.mt']] <- PercentageFeatureSet(so, pattern = '^mt-')


## Visualize QC metrics
VlnPlot(so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
qc_plot <- FeatureScatter(so, feature1 = "nCount_RNA", feature2 = "percent.mt") 
FeatureScatter(so, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave(qc_plot, height = 6, width = 6, filename = imageName)


## Save unfiltered data
saveRDS(so, file = 'PanCoV_SeuratData_Unfiltered.RDS')
so <- readRDS(file = 'PanCoV_SeuratData_Unfiltered.RDS')


## Filter data by feature number and mito %
so <- subset(so, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & percent.mt < 5)


## Preprocess data by sample
so.list <- SplitObject(so, split.by = "sample")
so.list <- lapply(so.list, function(x){
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000) 
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(so.list)
so.list <- lapply(so.list, function(x){
  x <- ScaleData(x, features = features) 
  x <- RunPCA(x, features = features)
})


## Integrate samples
immune.anchors <- FindIntegrationAnchors(object.list = so.list, anchor.features = features, reduction = 'rpca', k.anchor = 20)
so.int <- IntegrateData(anchorset = immune.anchors)
DefaultAssay(so.int) <- "integrated"


## Save integrated samples
saveRDS(so.int, file = 'PanCoV_SeuratData_vstIntegrated.RDS')

