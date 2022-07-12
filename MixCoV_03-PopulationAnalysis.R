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


## Load integrated data
so.int <- readRDS(file = 'PanCoV_SeuratData_vstIntegrated.RDS')


## Linear dimensional reduction
so.int <- ScaleData(so.int)
so.int <- RunPCA(so.int, npcs = 40)
ElbowPlot(so.int, ndims = 30, reduction = 'pca')


## UMAP dimensional reduction
DefaultAssay(so.int) <- 'integrated'
so.int <- RunUMAP(so.int, reduction = 'pca', dims = 1:12)
so.int <- FindNeighbors(so.int, reduction = 'pca', dims = 1:12)
so.int <- FindClusters(so.int, graph.name = 'integrated_snn', algorithm = 2, resolution = 0.31, random.seed = 1101)


## Normalize and Scale
DefaultAssay(so.int) <- 'RNA'
so.int <- NormalizeData(so.int, normalization.method = "LogNormalize", scale.factor = 10000) 
so.int <- FindVariableFeatures(so.int, selection.method = "vst", nfeatures = 2000)
so.int <- ScaleData(so.int, features = rownames(so.int))
DefaultAssay(so.int) <- 'integrated'


## Set population-specific marker list
marker.list <- list(Tcell = c('Cd3d','Cd8b1','Cd4','Ccr7','Cd44','Tcf7'), 
  CD4_Tcell = c('Cd3d','Cd4','Tbx21','Gata3','Foxp3','Il2ra','Gzmb','Ifng','Tgfb1','Il4'), 
  Plasma = c('Cd19','Ms4a1','Cr2','Cd9','Sdc1','H2-Aa','H2-Eb2','Prdm1','Pax5','Mki67','Irf4','Bach2'), 
  Bcell = c('Cd19','Ms4a1','Cr2','Cd27','Spn','H2-Eb2','H2-Aa','Ighm','Igha','Ighd','Cd24a','Fcer2a', 'Cd38', 'Pdcd1lg2', 'Cd80', 'Nt5e'), 
  Myeloid = c('Cd14','Fcgr3', 'Itgam', 'Itgax', 'Ncr1'),   
  Mono_Macro_DC_subsets = c('Itgam','Itgax','Cd24a','Bst2','Cd8b1','Csf1r','Cd14','Cd68','Adgre1','Ly6c1','Fcgr3','Fcgr1','Itgae','Sirpa', 'Siglech', 'Ccr2','Cx3cr1'))


## Compare marker expression of major immune cell types between clusters
GetDetectionMatrix <- function(markers, seurat.object = so.int, clusters = 'all'){
  if(clusters == 'all') clusters <- Idents(seurat.object) %>% sort() %>% unique() %>% as.character()
  cells <- CellsByIdentities(seurat.object, idents = clusters)
  data <- seurat.object$integrated@scale.data[markers, unlist(cells)] 
  m <- matrix(ncol = length(clusters), nrow = length(markers), dimnames = list(markers,clusters))
  for(marker in markers) m[marker,] <- sapply(1:length(cells), function(i) {(sum(data[marker, cells[[i]]] > .2))/length(cells[[i]])})
  return(m)
}
hm.data <- lapply(names(cluster.list), function(x) {GetDetectionMatrix(cluster.list[[x]], marker.list[[x]])})
for(i in 1:3) pheatmap(hm.data[[i]], scale = 'none', clustering_method = 'ward.D2',treeheight_row = 10, treeheight_col = 10, main = paste0(names(cluster.list)[i], ' cluster det rates'), cluster_rows = F,filename = paste0('umap_cluster_majorMarkers_', names(cluster.list)[i], '.pdf'), width = 3, height = 3)


## Subcluster populations 3 and 4
so <- FindSubCluster(so, cluster = '3', graph.name = 'integrated_snn', subcluster.name = 'intermediate', algorithm = 2, resolution = 0.08)
Idents(so) <- so[['intermediate']]
so <- FindSubCluster(so, cluster = '4', graph.name = 'integrated_snn', subcluster.name = 'intermediate', algorithm = 2, resolution = 0.1)
Idents(so) <- so[['intermediate']]
so[['celltype']] <- Idents(so)


## Label major immune cell types based on marker expression
so <- RenameIdents(so, '1' = 'Naive CD4 T', '2' = 'Naive CD8 T', '5' = 'CD8 T Effector', '8' = 'Activated CD4 T cells', '10' = 'CD8 TEM', '11' = 'Tcell-Bcell', '0' = 'Bcell', '3_0' = 'CD8 TCM', '3_1' = 'Bcell', '3_2' = 'Bcell', '4_0' = 'Bcell','4_1' = 'Naive CD8 T', '7' = 'Bcell', '13' = 'Plasma', '14' = 'Myeloid', '15' = 'Bcell','6' = 'NK cell', '9' = 'Myeloid', '12' = 'Myeloid', '16' = 'Plasma', '17' = 'Myeloid')
so[['celltype']] <- Idents(so)


## Save labeled data
saveRDS(so.int, file = 'PanCoV_SeuratData_vstIntegrated_v2.RDS')


