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


## Function compares gene detection rate across groups of cells (proportion of cells with scaled gene expr > 0.2) 
GetDetectionMatrix <- function(markers, seurat.object = so.int, clusters = 'all'){
  if(clusters == 'all') clusters <- Idents(seurat.object) %>% sort() %>% unique() %>% as.character()
  cells <- CellsByIdentities(seurat.object, idents = clusters)
  data <- seurat.object$integrated@scale.data[markers, unlist(cells)] 
  m <- matrix(ncol = length(clusters), nrow = length(markers), dimnames = list(markers,clusters))
  for(marker in markers) m[marker,] <- sapply(1:length(cells), function(i) {(sum(data[marker, cells[[i]]] > .2))/length(cells[[i]])})
  return(m)
}


## Load integrated data
so.int <- readRDS(file = 'PanCoV_SeuratData_vstIntegrated_v2.RDS')


## Set population-specific marker list
marker.list <- list(Tcell = c('Cd3d','Cd8b1','Cd4','Ccr7','Cd44','Tcf7'), 
  CD4_Tcell = c('Cd3d','Cd4','Tbx21','Gata3','Foxp3','Il2ra','Gzmb','Ifng','Tgfb1','Il4'), 
  Plasma = c('Cd19','Ms4a1','Cr2','Cd9','Sdc1','H2-Aa','H2-Eb2','Prdm1','Pax5','Mki67','Irf4','Bach2'), 
  Bcell = c('Cd19','Ms4a1','Cr2','Cd22','Fcer2a','Cd24a','Cd27','Cd38','Cd40','Spn','Cd69','Cd80','Cd83','Cd86','Cd93','Fas','H2-Aa','H2-Eb2', 'Irf4','Bcl6','Pax5','Prdm1','Spib','Ighm','Ighd','Igha','Ighg1','Ighg2b','Ighg3','Pdcd1lg2','Nt5e','Sdc1','Cd9'),
  Myeloid = c('Cd14','Fcgr3', 'Itgam', 'Itgax', 'Ncr1'),   
  Mono_Macro_DC_subsets = c('Itgam','Itgax','Cd24a','Bst2','Cd8b1','Csf1r','Cd14','Cd68','Adgre1','Ly6c1','Fcgr3','Fcgr1','Itgae','Sirpa', 'Siglech', 'Ccr2','Cx3cr1'))


## Get Myeloid subsets
so.sub <- SplitObject(so.int, split.by = 'cell.type')
so.sub <- so.sub[['Myeloid']] 
so.sub <- ScaleData(so.sub, vars.to.regress = 'percent.mt') 
so.sub <- RunPCA(so.sub, npcs = 40)
ElbowPlot(so.sub, ndims = 30, reduction = 'pca')


## visualize myeloid subclusters
so.sub <- RunUMAP(so.sub, reduction = 'pca', dims = 1:10)
so.sub <- FindNeighbors(so.sub, features = marker.list[['Mono_Macro_DC_subsets']], graph.name = 'markers')
so.sub <- FindClusters(so.sub, algorithm = 3, graph.name = 'markers', resolution = 0.1, random.seed = 1101)
DimPlot(so.sub, label = T, label.box = T)


## Identify myeloid subclusters
GetDetectionMatrix(as.character(0:6), marker.list[['Mono_Macro_DC_subsets']], so.sub)


## Relabel and visualize myeloid subclusters
so.sub <- RenameIdents(so.sub, '0' = 'Monocyte', '1' = 'cDC2', '2' = 'cDC1', '3' = 'Macrophage', '4' = 'pDC','5' = 'Macrophage', '6' = 'cDC1')
DimPlot(so.sub, label = T, label.box = T, repel = T) + theme(legend.position = "none")
ggsave(file = 'PanCoV_SeuratData_vstIntegrated_MyeloidSubcluster_UMAP.pdf', unit = 'in', width = 4, height = 3)
DotPlot(so.sub, features = marker.list[['Plasma']]) + coord_flip() + theme(axis.text.x = element_text(angle=90, hjust=1, size = 10)) 
ggsave(file = 'SubclusterMyeloid_Dotplot.pdf', unit = 'in', width = 4, height = 4)


## Save myeloid data
saveRDS(so.sub, file = 'PanCoV_SeuratData_vstIntegrated_MyeloidSubcluster.RDS')


## Rename myeloid clusters in so.int
cells <- CellsByIdentities(so.sub)
Idents(so.int) <- so.int[['cell.type']]
for(celltype in names(cells)) Idents(so.int, cells = cells[[celltype]]) <- celltype 
so.int[['cell.type']] <- Idents(so.int)


#####
## Get CD4 T cell subsets
so.sub <- SplitObject(so.int, split.by = 'cell.type')
so.sub <- so.sub[['Activated CD4 T cells']]
so.sub <- ScaleData(so.sub, vars.to.regress = 'percent.mt') 
so.sub <- RunPCA(so.sub, npcs = 40)
ElbowPlot(so.sub, ndims = 30, reduction = 'pca')


## Visualize CD4 T cell subclusters
so.sub <- RunUMAP(so.sub, reduction = 'pca', dims = 1:16)
so.sub <- FindNeighbors(so.sub, features = marker.list[['CD4_Tcell']], graph.name = 'markers')
so.sub <- FindClusters(so.sub, algorithm = 3, graph.name = 'markers', resolution = 0.25, random.seed = 3)
DimPlot(so.sub, label = T, label.box = T)


## Identify CD4 T cell subclusters
GetDetectionMatrix(as.character(0:(max(as.integer(Idents(so.sub))) - 1)), marker.list[['CD4_Tcell']], so.sub)

## Relabel and visualize CD4 T cell subclusters
so.sub <- RenameIdents(so.sub, '0' = 'Treg', '1' = 'Th2', '2' = 'Th1', '3' = 'Th2', '4' = 'Treg','5' = 'Th2')
DimPlot(so.sub, label = T, label.box = T, repel = T) + theme(legend.position = "none")
ggsave(file = 'PanCoV_SeuratData_vstIntegrated_CD4TcellSubcluster_UMAP.pdf', unit = 'in', width = 4, height = 3)
DotPlot(so.sub, features = marker.list[['CD4_Tcell']]) + coord_flip() + theme(axis.text.x = element_text(angle=90, hjust=1)) 
ggsave(file = 'SubclusterCd4Tcell_Dotplot.pdf', unit = 'in', width = 4, height = 3)


## Save CD4 T cell data
saveRDS(so.sub, file = 'PanCoV_SeuratData_vstIntegrated_CD4TcellSubcluster.RDS')


## Rename CD4 T cell clusters in so.int
cells <- CellsByIdentities(so.sub)
for(celltype in names(cells)) Idents(so.int, cells = cells[[celltype]]) <- celltype 
so.int[['cell.type']] <- Idents(so.int)


#####
## Get B cell subsets
so.sub <- SplitObject(so.int, split.by = 'cell.type')
so.sub <- so.sub[['Bcell']]
so.sub <- ScaleData(so.sub, vars.to.regress = 'percent.mt') 
so.sub <- RunPCA(so.sub, npcs = 40)
ElbowPlot(so.sub, ndims = 30, reduction = 'pca')


## Visualize B cell subclusters
so.sub <- RunUMAP(so.sub, reduction = 'pca', dims = 1:14)
so.sub <- FindNeighbors(so.sub, features = marker.list[['Bcell']], graph.name ='markers')
so.sub <- FindClusters(so.sub, algorithm = 3, graph.name = 'markers', resolution = 0.2, random.seed = 1101)
Idents(so.sub) <- so.sub[['markers_res.0.2']]


## Identify B cell subclusters
GetDetectionMatrix(as.character(0:(max(as.integer(Idents(so.sub))) - 1)), marker.list[['Bcell']], so.sub)


## Relabel and visualize B cell subclusters
so.sub <- RenameIdents(so.sub, '0' = 'Activated B cell', '1' = 'Naive B cell', '2' = 'Unswitched Memory B cell', '3' = 'Naive B cell', '4' = 'Switched Memory B cell', '5' = 'Plasmablast','6' = 'Activated B cell','7' = 'Activated B cell')
DimPlot(so.sub, label = T, label.box = T, repel = T) + theme(legend.position = "none")
ggsave(file = 'PanCoV_SeuratData_vstIntegrated_BcellSubcluster_UMAP_v2.pdf', unit = 'in', width = 4, height = 3)
DotPlot(so.sub, features = marker.subset[['Bcell']],cols = c('gray80','#EB8335'), cluster.idents = F)  + theme(axis.text.x = element_text(angle=90, hjust=1, size = 10)) 
ggsave(file = 'SubclusterBcell_Dotplot_v2.pdf', unit = 'in', width = 7, height = 2.6)


## Save B cell data
saveRDS(so.sub, file = 'PanCoV_SeuratData_vstIntegrated_BcellSubcluster_v2.RDS')


## Rename B cell clusters in so.int
cells <- CellsByIdentities(so.sub)
for(celltype in names(cells)) Idents(so.int, cells = cells[[celltype]]) <- celltype 
so.int[['cell.type']] <- Idents(so.int)


#####
## Get Plasma subsets
so.sub <- SplitObject(so.int, split.by = 'cell.type')
so.sub <- so.sub[['Plasma']]
so.sub <- ScaleData(so.sub, vars.to.regress = 'percent.mt') 
so.sub <- RunPCA(so.sub, npcs = 40)
ElbowPlot(so.sub, ndims = 30, reduction = 'pca')


## Visualize Plasma subclusters
so.sub <- RunUMAP(so.sub, reduction = 'pca', dims = 1:11)
so.sub <- FindNeighbors(so.sub, features = marker.list[['Plasma']], graph.name = 'markers')
so.sub <- FindClusters(so.sub, algorithm = 3, graph.name = 'markers', resolution = 0.1, random.seed = 1101)
DimPlot(so.sub, label = T, label.box = T)
DotPlot(so.sub, features = marker.list[['Plasma']]) + coord_flip() #+ theme(axis.text.x = element_text(angle=90, hjust=1)) 


## Identify Plasma subclusters
GetDetectionMatrix(as.character(0:(max(as.integer(Idents(so.sub))) - 1)), marker.list[['Plasma']], so.sub)


## Relabel and visualize Plasma subclusters
so.sub <- RenameIdents(so.sub, '0' = 'Plasmablast', '1' = 'Pre-plasmablast', '2' = 'Plasma cell', '3' = 'Pre-plasmablast')
DimPlot(so.sub, label = T, label.box = T, repel = T) + theme(legend.position = "none")
ggsave(file = 'PanCoV_SeuratData_vstIntegrated_PlasmaSubcluster_UMAP.pdf', unit = 'in', width = 4, height = 3)
DotPlot(so.sub, features = marker.list[['Plasma']]) + coord_flip() + theme(axis.text.x = element_text(angle=90, hjust=1, size = 10)) 
ggsave(file = 'SubclusterPlasma_Dotplot.pdf', unit = 'in', width = 4, height = 4)


## Save Plasma data
saveRDS(so.sub, file = 'PanCoV_SeuratData_vstIntegrated_PlasmaSubcluster.RDS')
so.sub  <- readRDS(file = 'PanCoV_SeuratData_vstIntegrated_PlasmaSubcluster.RDS')


## Rename Plasma clusters in so.int
cells <- CellsByIdentities(so.sub)
for(celltype in names(cells)) Idents(so.int, cells = cells[[celltype]]) <- celltype 
so.int[['cell.type']] <- Idents(so.int)


## Save labeled data
saveRDS(so.int, file = 'PanCoV_SeuratData_vstIntegrated_v3.RDS')


