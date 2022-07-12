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
so.int <- readRDS(file = 'PanCoV_SeuratData_vstIntegrated_v3.RDS')


## UMAP plot with updated celltype labels
ggcolor <- function(n) hcl(h = seq(15, 375, length = n + 1), l = 65, c = 100)[1:n]
dimplot.cols <- ggcolor(23)[sample(1:23, 22, replace = F)]
dimplot.cols <- c('#71B000','#FF6C92','#00C0B7','#CA9700','#00BB4B','#BE80FF','#F265E7','#00AEFA','#97A900','#F8766D','#00BF76', '#00BDD1','#2FB600','#8F91FF','#DD8D00','#FF64B3','#00C098','#B3A000','#EC823C','#DE71F9','#00B7E8','#00B0F6')
DimPlot(so, cols = dimplot.cols, label = T, label.box = T, repel = T, raster = T) + theme(legend.position = "none")
ggsave(file = 'PanCoV_SeuratData_vstIntegrated_UMAP_AllCelltypes.pdf', unit = 'in', width = 6, height = 4.5)


## UMAP plots, labeled by condition
DimPlot(so.int, label = T, group.by = 'treatment', label.box = T, repel = T, raster = T, pt.size = .01) + theme(legend.position = "none")
ggsave(file = 'PanCoV_SeuratData_vstIntegrated_UMAP_byTreatment.pdf', unit = 'in', width = 6, height = 4.5)
dimplot.cols <- c('#EC69EF', '#A58AFF', '#00C0B5', '#00C094', '#53B400', '#00BDD2', '#DA8F00', '#00B6EB', '#EB8335', '#FF6B96', '#F8766D', '#FF63B9', '#00BA38', '#00ABFD', '#A9A400', '#C49A00', '#D078FF', '#86AC00', '#00BE6D', '#619CFF', '#FB61D7')
DimPlot(so.int, cols = dimplot.cols, label = F, label.box = F, repel = F, raster = T, pt.size = .01) + theme(legend.position = "none")
ggsave(file = 'PanCoV_SeuratData_vstIntegrated_UMAP_AllCelltypes.pdf', unit = 'in', width = 6, height = 4.5)


## Heatmap of cluster-specific expression patterns
markers <- RunPrestoAll(so.int, assay = 'integrated',test.use = 'wilcox',slot = 'scale.data', max.cells.per.ident = 5000, random.seed = 1101)
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_diff)
p <- DoHeatmap(so.int, group.by = 'cell.type', slot = 'scale.data', features = unique(top10$gene))
ggsave(plot = p, filename = paste0('heatmap_byCelltype.pdf'), dpi = 300, width = 24, height = 6, units = 'in')
write.table(markers, paste0('allMarkers_wilcox.txt'), sep = '\t', quote = T, col.names = T, row.names = F)


## Featureplots of Celltype markers
celltypes <- list(CD4Tcell = c('Treg','Naive CD4 T','Th1','Th2'), CD8Tcell = c('Naive CD8 T','CD8 TCM', 'CD8TEM','CD8 T Effector'), Myeloid = c('NK cell','Macrophage','pDC','Monocyte','cDC1','cDC2'), Bcell = c('Activated B cell','Naive B cell', 'Unswitched B cell','Switched B cell'), Plasma = c('Plasmablast','Plasma cell','Pre-plasmablast'))
marker.list <- list(Tcell = c('Cd3d','Cd8b1','Cd4','Ccr7','Cd44','Tcf7'), 
  CD4_Tcell1 = c('Cd3d','Cd4','Tbx21','Gata3','Foxp3'), CD4_Tcell2 = c('Il2ra','Gzmb','Ifng','Tgfb1','Il4'), 
  Plasma1 = c('Cd19','Ms4a1','Cr2','Cd9','Sdc1','H2-Aa'), Plasma2 = c('H2-Eb2','Prdm1','Pax5','Mki67','Irf4','Bach2'), 
  Bcell1 = c('Cd19','Ms4a1','Cr2','Cd27','Spn','H2-Eb2'),  Bcell2 = c('H2-Aa','Ighm','Igha','Ighd','Cd24a','Fcer2a'), 
  Myeloid = c('Cd14','Fcgr3', 'Itgam', 'Itgax', 'Ncr1'), Mono_Macro_DC_subsets1 = c('Itgam','Itgax','Cd24a','Bst2','Cd8b1','Csf1r'), 
  Mono_Macro_DC_subsets2 = c('Cd14','Cd68','Adgre1','Ly6c1','Fcgr3','Fcgr1'), Mono_Macro_DC_subsets3 = c('Itgae','Sirpa', 'Siglech', 'Ccr2','Cx3cr1'))
for(i in names(marker.list)){
  FeaturePlot(so.int, pt.size = .05, raster = T, min.cutoff = 'q01', max.cutoff = 'q66', features = marker.list[[i]])
  ggsave(file = paste0('FeaturePlot_UMAP_',i,'.pdf'), unit = 'in', width = 5.5, height = 5)
}


## Barplot of celltype proportions between conditions
df <- so.int@meta.data[,c('sample','treatment','cell.type')]
df <- ddply(df, c('sample','treatment','cell.type'), summarise, cells = length(sample))
df.total <- ddply(df, 'sample', summarise, total = sum(cells))
df <- merge(df,df.total, by = 'sample')
df$percent <- 100*df$cells/df$total
cell.type.order <- sapply(unique(df$cell.type), function(x) {mean(df[(df$treatment == 'pbs') & (df$cell.type == x),'percent'])}) %>% structure(., names = as.character(unique(df$cell.type))) %>% sort(., decreasing = T) %>% names()
df$cell.type <- factor(df$cell.type, levels = cell.type.order)
df$treatment <- factor(df$treatment, levels = c('pbs','delta','panLow','panHigh'))
p <- ggplot(df, aes(fill=treatment, y= percent, x= cell.type)) + geom_boxplot(outlier.shape = NA, outlier.alpha = 0.1) +  
geom_point(position=position_dodge(width = 0.75), size = 0.2) + theme_classic() +  theme(axis.text.x = element_text(angle=90, hjust=1)) 
ggsave(plot = p, filename = 'Immune_Celltype_Comparisons.pdf', dpi = 300, width = 6, height = 3, units = 'in')
df <- so.int@meta.data[,c('sample','treatment','cell.type')]
df <- ddply(df, c('treatment','cell.type'), summarise, cells = length(cell.type))
df.total <- ddply(df, 'treatment', summarise, total = sum(cells))
df <- merge(df,df.total, by = 'treatment')
df$percent <- 100*df$cells/df$total
df$cell.type <- factor(df$cell.type, levels =(cell.type.order))
df$treatment <- factor(df$treatment, levels = rev(c('pbs','delta','panLow','panHigh')))
p <- ggplot(df, aes(fill=cell.type, y=percent, x=treatment)) + geom_bar(position = 'stack', stat = 'identity') +  coord_flip() + theme_classic() 
ggsave(plot = p, filename = 'Immune_Celltype_Comparisons_stackedBarplot.pdf', dpi = 300, width = 6, height = 3, units = 'in')


## Dotplots for celltype markers
celltype.list <- list(Bcell = c('Naive B cell','Activated B cell','Switched Memory B cell','Unswitched Memory B cell'),
  Plasma = c('Pre-plasmablast','Plasmablast','Plasma cell'), CD4Tcell = c('Naive CD4 T','Th1','Th2','Treg'), CD8Tcell = c('Naive CD8 T','CD8 T Effector','CD8 TEM','CD8 TCM'), 
  DC = c('cDC1','cDC2','pDC'), OtherMyeloid = c('Monocyte','Macrophage','NK cell'))
so.int[['cell.type']] <- factor(as.character(so.int[['cell.type']]), levels = unlist(celltype.list) %>% unique())
marker.list <- list(CD8Tcell = c('Cd3d','Cd8b1','Ccr7','Cd44','Tcf7','Gzmb','Ifng'), 
  CD4Tcell = c('Cd3d','Cd4','Tbx21','Gata3','Foxp3','Il2ra','Gzmb','Ifng','Tgfb1','Il4'), 
  Plasma = c('Cd19','Ms4a1','Cr2','Cd9','Sdc1','H2-Aa','H2-Eb2','Prdm1','Pax5','Mki67','Irf4','Bach2'), 
  Bcell = c('Fcer2a','Cd83','Cd69','H2-Aa','H2-Eb2','Ighd','Ighm','Ighg3','Spib','Cd38','Cd40','Sdc1'), 
  DC = c('Itgam','Itgax','Cd24a','Bst2','Cd8b1','Sirpa', 'Siglech','Itgae'),   
  OtherMyeloid = c('Itgam','Csf1r','Cd14','Cd68','Adgre1','Ly6c1','Fcgr3', 'Ccr2','Cx3cr1','Ncr1'))
DotPlot(so.int, features = marker.list[['OtherMyeloid']], idents = celltype.list[['OtherMyeloid']], cols = c('gray80','#009933'), cluster.idents = F) + theme(axis.text.x = element_text(angle=90, hjust=1, size = 10)) 
ggsave(filename = 'Feature_dotplot_myeloid.pdf', dpi = 300, width = 5.3, height = 4, units = 'in')
DotPlot(so.int, features = marker.list[['DC']], idents = celltype.list[['DC']], cols = c('gray80','#619CFF'), cluster.idents = F) + theme(axis.text.x = element_text(angle=90, hjust=1, size = 10)) 
ggsave(filename = 'Feature_dotplot_DC.pdf', dpi = 300, width = 4.5, height = 4, units = 'in')
DotPlot(so.int, features = marker.list[['CD4Tcell']], idents = celltype.list[['CD4Tcell']], cols = c('gray80','#F8766D'), cluster.idents = F) + theme(axis.text.x = element_text(angle=90, hjust=1, size = 10)) 
ggsave(filename = 'Feature_dotplot_CD4Tcell.pdf', dpi = 300, width = 5.2, height = 4, units = 'in')
DotPlot(so.int, features = marker.list[['CD8Tcell']], idents = celltype.list[['CD8Tcell']] %>% rev(), cols = c('gray80','#C49A00'), cluster.idents = F) + theme(axis.text.x = element_text(angle=90, hjust=1, size = 10)) 
ggsave(filename = 'Feature_dotplot_CD8Tcell.pdf', dpi = 300, width = 5.2, height = 4, units = 'in')
DotPlot(so.int, features = marker.list[['Bcell']], idents = celltype.list[['Bcell']], cols = c('gray80','#EB8335'), cluster.idents = F)  + theme(axis.text.x = element_text(angle=90, hjust=1, size = 10)) 
ggsave(filename = 'Feature_dotplot_Bcell.pdf', dpi = 300, width = 7, height = 4, units = 'in')
DotPlot(so.int, features = marker.list[['Plasma']], idents = celltype.list[['Plasma']], cols = c('gray80','#A58AFF'), cluster.idents = F) + theme(axis.text.x = element_text(angle=90, hjust=1, size = 10)) 
ggsave(filename = 'Feature_dotplot_Plasma.pdf', dpi = 300, width = 6, height = 4, units = 'in')


## Export metadata (nMolecules, nFeatures, sample, treatment, cell.type; UMAP embedding)
df <- so.int@meta.data[,c('seurat_clusters','lane','sample','treatment','nCount_RNA','nFeature_RNA','percent.mt','cell.type')]
write.table(df, paste0('so-int_metadata.txt'), sep = '\t', quote = T, col.names = T, row.names = T)
d <- so.int@meta.data[,c('sample','treatment','cell.type')]
df <- so.int@reductions$umap@cell.embeddings
df <- merge(df,d,by = 0)
write.table(df, paste0('so-int_umap-embeddings.txt'), sep = '\t', quote = F, col.names = T, row.names = T)

