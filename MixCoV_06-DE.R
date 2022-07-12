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


## Load data
so.int <- readRDS(file = 'PanCoV_SeuratData_vstIntegrated_v3.RDS')


## Get cellnames for following activated celltypes: 'CD4Tcell','CD8Tcell', and 'Bcell'
cells <- list(CD4Tcell = so.int@meta.data[which(so.int@meta.data$cell.type %in% c('Th1','Th2','Treg')),] %>% rownames(),
              CD8Tcell = so.int@meta.data[which(so.int@meta.data$cell.type %in% c('CD8 TEM','CD8 TCM','CD8 T Effector')),] %>% rownames(),
              Bcell = so.int@meta.data[which(so.int@meta.data$cell.type %in% c('Activated B cell','Unswitched Memory B cell','Switched Memory B cell')),] %>% rownames(),
              ActivatedBcell = so.int@meta.data[which(so.int@meta.data$cell.type == 'Activated B cell'),] %>% rownames(),
              UnswitchedMemoryBcell = so.int@meta.data[which(so.int@meta.data$cell.type == 'Unswitched Memory B cell'),] %>% rownames())


## DE analyses using EdgeR workflow, using cell detection rate and treatment as covariates
DefaultAssay(so.int) <- 'RNA'
so.int <- DietSeurat(so.int, assays = 'RNA', counts = T)
for(celltype in c('CD4Tcell','CD8Tcell','Bcell')){
  so <- subset(so.int, cells = cells[[celltype]])
  for(comp.cell1 in c('panHigh','panLow','delta')){
    for(comp.cell2 in c('delta','pbs')){
      if(comp.cell1 != comp.cell2){
        so.sub <- subset(so, cells = so@meta.data[which(so[['treatment']] ==comp.cell1 | so[['treatment']] == comp.cell2),] %>% rownames())
        data <- so.sub$RNA@counts
        gene.detection.rate <- apply(data,1, function(x) {sum(x > 0)/length(x)})
        data <- data[(gene.detection.rate > 0.05),]
        dge <- DGEList(data, group = factor(so.sub@meta.data$treatment))
        dge <- calcNormFactors(dge, method = 'TMM')
        treatment <- as.integer(so.sub@meta.data$treatment == comp.cell1)
        cell.detection.rate <- so.sub@meta.data$cell.detection.rate
        design <- model.matrix(~ cell.detection.rate + treatment)
        dge <- estimateDisp(dge, design = design)
        fit <- glmQLFit(dge, design = design)
        qlf <- glmQLFTest(fit, coef = 'treatment')
        tt <- topTags(qlf, n = Inf)$table
        filename <- paste0('DE_edgeR_',celltype,'_',comp.cell1,'-vs-',comp.cell2,'_results')
        write.table(tt, paste0(filename,'.txt'), sep = '\t', col.name = T, row.name = T, quote = F)
        save(dge,fit,qlf,tt, file = paste0(filename,'.Rdata'))
      }
    }
  }
}


##  Log-FC square plots
for(celltype in c('CD4Tcell','CD8Tcell','Bcell')){
  deComp <- lapply(c('delta','panHigh','panLow'), function(treatment) {
      read.delim(paste0('DE_edgeR_',celltype,'_',treatment,'-vs-pbs_results.txt'), header = T, sep = '\t')
  })
  names(deComp) <- c('delta','panHigh','panLow')
  SquarePlot(df = merge(deComp$panLow, deComp$panHigh, by = 0), lab.x = 'mixCoV-Lo', lab.y = 'mixCoV-Hi',
             filename = paste0('square2_mixCoVHi-vs-mixCoVLo_',celltype,'_v3.pdf'))
  SquarePlot(df = merge(deComp$delta, deComp$panLow, by = 0), lab.x = 'Delta', lab.y = 'mixCoV-Lo',
             filename = paste0('square2_mixCoVLo-vs-delta_',celltype,'_v3.pdf'))
  SquarePlot(df = merge(deComp$delta, deComp$panHigh, by = 0), lab.x = 'Delta', lab.y = 'mixCoV-Hi',
             filename = paste0('square2_mixCoVHi-vs-delta_',celltype,'_v3.pdf'))
}


