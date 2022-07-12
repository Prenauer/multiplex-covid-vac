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


## PW analysis
for(celltype in c('CD4Tcell','CD8Tcell','Bcell')){
  for(analysis in c('delta-vs-pbs','panHigh-vs-pbs','panLow-vs-pbs')){
    filename <- paste0('DE_edgeR_', celltype, '_',analysis,'_results.txt')
    if(file.exists(filename)){
      res <- read.delim(file = filename, sep = '\t', header = T)
      
      gostList <- res[(abs(res$logFC) > 0.5) & (res$FDR < .01),] %>% row.names()
      gostRes <- gprofiler2::gost(query = gostList, organism = "mmusculus",evcodes = T, ordered_query = T,significant = F,
                                  user_threshold = 1,correction_method = 'gSCS',domain_scope = 'known', sources = c('GO:BP'))
      gostRes$result$percent <-  gostRes$result$intersection_size/gostRes$result$term_size
      gScale <- structure(scale(res$logFC), names = rownames(res))
      gostRes$result$ActivationScore <- sapply(gostRes$result$intersection, function(x) {mean(gScale[str_split(x,',') %>% unlist()])})
      gostRes$result$parents <- sapply(gostRes$result$parents, function(x) {paste0(x, sep = '', collapse = ',')})
      write.table(gostRes$result,paste0('gse_gProfiler2_',analysis,'_',celltype,'_results.txt'), sep = '\t', col.names = T, row.names = F, quote = F)
      gostRes <- gostRes$result[(gostRes$result$term_size <= 500),]
      gostRes <- gostRes[abs(gostRes$ActivationScore) > .5,]
      gem <- data.frame(GO.ID = gostRes$term_id, Description = gostRes$term_name, p.Val = gostRes$p_value,FDR = gostRes$p_value, 
                        Phenotype = as.character(sign(gostRes$ActivationScore)), Genes = gostRes$intersection, Score = gostRes$ActivationScore, 
                        LogFDR = -log10(gostRes$p_value))
      write.table(gem,paste0('gse_gProfiler2_',analysis,'_',celltype,'_GEM.txt'), sep = '\t', col.names = T, row.names = F, quote = F)
    }
  }
}


## Network analysis
for(celltype in c('CD4Tcell','CD8Tcell','Bcell')){
  for(analysis in c('delta-vs-pbs','panHigh-vs-pbs','panLow-vs-pbs','panHigh-vs-delta','panLow-vs-delta')){
    gostRes <- read.delim(paste0('gse_gProfiler2_',analysis,'_',celltype,'_results.txt'), sep = '\t', header = T)
    gostRes <- gostRes[(gostRes$term_size <= 500) & (gostRes$intersection_size > 1)& (gostRes$p_value < 0.01),]
    gostRes <- gostRes[abs(gostRes$ActivationScore) > .5,]
    deRes <- read.delim(paste0('DE_edgeR_',celltype,'_',analysis,'_results.txt'), header = T, sep = '\t')
    v <- Pathway_Network_Plot(gostRes,clusters_to_label = 5, label.size = 4, filename = paste0('netplot2_',analysis,'_',celltype,'.pdf'))
    write.table(v,paste0('gse_gProfiler2_networkclustering_',analysis,'_',celltype,'.txt'), sep = '\t', col.names = T, row.names = F, quote = F)
    seeds <- v[which(nchar(v$label) > 0), 'name']
    MakeRidge(deRes = deRes, goRes = gostRes, goTerms = seeds, filename = paste0('ridgeplot2_gse_gProfiler_',analysis,'_',celltype,'.pdf'))
  }
}


## Get meta-PW data for all analyses
metaPW <- lapply(c('CD4Tcell','CD8Tcell','Bcell'), function(celltype){
  # for each analysis and each celltype, get the most signif GO term of the top 5 term-clusters (meta-PWs)
  data <- lapply(c('delta-vs-pbs','panHigh-vs-pbs','panLow-vs-pbs'), function(analysis){
    v <- read.delim(paste0('gse_gProfiler_networkclustering_',analysis,'_',celltype,'.txt'), sep = '\t', header = T)
    df <- data.frame()
    for(cluster in head(sort(unique(v$cluster)),5)) df <- rbind(df, v[which(v$cluster == cluster),] %>% .[with(., order(pval, decreasing = F)),c('name','description')] %>% head(.,1))
    return(df)
  })
  return(rbind(data[[1]],data[[2]],data[[3]]) %>% unique())
}) %>% structure(., names = c('CD4Tcell','CD8Tcell','Bcell'))


# Get all DE genes for each meta-PW, across all celltypes and DE-analyses
metaPWgenes <- sapply(c('CD4Tcell','CD8Tcell','Bcell'), function(celltype){
  sapply(metaPW[[celltype]]$name, function(go) {
    sapply(c('delta-vs-pbs','panHigh-vs-pbs','panLow-vs-pbs'), function(analysis){
      gostRes <- read.delim(paste0('gse_gProfiler_',analysis,'_',celltype,'_results.txt'), sep = '\t', header = T)
      return(gostRes[(gostRes$term_id == go),'intersection'] %>% str_split(., ',') %>% unlist())
    }) %>% unlist() %>% unique()
  })
})


# Make heatmaps of top 5 meta-PWs for each celltype
so.int  <- readRDS(file = 'PanCoV_SeuratData_vstIntegrated_v3.RDS')
DefaultAssay(so.int) <- 'RNA'
so.int <- DietSeurat(so.int, assays = 'RNA', counts = T)
cells <- list(CD4Tcell = so.int@meta.data[which(so.int@meta.data$cell.type %in% c('Th1','Th2','Treg')),] %>% rownames(),
              CD8Tcell = so.int@meta.data[which(so.int@meta.data$cell.type %in% c('CD8 TEM','CD8 TCM','CD8 T Effector')),] %>% rownames(),
              Bcell = so.int@meta.data[which(so.int@meta.data$cell.type %in% c('Activated B cell','Unswitched Memory B cell','Switched Memory B cell')),] %>% rownames())
for(celltype in c('CD4Tcell','CD8Tcell','Bcell')){
  so <- subset(so.int, cells = cells[[celltype]])
  DefaultAssay(so) <- 'RNA'
  so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 10000) 
  ## loop thru meta-PW, subset data by DE genes, scale, average expression by treatment, make heatmap
  for(i in 1:nrow(metaPW[[celltype]])){
    go_id <- metaPW[[celltype]]$name[i]
    go_desc <- metaPW[[celltype]]$description[i] 
    go_desc <- paste0(toupper(substr(go_desc,1,1)),substr(go_desc,2,nchar(go_desc)))
    glist <- metaPWgenes[[celltype]][[go_id]]
    df <- ScaleData(so, features = glist)
    df <- AverageExpression(df, assays='RNA',slot = 'scale.data',features = glist, group.by = 'treatment')$RNA
    fontsize_col <- 10
    if(nrow(df) > 22) fontsize_col <- floor(220/nrow(df))
    colnames(df) <- str_replace(colnames(df), 'panHigh', 'mixCoV-Hi') %>% str_replace(., 'panLow', 'mixCoV-Lo')  %>% str_replace(., 'delta', 'Delta')
    pheatmap::pheatmap(t(df), scale = 'none',clustering_method = 'ward.D2',fontsize_col = fontsize_col, treeheight_row = 10, treeheight_col = 10, cutree_cols = 2, 
                       cutree_rows = 2, main = go_desc, angle_col = '90', filename = paste0('heatmap2_',celltype,i,'_horiz.pdf'), height = 2,width = 4.5)
  }
}


##  Pathway bubble plots
for(celltype in c('CD4Tcell','CD8Tcell','Bcell')){
  deComp <- lapply(c('delta','panHigh','panLow'), function(treatment) {read.delim(paste0('DE_edgeR_',celltype,'_',treatment,'-vs-pbs_results.txt'), header = T, sep = '\t')})
  names(deComp) <- c('delta','panHigh','panLow')
  pwRes1 <- read.delim(paste0('gse_gProfiler_delta-vs-pbs_',celltype,'_results.txt'), header = T, sep = '\t')
  for(treatment in c('panHigh','panLow')){
    pwRes2 <- read.delim(paste0('gse_gProfiler_',treatment,'-vs-pbs_',celltype,'_results.txt'), header = T, sep = '\t')
    # merge df1 and pwRes2
    pwRes <- merge(pwRes1[,c('term_id','ActivationScore','p_value','term_name')],
                   pwRes2[,c('term_id','ActivationScore','p_value','term_name')], by = 'term_id')
    pwRes$LogFDR.x <- -log10(pwRes$p_value.x)
    pwRes$LogFDR.y <- -log10(pwRes$p_value.y)
    pwRes$ActivationRatio <- (pwRes$ActivationScore.y/pwRes$ActivationScore.x)
    pwRes <- pwRes[with(pwRes, order(LogFDR.y, decreasing = T)),]
    ## filter the pwRes.y for weak GO enrichment
    top <- list(all.up = which(pwRes$ActivationScore.y > (pwRes$ActivationScore.x + 1.2)), 
                all.dn = which(pwRes$ActivationScore.y < (pwRes$ActivationScore.x - 1.2)),
                top.up = pwRes[(pwRes$ActivationScore.y > (pwRes$ActivationScore.x + 1.2)) &  (pwRes$p_value.y < 0.01),] %>% head(.,4) %>% .[,'term_id'], 
                top.dn = pwRes[(pwRes$ActivationScore.y < (pwRes$ActivationScore.x - 1.2)) &  (pwRes$p_value.y < 0.01),] %>% head(.,4) %>% .[,'term_id']) 
    top[['label']] <- lapply(top[3:4], function(x) {
        sapply(x,function(go) {which(pwRes$term_id == go)})
    }) %>% unlist()%>% as.integer() %>% sort()
    pwRes$label <- sapply(pwRes$term_name.y, function(s) {paste0(toupper(substr(s,1,1)),substr(s,2,nchar(s)))}) %>% as.character()
    pwRes$label <- paste0(pwRes$label, ' (adj. p = ',format(pwRes$p_value.y,scientific = T,digits = 3),')')
    PathwayBubblePlot(data = pwRes,xlab = 'ActivationScore.x', ylab = 'ActivationScore.y', up = top$all.up, dn = top$all.dn, size = 'LogFDR.y', 
                      top.label = top$label, label = 'label',y_axis_label = paste0(str_replace(treatment, 'pan', 'mixCoV-'),' activ. score'),
                      title_label = paste0(str_replace(treatment, 'pan', 'mixCoV-'), 'Delta  vs PBS pathways'), 
                      filename = paste0('pwBubble2_gProfiler_',treatment,'-Delta_',celltype,'.pdf'))
  }
}


## Multi-ridge plots
for(celltype in c('CD4Tcell','CD8Tcell','Bcell')){
    deResList <- lapply(c('panHigh','panLow','delta'), function(treatment) {
      de <- read.delim(paste0('DE_edgeR_',celltype,'_',treatment,'-vs-pbs_results.txt'), header = T, sep = '\t')
      de$Gene <- rownames(de)
      return(de)
    })
    de <- merge(deResList[[1]][,c('Gene','logFC')],deResList[[2]][,c('Gene','logFC')], by = 'Gene')
    de <- merge(de,deResList[[3]][,c('Gene','logFC')], by = 'Gene')
    colnames(de) <- c('Gene','mixCoV-Hi','mixCoV-Lo','Delta')
    rownames(de) <- de$Gene
    de <- de[,-1]
    de <- reshape2::melt(as.matrix(de))
    colnames(de) <- c('Gene','Treatment','LogFC')
    ridgedata <- data.frame()
    for(n in names(metaPWgenes[[celltype]])) ridgedata <- rbind(ridgedata, data.frame(Term = n, de[which(de$Gene %in% metaPWgenes[[celltype]][[n]]),]))
    ridgedata$Treatment <- factor(ridgedata$Treatment, levels = c('Delta','mixCoV-Lo','mixCoV-Hi'))
    ridgedata$Term <- factor(ridgedata$Term, levels = rev(metaPW[[celltype]]$description))
  ggplot(ridgedata, aes(y = Term, x = LogFC, fill = Treatment)) +
    geom_density_ridges(alpha = 0.8) + 
    geom_vline(xintercept = 0, color = 'gray20', linetype = 2) +
    labs(x = 'Expr. log-FC', y = '') +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) + 
    theme_classic()
  ggsave(filename = paste0('multiridgeplot2_',celltype,'_pbs_v3.pdf'), dpi = 300, width = 5, height = 2.5, units = 'in')
}


