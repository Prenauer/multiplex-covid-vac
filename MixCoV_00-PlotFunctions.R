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


ggcolor <- function(n) hcl(h = seq(15, 375, length = n + 1), l = 65, c = 100)[1:n]

Pathway_Network_Plot <- function(gores, network_layout = layout_with_fr, clusters_to_label = 6, filename, width = 4, height = 6,label.size = 2.5, thresh.logfc = 0.5, thresh.fdr = 0.01){
  require(igraph)
  require(network)
  require(sna)
  require(ggplot2)
  require(ggrepel)
  require(stringr)
  
  ## Calculate similarity matrix
  SimMatrix <- function(res, overlap.pct = 0.5) # must have intersection and term_id columns
  { 
    require(stringr)
    # iterate rows, convert gene string to lists of character-vectors
    glists <- lapply(res$intersection, function(x) x %>% str_split(., ',') %>% unlist())
    gInter <- lapply(glists, function(g1) lapply(glists,function(g2) intersect(g1,g2) %>% length()) %>% unlist())
    gUnion <- lapply(glists, function(g1) lapply(glists,function(g2) union(g1,g2) %>% length()) %>% unlist())
    gMin <- lapply(glists, function(g1) lapply(glists,function(g2) min(length(g1),length(g2))) %>% unlist())
    m <- sapply(1:length(glists), function(i) {
      jaccard <- (gInter[[i]]/gUnion[[i]]) * (1 - overlap.pct)
      overlap <- (gInter[[i]]/gMin[[i]]) * overlap.pct
      return(jaccard + overlap)
    }) %>% as.matrix()
    dimnames(m) <- list(res$term_id,res$term_id)
    return(m)
  }
  ## Make edge matrix
  EdgeMatrix <- function(res, simThresh = 0.375, parallels = F, overlap.pct = 0.5){
    
    m <- SimMatrix(res, overlap.pct)
    ## filter by similarity score
    edges <- reshape2::melt(m, as.is=T)
    edges <- edges[(edges$Var1 != edges$Var2) & (edges$value >= simThresh),]
    if(!parallels){
      parallels <- apply(edges,1,function(x) sort(as.character(x[1:2])) %>% paste0(.,collapse = ','))
      edges <- edges[!duplicated(parallels),]
    }
    names(edges) <- c('from','to','similarity_coefficient')
    edges$name <- edges$from
    return(edges)
  }
  
  
  ## make iGraph
  e <- EdgeMatrix(gores)
  v <- gores[,c('term_id','term_name','term_size','percent','intersection','ActivationScore','p_value')]
  names(v) <- c('name','description','go_term_size','go_percent','genes','score','pval')
  net <- igraph::graph_from_data_frame(d = e,directed = F, vertices = v)
  
  ## cluster the nodes
  net <- set_edge_attr(net, 'weight', value =  edge_attr(net, 'similarity_coefficient'))
  l <- cluster_leiden(graph = net, objective_function = 'modularity', resolution_parameter = .2,n_iterations = 500)
  v$cluster <- membership(l)[v$name]
  
  ## choosing meta-pathways
  metaPW <- sapply(head(sort(unique(v$cluster)),clusters_to_label), function(cluster){
    goterms <- v[which(v$cluster == cluster),'name']
    df <- gores[which(gores$term_id %in% goterms),]
    return(df[with(df, order(p_value, decreasing = F)),'term_id'][1])
  })
  v$label <- ''
  v[which(v$name %in% metaPW),'label'] <- v[which(v$name %in% metaPW),'description']
  v$logFDR <- -log10(v$pval)
  
  ## prepare graphing attributes
  
  l <- network_layout(net) 
  if((max(l[,1]) - min(l[,1])) > (max(l[,2]) - min(l[,2]))) l <- l[,c(2,1)]
  v <- cbind(v,as.data.frame(l))
  v$cluster <- factor(v$cluster)
  v$label <- sapply(1:nrow(v), function(x) if(nchar(v$label[x]) > 0) paste0(toupper(substr(v$label[x],1,1)),substr(v$label[x],2,nchar(v$label[x])),' (adj.p = ',format(v$pval[x],scientific = T,digits = 3),')') else '')
  e <- igraph::as_data_frame(net, 'edges')
  e.names <- colnames(e)
  e <- merge(e, v[,c('name','V1','V2')], by.x = 'from', by.y = 'name')
  e <- merge(e, v[,c('name','V1','V2')], by.x = 'to', by.y = 'name')
  colnames(e) <- c(e.names, 'x0','y0','x1','y1')
  e$size <- e$similarity_coefficient *1
  v$size <- (v$logFDR/max(v$logFDR)) * 14
  
  ## plot the network                                                            
  p <- ggplot(e) +
    geom_segment(aes(x = x0, xend = x1, y = y0, yend = y1),size = e$size, color = 'gray80') +
    geom_point(data = v, aes(x = V1, y = V2, fill = cluster, size = size),colour = 'gray40',alpha = 0.8, shape = 21) +
    geom_label_repel(data = v, aes(x = V1, y = V2, label = str_wrap(label,40), fill = cluster, alpha = 0.3), size = label.size, min.segment.length = 1, fontface = 'bold', color = 'black', max.overlaps = 1000) + 
    theme_void() +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
          plot.background = element_rect(colour = "grey50")) +
    scale_fill_discrete(guide = 'none') + 
    scale_alpha(guide = 'none') + 
    scale_size_continuous(name = str_wrap('Adj. p (-log10)',8))
  
  ggsave(plot = p, filename, dpi = 300, height = 4, width = 4,units = 'in')
  return(v)
}

VolcanoPlot <- function(df, filename){
  #df <- df[with(df, order(logFC, decreasing = F)),]
  top <- list(up = df[(df$logFC > 0.5) & (df$LogFDR > 2),],
              dn = df[(df$logFC < -0.5) & (df$LogFDR > 2),])
  toptop <- list(up = head(df[(df$logFC > .8) & (df$LogFDR > 3),],10),
                 dn = head(df[(df$logFC < -.8) & (df$LogFDR > 3),],10))
  p <- ggplot(df, aes(x=logFC, y=LogFDR)) +
    geom_point(size = 1.5, color = 'gray10', fill = 'gray 50') +
    geom_segment(x = -0.5, xend = min(df$logFC), y = 2, yend = 2, color = 'gray50', linetype = 2) +
    geom_segment(x = 0.5, xend = max(df$logFC), y = 2, yend = 2, color = 'gray50', linetype = 2) +
    geom_segment(x = -0.5, xend = -0.5, y = 2, yend = max(df$LogFDR), color = 'gray50', linetype = 2) +
    geom_segment(x = 0.5, xend = 0.5, y = 2, yend = max(df$LogFDR), color = 'gray50', linetype = 2) +
    geom_point(data = top$up, aes(x=logFC, y=LogFDR), color = 'firebrick') +
    geom_point(data = top$dn, aes(x=logFC, y=LogFDR), color = 'steelblue') +
    geom_label_repel(data = toptop$up, aes(x=logFC, y=LogFDR, label = symbol_gene), alpha = .7, color = 'black', max.overlaps = 20) + 
    geom_label_repel(data = toptop$dn, aes(x=logFC, y=LogFDR, label = symbol_gene), alpha = .7, color = 'black', max.overlaps = 20) + 
    labs(x = 'Log fold-change', y = 'q value (-log10)') +
    theme_classic()
  ggsave(plot = p, filename = filename, dpi = 300, width = 4, height = 3, units = 'in')
}
                    
MakeRidge <- function(deRes, goRes, goTerms,thresh.logfc = 0.5, thresh.fdr = 0.01, width = 4, height = 3, filename = NULL){
  require(ggridges)
  j = lapply(goTerms, function(x) str_split(goRes[str_detect(x,goRes$term_id),'intersection'], pattern = ',') %>% unlist())
  tmp <-deRes[(abs(deRes$logFC) > thresh.logfc) & (deRes$FDR < thresh.fdr),]
  #tmp <- sapply(goTerms, function(x) tmp[which(row.names(tmp) %in% bp[[x]]$genes), 'logFC'])
  tmp <- lapply(j, function(x) tmp[x,'logFC'])
  tmp <- tmp[lapply(tmp,length) %>% unlist() > 2]
  goTerms <- goTerms[1:length(tmp)]
  goTermDesc <- sapply(goTerms, function(x) goRes[str_detect(x,goRes$term_id),'term_name'] %>% unique())
  goTermDesc <- sapply(goTermDesc, function(s) paste0(toupper(substr(s,1,1)),substr(s,2,nchar(s))))
  ridgedata <- data.frame()
  for(i in 1:length(tmp)) ridgedata <- suppressWarnings(rbind(ridgedata, data.frame(Term = goTermDesc[[i]], LogFC = tmp[[i]])))
  
  ridgedata$Term <- factor(ridgedata$Term, levels = rev(goTermDesc))
  col <- ggcolor(length(goTerms))
  ggplot((ridgedata), aes(y = (Term), x = (LogFC), fill = (Term))) +
    geom_density_ridges() + 
    scale_fill_manual(values = rev(col),guide = 'none') + 
    labs(x = 'Expr. log-FC', y = '') +
    geom_vline(xintercept = 0, linetype = 'dashed', color = 'gray50') +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) + 
    theme_classic()
  
  if(!is.null(filename))  ggsave(filename = filename, dpi = 300, width = width, height = height, units = 'in')
}

enrichResult <- function(result, gene = NULL, geneSets = NULL,gost = T, pvalueCutoff = 0.05, qvalueCutoff = 0.05, pAdjustMethod = 'BH'){
  if(gost){
    geneSets <- gprofiler2::gconvert(query = names(gene),organism = "mmusculus", target = "GO")
    geneSets <- geneSets[(geneSets$input == geneSets$name),]
    geneSets <- daply(geneSets, 'target', summarise, genes = paste0(name, sep = '', collapse = ','))
    geneSets <- lapply(geneSets, function(x) str_split(x, ','))
    result$ActivationScore <- sapply(gostRes$intersection, function(x) mean(gene[str_split(x,',') %>% unlist()]))
    result <- data.frame(ID = result$term_id, Description = result$term_name,
                         setSize = result$term_size, enrichmentScore = result$ActivationScore,
                         pvalue = result$p_value,
                         p.adjust = result$p_value, qvalues = result$p_value, 
                         core_enrichment = result$intersection)
  }
  if(is.null(geneSets)) geneSets <- list()
  x <- new("enrichResult", result = result, pvalueCutoff = pvalueCutoff, 
           pAdjustMethod = pAdjustMethod, qvalueCutoff = qvalueCutoff, 
           gene = as.character(gene), universe = 'UNKNOWN', geneSets = geneSets, 
           organism = "UNKNOWN", keytype = "UNKNOWN", ontology = "UNKNOWN", 
           readable = FALSE)
  return(x)
}
                                     
                                     
MultiRidge <- function(deResList, goResList, goTerms, filename){
  # get list of pooled term genes that were hits in the mixCov-vs-pbs analyses.
  for(i in 1:3) rownames(goResList[[i]]) <- goResList[[i]]$name
  goTermGenes <- lapply(goTerms, function(x) lapply(1:3, function(i) 
    if(sum(goResList[[i]]$name %in% x) > 0) 
      str_split(goResList[[i]][which(goResList[[i]]$name %in% x),'genes'],',') %>% unlist()) %>% unlist() %>% unique())
   ridgedataList <- lapply(1:length(deResList), function(i) {
    lapply(goTermGenes, function(x) {
      gList <- intersect(x,rownames(deResList[[i]]))
      return(deResList[[i]][gList,'logFC'])
    })
  })
  names(ridgedataList) <- names(deResList)
  goTermDesc <- sapply(goTerms, function(x) {
    df <- rbind(rbind(goResList[[1]],goResList[[2]]), goResList[[3]])
    df[str_detect(x,df$name),'description'] %>% unique()
  })
  goTermDesc <- sapply(goTermDesc, function(s) paste0(toupper(substr(s,1,1)),substr(s,2,nchar(s))))
  
  ridgedata <- data.frame()
  for(treatment in names(ridgedataList)){
    tmp <- data.frame()
    for(i in 1:length(goTerms)) tmp <- suppressWarnings(rbind(tmp, data.frame(Term = goTermDesc[i],Treatment = treatment, LogFC = ridgedataList[[treatment]][[i]])))
    ridgedata <- rbind(ridgedata, tmp)
  }
  ridgedata$Term <- factor(ridgedata$Term, levels = rev(goTermDesc))
  ridgedata$Treatment <- str_replace(ridgedata$Treatment, 'panHigh', 'mixCoV-Hi') %>% str_replace(., 'panLow', 'mixCoV-Lo')  %>% str_replace(., 'delta', 'Delta')
  ridgedata$Treatment <- factor(ridgedata$Treatment, levels = c('Delta','mixCoV-Lo','mixCoV-Hi'))
  
  ggplot(ridgedata, aes(y = Term, x = LogFC, fill = Treatment)) +
    geom_density_ridges(alpha = 0.8) + 
    geom_vline(xintercept = 0, color = 'gray20', linetype = 2) +
    labs(x = 'Expr. log-FC', y = '') +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) + 
    theme_classic()
  ggsave(filename = filename, dpi = 300, width = 5, height = 2.5, units = 'in')
}

SquarePlot <- function(df, logFC.x='logFC.x', logFC.y='logFC.y',FDR.y='FDR.y', filename, lab.x, lab.y){
  df[,'logFC.x'] <- df[,logFC.x]
  df[,'logFC.y'] <- df[,logFC.y]
  df[,'FDR.y'] <- df[,FDR.y]
  df$LogFDR <- -log10(df$FDR.y)
  df$LogFDR[which(df$LogFDR == Inf)] <- max(df$LogFDR[which(df$LogFDR != Inf)])
  df <- df[with(df, order(logFC.y, decreasing = T)),]
  top <- list(all.up = df[(df$logFC.y > 0.5),], 
              all.dn = df[(df$logFC.y < -0.5),],
              top.up = df[(df$logFC.y > 0.5) & (df$FDR.y < 0.01),] %>% head(.,8),
              top.dn = df[(df$logFC.y < -0.5) & (df$FDR.y < 0.01),] %>% tail(.,8))
  ggplot(df,aes(x = logFC.x, y = logFC.y, size = LogFDR)) +
    geom_point(color = 'gray40',  alpha = 0.2) +
    geom_point(data = top$all.up, color = 'firebrick', alpha = 0.6) +
    geom_point(data = top$all.dn, color = 'steelblue', alpha = 0.6) +
    scale_size(range = c(.2, 4), name= paste0(lab.y,' q val.')) +
    geom_vline(xintercept = c(-0.5,0.5), color = 'gray50', linetype = 2) +
    geom_hline(yintercept = c(-0.5,0.5), color = 'gray50', linetype = 2) +
    geom_label_repel(data = top$top.up, aes(x = logFC.x, y = logFC.y, label = Row.names), size = 3, fontface = 'italic', alpha = .7, color = 'black', max.overlaps = 20) + 
    geom_label_repel(data = top$top.dn, aes(x = logFC.x, y = logFC.y, label = Row.names), size = 3, fontface = 'italic', alpha = .7, color = 'black', max.overlaps = 20) +
    labs(x = paste0(lab.x,' log-FC'), y = paste0(lab.y,' log-FC'),title = paste0(lab.y,' DEG comparison')) +
    theme_classic()  
  ggsave(filename, dpi = 300, units = 'in', height = 3,width = 5)
}

PathwayBubblePlot <- function(data,xlab,ylab,up,dn,size,top.label = NULL, label,filename, x_axis_label = 'Delta activ. score',y_axis_label = 'mixCoV activ. score',title_label = '', height = 3, width = 4){
  data[,'xlab'] <- data[,xlab]
  data[,'ylab'] <- data[,ylab]
  data[,'size'] <- data[,size] + 0.5
  data[,'label'] <- data[,label]
  p <- ggplot(data, aes(x = xlab, y = ylab, size = size)) +
    geom_point(color = 'gray40',  alpha = 0.2) +
    scale_size(range = c(.2, 10), name= str_wrap('Adj. p (-log10)',8)) +
    geom_abline(slope = c(1,1), intercept = c(1.2,-1.2), color = 'gray50', linetype = 2) +
    geom_point(data = data[up,], color = 'firebrick', alpha = 0.6) +
    geom_point(data = data[dn,], color = 'steelblue', alpha = 0.6)
  if(!is.null(top.label)) p <- p + 
    geom_label_repel(data = data[top.label,],aes(x = xlab, y = ylab, label = str_wrap(label,25)),color = 'black',force = 1, size = 2, fontface = 'bold', max.overlaps = 800, alpha = 0.6)
    p <- p + labs(x = x_axis_label, y = y_axis_label,title = title_label) + theme_classic()  
  ggsave(plot = p, filename, dpi = 300, units = 'in', height = height,width = width)
}

