###Volcano Plot
library(ggplot2)
library(ggrepel)
library(org.Hs.eg.db)
library(AnnotationDbi)

###res_tableDE: object from running results in DESeq
###type: rnabulk/rnacyt/..etc 
###list_of_genes: particular group of genes to be highlighted, can leave it as c() if none highlighted
###label: label of particular group of genes to be highlighted, can leave it as '' if none highlighted
###folder_label: name of folder to save plots
###xlims: c(x_min, x_max); ylims: c(y_min, y_max)

plot_volcano <- function(res_tableDE, type, day,list_of_genes,label,folder_label,xlims,ylims){
  
  
  #  list_of_genes <- mapIds(org.Hs.eg.db,
  #                          keys = as.character(genes_table$X),
  #                          column = 'SYMBOL',
  #                          keytype ='ACCNUM',
  #                          multiVals = 'first')
  label <- 'Original'
  
  res_tableDE <- data.frame(res_tableDE)
  res_tableDE$SYMBOL <- mapIds(org.Hs.eg.db,
                               keys = as.character(rownames(res_tableDE)),
                               column = 'SYMBOL',
                               keytype ='ACCNUM',
                               multiVals = 'first')
  threshold_DE <- res_tableDE$padj < 0.05
  
  res_tableDE$threshold <- threshold_DE
  #res_tableDE$threshold[res_tableDE$SYMBOL %in% list_of_genes] <- 'highlight'
  
  res_tableDE$genelabels <- ""
  #res_tableDE$genelabels[res_tableDE$SYMBOL %in% list_of_genes] <- label
  #res_tableDE$genelabels[row.names(res_tableDE) %in% list_of_genes] <- label
  #View(res_tableDE)
  
  ### criteria for labelling
  #sub_dat <- subset(res_tableDE, SYMBOL %in% list_of_genes)
  right <- unname(quantile(res_tableDE$log2FoldChange,0.999))
  left <- unname(quantile(res_tableDE$log2FoldChange,0.001))
  up <- unname(quantile(-log10(res_tableDE$padj),0.998))
  
  ggplot(res_tableDE,aes(x = log2FoldChange, y = -log10(padj))) +
    
    ylim(ylims[[1]],ylims[[2]]) +
    xlim(xlims[[1]], xlims[[2]]) +
    
    #threshold lines
    geom_vline(xintercept = -1, linetype = 'dashed', color = 'brown',size = 0.7) +
    geom_vline(xintercept = 1, linetype = 'dashed', color = 'brown', size = 0.7) +
    geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'brown', size = 0.7) +
    
    #points
    geom_point(aes(colour = threshold),alpha = 0.5) +
    #geom_point(aes(colour = "#808080",alpha = 0.5),colour = "#A9A9A9",alpha = 0.6) +
    #geom_point(data = subset(res_tableDE, SYMBOL %in% list_of_genes),
               #geom_point(data = subset(res_tableDE, row.names(res_tableDE) %in% list_of_genes),          
    #           aes(fill = genelabels),colour = 'blue')+
    # geom_text_repel(aes(x = log2FoldChange, y = -log10(padj), 
    #                     label = ifelse(genelabels==TRUE,res_tableDE$SYMBOL,""),size = 3),
    #                 box.padding = unit(0.2,'lines'),
    #                 point.padding = unit(0.2, 'lines'))+
    
    scale_fill_manual(name = '',
                      values = c(label='blue'))+
    
    guides(colour = F)+
    
    ###add gene labels for particular group from list_of_genes
    geom_text_repel(#data = subset(res_tableDE, ((log2FoldChange > 1 | log2FoldChange < -1) & (FALSE))),
      #data = subset(res_tableDE, SYMBOL %in% c('Htra2','Tomm7','Sqstm1')), #for rna autophagy
      data = subset(res_tableDE,  (SYMBOL %in% c('AHSP','VEGFA','ACTB','ALAS2','SLC4A1','TMEM175') | (log2FoldChange>right) | (log2FoldChange<left) | (-log10(padj)>up))),
      #data = subset(res_tableDE, (log2FoldChange > 1.5)),
      aes(label = SYMBOL), size = 3,
      box.padding = unit(0.4, 'lines'),
      point.padding = unit(0.4, 'lines'),
      segment.size = 0.2, segment.colour = 'grey50')+
    
    theme_bw() + 
    
    ggtitle(paste(type,' ',label,sep = '')) +
    xlab(bquote(~log[2]~ "FC")) +
    ylab(bquote(~-log[10]~italic(p-adj)))+
    
    theme(legend.position = c(0.15,0.9),
          #legend.position = 'none',
          legend.title = element_blank(),
          legend.background = element_rect(color = 'black',size = 0.5,linetype= 'solid'),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = rel(1.5),hjust = 0.5),
          axis.title = element_text(size = rel(1.25)),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size=14))
  
  save_dir <- paste('plots/',label,'/',sep = '')
  
  
  #save_dir <- '/Users/Yue/Documents/HGLab/deseq_yue/Proper Version/volc_doublecheck/'
  #filepath <- paste(save_dir,folder_label,'/',type,toString(day),'.png',sep = '')
  filepath <- paste(save_dir,'/',type,' ',label,'.png',sep = '')
  ggsave(filepath,width = 9, height = 7)
  
}


plot_volcano(res_0_L24,'ER L24 T0',1,c(),'','',xlims=c(-5,5), ylims=c(0,10))

#plot_volcano(res_0_L24,'L24_0_label',1,c(subset(res_0_L24,log2FoldChange>1)),'','',xlims=c(-5,5), ylims=c(0,10))

plot_volcano(res_1_L24,'ER L24 T1',1,c(),'','',xlims=c(-5,5), ylims=c(0,10))
plot_volcano(res_2_L24,'ER L24 T2',1,c(),'','',xlims=c(-5,5), ylims=c(0,10))

plot_volcano(res_0_S15,'ER S15 T0',1,c(),'','',xlims=c(-5,5), ylims=c(0,26))
plot_volcano(res_1_S15,'ER S15 T1',1,c(),'','',xlims=c(-5,5), ylims=c(0,40))
plot_volcano(res_2_S15,'ER S15 T2',1,c(),'','',xlims=c(-5,5), ylims=c(0,20))

plot_volcano(res_0_S24,'ER S24 T0',1,c(),'','',xlims=c(-5,5), ylims=c(0,10))
plot_volcano(res_1_S24,'ER S24 T1',1,c(),'','',xlims=c(-5,5), ylims=c(0,10))
plot_volcano(res_2_S24,'ER S24 T2',1,c(),'','',xlims=c(-5,5), ylims=c(0,10))


plot_volcano(L24_0_S15,'L24 against S15 T0',1,c(),'','',xlims=c(-5,5),ylims=c(0,35))
plot_volcano(L24_1_S15,'L24 against S15 T1',1,c(),'','',xlims=c(-5,5),ylims=c(0,35))
plot_volcano(L24_2_S15,'L24 against S15 T2',1,c(),'','',xlims=c(-5,5),ylims=c(0,35))


plot_volcano(res_1_against_0_RNA,'RNA T1 against T0',1,c(),'','',xlims=c(-6,6), ylims=c(0,120))
plot_volcano(res_2_against_0_RNA,'RNA T2 against T0',1,c(),'','',xlims=c(-11,11), ylims=c(0,350))


plot_volcano(res_1_against_0_RPF,'RPF T1 against T0',1,c(),'','',xlims=c(-7,7), ylims=c(0,150))
plot_volcano(res_2_against_0_RPF,'RPF T2 against T0',1,c(),'','',xlims=c(-11,11), ylims=c(0,300))


plot_volcano(res_1_te_0,'TE T1 against T0',1,c(),'','',xlims=c(-5,5), ylims=c(0,12))
plot_volcano(res_2_te_0,'TE T2 against T0',1,c(),'','',xlims=c(-5,5), ylims=c(0,15))

