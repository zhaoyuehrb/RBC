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
  
  res_tableDE <- data.frame(res_tableDE)
  #res_tableDE$SYMBOL <- mapIds(org.Hs.eg.db,
  #                             keys = as.character(rownames(res_tableDE)),
  #                             column = 'SYMBOL',
  #                             keytype ='ACCNUM',
  #                             multiVals = 'first')
  threshold_DE <- (res_tableDE$padj < 0.05) & (res_tableDE$log2FoldChange >1)
  
  res_tableDE$threshold <- threshold_DE
  #res_tableDE$threshold[res_tableDE$SYMBOL %in% list_of_genes] <- 'highlight'
  
  #res_tableDE$genelabels <- ""
  #res_tableDE$genelabels[res_tableDE$SYMBOL %in% list_of_genes] <- label
  #View(res_tableDE)
  ggplot(res_tableDE,aes(x = log2FoldChange, y = -log10(padj))) +
    
    #ylim(ylims[[1]],ylims[[2]]) +
    ylim(0,10) +
    xlim(-2.5,2.5) +
    
    #threshold lines
    geom_vline(xintercept = -1, linetype = 'dashed', color = 'brown',size = 0.7) +
    geom_vline(xintercept = 1, linetype = 'dashed', color = 'brown', size = 0.7) +
    geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'brown', size = 0.7) +
    
    #points
    geom_point(data = res_tableDE[res_tableDE$threshold,] ,colour = 'red',alpha = 0.5,aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(data = res_tableDE[res_tableDE$threshold==FALSE,] ,colour = 'black',alpha = 0.5,aes(x = log2FoldChange, y = -log10(padj))) +
    #geom_point(alpha = 0.5) +
    scale_fill_manual(name = '',
                      values = c(label='blue'))+
    
    guides(colour = F)+
    
    ###add gene labels for particular group from list_of_genes
    #geom_text_repel(data = subset(res_tableDE, ((log2FoldChange > 1 | log2FoldChange < -1) & (FALSE))),
                    #data = subset(res_tableDE, SYMBOL %in% c('Htra2','Tomm7','Sqstm1')), #for rna autophagy
                    #data = subset(res_tableDE, SYMBOL %in% list_of_genes),
    #                aes(label = SYMBOL), size = 3,
    #                box.padding = unit(0.4, 'lines'),
    #                point.padding = unit(0.4, 'lines'),
    #                segment.size = 0.2, segment.colour = 'grey50')+
    
    theme_bw() + 
    
    ggtitle(paste(type,sep = '')) +
    xlab(bquote(~log[2]~ "FC")) +
    ylab(bquote(~-log[10]~italic(p-adj)))+
    
    theme(legend.position = c(0.15,0.9),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          #legend.position = 'none',
          legend.title = element_blank(),
          #legend.background = element_rect(color = 'black',size = 0.5,linetype= 'solid'),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = rel(1.5),hjust = 0.5),
          axis.title = element_text(size = rel(1.25)),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size=14))
  
  #filepath <- paste(save_dir,folder_label,'/',type,toString(day),'.png',sep = '')
  filepath <- paste('paperplots','/',type,'.png',sep = '')
  ggsave(filepath,width = 9, height = 7)
}


plot_volcano(res_0_L24,'ER L24 T0',1,c(),'','',xlims=c(-2.5,2.5), ylims=c(0,10))
plot_volcano(res_1_L24,'ER L24 T1',1,c(),'','',xlims=c(-2.5,2.5), ylims=c(0,10))
plot_volcano(res_2_L24,'ER L24 T2',1,c(),'','',xlims=c(-2.5,2.5), ylims=c(0,10))

plot_volcano(res_0_S15,'ER S15 T0',1,c(),'','',xlims=c(-2.5,2.5), ylims=c(0,10))
plot_volcano(res_1_S15,'ER S15 T1',1,c(),'','',xlims=c(-2.5,2.5), ylims=c(0,10))
plot_volcano(res_2_S15,'ER S15 T2',1,c(),'','',xlims=c(-2.5,2.5), ylims=c(0,10))

plot_volcano(res_0_S24,'ER S24 T0',1,c(),'','',xlims=c(-2.5,2.5), ylims=c(0,10))
plot_volcano(res_1_S24,'ER S24 T1',1,c(),'','',xlims=c(-2.5,2.5), ylims=c(0,10))
plot_volcano(res_2_S24,'ER S24 T2',1,c(),'','',xlims=c(-2.5,2.5), ylims=c(0,10))


plot_volcano(L24_0_S15,'L24 against S15 T0',1,c(),'','',xlims=c(-5,5),ylims=c(0,20))
plot_volcano(L24_1_S15,'L24 against S15 T1',1,c(),'','',xlims=c(-5,5),ylims=c(0,20))
plot_volcano(L24_2_S15,'L24 against S15 T2',1,c(),'','',xlims=c(-5,5),ylims=c(0,20))


plot_volcano(res_1_against_0_RNA,'RNA T1 against T0',1,c(),'','',xlims=c(-5,5), ylims=c(0,120))
plot_volcano(res_2_against_0_RNA,'RNA T2 against T0',1,c(),'','',xlims=c(-5,5), ylims=c(0,250))


plot_volcano(res_1_against_0_RPF,'RPF T1 against T0',1,c(),'','',xlims=c(-15,15), ylims=c(0,1200))
plot_volcano(res_2_against_0_RPF,'RPF T2 against T0',1,c(),'','',xlims=c(-15,15), ylims=c(0,1200))


plot_volcano(res_1_te_0,'TE T1 against T0',1,c(),'','',xlims=c(-5,5), ylims=c(0,12))
plot_volcano(res_2_te_0,'TE T2 against T0',1,c(),'','',xlims=c(-5,5), ylims=c(0,15))

