library(ggplot2)
library(ggrepel)
library(org.Hs.eg.db)







L24 <- read.table('../../Deseq2 Analysis/shift15/coding/csv_raw/ER_L24_t1.csv', sep = ',',header = T)
S15 <- read.table('../../Deseq2 Analysis/shift15/coding/csv_raw/ER_S15_t1.csv', sep = ',',header = T)
#TEfc <- read.table('csv_raw/TE_t2_t1.csv',sep = ',', header = T)
#RNAfc <- read.table('csv_raw/RNA_T2_against_T1.csv',sep = ',', header = T)

TIS <- read.table('../../Database/accnum2init_eff.csv',sep=',',header=T)
MFE_fpUTR <- read.table('../rnafold/utr_and_orf_minE_normbyLen.csv',sep=',',header=T)


accNums = L24$X

## filter unsure genes, and unused genes
TIS_notnan <- !is.nan(TIS$Efficiency)
TIS_id <- is.element(TIS$AccNum,accNums)
TIS_filtered <- TIS[TIS_id & TIS_notnan,]

MFE_filtered <- MFE_fpUTR[is.element(MFE_fpUTR$AccNum,TIS_filtered$AccNum),]

filteredGenetable <- read.table('../../Database/filteredGenesDetails_human_240118.txt',header=T)
ff <- is.element(filteredGenetable$AccNum,TIS_filtered$AccNum)
fftable <- filteredGenetable[ff,]


L24 <- L24[order(L24$log2FoldChange),]
L24 <- L24[L24$padj<0.05,]
L24up <- L24[L24$log2FoldChange>1,]
L24down <- head(L24,nrow(L24up))

## note: TIS_filtered and MFE_filtered have exact same accnums and order
L24_up_id = is.element(TIS_filtered$AccNum,L24up$X)
L24_down_id = is.element(TIS_filtered$AccNum,L24down$X)

S15 <- S15[order(S15$log2FoldChange),]
S15 <- S15[S15$padj<0.05,]
S15up <- S15[S15$log2FoldChange>1,]
S15down <- head(S15,nrow(S15up))

S15_up_id = is.element(TIS_filtered$AccNum,S15up$X)
S15_down_id = is.element(TIS_filtered$AccNum,S15down$X)

list_of_genes <- mapIds(org.Hs.eg.db,
                        keys = as.character(TIS_filtered$AccNum),
                        column = 'SYMBOL',
                        keytype ='ACCNUM',
                        multiVals = 'first')


group = rep('Others',length(list_of_genes))


group[L24_up_id] <- "L24 up"
group[L24_down_id] <- "L24 down"
group[S15_up_id] <- "S15 up"
group[S15_down_id] <- "S15 down"
group[L24_up_id & S15_up_id] <- "Both up"
group[L24_down_id & S15_down_id] <- "Both down"
group[L24_up_id & S15_down_id] <- "L24 up S15 down"
group[L24_down_id & S15_up_id] <- "L24 down S15 up"


df = data.frame(TIS = TIS_filtered$Efficiency,
                MFE = MFE_filtered$Energy,
                Symbol = list_of_genes,
                #X = TEfc$X,
                Group = group)
df$Group <- factor(df$Group, levels = c("L24 up","L24 down","S15 up","S15 down","Both up","Both down",'Others'))

ggplot(df, aes(x=MFE, y=TIS)) + 
  geom_point(data = df[df$Group == "Others",], alpha = 0.5, colour = "grey") + 
  geom_point(data = df[df$Group != "Others",], aes(colour = Group)) + 
  #geom_abline(slope=1, linetype = 'dashed', color = 'brown', size = 0.7) +
  #scale_fill_discrete(values = c("#d8b365", "#f5f5f5", "#5ab4ac")) +
  scale_color_manual(values=c('red','green','orange','blue','cyan','purple')) +
  #ylim(-1.5,2) + 
  #xlim(-4,8) +
  theme_bw() +

  ggtitle("TIS agst MFE_fpUTR_orf100_normByLen") +
  theme(legend.position = c(0.1,0.1),
        #legend.position = 'none',
        legend.title = element_blank(),
        legend.background = element_rect(color = 'black',size = 0.5,linetype= 'solid'),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = rel(1.5),hjust = 0.5),
        axis.title = element_text(size = rel(1.25)),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size=14))
ggsave('TIS agst MFE_fpUTR_orf100_norm_by_len.png',width = 9, height = 7)




make_plot <- function() {
  ggplot(df, aes(x=MFE, y=TIS)) + 
    geom_point(data = df[df$Group == "Others",], alpha = 0.5, colour = "grey") + 
    geom_point(data = df[df$Group != "Others",], aes(colour = Group)) + 
    #geom_abline(slope=1, linetype = 'dashed', color = 'brown', size = 0.7) +
    #scale_fill_discrete(values = c("#d8b365", "#f5f5f5", "#5ab4ac")) +
    scale_color_manual(values=c('red','green','orange','blue','cyan','purple')) +
    #ylim(-1.5,2) + 
    #xlim(-4,8) +
    theme_bw() +
    
    ggtitle("TIS agst MFE_fpUTR_orf100_normByLen") +
    theme(legend.position = c(0.2,0.92),
          #legend.position = 'none',
          legend.title = element_blank(),
          legend.background = element_rect(color = 'black',size = 0.5,linetype= 'solid'),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = rel(1.5),hjust = 0.5),
          axis.title = element_text(size = rel(1.25)),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size=14))
}


keyw = 'mitoComplexesAll'
genes_table <- read.table('../../Database/Gene Groups/mitoComplexAll.csv',sep=',',header=T)
View(genes_table)
genes_table = as.data.frame(sapply(genes_table, toupper))

group = rep('Others',length(list_of_genes))
igrp = is.element(fftable$AccNum,genes_table$X)
group[igrp] <- keyw
df = data.frame(TIS = TIS_filtered$Efficiency,
                MFE = MFE_filtered$Energy,
                Symbol = list_of_genes,
                #X = TEfc$X,
                Group = group)
df$Group <- factor(df$Group, levels = c(keyw,'Others'))
make_plot()
ggsave(paste('TIS agst MFE_fpUTR_orf100_norm_by_len_',keyw,'.png',sep=''),width = 9, height = 7)

