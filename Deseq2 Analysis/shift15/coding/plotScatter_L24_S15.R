library(ggplot2)
library(ggrepel)
library(org.Hs.eg.db)


L24 <- read.table('csv_raw/ER_L24_t1.csv', sep = ',',header = T)
S15 <- read.table('csv_raw/ER_S15_t1.csv', sep = ',',header = T)
TEfc <- read.table('csv_raw/TE_t2_t1.csv',sep = ',', header = T)
RNAfc <- read.table('csv_raw/RNA_T2_against_T1.csv',sep = ',', header = T)


L24 <- L24[order(L24$log2FoldChange),]
L24 <- L24[L24$padj<0.05,]
L24up <- L24[L24$log2FoldChange>1,]
L24down <- head(L24,nrow(L24up))
list_of_genes <- mapIds(org.Hs.eg.db,
                        keys = as.character(TEfc$X),
                        column = 'SYMBOL',
                        keytype ='ACCNUM',
                        multiVals = 'first')
up_id = is.element(TEfc$X,L24up$X)
down_id = is.element(TEfc$X,L24down$X)
group = rep('Others',7990)
group[up_id] <- "L24 up"
group[down_id] <- "L24 down"

df = data.frame(TEfc = TEfc$log2FoldChange,
                RNAfc = RNAfc$log2FoldChange,
                Symbol = list_of_genes,
                X = TEfc$X,
                Group = group)
df$Group <- factor(df$Group, levels = c("L24 up","L24 down","Others"))
ggplot(df, aes(x=RNAfc, y=TEfc)) + 
  geom_point(data = df[df$Group == "Others",], alpha = 0.5, colour = "grey") + 
  geom_point(data = df[df$Group != "Others",], aes(colour = Group)) + 
  geom_abline(slope=1, linetype = 'dashed', color = 'brown', size = 0.7) +
  #scale_fill_discrete(values = c("#d8b365", "#f5f5f5", "#5ab4ac")) +
  scale_color_manual(values=c('red','blue')) +
  #                    guide=TRUE) +
  
  ylim(-2,2) + 
  #ylim(min(df$RNAfc,df$TEfc),max(df$RNAfc,df$TEfc)) +
  xlim(min(df$RNAfc,df$TEfc),max(df$RNAfc,df$TEfc)) +
  theme_bw() +

   geom_text_repel(
     data = subset(df, (X %in% L24up$X | X %in% L24down$X) & (abs(RNAfc) > 1.8 | abs(TEfc > 1.4))),
     aes(label = Symbol), size = 3,
     box.padding = unit(0.4, 'lines'),
     point.padding = unit(0.4, 'lines'),
     force = 3,
     segment.size = 0.2, segment.colour = 'grey50'
   ) + 

  ggtitle("TE fold change against RNA fold change") +
  theme(legend.position = c(0.15,0.9),
        #legend.position = 'none',
        legend.title = element_blank(),
        legend.background = element_rect(color = 'black',size = 0.5,linetype= 'solid'),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = rel(1.5),hjust = 0.5),
        axis.title = element_text(size = rel(1.25)),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size=14))



#save_dir = 'compareBoth/'
#filename = paste(save_dir,file.label,'.png',sep = '')
ggsave('L24S15/L24 up down in TE and RNA.png',width = 9, height = 7)



S15 <- S15[order(S15$log2FoldChange),]
S15 <- S15[S15$padj<0.05,]
S15up <- S15[S15$log2FoldChange>1,]
S15down <- head(S15,nrow(S15up))

up_id = is.element(TEfc$X,S15up$X)
down_id = is.element(TEfc$X,S15down$X)
group = rep('Others',7990)
group[up_id] <- "S15 up"
group[down_id] <- "S15 down"

df = data.frame(TEfc = TEfc$log2FoldChange,
                RNAfc = RNAfc$log2FoldChange,
                Symbol = list_of_genes,
                X = TEfc$X,
                Group = group)
df$Group <- factor(df$Group, levels = c("S15 up","S15 down","Others"))
ggplot(df, aes(x=RNAfc, y=TEfc)) + 
  geom_point(data = df[df$Group == "Others",], alpha = 0.5, colour = "grey") + 
  geom_point(data = df[df$Group != "Others",], aes(colour = Group)) + 
  geom_abline(slope=1, linetype = 'dashed', color = 'brown', size = 0.7) +
  #scale_fill_discrete(values = c("#d8b365", "#f5f5f5", "#5ab4ac")) +
  scale_color_manual(values=c('red','blue')) +
  #                    guide=TRUE) +
  
  ylim(-2,2) +
  #ylim(min(df$RNAfc,df$TEfc),max(df$RNAfc,df$TEfc)) +
  xlim(min(df$RNAfc,df$TEfc),max(df$RNAfc,df$TEfc)) +
  theme_bw() +
  
  
  geom_text_repel(
    data = subset(df, (X %in% S15up$X | X %in% S15down$X) & (abs(RNAfc) > 2.3 | RNAfc < -1.85 | TEfc > 1.35 | TEfc < -0.9)),
    aes(label = Symbol), size = 3,
    box.padding = unit(0.4, 'lines'),
    point.padding = unit(0.4, 'lines'),
    force = 3,
    segment.size = 0.2, segment.colour = 'grey50'
  ) + 
ggtitle("TE fold change against RNA fold change") +
  theme(legend.position = c(0.15,0.9),
        #legend.position = 'none',
        legend.title = element_blank(),
        legend.background = element_rect(color = 'black',size = 0.5,linetype= 'solid'),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = rel(1.5),hjust = 0.5),
        axis.title = element_text(size = rel(1.25)),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size=14))
ggsave('L24S15/S15 up down in TE and RNA.png',width = 9, height = 7)










both_up <- is.element(TEfc$X,S15up$X) & is.element(TEfc$X,L24up$X)
group = rep('Others',7990)
group[is.element(TEfc$X,S15up$X)] <- "S15 up"
group[is.element(TEfc$X,L24up$X)] <- "L24 up"
group[both_up] <- "Both up"
df = data.frame(TEfc = TEfc$log2FoldChange,
                RNAfc = RNAfc$log2FoldChange,
                Symbol = list_of_genes,
                X = TEfc$X,
                Group = group)
df$Group <- factor(df$Group, levels = c("L24 up","S15 up","Both up","Others"))

ggplot(df, aes(x=RNAfc, y=TEfc)) + 
  geom_point(data = df[df$Group == "Others",], alpha = 0.5, colour = "grey") + 
  geom_point(data = df[df$Group != "Others",], aes(colour = Group)) + 
  geom_abline(slope=1, linetype = 'dashed', color = 'brown', size = 0.7) +
  #scale_fill_discrete(values = c("#d8b365", "#f5f5f5", "#5ab4ac")) +
  scale_color_manual(values=c('red','green','orange')) +
  #                    guide=TRUE) +
  
  ylim(-2,2) + 
  #ylim(min(df$RNAfc,df$TEfc),max(df$RNAfc,df$TEfc)) +
  xlim(min(df$RNAfc,df$TEfc),max(df$RNAfc,df$TEfc)) +
  theme_bw() +
  geom_text_repel(
    data = subset(df, (X %in% L24up$X | X %in% S15up$X) & (RNAfc > 2.3 | RNAfc < -1.7 | abs(TEfc > 1.6))),
    aes(label = Symbol), size = 3,
    box.padding = unit(0.4, 'lines'),
    point.padding = unit(0.4, 'lines'),
    force = 3,
    segment.size = 0.2, segment.colour = 'grey50'
  ) + 
  ggtitle("TE fold change against RNA fold change") +
    theme(legend.position = c(0.15,0.9),
          #legend.position = 'none',
          legend.title = element_blank(),
          legend.background = element_rect(color = 'black',size = 0.5,linetype= 'solid'),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = rel(1.5),hjust = 0.5),
          axis.title = element_text(size = rel(1.25)),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size=14))
ggsave('L24S15/L24 and S15 up in TE and RNA.png',width = 9, height = 7)

