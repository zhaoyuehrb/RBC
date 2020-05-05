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

L24_up_id = is.element(TEfc$X,L24up$X)
L24_down_id = is.element(TEfc$X,L24down$X)

S15 <- S15[order(S15$log2FoldChange),]
S15 <- S15[S15$padj<0.05,]
S15up <- S15[S15$log2FoldChange>1,]
S15down <- head(S15,nrow(S15up))

S15_up_id = is.element(TEfc$X,S15up$X)
S15_down_id = is.element(TEfc$X,S15down$X)


#L24 up and down
#anchor
group = rep('Others',7990)
group[L24_up_id] <- "L24 up"
group[L24_down_id] <- "L24 down"

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
  ylim(-2,2) + 
  xlim(-8,8) +
  theme_bw() +
  geom_text_repel(
    data = subset(df, (X %in% L24up$X | X %in% L24down$X) & (RNAfc > 2.3 | RNAfc < -1.3 | abs(TEfc > 1.4))),
    aes(label = Symbol), size = 5,
    box.padding = unit(0.4, 'lines'),
    point.padding = unit(0.4, 'lines'),
    force = 3,
    segment.size = 0.2, segment.colour = 'grey50'
  ) + 
  ggtitle("TE fold change against RNA fold change") +
  theme(#legend.position = c(0.8,0.9),
        #legend.position = 'none',
        #legend.title = element_blank(),
        #legend.background = element_rect(color = 'black',size = 0.5,linetype= 'solid'),
        #legend.text = element_text(size = 12),
        plot.title = element_text(size = rel(1.5),hjust = 0.5),
        axis.title = element_text(size = rel(2.25)),
        axis.text.y = element_text(size = 22),
        axis.text.x = element_text(size=24))
ggsave('scatter_L24_S15_feb/L24 up down t1 in TE and RNA t2_t1.png',width = 9, height = 7)

#S15 up and down
#anchor
group = rep('Others',7990)
group[S15_up_id] <- "S15 up"
group[S15_down_id] <- "S15 down"

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
  ylim(-2,2) + 
  xlim(-8,8) +
  theme_bw() +
  geom_text_repel(
    data = subset(df, (X %in% S15up$X | X %in% S15down$X) & (RNAfc > 2.5 | RNAfc < -1.5 | abs(TEfc > 1.6))),
    aes(label = Symbol), size = 5,
    box.padding = unit(0.4, 'lines'),
    point.padding = unit(0.4, 'lines'),
    force = 3,
    segment.size = 0.2, segment.colour = 'grey50'
  ) + 
  ggtitle("TE fold change against RNA fold change") +
  theme(#legend.position = c(0.8,0.9),
    #legend.position = 'none',
    #legend.title = element_blank(),
    #legend.background = element_rect(color = 'black',size = 0.5,linetype= 'solid'),
    #legend.text = element_text(size = 12),
    plot.title = element_text(size = rel(1.5),hjust = 0.5),
    axis.title = element_text(size = rel(2.25)),
    axis.text.y = element_text(size = 22),
    axis.text.x = element_text(size=24))
ggsave('scatter_L24_S15_feb/S15 up down t1 in TE and RNA t2_t1.png',width = 9, height = 7)














L24 <- read.table('csv_raw/ER_L24_t0.csv', sep = ',',header = T)
S15 <- read.table('csv_raw/ER_S15_t0.csv', sep = ',',header = T)



L24 <- L24[order(L24$log2FoldChange),]
L24 <- L24[L24$padj<0.05,]
L24up <- L24[L24$log2FoldChange>1,]
L24down <- head(L24,nrow(L24up))
list_of_genes <- mapIds(org.Hs.eg.db,
                        keys = as.character(TEfc$X),
                        column = 'SYMBOL',
                        keytype ='ACCNUM',
                        multiVals = 'first')

L24_up_id = is.element(TEfc$X,L24up$X)
L24_down_id = is.element(TEfc$X,L24down$X)

S15 <- S15[order(S15$log2FoldChange),]
S15 <- S15[S15$padj<0.05,]
S15up <- S15[S15$log2FoldChange>1,]
S15down <- head(S15,nrow(S15up))

S15_up_id = is.element(TEfc$X,S15up$X)
S15_down_id = is.element(TEfc$X,S15down$X)


#L24 up and down
#anchor
group = rep('Others',7990)
group[L24_up_id] <- "L24 up"
group[L24_down_id] <- "L24 down"

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
  ylim(-2,2) + 
  xlim(-8,8) +
  theme_bw() +
  geom_text_repel(
    data = subset(df, (X %in% L24up$X | X %in% L24down$X) & TRUE),
    aes(label = Symbol), size = 5,
    box.padding = unit(0.4, 'lines'),
    point.padding = unit(0.4, 'lines'),
    force = 3,
    segment.size = 0.2, segment.colour = 'grey50'
  ) + 
  ggtitle("TE fold change against RNA fold change") +
  theme(#legend.position = c(0.8,0.9),
    #legend.position = 'none',
    #legend.title = element_blank(),
    #legend.background = element_rect(color = 'black',size = 0.5,linetype= 'solid'),
    #legend.text = element_text(size = 12),
    plot.title = element_text(size = rel(1.5),hjust = 0.5),
    axis.title = element_text(size = rel(2.25)),
    axis.text.y = element_text(size = 22),
    axis.text.x = element_text(size=24))
ggsave('scatter_L24_S15_feb/L24 up down t1 in TE and RNA t2_t1.png',width = 9, height = 7)

#S15 up and down
#anchor
group = rep('Others',7990)
group[S15_up_id] <- "S15 up"
group[S15_down_id] <- "S15 down"

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
  ylim(-2,2) + 
  xlim(-8,8) +
  theme_bw() +
  geom_text_repel(
    data = subset(df, (X %in% S15up$X | X %in% S15down$X) & (RNAfc > 3 | RNAfc < -1.8 | TEfc > 1.5 | TEfc < -0.9)),
    aes(label = Symbol), size = 5,
    box.padding = unit(0.4, 'lines'),
    point.padding = unit(0.4, 'lines'),
    force = 3,
    segment.size = 0.2, segment.colour = 'grey50'
  ) + 
  ggtitle("TE fold change against RNA fold change") +
  theme(#legend.position = c(0.8,0.9),
    #legend.position = 'none',
    #legend.title = element_blank(),
    #legend.background = element_rect(color = 'black',size = 0.5,linetype= 'solid'),
    #legend.text = element_text(size = 12),
    plot.title = element_text(size = rel(1.5),hjust = 0.5),
    axis.title = element_text(size = rel(2.25)),
    axis.text.y = element_text(size = 22),
    axis.text.x = element_text(size=24))
ggsave('scatter_L24_S15_feb/S15 up down t1 in TE and RNA t2_t1.png',width = 9, height = 7)


































#L24 up and S15 up
#anchor
group = rep('Others',7990)
group[L24_up_id] <- "L24 up"
group[S15_up_id] <- "S15 up"

both_up <- L24_up_id & S15_up_id
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
  ylim(-1.5,2) + 
  xlim(-4,8) +
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
    theme(legend.position = c(0.8,0.9),
          #legend.position = 'none',
          legend.title = element_blank(),
          legend.background = element_rect(color = 'black',size = 0.5,linetype= 'solid'),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = rel(1.5),hjust = 0.5),
          axis.title = element_text(size = rel(1.25)),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size=14))
ggsave('scatter_L24_S15_feb/L24 and S15 up t1 in TE and RNA t2_t1.png',width = 9, height = 7)

#L24 down and S15 down
#anchor
group = rep('Others',7990)
group[L24_down_id] <- "L24 down"
group[S15_down_id] <- "S15 down"

both_down <- L24_down_id & S15_down_id
group[both_down] <- "Both down"

df = data.frame(TEfc = TEfc$log2FoldChange,
                RNAfc = RNAfc$log2FoldChange,
                Symbol = list_of_genes,
                X = TEfc$X,
                Group = group)
df$Group <- factor(df$Group, levels = c("L24 down","S15 down","Both down","Others"))

ggplot(df, aes(x=RNAfc, y=TEfc)) + 
  geom_point(data = df[df$Group == "Others",], alpha = 0.5, colour = "grey") + 
  geom_point(data = df[df$Group != "Others",], aes(colour = Group)) + 
  geom_abline(slope=1, linetype = 'dashed', color = 'brown', size = 0.7) +
  #scale_fill_discrete(values = c("#d8b365", "#f5f5f5", "#5ab4ac")) +
  scale_color_manual(values=c('red','green','orange')) +
  ylim(-1.5,2) + 
  xlim(-4,8) +
  theme_bw() +
  geom_text_repel(
    data = subset(df, (X %in% L24down$X | X %in% S15down$X) & (RNAfc > 2.3 | RNAfc < -1.7 | TEfc < -0.8 | TEfc > -0.1)),
    aes(label = Symbol), size = 3,
    box.padding = unit(0.4, 'lines'),
    point.padding = unit(0.4, 'lines'),
    force = 3,
    segment.size = 0.2, segment.colour = 'grey50'
  ) + 
  ggtitle("TE fold change against RNA fold change") +
  theme(legend.position = c(0.8,0.9),
        #legend.position = 'none',
        legend.title = element_blank(),
        legend.background = element_rect(color = 'black',size = 0.5,linetype= 'solid'),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = rel(1.5),hjust = 0.5),
        axis.title = element_text(size = rel(1.25)),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size=14))
ggsave('scatter_L24_S15_feb/L24 and S15 down t1 in TE and RNA t2_t1.png',width = 9, height = 7)




#L24 up and S15 down
#anchor
group = rep('Others',7990)
group[L24_up_id] <- "L24 up"
group[S15_down_id] <- "S15 down"

both <- L24_up_id & S15_down_id
group[both] <- "Both"

df = data.frame(TEfc = TEfc$log2FoldChange,
                RNAfc = RNAfc$log2FoldChange,
                Symbol = list_of_genes,
                X = TEfc$X,
                Group = group)
df$Group <- factor(df$Group, levels = c("L24 up","S15 down","Both","Others"))

ggplot(df, aes(x=RNAfc, y=TEfc)) + 
  geom_point(data = df[df$Group == "Others",], alpha = 0.5, colour = "grey") + 
  geom_point(data = df[df$Group != "Others",], aes(colour = Group)) + 
  geom_abline(slope=1, linetype = 'dashed', color = 'brown', size = 0.7) +
  #scale_fill_discrete(values = c("#d8b365", "#f5f5f5", "#5ab4ac")) +
  scale_color_manual(values=c('red','green','orange')) +
  ylim(-1.5,2) + 
  xlim(-4,8) +
  theme_bw() +
  geom_text_repel(
    data = subset(df, (X %in% L24up$X | X %in% S15down$X) & (RNAfc > 2.3 | RNAfc < -1.7 | TEfc < -0.9 | TEfc > 1)),
    aes(label = Symbol), size = 3,
    box.padding = unit(0.4, 'lines'),
    point.padding = unit(0.4, 'lines'),
    force = 3,
    segment.size = 0.2, segment.colour = 'grey50'
  ) + 
  ggtitle("TE fold change against RNA fold change") +
  theme(legend.position = c(0.8,0.9),
        #legend.position = 'none',
        legend.title = element_blank(),
        legend.background = element_rect(color = 'black',size = 0.5,linetype= 'solid'),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = rel(1.5),hjust = 0.5),
        axis.title = element_text(size = rel(1.25)),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size=14))
ggsave('scatter_L24_S15_feb/L24 up and S15 down t1 in TE and RNA t2_t1.png',width = 9, height = 7)



#L24 down and S15 up
#anchor
group = rep('Others',7990)
group[L24_down_id] <- "L24 down"
group[S15_up_id] <- "S15 up"

both <- L24_down_id & S15_up_id
group[both] <- "Both"

df = data.frame(TEfc = TEfc$log2FoldChange,
                RNAfc = RNAfc$log2FoldChange,
                Symbol = list_of_genes,
                X = TEfc$X,
                Group = group)
df$Group <- factor(df$Group, levels = c("L24 down","S15 up","Both","Others"))

ggplot(df, aes(x=RNAfc, y=TEfc)) + 
  geom_point(data = df[df$Group == "Others",], alpha = 0.5, colour = "grey") + 
  geom_point(data = df[df$Group != "Others",], aes(colour = Group)) + 
  geom_abline(slope=1, linetype = 'dashed', color = 'brown', size = 0.7) +
  #scale_fill_discrete(values = c("#d8b365", "#f5f5f5", "#5ab4ac")) +
  scale_color_manual(values=c('red','green','orange')) +
  ylim(-1.5,2) + 
  xlim(-4,8) +
  theme_bw() +
  geom_text_repel(
    data = subset(df, (X %in% L24down$X | X %in% S15up$X) & (RNAfc > 3 | RNAfc < -2 | TEfc < -0.7 | TEfc > 1.2)),
    aes(label = Symbol), size = 3,
    box.padding = unit(0.4, 'lines'),
    point.padding = unit(0.4, 'lines'),
    force = 3,
    segment.size = 0.2, segment.colour = 'grey50'
  ) + 
  ggtitle("TE fold change against RNA fold change") +
  theme(legend.position = c(0.8,0.9),
        #legend.position = 'none',
        legend.title = element_blank(),
        legend.background = element_rect(color = 'black',size = 0.5,linetype= 'solid'),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = rel(1.5),hjust = 0.5),
        axis.title = element_text(size = rel(1.25)),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size=14))
ggsave('scatter_L24_S15_feb/L24 down and S15 up t1 in TE and RNA t2_t1.png',width = 9, height = 7)




## anchor
L24_col = rep('noChange',7990)
L24_col[L24_down_id] <- "down"
L24_col[L24_up_id] <- "up"

S15_col = rep('noChange',7990)
S15_col[S15_down_id] <- "down"
S15_col[S15_up_id] <- "up"


df = data.frame(Gene = list_of_genes,
                TEfc = TEfc$log2FoldChange,
                RNAfc = RNAfc$log2FoldChange,
                L24 = L24_col,
                S15 = S15_col)

write.csv(df,"scatter_L24_S15_feb/tefc_rnafc_t2t1_full.csv", row.names = FALSE)

redundent = df$L24 == 'noChange' & df$S15 == 'noChange'
df_short = df[c(!redundent),]
write.csv(df,"scatter_L24_S15_feb/tefc_rnafc_t2t1_short.csv", row.names = FALSE)





# anchor

#L24 up and S15 up
#anchor
group = rep('Others',7990)

both_up <- L24_up_id & S15_up_id
group[both_up] <- "L24 up S15 up"

both_down <- L24_down_id & S15_down_id
group[both_down] <- "L24 down S15 down"

df = data.frame(TEfc = TEfc$log2FoldChange,
                RNAfc = RNAfc$log2FoldChange,
                Symbol = list_of_genes,
                X = TEfc$X,
                Group = group)
df$Group <- factor(df$Group, levels = c("L24 up S15 up","L24 down S15 down","Others"))

ggplot(df, aes(x=RNAfc, y=TEfc)) + 
  geom_point(data = df[df$Group == "Others",], alpha = 0.5, colour = "grey") + 
  geom_point(data = df[df$Group != "Others",], aes(colour = Group)) + 
  geom_abline(slope=1, linetype = 'dashed', color = 'brown', size = 0.7) +
  #scale_fill_discrete(values = c("#d8b365", "#f5f5f5", "#5ab4ac")) +
  scale_color_manual(values=c('red','green','orange')) +
  ylim(-1.5,2) + 
  xlim(-4,8) +
  theme_bw() +
  geom_text_repel(
    data = subset(df, (Group != 'Others') & (RNAfc > -2.3 | RNAfc < -1.7 | abs(TEfc > 1.6))),
    aes(label = Symbol), size = 3,
    box.padding = unit(0.4, 'lines'),
    point.padding = unit(0.4, 'lines'),
    force = 3,
    segment.size = 0.2, segment.colour = 'grey50'
  ) + 
  ggtitle("TE fold change against RNA fold change") +
  theme(legend.position = c(0.8,0.9),
        #legend.position = 'none',
        legend.title = element_blank(),
        legend.background = element_rect(color = 'black',size = 0.5,linetype= 'solid'),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = rel(1.5),hjust = 0.5),
        axis.title = element_text(size = rel(1.25)),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size=14))
ggsave('scatter_L24_S15_feb/Both up down in TE and RNA t2_t1.png',width = 9, height = 7)




#anchor

L24S15 <- read.table('csv_raw/L24_against_S15_t1.csv', sep = ',',header = T)

L24S15 <- L24S15[order(L24S15$log2FoldChange),]
L24S15 <- L24S15[L24S15$padj<0.05,]
L24S15down <- L24S15[L24S15$log2FoldChange< -1,]
L24S15up <- tail(L24S15,nrow(L24S15down))

L24_more_S15 = is.element(TEfc$X,L24S15up$X)
L24_less_S15 = is.element(TEfc$X,L24S15down$X)

group = rep('Others',7990)

both_up <- L24_up_id & S15_up_id
group[both_up] <- "L24 up S15 up"

both_down <- L24_down_id & S15_down_id
group[both_down] <- "L24 down S15 down"

group[L24_more_S15] <- "L24 > S15"
group[L24_less_S15] <- "L24 < S15"

df = data.frame(TEfc = TEfc$log2FoldChange,
                RNAfc = RNAfc$log2FoldChange,
                Symbol = list_of_genes,
                X = TEfc$X,
                Group = group)
df$Group <- factor(df$Group, levels = c("L24 up S15 up","L24 down S15 down","L24 > S15","L24 < S15","Others"))

ggplot(df, aes(x=RNAfc, y=TEfc)) + 
  geom_point(data = df[df$Group == "Others",], alpha = 0.5, colour = "grey") + 
  geom_point(data = df[df$Group != "Others",], aes(colour = Group)) + 
  geom_abline(slope=1, linetype = 'dashed', color = 'brown', size = 0.7) +
  #scale_fill_discrete(values = c("#d8b365", "#f5f5f5", "#5ab4ac")) +
  scale_color_manual(values=c('red','green','orange','blue')) +
  ylim(-1.5,2) + 
  xlim(-4,8) +
  theme_bw() +
  geom_text_repel(
    data = subset(df, (Group != 'Others') & (RNAfc > 3 | RNAfc < -0.5 | TEfc > 1.1 | TEfc < -0.7)),
    aes(label = Symbol), size = 3,
    box.padding = unit(0.4, 'lines'),
    point.padding = unit(0.4, 'lines'),
    force = 3,
    segment.size = 0.2, segment.colour = 'grey50'
  ) + 
  ggtitle("TE fold change against RNA fold change") +
  theme(legend.position = c(0.8,0.9),
        #legend.position = 'none',
        legend.title = element_blank(),
        legend.background = element_rect(color = 'black',size = 0.5,linetype= 'solid'),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = rel(1.5),hjust = 0.5),
        axis.title = element_text(size = rel(1.25)),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size=14))
ggsave('scatter_L24_S15_feb/TE and RNA t2_t1.png',width = 9, height = 7)



## anchor
L24_col = rep('noChange',7990)
L24_col[L24_down_id] <- "down"
L24_col[L24_up_id] <- "up"

S15_col = rep('noChange',7990)
S15_col[S15_down_id] <- "down"
S15_col[S15_up_id] <- "up"

L24S15_col = rep('noChange',7990)
L24S15_col[L24_more_S15] <- "up"
L24S15_col[L24_less_S15] <- "down"

df = data.frame(Gene = list_of_genes,
                TEfc = TEfc$log2FoldChange,
                RNAfc = RNAfc$log2FoldChange,
                L24 = L24_col,
                S15 = S15_col,
                L24overS15= L24S15_col)

write.csv(df,"scatter_L24_S15_feb/tefc_rnafc_t2t1_full.csv", row.names = FALSE)

redundent = df$L24 == 'noChange' & df$S15 == 'noChange' & df$L24overS15 == 'noChange'
df_short = df[c(!redundent),]
write.csv(df_short,"scatter_L24_S15_feb/tefc_rnafc_t2t1_short.csv", row.names = FALSE)
