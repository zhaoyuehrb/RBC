library(ggplot2)
library(ggrepel)
library(org.Hs.eg.db)

fp_files <- list.files(path="5pUTR/csv_raw/", pattern="*.csv", full.names=TRUE, recursive=FALSE)
cd_files <- list.files(path="coding/csv_raw/", pattern="*.csv", full.names=TRUE, recursive=FALSE)

for (i in c(1:21)) {
  
  file.name = strsplit(fp_files[i],'/')
  file.label = tail(file.name[[1]],1)
  file.label = strsplit(file.label,'.csv')
  
  f_fp <- read.table(fp_files[i], sep = ',',header = T) 
  f_cd <- read.table(cd_files[i], sep = ',',header = T) 
  
  commonGene = intersect(f_fp$X,f_cd$X)
  
  fp_id = is.element(f_fp$X,commonGene)
  cd_id = is.element(f_cd$X,commonGene)
  
  fp_specific = !fp_id
  
  fpGenes = f_fp$X[fp_specific]
  fpNames <- mapIds(org.Hs.eg.db,
                    keys = as.character(fpGenes),
                    column = 'SYMBOL',
                    keytype ='ACCNUM',
                    multiVals = 'first')
  
  df = data.frame(fpUTR.logFC = f_fp$log2FoldChange[fp_specific],
                  Symbol = fpNames)
  
  pos <- position_jitter(width = 0.01, seed = 2)
  
  right <- unname(quantile(df$fpUTR.logFC,0.98))
  left <- unname(quantile(df$fpUTR.logFC,0.02))
  
  ggplot(df, aes(x=fpUTR.logFC,y=0)) + 
    geom_jitter(
      data = subset(df, (fpUTR.logFC<=right) & (fpUTR.logFC>=left)),
      colour = "#878787",
    ) +
    geom_jitter(
      position=pos,
      data = subset(df, (fpUTR.logFC>right) | (fpUTR.logFC<left)),
      colour = "red",
    ) +
    theme_bw() +
    geom_text_repel(
      position=pos,
      data = subset(df, (fpUTR.logFC>right) | (fpUTR.logFC<left)),
      aes(label = Symbol), size = 4,
      box.padding = unit(0.4, 'lines'),
      point.padding = unit(0.4, 'lines'),
      segment.size = 0.2, segment.colour = 'grey50') +
    ggtitle(paste(file.label,"fpUTR logFC jitter plot",sep=" ")) +
    theme(legend.position = c(0.15,0.9),
          #legend.position = 'none',
          legend.title = element_blank(),
          legend.background = element_rect(color = 'black',size = 0.5,linetype= 'solid'),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = rel(1.5),hjust = 0.5),
          axis.title = element_text(size = rel(1.25)),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size=14))
  
  save_dir = 'fp specific/'
  filename = paste(save_dir,file.label,'.png',sep = '')
  ggsave(filename,width = 9, height = 7)
  
  
  cd_specific = !cd_id
  
  cdGenes = f_cd$X[cd_specific]
  cdNames <- mapIds(org.Hs.eg.db,
                    keys = as.character(cdGenes),
                    column = 'SYMBOL',
                    keytype ='ACCNUM',
                    multiVals = 'first')
  
  df = data.frame(cdUTR.logFC = f_cd$log2FoldChange[cd_specific],
                  Symbol = cdNames)
  
  pos <- position_jitter(width = 0.0001, seed = 2)
  
  right <- unname(quantile(df$cdUTR.logFC,0.995))
  left <- unname(quantile(df$cdUTR.logFC,0.005))
  
  ggplot(df, aes(x=cdUTR.logFC,y=0)) + 
    geom_jitter(
      data = subset(df, (cdUTR.logFC<=right) & (cdUTR.logFC>=left)),
      colour = "#878787",
    ) +
    geom_jitter(
      position=pos,
      data = subset(df, (cdUTR.logFC>right) | (cdUTR.logFC<left)),
      colour = "red",
    ) +
    theme_bw() +
    geom_text_repel(
      position=pos,
      data = subset(df, (cdUTR.logFC>right) | (cdUTR.logFC<left)),
      aes(label = Symbol), size = 4,
      box.padding = unit(0.4, 'lines'),
      point.padding = unit(0.4, 'lines'),
      segment.size = 0.2, segment.colour = 'grey50') +
    ggtitle(paste(file.label,"cd logFC jitter plot",sep=" ")) +
    theme(legend.position = c(0.15,0.9),
          #legend.position = 'none',
          legend.title = element_blank(),
          legend.background = element_rect(color = 'black',size = 0.5,linetype= 'solid'),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = rel(1.5),hjust = 0.5),
          axis.title = element_text(size = rel(1.25)),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size=14))
  
  save_dir = 'cd specific/'
  filename = paste(save_dir,file.label,'.png',sep = '')
  ggsave(filename,width = 9, height = 7)
}