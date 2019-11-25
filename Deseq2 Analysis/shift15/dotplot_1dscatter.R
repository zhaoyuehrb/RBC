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

  cd_id = is.element(f_cd$X,commonGene)
  
  fp_specific = !fp_id
  
  fpGenes = f_fp$X[fp_specific]
  fpNames <- mapIds(org.Hs.eg.db,
                          keys = as.character(fpGenes),
                          column = 'SYMBOL',
                          keytype ='ACCNUM',
                          multiVals = 'first')
  
  df = data.frame(fpUTR = f_fp$log2FoldChange[fp_specific],
                  Symbol = fpNames)
  
  right <- unname(quantile(df$fpUTR,0.99))
  left <- unname(quantile(df$fpUTR,0.01))
  
  ggplot(df, aes(x=fpUTR,y=0)) + 
    geom_dotplot(color='white',alpha = 0.6,binwidth = 0.05,stackratio = 5,dotsize = 0.5) +
    #geom_abline(slope=1, linetype = 'dashed', color = 'brown', size = 0.7) +
    #ylim(min(df$fpUTR,df$cd),max(df$fpUTR,df$cd)) +
    #xlim(min(df$fpUTR,df$cd),max(df$fpUTR,df$cd)) +
    theme_bw() +
    geom_text_repel(
      #hjust = T,
      force = 3,
      direction = 'both',
      data = subset(df, (fpUTR>right) | (fpUTR<left)),
      aes(label = Symbol), size = 4,
      box.padding = unit(0.4, 'lines'),
      point.padding = unit(0.4, 'lines'),
      segment.size = 0.2, segment.colour = 'grey50') +
    ggtitle(file.label) +
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
  
  df = data.frame(cd = f_cd$log2FoldChange[cd_specific],
                  Symbol = cdNames)
  
  right <- unname(quantile(df$cd,0.997))
  left <- unname(quantile(df$cd,0.003))
  ggplot(df, aes(x=cd,y=0)) + 
    geom_dotplot(color="white",alpha = 0.6,binwidth = 0.03,stackratio = 1.2,dotsize = 0.5) +
    #geom_abline(slope=1, linetype = 'dashed', color = 'brown', size = 0.7) +
    #ylim(min(df$fpUTR,df$cd),max(df$fpUTR,df$cd)) +
    #xlim(min(df$fpUTR,df$cd),max(df$fpUTR,df$cd)) +
    theme_bw() +
    geom_text_repel(
      #hjust = T,
      force = 3,
      data = subset(df, (cd>right) | (cd<left)),
      aes(label = Symbol), size = 4,
      box.padding = unit(0.4, 'lines'),
      point.padding = unit(0.4, 'lines'),
      segment.size = 0.2, segment.colour = 'grey50') +
    ggtitle(file.label) +
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







