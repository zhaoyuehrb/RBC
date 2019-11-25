library(ggplot2)
library(ggrepel)
library(org.Hs.eg.db)

fp_files <- list.files(path="5pUTR/csv_raw/", pattern="*.csv", full.names=TRUE, recursive=FALSE)
cd_files <- list.files(path="coding/csv_raw/", pattern="*.csv", full.names=TRUE, recursive=FALSE)

for (i in c(1:15)) {
  
  file.name = strsplit(fp_files[i],'/')
  file.label = tail(file.name[[1]],1)
  file.label = strsplit(file.label,'.csv')
  
  f_fp <- read.table(fp_files[i], sep = ',',header = T) 
  f_cd <- read.table(cd_files[i], sep = ',',header = T) 
  
  commonGene = intersect(f_fp$X,f_cd$X)
#  list_of_genes <- mapIds(org.Hs.eg.db,
#                          keys = as.character(commonGene),
#                          column = 'SYMBOL',
#                          keytype ='ACCNUM',
#                          multiVals = 'first')
  fp_id = is.element(f_fp$X,commonGene)
  cd_id = is.element(f_cd$X,commonGene)
  
  df = data.frame(fpUTR = f_fp$log2FoldChange[fp_id],
                  cd = f_cd$log2FoldChange[cd_id],
                  Symbol = list_of_genes)
  
  right <- unname(quantile(df$fpUTR,0.998))
  left <- unname(quantile(df$fpUTR,0.002))
  top <- unname(quantile(df$cd,0.998))
  bottom <- unname(quantile(df$cd,0.002))
  
  ggplot(df, aes(x=fpUTR, y=cd)) + 
    geom_point(colour = "cyan2",alpha = 0.6) +
    geom_abline(slope=1, linetype = 'dashed', color = 'brown', size = 0.7) +
    ylim(min(df$fpUTR,df$cd),max(df$fpUTR,df$cd)) +
    xlim(min(df$fpUTR,df$cd),max(df$fpUTR,df$cd)) +
    theme_bw() +
    geom_text_repel(
      data = subset(df, (fpUTR>right) | (fpUTR<left) | (cd>top) | (cd<bottom)),
      aes(label = Symbol), size = 3,
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
  
  save_dir = 'compareBoth/'
  filename = paste(save_dir,file.label,'.png',sep = '')
  ggsave(filename,width = 9, height = 7)
}
  
  
  
  
  
  
  
