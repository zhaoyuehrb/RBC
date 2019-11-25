library(dabestr)
library(ggplot2)
library(org.Hs.eg.db)
library(AnnotationDbi)


files <- list.files(path="csv_raw/", pattern="*.csv", full.names=TRUE, recursive=FALSE)

for (i in 1:21) {
  
  file.name = strsplit(files[i],'/')
  file.label = tail(file.name[[1]],1)
  file.label = strsplit(file.label,'.csv')
  
  
  deseq_out <- read.table(files[i], sep = ',',header = T) 
  deseq_out$GeneName <- mapIds(org.Hs.eg.db,
                          keys = as.character(deseq_out$X),
                          column = 'SYMBOL',
                          keytype ='ACCNUM',
                          multiVals = 'first')
  
  genes_table <- read.table('/Users/zhaoy/Documents/HGLab/Gene Groups/fpUTR_with_IREs_medium.csv',sep=',',header=T)
  idxs <-   deseq_out$GeneName %in% genes_table$Gene
  nidxs <-!idxs
  
  dat <- data.frame(matrix(ncol = 2, nrow = 5522))
  x <- c("Groups","logFC")
  colnames(dat) <- x
  dat$Groups[idxs] <- "fpUTR_with_IREs_medium"
  dat$Groups[!idxs] <- "Others"
  dat$logFC <- deseq_out$log2FoldChange
  two.group.unpaired <-
    dat %>%
    dabest(Groups,logFC,
           idx=c("Others","fpUTR_with_IREs_medium"),
           paired=FALSE,reps =10000)
  
  #two.group.unpaired
  
  
  dir.create(file.path(getwd(), "estimation stats/fpUTR_with_IREs_medium"), showWarnings = FALSE)
  outname = paste("estimation stats/fpUTR_with_IREs_medium/",file.label,"dabest.jpg")
  jpeg(outname)
  print(plot(two.group.unpaired,rawplot.ylim = c(-1.5,3)))
  dev.off()

  print(1+1)
}


for (i in 1:21) {
  
  file.name = strsplit(files[i],'/')
  file.label = tail(file.name[[1]],1)
  file.label = strsplit(file.label,'.csv')
  
  
  deseq_out <- read.table(files[i], sep = ',',header = T) 
  deseq_out$GeneName <- mapIds(org.Hs.eg.db,
                               keys = as.character(deseq_out$X),
                               column = 'SYMBOL',
                               keytype ='ACCNUM',
                               multiVals = 'first')
  
  genes_table <- read.table('/Users/zhaoy/Documents/HGLab/Gene Groups/tpUTR_with_IREs_medium.csv',sep=',',header=T)
  idxs <-   deseq_out$GeneName %in% genes_table$Gene
  nidxs <-!idxs
  
  dat <- data.frame(matrix(ncol = 2, nrow = 5522))
  x <- c("Groups","logFC")
  colnames(dat) <- x
  dat$Groups[idxs] <- "tpUTR_with_IREs_medium"
  dat$Groups[!idxs] <- "Others"
  dat$logFC <- deseq_out$log2FoldChange
  two.group.unpaired <-
    dat %>%
    dabest(Groups,logFC,
           idx=c("Others","tpUTR_with_IREs_medium"),
           paired=FALSE,reps =10000)
  
  #two.group.unpaired
  
  
  dir.create(file.path(getwd(), "estimation stats/tpUTR_with_IREs_medium"), showWarnings = FALSE)
  outname = paste("estimation stats/tpUTR_with_IREs_medium/",file.label,"dabest.jpg")
  jpeg(outname)
  print(plot(two.group.unpaired,rawplot.ylim = c(-1.5,3)))
  dev.off()
  
  print(1+1)
}