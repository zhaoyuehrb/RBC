genes_table <- read.table('Gene Groups/cytoskeleton.txt', sep = '\t',header = T) 
View(genes_table)

genes_table <- read.table('Gene Groups/CytoRPList.csv',sep=',',header=T)
View(genes_table)

genes_table <- read.table('Gene Groups/Erythroid_06032019List.txt',sep=',',header=T)
View(genes_table)
genes_table = as.data.frame(sapply(genes_table, toupper))


genes_table <- read.table('Gene Groups/HistoneList.csv',sep=',',header=T)
View(genes_table)

genes_table <- read.table('Gene Groups/lysosome.txt',sep='\t',header=T)
View(genes_table)

genes_table <- read.table('Gene Groups/mitoComplexAll.csv',sep='\t',header=T)
View(genes_table)

genes_table <- read.table('Gene Groups/MitoComplexI.txt',sep='\t',header=T)
View(genes_table)

genes_table <- read.table('Gene Groups/MitoComplexII.txt',sep='\t',header=T)
View(genes_table)

genes_table <- read.table('Gene Groups/MitoComplexIII.txt',sep='\t',header=T)
View(genes_table)

genes_table <- read.table('Gene Groups/MitoComplexIV.txt',sep='\t',header=T)
View(genes_table)

genes_table <- read.table('Gene Groups/MitoComplexV.txt',sep='\t',header=T)
View(genes_table)

genes_table <- read.table('Gene Groups/proteasome.txt',sep='\t',header=T)
View(genes_table)

genes_table <- read.table('Gene Groups/spliceosome1.csv',sep=',',header=T)
View(genes_table)

genes_table <- read.table('Gene Groups/autophagyRelated.csv',sep=',',header=T)
View(genes_table)

genes_table <- read.table('Gene Groups/MitoRPList.csv',sep=',',header=T)
View(genes_table)

genes_table <- read.table('Gene Groups/Erythroid_29102019.csv',sep=',',header=T)
View(genes_table)

genes_table <- read.table('Gene Groups/fpUTR_with_IREs.csv',sep=',',header=T)
View(genes_table)

genes_table <- read.table('Gene Groups/tpUTR_with_IREs.csv',sep=',',header=T)
View(genes_table)





genes_table <- read.table('csv_filtered/RNA_T2_against_T0_Up.csv',sep=',',header=T)
View(genes_table)
genes_table <- genes_table[order(-genes_table$log2FoldChange),]
genes_table <- genes_table[1:801,]

genes_table <- read.table('Gene Groups/custom/top10RPFt2.csv',sep=',',header=T)
View(genes_table)


