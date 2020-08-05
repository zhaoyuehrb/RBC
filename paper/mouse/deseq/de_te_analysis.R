library(DESeq2)
library(IHW)

#load read count table
#df_te <- read.table('/dir', sep = '\t',header = T) # dir = path where df_te is stored
df_te <- read.table('../collate_rpkms.csv', sep = ',',header = T) 
#df_te <- read.table('../collate_rpkms_noFil.csv', sep = ',',header = T) 
rownames(df_te) <- df_te[,1]
df_te[,1]<-NULL
df_te <- df_te*10


###call design matrix
colData = read.table('combined_design.txt',sep=',',header=TRUE) # dir = path where de_te_design is stored

###deseq analysis
###use a design formula that models the assay different at t0, the diff over time,
###and any assay-specific differences over time (the interaction term assay:time)
design = ~ assay
colData$assay <- factor(colData$assay, levels = c("rpl10a_total","rpl25_total","rpl22_total","rpl10a","rpl25","rpl22"))
dds <- DESeqDataSetFromMatrix(countData = df_te, colData = colData,design = design)

#run DESeq
dds <- DESeq(dds)

alpha <- 0.05
res <- results(dds,alpha = alpha)
resultsNames(dds)
test_design <- model.matrix(design, colData)

###function to compare different combinations of factors
filter_sig_res <- function(dds, factors) {
  res_df <- results(dds, contrast = (factors),test = "Wald",alpha = alpha,cooksCutoff = FALSE,independentFiltering = FALSE)
  
  summary(res_df)
  res_df <- na.omit(res_df)
}
#### ER
res_rpl10 <- filter_sig_res(dds,c("assay","rpl10a","rpl10a_total"))


###########################
colData$assay <- factor(colData$assay, levels = c("rpl25_total","rpl10a_total","rpl22_total","rpl10a","rpl25","rpl22"))
dds <- DESeqDataSetFromMatrix(countData = df_te, colData = colData,design = design)

#run DESeq
dds <- DESeq(dds)

alpha <- 0.05
res <- results(dds,alpha = alpha)
resultsNames(dds)
test_design <- model.matrix(design, colData)

###function to compare different combinations of factors
filter_sig_res <- function(dds, factors) {
  res_df <- results(dds, contrast = (factors),test = "Wald",alpha = alpha,cooksCutoff = FALSE,independentFiltering = FALSE)
  
  summary(res_df)
  res_df <- na.omit(res_df)
}


#### ER
res_rpl25 <- filter_sig_res(dds,c("assay","rpl25","rpl25_total"))




###########################
colData$assay <- factor(colData$assay, levels = c("rpl22_total","rpl25_total","rpl10a_total","rpl10a","rpl25","rpl22"))
dds <- DESeqDataSetFromMatrix(countData = df_te, colData = colData,design = design)

#run DESeq
dds <- DESeq(dds)

alpha <- 0.05
res <- results(dds,alpha = alpha)
resultsNames(dds)
test_design <- model.matrix(design, colData)

###function to compare different combinations of factors
filter_sig_res <- function(dds, factors) {
  res_df <- results(dds, contrast = (factors),test = "Wald",alpha = alpha,cooksCutoff = FALSE,independentFiltering = FALSE)
  
  summary(res_df)
  res_df <- na.omit(res_df)
}


#### ER
res_rpl22 <- filter_sig_res(dds,c("assay","rpl22","rpl22_total"))










write.csv(res_rpl10,"res_rpl10.csv")
write.csv(res_rpl25,"res_rpl25.csv")
write.csv(res_rpl22,"res_rpl22.csv")

