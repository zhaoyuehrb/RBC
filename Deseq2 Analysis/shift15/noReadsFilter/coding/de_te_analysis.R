library(DESeq2)
library(IHW)

#load read count table
#df_te <- read.table('/dir', sep = '\t',header = T) # dir = path where df_te is stored
df_te <- read.table('cdReadsNotFiltered.csv', sep = ',',header = T) 

rownames(df_te) <- df_te[,1]
df_te[,1]<-NULL

###call design matrix
colData = read.table('combined_design.txt',sep='\t',header=TRUE) # dir = path where de_te_design is stored

###deseq analysis
###use a design formula that models the assay different at t0, the diff over time,
###and any assay-specific differences over time (the interaction term assay:time)
design = ~ combined
#colData$assay <- factor(colData$assay, levels = c("RPF0")) #"RNA"))
dds <- DESeqDataSetFromMatrix(countData = df_te, colData = colData,design = design)

# ###run analysis to produce genes which at one or more time points after time 0 showed a TE effect
# reduced = ~ assay + time
# dds <- DESeq(dds, test='LRT', reduced = reduced)

#run DESeq
dds <- DESeq(dds)

alpha <- 0.05
res <- results(dds,alpha = alpha)
#res <- results(dds,alpha = alpha, filterFun = IHW)
resultsNames(dds)
test_design <- model.matrix(design, colData)

###function to compare different combinations of factors
filter_sig_res <- function(dds, factors) {
  res_df <- results(dds, contrast = (factors),test = "Wald",alpha = alpha,cooksCutoff = FALSE,independentFiltering = FALSE)
  
  summary(res_df)
  res_df <- na.omit(res_df)
  # #
  # res_df$SYMBOL <- mapIds(org.Mm.eg.db,
  #                         keys = as.character(rownames(res_df)),
  #                         column = 'SYMBOL',
  #                         keytype = 'ACCNUM',
  #                         multiVals = 'first')
  # ###remove statistically insignificant genes
  # res_df <- res_df[res_df$padj<alpha,]
  # res_df <- data.frame(ACCNUM=row.names(data.frame(res_df)), res_df,row.names = NULL)
  # # #
  # ###select columns Accnum and Log2FC
  # res_df <- res_df[,c("ACCNUM","log2FoldChange")]
  # res_df <- data.frame(res_df, row.names=NULL)
  # # 
  # # ###select AccNum, log2FC and padj
  # # res_df <- data.frame(ACCNUM = row.names(data.frame(res_df)), res_df, row.names = NULL)
  # # res_df <- res_df[,c("ACCNUM","log2FoldChange","padj")]
  
  ###add symbol and select padj
  # res_df <- data.frame(res_df)
  # rownames(res_df) <- NULL
  # res_df <- res_df[,c("SYMBOL","log2FoldChange")]
}


#### ER
res_0_L24 <- filter_sig_res(dds,c("combined","L240","RPF0"))

#res_0_L24_padj <- subset(res_0_L24,padj>0.05)

res_1_L24 <- filter_sig_res(dds,c("combined","L241","RPF1"))
res_2_L24 <- filter_sig_res(dds,c("combined","L242","RPF2"))

res_0_S15 <- filter_sig_res(dds,c("combined","S150","RPF0"))
res_1_S15 <- filter_sig_res(dds,c("combined","S151","RPF1"))
res_2_S15 <- filter_sig_res(dds,c("combined","S152","RPF2"))

res_0_S24 <- filter_sig_res(dds,c("combined","S240","RPF0"))
res_1_S24 <- filter_sig_res(dds,c("combined","S241","RPF1"))
res_2_S24 <- filter_sig_res(dds,c("combined","S242","RPF2"))

L24_0_S15 <- filter_sig_res(dds,c("combined","L240","S150"))
L24_1_S15 <- filter_sig_res(dds,c("combined","L241","S151"))
L24_2_S15 <- filter_sig_res(dds,c("combined","L242","S152"))






### RNA
res_1_against_0_RNA <- filter_sig_res(dds,c("combined","RNA1","RNA0"))
res_2_against_0_RNA <- filter_sig_res(dds,c("combined","RNA2","RNA0"))
res_2_against_1_RNA <- filter_sig_res(dds,c("combined","RNA2","RNA1"))


### RPF

res_1_against_0_RPF <- filter_sig_res(dds,c("combined","RPF1","RPF0"))
res_2_against_0_RPF <- filter_sig_res(dds,c("combined","RPF2","RPF0"))


## TE

res_0_TE <- filter_sig_res(dds,c("combined","RPF0","RNA0"))
res_1_TE <- filter_sig_res(dds,c("combined","RPF1","RNA1"))
res_2_TE <- filter_sig_res(dds,c("combined","RPF2","RNA2"))





write.csv(res_1_against_0_RNA,"csv_raw/RNA_T1_against_T0.csv")
write.csv(res_2_against_0_RNA,"csv_raw/RNA_T2_against_T0.csv")
write.csv(res_2_against_1_RNA,"csv_raw/RNA_T2_against_T1.csv")
write.csv(res_0_L24,"csv_raw/ER_L24_t0.csv")
write.csv(res_1_L24,"csv_raw/ER_L24_t1.csv")
write.csv(res_2_L24,"csv_raw/ER_L24_t2.csv")

write.csv(res_0_S15,"csv_raw/ER_S15_t0.csv")
write.csv(res_1_S15,"csv_raw/ER_S15_t1.csv")
write.csv(res_2_S15,"csv_raw/ER_S15_t2.csv")

write.csv(res_0_S24,"csv_raw/ER_S24_t0.csv")
write.csv(res_1_S24,"csv_raw/ER_S24_t1.csv")
write.csv(res_2_S24,"csv_raw/ER_S24_t2.csv")

write.csv(L24_0_S15,"csv_raw/L24_against_S15_t0.csv")
write.csv(L24_1_S15,"csv_raw/L24_against_S15_t1.csv")
write.csv(L24_2_S15,"csv_raw/L24_against_S15_t2.csv")
write.csv(res_1_against_0_RPF,"csv_raw/RPF_T1_against_T0.csv")
write.csv(res_2_against_0_RPF,"csv_raw/RPF_T2_against_T0.csv")
write.csv(res_0_TE,"csv_raw/TE_t0.csv")
write.csv(res_1_TE,"csv_raw/TE_t1.csv")
write.csv(res_2_TE,"csv_raw/TE_t2.csv")


### count stuff
#length(which(res_1_against_0_RPF$log2FoldChange> 1))
#length(which(res_1_against_0_RPF$log2FoldChange< -1))
#length(which(res_1_against_0_RPF$log2FoldChange>= -1 & res_1_against_0_RPF$log2FoldChange <=1))

