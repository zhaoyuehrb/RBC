
###deseq analysis
###use a design formula that models the assay different at t0, the diff over time,
###and any assay-specific differences over time (the interaction term assay:time)
design = ~ assay + time + assay:time
colData$assay <- factor(colData$assay, levels = c("RNA","RPF","L24","S15","S24"))
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

###function to compare different combinations of factors
filter_sig_res <- function(dds, factors) {
  res_df <- results(dds, contrast = list(factors),test = "Wald",alpha = alpha,cooksCutoff = FALSE,independentFiltering = FALSE)
  
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

###deseq at individual time points for TE, controlling for the baseline
res_1_te_0 <- filter_sig_res(dds, c("assayRPF.timet1"))
res_2_te_0 <- filter_sig_res(dds, c("assayRPF.timet2"))

write.csv(res_1_te_0,"csv_raw/TE_t1_t0.csv")
write.csv(res_2_te_0,"csv_raw/TE_t2_t0.csv")
