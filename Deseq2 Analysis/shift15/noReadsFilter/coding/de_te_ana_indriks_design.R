
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
 
}

###deseq at individual time points for TE, controlling for the baseline
res_1_te_0 <- filter_sig_res(dds, c("assayRPF.timet1"))
res_2_te_0 <- filter_sig_res(dds, c("assayRPF.timet2"))

write.csv(res_1_te_0,"csv_raw/TE_t1_t0.csv")
write.csv(res_2_te_0,"csv_raw/TE_t2_t0.csv")








colData$time <- factor(colData$time, levels = c("t1","t2","t0"))
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
  
}

###deseq at individual time points for TE, controlling for the baseline
res_2_te_1 <- filter_sig_res(dds, c("assayRPF.timet2"))

write.csv(res_2_te_1,"csv_raw/TE_t2_t1.csv")


