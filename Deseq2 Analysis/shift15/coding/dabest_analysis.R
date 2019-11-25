library(dabestr)
library(ggplot2)
dat = read.table('csv_dabest/ER_S15_t0_erythroid.csv',sep=',',header = T)
two.group.unpaired <-
  dat %>%
  dabest(Groups,logFC,
         idx=c("non-Erythroids","Erythroids"),
         paired=FALSE,reps =10000)

two.group.unpaired

jpeg("csv_dabest/ER_S15_t0_erythroid_dabest.jpg")
plot(two.group.unpaired,rawplot.ylim = c(-1.5,3))
dev.off()

dat = read.table('csv_dabest/ER_S15_t1_erythroid.csv',sep=',',header = T)
two.group.unpaired <-
  dat %>%
  dabest(Groups,logFC,
         idx=c("non-Erythroids","Erythroids"),
         paired=FALSE,reps =10000)

two.group.unpaired

jpeg("csv_dabest/ER_S15_t1_erythroid_dabest.jpg")
plot(two.group.unpaired,rawplot.ylim = c(-1.5,3))
dev.off()

dat = read.table('csv_dabest/ER_S15_t2_erythroid.csv',sep=',',header = T)
two.group.unpaired <-
  dat %>%
  dabest(Groups,logFC,
         idx=c("non-Erythroids","Erythroids"),
         paired=FALSE,reps =10000)

two.group.unpaired

jpeg("csv_dabest/ER_S15_t2_erythroid_dabest.jpg")
plot(two.group.unpaired,rawplot.ylim = c(-1.5,3))
dev.off()