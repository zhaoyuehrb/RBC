library(DiffLogo)
library(seqLogo)
library(MotifDb)

er_dir <- 'deseq_wo_replace/pwms/'
cluster_dir <- 'kmeans_jm/pwm/'
createPlot <- function(pwm1,pwm2,cluster) {
  diffLogoObj = createDiffLogoObject(pwm1 = pwm1, pwm2 = pwm2)
  #diffLogoObj$ylim.negMax=-0.05
  #diffLogoObj$ylim.posMax=0.05
  widthToHeightRatio = 16/10
  size = length(pw1) * 2 / 4
  resolution = 200
  width = size * widthToHeightRatio
  height = size
  fileName = paste('kmeans_jm/',batch,'/',batch,' vs cluster',cluster,'.png',sep='')
  png(filename = fileName, res = resolution,width = width * resolution, height = height * resolution)
  #diffLogoFromPwm(pwm1 = pwm1, pwm2 = pwm2)
  diffLogo(diffLogoObj)
  dev.off()
}
batches <- clist <- c("ER_L24_t1", "ER_S15_t0","ER_S15_t1","ER_S15_t2")

for (batch in batches) {
  for (cluster in c(1:8)) {
  #dir.create(file.path(getwd(),'difflogo/',batch), showWarnings = FALSE)
  
  
  
  pwm_en = read.csv(paste(er_dir,batch,'/',batch,'_enriched.csv',sep=''),row.names=1)
  pwm_cluster = read.csv(paste(cluster_dir,'cluster',cluster,'_pwm.csv',sep=''),row.names=1)
  
  
  createPlot(pwm1 = pwm_en,pwm2 = pwm_cluster,cluster = cluster)
}
}



diffLogoObj = createDiffLogoObject(pwm1 = pw1, pwm2 = pw2)
diffLogoObj$ylim.negMax=-0.05
diffLogoObj$ylim.posMax=0.05
widthToHeightRatio = 16/10
size = length(pw1) * 2 / 4
resolution = 200
width = size * widthToHeightRatio
height = size
fileName = "test.png"
png(filename = fileName, res = resolution,width = width * resolution, height = height * resolution)
diffLogoFromPwm(pwm1 = pw1, pwm2 = pw2)
dev.off()

