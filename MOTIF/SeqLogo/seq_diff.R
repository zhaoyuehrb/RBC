library(DiffLogo)
library(seqLogo)
library(MotifDb)

dir <- 'deseq_wo_replace/pwms/'
createPlot <- function(pwm1,pwm2,compare) {
  diffLogoObj = createDiffLogoObject(pwm1 = pwm1, pwm2 = pwm2)
  #diffLogoObj$ylim.negMax=-0.05
  #diffLogoObj$ylim.posMax=0.05
  widthToHeightRatio = 16/10
  size = length(pw1) * 2 / 4
  resolution = 200
  width = size * widthToHeightRatio
  height = size
  fileName = paste('difflogo/',batch,'/',batch,'_',compare,'.png',sep='')
  png(filename = fileName, res = resolution,width = width * resolution, height = height * resolution)
  #diffLogoFromPwm(pwm1 = pwm1, pwm2 = pwm2)
  diffLogo(diffLogoObj)
  dev.off()
}
batches <- clist <- c("ER_L24_t1", "ER_S15_t0","ER_S15_t1","ER_S15_t2")

for (batch in batches) {
  dir.create(file.path(getwd(),'difflogo/',batch), showWarnings = FALSE)
  
  
  
  pwm_en = read.csv(paste(dir,batch,'/',batch,'_enriched.csv',sep=''),row.names=1)
  pwm_de = read.csv(paste(dir,batch,'/',batch,'_depleted.csv',sep=''),row.names=1)
  pwm_bg = read.csv(paste(dir,batch,'/bg.csv',sep=''),row.names=1)
  pwm_left = read.csv(paste(dir,batch,'/',batch,'_leftstd_comb50.csv',sep=''),row.names=1)
  pwm_mid = read.csv(paste(dir,batch,'/',batch,'_median_comb50.csv',sep=''),row.names=1)
  pwm_right = read.csv(paste(dir,batch,'/',batch,'_rightstd_comb50.csv',sep=''),row.names=1)
  

  
  createPlot(pwm1 = pwm_en,pwm2 = pwm_de,compare = 'enrich vs deplete')
  createPlot(pwm1 = pwm_en,pwm2 = pwm_bg,compare = 'enrich vs background')
  createPlot(pwm1 = pwm_en,pwm2 = pwm_left,compare = 'enrich vs left')
  createPlot(pwm1 = pwm_en,pwm2 = pwm_mid,compare = 'enrich vs mode')
  createPlot(pwm1 = pwm_en,pwm2 = pwm_right,compare = 'enrich vs right')
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

