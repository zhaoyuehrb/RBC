library(DiffLogo)
library(seqLogo)
library(MotifDb)
## import motifs
hitIndeces = grep ('CTCF', values (MotifDb)$geneSymbol, ignore.case=TRUE)
list = as.list(MotifDb[hitIndeces])
sequenceCounts = as.numeric(values (MotifDb)$sequenceCount[hitIndeces])
names(sequenceCounts) = names(list)
pwm1 = reverseComplement(list$"Hsapiens-JASPAR_CORE-CTCF-MA0139.1"[, 2:18]) ## > pwm2 = list$"Hsapiens-jolma2013-CTCF"
n1 = sequenceCounts["Hsapiens-JASPAR_CORE-CTCF-MA0139.1"]
n2 = sequenceCounts["Hsapiens-jolma2013-CTCF"]
## DiffLogo can also handle motifs of different length
pwm_long = reverseComplement(list$"Hsapiens-JASPAR_CORE-CTCF-MA0139.1") ## reverse > pwm_short = list$"Hsapiens-jolma2013-CTCF"
pwm2 = list$"Hsapiens-jolma2013-CTCF"

## import five DNA motifs from matrix
motif_folder = "extdata/pwm"
motif_names_dna = c("H1-hESC", "MCF7", "HeLa-S3", "HepG2", "HUVEC")
motifs_dna = list()
for (name in motif_names_dna) {
fileName = paste(motif_folder,"/",name,".pwm",sep="")
file = system.file(fileName, package = "DiffLogo")
motifs_dna[[name]] = getPwmFromPwmFile(file)
}
sampleSizes_dna = c("H1-hESC"=100, "MCF7"=100, "HeLa-S3"=100, "HepG2"=100, "HUVEC"=100)
## import three DNA motifs from table
motif_folder = "extdata/alignments"
motif_names_dna2 = c("Mad", "Max", "Myc")
motifs_dna2 = list()
for (name in motif_names_dna2) {

fileName = paste(motif_folder,"/",name,".txt",sep="")
file = system.file(fileName, package = "DiffLogo")
motifs_dna2[[name]] = getPwmFromAlignmentFile(file)
}
> ## import three ASN motifs from fasta files
motif_folder = "extdata/alignments"
motif_names_asn = c("orderByRPF/seqs/RPF0_top.seq")
motifs_asn = list()
for (name in motif_names_asn) {
 fileName = paste(name,".fa",sep="")
 file = system.file(fileName, package = "DiffLogo")
 motifs_asn[[name]] = getPwmFromFastaFile(file, FULL_ALPHABET)
}



## plot custom sequence logo
par(mfrow=c(2,1), pin=c(3, 1), mar = c(2, 4, 1, 1))
DiffLogo::seqLogo(pwm = pwm1)
DiffLogo::seqLogo(pwm = pwm2, stackHeight = sumProbabilities)
par(mfrow=c(1,1), pin=c(1, 1), mar=c(5.1, 4.1, 4.1, 2.1))


## plot DiffLogo
diffLogoFromPwm(pwm1 = pwm1, pwm2 = pwm2)

## diffLogoFromPwm is a convenience function for
diffLogoObj = createDiffLogoObject(pwm1 = pwm1, pwm2 = pwm2)
diffLogo(diffLogoObj)
## mark symbol stacks with significant changes
diffLogoObj = enrichDiffLogoObjectWithPvalues(diffLogoObj, n1, n2)
diffLogo(diffLogoObj)
## plot DiffLogo for PWMs of different length
diffLogoFromPwm(pwm1 = pwm_long, pwm2 = pwm_short, align_pwms=TRUE)
