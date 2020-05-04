import loadCollectionsModule_122710 as lcm
import utilityModule as um
import sys, time
################################################################################
def spliceInterval(locus, gene):
    intervalLength = locus.len() - 1
    for intron in gene.introns():
        if locus.contains(intron):
            intervalLength -= intron.len()
    return intervalLength
############################################################################
def spliceCDExons(gene):
    cdExons = gene.cdExons()

    sumExonLength = 0
    for exon in cdExons:
        exonLength = exon.len()
        sumExonLength += exonLength
    return sumExonLength
#################################################################################################
#########################################################################################
def getAsitePosition_genomicHit(hit, gene, exptType):
    if gene.sense() == '+':
        interval = um.Locus(gene.chr(), gene.cdLocus().start(), hit.start(), gene.sense())
        if hit.start() >= gene.cdLocus().start():
            hitStartPosition = spliceInterval(interval, gene)
        else:
            hitStartPosition = -(spliceInterval(interval, gene))
    else:
        interval = um.Locus(gene.chr(), hit.end(), gene.cdLocus().end(), gene.sense())
        if hit.end() <= gene.cdLocus().end():
            hitStartPosition = spliceInterval(interval, gene)
        else:
            hitStartPosition = -(spliceInterval(interval, gene))
    if exptType == 'RPF':
        shiftedHitStartPosition = hitStartPosition + 15
    elif exptType == 'RNA':
        shiftedHitStartPosition = hitStartPosition
    return shiftedHitStartPosition
##################################################################################################
def getAsitePosition_refseqHit(refseqHit, gene, exptType):
    if gene.sense() == '+':
        interval = um.Locus(gene.chr(), gene.cdLocus().start(), refseqHit.genomeStart(), gene.sense())
        if refseqHit.genomeStart() >= gene.cdLocus().start():
            hitStartPosition = spliceInterval(interval, gene)
        else:
            hitStartPosition = -(spliceInterval(interval, gene))
    else:
        interval = um.Locus(gene.chr(), refseqHit.genomeStart(), gene.cdLocus().end(), gene.sense())
        if refseqHit.genomeStart() <= gene.cdLocus().end():
            hitStartPosition = spliceInterval(interval, gene)
        else:
            hitStartPosition = -(spliceInterval(interval, gene))
    if exptType == 'RPF':
        shiftedHitStartPosition = hitStartPosition + 15
    elif exptType == 'RNA':
        shiftedHitStartPosition = hitStartPosition
    return shiftedHitStartPosition
#################################################################################
#############################################################################
def loadAccNumToRawReads(filePath, threshold):
    infile = open(filePath)
    infileLines = infile.readlines()
    infile.close()

    accNumToRawReads = {}
    for line in infileLines[1:-1]:
        line = line.rstrip('\n')
        cols = line.split('\t')
        accNum = cols[0]
        WTReads = float(cols[2])
        KOReads = float(cols[1])
        if KOReads >= threshold:
            accNumToRawReads[accNum] = (KOReads, WTReads)
    totalMappedReads_KO, totalMappedReads_WT = map(float, infileLines[-1].split('\t'))

    print 'Number of genes loaded: ', len(accNumToRawReads)
    print 'TotalMappedReads_KO for this set: ', totalMappedReads_KO
    print 'TotalMappedReads_WT for this set: ', totalMappedReads_WT
    return accNumToRawReads, totalMappedReads_KO, totalMappedReads_WT
################################################################################
rootdirectory = sys.argv[1]
lenRange = '18,37'
barcode = rootdirectory[:-1].split('/')[-1]
#exptType = transfection.split('_')[-1]
refseqDirectory = rootdirectory + barcode + '_seed29_refseq/'
##mir = sys.argv[2]
##footprint_threshold = int(sys.argv[2])
##rnaseq_threshold = int(sys.argv[4])
exptType = sys.argv[2]

##########################################################################
refPath = '/home/hguo/Documents/annotations/hg19/refFlat/refFlat_240118.txt'
splicedGeneORFpath = '/home/hguo/Documents/annotations/hg19/refFlat/splicedGenes_ORF_240118.txt'
chrList = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']
filteredGeneList = lcm.getFilteredGeneList_human(refPath, splicedGeneORFpath, chrList)
##############################################################################

##refseqDirectoryCore = (rootdirectory[:-1].split('/')[-1])+'_refseq/'
genomeDirectory = rootdirectory#+genomeDirectoryCore
##refseqDirectory = rootdirectory+refseqDirectoryCore

#chrList = ['12']#, '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']
chrList = map(lambda i: 'chr'+i, chrList)   ## no chrM !! 

chromToGenes = {}
notOnChrom = []

## populate chromToGenes dict
for chrom in chrList:
    chromToGenes[chrom] = []
for gene in filteredGeneList:
    chrom = gene.chr()
    if chrom in chromToGenes:
        chromToGenes[chrom].append(gene)
    else:
        notOnChrom.append(chrom)
        
numGenesInDict = 0
for chrom, genes in chromToGenes.iteritems():
    numGenesInDict += len(genes)
    #print chrom, len(genes)
print 'Number of genes in chromToGenes: ', numGenesInDict
print 'Number of entries in chromToGenes: ', len(chromToGenes)

## load the refseqHits (genomic hits will be loaded by chrom)
uniqRefseqHits = lcm.getAllUniqRefseqHits(refseqDirectory, lenRange)

#############################################################################
################################################################################
#############################################################################
## initialize frameToReads dict
frameToReads = {}
for i in xrange(3):
    frameToReads[i] = 0
#############################################################################

numGenesContribute = 0
for chrom in chrList:
    print '=========== '+chrom+' ===================', time.asctime()
    uniqHitsCollection = lcm.getUniqGenomicHitsByChrom(chrom, genomeDirectory, lenRange)
    print 'Number of uniq hits that map to '+chrom+': ', len(uniqHitsCollection)

    exonHitsCollection, intronHitsCollection = lcm.getUniqTXExonIntronHits(uniqHitsCollection, filteredGeneList)    

    ##################################################################
    ##################################################################
    for gene in chromToGenes[chrom]:
        accNum = gene.accNum()

        ORF_length = spliceCDExons(gene)
        numGenesContribute += 1
        if accNum in exonHitsCollection:                    
            exonHits = exonHitsCollection[accNum]
            for exonHit in exonHits:
                shiftedHitStart = getAsitePosition_genomicHit(exonHit, gene, exptType)
                fraction = shiftedHitStart/float(ORF_length)
                if 0 <= fraction < 1:
                    ## this hit is within the coding region
                    frame = shiftedHitStart%3
                    if frame == 0:
                        frameToReads[0] += exonHit.readNum()
                    elif frame == 1:
                        frameToReads[1] += exonHit.readNum()
                    elif frame == 2:
                        frameToReads[2] += exonHit.readNum()

        if accNum in uniqRefseqHits:
            refseqHits = uniqRefseqHits[accNum]
            for refseqHit in refseqHits:
                shiftedHitStart = getAsitePosition_refseqHit(refseqHit, gene, exptType)
                fraction = shiftedHitStart/float(ORF_length)
                if 0 <= fraction < 1:
                ## this hit is within the coding region
                    frame = shiftedHitStart%3
                    if frame == 0:
                        frameToReads[0] += refseqHit.readNum()
                    elif frame == 1:
                        frameToReads[1] += refseqHit.readNum()
                    elif frame == 2:
                        frameToReads[2] += refseqHit.readNum()

#################################################################################
print 'Number of genes that contributed: ', numGenesContribute
################################################################################
####
## whether divide by totalMappedReads or not,
## fraction for each frame is going to be the same
####
###########################################################################
print '============================================'
for i in xrange(3):
    print i, frameToReads[i]
    
import csv
key = rootdirectory[:-1].split('/')[-2] + '_' + barcode
fname = key + '_frame.csv'
with open(fname,"w") as f:
    cr = csv.writer(f,delimiter=",",lineterminator="\n")
    for i in xrange(3):
        cr.writerow([i,frameToReads[i]])
