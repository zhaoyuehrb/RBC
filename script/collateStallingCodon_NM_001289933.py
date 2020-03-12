import loadCollectionsModule_122710 as lcm
import utilityModule as um
import sys, time
import os

#script_dir=os.path.dirname(__file__)
genomeDirectory =sys.argv[1]
outDirectory = sys.argv[2]  # '../../codonFreqWtAcc/'
lenRange = sys.argv[3]      #'30,31'
cutRange = sys.argv[4]      #'15,18' for shift 15 bits to the right
cutStr = cutRange.split(',')
left=int(cutStr[0])
right=int(cutStr[1])


################################################################################
def spliceInterval(locus, gene):
    intervalLength = locus.len() - 1
    for intron in gene.introns():
        if locus.contains(intron):
            intervalLength -= intron.len()
    return intervalLength
############################################################################

############################################################################
def spliceCDExons(gene):
    cdExons = gene.cdExons()

    sumExonLength = 0
    for exon in cdExons:
        exonLength = exon.len()
        sumExonLength += exonLength
    return sumExonLength



############################################################################
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





##########################################################################
refPath = '/home/hguo/Documents/annotations/hg19/refFlat/refFlat_240118.txt'
splicedGeneORFpath = '/home/hguo/Documents/annotations/hg19/refFlat/splicedGenes_ORF_240118.txt'
chrList = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']
## all genes
filteredGeneList = lcm.getFilteredGeneList_human(refPath, splicedGeneORFpath, chrList)
##############################################################################

f = open('background geneList.txt','r')
bgFile = f.readlines()
f.close()
selectedGenes = str(bgFile[0]).split('\r')
selectedGenes = ['NM_001289933']

selectGeneList = []
for gene in filteredGeneList:
	if gene._accNum in selectedGenes:
		selectGeneList.append(gene)
print 'Select Genes: ' , len(selectGeneList)

chrList = map(lambda i: 'chr'+i, chrList)
chromToGenes = {}
notOnChrom = []
## populate chromToGenes dict
for chrom in chrList:
    chromToGenes[chrom] = []
for gene in selectGeneList:
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
print 'Number of lines in notOnChrom: ', len(notOnChrom)


geneToExonRefseqIntronHits = {}
for chrom in chrList:
    for gene in chromToGenes[chrom]:
        accNum = gene.accNum()
        if not geneToExonRefseqIntronHits.has_key(accNum):
            geneToExonRefseqIntronHits[gene] = [0, 0, 0, 0, 0, 0, accNum]



outfileString = genomeDirectory+outDirectory+genomeDirectory[-13:-1]+' codon cut at '+cutStr[0]+cutStr[1]+'.txt'
outfile = open(outfileString, 'w')
firstLineList = ['Codon', 'Reads','AccNum']
firstLineString = ','.join(firstLineList)
outfile.write(firstLineString+'\n')
lineWritten = 1  

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
        ## exon column
        if accNum in exonHitsCollection:
            exonHits = exonHitsCollection[accNum]
            for exonHit in exonHits:
                geneToExonRefseqIntronHits[gene][0] += exonHit.readNum()

                shiftedHitStart = getAsitePosition_genomicHit(exonHit, gene, 'RPF')
                fraction = shiftedHitStart/float(ORF_length)
                #if 0 <= fraction < 1:
                if 0 <= fraction < 1 and shiftedHitStart >= 27:# and ORF_length - shiftedHitStart > 17:
                    ## this hit is within the coding region
					#print exonHit.readNum()
					#print exonHit.seq()
					outList = [exonHit.seq()[left:right], str(exonHit.readNum()), accNum]
					outListString = ','.join(outList)
					outfile.write(outListString+'\n')
					lineWritten +=1

outfile.close()
print "Lines written: ", lineWritten









