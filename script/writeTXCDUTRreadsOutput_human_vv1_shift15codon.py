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

###############################################################################
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
################################################################################
################################################################################
##########################################################################
rootdirectory = sys.argv[1]
lenRange = sys.argv[2]
exptType = sys.argv[3]
refdir = sys.argv[4]

##########################################################################
refPath = '/home/hguo/Documents/annotations/hg19/refFlat/refFlat_240118.txt'
splicedGeneORFpath = '/home/hguo/Documents/annotations/hg19/refFlat/splicedGenes_ORF_240118.txt'
chrList = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'M']
filteredGeneList = lcm.getFilteredGeneList_human(refPath, splicedGeneORFpath, chrList)
##############################################################################

#genomeDirectoryCore = (rootdirectory[:-1].split('/')[-1])+'_genome/'
#refseqDirectoryCore = (rootdirectory[:-1].split('/')[-1])+'_refseq/'
genomeDirectory = rootdirectory#+genomeDirectoryCore
refseqDirectory = rootdirectory+refdir#refseqDirectoryCore

chrList = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']
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
## first, populate the dict with all available genes
geneToExonRefseqIntronHits = {}
for chrom in chrList:
    for gene in chromToGenes[chrom]:
        accNum = gene.accNum()
        if not geneToExonRefseqIntronHits.has_key(accNum):
            geneToExonRefseqIntronHits[gene] = [0, 0, 0, 0, 0, 0, accNum]
            
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

                shiftedHitStart = getAsitePosition_genomicHit(exonHit, gene, exptType)
                fraction = shiftedHitStart/float(ORF_length)
                if 0 <= fraction < 1 and shiftedHitStart >=45:
                    ## this hit is within the coding region
                    geneToExonRefseqIntronHits[gene][3] += exonHit.readNum()
                elif fraction < 0:
                    ## this hit is in the fpUTR
                    geneToExonRefseqIntronHits[gene][4] += exonHit.readNum()
                elif fraction >= 1:
                    ## this hit is in the tpURE
                    geneToExonRefseqIntronHits[gene][5] += exonHit.readNum()
                    
                
        ## refseq column
        if accNum in uniqRefseqHits:
            refseqHits = uniqRefseqHits[accNum]
            for refseqHit in refseqHits:
                geneToExonRefseqIntronHits[gene][1] += refseqHit.readNum()

                shiftedHitStart = getAsitePosition_refseqHit(refseqHit, gene, exptType)
                fraction = shiftedHitStart/float(ORF_length)
                if 0 <= fraction < 1 and shiftedHitStart >=45:
                    ## this hit is within the coding region
                    geneToExonRefseqIntronHits[gene][3] += refseqHit.readNum()
                elif fraction < 0:
                    ## this hit is in the fpUTR
                    geneToExonRefseqIntronHits[gene][4] += refseqHit.readNum()
                elif fraction >= 1:
                    ## this hit is in the tpURE
                    geneToExonRefseqIntronHits[gene][5] += refseqHit.readNum()

                    
        ## intron column
        if accNum in intronHitsCollection:
            intronHits = intronHitsCollection[accNum]
            for intronHit in intronHits:
                geneToExonRefseqIntronHits[gene][2] += intronHit.readNum()


## after going through all chr's, all FOUR columns shld be populated
## G+R TX reads need to be summed here, but G+R CD reads shld already be summed
date = rootdirectory[:-1].split('/')[-2]
laneIdentity = rootdirectory[:-1].split('/')[-1]
##outDirectory = '/lab/bartel4_ata/huili/122610/mouse_querySeq/txCDutrReadsOutput/'
outDirectory = rootdirectory + '../'
outfileString = outDirectory+'gene_TXCDUTR_ReadOutput_'+date+'_'+laneIdentity+'_shift15.txt'
outfile = open(outfileString, 'w')
firstLineList = ['AccNum', 'ExonReads', 'RefseqReads', 'IntronReads', 'GenomeRefseqReads', 'CD_GenomeRefseqReads', 'fpUTR_reads', 'tpUTR_reads']
firstLineString = '\t'.join(firstLineList)
outfile.write(firstLineString+'\n')

linesWritten = 1
for gene, readsList in geneToExonRefseqIntronHits.iteritems():
    outList = [readsList[-1]]
    outList.extend(readsList[:3])
    outList.append(readsList[0]+readsList[1])
    outList.extend(readsList[3:-1])
    outList = map(lambda i: str(i), outList)
    outListString = '\t'.join(outList)
    outfile.write(outListString+'\n')
    linesWritten += 1
outfile.close()
print 'Number of lines written: ', linesWritten
