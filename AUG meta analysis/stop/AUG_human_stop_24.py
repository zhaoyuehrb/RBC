import loadCollectionsModule_122710 as lcm
import utilityModule as um
import sys, time,csv
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
################################################################################
#########################################################################################
def getAsitePosition_genomicHit(hit, gene, exptType):
    if gene.sense() == '+':
        interval = um.Locus(gene.chr(), gene.cdLocus().end(), hit.start(), gene.sense())
        if hit.start() >= gene.cdLocus().end():
            hitStartPosition = spliceInterval(interval, gene)
        else:
            hitStartPosition = -(spliceInterval(interval, gene))
    else:
        interval = um.Locus(gene.chr(), hit.end(), gene.cdLocus().start(), gene.sense())
        if hit.end() <= gene.cdLocus().start():
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
        interval = um.Locus(gene.chr(), gene.cdLocus().end(), refseqHit.genomeStart(), gene.sense())
        if refseqHit.genomeStart() >= gene.cdLocus().end():
            hitStartPosition = spliceInterval(interval, gene)
        else:
            hitStartPosition = -(spliceInterval(interval, gene))
    else:
        interval = um.Locus(gene.chr(), refseqHit.genomeStart(), gene.cdLocus().start(), gene.sense())
        if refseqHit.genomeStart() <= gene.cdLocus().start():
            hitStartPosition = spliceInterval(interval, gene)
        else:
            hitStartPosition = -(spliceInterval(interval, gene))
    if exptType == 'RPF':
        shiftedHitStartPosition = hitStartPosition + 15
    elif exptType == 'RNA':
        shiftedHitStartPosition = hitStartPosition
    return shiftedHitStartPosition
################################################################################
##################################################################################
def loadAccNumToLengths(): 
    #infile = open('/lab/bartel4_ata/huili/122610/mouse_querySeq/filteredGenesDetails_mouse_122610.txt')
    infile = open('/home/hguo/Documents/processing/filteredGenesDetails_human_240118.txt')
    infileLines = infile.readlines()
    infile.close()

    accNumToLengths = {}
    for line in infileLines[1:]:
        line = line.rstrip('\n')
        cols = line.split('\t')
        accNum = cols[1]
        ORFlength = int(cols[4])
        fpUTR = int(cols[5])
        tpUTR = int(cols[6])
        introns = int(cols[7])
        accNumToLengths[accNum] = [ORFlength, introns, fpUTR, tpUTR]
    print 'Number of genes loaded in accNumToLengths: ', len(accNumToLengths)
    return accNumToLengths
###############################################################################
###############################################################################
def loadAccNumList(infileString):
    infile = open(infileString)
    infileLines = infile.readlines()
    infile.close()

    accNumList = []
    for line in infileLines[1:]:
        line = line.rstrip('\n')
        accNumList.append(line)
    print 'Number of accNums in accNumList: ', len(accNumList)
    return accNumList
################################################################################
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
transfection = rootdirectory[:-1].split('/')[-1]
date = rootdirectory[:-1].split('/')[-2]
#exptType = transfection.split('_')[-1]
#mir = sys.argv[2]
#footprint_threshold = 10#int(sys.argv[3])
#rnaseq_threshold = 10#int(sys.argv[4])
exptType = sys.argv[2]
barcode = rootdirectory[:-1].split('/')[-1]
out_dir = sys.argv[3]
in_dir = rootdirectory[:-1].split('/')[-2]
print out_dir + 'AUGbaseToReads_' + in_dir +'_' + barcode + '.txt'

##accNumFile = sys.argv[5]
##accNumRemark = accNumFile.split('.')[0]

#################################################################################
## this part is only used for totalMappedReads
##if transfection.find('s_4') or transfection.find('s_6')!= -1:
##    mockOrMir = 'WT'
##elif transfection.find('s_5') or transfection.find('s_7') != -1:
##    mockOrMir = 'KO'

# if transfection.find('s_1') != -1 or transfection.find('WT') != -1:
#     mockOrMir = 'WT'
# elif transfection.find('s_2') != -1 or transfection.find('KO') != -1:
#     mockOrMir = 'KO'

# inDirectory = '/lab/bartel4_ata/huili/122610/mir155mouse/densityGroups/'
# print '======== load footprint genes that pass threshold ...', footprint_threshold, '================='
# accNumToRawReads_footprint, totalMappedReads_KO_footprint, totalMappedReads_WT_footprint = loadAccNumToRawReads(inDirectory+mir+'_RPF_filteredRawReads.txt', footprint_threshold)
# print '======== load RNAseq genes that pass threshold ...', rnaseq_threshold, '================='
# accNumToRawReads_rnaseq, totalMappedReads_KO_rnaseq, totalMappedReads_WT_rnaseq = loadAccNumToRawReads(inDirectory+mir+'_RNA_filteredRawReads.txt', rnaseq_threshold)

# commonGeneToRawReads = {}
# for gene, rawReads in accNumToRawReads_footprint.iteritems():
#     if gene in accNumToRawReads_rnaseq:
#         ## this gene is common, passes both thresholds
#         commonGeneToRawReads[gene] = (rawReads[0], rawReads[1], accNumToRawReads_rnaseq[gene][0], accNumToRawReads_rnaseq[gene][1])
        
# if exptType == 'RPF' and mockOrMir == 'WT':
#     totalMappedReads = totalMappedReads_WT_footprint
# elif exptType == 'RNA' and mockOrMir == 'WT':
#     totalMappedReads = totalMappedReads_WT_rnaseq
# elif exptType == 'RPF' and mockOrMir == 'KO':
#     totalMappedReads = totalMappedReads_KO_footprint
# elif exptType == 'RNA' and mockOrMir == 'KO':
#     totalMappedReads = totalMappedReads_KO_rnaseq
###################################################################################
accNumToLengths = loadAccNumToLengths()
#accNumList = loadAccNumList(accNumFile)

########################################################################
##########################################################################
# refPath = '/lab/bartel4_ata/huili/refFlats/mouse_refFlats/refFlat_122610.txt'
# splicedGeneORFpath = '/lab/bartel4_ata/huili/refFlats/mouse_refFlats/splicedGenes_ORF_122610.txt'
# chrList = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', 'X', 'Y', 'M']
# filteredGeneList = lcm.getFilteredGeneList_mouse(refPath, splicedGeneORFpath, chrList)
refPath = '/home/hguo/Documents/annotations/hg19/refFlat/refFlat_240118.txt'
splicedGeneORFpath = '/home/hguo/Documents/annotations/hg19/refFlat/splicedGenes_ORF_240118.txt'
chrList = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'M']
filteredGeneList = lcm.getFilteredGeneList_human(refPath, splicedGeneORFpath, chrList)
##############################################################################

#genomeDirectoryCore = (rootdirectory[:-1].split('/')[-1])+'_genome/'

#refseqDirectoryCore = (rootdirectory[:-1].split('/')[-1])+'_refseq/'
genomeDirectory = rootdirectory#+genomeDirectoryCore
refseqDirectory = rootdirectory+barcode+'_seed29_refseq/'


totalMappedReads = 0
rpkmfile = rootdirectory+'../shift15_geneTXCD_RPKMoutput_'+rootdirectory[:-1].split('/')[-1]+'_shift15.txt'
with open(rpkmfile) as fin:
    headerline = fin.next()
    for row in csv.reader(fin):
        str1 = ''.join(row)
        #print str1
        eles = str1.split('\t')
        txR = eles[1]
        totalMappedReads += float(txR)

print totalMappedReads


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
print 'uniqRefseqHits len: ' ,len(uniqRefseqHits)

#############################################################################
##################################################################################
bgName = sys.argv[4] # 'background gene list.txt'
f = open(bgName,'r')
bgFile = f.readlines()
f.close()
bgFile = [x.strip() for x in bgFile]
#selectedGenes = str(bgFile[0]).split(',')
print 'selectedGenes list len: ' , len(bgFile)
#selectedGenes = ['NM_001256799']

selectGeneAccDict = {}
for gene in filteredGeneList:
    if gene._accNum in bgFile:
        selectGeneAccDict[gene._accNum] = 1
print 'Select Genes: ' , len(selectGeneAccDict)


## initialize baseToReads dict
got = 0
## initialize codonToReads dict
baseToReads = {}
codonToReads = {}
for i in xrange(-102, 303, 1):
    baseToReads[i] = 0
    codon = i/3
    codonToReads[codon] = 0
print 'Number of bases in baseToReads: ', len(baseToReads)
print 'Number of codons in codonToReads: ', len(codonToReads)
##############################################################################
histoneAccNums = ['NM_206882', 'NM_030082', 'NM_175652', 'NM_178218', 'NM_178216', 'NM_054045', 'NM_175666', 'NM_178214', 'NM_175659', 'NM_178182', 'NM_033596', 'NM_013549', 'NM_178212', 'NM_145713', 'NM_178194', 'NM_178195', 'NM_175663', 'NM_175664', 'NM_023422', 'NM_178199', 'NM_178200', 'NM_178197', 'NM_178198', 'NM_175665', 'NM_178202', 'NM_178201', 'NM_178196', 'NM_010377', 'NM_015787', 'NM_030609', 'NM_020034', 'NM_015786', 'NM_178203', 'NM_175653', 'NM_013550', 'NM_013548', 'NM_145073', 'NM_178204', 'NM_178205', 'NM_178206', 'NM_178207', 'NM_178215', 'NM_178213', 'NM_175662', 'NM_178185', 'NM_178184', 'NM_178183', 'NM_178187', 'NM_178188', 'NM_178186', 'NM_175661', 'NM_175658', 'NM_178189', 'NM_175660', 'NM_178208', 'NM_178193', 'NM_178192', 'NM_175655', 'NM_175654', 'NM_178211', 'NM_178210', 'NM_175656', 'NM_153173', 'NM_175657']

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
        ORFlength = accNumToLengths[accNum][0]
        if ORFlength >= 400 and accNum in selectGeneAccDict:
        #if accNum in commonGeneToRawReads and ORFlength >= 600:
        #if accNum in commonGeneToRawReads:
        #if ORFlength >= 600:
            numGenesContribute += 1

            if accNum in exonHitsCollection:
                exonHits = exonHitsCollection[accNum]
                for exonHit in exonHits:
                    shiftedHitStart = getAsitePosition_genomicHit(exonHit, gene, exptType)
                    if -102 <= shiftedHitStart < 303:
                        ## record, shifted position relative to cdStart
                        baseToReads[shiftedHitStart] += exonHit.readNum()
                        ## now figure out codon
                        codon = shiftedHitStart/3
                        codonToReads[codon] += exonHit.readNum()
                        
            if accNum in uniqRefseqHits:
                refseqHits = uniqRefseqHits[accNum]
                for refseqHit in refseqHits:
                    shiftedHitStart = getAsitePosition_refseqHit(refseqHit, gene, exptType)
                    if -102 <= shiftedHitStart < 303:
                        ## record, shifted position relative to cdStart
                        baseToReads[shiftedHitStart] += refseqHit.readNum()
                        ## now figure out codon
                        codon = shiftedHitStart/3
                        codonToReads[codon] += refseqHit.readNum()
                        if got == 0:
                            print 'reached refseq reads'
                            got = 1
##########################################################################################
print 'Number of genes that contributed: ', numGenesContribute
## this shld be equal to number of genes in accNumList
#############################################################################                      

## output baseToReads
## output codonToReads
in_dir = rootdirectory[:-1].split('/')[-2]

base_outfileString = out_dir + 'AUGbaseToReads_' + in_dir + '_' + barcode + '.txt'
base_outfile = open(base_outfileString, 'w')
base_linesWritten = 0
codon_outfileString = out_dir + 'AUGcodonToReads_' + in_dir + '_' + barcode + '.txt'
codon_outfile = open(codon_outfileString, 'w')
codon_linesWritten = 0





for base in xrange(-102, 303, 1):
    base_outList = [base, ((((baseToReads[base])*1000000)/totalMappedReads)/float(numGenesContribute))]
    base_strOutList = map(str, base_outList)
    base_outListString = '\t'.join(base_strOutList)
    base_outfile.write(base_outListString+'\n')
    base_linesWritten += 1

    if base%3 == 0:
        codon = base/3
        codon_outList = [codon, ((((codonToReads[codon])*1000000)/totalMappedReads)/float(numGenesContribute))]
        codon_strOutList = map(str, codon_outList)
        codon_outListString = '\t'.join(codon_strOutList)
        codon_outfile.write(codon_outListString+'\n')
        codon_linesWritten += 1
base_outfile.write(str(numGenesContribute)+'\n')
base_linesWritten += 1
base_outfile.close()
codon_outfile.write(str(numGenesContribute)+'\n')
codon_linesWritten += 1
codon_outfile.close()
print 'Nunber of lines written for baseToReads: ', base_linesWritten
print 'Number of lines written for codonToReads: ', codon_linesWritten
print '========= Files written =============='
