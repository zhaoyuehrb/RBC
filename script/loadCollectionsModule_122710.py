import utilityModule as um
import os
###############################################################################
def filterAA(hit):
     seq = hit.seq()
     if seq.find('AAAAAAAAAAAAAAAAAA') != -1:
          return False
     return True
############################################################################################################################################################
def spliceTXExons(gene):
    txExons = gene.txExons()

    sumExonLength = 0
    for exon in txExons:
        exonLength = exon.len()
        sumExonLength += exonLength
    return sumExonLength
#####################################################################################################
def spliceCDExons(gene):
    cdExons = gene.cdExons()

    sumExonLength = 0
    for exon in cdExons:
        exonLength = exon.len()
        sumExonLength += exonLength
    return sumExonLength
#######################################################################################
def possibleNMDsubstrate(gene):
    if gene.sense() == '+':
        lastExon = gene.txExons()[-1]
        if lastExon.start() <= gene.cdLocus().end() <= lastExon.end():
            # last base in last exon, ok
            return False
        else:
            # last base is not in last exon, see if in second last exon
            secondLastExon = gene.txExons()[-2]
            if secondLastExon.start() <= gene.cdLocus().end() <= secondLastExon.end():
                # last base is in second last exon, see if more than 50nt away
                if (secondLastExon.end() - gene.cdLocus().end()) > 50:
                    # more than 50nt and in second last exon, NMD substrate
                    return True
                else:
                    # less than 50nt away though in second last exon, safe
                    return False
            else:
                # not in second last and not in last exon, NMD substrate
                return True
    else:
        lastExon = gene.txExons()[0]
        if lastExon.start() <= gene.cdLocus().start() <= lastExon.end():
            # last base in last exon, ok
            return False
        else:
            # last base is not in last exon, see if in second last exon
            secondLastExon = gene.txExons()[1]
            if secondLastExon.start() <= gene.cdLocus().start() <= secondLastExon.end():
                # last base is in second last exon, see if more than 50nt away
                if (gene.cdLocus().start() - secondLastExon.start()) > 50:
                    # more than 50nt and in second last exon, NMD substrate
                    return True
                else:
                    # less than 50nt away though in second last exon, safe
                    return False
            else:
                # not in second last and not in last exon, NMD substrate
                return True
########################################################################################
def prematureStop(ORFseq):
    stopCodons = ['TGA', 'TAG', 'TAA']
    for i in range(0, len(ORFseq)-3, 3):
        codon = ORFseq[i:i+3]
        if codon in stopCodons:
            return True
    return False
#################################################################################################
def getFilteredGeneList_mouse(refPath, splicedGeneORFpath, chrList):
    chrList = map(lambda i: 'chr'+i, chrList)
    stopCodons = ['TGA', 'TAG', 'TAA']
    
    genesList = um.makeGenes(refPath)
    ORFseqs = um.get_fa(splicedGeneORFpath)
    print 'Number of entries in ORFseqs:', len(ORFseqs)

    genesDict = {}
    for gene in genesList:
        chrom = gene.chr()
        if chrom in chrList:
            if gene.accNum()[:2] == 'NM':
                if not possibleNMDsubstrate(gene):
                    ORF_length = spliceCDExons(gene)
                    if ORF_length%3 == 0:
                        cdLocusString = str(gene.cdLocus())
                        key = gene.accNum()+'_'+cdLocusString
                        ORFseq = ORFseqs[key]
                        if (ORFseq[-3:] in stopCodons):
                            if not prematureStop(ORFseq):
                                if (not genesDict.has_key(gene.name())):
                                    genesDict[gene.name()] = [gene]
                                else:
                                    genesDict[gene.name()].append(gene)
             
    print 'Length of genesDict: ', len(genesDict)

    filteredGeneList = []
    for key in genesDict.keys():
        if len(genesDict[key]) == 1:
            filteredGeneList.append(genesDict[key][0])
        else:
            lengthLoci = []
            for gene in genesDict[key]:
                mrnaLength = spliceTXExons(gene)
                #txLocus = gene.txLocus()
                #length = txLocus.len()
                lengthLoci.append(mrnaLength)
            maxLength = max(lengthLoci)
            index_max = lengthLoci.index(maxLength)
            filteredGeneList.append(genesDict[key][index_max])
    print 'Length of filtered gene list: ', len(filteredGeneList)
    return filteredGeneList
##################################################################################################
def getFilteredGeneList_human(refPath, splicedGeneORFpath, chrList):
    chrList = map(lambda i: 'chr'+i, chrList)
    stopCodons = ['TGA', 'TAG', 'TAA']
    
    genesList = um.makeGenes(refPath)
    ORFseqs = um.get_fa(splicedGeneORFpath)
    print 'Number of entries in ORFseqs:', len(ORFseqs)

    genesDict = {}
    for gene in genesList:
        chrom = gene.chr()
        if chrom in chrList:
            if gene.accNum()[:2] == 'NM':# and (gene.accNum() != 'NM_057179'):
                if not possibleNMDsubstrate(gene):
                    ORF_length = spliceCDExons(gene)
                    if ORF_length%3 == 0:
                        cdLocusString = str(gene.cdLocus())
                        key = gene.accNum()+'_'+cdLocusString
                        ORFseq = ORFseqs[key]
                        if (ORFseq[-3:] in stopCodons):
                            if not prematureStop(ORFseq):
                                if (not genesDict.has_key(gene.name())):
                                    genesDict[gene.name()] = [gene]
                                else:
                                    genesDict[gene.name()].append(gene)
             
    print 'Length of genesDict: ', len(genesDict)

    filteredGeneList = []
    for key in genesDict.keys():
        if len(genesDict[key]) == 1:
            filteredGeneList.append(genesDict[key][0])
        else:
            lengthLoci = []
            for gene in genesDict[key]:
                mrnaLength = spliceTXExons(gene)
                #txLocus = gene.txLocus()
                #length = txLocus.len()
                lengthLoci.append(mrnaLength)
            maxLength = max(lengthLoci)
            index_max = lengthLoci.index(maxLength)
            filteredGeneList.append(genesDict[key][index_max])
    print 'Length of filtered gene list: ', len(filteredGeneList)
    return filteredGeneList
##################################################################################################
##################################################################################################
def linefilter(splitline, chrom):
    chromlabel = splitline[1]
    rrnalabel = splitline[10]
    ncrnalabel = splitline[11]
    return chromlabel == chrom and rrnalabel == 'NA' and ncrnalabel == 'NA'
##############################################################################################
def getUniqGenomicHitsByChrom(chrom, directory, lenRange):
    lenRange = map(int, lenRange.split(','))
    shortest, longest = str(lenRange[0]), str(lenRange[1] - 1)
    datasets = [directory]
    print 'Dataset is from: ', directory
    print '======== Loading hits from '+chrom+' only ...... ========='
    print '======== EXCLUDING rRNA and other RNA! ==================='
    requirementTrue = lambda x: True
    genomic_hits = um.getAllGenomicHits(datasets, lenRange, lambda line: linefilter(line, chrom), requirementTrue)
    print 'Number of such hits in list for '+shortest+'-'+longest+' nt: ', len(genomic_hits)

    ## filter the list for AA reads and get hitsCollection
    filteredList = filter(filterAA, genomic_hits)
    print
    print '=============== AFTER FILTERING 18xA ===================================='
    print 'Number of hits in filtered list for '+shortest+'-'+longest+' nt: ', len(filteredList)
    uniqGenomicHits = um.LocusCollection(filteredList, 1000)
    print 'Number of hits in collection for '+shortest+'-'+longest+' nt: ', len(uniqGenomicHits)
    return uniqGenomicHits
####################################################################################################
def getTXhitsCollection(uniqHitsCollection, fgList):
    ## genomicHits is a locusCollection, all hits are UNIQUE already
    gene_hitsDict = {}   # Dict A
    gene_repeatHitsDict = {}     # Dict B

    for gene in fgList:
        for txExon in gene.txExons():
            hitsInExon = uniqHitsCollection.getContained(txExon, sense = 'sense')
            for hit in hitsInExon:
                if (not gene_hitsDict.has_key(hit)):
                    gene_hitsDict[hit] = None
                else:
                    # this hit was grabbed before, put in Dict B
                    if (not gene_repeatHitsDict.has_key(hit)):
                        gene_repeatHitsDict[hit] = None
                    else:
                        # dict B already has the key, this hit is grabbed at least the third time
                        # but it's ok becos already listed as a repeat
                        pass
    unambiExonHits = {}    # Dict C
    for hit in gene_hitsDict:
        if hit not in gene_repeatHitsDict:
            # grabbed once only, add to dict C
            unambiExonHits[hit] = None
        else:
            # this hit is in Dict B, ie it was grabbed more than once
            # don't add to dict C
            pass
    ## Dict C is now filled with hits from unambiguous TX exons
    txHitsCollection = um.LocusCollection(unambiExonHits.keys(), 1000)
    print 'Number of hits in txHitsCollection: ', len(txHitsCollection)
    return txHitsCollection
#######################################################################################################
def getTXhitsCollectionOverlap(uniqHitsCollection, fgList):
    ## genomicHits is a locusCollection, all hits are UNIQUE already
    gene_hitsDict = {}   # Dict A
    gene_repeatHitsDict = {}     # Dict B

    for gene in fgList:
        for txExon in gene.txExons():
            hitsInExon = uniqHitsCollection.getOverlap(txExon, sense = 'sense')
            for hit in hitsInExon:
                if (not gene_hitsDict.has_key(hit)):
                    gene_hitsDict[hit] = None
                else:
                    # this hit was grabbed before, put in Dict B
                    if (not gene_repeatHitsDict.has_key(hit)):
                        gene_repeatHitsDict[hit] = None
                    else:
                        # dict B already has the key, this hit is grabbed at least the third time
                        # but it's ok becos already listed as a repeat
                        pass
    unambiExonHits = {}    # Dict C
    for hit in gene_hitsDict:
        if hit not in gene_repeatHitsDict:
            # grabbed once only, add to dict C
            unambiExonHits[hit] = None
        else:
            # this hit is in Dict B, ie it was grabbed more than once
            # don't add to dict C
            pass
    ## Dict C is now filled with hits from unambiguous TX exons
    txHitsCollection = um.LocusCollection(unambiExonHits.keys(), 1000)
    print 'Number of hits in txHitsCollection: ', len(txHitsCollection)
    return txHitsCollection
#######################################################################################################
##################################################################################################
def minimalfilter(splitline):
    chromlabel = splitline[1]
    rrnalabel = splitline[10]
    ncrnalabel = splitline[11]
    return rrnalabel == 'NA' and ncrnalabel == 'NA'
##########################################################################################################
def getAllUniqGenomicHits(directory, lenRange):
    lenRange = map(int, lenRange.split(','))
    shortest, longest = str(lenRange[0]), str(lenRange[1] - 1)
    datasets = [directory]
    print 'Dataset is from: ', directory
    print '======== Loading ALL QUALIFIED hits ========='
    print '======== EXCLUDING rRNA and other RNA! ==================='
    requirementTrue = lambda x: True
    genomic_hits = um.getAllGenomicHits(datasets, lenRange, lambda line: minimalfilter(line), requirementTrue)
    print 'Number of such hits in list for '+shortest+'-'+longest+' nt: ', len(genomic_hits)

    ## filter the list for AA reads and get hitsCollection
    filteredList = filter(filterAA, genomic_hits)
    print
    print '=============== AFTER FILTERING 18xA ===================================='
    print 'Number of hits in filtered list for '+shortest+'-'+longest+' nt: ', len(filteredList)
    uniqGenomicHits = um.LocusCollection(filteredList, 1000)
    print 'Number of hits in collection for '+shortest+'-'+longest+' nt: ', len(uniqGenomicHits)
    return uniqGenomicHits
###################################################################################################
######################################################################################################
class GeneToRank:
    def __init__(self, name, coverage, readNum, numHits, fgInstance):
        self.name = name
        self.coverage = coverage
        self.readNum = readNum
        self.numHits = numHits
        self.fgInstance = fgInstance
    def getName(self):
        return self.name
    def getCoverage(self):
        return self.coverage
    def getReadNum(self):
        return self.readNum
    def getNumHits(self):
        return self.numHits
    def getFGinstance(self):
        return self.fgInstance
    def __cmp__(self, other):
        if self.coverage > other.coverage:
            return -1
        elif self.coverage < other.coverage:
            return 1
        else:
            if self.readNum > other.readNum:
                return -1
            elif self.readNum < other.readNum:
                return 1
            else:
                return 0
##############################################################################################################
def getGenomeQuantifiedGenes(txHits, fgList):
    quantifiedGenes = []
    totalExonReadNum = 0
    for gene in fgList:
        splicedLength = spliceTXExons(gene)
        readNumExons = 0
        numHitsExons = 0
        for txExon in gene.txExons():
            hitsInExon = txHits.getOverlap(txExon, sense = 'sense')
            numHitsExons += len(hitsInExon)
            for hit in hitsInExon:
                readNumExons += hit.readNum()
        mRNACoverage = numHitsExons/float(splicedLength)     #### this may be changed
        if readNumExons != 0: ## has at least one read to this gene's txExons
            geneToRank = GeneToRank(gene.name(), mRNACoverage, readNumExons, numHitsExons, gene)
            quantifiedGenes.append(geneToRank)
        totalExonReadNum += readNumExons
    print 'Number of quantified genes: ', len(quantifiedGenes)
    print 'TotalExonReadNum: ', totalExonReadNum
    return quantifiedGenes
##################################################################################
def convertReadsToRPKM(gene, readNum, totalTXexonreadNum):
    sumLengthExons = spliceTXExons(gene)
    rpkm = (1000000000*readNum)/float(totalTXexonreadNum*sumLengthExons)
    return rpkm
##################################################################################################
def getRPKMdict(quantifiedGenes, totalTXreadNum):
     rpkmDict = {}
     for qGene in quantifiedGenes:
        fgInstance = qGene.getFGinstance()
        readNum = qGene.getReadNum()
        rpkm = convertReadsToRPKM(fgInstance, readNum, totalTXreadNum)
        rpkmDict[fgInstance.accNum()] = (readNum, rpkm)
     print 'Number of genes in rpkmDict: ', len(rpkmDict)
     return rpkmDict
##################################################################################
#################################################################################
class RefseqHit:
    def __init__(self, entryID, chrom, accNum, genomeStart, length, sense, blockStarts, blockEnds, seq):
        self.entryID = entryID
        self.chrom = chrom
        self._accNum = accNum
        self._genomeStart = int(genomeStart)
        self.length = int(length)
        self._sense = sense
        self._blockStarts = blockStarts
        self._blockEnds = blockEnds
        self._seq = seq

    def getEntryID(self):
        return self.entryID
    def chr(self):
        return self.chrom
    def accNum(self):
        return self._accNum
    def genomeStart(self):
        return self._genomeStart
    def len(self):
        return self.length
    def sense(self):
        return self._sense
    def blockStarts(self):
        return self._blockStarts
    def blockEnds(self):
        return self._blockEnds

    def readNum(self):
        readNum = int(self.entryID.split('_')[1])
        return readNum
    def seq(self):
        return self._seq

    def firstBlockStart(self):
        return int(self.blockStarts().split(',')[0])
    def secondBlockStart(self):
        return int(self.blockStarts().split(',')[1])
    def firstBlockEnd(self):
        return int(self.blockEnds().split(',')[0])
    def secondBlockEnd(self):
        return int(self.blockEnds().split(',')[1])

    def firstBlockCoords(self):
        firstBlockCoords = [self.firstBlockStart(), self.firstBlockEnd()]
        firstBlockCoords.sort()
        return firstBlockCoords
    def secondBlockCoords(self):
        secondBlockCoords = [self.secondBlockStart(), self.secondBlockEnd()]
        secondBlockCoords.sort()
        return secondBlockCoords       

    def firstBlockSize(self):
        firstBlockCoords = self.firstBlockCoords()
        return firstBlockCoords[1] - firstBlockCoords[0] + 1
    def secondBlockSize(self):
        secondBlockCoords = self.secondBlockCoords()
        return secondBlockCoords[1] - secondBlockCoords[0] + 1

    def genomeEnd(self):
        return int(self.blockEnds().split(',')[1])

    def span(self):
        actualCoords = [self.genomeStart(), self.genomeEnd()]
        actualCoords.sort()
        return actualCoords[1] - actualCoords[0] + 1        

    def __str__(self):
        firstBlockString = '-'.join(map(lambda i: str(i), self.firstBlockCoords()))
        secondBlockString = '-'.join(map(lambda i: str(i), self.secondBlockCoords()))
        return self.chr()+'('+self.sense()+'):'+firstBlockString+'_'+secondBlockString
        
##################################################################################
def getAllUniqRefseqHits(directory, lengthTuple):
    corefilename = directory[:-1].split('/')[-1]

    lenTuple = map(int, lengthTuple.split(','))
    lengths = []
    for i in range(lenTuple[0], lenTuple[1]):
        lengths.append(i)
    dirfiles = os.listdir(directory)
    avaiLengths = []
    for dirfile in dirfiles:
        if dirfile[-3:] == '.bl':
            if int(dirfile.split('.')[-2]) in lengths:
                avaiLengths.append(int(dirfile.split('.')[-2]))
##    shortest = int(lengthTuple.split(',')[0])
##    longest = int(lengthTuple.split(',')[1])
##    lengths = []
##    #print shortest, longest
##    #excludedLengths = [31, 32, 33, 34, 35]
##    excludedLengths = [19, 31, 32, 33, 34, 35]
##    for i in range(shortest, longest, 1):
##        if i not in excludedLengths:
##            lengths.append(i)
##    print lengths

    geneToRefseqHits = {}
    seqToRefseqHits  = {}
    hitsDict = {}

    for length in avaiLengths:
        length = str(length)
        infileString = directory+corefilename+'.'+length+'.bl'
        infile = open(infileString)
        infileLines = infile.readlines()
        infile.close()

        #print 'Number of bl entries: ', len(infileLines[:-1])
        faEntries = um.get_fa(directory+corefilename+'.'+length+'.fa')  
        #print 'Number of faEntries: ', len(faEntries)

        for line in infileLines[:-1]:
            line = line.rstrip('\n')
            cols = line.split('\t')
            entryID = cols[0]
            chrom = cols[1]
            accNum = cols[2]
            readlength = int(cols[3])
            genomeStart = int(cols[8])
            surrogateEnd = int(cols[9])
            blockStarts = cols[10]
            blockEnds = cols[11]
            if genomeStart < surrogateEnd:
                sense = '+'
            else:
                sense = '-'
            seq = faEntries[entryID]
            refseqHitInstance = RefseqHit(entryID, chrom, accNum, genomeStart, readlength, sense, blockStarts, blockEnds, seq)
            if not hitsDict.has_key(refseqHitInstance):
                hitsDict[refseqHitInstance] = None
                
##            if not geneToRefseqHits.has_key(accNum):
##                geneToRefseqHits[accNum] = [refseqHitInstance]
##            else:
##                geneToRefseqHits[accNum].append(refseqHitInstance)
            if not seqToRefseqHits.has_key(seq):
                seqToRefseqHits[seq] = [refseqHitInstance]
            else:
                seqToRefseqHits[seq].append(refseqHitInstance)

    print directory, len(seqToRefseqHits)
    print 'Number of hits in refseq hits list: ', len(hitsDict) 
    filteredRefseqHits = filter(filterAA, hitsDict.keys())
    print '=============== AFTER FILTERING 18xA, refseq hits ============================='
    print 'Number of hits in filtered refseq hits list: ', len(filteredRefseqHits)

    for hit in filteredRefseqHits:
        accNum = hit.accNum()
        if not geneToRefseqHits.has_key(accNum):
            geneToRefseqHits[accNum] = [hit]
        else:
            geneToRefseqHits[accNum].append(hit)            
			
    print 'Number of genes quantified by refseq hits: ', len(geneToRefseqHits)
    return geneToRefseqHits
################################################################################
def getQuantifiedGenes(txHits, geneToRefseqHits, fgList):
    quantifiedGenes = []
    totalExonReadNum = 0
    for gene in fgList:
        splicedLength = spliceTXExons(gene)
        readNumExons = 0
        numHitsExons = 0
        for txExon in gene.txExons():
            hitsInExon = txHits.getOverlap(txExon, sense = 'sense')
            numHitsExons += len(hitsInExon)
            for hit in hitsInExon:
                readNumExons += hit.readNum()
        ## at this point, readNumExons and numHitsExons filled with genome hits
        ## now fill refseq hits
        accNum = gene.accNum()
        if accNum in geneToRefseqHits:
            hitsInGene = geneToRefseqHits[accNum]
            numHitsExons += len(hitsInGene)
            for hit in hitsInGene:
                readNumExons += hit.readNum()
        ## now shld be filled with both genome hits and refseq hits
        mRNACoverage = numHitsExons/float(splicedLength)
        if readNumExons != 0: ## has at least one read to this gene's txExons
            geneToRank = GeneToRank(gene.name(), mRNACoverage, readNumExons, numHitsExons, gene)
            quantifiedGenes.append(geneToRank)
        totalExonReadNum += readNumExons
    print 'Number of quantified genes: ', len(quantifiedGenes)
    print 'TotalExonReadNum: ', totalExonReadNum
    return quantifiedGenes
############################################################################
def getUniqRefseqHitsByChrom(directory, lengthTuple, chromosome):
    corefilename = directory[:-1].split('/')[-1]

    lenTuple = map(int, lengthTuple.split(','))
    lengths = []
    for i in range(lenTuple[0], lenTuple[1]):
        lengths.append(i)
    dirfiles = os.listdir(directory)
    avaiLengths = []
    for dirfile in dirfiles:
        if dirfile[-3:] == '.bl':
            if int(dirfile.split('.')[-2]) in lengths:
                avaiLengths.append(int(dirfile.split('.')[-2]))
##    shortest = int(lengthTuple.split(',')[0])
##    longest = int(lengthTuple.split(',')[1])
##    lengths = []
##    #print shortest, longest
##    #excludedLengths = [31, 32, 33, 34, 35]
##    excludedLengths = [19, 31, 32, 33, 34, 35]
##    for i in range(shortest, longest, 1):
##        if i not in excludedLengths:
##            lengths.append(i)
##    print lengths

    geneToRefseqHits = {}
    seqToRefseqHits  = {}
    hitsDict = {}

    for length in avaiLengths:
        length = str(length)
        infileString = directory+corefilename+'.'+length+'.bl'
        infile = open(infileString)
        infileLines = infile.readlines()
        infile.close()

        #print 'Number of bl entries: ', len(infileLines[:-1])
        faEntries = um.get_fa(directory+corefilename+'.'+length+'.fa')  
        #print 'Number of faEntries: ', len(faEntries)

        for line in infileLines[:-1]:
            line = line.rstrip('\n')
            cols = line.split('\t')
            entryID = cols[0]
            chrom = cols[1]
            accNum = cols[2]
            readlength = int(cols[3])
            genomeStart = int(cols[8])
            surrogateEnd = int(cols[9])
            blockStarts = cols[10]
            blockEnds = cols[11]
            if genomeStart < surrogateEnd:
                sense = '+'
            else:
                sense = '-'
            seq = faEntries[entryID]
            if chrom == chromosome:
                refseqHitInstance = RefseqHit(entryID, chrom, accNum, genomeStart, readlength, sense, blockStarts, blockEnds, seq)
                if not hitsDict.has_key(refseqHitInstance):
                    hitsDict[refseqHitInstance] = None
                
##            if not geneToRefseqHits.has_key(accNum):
##                geneToRefseqHits[accNum] = [refseqHitInstance]
##            else:
##                geneToRefseqHits[accNum].append(refseqHitInstance)
                if not seqToRefseqHits.has_key(seq):
                    seqToRefseqHits[seq] = [refseqHitInstance]
                else:
                    seqToRefseqHits[seq].append(refseqHitInstance)

    print directory, len(seqToRefseqHits)
    print 'Number of hits in refseq hits list for ',chromosome,' : ', len(hitsDict) 
    filteredRefseqHits = filter(filterAA, hitsDict.keys())
    print '=============== AFTER FILTERING 18xA, refseq hits ============================='
    print 'Number of hits in filtered refseq hits list ',chromosome,' : ', len(filteredRefseqHits)

    for hit in filteredRefseqHits:
        accNum = hit.accNum()
        if not geneToRefseqHits.has_key(accNum):
            geneToRefseqHits[accNum] = [hit]
        else:
            geneToRefseqHits[accNum].append(hit)            
			
    print 'Number of genes on quantified by refseq hits, on ',chromosome,' : ', len(geneToRefseqHits)
    return geneToRefseqHits
#######################################################################################
##
## added on 6/25/09, intronHitsCollection
#######################################################################################################
def getIntronHitsCollectionOverlap(uniqHitsCollection, fgList):
    ## genomicHits is a locusCollection, all hits are UNIQUE already
    gene_hitsDict = {}   # Dict A
    gene_repeatHitsDict = {}     # Dict B

    for gene in fgList:
        for intron in gene.introns():
            hitsInIntron = uniqHitsCollection.getOverlap(intron, sense = 'sense')
            for hit in hitsInIntron:
                if (not gene_hitsDict.has_key(hit)):
                    gene_hitsDict[hit] = None
                else:
                    # this hit was grabbed before, put in Dict B
                    if (not gene_repeatHitsDict.has_key(hit)):
                        gene_repeatHitsDict[hit] = None
                    else:
                        # dict B already has the key, this hit is grabbed at least the third time
                        # but it's ok becos already listed as a repeat
                        pass
    unambiIntronHits = {}    # Dict C
    for hit in gene_hitsDict:
        if hit not in gene_repeatHitsDict:
            # grabbed once only, add to dict C
            unambiIntronHits[hit] = None
        else:
            # this hit is in Dict B, ie it was grabbed more than once
            # don't add to dict C
            pass
    ## Dict C is now filled with hits from unambiguous TX exons
    intronHitsCollection = um.LocusCollection(unambiIntronHits.keys(), 1000)
    print 'Number of hits in intronHitsCollection: ', len(intronHitsCollection)
    return intronHitsCollection
#######################################################################################################
def passBound_single(exon, hit):
    if exon.sense() == '+':
        lowerbd = 0 ## zero-based start of exon
        upperbdRaw = exon.end() - 14
        upperbd = upperbdRaw - exon.start() ## make upperbd relative to exonStart
        hitStart = hit.start() - exon.start()
        #print hit.seq()
        #print 'Raw hit start: ', hit.start()
        #print 'Relative hit start: ', hitStart

    else:
        lowerbd = 0 ## zero-based start of exon
        upperbdRaw = exon.start() + 14
        upperbd = exon.end() - upperbdRaw ## make upperbd relative to exonStart
        hitStart = exon.end() - hit.end()
        #print hit.seq()
        #print 'Raw hit start: ', hit.end()
        #print 'Relative hit start: ', hitStart
    if lowerbd <= hitStart <= upperbd:
        return True
    else:
        return False
################################################################################
def passBound_multiFirst(exon, hit):
    if exon.sense() == '+':
        lowerbd = 0 ## zero-based start of exon
        upperbdRaw = exon.end() - 18
        upperbd = upperbdRaw - exon.start() ## make upperbd relative to exonStart
        hitStart = hit.start() - exon.start()
        #print hit.seq()
        #print 'Raw hit start: ', hit.start()
        #print 'Relative hit start: ', hitStart

    else:
        lowerbd = 0 ## zero-based start of exon
        upperbdRaw = exon.start() + 18
        upperbd = exon.end() - upperbdRaw ## make upperbd relative to exonStart
        hitStart = exon.end() - hit.end()
        #print hit.seq()
        #print 'Raw hit start: ', hit.end()
        #print 'Relative hit start: ', hitStart
    if lowerbd <= hitStart <= upperbd:
        return True
    else:
        return False
##################################################################################
def passBound_multiLast(exon, hit):
    if exon.sense() == '+':
        lowerbd = -8 ## zero-based start of exon
        upperbdRaw = exon.end() - 14
        upperbd = upperbdRaw - exon.start() ## make upperbd relative to exonStart
        hitStart = hit.start() - exon.start()
        #print hit.seq()
        #print 'Raw hit start: ', hit.start()
        #print 'Relative hit start: ', hitStart

    else:
        lowerbd = -8 ## zero-based start of exon
        upperbdRaw = exon.start() + 14
        upperbd = exon.end() - upperbdRaw ## make upperbd relative to exonStart
        hitStart = exon.end() - hit.end()
        #print hit.seq()
        #print 'Raw hit start: ', hit.end()
        #print 'Relative hit start: ', hitStart
    if lowerbd <= hitStart <= upperbd:
        return True
    else:
        return False
################################################################################
def passBound_multiMiddle(exon, hit):
    if exon.sense() == '+':
        lowerbd = -8 ## zero-based start of exon
        upperbdRaw = exon.end() - 18
        upperbd = upperbdRaw - exon.start() ## make upperbd relative to exonStart
        hitStart = hit.start() - exon.start()
        #print hit.seq()
        #print 'Raw hit start: ', hit.start()
        #print 'Relative hit start: ', hitStart

    else:
        lowerbd = -8 ## zero-based start of exon
        upperbdRaw = exon.start() + 18
        upperbd = exon.end() - upperbdRaw ## make upperbd relative to exonStart
        hitStart = exon.end() - hit.end()
        #print hit.seq()
        #print 'Raw hit start: ', hit.end()
        #print 'Relative hit start: ', hitStart
    if lowerbd <= hitStart <= upperbd:
        return True
    else:
        return False
##################################################################################
def passBound_intron(intron, hit):
    if intron.sense() == '+':
        lowerbd = 0 ## zero-based start of intron
        upperbdRaw = intron.end() - 8
        upperbd = upperbdRaw - intron.start() ## make upperbd relative to intronStart
        hitStart = hit.start() - intron.start()
        #print hit.seq()
        #print 'Raw hit start: ', hit.start()
        #print 'Relative hit start: ', hitStart

    else:
        lowerbd = 0 ## zero-based start of intron
        upperbdRaw = intron.start() + 8
        upperbd = intron.end() - upperbdRaw ## make upperbd relative to intronStart
        hitStart = intron.end() - hit.end()
        #print hit.seq()
        #print 'Raw hit start: ', hit.end()
        #print 'Relative hit start: ', hitStart
    if lowerbd <= hitStart <= upperbd:
        return True
    else:
        return False
###########################################################################################
def getUniqTXExonIntronHits(uniqHitsCollection, fgList):
    gene_hitsDict = {}   # Dict A
    gene_repeatHitsDict = {}     # Dict B

    intron_hitsDict = {}   # Dict X
    intron_repeatHitsDict = {} # Dict Y

    for gene in fgList:
        numExons = len(gene.txExons())
        if numExons == 1:
            exon = gene.txExons()[0]
            hitsOverlapExon = uniqHitsCollection.getOverlap(exon, sense = 'sense')
            for hit in hitsOverlapExon:
                if passBound_single(exon, hit):
                    if (not gene_hitsDict.has_key(hit)):
                        gene_hitsDict[hit] = gene.accNum()
                    else:
                        ## already in dictA, this is a repeat
                        ## put in dictB
                        if (not gene_repeatHitsDict.has_key(hit)):
                            gene_repeatHitsDict[hit] = gene.accNum()
                        else:
                            ## already in dictB
                            ## at least grabbed the third time
                            pass
                else:
                    ## hit does not start within acceptable boundaries
                    ## even though overlap this exon
                    pass
               
        elif numExons > 1:
            if gene.sense() == '+':
                firstExon = gene.txExons()[0]
                lastExon = gene.txExons()[-1]
            else:
                firstExon = gene.txExons()[-1]
                lastExon = gene.txExons()[0]

            ## first exon #################
            hitsOverlapFirstExon = uniqHitsCollection.getOverlap(firstExon, sense = 'sense')
            for hit in hitsOverlapFirstExon:
                if passBound_multiFirst(firstExon, hit):
                    if (not gene_hitsDict.has_key(hit)):
                        gene_hitsDict[hit] = gene.accNum()
                    else:
                        ## already in dictA, this is a repeat
                        ## put in dictB
                        if (not gene_repeatHitsDict.has_key(hit)):
                            gene_repeatHitsDict[hit] = gene.accNum()
                        else:
                            ## already in dictB
                            ## at least grabbed the third time
                            pass
                else:
                    ## hit does not start within acceptable boundaries
                    ## even though overlap this exon
                    pass

            # last exon #######################
            hitsOverlapLastExon = uniqHitsCollection.getOverlap(lastExon, sense = 'sense')
            for hit in hitsOverlapLastExon:
                if passBound_multiLast(lastExon, hit):
                    if (not gene_hitsDict.has_key(hit)):
                        gene_hitsDict[hit] = gene.accNum()
                    else:
                        ## already in dictA, this is a repeat
                        ## put in dictB
                        if (not gene_repeatHitsDict.has_key(hit)):
                            gene_repeatHitsDict[hit] = gene.accNum()
                        else:
                            ## already in dictB
                            ## at least grabbed the third time
                            pass
                else:
                    ## hit does not start within acceptable boundaries
                    ## even though overlap this exon
                    pass

            ## middle exons ###################
            for exon in gene.txExons()[1:-1]:
                hitsOverlapExon = uniqHitsCollection.getOverlap(exon, sense = 'sense')
                for hit in hitsOverlapExon:
                    if passBound_multiMiddle(exon, hit):
                        if (not gene_hitsDict.has_key(hit)):
                            gene_hitsDict[hit] = gene.accNum()
                        else:
                            ## already in dictA, this is a repeat
                            ## put in dictB
                            if (not gene_repeatHitsDict.has_key(hit)):
                                gene_repeatHitsDict[hit] = gene.accNum()
                            else:
                                ## already in dictB
                                ## at least grabbed the third time
                                pass
                    else:
                        ## hit does not start within acceptable boundaries
                        ## even though overlap this exon
                        pass

            ## introns #######################
            for intron in gene.introns():
                hitsOverlapIntron = uniqHitsCollection.getOverlap(intron, sense = 'sense')
                for hit in hitsOverlapIntron:
                    if passBound_intron(intron, hit):
                        if (not intron_hitsDict.has_key(hit)):
                            intron_hitsDict[hit] = gene.accNum()
                        else:
                            ## already in dictX, this is a repeat
                            ## put in dictY
                            if (not intron_repeatHitsDict.has_key(hit)):
                                intron_repeatHitsDict[hit] = gene.accNum()
                            else:
                                ## already in dictY
                                ## at least grabbed the third time
                               pass
                    else:
                        ## hit does not start within acceptable boundaries
                        ## even though overlap this exon
                        pass

    print 'Number of hits in dictA: ', len(gene_hitsDict)
    print 'Number of hits in dictB: ', len(gene_repeatHitsDict)

    print 'Number of hits in dictX: ', len(intron_hitsDict)
    print 'Number of hits in dictY: ', len(intron_repeatHitsDict)

    ################# unambi exons ##########################################

    unambiExonHits = {} # Dict C
    for hit in gene_hitsDict:
        if hit not in gene_repeatHitsDict:
            # grabbed once only, add to dict C
            unambiExonHits[hit] = gene_hitsDict[hit]
        else:
            # this hit is in Dict B, ie it was grabbed more than once
            # don't add to dict C
            pass

    ## Dict C is now filled with hits from unambiguous TX exons
    print 'Number of hits in dictC, chr1: ', len(unambiExonHits)
    txHitsCollection = um.LocusCollection(unambiExonHits.keys(), 1000)
    print 'Number of hits in txHitsCollection: ', len(txHitsCollection)

    ########### unambi introns #############################################

    unambiIntronHits = {} # Dict Z
    for hit in intron_hitsDict:
        if hit not in intron_repeatHitsDict:
            # grabbed once only, add to dict Z
            unambiIntronHits[hit] = intron_hitsDict[hit]
        else:
            # this hit is in Dict Y, ie it was grabbed more than once
            # don't add to dict Z
            pass

    ## Dict Z is now filled with hits from unambiguous introns
    print 'Number of hits in dictZ, chr1: ', len(unambiIntronHits)
    intronHitsCollection = um.LocusCollection(unambiIntronHits.keys(), 1000)
    print 'Number of hits in intronHitsCollection: ', len(intronHitsCollection)

    ############### unambi overall ##############################################

    unambiguousExonHits = {}
    unambiguousIntronHits = {}

    for hit in unambiIntronHits:
        if hit not in unambiExonHits:
            ## this hit was not grabbed by any exon
            ## therefore unambiguously intronic
            unambiguousIntronHits[hit] = unambiIntronHits[hit]
        else:
            ## this hit was grabbed by an exon before
            ## do not add to unambiguousIntronHits dict
            pass

    print 'Number of unambiguously intronic hits, chr1: ', len(unambiguousIntronHits)

    for hit in unambiExonHits:
        if hit not in unambiIntronHits:
            ## this hit was not grabbed by any intron
            ## therefore unambiguously exonic
            unambiguousExonHits[hit] = unambiExonHits[hit]
        else:
            ## this hit was grabbed by an intron before
            ## do not add to unambiguousExonHits dict
            #print hit
            pass
    

    print 'Number of unambiguously exonic hits, chr1: ', len(unambiguousExonHits)

    ############ organise by genes #######################################
    geneToIntronHits = {}
    for intronHit, accNum in unambiguousIntronHits.iteritems():
        if (not geneToIntronHits.has_key(accNum)):
            geneToIntronHits[accNum] = [intronHit]
        else:
            geneToIntronHits[accNum].append(intronHit)
    print 'Number of genes in geneToIntronHits: ', len(geneToIntronHits)
    totalIntronHits = 0
    totalIntronReadNum = 0
    for gene in geneToIntronHits:
        totalIntronHits += len(geneToIntronHits[gene])
        for intronHit in geneToIntronHits[gene]:
           totalIntronReadNum += intronHit.readNum()
    print 'Total number of intron hits by geneToIntronHits: ', totalIntronHits
    print 'Total intron readNum: ', totalIntronReadNum

    geneToExonHits = {}
    for exonHit, accNum in unambiguousExonHits.iteritems():
        if (not geneToExonHits.has_key(accNum)):
            geneToExonHits[accNum] = [exonHit]
        else:
            geneToExonHits[accNum].append(exonHit)

    print 'Number of genes in geneToExonHits: ', len(geneToExonHits)
    totalExonHits = 0
    totalExonReadNum = 0
    for gene in geneToExonHits:
        totalExonHits += len(geneToExonHits[gene])
        for exonHit in geneToExonHits[gene]:
            totalExonReadNum += exonHit.readNum()
    print 'Total number of exon hits by geneToExonHits: ', totalExonHits
    print 'Total exon readNum: ', totalExonReadNum

    return geneToExonHits, geneToIntronHits
#########################################################################################
    
