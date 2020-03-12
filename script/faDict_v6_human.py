import sys, time

######################################################################################################
def getChromCollatedKeys(chrom, directory,corefilename):
    hl_path = '/home/hguo/Documents/indexes/hg19/refseq/chromGenes_240118_human/'
    fafile = hl_path+chrom+'_genes.fa'
    rfCoordfile = hl_path+chrom+'_refseqCoords.txt'

    ## populate chromBaseToGenomeCoord dict with fafile

    chromBaseToGenomeCoord = {}
    infile = open(fafile)
    header = infile.readline()
    infileBases = infile.read()
    infile.close()

    for i in range(len(infileBases)):
        chromBaseToGenomeCoord[i] = None
    print 'Number of entries in chromBaseToGenomeCoord for '+chrom+' : ', len(chromBaseToGenomeCoord)

    inRFfile = open(rfCoordfile)
    inRFlines = inRFfile.readlines()
    inRFfile.close()

    annotatedGenes = 0
    for line in inRFlines:
        line = line.rstrip('\n')
        cols = line.split('\t')
        #accNum = cols[0]
        accNum = ''.join(cols[0].split('_'))    # remove the underscore to prevent complications from parsing collated key later
        txStart = int(cols[1])
        txEnd = int(cols[2])
        rStart = int(cols[3])
        rEnd = int(cols[4])
        exStarts = map(lambda i: int(i), cols[5].split(','))
        exEnds = map(lambda i: int(i), cols[6].split(','))
        numExons = len(exStarts)
        refStarts = map(lambda i: int(i), cols[7].split(','))
        refEnds = map(lambda i: int(i), cols[8].split(','))

        if txStart < txEnd:
            sense = '+'
        else:
            sense = '-'
        if sense == '+':
            for exon in range(numExons):
                currentGenomeCoord = exStarts[exon]
                for base in range(refStarts[exon], refEnds[exon]+1):
                    chromBaseToGenomeCoord[base] = (currentGenomeCoord, sense, accNum)
                    currentGenomeCoord += 1
        elif sense == '-':
            for exon in range(numExons):
                currentGenomeCoord = exStarts[exon]
                for base in range(refStarts[exon], refEnds[exon]+1):
                    chromBaseToGenomeCoord[base] = (currentGenomeCoord, sense, accNum)
                    currentGenomeCoord -= 1
        annotatedGenes += 1

    annotated = 0
    baselist = []
    for base, genomeCoord in chromBaseToGenomeCoord.iteritems():
        if genomeCoord != None:
            annotated += 1
            baselist.append((base, genomeCoord))
    print 'Number of annotated bases: ', annotated, len(baselist)
    print 'Number of annotated genes: ', annotatedGenes

    ## converting alignments in .map file to collatedKeys
    #corefilename = directory[:-1].split('/')[-1]
    infilePaths = [directory+corefilename+'_zeroMM.map', directory+corefilename+'_oneMM.map', directory+corefilename+'_twoMM.map']

    inMapfileLines = []
    for i in range(len(infilePaths)):
        infile = open(infilePaths[i])
        partfileLines = infile.readlines()
        infile.close()
        print 'Number of lines in part '+str(i)+': ', len(partfileLines)
        inMapfileLines.extend(partfileLines)
    print 'Number of lines in total: ', len(inMapfileLines)

    collatedKeyToReadNum = {}
    complement = []
    chromAlign = 0
    for mapline in inMapfileLines:
        mapline = mapline.rstrip('\n')
        cols = mapline.split('\t')
        refsense = cols[1]
        chromosome = cols[2].split('_')[0]
        refOffset = int(cols[3])
        read = cols[4]
        if refsense == '+':
            if chromosome == chrom: ## take the alignments from this chrom only
                # convert refOffset to genomeCoords
                chromAlign += 1
                genomeRawStart = chromBaseToGenomeCoord[refOffset][0]
                sense = chromBaseToGenomeCoord[refOffset][1]
                accNum = chromBaseToGenomeCoord[refOffset][2]
                key = chrom+'('+sense+')'+accNum+':'+str(genomeRawStart)+'_'+str(len(read))

                if not collatedKeyToReadNum.has_key(key):
                    collatedKeyToReadNum[key] = 1
                else:
                    collatedKeyToReadNum[key] += 1
            else: # alignments for other chromosomes
                pass
        else:
            complement.append(read)
            
    print 'Number of alignments to complement of transcripts: ', len(complement)
    print 'Number of alignments that match to transcripts: ', chromAlign
    print complement[:5]
    print 'Length of collatedKeyToReadNum before output length files: ', len(collatedKeyToReadNum), time.asctime()

    return collatedKeyToReadNum, chromAlign
#######################################################################################################
directory = sys.argv[1]
#corefilename = directory[:-1].split('/')[-1]
corefilename = sys.argv[2]
chromosomes = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']
chromosomes = map(lambda i: 'chr'+i, chromosomes)
#chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7']
lengthToKeys = {}
chromAligns = 0
chromToCollatedKeyToReadNum = {}
for chrom in chromosomes:
    collatedKeyToReadNum, chromAlign = getChromCollatedKeys(chrom, directory,corefilename)
    chromAligns += chromAlign
    chromToCollatedKeyToReadNum[chrom] = collatedKeyToReadNum
    for entry in collatedKeyToReadNum:
        length = entry.split('_')[1]
        if not lengthToKeys.has_key(length):
            lengthToKeys[length] = [entry]
        else:
            # already has this length, just append
            lengthToKeys[length].append(entry)
    print 'Number of lengths in lengthToKeys: ', len(lengthToKeys)

print 'Total number of alignments to chromGenes, sense: ',  chromAligns
print 'Number of chromosome entries in chromToCollatedKeyToReadNum: ', len(chromToCollatedKeyToReadNum)

print lengthToKeys.keys()
totalAligns = 0
for length, details in lengthToKeys.iteritems():
    totalAligns += len(details)
print 'Number of alignments in lengthToKeys: ', totalAligns

totalReads = 0
uniqueID = 0
for length in lengthToKeys:
    outfilename = '.'.join([directory+corefilename, length, 'bl'])
    outfile = open(outfilename, 'w')
    linesWritten = 0
    for key in lengthToKeys[length]:
        uniqueID += 1
        # fill in to outString
        chrom = key.split('(')[0]
        readNum = chromToCollatedKeyToReadNum[chrom][key]
        sense = key.split('(')[1][0]
        squishedAccNum = key.split(':')[0].split(')')[1]
        accNum = 'NM_'+squishedAccNum.split('NM')[1]
        rawStart = key.split('_')[0].split(':')[1]
        if sense == '+':
            rawEnd = int(rawStart) + int(length) - 1
        elif sense == '-':
            rawEnd = int(rawStart) - int(length) + 1
        newKey = str(uniqueID)+'_'+str(readNum)
        outList = [newKey, chrom, accNum, length, '0', '0', '1', length, rawStart, str(rawEnd), 'NA', 'NA']
        outString = '\t'.join(outList)
        outfile.write(outString+'\n')
        linesWritten += 1
        totalReads += readNum
    print 'Number of linesWritten to file '+length+': ', linesWritten
    outfile.write('# blasting completed')
    print '====== Last line written ======'
    outfile.close()
    print '=============================================='

print 'Total reads: ', totalReads



