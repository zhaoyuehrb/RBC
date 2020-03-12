import sys, os
########################################################################
chromosomes = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']
hl_path = '/home/hguo/Documents/indexes/hg19/refseq/chromGenes_240118_human/'
rfCoordfiles = map(lambda i: hl_path+'/chr'+i+'_refseqCoords.txt', chromosomes)

inRFlines = []
for rfCoordfile in rfCoordfiles:
    inRFfile = open(rfCoordfile)
    partfileLines = inRFfile.readlines()
    inRFlines.extend(partfileLines)
    inRFfile.close()

geneToExStartsEnds = {}
for line in inRFlines:
    line = line.rstrip('\n')
    cols = line.split('\t')
    accNum = cols[0]
    exStarts = map(lambda i: int(i), cols[5].split(','))
    exEnds = map(lambda i: int(i), cols[6].split(','))
    geneToExStartsEnds[accNum] = (exStarts, exEnds)

print 'Number of entries in geneToExStartsEnds: ', len(geneToExStartsEnds)

directory = sys.argv[1]
#corefilename = directory[:-1].split('/')[-1]
corefilename = sys.argv[2]

lenTuple = (18,37)
lengths = []
for i in range(lenTuple[0], lenTuple[1]):
    lengths.append(i)
dirfiles = os.listdir(directory)
avaiLengths = []
for dirfile in dirfiles:
    if dirfile[-3:] == '.bl':
        if int(dirfile.split('.')[-2]) in lengths:
            avaiLengths.append(int(dirfile.split('.')[-2]))
#lengths = [18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 36]
#lengths = [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 36]
for length in avaiLengths:
    infilePath = directory+corefilename+'.'+str(length)+'.bl'
    infile = open(infilePath)
    infileLines = infile.readlines()
    infile.close()

    outfilePath = directory+corefilename+'.'+str(length)+'.blbed2'
    outfile = open(outfilePath, 'w')
    linesWritten = 0
    for line in infileLines[:-1]:
        line = line.rstrip('\n')
        cols = line.split('\t')
        accNum = cols[2]
        genomeStart = int(cols[8])
        length = int(cols[7])
        surrogateEnd = int(cols[9])
        exStarts = geneToExStartsEnds[accNum][0]
        exEnds = geneToExStartsEnds[accNum][1]
        if genomeStart < surrogateEnd:
            sense = '+'
        else:
            sense = '-'
        if sense == '+':
            blockStarts = [genomeStart]
            #exonNum = 0
            for i in range(len(exStarts)):
            #for exStart in exStarts:
                #exonNum += 1
                if genomeStart < exStarts[i]:
                    # hit started in previous exon
                    # therefore 2nd block starts with current exStart
                    blockStarts.append(exStarts[i])
                    blockEnds = [exEnds[i-1]]
                    firstBlockSize = blockEnds[0] - blockStarts[0] + 1
                    secondBlockSize = length - firstBlockSize
                    secondEnd = secondBlockSize - 1 + blockStarts[1]
                    blockEnds.append(secondEnd)
                    break
        elif sense == '-':
            blockStarts = [genomeStart]
            #exonNum = 0
            for i in range(len(exStarts)):
            #for exStart in exStarts:
                #exonNum += 1
                if genomeStart > exStarts[i]:
                    # hits started in previous exon
                    # therefore 2nd block starts with current exStart
                    blockStarts.append(exStarts[i])
                    blockEnds = [exEnds[i-1]]
                    firstBlockSize = blockStarts[0] - blockEnds[0] + 1
                    secondBlockSize = length - firstBlockSize
                    secondEnd = blockStarts[1] + 1 - secondBlockSize
                    blockEnds.append(secondEnd)
                    break

        ## now blockStarts and blockEnds lists filled
        ## modify cols[10] and cols[11]
        blockStarts = map(lambda i: str(i), blockStarts)
        blockEnds = map(lambda i: str(i), blockEnds)
        cols[10] = ','.join(blockStarts)
        cols[11] = ','.join(blockEnds)
        outString = '\t'.join(cols)
        outfile.write(outString+'\n')
        linesWritten += 1
    print 'Number of linesWritten to file '+str(length)+': ', linesWritten
    outfile.write('# blasting completed')
    print '====== Last line written ======'
    outfile.close()
    print '=============================================='
    
    
    

    
