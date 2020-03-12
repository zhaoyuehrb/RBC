import sys
import os
## input: inDirectory1, inDirectory2, outDirectory, nextIdNum
####################################################################
def loadCoordsToKeys(infileString):
    infile = open(infileString)
    infileLines = infile.readlines()
    infile.close()

    #print infileLines[-2]
    coordsToKey = {}
    coordsToLabels = {}
    for line in infileLines[:-1]:
        line = line.rstrip('\n')
        cols = line.split('\t')
        key = cols[0]
        chrom = cols[1]
        length = int(cols[3])
        rawStart = int(cols[8])
        rawEnd = int(cols[9])
        rrnaLabel = cols[10]
        otherrnaLabel = cols[11]
        if rawStart < rawEnd:
            sense = '+'
        else:
            sense = '-'
        coords = chrom+'('+sense+'):'+str(rawStart)+'_'+str(length)
        if not coords in coordsToKey:
            coordsToKey[coords] = key
        if not coords in coordsToLabels:
            coordsToLabels[coords] = (rrnaLabel, otherrnaLabel)
    print 'Number of coords loaded in coordsToKey: ', len(coordsToKey)
    print 'Number of coords loaded in coordsToLabels: ', len(coordsToLabels)
    return coordsToKey, coordsToLabels
###########################################################################

inDirectory1, inDirectory2 = sys.argv[1], sys.argv[2]
outDirectory = sys.argv[3]

allLengths = [18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 36]

dirfiles1 = os.listdir(inDirectory1)
avaiLengths1 = []
for dirfile1 in dirfiles1:
    if dirfile1[-3:] == '.bl':
        if int(dirfile1.split('.')[-2]) in allLengths:
            avaiLengths1.append(int(dirfile1.split('.')[-2]))

dirfiles2 = os.listdir(inDirectory2)
avaiLengths2 = []
for dirfile2 in dirfiles2:
    if dirfile2[-3:] == '.bl':
        if int(dirfile2.split('.')[-2]) in allLengths:
            avaiLengths2.append(int(dirfile2.split('.')[-2]))

lengths = []
onlyInDirectory1 = []
onlyInDirectory2 = []
for avaiLength in avaiLengths1:
    if avaiLength in avaiLengths2:
        lengths.append(avaiLength)
    else:
        onlyInDirectory1.append(avaiLength)
for avaiLength in avaiLengths2:
    if avaiLength not in lengths:
        onlyInDirectory2.append(avaiLength)

print lengths
print onlyInDirectory1
print onlyInDirectory2

nextIdNum = 1
for i in lengths:
    corefilename1 = inDirectory1[:-1].split('/')[-1]
    corefilename2 = inDirectory2[:-1].split('/')[-1]
    outCorefilename = outDirectory[:-1].split('/')[-1]

    infileString1 = inDirectory1+corefilename1+'.'+str(i)+'.bl'
    infileString2 = inDirectory2+corefilename2+'.'+str(i)+'.bl'
    infileString1 = inDirectory1+corefilename1+'_seed29_genome.'+str(i)+'.bl'
    infileString2 = inDirectory2+corefilename2+'_seed29_genome.'+str(i)+'.bl'


    coordsToKey1, coordsToLabels1 = loadCoordsToKeys(infileString1)
    coordsToKey2, coordsToLabels2 = loadCoordsToKeys(infileString2)

    ## get all keys
    coordsToKey = {}
    for coord in coordsToKey1:
        coordsToKey[coord] = [coordsToKey1[coord]]
        if coord in coordsToKey2:
            coordsToKey[coord].append(coordsToKey2[coord])

    for coord in coordsToKey2:
        if coord not in coordsToKey:
            ## coord was not in first .bl file
            coordsToKey[coord] = [coordsToKey2[coord]]

    print 'Number of coords in combined coordsToKey: ', len(coordsToKey)

    ## get all labels
    coordsToLabels = {}
    for coord in coordsToLabels1:
        coordsToLabels[coord] = [coordsToLabels1[coord]]
        if coord in coordsToLabels2:
            coordsToLabels[coord].append(coordsToLabels2[coord])

    for coord in coordsToLabels2:
        if coord not in coordsToLabels:
            ## coord was not in first .bl file
            coordsToLabels[coord] = [coordsToLabels2[coord]]

    print 'Number of coords in combined coordsToLabels: ', len(coordsToLabels)

    outfileString = outDirectory+outCorefilename+'.'+str(i)+'.bl'
    outfile = open(outfileString, 'w')

    linesWritten= 0
    for coord in coordsToKey:
        readNum1 = int(coordsToKey[coord][0].split('_')[1])
        if len(coordsToKey[coord]) > 1:
            readNum2 = int(coordsToKey[coord][1].split('_')[1])
        else:
            readNum2 = 0
        sumReads = readNum1 + readNum2
        chrom = coord.split('(')[0]
        sense = coord.split('(')[1][0]
        rawStart = coord.split('_')[0].split(':')[1]
        length = coord.split('_')[1]
        if sense == '+':
            rawEnd = int(rawStart) + int(length) - 1
        elif sense == '-':
            rawEnd = int(rawStart) - int(length) + 1
        newId = str(nextIdNum)+'_'+str(sumReads)
        outList = [newId, chrom, '-', length, '0', '0', '1', length, rawStart, str(rawEnd)]
        labels = [coordsToLabels[coord][0][0], coordsToLabels[coord][0][1]]
        outList.extend(labels)
        outListString = '\t'.join(outList)
        outfile.write(outListString+'\n')
        linesWritten += 1
        nextIdNum += 1
    print 'Number of linesWritten to file '+str(i)+': ', linesWritten
    outfile.write('# blasting completed')
    print '====== Last line written ======'
    outfile.close()
    print '=============================================='  

for i in onlyInDirectory1:
    corefilename1 = inDirectory1[:-1].split('/')[-1]
    infileString1 = inDirectory1+corefilename1+'_seed29_genome.'+str(i)+'.bl'
    outCorefilename = outDirectory[:-1].split('/')[-1]

    infile1 = open(infileString1)
    infileLines1 = infile1.readlines()
    infile1.close()

    outfileString = outDirectory+outCorefilename+'.'+str(i)+'.bl'
    outfile = open(outfileString, 'w')
    linesWritten= 0
    for line in infileLines1[:-1]:
        line = line.rstrip('\n')
        cols = line.split('\t')
        oldId = cols[0]
        readNum = oldId.split('_')[1]
        newId = str(nextIdNum)+'_'+readNum
        details = cols[1:]
        outList = [newId]
        outList.extend(details)
        outListString = '\t'.join(outList)
        outfile.write(outListString+'\n')
        linesWritten += 1
        nextIdNum += 1
    print 'Number of linesWritten to file '+str(i)+': ', linesWritten
    outfile.write('# blasting completed')
    print '====== Last line written ======'
    outfile.close()
    print '=============================================='

for i in onlyInDirectory2:
    corefilename2 = inDirectory2[:-1].split('/')[-1]
    infileString2 = inDirectory2+corefilename2+'_seed29_genome.'+str(i)+'.bl'
    outCorefilename = outDirectory[:-1].split('/')[-1]

    infile2 = open(infileString2)
    infileLines2 = infile2.readlines()
    infile2.close()
    linesWritten= 0
    outfileString = outDirectory+outCorefilename+'.'+str(i)+'.bl'
    outfile = open(outfileString, 'w')

    for line in infileLines2[:-1]:
        line = line.rstrip('\n')
        cols = line.split('\t')
        oldId = cols[0]
        readNum = oldId.split('_')[1]
        newId = str(nextIdNum)+'_'+readNum
        details = cols[1:]
        outList = [newId]
        outList.extend(details)
        outListString = '\t'.join(outList)
        outfile.write(outListString+'\n')
        linesWritten += 1
        nextIdNum += 1
    print 'Number of linesWritten to file '+str(i)+': ', linesWritten
    outfile.write('# blasting completed')
    print '====== Last line written ======'
    outfile.close()
    print '=============================================='
