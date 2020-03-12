import sys
import os
#######################################################################
directory = sys.argv[1]
dat_set = sys.argv[2]
outPath = sys.argv[3] #'../wigFilesAHSP/'
if30 = sys.argv[4]
#totalMappedReads = float(sys.argv[2])
#lengths = [int(sys.argv[2])]

#len_str = sys.argv[2]   # format as [2,3,4]
#lengths = map(int,len_str.strip('[]').split(','))

totalReads = 0

corefilename = directory[:-1].split('/')[-1]
print(corefilename)
#chrList = '['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'M']
chrList = ['1','9']
#chrList = ['M']

chrList = map(lambda i: 'chr'+i, chrList)
chromToWigStartToScore = {}
for chrom in chrList:
    chromToWigStartToScore[chrom] = {}

if if30 == '30':
    lengths = [30]
else:    
    lengths = [18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 36]

for length in lengths:
    infileString = directory+corefilename+'_seed29_genome.'+str(length)+'.bl'
    infile = open(infileString)
    infileLines = infile.readlines()
    infile.close()

    for line in infileLines[:-1]:
        line = line.rstrip('\n')
        cols = line.split('\t')
        readNum = int(cols[0].split('_')[1])
        chrom = cols[1]
        rawStart = int(cols[8])
        rawEnd = int(cols[9])
        
        entryID = cols[0]
        readNum = int(entryID.split('_')[1])
        totalReads += readNum
        if chrom in chrList:
            if rawStart < rawEnd:
                # sense is +
                for position in xrange(rawStart, rawEnd+1):
                    # label, make one-based!
                    wigStart = position + 1
                    if not chromToWigStartToScore[chrom].has_key(wigStart):
                       chromToWigStartToScore[chrom][wigStart] = readNum
                    else:
                       chromToWigStartToScore[chrom][wigStart] += readNum


            elif rawStart > rawEnd:
                # sense is -
                for position in xrange(rawStart, rawEnd-1, -1):
                    ## label, make one-based!
                    wigStart = position + 1
                    if not chromToWigStartToScore[chrom].has_key(wigStart):
                        chromToWigStartToScore[chrom][wigStart] = readNum
                    else:
                        chromToWigStartToScore[chrom][wigStart] += readNum

                ## change to only record start/end
## refseq
ref_dir = directory+corefilename+'_seed29_refseq/'
blFiles = []
for file in os.listdir(ref_dir):
    if file.endswith('.bl'):
        blFiles.append(file)

for fileName in blFiles:
    
    infileString = ref_dir + fileName
    infile = open(infileString)
    infileLines = infile.readlines()
    infile.close()

    for line in infileLines[:-1]:
        line = line.rstrip('\n')
        cols = line.split('\t')
        readNum = int(cols[0].split('_')[1])
        chrom = cols[1]
        rawStarts = (cols[10]).split(',')
        rawEnds = (cols[11]).split(',')
        
        entryID = cols[0]
        readNum = int(entryID.split('_')[1])
        totalReads += readNum
        if chrom in chrList:
            for idx in range(0,2):
                rawStart = int(rawStarts[idx])
                rawEnd = int(rawEnds[idx])
                if rawStart < rawEnd:
                    # sense is +
                    for position in xrange(rawStart, rawEnd+1):
                        # label, make one-based!
                        wigStart = position + 1
                        if not chromToWigStartToScore[chrom].has_key(wigStart):
                           chromToWigStartToScore[chrom][wigStart] = readNum
                        else:
                           chromToWigStartToScore[chrom][wigStart] += readNum


                elif rawStart > rawEnd:
                    # sense is -
                    for position in xrange(rawStart, rawEnd-1, -1):
                        ## label, make one-based!
                        wigStart = position + 1
                        if not chromToWigStartToScore[chrom].has_key(wigStart):
                            chromToWigStartToScore[chrom][wigStart] = readNum
                        else:
                            chromToWigStartToScore[chrom][wigStart] += readNum

                ## change to only record start/end




## output

outfileString = outPath+'bigWiggle_all_'+corefilename+'.wig'
outfile = open(outfileString, 'w')
#header = 'track type=wiggle_0 name=mir155KORNA description=m155_32_RNA visibility=full color=250,150,255 altColor=250,150,255 priority=10'
#secondHeader = 'variableStep chrom='+chromosome+' span=1'
#outfile.write(header+'\n')
#outfile.write(secondHeader+'\n')

#len_str = len_str.strip('[]').replace(', ','')
outfile.write('browser position chr16:31,539,221-31,540,124\n')
outfile.write('track type=wiggle_0 name='+dat_set+' description='+dat_set+'\n')

for chrom in chromToWigStartToScore:
    secondHeader = 'variableStep chrom='+chrom+' span=1'
    outfile.write(secondHeader+'\n')
    linesWritten = 1
    baseKeys = chromToWigStartToScore[chrom].keys()[:]
    baseKeys.sort()
    for base in baseKeys:
        outList = [str(base)]
        rawScore = chromToWigStartToScore[chrom][base]
        #normScore = (rawScore/float(totalMappedReads))*1000000
        #normScore = (rawScore/float(totalReads))*1000000
        normScore = rawScore
        outList.append(str(normScore))
        outListString = '\t'.join(outList)
        outfile.write(outListString+'\n')
        linesWritten += 1
    print 'Number of lines written for '+chrom+': ', linesWritten
outfile.close()
                    
                

                
