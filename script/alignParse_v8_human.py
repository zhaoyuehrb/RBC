## the three .map files from bowtieFlow shld be in the same directory
## provide directory as argument
## eg s_8_old_seed21.0.map

import sys
import os



in_type = sys.argv[1] #RNAseq or RPF


# eg directory = bartel_lab.3\huili\bowtie\bowtie-0.9.9\datasets\s_8_seed21/
#corefilename = s_8_seed21
header = sys.argv[2] # folder_name and filename b4 "_seed..." e.g. CGATGT-s_6_1
directory = sys.argv[3]
corefilename =header+"_seed29_genome"
directory = directory+header+"/"
# eg .map file = s_8_seed21_zeroMM.map
#corefilename = directory[:-1].split('/')[-1]

infilePaths = [directory+corefilename+'_zeroMM.map', directory+corefilename+'_oneMM.map', directory+corefilename+'_twoMM.map']
#infilePaths = [directory+corefilename+'_zeroMM.map', directory+corefilename+'_oneMM.map']

infileLines = []
for i in range(len(infilePaths)):
    infile = open(infilePaths[i])
    partfileLines = infile.readlines()
    infile.close()
    print 'Number of lines in part '+str(i)+': ', len(partfileLines)
    infileLines.extend(partfileLines)
print 'Number of lines in total: ', len(infileLines)

###############################################
### check if lines include linebreaks!!
###############################################

## fill in the rest from alignParse_v2.py

##chrDict = {'000067': 'chr1', '000068': 'chr2', '000069': 'chr3', '000070': 'chr4', \
##           '000071': 'chr5', '000072': 'chr6', '000073': 'chr7', '000074': 'chr8', \
##           '000075': 'chr9', '000076': 'chr10', '000077': 'chr11', '000078': 'chr12', \
##           '000079': 'chr13', '000080': 'chr14', '000081': 'chr15', '000082': 'chr16', \
##           '000083': 'chr17', '000084': 'chr18', '000085': 'chr19', '000086': 'chrX', \
##           '000087': 'chrY', '005089': 'chrM'} 


collatedKeyToReadNum = {}
uniqueID = 0

for line in infileLines:
    line = line.rstrip('\n')
    cols = line.split('\t')
    sense = cols[1]
##    chromKey = cols[2].split('ref')[1].split('.')[0].split('_')[1]
##    chrom = chrDict[chromKey]
    chrom = cols[2]
    offset = int(cols[3])
    read = cols[4]
    if sense == '+':
        realOffset = offset
    elif sense == '-':
        realOffset = offset + len(read) - 1
    key = chrom+'('+sense+'):'+str(realOffset)+'_'+str(len(read))

    if not collatedKeyToReadNum.has_key(key):
        collatedKeyToReadNum[key] = 1
    else: # already has this key
        collatedKeyToReadNum[key] += 1

print 'Length of collatedKeyToReadNum before output length files: ', len(collatedKeyToReadNum)
##count = 0
##for key in collatedKeyToReadNum:
##    print key, collatedKeyToReadNum[key]
##    count += 1
##    if count == 20:
##        break

lengthToKeys = {}
for entry in collatedKeyToReadNum:
    length = entry.split('_')[1]
    if not lengthToKeys.has_key(length):
        lengthToKeys[length] = [entry]
    else:
        # already has this length, just append
        lengthToKeys[length].append(entry)

print 'Number of lengths in lengthToKeys: ', len(lengthToKeys)

for length in lengthToKeys:
    outfilename = '.'.join([directory+corefilename, length, 'bl'])
    outfile = open(outfilename, 'w')
    linesWritten = 0
    for key in lengthToKeys[length]:
        uniqueID += 1
        # fill in to outString
        readNum = collatedKeyToReadNum[key]
        chrom = key.split('(')[0]
        sense = key.split('(')[1][0]
        rawStart = key.split('_')[0].split(':')[1]
        if sense == '+':
            rawEnd = int(rawStart) + int(length) - 1
        elif sense == '-':
            rawEnd = int(rawStart) - int(length) + 1
        newKey = str(uniqueID)+'_'+str(readNum)
        outList = [newKey, chrom, '-', length, '0', '0', '1', length, rawStart, str(rawEnd), 'NA', 'NA']
        outString = '\t'.join(outList)
        outfile.write(outString+'\n')
        linesWritten += 1

    print 'Number of linesWritten to file '+length+': ', linesWritten
    outfile.write('# blasting completed')
    print '====== Last line written ======'
    outfile.close()
    print '=============================================='
