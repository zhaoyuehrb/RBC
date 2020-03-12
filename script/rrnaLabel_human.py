## inital file
## load the length bl file as lines, not into a hitsCollection
## supply directory

## call next file, providing directory and length
## wait for keysList to be returned

## with keysList, go through lines and modify cols[10] when key is in keysList

import sys
import os
import ncRNAoverlapKeys_human as ro
############################################################################

#directory = sys.argv[1]
length = sys.argv[1]#2]
rnaType = sys.argv[2]#3]
#corefilename = directory[:-1].split('/')[-1]
prefix = sys.argv[3]
directory = sys.argv[4]


corefilename = prefix+"_seed29_genome"
infileString = directory+corefilename+'.'+length+'.bl'
infile = open(infileString)
infileLines = infile.readlines()
infile.close()

print '======== '+corefilename+'.'+length+'.bl =========='
print 'Number of lines in infile: ', len(infileLines)

rrnaKeys = ro.getKeys(directory, length, rnaType)
print 'Number of '+rnaType+' keys: ', len(rrnaKeys)
#print rrnaKeys
outfileString = directory+corefilename+'.'+length+'.rll'
outfile = open(outfileString, 'w')

linesWritten = 0
for line in infileLines[:-1]:
    line = line.rstrip('\n')
    cols = line.split('\t')
    chrom = cols[1]
    length = int(cols[3])
    rawStart = int(cols[8])
    rawEnd = int(cols[9])
    if rawStart < rawEnd:
        sense = '+'
    elif rawStart > rawEnd:
        sense = '-'
    collatedKey = chrom+'('+sense+'):'+str(rawStart)+'_'+str(length)
    if collatedKey in rrnaKeys:
        cols[10] = rnaType
    outString = '\t'.join(cols)
    outfile.write(outString+'\n')
    linesWritten += 1
        


##linesWritten = 0
##for line in infileLines[:-1]:
##    line = line.rstrip('\n')
##    cols = line.split('\t')
##    entryID = cols[0]
##    if entryID in rrnaKeys:
##        cols[10] = rnaType
##    outString = '\t'.join(cols)
##    outfile.write(outString+'\n')
##    linesWritten += 1

print 'Number of lines written: ', linesWritten
outfile.write('# blasting completed')
print '============= last line written ================'
outfile.close()



