import sys
import math
import os
#############################################################################
def loadAllAccNumToReads(filePath):
    infile = open(filePath)
    infileLines = infile.readlines()
    infile.close()

    accNumToReads = {}
    totalMappedReads = 0
    for line in infileLines[1:]:
        line = line.rstrip('\n')
        cols = line.split('\t')
        
        accNum = cols[0]
        txReads = float(cols[4])
        cdReads = float(cols[5])
        intronReads = float(cols[3])
        totalMappedReads += txReads
        totalMappedReads += intronReads
        accNumToReads[accNum] = (txReads, cdReads)

    print 'Number of genes loaded in accNumToReads: ', len(accNumToReads)
    print 'Total number of mapped reads: ', totalMappedReads
    return accNumToReads, totalMappedReads
############################################################################
def loadAccNumToLengths():
    infile = open('/home/hguo/Documents/processing/filteredGenesDetails_human_240118.txt')
    infileLines = infile.readlines()
    infile.close()

    accNumToLengths = {}
    for line in infileLines[1:]:
        line = line.rstrip('\n')
        cols = line.split('\t')
        accNum = cols[1]
        mRNAlength = int(cols[3])
        ORFlength = int(cols[4])
        accNumToLengths[accNum] = (mRNAlength, ORFlength)
    print 'Number of genes loaded in accNumToLengths: ', len(accNumToLengths)
    return accNumToLengths
###############################################################################
def calculateRPKM(readNum, mRNAlength, totalReads):
    numerator = readNum*(10**9)
    denominator = totalReads*mRNAlength
    rpkm = numerator/float(denominator)
    return rpkm
##################################################################################

target_directory = sys.argv[1]

filename = target_directory+sys.argv[2]
date = filename.split('_')[-4]
lane = '_'.join(filename.split('_')[-3:])
outfileString = target_directory+'geneTXCD_RPKMoutput_'+date+'_'+lane

## load reads, including genes with zero read count
accNumToReads, totalReads = loadAllAccNumToReads(filename)
## load lengths
accNumToLengths = loadAccNumToLengths()

## calculate rpkm
outfile = open(outfileString, 'w')
firstLineList = ['AccNum', 'txReads', 'txRPKM', 'cdReads', 'cdRPKM']
firstLineString = '\t'.join(firstLineList)
outfile.write(firstLineString+'\n')

linesWritten = 1
for accNum in accNumToReads:
    txReads, cdReads = accNumToReads[accNum][0], accNumToReads[accNum][1]
    if txReads != 0:
        tx_rpkm = calculateRPKM(txReads, accNumToLengths[accNum][0], totalReads)
    else:
        tx_rpkm = '-'
    if cdReads != 0:
        cd_rpkm = calculateRPKM(cdReads, accNumToLengths[accNum][1], totalReads)
    else:
        cd_rpkm = '-'

    outList = [accNum, txReads, tx_rpkm, cdReads, cd_rpkm]
    strOutList = map(lambda i: str(i), outList)
    outListString = '\t'.join(strOutList)
    outfile.write(outListString+'\n')
    linesWritten += 1
outfile.close()
print 'Number of lines written: ', linesWritten




