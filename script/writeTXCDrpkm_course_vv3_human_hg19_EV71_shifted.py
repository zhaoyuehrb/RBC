import sys
import math
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
def loadAccNumToGene():
    infile = open('/home/hguo/Documents/processing/filteredGenesDetails_human_240118.txt')
    infileLines = infile.readlines()
    infile.close()

    accNumToGene = {}
    for line in infileLines[1:]:
        line = line.rstrip('\n')
        cols = line.split('\t')
        name = cols[0]
        accNum = cols[1]
        accNumToGene[accNum] = name
    print 'Number of genes loaded in accNumToGene: ', len(accNumToGene)
    return accNumToGene
##########################################################################################

folder = sys.argv[1]
shift_len = sys.argv[2] #i.e. 'shift34'
barcode = folder.split('/')[-2]
filename = folder + '../gene_TXCDUTR_ReadOutput_'+barcode+'_'+shift_len+'.txt'
#chrE_readNum = int(sys.argv[2])
with open(folder + 'blDetails.txt','r') as file:
    for line in file:
        if 'Combined chrE reads:' in line:
            break

chrE_readNum = int(line.split(' ')[-1])




date = filename.split('_')[-6]
lane = '_'.join(filename.split('_')[-5:])
outfileString =folder+'../'+shift_len+'_'+'geneTXCD_RPKMoutput_EVadjusted_'+barcode+'.txt'

## load reads, including genes with zero read count
accNumToReads, totalReads = loadAllAccNumToReads(filename)
## load lengths
accNumToLengths = loadAccNumToLengths()
## load gene names
accNumToGene = loadAccNumToGene()

## calculate rpkm
outfile = open(outfileString, 'w')
firstLineList = ['AccNum', 'GeneName', 'txReads', 'txRPKM', 'cdReads', 'cdRPKM']
firstLineString = '\t'.join(firstLineList)
outfile.write(firstLineString+'\n')

linesWritten = 1
for accNum in accNumToReads:
    txReads, cdReads = accNumToReads[accNum][0], accNumToReads[accNum][1]
    if txReads != 0:
        tx_rpkm = calculateRPKM(txReads, accNumToLengths[accNum][0], totalReads + chrE_readNum)
    else:
        tx_rpkm = '-'
    if cdReads != 0:
        cd_rpkm = calculateRPKM(cdReads, accNumToLengths[accNum][1], totalReads + chrE_readNum)
    else:
        cd_rpkm = '-'

    outList = [accNum, accNumToGene[accNum], txReads, tx_rpkm, cdReads, cd_rpkm]
    strOutList = map(lambda i: str(i), outList)
    outListString = '\t'.join(strOutList)
    outfile.write(outListString+'\n')
    linesWritten += 1
outfile.close()
print 'Number of lines written: ', linesWritten




