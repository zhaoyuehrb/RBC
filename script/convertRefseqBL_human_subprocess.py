#############################################################
## need to change this file, to only take in one argument
## the bl file, and nothing else
## length etc, move to one file above
################################################################

import os, sys, time
import utilityModule as um
import subprocess
################################################################################3
def makeSplitfiles(infileLines, outfileString, numSplitfiles, readLength):

    splitfileNames = []
    for i in range(numSplitfiles):
        splitfileString = outfileString[:-5]+str(i)+'.'+readLength+'.bl'
        splitfileNames.append(splitfileString)

    splitfile = open(splitfileNames[0], 'w')
    index = 0
    count = 0
    for line in infileLines:
        splitfile.write(line)
        count += 1
        if count == 20000:
            splitfile.write('# blasting completed')
            splitfile.close()
            index += 1
            splitfile = open(splitfileNames[index], 'w')
            count = 0
    splitfile.write('# blasting completed')
    splitfile.close()
    return splitfileNames
################################################################################

## lines per file: 20000

## argument: queryBLfilename

queryFilePath = sys.argv[1] ## queryFilePath includes directory
#directory = '/lab/bartel/huili/Databases/hg18/'
directory = '/home/hguo/Downloads/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes/nibDirectory/'

lastSlash = queryFilePath.rfind('/')
inDirectory = queryFilePath[:lastSlash]

infile = open(queryFilePath)
infileLines = infile.readlines()
infile.close()

outfileString = queryFilePath[:-3]+'.fa'
outfile = open(outfileString, 'w')

length = queryFilePath.split('.')[-2]

numLines = len(infileLines[:-1])
if numLines <= 20000:
    needToSplit = False
else:
    needToSplit = True

if not needToSplit:
    linesWritten = 0
    for line in infileLines[:-1]:
        line = line.rstrip('\n')
        cols = line.split('\t')
        key = cols[0]
        chrom = cols[1]
        accNum = cols[2]
        genomeStart = int(cols[8])
        surrogateEnd = int(cols[9])
        if genomeStart < surrogateEnd:
            sense = '+'
        else:
            sense = '-'
        blockStarts = map(lambda i: int(i), cols[10].split(','))
        blockEnds = map(lambda i: int(i), cols[11].split(','))

        blockOne = um.Locus(chrom, blockStarts[0], blockEnds[0], sense)
        blockOneSeq = um.getSequence(blockOne, directory)
        blockTwo = um.Locus(chrom, blockStarts[1], blockEnds[1], sense)
        blockTwoSeq = um.getSequence(blockTwo, directory)
        actualSeq = blockOneSeq+blockTwoSeq

        header = '>'+key
        package = header+'\n'+actualSeq
        outfile.write(package+'\n')
        linesWritten += 2
    outfile.close()
    print '======= '+outfileString+' ============'
    print 'Number of linesWritten: ', linesWritten
    print
    trackfilePath = queryFilePath.rstrip('.bl') + '.track'
    trackfile = open(trackfilePath, 'w')
    trackfile.write('this file is done')
    trackfile.close()
    
else:
    print '======= '+outfileString+' ============'
    print '======= needs splitting =============='
    numSplitfiles = numLines/20000
    if numLines%20000 != 0:
        ## there is a remainder
        numSplitfiles += 1
    splitfileList = makeSplitfiles(infileLines[:-1], outfileString, numSplitfiles, length)
    print 'Number of splitBLfiles made: ', len(splitfileList)
    ## return a list of filenames, files already split

    ## if need to split, then need to make trackfiles
    ## for each splitfile, run same command, send output to /dev/null
    ## when each one splitfile completes, make a trackfile
    
    for sf in splitfileList:
        command = 'python convertRefseqBL_human_subprocess.py '+sf+' > /dev/null'
        output = subprocess.check_output(['bash', '-c', command])
    ## every 10 seconds, check number of trackfiles
    ## if number of trackfiles == number of splitfiles
    ## then all are done, then trigger combine
        
    triggerCombine = False
    while triggerCombine == False:
        time.sleep(10)
        trackCommand = 'ls '+inDirectory+'/*'+length+'.track -l | wc -l'
        trackOutput = subprocess.check_output(['bash', '-c', trackCommand])
        if int(trackOutput) == len(splitfileList):
            ## number of trackfiles = number of splitfiles
            ## can trigger combine
            triggerCombine = True
        else:
	    print "notequal"
        
    if triggerCombine:
        ## ie all the splitBLfiles have been converted to .fa files
        ## combfileString = outfileString from beginning

        resultfileList = map(lambda i: i[:-3]+'.fa', splitfileList)
        outlines = 0
        for rf in resultfileList:
            resultSplitfile = open(rf)
            resultLines = resultSplitfile.readlines()
            resultSplitfile.close()

            for rline in resultLines:
                outfile.write(rline)
                outlines += 1
        outfile.close()
        print '======= '+outfileString+' ============'
        print 'Number of outlinesWritten: ', outlines
        print

        ## delete resultsplitfiles
        for rf in resultfileList:
            os.remove(rf)
        ## delete splitfiles
        for sf in splitfileList:
            os.remove(sf)
        ## delete trackfiles
        trackfileList = map(lambda i: i[:-3]+'.track', splitfileList)
        for tf in trackfileList:
            os.remove(tf)

            
            
        
        
         
        
        
              
        



                
            
            
        
        

    
    
