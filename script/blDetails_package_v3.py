import sys, os
#################################################################################################
def getDetails(infileString):
    infile = open(infileString)
    infileLines = infile.readlines()
    infile.close()

    totalReads = 0   ## this includes chrM and chrY !!!!
    rRNAreads = 0
    otherRNAreads = 0

    chrMreads = 0
    chrYreads = 0

    potentialGeneReads = 0
    
    ## mitochondria-related
    mt_rRNA_reads = 0
    mt_tRNA_reads = 0
    mt_ORF_reads = 0    

    totalHits = 0   ## this includes chrM and chrY !!!!
    rRNAhits = 0
    otherRNAhits = 0

    chrMhits = 0
    chrYhits = 0 

    potentialGeneHits = 0
    
    ## mitochondria-related
    mt_rRNA_hits = 0
    mt_tRNA_hits = 0
    mt_ORF_hits = 0

    for line in infileLines[:-1]:
        line = line.rstrip('\n')
        cols = line.split('\t')

        entryID = cols[0]
        chrom = cols[1]
        rRNAlabel = cols[10]
        otherlabel = cols[11]
        readNum = int(entryID.split('_')[1])

        totalReads += readNum
        totalHits += 1
        if rRNAlabel == 'rRNA':
            rRNAreads += readNum
            rRNAhits += 1
        if otherlabel != 'NA':
            otherRNAreads += readNum
            otherRNAhits +=1

        if chrom == 'chrM':
            chrMreads += readNum
            chrMhits += 1
            if rRNAlabel == 'mt-rRNA':
                mt_rRNA_reads += readNum
		mt_rRNA_hits += 1
	    elif rRNAlabel == 'mt-tRNA':
		mt_tRNA_reads += readNum
		mt_tRNA_hits += 1
	    elif rRNAlabel == 'NA':
		mt_ORF_reads += readNum
		mt_ORF_hits += 1
				
        if chrom == 'chrY':
            chrYreads += readNum
            chrYhits += 1
        if rRNAlabel == 'NA' and otherlabel == 'NA' and chrom != 'chrM' and chrom != 'chrY':
            ## this hit could potentially be taken for a gene
            potentialGeneReads += readNum
            potentialGeneHits += 1
        
    return totalReads, rRNAreads, otherRNAreads, chrMreads, chrYreads, potentialGeneReads, mt_rRNA_reads, mt_tRNA_reads, mt_ORF_reads, totalHits, rRNAhits, otherRNAhits, chrMhits, chrYhits, potentialGeneHits, mt_rRNA_hits, mt_tRNA_hits, mt_ORF_hits
#########################################################################################################################


#target_dir = os.path.join(script_dir,rel_path)

rootdirectory = sys.argv[1]
genomeDirectoryCore = (rootdirectory[:-1].split('/')[-1])
genomeDirectory = rootdirectory
corefilename_genome = genomeDirectory[:-1].split('/')[-1] + "_seed29_genome"


# genomeDirectoryCore = (rootdirectory[:-1].split('/')[-1])+'_genome/'
# genomeDirectory = os.path.join(script_dir,rel_path)#rootdirectory+genomeDirectoryCore
# corefilename_genome = "CGATGT-s_6_1_seed29_genome"#genomeDirectory[:-1].split('/')[-1]


lenRange = (18,37)
lengths = []
for i in range(lenRange[0], lenRange[1]):
	lengths.append(i)
	
#lengths = [18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 36]

dirfiles = os.listdir(genomeDirectory)
avaiLengths = []
for dirfile in dirfiles:
    if dirfile[-3:] == '.bl':
        if int(dirfile.split('.')[-2]) in lengths:
            avaiLengths.append(int(dirfile.split('.')[-2]))

combineReads = 0
combineHits = 0

combineYreads = 0
combineMreads = 0
combinePotentialGeneReads = 0

combineRRNA = 0
combineOTHER = 0

combine_mt_rRNA = 0
combine_mt_tRNA = 0
combine_mt_ORF = 0

lengthToTotal = {}
lengthTorRNA = {}
lengthToGenes = {}

lengthToCHROMMtotal = {}
lengthToCHROMMrRNA = {}
lengthToCHROMMtRNA = {}
lengthToCHROMM_ORF = {}

## check total number of reads, and total number of uniq reads
for avaiLength in avaiLengths:
    print '============ '+str(avaiLength)+ ' =============='
    ## open file and count readNum, numHits for rRNA
    ## return and add to dict
    filePath = genomeDirectory+corefilename_genome+'.'+str(avaiLength)+'.bl'
    totalReads, rRNAreads, otherRNAreads, chrMreads, chrYreads, potentialGeneReads, mt_rRNA_reads, mt_tRNA_reads, mt_ORF_reads, totalHits, rRNAhits, otherRNAhits, chrMhits, chrYhits, potentialGeneHits, mt_rRNA_hits, mt_tRNA_hits, mt_ORF_hits = getDetails(filePath)
    print 'Reads: ', totalReads, rRNAreads, otherRNAreads, chrMreads, chrYreads, potentialGeneReads, mt_rRNA_reads, mt_tRNA_reads, mt_ORF_reads
    print 'Hits: ', totalHits, rRNAhits, otherRNAhits, chrMhits, chrYhits, potentialGeneHits, mt_rRNA_hits, mt_tRNA_hits, mt_ORF_hits
    combineReads += totalReads
    combineHits += totalHits
    combineYreads += chrYreads
    combineMreads += chrMreads
    combinePotentialGeneReads += potentialGeneReads
    combineRRNA += rRNAreads
    combineOTHER += otherRNAreads
    
    ## mitochondria-related
    combine_mt_rRNA += mt_rRNA_reads
    combine_mt_tRNA += mt_tRNA_reads
    combine_mt_ORF += mt_ORF_reads
    
    lengthToTotal[avaiLength] = totalReads
    lengthTorRNA[avaiLength] = rRNAreads
    lengthToGenes[avaiLength] = potentialGeneReads
    
    ## mitochondria-related
    lengthToCHROMMtotal[avaiLength] = chrMreads
    lengthToCHROMMrRNA[avaiLength] = mt_rRNA_reads
    lengthToCHROMMtRNA[avaiLength] = mt_tRNA_reads
    lengthToCHROMM_ORF[avaiLength] = mt_ORF_reads

print
print '======== combined ====================='
print 'Combined number of reads ', combineReads
print 'Combined number of hits: ', combineHits
print
print 'Combined chrY reads: ', combineYreads
print 'Combined chrM reads: ', combineMreads
print 'Combined potential gene reads: ', combinePotentialGeneReads
print
print 'Combined rRNA reads ', combineRRNA
print 'Combined otherRNA reads: ', combineOTHER
print
print 'Combined mito rRNA reads: ', combine_mt_rRNA
print 'Combined mito tRNA reads: ', combine_mt_tRNA
print 'Combined mito ORF reads: ', combine_mt_ORF

for length in lengthToTotal:
    print length, lengthToTotal[length], lengthTorRNA[length], lengthToGenes[length], lengthToCHROMMtotal[length], lengthToCHROMMrRNA[length], lengthToCHROMMtRNA[length], lengthToCHROMM_ORF[length]

totalReadsList = []
rrnaReadsList = []
genesReadsList = []

## mitochondria-related
chromM_total_readsList = []
chromM_rRNA_readsList = []
chromM_tRNA_readsList = []
chromM_ORF_readsList = []

avaiLengths.sort()  
for avaiLength in avaiLengths:
	totalReadsList.append(lengthToTotal[avaiLength])
	rrnaReadsList.append(lengthTorRNA[avaiLength])
	genesReadsList.append(lengthToGenes[avaiLength])
	
	chromM_total_readsList.append(lengthToCHROMMtotal[avaiLength])
	chromM_rRNA_readsList.append(lengthToCHROMMrRNA[avaiLength])
	chromM_tRNA_readsList.append(lengthToCHROMMtRNA[avaiLength])
	chromM_ORF_readsList.append(lengthToCHROMM_ORF[avaiLength])
	
print avaiLengths
print '%%%%%%%%%%%%%%%%%% TOTAL %%%%%%%%%%%%%%%%%%%%'
print totalReadsList
print
print '%%%%%%%%%%%%%%%%%% rRNA %%%%%%%%%%%%%%%%%%%%'
print rrnaReadsList
print
print '%%%%%%%%%%%%%%%%%% potential genes %%%%%%%%%%%%%%%%%%%%'
print genesReadsList
print
print '%%%%%%%%%%%%%%%%%% chrM TOTAL %%%%%%%%%%%%%%%%%%%%'
print chromM_total_readsList
print 
print '%%%%%%%%%%%%%%%%%% chrM rRNA TOTAL %%%%%%%%%%%%%%%%%%%%'
print chromM_rRNA_readsList
print 
print '%%%%%%%%%%%%%%%%%% chrM tRNA TOTAL %%%%%%%%%%%%%%%%%%%%'
print chromM_tRNA_readsList
print 
print '%%%%%%%%%%%%%%%%%% chrM ORF TOTAL %%%%%%%%%%%%%%%%%%%%'
print chromM_ORF_readsList
print 

    
    

    
    
    
