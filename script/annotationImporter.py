# this module defines functions for importing annotation sets.
# these functions may contain information for the reading of
# specific datasets, including directory paths.


import utilityModule, os


class RepeatLocus(utilityModule.Locus):
    # singleton model for names
    __classifyD = dict()
    def __init__(self,locus,rName,rClass,rFamily,rStart,rEnd):
        utilityModule.Locus.__init__(self,locus.chr(),locus.start(),locus.end(),locus.sense())
        if not(self.__classifyD.has_key(rName)): self.__classifyD[rName] = rName
        if not(self.__classifyD.has_key(rClass)): self.__classifyD[rClass] = rClass
        if not(self.__classifyD.has_key(rFamily)): self.__classifyD[rFamily] = rFamily
        self.__repName = self.__classifyD[rName]
        self.__repClass = self.__classifyD[rClass]
        self.__repFamily = self.__classifyD[rFamily]
        self.__repStart = rStart
        self.__repEnd = rEnd
        self.checkRep()
    def checkRep(self):
        utilityModule.Locus.checkRep(self)
        # repStart must be a smaller number than repEnd
        if self.__repStart > self.__repEnd:
            raise ValueError("repStart ("+str(self.__repStart)+") is too big for repEnd ("+str(self.__repEnd)+").")
    def repName(self): return self.__repName
    def repClass(self): return self.__repClass
    def repFamily(self): return self.__repFamily
    def repStart(self): return self.__repStart
    def repEnd(self): return self.__repEnd
    def repCoords(self): return [self.__repStart,self.__repEnd]


# myFilter is a lambda that will allow me to limit the set of repeat loci returned
# to a particular subset of loci.  It must take a single arguement which is an instance
# of RepeatLocus.
# the information from each column of the input files is summarized below.
# 0: ?
# 1: ?
# 2: ?
# 3: ?
# 4: ?
# 5: chr name
# 6: start coord
# 7: end coord
# 8: ?
# 9: strand
# 10: repeat name
# 11: repeat class
# 12: repeat family
# 13: rep. start
# 14: rep. end
# 15: ?
# 16: ?
# assembly is like 'mm8', this is a check to make sure that the correct
#   annotations are being used.
def getRmskRepeats(assembly,myFilter = lambda rl: True, winSize = 10000):
    if assembly=='mm8': repeatDirectory = '/lab/bartel/ruby/Solexa/databases/mm8/annotations/Rmsk/'
    #elif assembly=='hg18': repeatDirectory = '/lab/bartel1_ata/ruby/Solexa/databases/hg18/Annotation/Rmsk/'
    elif assembly=='hg19': repeatDirectory = '/home/hguo/Documents/annotations/hg19/Rmsk/'
    elif assembly=='ce4': repeatDirectory = '/lab/bartel/ruby/Solexa/databases/ce4/Annotations/Rmsk/'
    elif assembly == 'mm9': repeatDirectory = '/home/hguo/Documents/annotations/mm9/Rmsk/'
    else: raise ValueError("assembly '"+assembly+"' isn't valid.")
    allRmskFiles = map(lambda f: repeatDirectory+f, os.listdir(repeatDirectory))
    
    repeatCollection = utilityModule.LocusCollection([],winSize)
    for rmskFilename in allRmskFiles:
        print rmskFilename
        rmskFile = open(rmskFilename)
        lineNum = 0
        for line in rmskFile.readlines():
            lineNum += 1
            if line[-1]=='\n': line = line[:-1]
            cols = line.split()
            # check for the correct number of columns
            if len(cols)!=17: raise ValueError("Line "+str(lineNum)+" in "+rmskFilename+" is incorrect:\n'"+line+"'")
            localLocus = utilityModule.Locus(cols[5],int(cols[6]),(int(cols[7]) - 1),cols[9])
            localRepeat = RepeatLocus(localLocus,cols[10],cols[11],cols[12],int(cols[13]),int(cols[14]))
            if myFilter(localRepeat): repeatCollection.append(localRepeat)
        rmskFile.close()

    return repeatCollection


# this forwards to getSangerGenes so that i don't have to fix old scripts
def getElegansCodingGenes(genomeName,test = lambda g: True):
    return getSangerCodingGenes(genomeName,test)


# this gets all of the annotated protein-coding genes for nematodes
# the test lambda takes one Gene instance as an argument
def getSangerCodingGenes(genomeName,test = lambda g: True):
    genomeToFile = dict()
    genomeToFile['ce2'] = '/lab/bartel/ruby/454/C_elegans/genome/sangerGene.txt'
    genomeToFile['ce4'] = '/lab/bartel/ruby/Solexa/databases/ce4/Annotations/sangerGene.txt'
#    genomeToFile['caePb1'] = 
#    genomeToFile['caeRem2'] = 
#    genomeToFile['cb3'] = 
#    genomeToFile['priPac1'] = 
    filename = genomeToFile[genomeName]
    geneList = []
    f = open(filename)
    for line in f.readlines():
        raw_line = line.split('\t')
        name = raw_line[0]
        chr = raw_line[1]
        sense = raw_line[2]
        txCoords = map(int,raw_line[3:5])
        cdCoords = map(int,raw_line[5:7])
        exStarts = map(int,raw_line[8].split(',')[:-1])
        exEnds = map(int,raw_line[9].split(',')[:-1])
        # fix the edges to work with the Gene specs
        exEnds = map(lambda n: n-1, exEnds)
        txCoords = [txCoords[0], txCoords[1] - 1]
        cdCoords = [cdCoords[0], cdCoords[1] - 1]
        if cdCoords[0]!=0 or cdCoords[1]!=0:
            newGene = utilityModule.Gene(name,chr,sense,txCoords,cdCoords,exStarts,exEnds)
            if test(newGene): geneList.append(newGene)
    f.close()
    return geneList



# this forwards to getSangerGenes so that i don't have to fix old scripts
def getElegansNcGenes(genomeName,test = lambda g: True):
    return getSangerNcGenes(genomeName,test)


# this gets all of the annotated non-protein-coding genes for ce2
# the test lambda takes one Gene instance as an argument
def getSangerNcGenes(genomeName,test = lambda g: True):
    genomeToFile = dict()
    genomeToFile['ce2'] = '/lab/bartel/ruby/454/C_elegans/genome/sangerGene.txt'
    genomeToFile['ce4'] = '/lab/bartel/ruby/Solexa/databases/ce4/Annotations/sangerGene.txt'
    filename = genomeToFile[genomeName]
    geneList = []
    f = open(filename)
    for line in f.readlines():
        raw_line = line.split('\t')
        name = raw_line[0]
        chr = raw_line[1]
        sense = raw_line[2]
        txCoords = map(int,raw_line[3:5])
        cdCoords = map(int,raw_line[5:7])
        exStarts = map(int,raw_line[8].split(',')[:-1])
        exEnds = map(int,raw_line[9].split(',')[:-1])
        # fix the edges to work with the Gene specs
        exEnds = map(lambda n: n-1, exEnds)
        txCoords = [txCoords[0], txCoords[1] - 1]
        cdCoords = [cdCoords[0], cdCoords[1] - 1]
        if cdCoords[0]==0 and cdCoords[1]==0:
            newGene = utilityModule.Gene(name,chr,sense,txCoords,cdCoords,exStarts,exEnds)
            if test(newGene): geneList.append(newGene)
    f.close()
    return geneList



# mm8 refGenes
def getMouseRefGenes(test = lambda g: True):
    filename = '/lab/bartel/ruby/Solexa/databases/mm8/annotations/refGene.txt'
    geneList = []
    f = open(filename)
    for line in f.readlines():
        raw_line = line.split('\t')
        name = raw_line[1]
        chr = raw_line[2]
        sense = raw_line[3]
        txCoords = map(int,raw_line[4:6])
        cdCoords = map(int,raw_line[6:8])
        exStarts = map(int,raw_line[9].split(',')[:-1])
        exEnds = map(int,raw_line[10].split(',')[:-1])
        # fix the edges to work with the Gene specs
        exEnds = map(lambda n: n-1, exEnds)
        txCoords = [txCoords[0], txCoords[1] - 1]
        cdCoords = [cdCoords[0], cdCoords[1] - 1]
        if cdCoords[0]!=0 or cdCoords[1]!=0:
            newGene = utilityModule.Gene(name,chr,sense,txCoords,cdCoords,exStarts,exEnds)
            if test(newGene): geneList.append(newGene)
    f.close()
    return geneList




# gets the list of genes corresponding to each mountain from Kim et al 2001
# (C. elegans gene sets).  returns a dictionary whose keys are the names of
# gene collections and whose values are lists of the genes in that collection.
def getKimMountains():
    kimDir = '/lab/bartel/ruby/Solexa/databases/KimMountains/'
    allMountainFiles = map(lambda f: kimDir+f,
                           filter(lambda f: f.find('README')==-1, os.listdir(kimDir)))
    mountToGenes = dict()
    for mfName in allMountainFiles:
        mount = mfName.split('/')[-1].split('.')[0]
        genes = []
        f = open(mfName)
        for i in f.readlines():
            i = i.split('#')[0].strip()
            genes.append(i)
        mountToGenes[mount] = genes
    return mountToGenes
            
