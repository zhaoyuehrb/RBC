
# this is a module with functions for dealing with sequence data in general.
# it was originally built from compoents of 454/DM_fff/fofifo.py.

import string, os


############################ DATA ABSTRACTIONS ##############################


    

    

# species

class Species:
    # seq can be either a string or a DenseSeq instance
    def __init__(self,seq,readNum):
        self._seq = seq
        self._readNum = readNum
    def seq(self): return self._seq
    def readNum(self): return self._readNum
    def __hash__(self): return self._seq.__hash__()



# genomic species

class GenomicSpecies:
    __uninitiatedString = 'uninitiated'
    def __init__(self,species,loci):
        self._species = species
        # make loci non-redundant
        tempLoci = dict()
        for i in loci: tempLoci[i] = None
        self._hits = map(lambda i: GenomicHit(self,i), tempLoci.keys())
    def seq(self): return self._species.seq()
    def denseSeq(self): return self._species.denseSeq()
    def species(self): return self._species
    def readNum(self): return self._species.readNum()
    def hits(self):
        if self._hits==self.__uninitiatedString: self._hits = map(lambda i: GenomicHit(self,i), self._loci)
        return map(lambda i: i, self._hits)
    def numLoci(self): return len(self._hits)
    
    def __hash__(self): return self._species.__hash__()
    




# locus

class Locus:
    # this may save some space by reducing the number of chromosome strings
    # that are associated with Locus instances (see __init__).
    __chrDict = dict()
    __senseDict = {'+':'+', '-':'-', '.':'.'}
    # chr = chromosome name (string)
    # sense = '+' or '-' (or '.' for an ambidexterous locus)
    # start,end = ints of the start and end coords of the locus;
    #      end coord is the coord of the last nucleotide.
    def __init__(self,chr,start,end,sense):
        coords = [start,end]
        coords.sort()
        # this method for assigning chromosome should help avoid storage of
        # redundant strings.
        if not(self.__chrDict.has_key(chr)): self.__chrDict[chr] = chr
        self._chr = self.__chrDict[chr]
        self._sense = self.__senseDict[sense]
        self._start = int(coords[0])
        self._end = int(coords[1])
    def chr(self): return self._chr
    def start(self): return self._start  ## returns the smallest coordinate
    def end(self): return self._end   ## returns the biggest coordinate
    def len(self): return self._end - self._start + 1
    def getAntisenseLocus(self):
        if self._sense=='.': return self
        else:
            switch = {'+':'-', '-':'+'}
            return Locus(self._chr,self._start,self._end,switch[self._sense])
    def coords(self): return [self._start,self._end]  ## returns a sorted list of the coordinates
    def sense(self): return self._sense
    # returns boolean; True if two loci share any coordinates in common
    def overlaps(self,otherLocus):
        if self.chr()!=otherLocus.chr(): return False
        elif not(self._sense=='.' or \
                 otherLocus.sense()=='.' or \
                 self.sense()==otherLocus.sense()): return False
        elif self.start() > otherLocus.end() or otherLocus.start() > self.end(): return False
        else: return True
        
    # returns boolean; True if all the nucleotides of the given locus overlap
    #      with the self locus
    def contains(self,otherLocus):
        if self.chr()!=otherLocus.chr(): return False
        elif not(self._sense=='.' or \
                 otherLocus.sense()=='.' or \
                 self.sense()==otherLocus.sense()): return False
        elif self.start() > otherLocus.start() or otherLocus.end() > self.end(): return False
        else: return True
        
    # same as overlaps, but considers the opposite strand
    def overlapsAntisense(self,otherLocus):
        return self.getAntisenseLocus().overlaps(otherLocus)
    # same as contains, but considers the opposite strand
    def containsAntisense(self,otherLocus):
        return self.getAntisenseLocus().contains(otherLocus)
    def __hash__(self): return self._start + self._end
    def __eq__(self,other):
#        if self.__class__ != other.__class__: return False
        if self.chr()!=other.chr(): return False
        if self.start()!=other.start(): return False
        if self.end()!=other.end(): return False
        if self.sense()!=other.sense(): return False
        return True
    def __ne__(self,other): return not(self.__eq__(other))
    def __str__(self): return self.chr()+'('+self.sense()+'):'+'-'.join(map(str,self.coords()))
    def checkRep(self):
        pass



# genomic hit

class GenomicHit(Locus):
    def __init__(self,genomicSpecies,locus):
        Locus.__init__(self,locus.chr(),locus.start(),locus.end(),locus.sense())
        self._gSpecies = genomicSpecies
    # in case the non-normalized number of reads is required, 'False' can be provided
    # as an argument.
    def readNum(self, normalized=True):
        if normalized:
            return float(self._gSpecies.readNum()) / self._gSpecies.numLoci()
        else:
            return float(self._gSpecies.readNum())

    # gets the fraction of reads that should be counted from this hit as overlapping the
    # provided locus.  ie if half the hit overlaps the locus, then only half of the reads
    # are counted.  NOTE: THIS IS A SENSE-INDEPENDENT DEFINITION OF OVERLAP!
    def overlapReads(self,locus,normalized=True):
        if locus.contains(self) or locus.containsAntisense(self): return self.readNum(normalized)
        elif not(locus.overlaps(self) or locus.overlapsAntisense(self)): return 0.0
        else:
            overlapLen = min([locus.end(),self.end()]) - max([locus.start(),self.start()]) + 1
            readFraction = float(overlapLen) / (self.end() - self.start() + 1)
            return self.readNum(normalized) * readFraction
    def seq(self): return self._gSpecies.seq()
    def fiveEndLocus(self):
        if self.sense()=='+': return Locus(self.chr(),self.start(),self.start(),self.sense())
        elif self.sense()=='-': return Locus(self.chr(),self.end(),self.end(),self.sense())
        else: raise ValueError("sense of '"+str(self.sense())+"' is not a valid sense entry for a hit.")
    def threeEndLocus(self):
        if self.sense()=='+': return Locus(self.chr(),self.end(),self.end(),self.sense())
        elif self.sense()=='-': return Locus(self.chr(),self.start(),self.start(),self.sense())
        else: raise ValueError("sense of '"+str(self.sense())+"' is not a valid sense entry for a hit.")
    def getGenomicSpecies(self): return self._gSpecies





# this will allow me to optimize searching/extraction for speed and abstract
# away optimized organization

class LocusCollection:
    def __init__(self,loci,windowSize):
        ### top-level keys are chr, then strand, no space
        self.__chrToCoordToLoci = dict()
        self.__loci = dict()
        self.__winSize = windowSize
        for lcs in loci: self.__addLocus(lcs)

    def __addLocus(self,lcs):
        if not(self.__loci.has_key(lcs)):
            self.__loci[lcs] = None
            if lcs.sense()=='.': chrKeyList = [lcs.chr()+'+', lcs.chr()+'-']
            else: chrKeyList = [lcs.chr()+lcs.sense()]
            for chrKey in chrKeyList:
                if not(self.__chrToCoordToLoci.has_key(chrKey)): self.__chrToCoordToLoci[chrKey] = dict()
                for n in self.__getKeyRange(lcs):
                    if not(self.__chrToCoordToLoci[chrKey].has_key(n)): self.__chrToCoordToLoci[chrKey][n] = []
                    self.__chrToCoordToLoci[chrKey][n].append(lcs)

    def __getKeyRange(self,locus):
        start = locus.start() / self.__winSize
        end = locus.end() / self.__winSize + 1 ## add 1 because of the range
        return range(start,end)

    def __len__(self): return len(self.__loci)
        
    def append(self,new): self.__addLocus(new)
    def extend(self,newList):
        for lcs in newList: self.__addLocus(lcs)
    def hasLocus(self,locus):
        return self.__loci.has_key(locus)
    def remove(self,old):
        if not(self.__loci.has_key(old)): raise ValueError("requested locus isn't in collection")
        del self.__loci[old]
        if old.sense()=='.': senseList = ['+','-']
        else: senseList = [old.sense()]
        for k in self.__getKeyRange(old):
            for sense in senseList:
                self.__chrToCoordToLoci[old.chr()+sense][k].remove(old)

    def getWindowSize(self): return self.__winSize
    def getLoci(self): return self.__loci.keys()
    def getChrList(self):
        # i need to remove the strand info from the chromosome keys and make
        # them non-redundant.
        tempKeys = dict()
        for k in self.__chrToCoordToLoci.keys(): tempKeys[k[:-1]] = None
        return tempKeys.keys()
            
    def __subsetHelper(self,locus,sense):
        sense = sense.lower()
        if ['sense','antisense','both'].count(sense)!=1:
            raise ValueError("sense command invalid: '"+sense+"'.")
        matches = dict()
        senses = ['+','-']
        if locus.sense()=='.' or sense=='both': lamb = lambda s: True
        elif sense=='sense': lamb = lambda s: s==locus.sense()
        elif sense=='antisense': lamb = lambda s: s!=locus.sense()
        else: raise ValueError("sense value was inappropriate: '"+sense+"'.")
        for s in filter(lamb, senses):
            chrKey = locus.chr()+s
            if self.__chrToCoordToLoci.has_key(chrKey):
                for n in self.__getKeyRange(locus):
                    if self.__chrToCoordToLoci[chrKey].has_key(n):
                        for lcs in self.__chrToCoordToLoci[chrKey][n]:
                            matches[lcs] = None
        return matches.keys()
        
    # sense can be 'sense' (default), 'antisense', or 'both'
    # returns all members of the collection that overlap the locus
    def getOverlap(self,locus,sense='sense'):
        matches = self.__subsetHelper(locus,sense)
        ### now, get rid of the ones that don't really overlap
        realMatches = dict()
        if sense=='sense' or sense=='both':
            for i in filter(lambda lcs: lcs.overlaps(locus), matches):
                realMatches[i] = None
        if sense=='antisense' or sense=='both':
            for i in filter(lambda lcs: lcs.overlapsAntisense(locus), matches):
                realMatches[i] = None 
        return realMatches.keys()

    # sense can be 'sense' (default), 'antisense', or 'both'
    # returns all members of the collection that are contained by the locus
    def getContained(self,locus,sense='sense'):
        matches = self.__subsetHelper(locus,sense)
        ### now, get rid of the ones that don't really overlap
        realMatches = dict()
        if sense=='sense' or sense=='both':
            for i in filter(lambda lcs: locus.contains(lcs), matches):
                realMatches[i] = None
        if sense=='antisense' or sense=='both':
            for i in filter(lambda lcs: locus.containsAntisense(lcs), matches):
                realMatches[i] = None
        return realMatches.keys()

    # sense can be 'sense' (default), 'antisense', or 'both'
    # returns all members of the collection that contain the locus
    def getContainers(self,locus,sense='sense'):
        matches = self.__subsetHelper(locus,sense)
        ### now, get rid of the ones that don't really overlap
        realMatches = dict()
        if sense=='sense' or sense=='both':
            for i in filter(lambda lcs: lcs.contains(locus), matches):
                realMatches[i] = None
        if sense=='antisense' or sense=='both':
            for i in filter(lambda lcs: lcs.containsAntisense(locus), matches):
                realMatches[i] = None
        return realMatches.keys()




# this is a gene object.  unlike the previous gene_object, this actually represents
# a gene, as opposed to a whole set of genes, which was a poor design in the first place.
class Gene:
    # name = name of the gene (string)
    # txCoords = list of coords defining the boundaries of the transcipt
    # cdCoords = list of coords defining the beginning and end of the coding region
    # exStarts = list of coords marking the beginning of each exon
    # exEnds = list of coords marking the end of each exon
    # IF THIS IS A NON-CODING GENE, cdCoords => [0,0]
#    def __init__(self,name,chr,sense,txCoords,cdCoords,exStarts,exEnds):
#        self._name = name
#        self._txLocus = Locus(chr,min(txCoords),max(txCoords),sense)
#        self._cdLocus = Locus(chr,min(cdCoords),max(cdCoords),sense)
#
#        exStarts = map(lambda i: i, exStarts)
#        exEnds = map(lambda i: i, exEnds)
#        exStarts.sort()
#        exEnds.sort()
#        
#        self._txExons = []
#        self._cdExons = []
#        self._introns = []
#        
#        for n in range(len(exStarts)):
#            if n==0:
#                self._txExons.append(Locus(chr,txCoords[0],exEnds[n]-1,sense))
#                self._cdExons.append(Locus(chr,cdCoords[0],exEnds[n]-1,sense))
#            elif n==len(exStarts)-1:
#                self._txExons.append(Locus(chr,txCoords[0],txCoords[1],sense))
#                self._cdExons.append(Locus(chr,cdCoords[0],cdCoords[1],sense))
#            else:
#                newExon = Locus(chr,exStarts[n],exEnds[n]-1,sense)
#                self._txExons.append(newExon)
#                self._cdExons.append(newExon)
#            if n < len(exStarts)-1: self._introns.append(Locus(chr,exEnds[n],exStarts[n+1]-1,sense))
#
#        if sense=='+':
#            self._fpUtr = Locus(chr,txCoords[0],cdCoords[0]-1,sense)
#            self._tpUtr = Locus(chr,cdCoords[1]+1,txCoords[1],sense)
#        elif sense=='-':
#            self._fpUtr = Locus(chr,cdCoords[1]+1,txCoords[1],sense)
#            self._tpUtr = Locus(chr,txCoords[0],cdCoords[0]-1,sense)
    def __init__(self,name,accNum,chr,sense,txCoords,cdCoords,exStarts,exEnds):
         self._name = name
         self._accNum = accNum
         self._txLocus = Locus(chr,min(txCoords),max(txCoords) - 1,sense)
         if cdCoords == None:
             self._cdLocus = None
         else:
             self._cdLocus = Locus(chr,min(cdCoords),max(cdCoords) - 1,sense)

         exStarts = map(lambda i: i, exStarts)
         exEnds = map(lambda i: i, exEnds)
         exStarts.sort()
         exEnds.sort()

         self._txExons = []
         self._cdExons = []
         self._introns = []

         cd_exon_count = 0

         for n in range(len(exStarts)):
             first_locus = Locus(chr,exStarts[n],exStarts[n],sense)
             second_locus = Locus(chr,exEnds[n]-1,exEnds[n]-1,sense)

             # Add the transcription unit exon
             tx_exon = Locus(chr,exStarts[n],exEnds[n]-1,sense)

             self._txExons.append(tx_exon)

             # Add Coding Exons
             # Need to make sure that the current exon is actually in the coding region of the gene first
             if self.isCoding() and tx_exon.overlaps(self._cdLocus):
                 if not first_locus.overlaps(self._cdLocus):
                     first_coord = min(cdCoords)
                 else:
                     first_coord = exStarts[n]

                 if not second_locus.overlaps(self._cdLocus):
                     second_coord = max(cdCoords) - 1
                 else:
                     second_coord = exEnds[n]-1

                 new_cd_exon = Locus(chr,first_coord,second_coord,sense)
                 self._cdExons.append(new_cd_exon)

             # Add Introns
             if n < len(exStarts)-1:
                 self._introns.append(Locus(chr,exEnds[n]+1-1,exStarts[n +1]-1,sense))


##         if sense == '+':
##             if min(txCoords) != min(cdCoords):
##                 # txStart not equal to cdStart, calculate fpUTR the usual way
##                 self._fpUTR = Locus(chr,min(txCoords),min(cdCoords)-1,sense)
##             else:
##                 # txStart = cdStart, no fpUTR
##                 self._fpUTR = None
##
##             if max(txCoords) != max(cdCoords):
##                 # txEnd not equal to cdEnd, calc tpUTR the usual way
##                 self._tpUTR = Locus(chr,max(cdCoords)+1-1,max(txCoords)-1,sense)
##             else:
##                 # txEnd = cdEnd, no tpUTR
##                 self._tpUTR = None
##                 
##         elif sense == '-':
##             if max(txCoords) != max(cdCoords):
##                 # txStart not equal to cdStart, calc fpUTR the usual way
##                 self._fpUTR = Locus(chr,max(cdCoords)+1-1,max(txCoords)-1,sense)
##             else:
##                 # txStart = cdStart, no fpUTR
##                 self._fpUTR = None
##
##             if min(txCoords) != min(cdCoords):
##                 # txEnd not equal to cdEnd, calc tpUTR the usual way
##                 self._tpUTR = Locus(chr,min(txCoords),min(cdCoords)-1,sense)
##             else:
##                 # txEnd = cdEnd, no tpUTR
##                 self._tpUTR = None
                 
                 
                 
                 
                 

         if self.isCoding():
             if sense=='+':
##                 if min(txCoords) != min(cdCoords):
##                     self._fpUTR = Locus(chr,min(txCoords),min(cdCoords)-1,sense)
##                 else:
##                     self._fpUTR = Locus(chr,min(txCoords),min(txCoords)-1, sense)
##                 if max(txCoords) != max(cdCoords):
##                     self._tpUTR = Locus(chr,max(cdCoords)+1-1,max(txCoords)-1,sense)
##                 else:
##                     self._tpUTR = Locus(chr,max(txCoords)-1,max(txCoords)-1-1,sense)
                 self._fpUTR = Locus(chr,min(txCoords),min(cdCoords)-1,sense)
                 self._tpUTR = Locus(chr,max(cdCoords)+1-1,max(txCoords)-1,sense)
             elif sense=='-':
                 
##                 if max(txCoords) != max(cdCoords):
##                     self._fpUTR = Locus(chr,max(cdCoords)+1-1,max(txCoords)-1,sense)
##                 else:
##                     self._fpUTR = Locus(chr,max(txCoords)-1,max(txCoords)-1-1,sense)
##                 if min(txCoords) != min(cdCoords):
##                     self._tpUTR = Locus(chr,min(txCoords),min(cdCoords)-1,sense)
##                 else:
##                     self._tpUTR = Locus(chr,min(txCoords),min(txCoords)-1,sense)
                 self._fpUTR = Locus(chr,max(cdCoords)+1-1,max(txCoords)-1,sense)
                 self._tpUTR = Locus(chr,min(txCoords),min(cdCoords)-1,sense)
           
            
         else:
             self._fpUTR = None
             self._tpUTR = None
            

    def name(self): return self._name
    def accNum(self): return self._accNum
    def chr(self): return self._txLocus.chr()
    def sense(self): return self._txLocus.sense()
    def txLocus(self): return self._txLocus   ## locus of full transcript
    def cdLocus(self): return self._cdLocus   ## locus from start codon to end codon
    def txExons(self): return map(lambda i: i, self._txExons)  ## list of loci
    def cdExons(self): return map(lambda i: i, self._cdExons)  ## list of loci
    def introns(self): return map(lambda i: i, self._introns)  ## list of loci
    def fpUtr(self): return self._fpUTR  ## locus
    def tpUtr(self): return self._tpUTR  ## locus
    def isCoding(self): return not(self._cdLocus.start()==0 and self._cdLocus.end()==0)  # boolean; is this gene protein-coding?
    def __hash__(self): return self._txLocus.__hash__()






# getEnvironment
# returns the environment of the indicated file after interpretation
# by the python interpreter
def getEnvironment(filename):
    sandbox = dict()
    file = open(filename)
    text = file.read()
    file.close()
    exec text in sandbox
    return sandbox
    



# getChrSize
# ------------------------------------------------------------------------------
# args: nib_file_with_path: the name of the .nib file for the chromosome whose size
#            you want to determine, complete with path
# returns: chrsize: the integer length of the chromosome (# nucs)
# ------------------------------------------------------------------------------
def getChrSize(nib_file_with_path):
    sri,size_report = os.popen4('nibFrag '+nib_file_with_path+' 0 10000000000000 + stdout')
    sri.close()
    chrsize = int(size_report.read().split(' ')[7][:-1])
    size_report.close()
    return chrsize




# get_fa
# ------------------------------------------------------------------------------
# args: filename: name of fasta sequence file
#       onlyKeys: (optional) a list of keys; only get the entries for these keys.
# returns: d: dictionary whose keys are title lines of .fa file entries, values
#               are sequences
# ------------------------------------------------------------------------------
# this function opens a .fa file and returns a dictionary of its contents
def get_fa(filename,onlyKeys = 'All'):
    keyD = dict()
    if onlyKeys!='All':
        for k in onlyKeys: keyD[k] = None
        
    def addEntry(current,d):
        entry = current.split('\n')
        if entry[0]!='' and (onlyKeys=='All' or keyD.has_key(entry[0])):
            d[entry[0]] = ''.join(entry[1:])

    d = dict()
    bf = open(filename)
    current = ['\n']
    for line in bf:
        if line[0]=='>':
            addEntry(''.join(current),d)
            if len(line)>1: current = [line[1:]]
            else: current = ['\n']
        else: current.append(line)
    addEntry(''.join(current),d)
    bf.close()
    return d


# getFolds
# ------------------------------------------------------------------------------
# args: sl: a list of sequences
# returns: ind_folds: a list of bracket-notation folds corresponding to the mfe
#               fold for the sequence of the same index in the sl list
# ------------------------------------------------------------------------------
# uses RNAfold (zuker algorithm) to generate an mfe fold for each of the sequences
# in sl
def getFolds(sl):
    numr = range(len(sl))
    ind_folds = []    
    all_s = '\n'.join(sl)
    fi,fo = os.popen2('RNAfold')
    fi.write(all_s)
    fi.close()
    b = fo.read().split('\n')
    fo.close()
    missed = 0
    for n in numr:
        if b[n*2]=='':
            missed+=1
            c = ['.'*len(sl[n])]
        else: c = b[n*2+1-missed].split(' ')
        ind_folds.append(c[0])
    return ind_folds

def getFold(seq):
    return getFolds([seq])[0]


# getSequence
# ------------------------------------------------------------------------------
# args: chr: a string with a chromosome name, like 'chr4' or 'chrUn_random'
#       start: an integer coordinate for the start position of the fetched sequence
#       finish: an integer coordinate for the end position of the fetched sequence
#       side: a string '+' or a '-', indicating strand to fetch from
#       f: a directory pathway to the nibfile of interest
# returns:
#       a string which is the fetched sequence
# ------------------------------------------------------------------------------
# used by get_queries to nibFrag the sequences to be scored.  uses jim kent's nibFrag
# to fetch sequences from a .nib indexed sequence file.
def getSequence(locus,f):
    if f[-1]=='/': f = f[:-1]
    chr = locus.chr()
    start = locus.start()
    finish = locus.end() + 1
    side = locus.sense()
    a = os.popen('nibFrag '+f+'/'+chr+'.nib '+str(start)+' '+str(finish)+' '+str(side)+' stdout','r')
    b = a.read().split('\n')
    a.close()
    return ''.join(b[1:]).upper()



# reverseComp
# ------------------------------------------------------------------------------
# args: st: a string
# returns: the reverse compliment of the string
# ------------------------------------------------------------------------------
# makes reverse compliment of a string
def reverseComp(st):
    comp = string.maketrans('ATCG','TAGC')
    return reverseString(st).translate(comp)


# reverseString
# ------------------------------------------------------------------------------
# args: st: a string
# returns: the same string with all characters appearing in the reverse order
# ------------------------------------------------------------------------------
# reverses the order of characters in a string
def reverseString(st):
    li = []
    for i in st: li.append(i)
    li.reverse()
    return ''.join(li)


# gets the indicated percentile value from the indicated numList
def getPercentile(numList,percentile):
    if len(numList)==0:
        return 'NA'
    else:
        spot = len(numList) * float(percentile)
        baseBin = int(spot)
        fraction = 1.0 - (spot % 1)
        numList = map(lambda i: i, numList)
        numList.sort()
        numList.append(numList[-1])
        return numList[baseBin]*fraction + numList[baseBin+1]*(1.0 - fraction)


################### DEALING WITH SOLEXA PROCESSING STRUCTURES ####################


# makeGenomicHitsCollection
# this is the function that is used over and over to get genomic hits
# in the form of a LocusCollection.
# ------------------------------------------------------------------------------
# listOfDatasets: a list of strings, each string being the complete
#      pathway to a directory that contains a dataset (with .solexa file,
#      .fa and .bl files, and a masterFile.py).
# lenRange: a pair of integers indicating the shortest and longest+1 read
#      length that will be returned.
# winsize: an integer indicating the window size for the LocusCollection
# requirementGS: a lambda that takes a GenomicSpecies as an arg and returns
#      a boolean.  decides which genomicSpeceis will be included.
# requirementGH: a lambda that takes a GenomicHit as an arg and returns
#      a boolean.  decides which genomicHits will be included.
# ------------------------------------------------------------------------------
def makeGenomicHitsCollection(listOfDatasets,lenRange,winsize,requirementGS,requirementGH):
    genomicSpecies = getAllGenomicSpecies(listOfDatasets,lenRange,requirementGS)
    genomicHits = makeGenomicHits(genomicSpecies,requirementGH)
    ghc = LocusCollection(genomicHits,winsize)
    return ghc

# or if you just want the hits as a list
def getAllGenomicHits(listOfDatasets,lenRange,requirementGS,requirementGH):
    genomicSpecies = getAllGenomicSpecies(listOfDatasets,lenRange,requirementGS)
    genomicHits = makeGenomicHits(genomicSpecies,requirementGH)
    return genomicHits



# parseBlasts
# ------------------------------------------------------------------------------
# args: filename: name of fasta sequence file
# returns: a list of genomicSpecies
#               are lists of loci (instances of Locus)
# ------------------------------------------------------------------------------
# this function opens a .bl file and returns a dictionary of its contents
def parseBlasts(filename,requirement):
    # assumes that the query file has the same name, but with '.fa' instead of '.bl'
    hits = dict()
    lastLine = ''  ## i will use this to check that the last line indicates the
                   ## .bl file was completed.
    bf = open(filename)
    for b in bf.readlines(): ## iterates through lines
        b = b.rstrip()
        lastLine = b
        if len(b)>0 and b[0]!='#':
            entry = b.split('\t')
            if len(entry)<3: print entry, filename
            #for n in [2]: entry[n] = float(entry[n])  # if you want to use scores, make 10 and 11 floats as well
            for n in range(3,10):
                if entry[n]!='NA': entry[n] = int(entry[n])
            if entry[8]<entry[9]: strand = '+'
            else: strand = '-'
            if requirement(entry):
                newLocus = Locus(entry[1],entry[8],entry[9],strand)
                if not(hits.has_key(entry[0])): hits[entry[0]] = []
                hits[entry[0]].append(newLocus)

    bf.close()
    if lastLine.find('# blasting completed')!=0:
        raise ValueError("blasting was not completed; last line = '"+lastLine+"'")

    # this should save space by only collecting fasta entries for sequences that
    # actually mapped to the genome.
    queries = get_fa(filename[:-3]+'.fa', hits.keys())

    gsList = []
    for k in hits.keys():
        sp = parseSpecies(k,queries[k])
        gs = GenomicSpecies(sp,hits[k])
        gsList.append(gs)
        
        ## empty the source dictionaries as their contents is used
        ## to save memory.
        del hits[k]
        del queries[k]
        
    return gsList

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### parseBlasts
### ------------------------------------------------------------------------------
### args: filename: name of fasta sequence file
### returns: a list of genomicSpecies
###               are lists of loci (instances of Locus)
### ------------------------------------------------------------------------------
### this function opens a .bl file and returns a dictionary of its contents
##def parseBlasts(filename,requirement):
##    # assumes that the query file has the same name, but with '.fa' instead of '.bl'
##    hits = dict()
##    lastLine = ''  ## i will use this to check that the last line indicates the
##                   ## .bl file was completed.
##    bf = open(filename)
##    for b in bf.readlines(): ## iterates through lines
##        b = b.rstrip()
##        lastLine = b
##        if len(b)>0 and b[0]!='#':
##            entry = b[:-1].split('\t')
##            if len(entry)<3: print entry, filename
##            for n in [2]: entry[n] = float(entry[n])  # if you want to use scores, make 10 and 11 floats as well
##            for n in range(3,10):
##                if entry[n]!='NA': entry[n] = int(entry[n])
##            if entry[8]<entry[9]: strand = '+'
##            else: strand = '-'
##            newLocus = Locus(entry[1],entry[8],entry[9],strand)
##            if not(hits.has_key(entry[0])): hits[entry[0]] = []
##            hits[entry[0]].append(newLocus)
##
##    bf.close()
##    if lastLine.find('# blasting completed')!=0:
##        raise ValueError("blasting was not completed; last line = '"+lastLine+"'")
##
##    # this should save space by only collecting fasta entries for sequences that
##    # actually mapped to the genome.
##    queries = get_fa(filename[:-3]+'.fa', hits.keys())
##
##    gsList = []
##    for k in hits.keys():
##        sp = parseSpecies(k,queries[k])
##        gs = GenomicSpecies(sp,hits[k])
##        if requirement(gs): gsList.append(gs)
##        
##        ## empty the source dictionaries as their contents is used
##        ## to save memory.
##        del hits[k]
##        del queries[k]
##        
##    return gsList

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



# import the contents of the masterFile
def getMaster(directory):
    if directory[-1]!='/': directory += '/'
    masterFileName = directory + 'masterFile.py'
    return getEnvironment(masterFileName)



# gets all the species from a dataset
# requires: there is only one fasta entry for any given sequence (the directory contains
# only one dataset).
def getSpecies(directory,lenRange):
    if directory[-1]!='/': directory += '/'

    # get the fasta query files
    seqs = dict()
    fastaFilenames = getBlastQueryFiles(directory,lenRange)
    for sf in fastaFilenames: seqs.update(get_fa(sf))

    speciesList = []
    for k in seqs.keys():
        speciesList.append(parseSpecies(k,seqs[k]))
        
    return speciesList




def parseSpecies(key,seq):
    # the number of reads is stored in the key, after the underscore
    readNum = int(key.split('_')[1])
    return Species(seq,readNum)
    



# gets all of the species from a list of directories,
# combining read counts from all represented datasets
def getAllSpecies(directoryList,lenRange):
    seqToSpecies = dict()
    for directory in directoryList:
        for species in getSpecies(directory,lenRange):
            seq = species.seq()
            if seqToSpecies.has_key(seq):
                seqToSpecies[seq] = Species(seq, species.readNum() + seqToSpecies[seq].readNum())
            else:
                seqToSpecies[seq] = species
    return seqToSpecies.values()
            

# this function gets the primary query files (one file per directory per query length)
def getBlastQueryFiles(directoryName, rangeTuple=False):
    # import the contents of the masterFile
    master = getMaster(directoryName)

    # get the names of all the query files
    coreFilename = '.'.join(master['solexaFile'].split('.')[:-1])
    coreFilename = coreFilename.split('/')[-1]
    allFilenames = os.listdir(directoryName)

    queryFiles = []
    for fn in allFilenames:
        fnSplit = fn.split('.')
        # since I don't know that the input will have the appropriate length or
        # that the portion of the name in the length position will resolve to
        # an integer, I put those checks within a try statement, with addition to
        # the query list at the end.
        try:
            ### first, check the name format (corename . integer . fa)
            if fnSplit[-1]=='fa' and \
               '.'.join(fnSplit[:-2])==coreFilename and \
               int(fnSplit[-2]) > 0:
                ### if appropriate, check the rangeTuple
                if not(rangeTuple) or rangeTuple[0] <= int(fnSplit[-2]) < rangeTuple[1]:
                    queryFiles.append(fn)
        except:
            pass

    queryFiles = map(lambda i: directoryName+i, queryFiles)
    return queryFiles


# this function gets the blast result files from a directory without getting
# any of the supporting blast result files.
def getBlastResultFiles(directoryName,lenRange):
    allFilenames = os.listdir(directoryName)
    queryFiles = getBlastQueryFiles(directoryName,lenRange)
    blastResultFiles = []
    for qf in queryFiles:
        noDirQf = qf.split('/')[-1]
        noDirRf = noDirQf[:-3]+'.bl'
        if allFilenames.count(noDirRf)==1: blastResultFiles.append(qf[:-3]+'.bl')
    return blastResultFiles


# gets all of the blast results from a directory.  returns a dictionary
# whose keys are sequences with blast matches to the genome, and values
# are lists of loci.
def getGenomicSpecies(directory,lenRange,requirement):
    if directory[-1]!='/': directory += '/'
    # raises an exception if blasting is not complete for a directory
    ### directoryBlastValidate(directory,lenRange)
    ### i put this test in the parseBlast function, which will
    ### raise an exception if any of the files aren't there or aren't done.

    gsList = []
    # get the names of all the blast files
    blastFilenames = getBlastResultFiles(directory,lenRange)
    for bf in blastFilenames: gsList.extend(parseBlasts(bf,requirement))
    
    return gsList



# gets all of the Genomic species from a list of directories,
# combining read counts from all represented datasets.
# note that some genomic species may have >50 hits associated with them
# if multiple directories are provided.
#
# I need to filter things out that are useless as i load them.
# I therefore allow a lambda to be provided as an
# optional arg.  the lambda must take a GenomicSpecies as its only
# arg and return a boolean (True -> keep the gs instance).
# NOTE: the requirements will be imposed after EACH dataset is loaded,
# so readNum (which may increase for a particular abstract species
# as more datasets are loaded) should NOT be considered.
#
# note that the number of hits corresponding to a particular sequence here is underspecified;
# the max. number of hits returned for a sequence may be >500 if multiple directories are provided
#
def getAllGenomicSpecies(directoryList,lenRange,requirement = lambda gs: True):
    seqToGs = dict()
    for directory in directoryList:
        localGS = getGenomicSpecies(directory,lenRange,requirement)
        for gs in localGS:
            seq = gs.seq()
            if seqToGs.has_key(seq):
                readNum = gs.readNum() + seqToGs[seq].readNum()
                loci = gs.hits()
                loci.extend(seqToGs[seq].hits())
                seqToGs[seq] = GenomicSpecies(Species(seq, readNum), loci)
            else:
                seqToGs[seq] = gs
        print directory, len(seqToGs)
    return seqToGs.values()






def getGenomicHits(directory,lenRange):
    gsList = getGenomicSpecies(directory,lenRange)
    hitList = []
    for gs in gsList: hitList.extend(gs.hits())
    return hitList


def makeGenomicHits(gsList,requirement = lambda gh: True):
    hitList = []
    for gs in gsList:
        for gh in gs.hits():
            if requirement(gh): hitList.append(gh)
    return hitList






# raises an exception if blasting is not complete for a directory
def directoryBlastValidate(directoryName,lenRange):
    # import the contents of the masterFile
    master = getMaster(directoryName)
    # get the names of all the query files
    allFilenames = os.listdir(directoryName)
    queryFiles = getBlastQueryFiles(directoryName,lenRange)

    allDone = []
    notDone = []
    for qf in queryFiles:
        qfResultFile = qf[:-3]+'.bl'
        if blastIsComplete(qfResultFile): allDone.append(qf)
        else: notDone.append(qf)

    if len(notDone)!=0: raise ValueError("Blasting of the following files is not complete:\n" + \
                                         '\n'.join(notDone))




# figures out if the blast job for the given query file is complet
# requires: resultFilename has a path to find the file from the current directory
def blastIsComplete(resultFilename):
    try:
        lastLine = ''
        blFile = open(resultFilename)
        for i in blFile: lastLine = i
        blFile.close()
        if lastLine.find('# blasting completed')==0: return True
        else: return False
    except IOError:
        return False






# makeGenes
# ------------------------------------------------------------------------------
# args: filename: string; name of the file from which genes will be extracted
# returns: list of Gene objects
# ------------------------------------------------------------------------------
def makeGenes(geneFileName):    
    output = []
    f = open(geneFileName)
    for line in f.readlines():
        raw_line = line.split('\t')
        name = raw_line[0]
        accNum = raw_line[1]
        chr = raw_line[2]
        sense = raw_line[3]
        txCoords = map(int,raw_line[4:6])
        cdCoords = map(int,raw_line[6:8])
        exStarts = map(int, filter(lambda i: i!='', raw_line[9].split(',')))
        exEnds = map(int, filter(lambda i: i!='' and i!= '\n', raw_line[10].split(',')))

        new_gene = Gene(name,accNum,chr,sense,txCoords,cdCoords,exStarts,exEnds)
        # this check is because drosophila has six genes whose strand is designated '.';
        # they are all variants of CG32491.
        if new_gene.sense()=='+' or new_gene.sense()=='-': output.append(new_gene)
        
    f.close()
    return output


# makeNoncodingGenes
# ------------------------------------------------------------------------------
# args: filename: string; name of the file from which genes will be extracted
# returns: list of Gene objects
# ------------------------------------------------------------------------------
def makeNoncodingGenes(directory):
    if directory[-1]!='/': directory += '/'
    master = getMaster(directory)
    
    output = []
    
    f = open(master['noncodingGeneFile'])
    for line in f.readlines():
        raw_line = line.split('\t')
        name = raw_line[4]
        chr = raw_line[1]
        sense = raw_line[6]
        txCoords = map(int,raw_line[2:4])
        cdCoords = [0,0]
        exStarts = map(int, filter(lambda i: i!='', raw_line[2].split(',')))
        exEnds = map(int, filter(lambda i: i!='', raw_line[3].split(',')))

        new_gene = Gene(name,chr,sense,txCoords,cdCoords,exStarts,exEnds)
        # this check is because drosophila has six genes whose strand is designated '.';
        # they are all variants of CG32491.
        if new_gene.sense()=='+' or new_gene.sense()=='-': output.append(new_gene)

    f.close()
    return output



# gets the nibDirectory for a series of datasets, and confirms that
# the nibDirectory is common to the coordinates that have been assigned
# for all of the indicated datasets.
def getNibDirectory(directoryList):
    if len(directoryList)==0:
        raise ValueError("directoryList must contain at least one directory to be meaningful.")
    nibDirDict = dict()
    for d in directoryList:
        master = getMaster(d)
        if not(master.has_key('nibDirectory')):
            raise ValueError("Directory '"+d+"' doesn't have a nibDirectory")
        nibDirDict[master['nibDirectory']] = None
    nibDirs = nibDirDict.keys()
    if len(nibDirs)!=1:
        raise ValueError("directories refer to multiple nibDirectories: "+str(nibDirs))
    return nibDirs[0]
        




# RNA(seq) will convert DNA sequence seq to RNA by changing the T's to U's
def RNA(seq):
    tr = string.maketrans('Tt','Uu')
    return seq.translate(tr,'')

# DNA(seq) will convert RNA sequence seq to DNA by changing the U's to T's
def DNA(seq):
    tr = string.maketrans('Uu','Tt')
    return seq.translate(tr,'')




