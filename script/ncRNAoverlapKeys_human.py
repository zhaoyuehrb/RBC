
## called file, from initial file
## load rRNA loci from ai, Rmsk
## load length file with datasets --> directory, length range tuple

import utilityModule as um
import annotationImporter as ai
import time
###############################################################################
def findrRNA(repeatLocus, rnaType):
    return repeatLocus.repClass() == rnaType
#####################################################################
def getKeys(directory, wantedLength, rnaType):
    print time.asctime()
    rrnaCollection = ai.getRmskRepeats('hg19', lambda repeatLocus: findrRNA(repeatLocus, rnaType))
    print time.asctime()
    print 'Number of loci in '+rnaType+' collection: ', len(rrnaCollection)
    print
    
    length = int(wantedLength)
    lenTuple = (length, length+1)
    datasets = [directory]
    requirementTrue = lambda i: True

    genomic_hits = um.getAllGenomicHits(datasets, lenTuple, requirementTrue, requirementTrue)
    print 'Number of genomic hits loaded: ', len(genomic_hits)
    genomicHits = um.LocusCollection(genomic_hits, 1000)
    print 'Number of hits in genomicHits collection: ', len(genomicHits)
    print '========================================'

    rrna_hits = []
    for locus in rrnaCollection.getLoci():
        hits = genomicHits.getOverlap(locus, sense = 'sense')
        rrna_hits.extend(hits)
    rrnaHitsCollection = um.LocusCollection(rrna_hits, 1000)
    print 'Number of '+rnaType+' hits in collection: ', len(rrnaHitsCollection)
    print '========================================'

    rrnaKeys = []
    readsFromrrna = 0    
    for hit in rrnaHitsCollection.getLoci():
        chrom = hit.chr()
        length = hit.len()
        sense = hit.sense()
        if sense == '+':
            rawStart = hit.start()
        elif sense == '-':
            rawStart = hit.end()
        collatedKey = chrom+'('+sense+'):'+str(rawStart)+'_'+str(length)
        rrnaKeys.append(collatedKey)
        readsFromrrna += hit.readNum()

    print '===== Finished cycling through hitsCollection =========', time.asctime()

##    corefilename = directory[:-1].split('/')[-1]
##    faEntries = um.get_fa(directory+corefilename+'.'+wantedLength+'.fa')
##    print 'Number of entries in .fa file: ', len(faEntries)
##    seqToKey = {}
##    for key in faEntries.keys():
##        seq = faEntries[key].upper()
##        seqToKey[seq] = key
##    print 'Number of entries in seqToKey: ', len(seqToKey)
##    print time.asctime()
##
##    rrnaKeys = []
##    readsFromrrna = 0
##    for hit in rrnaHitsCollection.getLoci():
##        seq = hit.seq()
##        key = seqToKey[seq]
##        rrnaKeys.append(key)
##        readsFromrrna += hit.readNum()
    print 'Number of reads from '+rnaType+' hits: ', readsFromrrna
    return rrnaKeys
    
