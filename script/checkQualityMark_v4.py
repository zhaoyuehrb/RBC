## open fastq file
## for each package, check is last character of name is 1 or 0
## if 1, carry on
## if 0, triage
## check if number of packages that pass is equal to number passing chastity filter

import sys
import os

in_dir = sys.argv[1] #RNAseq or RPF
prefix = sys.argv[2]
out_dir = sys.argv[3]


script_dir=os.path.dirname(__file__)

infileString = in_dir + "/" + prefix+ '_sequence.txt'
infile = open(infileString)

triagefileString = out_dir  + prefix +'/triage_'+prefix+'_sequence.txt'
triagefile = open(triagefileString, 'w+')

outfileString = out_dir + prefix + '/adapTrim_'+prefix+'_sequence.txt'
outfile = open(outfileString, 'w+')

triagePackages = 0
outPackages = 0

lines = 0
packages = 0
passQual = 0
lastLine = False
while lastLine == False:
    package = []
    for i in range(0, 4):
        infileLine = infile.readline()
        if infileLine == '':
            lastLine = True
        else:
            lines += 1
        ## store line in package
        package.append(infileLine)
    if not lastLine:
        package = map(lambda i: i.rstrip('\n'), package)
        name = package[0]
        seq = package[1]
        rename = package[2]
        qual = package[3]
        if package[0].split(';')[1] == '1':
            passQual += 1
            ## if passQual
            ## do the same things as fastqParse did

##            seq = seq[:36]    ## shld not change 36bp read, but will trim 54bp read to 36bp
##            qual = qual[:36]
            ##########################################################################
            ## keep length of read to find TCGTAT in
            ## if TCGTAT not found, read will be 40 nt
            ## --> then, trim it to 36 nt for mapping
            ## --> else, trim to the position where adaptor is found for mapping
            ##########################################################################
            polyA = seq.upper().find('A'*18)
            if polyA == -1:
                # did not find polyA, continue
                adaptorMatch = seq.upper().rfind('TCGTAT')
                #if adaptorMatch == -1:
                # modified 16 September to let sequencing with length 50 behave same as length 40
                if adaptorMatch == -1 or adaptorMatch > 34: 
                    # did not find adaptor, output whole package
                    ###############################
                    seq = seq[:36]
                    qual = qual[:36]
                    ###############################
                    packageList = [name, seq.upper(), rename, qual]
                    outPackage = '\n'.join(packageList)
                    outfile.write(outPackage+'\n')
                    outPackages += 1
                else:
                    # found adaptor, trim
                    seq = seq[:adaptorMatch].upper()
                    qual = qual[:adaptorMatch]
                    if len(seq) >= 18:
                        # at least 18nt after trimming, output trimmed package
                        packageList = [name, seq, rename, qual]
                        outPackage = '\n'.join(packageList)
                        outfile.write(outPackage+'\n')
                        outPackages += 1
                    else:
                        # less than 18nt after trimming, triage
                        packageList = [name, seq, rename, qual]
                        triagePackage = '\n'.join(packageList)
                        triagefile.write(triagePackage+'\n')
                        triagePackages += 1
            else:
                # found polyA, triage
                packageList = [name, seq.upper(), rename, qual]
                triagePackage = '\n'.join(packageList)
                triagefile.write(triagePackage+'\n')
                triagePackages += 1
        packages += 1
infile.close()
triagefile.close()
outfile.close()

print
print in_dir, prefix
print
print 'Number of lines read: ', lines
print 'Number of packages: ', packages
print 'Number of packages that would passQual: ', passQual
print
print 'Number of packages in outfile: ', outPackages
print 'Number of packages in triagefile: ', triagePackages
