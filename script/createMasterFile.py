import sys

#solexaFile = '/mnt/31bba345-2359-4c15-97af-4bf5db60b075/zhaoy/preprocessing/processedData/TTAGGC-s_7_1_seed29_genome.solexa'
#nibDirectory = '/home/hguo/Downloads/Mus_musculus/UCSC/mm9/Sequence/Chromosomes/nibDirectory/'

directory = sys.argv[1]
barcode = directory.split('/')[-2]
solexaFile = directory + barcode + '_seed29_genome.solexa'
nibDirectory = '/home/hguo/Downloads/Mus_musculus/UCSC/mm9/Sequence/Chromosomes/nibDirectory/'

with open(directory+'masterFile.py', 'w+') as the_file:
    the_file.write('solexaFile = \''+solexaFile+'\'\n')
    the_file.write('nibDirectory = \''+nibDirectory+'\'\n')
    
