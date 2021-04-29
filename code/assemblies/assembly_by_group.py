#################################
# coassembly_by_group_metaspades.py
# Benjamin D. Peterson

# This script will read in a csv file with the
# group information for coassembling metagenomes
# and start metaSPADes on these groups.
#################################


#######################################
# Load needed python packages
#######################################
import os
import sys
import numpy
import pandas as pd


#######################################
# Read in the needed inputs from the command line
#######################################
group = sys.argv[1]
clusterSpreadSheetName = sys.argv[2]
readLocation = sys.argv[3]
output = sys.argv[4]


#######################################
# Read in the grouping file
#######################################
group2metagenome = pd.read_csv(clusterSpreadSheetName)


#######################################
# Select the metagenomes that correspond to the selected group
#######################################
groupMetagenomes = group2metagenome.loc[group2metagenome['groupID'] == group]
neededMetagenomes = groupMetagenomes['metagenomeID']
print("Assembling " + group + " using " + neededMetagenomes)


#######################################
# Loop over the metagenome names, adding them to the metaSPADes command
#######################################
metaSPadesCommand = 'metaspades.py -t 24 -m 1000 -k 21,33,55,77,99,127 '
for metagenome in neededMetagenomes:
    metaSPadesCommand = metaSPadesCommand + '--pe-1 1 ' + readLocation + '/' + metagenome + '_R1.fastq.gz '
    metaSPadesCommand = metaSPadesCommand + '--pe-2 1 ' + readLocation + '/' + metagenome + '_R2.fastq.gz '
    metaSPadesCommand = metaSPadesCommand + '--pe-s 1 ' + readLocation + '/' + metagenome + '_single.fastq.gz '
    metaSPadesCommand = metaSPadesCommand + '--pe-m 1 ' + readLocation + '/' + metagenome + '_merged.fastq.gz '

metaSPadesCommand = metaSPadesCommand + '-o ' + output
print(metaSPadesCommand)


#######################################
# Run metaSPADes
#######################################
os.system(metaSPadesCommand)
