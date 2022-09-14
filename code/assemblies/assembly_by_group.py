#################################
# coassembly_by_group.py
# Benjamin D. Peterson

# This script will read in a csv file with the
# group information for coassembling metagenomes
# and start metaSPADes on these groups.
#################################


####---------------------------------####
# Load needed python packages
####---------------------------------####
import os
import sys
import numpy
import argparse
import pandas as pd


####---------------------------------####
# Set up parser
####---------------------------------####
# Set up an argument parser
parser = argparse.ArgumentParser()
parser.add_argument('--groupName')
parser.add_argument('--groupSpreadsheetName')
parser.add_argument('--readLocation')
parser.add_argument('--mergedReads')
parser.add_argument('--output')
# Settings for metaSPADes
parser.add_argument('--threads')
parser.add_argument('--memory')


####---------------------------------####
# Parse input variables
####---------------------------------####
# Parse names from argument
inputs = parser.parse_args()
GROUP_NAME = inputs.groupName
GROUP_SPREADSHEET_NAME = inputs.groupSpreadsheetName
READ_LOCATION = inputs.readLocation
MERGED_READS = inputs.mergedReads
OUTPUT = inputs.output

# Settings for metaSPADes
THREADS = inputs.threads
if THREADS is None:
    THREADS = 6
MEMORY = inputs.memory
if MEMORY is None:
    THREADS = 100
print("Running with " + THREADS + " and " + MEMORY + " Gbs of memory.")

# Check if we're using merged reads
if MERGED_READS == "yes":
    print("Using merged reads")
elif MERGED_READS is None:
    print("No merged reads, since you didn't specify")
else:
    print("No merged reads. If you want merged reads,")
    print("you need to supply '--mergedReads yes'")


####---------------------------------####
# Read in the grouping file
####---------------------------------####
group2mg = pd.read_csv(GROUP_SPREADSHEET_NAME)


####---------------------------------####
# Select the metagenomes that correspond to the selected group
####---------------------------------####
groupMetagenomes = group2mg.loc[group2mg['groupID'] == group]
neededMetagenomes = groupMetagenomes['metagenomeID']
print("Assembling " + group + " using " + neededMetagenomes)


####---------------------------------####
# Loop over the metagenome names, adding them to the metaSPADes command
####---------------------------------####
metaSPadesCommand = 'metaspades.py -t ' + THREADS + ' -m ' + MEMORY + ' -k 21,33,55,77,99,127 '
for metagenome in neededMetagenomes:
    metaSPadesCommand = metaSPadesCommand + '--pe-1 1 ' + readLocation + '/' + metagenome + '_R1.fastq.gz '
    metaSPadesCommand = metaSPadesCommand + '--pe-2 1 ' + readLocation + '/' + metagenome + '_R2.fastq.gz '
    metaSPadesCommand = metaSPadesCommand + '--pe-s 1 ' + readLocation + '/' + metagenome + '_single.fastq.gz '
    if MERGED_READS == 'yes':
        metaSPadesCommand = metaSPadesCommand + '--pe-m 1 ' + readLocation + '/' + metagenome + '_merged.fastq.gz '

metaSPadesCommand = metaSPadesCommand + '-o ' + output
print(metaSPadesCommand)


#######################################
# Run metaSPADes
#######################################
os.system(metaSPadesCommand)
