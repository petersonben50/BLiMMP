
#####################################
# metagenomeAssemblyAnalysis/calculate_depth_contigs.py
# Benjamin D. Peterson

# This reads in a depth file (generated by
# the depth function in samtools) and
# averages the depth of each contig.

# It requires only that you input the
# depth file.
#####################################

#####################################
# Load up libraries
#####################################
import os
import sys
import pandas as pd

#####################################
# Read in input from command line
#####################################
depth = sys.argv[1]
filteredLength = int(sys.argv[2])
outputFileName = sys.argv[3]

#####################################
# Open up depth file
#####################################
depthTable = pd.read_table(depth, names = ['contigs', 'locus', 'depth'])

#####################################
# Filter out residues at start of contig
#####################################
depthTable = depthTable[depthTable['locus'] >= filteredLength]


#####################################
# Filter out residues at end of contig
#####################################

# First, find the max contig length
lengthOfContig = depthTable[['contigs', 'locus']].groupby('contigs').max()
lengthOfContig.rename(columns = {'locus':'lengthOfContig'}, inplace = True)

# Substract the length to filter out
lengthOfContig['lengthOfContig'] = lengthOfContig['lengthOfContig'] - filteredLength

# Then, join max contig length DF with depth table
depthTable = pd.merge(depthTable, lengthOfContig, on='contigs', how='outer')

# Filter out end of contig
depthTable = depthTable[depthTable['locus'] <= depthTable['lengthOfContig']]

depthTable = depthTable[['contigs', 'depth']]

#####################################
# Aggregate depth by contig
#####################################
# Find average coverage across contig
depthTable = depthTable.groupby('contigs').mean()
# Read out data
depthTable.to_csv(outputFileName, sep='\t', header = False)