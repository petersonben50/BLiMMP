# First source the python viz venv:
# source ~/virtual-envs/py_viz/bin/activate

import os
import sys
import pandas as pd
from Bio import SeqIO

#gffName = 'Bacteroidetes_0005.annot.gff'
#fnaFile = 'Bacteroidetes_0005.fna'
#central_gene = 'Bacteroidetes_0005_113_7'
#distanceIncluded = 20000
#output = 'Bacteroidetes_0005_113_7'

gffName = sys.argv[1]
fnaFile = sys.argv[2]
central_gene = sys.argv[3]
distanceIncluded = sys.argv[4]
distanceIncluded = int(distanceIncluded)
output = sys.argv[5]

listOfNames = ['sequence', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phrase', 'attributes']
df = pd.read_csv(gffName, sep = '\t', names = listOfNames)

def get_element(my_list, position):
    return my_list[position]

def reversed_string(a_string):
    return a_string[::-1]

# Keep only gene sequences from this scaffold
scaffold = df.loc[df['attributes'].str.contains(central_gene),'sequence'].item()
df = df.loc[df['sequence'] == scaffold]

# Read in NA sequence, and set length
for seq_record in SeqIO.parse(fnaFile, "fasta"):
    if seq_record.id == scaffold:
        print(seq_record.id + " matches")
        fastaSequence = str(seq_record.seq)
        fastaLength = len(fastaSequence)

# Calculate start coordinates of gene of interest
df['start'] = df['start'].astype(int)
df['end'] = df['end'].astype(int)
startcoord_ref = int(df.loc[df['attributes'].str.contains(central_gene),'start'])
endcoord_ref = int(df.loc[df['attributes'].str.contains(central_gene),'end'])

# We want to use this to line up sequences in Geneious, so
# if we're looking for more length than we have, we'll fill
# in the sequence with N's.
front_end_padding = (startcoord_ref - distanceIncluded)*-1
if front_end_padding < 0:
    front_end_padding = 0

back_end_padding = ((fastaLength - (endcoord_ref + distanceIncluded))*-1 - 1)
if back_end_padding < 0:
    back_end_padding = 0

# Figure out direction of gene
sign = df.loc[df['attributes'].str.contains(central_gene),'strand'].item()
print(central_gene + " is on " + sign + "strand")

if sign == "+": # Subtract all coords from initial start coord
    strand_key = {"+" : "+", "-" : "-"}
    startVector = (df['start'] - startcoord_ref) + distanceIncluded
    endVector = (df['end'] - startcoord_ref) + distanceIncluded
elif sign == "-": # -1 * Subtract all coords from initial end (true start)
    vector = -1.0
    strand_key = {"-" : "+", "+" : "-"} # Flip strands
    startVector = (vector*(df['end'] - endcoord_ref)) + distanceIncluded
    endVector = (vector*(df['start'] - endcoord_ref)) + distanceIncluded

df['start'] = startVector.apply(int)
df['end'] = endVector.apply(int)
df['strand'] = df['strand'].map(strand_key)


# Filter by columns
df = df.loc[df['start'] > 0]
df = df.loc[df['end'] < (distanceIncluded*2 + (endcoord_ref - startcoord_ref))]


startFastaPosition = startcoord_ref - distanceIncluded - 1
if startFastaPosition < 0:
    startFastaPosition = 0

endFastaPosition = endcoord_ref + distanceIncluded - 1
if endFastaPosition > len(fastaSequence):
    endFastaPosition = len(fastaSequence)

fastaSequence = fastaSequence[startFastaPosition:endFastaPosition]

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
if sign == "-":
    listBases = list(fastaSequence)
    reversedBases = [complement[base] for base in listBases]
    reversedBases = ''.join(reversedBases)
    fastaSequence = reversed_string(reversedBases)
    fastaSequence = ("N" * back_end_padding) + fastaSequence + ("N" * front_end_padding)

if sign == "+":
    fastaSequence = ("N" * front_end_padding) + fastaSequence + ("N" * back_end_padding)


# Write out files
gffOutput = output + "_neighborhood.gff"
df.to_csv(gffOutput, sep = '\t', index = False, header = False)

fastaOutput = output + "_neighborhood.fna"
with open(fastaOutput, 'w') as outFile:
    outFile.write('>' + scaffold + '\n' + fastaSequence + '\n')
