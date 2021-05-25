

# First source the python viz venv:
# source ~/virtual-envs/py_viz/bin/activate

import os
import sys
import pandas as pd
from Bio import SeqIO

geneListName = sys.argv[1]
gffName = sys.argv[2]
output = sys.argv[3]

geneListName = 'identification/hgcA_good.txt'
gffName = 'scaffolds/hgcA_scaffolds.gff'
output = 'hgcB/downstream_gene_list.txt'

# Read in the GFF file
listOfNames = ['sequence', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phrase', 'attributes']
df = pd.read_csv(gffName, sep = '\t', names = listOfNames)

# Find names
with open(output, 'w') as outFile:
    for idRawNewline in open(geneListName).readlines():
        idRaw = idRawNewline.rsplit("\n")[0]
        scaffold = idRaw.rsplit("_", 1)[0]
        dfScaffold = df[df['sequence'].str.match(scaffold)]
        print("Working on " + scaffold)
        # Find the gene ID and orientation of hgcA
        geneNumber = idRaw.split("_")[3]
        geneID = idRaw.split("_")[2].lstrip('0') + "_" + geneNumber + ";"
        strandvalue = dfScaffold.loc[dfScaffold['attributes'].str.contains(geneID), 'strand'].item()
        # Retrieve downstream gene ID
        if strandvalue == "-":
            print(idRaw + " is on the reverse strand")
            hgcbLocation = int(geneNumber) - 1
        else:
            print(idRaw + " is on the forward strand")
            hgcbLocation = int(geneNumber) + 1
        # Get downstream gene name
        hgcBgeneID = scaffold + "_" + str(hgcbLocation)
        print(hgcBgeneID)
        outFile.write(hgcBgeneID + '\n')
