#!/bin/sh

#########################
# code/metatranscriptomes/pseudomapping_mRNA reads.sh
# Benjamin D. Peterson
#########################


##################################################
##################################################
# Pseudoalign RNA reads to references
##################################################
##################################################
mkdir ~/BLiMMP/dataEdited/metatranscriptomes/alignment
mkdir ~/BLiMMP/dataEdited/metatranscriptomes/alignment/index
mkdir ~/BLiMMP/dataEdited/metatranscriptomes/alignment/reports

# Index nucleic acid ORFs for kallisto
screen -S kallisto_indexing
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
PYTHONPATH=""
referenceDirectory=/home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/assemblies/ORFs
conda activate kallisto
cd ~/BLiMMP/dataEdited/metatranscriptomes/alignment/index
awk -F '\t' '{ print $2 }' ~/BLiMMP/dataEdited/metatranscriptomes/reports/pseudomapping_key.tsv | \
  sort | uniq | \
  while read assemblyID
  do
    echo "Indexing" $assemblyID
    kallisto index -i $assemblyID.idx $referenceDirectory/$assemblyID.fna
  done

# Submit kallisto analyses
chmod +x /home/GLBRCORG/bpeterson26/BLiMMP/code/kallisto_pseudoalignment.sh
condor_submit /home/GLBRCORG/bpeterson26/BLiMMP/code/kallisto_pseudoalignment.sub
