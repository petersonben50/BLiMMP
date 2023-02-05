#!/bin/sh


##############################
# code/assemblies/ORF_prediction_assemblies.sh
# Benjamin D. Peterson

# This script uses Prodigal to predict
# the open reading frames of the assemblies
# from the BLiMMP metagenomes.
# Submission file loops over all the assemblies.
##############################


##############################
# Set up the environment
##############################
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
cd $workingDirectory


##############################
# Run ORF prediction
##############################
if [ -e scaffolds/$assembly\_assembly.fna ]; then
  echo $assembly "has been cleaned. Let's predict ORFs."
  if [ -e ORFs/$assembly.faa ]; then
    echo "You already predicted ORFs for" $assembly". Relax."
  else
    echo "Time to predict some proteins for" $assembly
    prodigal -i scaffolds/$assembly\_assembly.fna \
              -o ORFs/$assembly.gff \
              -f gff \
              -a ORFs/$assembly.faa \
              -d ORFs/$assembly.fna \
              -p meta

    python /home/GLBRCORG/bpeterson26/BLiMMP/code/cleanFASTA.py ORFs/$assembly.fna
    mv -f ORFs/$assembly.fna_temp.fasta ORFs/$assembly.fna
    python /home/GLBRCORG/bpeterson26/BLiMMP/code/cleanFASTA.py ORFs/$assembly.faa
    mv -f ORFs/$assembly.faa_temp.fasta ORFs/$assembly.faa
  fi
else
  echo $assembly "has not yet been cleaned (possibly not assembled)."
fi
