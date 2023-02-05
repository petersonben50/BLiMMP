#!/bin/sh


##############################
# code/metatranscriptomes/sortmerna_featureOfInterest.sh
# Benjamin D. Peterson

# This script runs sortmerna on a metagenome,
# using a specific reference database.
##############################

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate sortmerna
PYTHONPATH=""

# Get set up
cd $workingDirectory
echo "Mapping against these databases:"
ls -l $reference

# Clean up any old files
rm -rf sortmernaTempIS_$MT_subset
rm -f nonRNA_$featureName
rm -f nonRNA_non$featureName

# Run job
sortmerna --ref $reference \
          --reads $MT_subset \
          --fastx \
          --paired_in \
          --workdir  sortmernaTempIS_$MT_subset \
          --aligned nonRNA_$featureName \
          --other nonRNA_non$featureName \
          --threads 7

# Clean up temp file
rm -rf sortmernaTempIS_$MT_subset
