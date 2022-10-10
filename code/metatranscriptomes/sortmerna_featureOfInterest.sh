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

cd $workingDirectory
echo "Mapping against these databases:"
ls -l $reference

sortmerna --ref $reference \
          --reads $MT_subset \
          --fastx \
          --paired_in \
          --workdir  sortmernaTempIS_MT_subset_$MT_subset \
          --aligned $MT_subset\_nonRNA_$featureName \
          --other $MT_subset\_nonRNA_non$featureName \
          --threads 3

rm -rf sortmernaTempIS_MT_subset_$MT_subset
