#!/bin/sh

##############################
# code/me/anvioDB_generation.sh
# Benjamin D. Peterson

#
##############################

echo "Working in the sortmerna conda environment"
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate sortmerna
PYTHONPATH=""
cd $workingDirectory

echo ""
echo "Identifying rRNA reads in" $MT_subset
echo ""
echo "Mapping against these databases:"
ls -l $referenceDB

sortmerna --ref $referenceDB/rfam-5.8s-database-id98.fasta \
              --ref $referenceDB/rfam-5s-database-id98.fasta \
              --ref $referenceDB/silva-arc-16s-id95.fasta \
              --ref $referenceDB/silva-arc-23s-id98.fasta \
              --ref $referenceDB/silva-bac-16s-id90.fasta \
              --ref $referenceDB/silva-bac-23s-id98.fasta \
              --ref $referenceDB/silva-euk-18s-id95.fasta \
              --ref $referenceDB/silva-euk-28s-id98.fasta \
          --reads $MT_subset \
          --fastx \
          --paired_in \
          --workdir  sortmernaTemp_MT_subset_$MT_subset \
          --aligned rRNA \
          --other nonRNA \
          --threads 3

echo ""
echo "Removing working directory"
rm -rf sortmernaTemp_MT_subset_$MT_subset
