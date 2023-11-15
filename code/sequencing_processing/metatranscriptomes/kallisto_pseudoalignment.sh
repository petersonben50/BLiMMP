#!/bin/sh

##############################
# code/executables/kallisto_pseudoalignment.sh
# Benjamin D. Peterson

# This is a script that will pseudoalign
# reads from a transcriptome to a specified
# set of ORFs that have been indexed.
##############################

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
PYTHONPATH=""
conda activate kallisto

cd $outputDirectory

echo "Pseudoaligning" $mtID "to" $assemblyID
kallisto quant \
          -i $referenceDirectory/$assemblyID.idx \
          -o $mtID\_to_$assemblyID\_output \
          -b 100 \
          -t 7 \
          $readsLocation/$mtID\_mRNA_R1.fastq \
          $readsLocation/$mtID\_mRNA_R2.fastq

# Gather and rename abundance files
cp $mtID\_to_$assemblyID\_output/abundance.tsv ./$mtID\_to_$assemblyID\_kallisto.tsv
cp $mtID\_to_$assemblyID\_output/run_info.json reports/$mtID\_to_$assemblyID\_run_info.json
rm -rf $mtID\_to_$assemblyID\_output
