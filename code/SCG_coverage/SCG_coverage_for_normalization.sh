#!/bin/sh

######################
# code/SCG_coverage_for_normalization.sh
# Benjamin D. Peterson
######################

######################
# Get set up
######################
mkdir /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/scg_abundance
mkdir /home/GLBRCORG/bpeterson26/BLiMMP/reports/scgAbund
cd /home/GLBRCORG/bpeterson26/BLiMMP/reports/scgAbund
mkdir outs errs logs


######################
# Submit jobs
######################
chmod +x /home/GLBRCORG/bpeterson26/BLiMMP/code/SCG_abundance_in_assemblies.sh
condor_submit /home/GLBRCORG/bpeterson26/BLiMMP/code/SCG_abundance_in_assemblies.sub


######################
# Concatenate files
######################
cd ~/BLiMMP/dataEdited/scg_abundance
cat *_scg_coverage.tsv > scg_coverage.tsv
rm -rf working_directory_*
