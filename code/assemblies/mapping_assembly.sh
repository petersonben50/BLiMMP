#!/bin/sh

######################
# code/assemblies/mapping_assembly.sh
# Benjamin D. Peterson
######################

############################################
############################################
# Generate cleaned metagenomes for mapping
############################################
############################################
# For mapping purposes, we don't want the paired reads
# to be merged.
mkdir ~/BLiMMP/dataEdited/mapping/
mkdir ~/BLiMMP/dataEdited/mapping/reads_for_mapping
grep 'BLI2' ~/BLiMMP/dataEdited/metagenomes/reports/naming_key.tsv | \
  awk '{ print $2 }' \
  > ~/BLiMMP/dataEdited/metagenomes/reports/metagenome_list.txt
chmod +x /home/GLBRCORG/bpeterson26/BLiMMP/code/mapping_clean_metagenome.sh
condor_submit /home/GLBRCORG/bpeterson26/BLiMMP/code/mapping_clean_metagenome.sub
