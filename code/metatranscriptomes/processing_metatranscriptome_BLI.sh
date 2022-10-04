#!/bin/sh

#########################
# code/metatranscriptomes/processing_metatranscriptomes_BLI.sh
# Benjamin D. Peterson
#########################

# Upload metatranscriptomes_list.txt to ~/BLiMMP/metadata
#########################
# Set up needed directories
#########################
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
mkdir ~/BLiMMP/dataEdited/metatranscriptomes
mkdir ~/BLiMMP/dataEdited/metatranscriptomes/workingDirectory
mkdir ~/BLiMMP/dataEdited/metatranscriptomes/reports
mkdir ~/BLiMMP/reports/metatranscriptomes
cd ~/BLiMMP/reports/metatranscriptomes
mkdir outs errs logs


#########################
# Download fastp
#########################
cd ~/BLiMMP/code
wget http://opengene.org/fastp/fastp
chmod +x ./fastp


#########################
# Trim metatranscriptomes with fastp
# This saved out interleaved files that
# we immediately split up.
#########################
condor_status -avail -long -attributes Name,Cpus,Memory
chmod +x /home/GLBRCORG/bpeterson26/BLiMMP/code/trimming_MT_fastp.sh
condor_submit /home/GLBRCORG/bpeterson26/BLiMMP/code/trimming_MT_fastp.sub

transcriptomeDirectory=/home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/metatranscriptomes
cd $transcriptomeDirectory/workingDirectory
ls *_splitFiles_* > split_file_list.txt
wc -l split_file_list.txt


#########################
# Run sortmerna on metatranscriptomes
#########################
chmod +x /home/GLBRCORG/bpeterson26/BLiMMP/code/sortmerna_rRNA_silva.sh
condor_submit /home/GLBRCORG/bpeterson26/BLiMMP/code/sortmerna_rRNA_silva.sub
