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
referenceDB=/home/GLBRCORG/bpeterson26/references/sortmerna_dbs
originalReadsDirectory=/home/GLBRCORG/bpeterson26/BLiMMP/dataRaw/metatranscriptomes
transcriptomeDirectory=/home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/metatranscriptomes


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
screen -S MT_processing
originalReadsDirectory=/home/GLBRCORG/bpeterson26/BLiMMP/dataRaw/metatranscriptomes
transcriptomeDirectory=/home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/metatranscriptomes
fastp=~/BLiMMP/code/fastp
linesToCut=160000
cd $transcriptomeDirectory/workingDirectory

condor_status -avail -long -attributes Name,Cpus,Memory
chmod +x /home/GLBRCORG/bpeterson26/BLiMMP/code/trimming_MT_fastp.sh
condor_submit /home/GLBRCORG/bpeterson26/BLiMMP/code/trimming_MT_fastp.sub

transcriptomeDirectory=/home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/metatranscriptomes
cd $transcriptomeDirectory/workingDirectory
ls *_splitFiles_* > split_file_list.txt
