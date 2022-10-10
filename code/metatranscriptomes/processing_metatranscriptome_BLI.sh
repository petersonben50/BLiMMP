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

# Check which ones are held up:
cd ~/BLiMMP/reports/metatranscriptomes/outs
ls sortmerna_subset_*.out | wc -l
# 7719 total, which makes the original number of jobs
tail -q -n 1 sortmerna_subset_*.out | grep -v 'Removing working directory'
# [writeLog:898] Using Log file: "/mnt/bigdata/linuxhome/bpeterson26/BLiMMP/dataEdited/metatranscriptomes/workingDirectory/BLI20_MT_007_splitFiles_hz_rRNA.log"
grep -L 'Removing working directory' sortmerna_subset_*.out
# sortmerna_subset_BLI20_MT_007_splitFiles_hz.out
# sortmerna_subset_BLI21_MT_009_splitFiles_ow.out
condor_rm 321179
# Need to manually run these, there are a bunch of intermediate files left.
# The files of interest look intact, but best just to re-run.
cd ~/BLiMMP/dataEdited/metatranscriptomes/workingDirectory
rm -rf BLI20_MT_007_splitFiles_hz_nonRNA.BLI20_MT_007_splitFiles_hz BLI20_MT_007_splitFiles_hz_rRNA.BLI20_MT_007_splitFiles_hz BLI20_MT_007_splitFiles_hz_rRNA.log sortmernaTemp_MT_subset_BLI20_MT_007_splitFiles_hz
# In screen, set MT_subset to BLI20_MT_007_splitFiles_hz and run the scripts in sortmerna_rRNA_silva.sh
rm -rf  sortmernaTemp_MT_subset_BLI21_MT_009_splitFiles_owsortmerna_keys_2* BLI21_MT_009_splitFiles_ow_nonRNA.BLI21_MT_009_splitFiles_ow BLI21_MT_009_splitFiles_ow_rRNA.BLI21_MT_009_splitFiles_ow sortmernaTemp_MT_subset_BLI21_MT_009_splitFiles_ow
# In screen, set MT_subset to BLI21_MT_009_splitFiles_ow and run the scripts in sortmerna_rRNA_silva.sh



#########################
# Combine and count non-rRNA and rRNA reads from each sample
#########################
transcriptomeDirectory=/home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/metatranscriptomes
mkdir $transcriptomeDirectory/rRNA_reads
mkdir $transcriptomeDirectory/workingDirectory_IS
echo -e 'mtID\trRNA_reads\tnonrRNA_reads' > $transcriptomeDirectory/reports/mt_read_counts_rRNA.tsv
cat /home/GLBRCORG/bpeterson26/BLiMMP/metadata/metatranscriptome_list.txt | while read mtID
do
  echo "Counting rRNA reads for" $mtID
  cat rRNA.$mtID\_splitFiles_* > $transcriptomeDirectory/rRNA_reads/$mtID\_rRNA.fastq
  rRNA_lines=`cat $transcriptomeDirectory/rRNA_reads/$mtID\_rRNA.fastq | wc -l`
  rRNA_counts=`expr $rRNA_lines / 4`

  echo "Counting non-rRNA reads for" $mtID
  cat nonRNA.$mtID\_splitFiles_* > $transcriptomeDirectory/workingDirectory_IS/$mtID\_nonRNA.fastq
  nonrRNA_lines=`cat $transcriptomeDirectory/workingDirectory_IS/$mtID\_nonrRNA.fastq | wc -l`
  nonrRNA_counts=`expr $nonrRNA_lines / 4`

  echo -e "$mtID\t$rRNA_counts\t$nonrRNA_counts" >> $transcriptomeDirectory/reports/mt_read_counts_rRNA.tsv
done
cd ~/BLiMMP/dataEdited/metatranscriptomes
rm -rf workingDirectory


#########################
# Pull out reads mapping to internal standard
#########################
cd ~/BLiMMP/dataEdited/metatranscriptomes/workingDirectory_IS
# Split them up first
linesToCut=800000
cat /home/GLBRCORG/bpeterson26/BLiMMP/metadata/metatranscriptome_list.txt | while read mtID
do
  echo "Splitting up" $mtID
  split -l $linesToCut $mtID\_nonrRNA.fastq $mtID\_nonrRNA_splitFiles_
done
ls *_nonrRNA_splitFiles_* > nonrRNA_splitFiles_list.txt

# Run submission file to identify IS reads
chmod +x /home/GLBRCORG/bpeterson26/BLiMMP/code/sortmerna_featureOfInterest.sh
condor_submit /home/GLBRCORG/bpeterson26/BLiMMP/code/sortmerna_featureOfInterest_internalStandard.sub
