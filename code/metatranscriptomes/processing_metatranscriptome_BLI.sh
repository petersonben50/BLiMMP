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

# Check to make sure they're done:
cd ~/BLiMMP/dataEdited/metatranscriptomes/workingDirectory_IS
find * -empty | wc -l
# We have 1330 files that are empty.
# Seven of these are sortmernaTempIS files. The jobs they correspond to,
# some appear to have been completed, others just have empty IS and nonIS
# files.
# Need to look at this more closely:
mkdir reads_IS reads_nonIS logs
#mv *_nonRNA_IS.BLI* reads_IS
#mv *_nonRNA_nonIS.BLI* reads_nonIS
#mv *.log logs
# When I tried to do this, it reported that one of the file was busy.
# This file was one of the ones listed as empty. Perhaps they're still
# being written even through the job is done? Let's look:
find * -empty | head
# First one is BLI20_MT_001_nonrRNA_splitFiles_di
# The log and both split files are empty.
cd ~/BLiMMP/reports/metatranscriptomes
ls */*BLI20_MT_001_nonrRNA_splitFiles_di*
cat errs/sortmerna_IS_BLI20_MT_001_nonrRNA_splitFiles_di.err
# Got it. The submission scripts generated an error
# that a folder was already in existance. Two potential
# explanations:
# 1. I started the run, stopped it, then didn't delete all
#    the old files. This is unlikely, I usually reset the
#    the directory when I start anew. But possible.
# 2. I noticed that at times, it seems that jobs were sent
#    back, in that the jobs that were idle would increase,
#    usually only by 2 or so. So, a job may start, reconsider,
#    then stop, but already have started creating files.
#    If this is the case, at the start of the submission
#    file we just need to add a script to delete that folder
#    if it exists.
# Interesting that those folders do not exist now.

# Looks like if things went well, there are no error outputs:
cd ~/BLiMMP/reports/metatranscriptomes/errs
find * -empty | wc -l
find sortmerna_IS_BLI2* -empty | wc -l
# Damn, only 4921.
cd ~/BLiMMP/dataEdited/metatranscriptomes/workingDirectory_IS
wc -l nonrRNA_splitFiles_list.txt
# Out of 5621. That's a lot that failed.

# Let's try again. Reset, new folder:
cd ~/BLiMMP/dataEdited/metatranscriptomes
mv workingDirectory_IS original_workingDirectory_IS
mkdir workingDirectory_IS temp_nonrRNA
mv original_workingDirectory_IS/*fastq temp_nonrRNA
cat original_workingDirectory_IS/nonrRNA_splitFiles_list.txt | while read subset
do
  mv original_workingDirectory_IS/$subset workingDirectory_IS/$subset
done
mv original_workingDirectory_IS/nonrRNA_splitFiles_list.txt workingDirectory_IS/nonrRNA_splitFiles_list.txt
condor_submit /home/GLBRCORG/bpeterson26/BLiMMP/code/sortmerna_featureOfInterest_internalStandard.sub
