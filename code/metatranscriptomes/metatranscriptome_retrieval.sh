#!/bin/sh

#########################
# code/metatranscriptomes/metatranscriptome_retrieval.sh
# Benjamin D. Peterson
#########################


##################################################
##################################################
# Transfer data
##################################################
##################################################

#########################
# Transfer 2021 data
#########################
screen -S BLI21_MT_transfer
mkdir ~/BLiMMP/dataRaw/metatranscriptomes
cd ~/BLiMMP/dataRaw/metatranscriptomes

#source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate lftp
lftp -c 'set ssl:verify-certificate no set ftp:ssl-protect-data true set ftp:ssl-force true; open -u n220930_150PE_NVS1A_S4_L1_Peterson,ui8Phooquushae7 -e "mirror -c; quit" ftp://gslanalyzer.qb3.berkeley.edu:990'

# Compare md5sum values of downloaded files to the ones they sent
mkdir reports
md5sum KMBP020*gz > reports/md5sum_post_download.txt
mv md5sum.txt laneBarcode.html reports
cd reports
cat md5sum_post_download.txt
cat md5sum.txt
# They match up!


##################################################
##################################################
# Generate and upload the naming file
##################################################
##################################################
# Upload MT_naming_key.tsv to the reports file
# ~/BLiMMP/dataRaw/metatranscriptomes/reports/naming_key.tsv
# This naming_key.tsv file is a tab-separated file with two columns.
# The first has the QB3 name. The second has our metatranscriptomes ID.
rawReads=~/BLiMMP/dataRaw/metatranscriptomes
ancillary_info=~/BLiMMP/dataRaw/metatranscriptomes/reports

# Rename the metagenomes
grep 'BLI2' $ancillary_info/MT_naming_key.tsv | while read line
do
  sequencingID=`echo $line | awk '{ print $1 }'`
  metagenomeID=`echo $line | awk '{ print $2 }'`
  echo "Moving" $rawReads/$sequencingID\_*R1*fastq.gz "to" $rawReads/$metagenomeID\_R1.fastq.gz
  mv $rawReads/$sequencingID\_*R1*fastq.gz $rawReads/$metagenomeID\_R1.fastq.gz
  mv $rawReads/$sequencingID\_*R2*fastq.gz $rawReads/$metagenomeID\_R2.fastq.gz
done

# Clean up
rm -rf reports
