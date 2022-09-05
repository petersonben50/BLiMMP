#!/bin/sh

#########################
# code/metagenomes/metagenome_retrieval.sh
# Benjamin D. Peterson
#########################


##################################################
##################################################
# Transfer data
##################################################
##################################################

#########################
# Transfer 2020 data
#########################
screen -S BLI20_MG_transfer
mkdir ~/BLiMMP
mkdir ~/BLiMMP/dataRaw
mkdir ~/BLiMMP/dataRaw/metagenomes
mkdir ~/BLiMMP/dataRaw/metagenomes/2020_data
cd ~/BLiMMP/dataRaw/metagenomes/2020_data

#source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
#conda activate lftp
#lftp -c 'set ssl:verify-certificate no set ftp:ssl-protect-data true set ftp:ssl-force true; open -u n210421_150PE_NVS1B_S1_Peterson,Shiezeicie7Eab5 -e "mirror -c; quit" ftp://gslanalyzer.qb3.berkeley.edu:990'
mv BP_* ../
# Check md5sum


#########################
# Transfer 2021 data
#########################
screen -S BLI21_MG_transfer
cd ~/BLiMMP/dataRaw/metagenomes
mkdir 2021_data
cd 2021_data
#source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
#conda activate lftp
#lftp -c 'set ssl:verify-certificate no set ftp:ssl-protect-data true set ftp:ssl-force true; open -u n220708_150PE_NVS1A_S4_Peterson_M001705,teoDiishech4chi -e "mirror -c; quit" ftp://gslanalyzer.qb3.berkeley.edu:990'
# ls -alh
# Stash the L001 files
mkdir 2021_L001
mv *_L001_* 2021_L001

# Compare md5sum values of downloaded files to the ones they sent
#md5sum KMBP011*gz > post_download_md5sum.txt
#cat post_download_md5sum.txt
#cat md5sum.txt
# They match up!
mv KMBP011* ../



##################################################
##################################################
# Generate and upload the naming file
##################################################
##################################################
mkdir ~/BLiMMP/metadata
# Upload a list of the metagenomes:
# ~/BLiMMP/dataEdited/metagenomes/reports/naming_key.tsv
rawReads=~/BLiMMP/dataRaw/metagenomes
ancillary_info=~/BLiMMP/dataEdited/metagenomes/reports

# Rename the metagenomes
grep 'BLI2' $ancillary_info/naming_key.tsv | while read line
do
  sequencingID=`echo $line | awk '{ print $1 }'`
  metagenomeID=`echo $line | awk '{ print $2 }'`
  echo "Moving" $rawReads/$sequencingID\_*R1*fastq.gz "to" $rawReads/$metagenomeID\_R1.fastq.gz
  mv $rawReads/$sequencingID\_*R1*fastq.gz $rawReads/$metagenomeID\_R1.fastq.gz
  mv $rawReads/$sequencingID\_*R2*fastq.gz $rawReads/$metagenomeID\_R2.fastq.gz
done


##################################################
##################################################
# Clean up!
##################################################
##################################################
# Only want our metagenomes of interest in here.
# Download 2020_data and 2021_data folders to local repository: dataRaw/metagenomes
rm 2020_data 2021_data
# Move other people's metagenomes
mkdir ~/collaborators/2020_other_metagenomes
mv BP_* ~/collaborators/2020_other_metagenomes
