#!/bin/sh

#########################
# code/metagenomes/metagenome_processing.sh
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


##################################################
##################################################
# Trim metagenomes
##################################################
##################################################
screen -S BLiMMP_metagenome_trimming
mkdir ~/BLiMMP/dataEdited
mkdir ~/BLiMMP/dataEdited/metagenomes
mkdir ~/BLiMMP/dataEdited/metagenomes/reports
# Upload naming_key.tsv to the reports file
# This naming_key.tsv file is a tab-separated file with two columns.
# The first has the QB3 name. The second has our metagenome ID.
# Data was originally found in our dilution table spreadsheets:
# 2020 data: /Users/benjaminpeterson/Documents/research/BLiMMP/dataEdited/dnaSequencing/2020/samplePrep/KMBP010_dilutions.xlsx
# 2021 data: /Users/benjaminpeterson/Documents/research/BLiMMP/dataEdited/dnaSequencing/2021/samplePrep/BLI21_MG_dilutions.xlsx

# Download fastp
#mkdir ~/BLiMMP/code
#cd ~/BLiMMP/code
#wget http://opengene.org/fastp/fastp
#chmod a+x ./fastp
fastp=~/BLiMMP/code/fastp

cd ~/BLiMMP/dataRaw/metagenomes
read_storage=~/BLiMMP/dataEdited/metagenomes
ancillary_info=~/BLiMMP/dataEdited/metagenomes/reports
rawReads=~/BLiMMP/dataRaw/metagenomes

echo -e '\n\n' >> $ancillary_info/metagenome_run_info.txt
date >> $ancillary_info/metagenome_run_info.txt
grep 'BLI2' $ancillary_info/naming_key.tsv | while read line
do
  sequencingID=`echo $line | awk '{ print $1 }'`
  metagenomeID=`echo $line | awk '{ print $2 }'`
  if [ ! -e $read_storage/$metagenomeID\_R1.fastq.gz ]; then
    echo "Processing" >> $ancillary_info/metagenome_run_info.txt
    ls $rawReads/$sequencingID\_*R1*fastq.gz >> $ancillary_info/metagenome_run_info.txt
    ls $rawReads/$sequencingID\_*R2*fastq.gz >> $ancillary_info/metagenome_run_info.txt
    echo "into" $metagenomeID >> $ancillary_info/metagenome_run_info.txt
    $fastp --in1 $rawReads/$sequencingID\_*R1*fastq.gz \
            --in2 $rawReads/$sequencingID\_*R2*fastq.gz \
            --out1 $read_storage/$metagenomeID\_R1.fastq.gz \
            --out2 $read_storage/$metagenomeID\_R2.fastq.gz \
            --unpaired1 $read_storage/$metagenomeID\_single.fastq.gz \
            --unpaired2 $read_storage/$metagenomeID\_single.fastq.gz \
            --merge \
            --merged_out $read_storage/$metagenomeID\_merged.fastq.gz \
            --failed_out $ancillary_info/$metagenomeID\_failed.fastq.gz \
            --html $ancillary_info/$metagenomeID\_report.html \
            --detect_adapter_for_pe \
            --cut_tail \
            --cut_tail_window_size 10 \
            --cut_tail_mean_quality 20 \
            --length_required 100
  else
    echo "Already processed" $metagenomeID >> $ancillary_info/metagenome_run_info.txt
  fi
done





############################################
############################################
# Check size of metagenomes
############################################
############################################

# Count reads in metagenome pre-trimming
screen -S BLiMMP_metagenome_read_counting
cd ~/BLiMMP/dataRaw/metagenomes
rawReads=~/BLiMMP/dataRaw/metagenomes
ancillary_info=~/BLiMMP/dataEdited/metagenomes/reports
echo -e "metagenomeID\tforwardReads\treverseReads" > $ancillary_info/metagenome_read_count_pre_trimming.tsv
grep 'BLI2' $ancillary_info/naming_key.tsv | while read line
do
  sequencingID=`echo $line | awk '{ print $1 }'`
  metagenomeID=`echo $line | awk '{ print $2 }'`
  echo "Working on" $metagenomeID
  forwardCount=$(zgrep -c "^@" $rawReads/$sequencingID\_*R1*fastq.gz)
  reverseCount=$(zgrep -c "^@" $rawReads/$sequencingID\_*R1*fastq.gz)
  echo -e $metagenomeID"\t"$forwardCount"\t"$reverseCount >> $ancillary_info/metagenome_read_count_pre_trimming.tsv
done


# Count reads in metagenome post-trimming
screen -S BLiMMP_metagenome_read_counting_post
cd ~/BLiMMP/dataRaw/metagenomes
read_storage=~/BLiMMP/dataEdited/metagenomes
ancillary_info=~/BLiMMP/dataEdited/metagenomes/reports
echo -e "metagenomeID\tforwardReads\treverseReads\tsingleReads\tmergedReads" > $ancillary_info/metagenome_read_count.tsv
grep 'BLI2' $ancillary_info/naming_key.tsv | while read line
do
  metagenomeID=`echo $line | awk '{ print $2 }'`
  echo "Working on" $metagenomeID
  forwardCount=$(zgrep -c "^@" $read_storage/$metagenomeID\_R1.fastq.gz)
  reverseCount=$(zgrep -c "^@" $read_storage/$metagenomeID\_R2.fastq.gz)
  singleCount=$(zgrep -c "^@" $read_storage/$metagenomeID\_single.fastq.gz)
  mergeCount=$(zgrep -c "^@" $read_storage/$metagenomeID\_merged.fastq.gz)
  echo -e $metagenomeID"\t"$forwardCount"\t"$reverseCount"\t"$singleCount"\t"$mergeCount >> $ancillary_info/metagenome_read_count.tsv
done


# Coverage pre-trimming
screen -S EG_MG_coverage_counting_pre
cd ~/BLiMMP/dataEdited/metagenomes/reports
code=~/BLiMMP/code/readfq-master
read_storage=~/BLiMMP/dataRaw/metagenomes
ancillary_info=~/BLiMMP/dataEdited/metagenomes/reports
echo -e "sequencingID\tmetagenomeID\tR1\tR2" > metagenome_coverage_pre_trimming.tsv
grep 'BLI2' $ancillary_info/naming_key.tsv | while read line
do
  sequencingID=`echo $line | awk '{ print $1 }'`
  metagenomeID=`echo $line | awk '{ print $2 }'`
  echo "Counting coverage in" $sequencingID
  R1_count=$($code/kseq_fastq_base $read_storage/$sequencingID*_R1*.fastq.gz | \
                awk -F " " '{ print $5 }')
  R2_count=$($code/kseq_fastq_base $read_storage/$sequencingID*_R2*.fastq.gz | \
                awk -F " " '{ print $5 }')
  echo -e $sequencingID"\t"$metagenomeID"\t"$R1_count"\t"$R2_count >> metagenome_coverage_pre_trimming.tsv
done


# Coverage post-trimming
screen -S EG_MG_coverage_counting_post
cd ~/BLiMMP/dataEdited/metagenomes/reports
code=~/BLiMMP/code/readfq-master
read_storage=~/BLiMMP/dataEdited/metagenomes
ancillary_info=~/BLiMMP/dataEdited/metagenomes/reports
IFS=$'\n'
echo -e "metagenomeID\tR1\tR2\tsingle\tmerged" > metagenome_coverage.tsv
grep 'BLI2' $ancillary_info/naming_key.tsv | while read line
do
  metagenomeID=`echo $line | awk '{ print $2 }'`
  echo "Counting coverage in" $metagenomeID
  R1_count=$($code/kseq_fastq_base $read_storage/$metagenomeID\_R1.fastq.gz | \
                awk -F " " '{ print $5 }')
  R2_count=$($code/kseq_fastq_base $read_storage/$metagenomeID\_R2.fastq.gz | \
                awk -F " " '{ print $5 }')
  single_count=$($code/kseq_fastq_base $read_storage/$metagenomeID\_single.fastq.gz | \
                awk -F " " '{ print $5 }')
  merged_count=$($code/kseq_fastq_base $read_storage/$metagenomeID\_merged.fastq.gz | \
                awk -F " " '{ print $5 }')
  echo -e $metagenomeID"\t"$R1_count"\t"$R2_count"\t"$single_count"\t"$merged_count >> metagenome_coverage.tsv
done




############################################
############################################
# Run Mash on metagenomes
############################################
############################################

screen -S BLI_mash
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics


# Generate sketches for each metagenome

mkdir ~/BLiMMP/dataEdited/mash_data
mkdir ~/BLiMMP/dataEdited/mash_data/temp_MG_files
mkdir ~/BLiMMP/dataEdited/mash_data/sketch_files
cd ~/BLiMMP/dataEdited/mash_data
read_storage=~/BLiMMP/dataEdited/metagenomes


cat ~/HellsCanyon/metadata/lists/MG_list_all.txt | while read metagenome
do
  if [ ! -e sketch_files/$metagenome.msh ]; then
    cat $read_storage/$metagenome*.fastq.gz > temp_MG_files/$metagenome.fastq.gz
    mash sketch -S 50 \
                -r \
                -m 2 \
                -k 21 \
                -s 100000 \
                -o sketch_files/$metagenome \
                temp_MG_files/$metagenome.fastq.gz
    rm -f temp_MG_files/$metagenome.fastq.gz
  else
    echo "Already sketched" $metagenome
  fi
done




# Make sure the key linking the metagenomes to the
# assemblies is up-to-date:
# ~/BLiMMP/metadata/assembly_key.tsv

# Also upload a list of the assemblies:
# ~/BLiMMP/metadata/assembly_list.txt
