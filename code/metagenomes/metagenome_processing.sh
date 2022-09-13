#!/bin/sh

#########################
# code/metagenomes/metagenome_processing.sh
# Benjamin D. Peterson
#########################

#########################
# Requirements:
# 1. Raw metagenomes. Can be downloaded from NCBI.



##################################################
##################################################
# Trim metagenomes
##################################################
##################################################
screen -S BLiMMP_metagenome_trimming
mkdir ~/BLiMMP/dataEdited
mkdir ~/BLiMMP/dataEdited/metagenomes
mkdir ~/BLiMMP/dataEdited/metagenomes/reports

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
  metagenomeID=`echo $line | awk '{ print $2 }'`
  if [ ! -e $read_storage/$metagenomeID\_R1.fastq.gz ]; then
    echo "Processing" >> $ancillary_info/metagenome_run_info.txt
    ls $rawReads/$metagenomeID\_*R1*fastq.gz >> $ancillary_info/metagenome_run_info.txt
    ls $rawReads/$metagenomeID\_*R2*fastq.gz >> $ancillary_info/metagenome_run_info.txt
    $fastp --in1 $rawReads/$metagenomeID\_*R1*fastq.gz \
            --in2 $rawReads/$metagenomeID\_*R2*fastq.gz \
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
  metagenomeID=`echo $line | awk '{ print $2 }'`
  echo "Counting pre-trimming reads in" $metagenomeID >> $ancillary_info/metagenome_run_info.txt
  forwardCount=$(zgrep -c "^@" $rawReads/$metagenomeID\_*R1*fastq.gz)
  reverseCount=$(zgrep -c "^@" $rawReads/$metagenomeID\_*R1*fastq.gz)
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
  echo "Counting post-trimming reads in" $metagenomeID >> $ancillary_info/metagenome_run_info.txt
  forwardCount=$(zgrep -c "^@" $read_storage/$metagenomeID\_R1.fastq.gz)
  reverseCount=$(zgrep -c "^@" $read_storage/$metagenomeID\_R2.fastq.gz)
  singleCount=$(zgrep -c "^@" $read_storage/$metagenomeID\_single.fastq.gz)
  mergeCount=$(zgrep -c "^@" $read_storage/$metagenomeID\_merged.fastq.gz)
  echo -e $metagenomeID"\t"$forwardCount"\t"$reverseCount"\t"$singleCount"\t"$mergeCount >> $ancillary_info/metagenome_read_count.tsv
done



############################################
############################################
# Run Mash on metagenomes
############################################
############################################

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics



mkdir ~/BLiMMP/dataEdited/mash_data
mkdir ~/BLiMMP/dataEdited/mash_data/sketch_files
cd ~/BLiMMP/dataEdited/mash_data
read_storage=~/BLiMMP/dataEdited/metagenomes
ancillary_info=~/BLiMMP/dataEdited/metagenomes/reports

# Generate sketches for each metagenome
grep 'BLI2' $ancillary_info/naming_key.tsv | while read line
do
  metagenomeID=`echo $line | awk '{ print $2 }'`
  if [ ! -e sketch_files/$metagenomeID.msh ]; then
    echo "Generating mash sketch for" $metagenomeID >> $ancillary_info/metagenome_run_info.txt
    mash sketch -p 12 \
                -o sketch_files/$metagenomeID \
                -k 21 \
                -s 1000000 \
                -S 50 \
                -r \
                -m 2 \
                $read_storage/$metagenomeID\_R1.fastq.gz
  else
    echo "Already sketched" $metagenomeID >> $ancillary_info/metagenome_run_info.txt
  fi
done

# Combine sketches
mash paste BLI_MG_sketches sketch_files/BLI*.msh

# Calculate distances between metagenomes
mash dist -S 50 \
          BLI_MG_sketches.msh \
          BLI_MG_sketches.msh \
          > BLI_MG_sketches.dist



##########################################################
##########################################################
# Assemble needed sequences using EMIRGE
##########################################################
##########################################################
mkdir ~/BLiMMP/dataEdited/16S_from_MG
mkdir ~/BLiMMP/dataEdited/16S_from_MG/emirge_output
mkdir ~/BLiMMP/reports/emirge

######################
# Run EMIRGE on all samples
######################
cd ~/BLiMMP/dataEdited/16S_from_MG
emirge_output=~/BLiMMP/dataEdited/16S_from_MG/emirge_output
rm -f MGs_needing_emirge.txt
grep 'BLI2' $ancillary_info/naming_key.tsv | while read line
do
  metagenomeID=`echo $line | awk '{ print $2 }'`
  if [ ! -d $emirge_output/$metagenomeID ]; then
    echo "EMIRGE has not been run on" $metagenomeID
    echo $metagenomeID >> MGs_needing_emirge.txt
  else
    echo "EMIRGE was already run on" $metagenomeID
  fi
done

chmod +x /home/GLBRCORG/bpeterson26/HellsCanyon/code/executables/emirge_run.sh
condor_submit /home/GLBRCORG/bpeterson26/BLiMMP/code/emirge_run.sub
