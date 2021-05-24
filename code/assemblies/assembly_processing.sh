#!/bin/sh

#########################
# code/assemblies/assembly_processing.sh
# Benjamin D. Peterson
#########################


##################################################
##################################################
# Upload list files
##################################################
##################################################
mkdir ~/BLiMMP/metadata
# Upload a list of the metagenomes:
# ~/BLiMMP/metadata/metagenome_list.txt

# Make sure the key linking the metagenomes to the
# assemblies is up-to-date:
# ~/BLiMMP/metadata/assembly_key.tsv

# Also upload a list of the assemblies:
# ~/BLiMMP/metadata/assembly_list.txt





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
cd ~/BLiMMP/dataRaw/metagenomes

#source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
#conda activate lftp
#lftp -c 'set ssl:verify-certificate no set ftp:ssl-protect-data true set ftp:ssl-force true; open -u n210421_150PE_NVS1B_S1_Peterson,Shiezeicie7Eab5 -e "mirror -c; quit" ftp://gslanalyzer.qb3.berkeley.edu:990'





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

# Downlaode fastp
#mkdir ~/BLiMMP/code
#cd ~/BLiMMP/code
#wget http://opengene.org/fastp/fastp
#chmod a+x ./fastp
fastp=~/BLiMMP/code/fastp

cd ~/BLiMMP/dataRaw/metagenomes
read_storage=~/BLiMMP/dataEdited/metagenomes
ancillary_info=~/BLiMMP/dataEdited/metagenomes/reports
rawReads=~/BLiMMP/dataRaw/metagenomes

grep 'BLI20_MG_00' $ancillary_info/naming_key.tsv | while read line
do
  sequencingID=`echo $line | awk '{ print $1 }'`
  metagenomeID=`echo $line | awk '{ print $2 }'`
  echo "Cleaning up:"
  ls $rawReads/$sequencingID\_*R1*fastq.gz
  ls $rawReads/$sequencingID\_*R2*fastq.gz
  echo "into" $metagenomeID
  if [ ! -e $read_storage/$metagenomeID\_R1.fastq.gz ]; then
    echo "Processing" $metagenomeID
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
    echo "Already processed" $metagenome
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
grep 'BLI20_MG_00' $ancillary_info/naming_key.tsv | while read line
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
awk -F '\t' '{ print $1 }' ~/BLiMMP/metadata/metagenome_list.txt | while read metagenomeID
do
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
grep 'BLI20_MG_00' $ancillary_info/naming_key.tsv | while read line
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
IFS=$'\n'
echo -e "metagenomeID\tR1\tR2\tsingle\tmerged" > metagenome_coverage.tsv
cat ~/BLiMMP/metadata/metagenome_list.txt) | while read metagenome
do
  echo "Counting coverage in" $metagenome
  R1_count=$($code/kseq_fastq_base $read_storage/$metagenome\_R1.fastq.gz | \
                awk -F " " '{ print $5 }')
  R2_count=$($code/kseq_fastq_base $read_storage/$metagenome\_R2.fastq.gz | \
                awk -F " " '{ print $5 }')
  single_count=$($code/kseq_fastq_base $read_storage/$metagenome\_single.fastq.gz | \
                awk -F " " '{ print $5 }')
  merged_count=$($code/kseq_fastq_base $read_storage/$metagenome\_merged.fastq.gz | \
                awk -F " " '{ print $5 }')
  echo -e $metagenome"\t"$R1_count"\t"$R2_count"\t"$single_count"\t"$merged_count >> metagenome_coverage.tsv
done





############################################
############################################
# Metagenome assembly
############################################
############################################

######################
# Assemblies by metaSPAdes
######################
mkdir ~/BLiMMP/dataEdited/assemblies
mkdir ~/BLiMMP/dataEdited/assemblies/assembly_files

screen -S BLiMMP_metagenome_coassembly
cd ~/BLiMMP/dataEdited/assemblies
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=""
PERL5LIB=""

code=~/BLiMMP/code
assembly_grouping=~/BLiMMP/metadata/assembly_key.csv
read_storage=~/BLiMMP/dataEdited/metagenomes
output=~/BLiMMP/dataEdited/assemblies/assembly_files

#chmod +x $code/assembly_by_group.py

awk -F '\t' '{ print $1 }' ~/BLiMMP/metadata/assembly_list.txt | while read assembly
do
  if [ ! -d $output/$assembly ]; then
    mkdir $output/$assembly
  fi
  if [ ! -e $output/$assembly/scaffolds.fasta ]; then
    echo "Assembling" $assembly
    python $code/assembly_by_group.py $assembly \
                                      $assembly_grouping \
                                      $read_storage \
                                      $output/$assembly
    # To continue a paused run:
    # assembly=XXXXXX
    # metaspades.py --continue -o assembly
  else
    echo $assembly "already assembled"
  fi
done


######################
# Clean up assemblies
######################
screen -S HCC_clean_metagenome_assembly
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate anvio6.2
PYTHONPATH=''
cd ~/BLiMMP/dataEdited/assemblies
mkdir scaffolds
mkdir scaffolds/renaming_reports

awk -F '\t' '{ print $1 }' ~/BLiMMP/metadata/assembly_list.txt | while read assembly
do
  if [ -e assembly_files/$assembly/scaffolds.fasta ]; then
    echo $assembly "is done assembling. Let's clean it."
    if [ -e scaffolds/$assembly\_assembly.fna ]; then
      echo "Dude, you already cleaned the assembly for" $assembly". Relax."
    else
      echo "Cleaning the assembly for" $assembly
      anvi-script-reformat-fasta assembly_files/$assembly/scaffolds.fasta \
                                  -o scaffolds/$assembly\_assembly.fna \
                                  -l 1000 \
                                  --simplify-names \
                                  --prefix $assembly \
                                  --report-file scaffolds/renaming_reports/$assembly\_report_file.txt
    fi
  else
    echo "Yo, you gotta go back and assemble" $assembly
  fi
done
conda deactivate


############################################
############################################
# Predict ORFs
############################################
############################################

# Get set up
#cd /home/GLBRCORG/bpeterson26/BLiMMP/
#cp /home/GLBRCORG/bpeterson26/HellsCanyon/code/generalUse/cleanFASTA.py code/cleanFASTA.py
#mkdir reports
#cd reports
#mkdir outs errs logs


mkdir /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/assemblies/ORFs
chmod +x /home/GLBRCORG/bpeterson26/BLiMMP/code/ORF_prediction_assemblies.sh
condor_submit /home/GLBRCORG/bpeterson26/BLiMMP/code/ORF_prediction_assemblies.sub





############################################
############################################
# Get assembly stats
############################################
############################################

screen -S BLI_assembly_stats
cd ~/BLiMMP/dataEdited/assemblies
mkdir reports
code=~/BLiMMP/code
cp ~/HellsCanyon/code/assembly_processing/abyss-fac.pl $code/

# Use abyss-fac to check stats
awk -F '\t' '{ print $1 }' ~/BLiMMP/metadata/assembly_list.txt | while read assembly
do
  if [ -e scaffolds/$assembly\_assembly.fna ]; then
    echo $assembly "has been cleaned. Let's get some stats on it."
    if [ ! -e reports/$assembly\_report.txt ]; then
      perl $code/abyss-fac.pl scaffolds/$assembly\_assembly.fna \
          > reports/$assembly\_report.txt
    else
      echo "Already checked" $assembly
    fi
  else
    echo $assembly "has not been cleaned. Do that first."
  fi
done

# Aggregate stats
cd ~/BLiMMP/dataEdited/assemblies/reports
echo -e "n\tn:200\tL50\tmin\tN80\tN50\tN20\tmax\tsum\tassemblyID" > all_assemblies_stats.txt
for file in *report.txt
do
  tail -n +2 $file >> all_assemblies_stats.txt
done
rm BLI*_report.txt

# Clean up the report file
awk -F '\t' '{ print $1 }' ~/BLiMMP/metadata/assembly_list.txt | while read assembly
do
  sed "s/scaffolds\/$assembly\_assembly.fna/$assembly/" all_assemblies_stats.txt \
    > all_assemblies_stats.txt_edited
  mv -f all_assemblies_stats.txt_edited all_assemblies_stats.txt
done


# Count ORFs
cd ~/BLiMMP/dataEdited/assemblies
echo -e 'assemblyID\tORF_count' > reports/ORF_counts.tsv
awk -F '\t' '{ print $1 }' ~/BLiMMP/metadata/assembly_list.txt | while read assembly
do
  if [ -e ORFs/$assembly.faa ]; then
    echo "Count ORFs in" $assembly
    ORF_count=$(grep '>' ORFs/$assembly.faa | wc -l)
    echo -e "$assembly\t$ORF_count" >> reports/ORF_counts.tsv
  fi
done
