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


############################################
############################################
# Calculate total coverage of each metagenome
############################################
############################################

screen -S BLI_MG_coverage_counting
cd ~/BLiMMP/dataEdited/mapping/
code=~/BLiMMP/code/readfq-master
read_storage=~/BLiMMP/dataEdited/metagenomes

echo -e "metagenomeID\tR1\tR2\tsingle" > metagenome_coverage.tsv
cat /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/metagenomes/reports/metagenome_list.txt | while read metagenome
do
  echo "Counting coverage in" $metagenome
  R1_count=$($code/kseq_fastq_base $read_storage/$metagenome\_R1.fastq.gz | \
                awk -F " " '{ print $5 }')
  R2_count=$($code/kseq_fastq_base $read_storage/$metagenome\_R2.fastq.gz | \
                awk -F " " '{ print $5 }')
  single_count=$($code/kseq_fastq_base $read_storage/$metagenome\_R2.fastq.gz | \
                    awk -F " " '{ print $5 }')
  echo -e $metagenome"\t"$R1_count"\t"$R2_count"\t"$single_count >> metagenome_coverage.tsv
done


############################################
############################################
# Mapping reads to scaffolds
############################################
############################################

######################
# Generate indices
######################
screen -S BLI_assembly_indexing
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PERL5LIB=""
mkdir ~/BLiMMP/dataEdited/mapping/indices
cd ~/BLiMMP/dataEdited/mapping/indices
scaffolds=/home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/assemblies/scaffolds

cat ~/BLiMMP/metadata/assembly_list.txt | while read assembly
do
  if [ -e $scaffolds/$assembly\_assembly.fna ]; then
    if [ ! -e bowtie_index_$assembly.1.bt2 ]; then
      echo "Building index for" $assembly
      bowtie2-build $scaffolds/$assembly\_assembly.fna \
                    bowtie_index_$assembly
    else
      echo "Already built index for" $assembly
    fi
  else
    echo $assembly "has not been cleaned yet"
  fi
done


######################
# Map reads and process output
######################
# Set up for mapping
mkdir /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/mapping/reports
# Upload the mapping_key.tsv file here
mkdir /home/GLBRCORG/bpeterson26/BLiMMP/reports/readMappingAssemblies/
cd /home/GLBRCORG/bpeterson26/BLiMMP/reports/readMappingAssemblies/
mkdir errs outs logs
chmod +x /home/GLBRCORG/bpeterson26/BLiMMP/code/read_mapping_assemblies.sh

# Run mapping
condor_submit /home/GLBRCORG/bpeterson26/BLiMMP/code/read_mapping_assemblies.sub
