#!/bin/sh

######################
# code/assemblies/mapping_assemblies.sh
# Benjamin D. Peterson
######################

############################################
############################################
# Calculate total coverage of each metagenome
############################################
############################################

screen -S BLI_MG_coverage_counting
mkdir ~/BLiMMP/dataEdited/mapping/
cd ~/BLiMMP/dataEdited/mapping/
code=~/BLiMMP/code/readfq-master
read_storage=~/BLiMMP/dataEdited/metagenomes

echo -e "metagenomeID\tR1\tR2\tsingle\tmerged" > metagenome_coverage.tsv
awk -F '\t' '{ print $1 }' ~/BLiMMP/metadata/metagenome_list.txt | while read metagenome
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
# Mapping reads to scaffolds
############################################
############################################

mkdir ~/BLiMMP/dataEdited/mapping

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

awk -F '\t' '{ print $1 }' ~/BLiMMP/metadata/assembly_list.txt | while read assembly
do
  if [ -e $scaffolds/$assembly\_assembly.fna ]; then
    echo $assembly "has been cleaned, let's map to it"
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

# Identify mapping pairs that have not been completed
cd ~/BLiMMP/dataEdited/mapping/
rm -f ~/BLiMMP/metadata/mapping_key_pairsRemaining.tsv
touch ~/BLiMMP/metadata/mapping_key_pairsRemaining.tsv
cat ~/BLiMMP/metadata/mapping_key.tsv | while read -r line
do
  metagenome=`echo $line | awk -F '\t' '{ print $1 }'`
  assembly=`echo $line | awk -F '\t' '{ print $2 }'`
  if [ -e $metagenome\_to_$assembly.bam ]; then
    echo $metagenome "mapping to" $assembly "already done"
  else
    echo $metagenome "mapping to" $assembly "still needs to be done"
    echo $metagenome", "$assembly >> ~/BLiMMP/metadata/mapping_key_pairsRemaining.tsv
  fi
done


# Set up for mapping
mkdir /home/GLBRCORG/bpeterson26/BLiMMP/reports/readMappingAssemblies/
cd /home/GLBRCORG/bpeterson26/BLiMMP/reports/readMappingAssemblies/
mkdir errs outs logs
chmod +x /home/GLBRCORG/bpeterson26/BLiMMP/code/read_mapping_assemblies.sh

# Run mapping
condor_submit /home/GLBRCORG/bpeterson26/BLiMMP/code/read_mapping_assemblies.sub
