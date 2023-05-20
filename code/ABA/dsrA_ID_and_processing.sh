#!/bin/sh

######################
# code/ABA/dsrA_ID_and_processing.sh
# Benjamin D. Peterson

# This set of scripts will contain the bash
# code we need to pull out the hgcA sequences
# from our assembly and process them.
######################


############################################
############################################
# Identify dsrA sequences with GID
############################################
############################################

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PATH="/home/GLBRCORG/bpeterson26/BLiMMP/code/HomeBio/bin:$PATH"

#cd ~/BLiMMP/code
#git clone https://github.com/petersonben50/HomeBio
#cd ~/BLiMMP/code/HomeBio
#git pull
HomeBio=~/BLiMMP/code/HomeBio

# Set up dataset
# rm -fr ~/BLiMMP/dataEdited/ABA/dsrA
python3 $HomeBio/bin/ABA_GID.py --orf_folder ~/BLiMMP/dataEdited/assemblies/ORFs/ \
                  --hmm $HomeBio/reference_data/HMMs/hmm_folder/TIGR02064.HMM \
                  --output_location ~/BLiMMP/dataEdited/ABA/dsrA \
                  --output_prefix dsrA \
                  --metagenome_list /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/metagenomes/reports/metagenome_list.txt \
                  --metagenomes_location /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/mapping/ \
                  --reference_aa_dataset $HomeBio/reference_data/sequence_databases/dsrA/muller_DsrAB_dataset_final.faa \
                  --number_threads 30
              > ~/BLiMMP/dataEdited/ABA/GID_log_dsrA.txt

