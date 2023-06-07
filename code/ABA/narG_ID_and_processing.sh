#!/bin/sh

######################
# code/ABA/narG_ID_and_processing.sh
# Benjamin D. Peterson
######################


############################################
############################################
# Identify narG sequences with GID
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
# rm -fr ~/BLiMMP/dataEdited/ABA/narG
python3 $HomeBio/bin/ABA_GID.py --orf_folder ~/BLiMMP/dataEdited/assemblies/ORFs/ \
                  --hmm $HomeBio/reference_data/HMMs/hmm_folder/TIGR01580.HMM \
                  --output_location ~/BLiMMP/dataEdited/ABA/narG \
                  --output_prefix narG \
                  --metagenome_list /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/metagenomes/reports/metagenome_list.txt \
                  --metagenomes_location /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/mapping/ \
                  --metatranscriptome_location /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/metatranscriptomes/alignment \
                  --reference_aa_dataset $HomeBio/reference_data/sequence_databases/narG/NarG_luke_database.faa \
                  --number_threads 30 \
                  --cluster_cutoff 0.8 \
                  > ~/BLiMMP/dataEdited/ABA/GID_log_narG.txt

############################################
############################################
# Cluster seqs at defined cut-off
############################################
############################################
cd-hit -g 1 \
        -i ~/BLiMMP/dataEdited/ABA/narG/narG.faa \
        -o ~/BLiMMP/dataEdited/ABA/narG/working_directory/narG_derep_99.faa \
        -c 0.99 \
        -n 5 \
        -d 0
python $HomeBio/bin/FM_CDHIT_parsing.py --clstr_in  ~/BLiMMP/dataEdited/ABA/narG/working_directory/narG_derep_99.faa.clstr \
                                        --clstr_out ~/BLiMMP/dataEdited/ABA/narG/narG_derep_99.tsv