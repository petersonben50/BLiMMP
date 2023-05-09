#!/bin/sh

#########################
# SCG_abundance_in_assemblies.sh
# Benjamin D. Peterson

# This is the executable script to
# run SiCoGeCo on a range 
#########################


# Will read in these variables from
# submission file:
# geneName


#########################
# Activate environment
#########################
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PATH="/home/GLBRCORG/bpeterson26/BLiMMP/code/HomeBio/bin:$PATH"
HomeBio=~/BLiMMP/code/HomeBio


#########################
# Run ABA_SiCoGeCo.py script to pull out rp16 gene depths
#########################
python3 $HomeBio/bin/ABA_SiCoGeCo.py --gene_name $geneName \
                                      --scg_hmms_location $HomeBio/reference_data/HMMs/hmm_folder \
                                      --scg_hmms_key $HomeBio/reference_data/HMMs/rp16_key.csv \
                                      --assembly_location /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/assemblies/ORFs \
                                      --metagenome_list /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/metagenomes/reports/metagenome_list.txt \
                                      --mapping_location /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/mapping \
                                      --output_directory /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/scg_coverage/processing_directory \
                                      --homebio_location /home/GLBRCORG/bpeterson26/BLiMMP/code/HomeBio \
                                      --number_threads 7 \
                                      --length_to_trim 150 \
                                      --use_na
