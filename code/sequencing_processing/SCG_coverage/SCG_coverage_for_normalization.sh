#!/bin/sh

######################
# code/SCG_coverage_for_normalization.sh
# Benjamin D. Peterson
######################

######################
# Get set up
######################
#cd BLiMMP/code/HomeBio/
#git pull
#rm -rf /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/scg_coverage
#rm -rf /home/GLBRCORG/bpeterson26/BLiMMP/reports/scg_coverage
mkdir /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/scg_coverage
mkdir /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/scg_coverage/processing_directory
mkdir /home/GLBRCORG/bpeterson26/BLiMMP/reports/scg_coverage
cd /home/GLBRCORG/bpeterson26/BLiMMP/reports/scg_coverage
mkdir outs errs logs

######################
# Gather coverage information
######################
chmod +x /home/GLBRCORG/bpeterson26/BLiMMP/code/execute_SCG_abundance_in_assemblies.sh
condor_submit /home/GLBRCORG/bpeterson26/BLiMMP/code/submit_SCG_abundance_in_assemblies.sub
cd ~/BLiMMP/dataEdited/scg_coverage
cp */*/*_coverage.tsv  */*/*_G2A.tsv */*/*_cluster_data.tsv ./


######################
# Cluster sequences
######################
cd ~/BLiMMP/dataEdited/scg_coverage/processing_directory
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
HomeBio=/home/GLBRCORG/bpeterson26/BLiMMP/code/HomeBio

cat /home/GLBRCORG/bpeterson26/BLiMMP/references/scg_list.txt | while read gene_name
do
    cd-hit -g 1 \
            -i $gene_name/$gene_name.fna \
            -o $gene_name/$gene_name\_drep_take2.fna \
            -c 0.995 \
            -n 5 \
            -d 0
    python $HomeBio/bin/FM_CDHIT_parsing.py --clstr_in $gene_name/$gene_name\_drep_take2.fna.clstr \
                                            --clstr_out /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/scg_coverage/$gene_name\_finalCluster.tsv
done