screen -S dsrA
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate batch_HMMs
PYTHONPATH=''
PERL5LIB=''
cd ~/BLiMMP/dataEdited/bin_based_analyses/hgcA_bins/metabolism/batch_HMMs/
mkdir manual_confirmation
HomeBio=~/BLiMMP/code/HomeBio
$HomeBio/reference_data/sequence_databases/dsrA/muller_DsrA_dataset_final.faa
