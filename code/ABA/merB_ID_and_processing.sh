#!/bin/sh

######################
# code/ABA/merB_ID_and_processing.sh
# Benjamin D. Peterson

# This set of scripts will contain the bash
# code we need to pull out the merB sequences
# from our assembly and process them.
######################


############################################
############################################
# Identify merB sequences with GID
############################################
############################################

screen -S merB_GID
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PATH="/home/GLBRCORG/bpeterson26/BLiMMP/code/HomeBio/bin:$PATH"

#cd ~/BLiMMP/code
#git clone https://github.com/petersonben50/HomeBio
#cd ~/BLiMMP/code/HomeBio
#git pull

# Prep references
python ~/BLiMMP/code/cleanFASTA.py ~/BLiMMP/references/merB/Supplementary_Dataset_4_MerB_homologs.faa
sed 's/-//g' /home/GLBRCORG/bpeterson26/BLiMMP/references/merB/Supplementary_Dataset_4_MerB_homologs.faa_temp.fasta \
    > /home/GLBRCORG/bpeterson26/BLiMMP/references/merB/merB_refs.faa
# Cluster at 70% identity
cd-hit -g 1 \
        -i /home/GLBRCORG/bpeterson26/BLiMMP/references/merB/merB_refs.faa \
        -o /home/GLBRCORG/bpeterson26/BLiMMP/references/merB/merB_refs_drep.faa \
        -c 0.70 \
        -n 5 \
        -d 0

# Set up dataset
# rm -fr ~/BLiMMP/dataEdited/ABA/merB
python3 ~/BLiMMP/code/ABA_GID.py --orf_folder ~/BLiMMP/dataEdited/assemblies/ORFs/ \
                  --hmm ~/references/merB/merB_christakis.hmm \
                  --output_location ~/BLiMMP/dataEdited/ABA/merB \
                  --output_prefix merB \
                  --metagenome_list /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/metagenomes/reports/metagenome_list.txt \
                  --metagenomes_location /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/mapping/ \
                  --metatranscriptome_location /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/metatranscriptomes/alignment \
                  --reference_aa_dataset /home/GLBRCORG/bpeterson26/BLiMMP/references/merB/merB_refs_drep.faa \
                  --number_threads 30 \
                  > ~/BLiMMP/dataEdited/ABA/GID_log_merB.txt

rm -rf ~/BLiMMP/dataEdited/ABA/merB/working_directory/


######################
# Concatenate merB faa files with references
######################
cat /home/GLBRCORG/bpeterson26/BLiMMP/references/merB/merB_refs_drep.faa \
    /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/ABA/merB/merB.faa \
    > /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/ABA/merB/merB_with_refs.faa
sed 's/assembly/a/' /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/ABA/merB/merB_with_refs.faa > /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/ABA/merB/merB_with_refs_short_names.faa


# Download the merB_with_refs.faa file to your local machine
# Upload to http://prodata.swmed.edu/promals3d/promals3d.php
# Enter 3F0P for the PDB ID.
# Run with default settings.
# Download the alignment and save as merB_with_refs_aligned.faa
# Results: http://prodata.swmed.edu/promals3d/showResult.php?name=QUERY_Ryzrl7&email=&target_name=QUERY_Ryzrl7&jobname=merB_alignment_with_refs
# Save to local machine as merB_with_refs_short_names.afa


######################
# First pull out merB+ scaffolds
######################
screen -S BLI_merB_context
cd /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/ABA/merB
# rm -r scaffolds
mkdir working_directory/scaffolds
# Get list of gene names:
grep '>' merB.faa | sed 's/>//' | awk -F ' ' '{print $1}' > merB_genes.txt

# Pull out scaffold FNA seqences and GFF entries
rm -f working_directory/scaffolds/merB_scaffolds.fna working_directory/scaffolds/merB_scaffolds.gff
awk -F '_' '{ print $1"_"$2"_"$3 }' merB_genes.txt | while read scaffold
do
  assemblyName=$(echo $scaffold | awk -F '_' '{ print $1"_"$2 }')
  echo "Pulling out" $scaffold "from" $assemblyName
  grep -A 1 -m 1 $scaffold\$ ~/BLiMMP/dataEdited/assemblies/scaffolds/$assemblyName\_assembly.fna \
      >> working_directory/scaffolds/merB_scaffolds.fna
  echo "Pulling out" $scaffold "GFF entries"
  awk -v scaffold="$scaffold" '{ if ($1 == scaffold) print }' ~/BLiMMP/dataEdited/assemblies/ORFs/$assemblyName.gff \
      >> working_directory/scaffolds/merB_scaffolds.gff
done















######################
# Isolate gene neighborhoods
######################

screen -S BLI_merB_gene_neighborhood
cd ~/BLiMMP/dataEdited/ABA/merB
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate py_viz
PYTHONPATH=''
scripts=~/BLiMMP/code

awk -F '\t' '{ print $1 }' merB_G2A.tsv | while read merB_id
do
  scaffold_id=$(echo $merB_id | cut -d '_' -f 1-3)
  gene_id=$(echo $merB_id | \
              cut -d '_' -f 3-4 | \
              sed 's/^0*//g')
  if [ ! -f working_directory/scaffolds/temp_$gene_id\_neighborhood.gff ]
  then
    echo "Working on" $merB_id", on scaffold" $scaffold_id
    awk -F '\t' -v scaffold_id="$scaffold_id" '$1 == scaffold_id { print $0 }' working_directory/scaffolds/merB_scaffolds.gff > working_directory/scaffolds/temp_scaffolds.gff
    echo "Searching for" $gene_id
    python $scripts/gene_neighborhood_extraction.py working_directory/scaffolds/temp_scaffolds.gff \
                                                    working_directory/scaffolds/merB_scaffolds.fna \
                                                    $gene_id \
                                                    5000 \
                                                    working_directory/scaffolds/temp_$gene_id
    rm -f working_directory/scaffolds/temp_scaffolds.gff
  else
    echo $merB_id "already processed"
  fi
done
cat working_directory/scaffolds/temp_*.gff > merB_geneNeighborhood_raw.gff
cat working_directory/scaffolds/temp_*.fna > merB_geneNeighborhood_raw.fna
rm -f working_directory/scaffolds/temp_*.gff working_directory/scaffolds/temp_*.fna

# Pull out ORF files
sed 's/\tProdigal.*ID=[0-9]*_/_/' merB_geneNeighborhood_raw.gff | awk -F ';' '{ print $1 }' | while read orf_id
do
  assembly_id=$(echo $orf_id | cut -d '_' -f 1-2)
  echo "Pulling out" $orf_id "from" $assembly_id
  grep -A 1 -m 1 $orf_id$ ~/BLiMMP/dataEdited/assemblies/ORFs/$assembly_id.faa \
      >> merB_geneNeighborhood_ORFs.faa
done





# Run KOFAMscan on neighborhood ORFs
# See instructions for ~/BLiMMP/code/kofam_scan-1.3.0/config.yml file setup in merB_bin_analysis.sh
screen -S kofamscan
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate kofamscan
PYTHONPATH=""
PERL5LIB=""

merB_folder=/home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/ABA/merB
cd /home/GLBRCORG/bpeterson26/BLiMMP/code/kofam_scan-1.3.0
./exec_annotation -f detail-tsv \
                  -o $merB_folder/working_directory/merB_GN_kofamscan.tsv \
                  $merB_folder/merB_geneNeighborhood_ORFs.faa
head -n 1 $merB_folder/working_directory/merB_GN_kofamscan.tsv > $merB_folder/working_directory/merB_GN_kofamscan_qualityHits.tsv
grep '*' $merB_folder/working_directory/merB_GN_kofamscan.tsv >> $merB_folder/working_directory/merB_GN_kofamscan_qualityHits.tsv
# Remove the first column of the file using awk, keeping the columns separated by tabs
awk -F '\t' '{ print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7 }' $merB_folder/working_directory/merB_GN_kofamscan_qualityHits.tsv > $merB_folder/merB_GN_kofamscan_qualityHits_clean.tsv
exit