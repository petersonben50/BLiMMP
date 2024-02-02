#!/bin/sh

######################
# code/hgcA_analysis/hgcA_analysis.sh
# Benjamin D. Peterson

# This set of scripts will contain the bash
# code we need to pull out the hgcA sequences
# from our assembly and process them.
######################


############################################
############################################
# Identify hgcA sequences with GID
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
HomeBio=~/BLiMMP/code/HomeBio/bin

# Set up dataset
# rm -fr ~/BLiMMP/dataEdited/ABA/hgcA
mkdir ~/BLiMMP/dataEdited/ABA
python3 $HomeBio/ABA_GID.py --orf_folder ~/BLiMMP/dataEdited/assemblies/ORFs/ \
                  --hmm ~/references/hgcA/hgcA.hmm \
                  --output_location ~/BLiMMP/dataEdited/ABA/hgcA \
                  --output_prefix hgcA \
                  --metagenome_list /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/metagenomes/reports/metagenome_list.txt \
                  --metagenomes_location /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/mapping/ \
                  --reference_aa_dataset ~/references/13105370/Hg-MATE-Db.v1.01142021_ISOCELMAG_HgcA.fas \
              > ~/BLiMMP/dataEdited/ABA/GID_log_hgcA.txt




############################################
############################################
# Generate tree for classification
############################################
############################################
cd ~/BLiMMP/dataEdited/ABA/hgcA
sed -e '/^[A-Z]/s/-//g' ~/BLiMMP/references/hgcA/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.refpkg/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full_align-bmge.fasta > working_directory/HgMATE_rough_phy_db_all.faa
sed -i -e '/^-/s/-//g' working_directory/HgMATE_rough_phy_db_all.faa
cd-hit -g 1 \
        -i working_directory/HgMATE_rough_phy_db_all.faa \
        -o working_directory/HgMATE_rough_phy_db.faa \
        -c 0.70 \
        -n 5 \
        -d 0
cat hgcA_for_phylogeny.faa \
    working_directory/HgMATE_rough_phy_db.faa \
    > working_directory/hgcA_for_phylogeny_with_refs.faa
muscle -align working_directory/hgcA_for_phylogeny_with_refs.faa \
        -output working_directory/hgcA_for_phylogeny_with_refs.afa
trimal -in working_directory/hgcA_for_phylogeny_with_refs.afa \
        -out working_directory/hgcA_for_phylogeny_with_refs_cleaned.afa \
        -gt 0.5
FastTree working_directory/hgcA_for_phylogeny_with_refs_cleaned.afa \
    > hgcA_phylogeny.tree
# Download this to my local computer.
conda deactivate


############################################
############################################
# Classify hgcA seqs with pplacer workflow
############################################
############################################
conda activate hgcA_classifier
PYTHONPATH=""
PERL5LIB=""
python $HomeBio/HG_hgcA_AutoClass.py --fasta_file /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/ABA/hgcA/hgcA.faa \
                                    --ref_package /home/GLBRCORG/bpeterson26/BLiMMP/references/hgcA/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.refpkg \
                                    --output_name hgcA_for_autoClass \
                                    --output_location /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/ABA/hgcA


############################################
############################################
# Genomic context for hgcA
############################################
############################################

######################
# First pull out hgcA+ scaffolds
######################
screen -S BLI_hgcA_context
cd /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/ABA/hgcA
# rm -r scaffolds
mkdir working_directory/scaffolds
# Get list of gene names:
grep '>' hgcA.faa | sed 's/>//' | awk -F ' ' '{print $1}' > hgcA_genes.txt

# Pull out scaffold FNA seqences and GFF entries
rm -f working_directory/scaffolds/hgcA_scaffolds.fna working_directory/scaffolds/hgcA_scaffolds.gff
awk -F '_' '{ print $1"_"$2"_"$3 }' hgcA_genes.txt | while read scaffold
do
  assemblyName=$(echo $scaffold | awk -F '_' '{ print $1"_"$2 }')
  echo "Pulling out" $scaffold "from" $assemblyName
  grep -A 1 -m 1 $scaffold\$ ~/BLiMMP/dataEdited/assemblies/scaffolds/$assemblyName\_assembly.fna \
      >> working_directory/scaffolds/hgcA_scaffolds.fna
  echo "Pulling out" $scaffold "GFF entries"
  awk -v scaffold="$scaffold" '{ if ($1 == scaffold) print }' ~/BLiMMP/dataEdited/assemblies/ORFs/$assemblyName.gff \
      >> working_directory/scaffolds/hgcA_scaffolds.gff
done


######################
# Search downstream genes for hgcB
######################
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate py_viz
PYTHONPATH=''
#IFS=$'\n'
scripts=~/BLiMMP/code
ORFs=~/BLiMMP/dataEdited/assemblies/ORFs
cd ~/BLiMMP/dataEdited/ABA/hgcA
mkdir working_directory/hgcB
rm -f working_directory/hgcB/downstream_gene_list.txt
python $scripts/retrieve_downstream_gene_name.py \
          hgcA_genes.txt \
          working_directory/scaffolds/hgcA_scaffolds.gff \
          working_directory/hgcB/downstream_gene_list.txt
conda deactivate

# Extract downstream amino acid sequences
cd working_directory/hgcB
rm -f downstream_genes.faa
cat downstream_gene_list.txt | while read gene
do
  assemblyName=$(echo $gene | cut -d "_" -f1,2)
  echo "Pulling out" $gene "faa entry from" $assemblyName
  grep -A 1 -m 1 \
        $gene$ \
        $ORFs/$assemblyName.faa \
      >> downstream_genes.faa
done
grep '>' downstream_genes.faa | \
  sed 's/>//' > downstream_genes_present.txt

# Search adjacent genes with hgcB HMM
cd ~/BLiMMP/dataEdited/ABA/hgcA/working_directory/hgcB
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''
scripts=~/BLiMMP/code/HomeBio/fasta_manipulation
hmmsearch --tblout hgcB.out \
          --cpu 4 \
          -T 30 \
          ~/references/hgcA/hgcB_5M.HMM \
          downstream_genes.faa \
          > hgcB_report.txt

# Check the hgcB hits
rm -f hgcB.faa
grep -v "#" hgcB.out | awk '{ print $1 }' | while read geneID
do
  grep -A 1 $geneID$ downstream_genes.faa >> hgcB.faa
done
sed -i 's/\*//' hgcB.faa
grep '>' hgcB.faa | wc -l
# Align sequences to HMM
hmmalign --amino \
         -o hgcB.sto \
         ~/references/hgcA/hgcB_5M.HMM \
         hgcB.faa
# Convert alignment to fasta format
python ~/BLiMMP/code/convert_stockhold_to_fasta.py hgcB.sto
python ~/BLiMMP/code/cleanFASTA.py hgcB.afa
mv -f hgcB.afa_temp.fasta hgcB.afa
cat hgcB.afa | while read line; do echo ${line^^} >> hgcB_clean.afa; done
sed -i 's/_ASSEMBLY/_assembly/' hgcB_clean.afa
sed -i 's/_COASSEMBLY/_coassembly/' hgcB_clean.afa
cp hgcB_clean.afa ../../
# Check sequences in hgcA_data_processing.Rmd, on my local computer.

# Check GFF entries for cut hgcB seqs
cd ~/BLiMMP/dataEdited/ABA/hgcA/working_directory/scaffolds
declare -a hgcB_trunc=("BLI20_assembly005_000000188864_1" "BLI21_coassembly_000000421460_3" "BLI21_coassembly_000000483203_1" "BLI20_assembly001_000000070404_2" "BLI21_coassembly_000001026396_1" "BLI20_coassembly_000000275727_4")
for gene_id in "${hgcB_trunc[@]}"
do
  scaffold_id=$(echo $gene_id | cut -d '_' -f 1-3)
  grep $scaffold_id hgcA_scaffolds.gff
done
# They all check out, so I'll just use the hgcB.faa
# file as our final set of hgcB seqs.


######################
# Isolate gene neighborhoods
######################

screen -S BLI_hgcA_gene_neighborhood
cd ~/BLiMMP/dataEdited/ABA/hgcA
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate py_viz
PYTHONPATH=''
scripts=~/BLiMMP/code

cat ~/BLiMMP/dataEdited/ABA/hgcA/hgcA_genes.txt | while read hgcA_id
do
  scaffold_id=$(echo $hgcA_id | cut -d '_' -f 1-3)
  gene_id=$(echo $hgcA_id | \
              cut -d '_' -f 3-4 | \
              sed 's/^0*//g')
  if [ ! -f working_directory/scaffolds/temp_$gene_id\_neighborhood.gff ]
  then
    echo "Working on" $hgcA_id", on scaffold" $scaffold_id
    awk -F '\t' -v scaffold_id="$scaffold_id" '$1 == scaffold_id { print $0 }' working_directory/scaffolds/hgcA_scaffolds.gff > working_directory/scaffolds/temp_scaffolds.gff
    echo "Searching for" $gene_id
    python $scripts/gene_neighborhood_extraction.py working_directory/scaffolds/temp_scaffolds.gff \
                                                    working_directory/scaffolds/hgcA_scaffolds.fna \
                                                    $gene_id \
                                                    5000 \
                                                    working_directory/scaffolds/temp_$gene_id
    rm -f working_directory/scaffolds/temp_scaffolds.gff
  else
    echo $hgcA_id "already processed"
  fi
done
cat working_directory/scaffolds/temp_*.gff > hgcA_geneNeighborhood_raw.gff
cat working_directory/scaffolds/temp_*.fna > hgcA_geneNeighborhood_raw.fna
rm -f working_directory/scaffolds/temp_*.gff working_directory/scaffolds/temp_*.fna

# Pull out ORF files
sed 's/\tProdigal.*ID=[0-9]*_/_/' hgcA_geneNeighborhood_raw.gff | awk -F ';' '{ print $1 }' | while read orf_id
do
  assembly_id=$(echo $orf_id | cut -d '_' -f 1-2)
  echo "Pulling out" $orf_id "from" $assembly_id
  grep -A 1 -m 1 $orf_id$ ~/BLiMMP/dataEdited/assemblies/ORFs/$assembly_id.faa \
      >> hgcA_geneNeighborhood_ORFs.faa
done

# Run KOFAMscan on neighborhood ORFs
# See instructions for ~/BLiMMP/code/kofam_scan-1.3.0/config.yml file setup in hgcA_bin_analysis.sh
screen -S kofamscan
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate kofamscan
PYTHONPATH=""
PERL5LIB=""

hgcA_folder=/home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/ABA/hgcA
cd /home/GLBRCORG/bpeterson26/BLiMMP/code/kofam_scan-1.3.0
./exec_annotation -f detail-tsv \
                  -o $hgcA_folder/working_directory/hgcA_GN_kofamscan.tsv \
                  $hgcA_folder/hgcA_geneNeighborhood_ORFs.faa
head -n 1 $hgcA_folder/working_directory/hgcA_GN_kofamscan.tsv > $hgcA_folder/working_directory/hgcA_GN_kofamscan_qualityHits.tsv
grep '*' $hgcA_folder/working_directory/hgcA_GN_kofamscan.tsv >> $hgcA_folder/working_directory/hgcA_GN_kofamscan_qualityHits.tsv
# Remove the first column of the file using awk, keeping the columns separated by tabs
awk -F '\t' '{ print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7 }' $hgcA_folder/working_directory/hgcA_GN_kofamscan_qualityHits.tsv > $hgcA_folder/hgcA_GN_kofamscan_qualityHits_clean.tsv

conda deactivate

# Find references from RefSeq
screen -S blast
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=""
PERL5LIB=""
hgcA_folder=/home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/ABA/hgcA
head -n 2 $hgcA_folder/hgcA_geneNeighborhood_ORFs.faa > $hgcA_folder/hgcA_geneNeighborhood_ORFs_TESTONLY.faa
/opt/bifxapps/ncbi-blast-2.6.0+/bin/blastp -query $hgcA_folder/hgcA_geneNeighborhood_ORFs.faa \
                                            -db /opt/bifxapps/ncbi-blastdb/refseq_protein \
                                            -evalue 0.0001 \
                                            -outfmt '6 qseqid evalue sseqid sseq staxids' \
                                            -max_target_seqs 5 \
                                            -num_threads 50 \
                                            -out $hgcA_folder/working_directory/hgcA_GN_blast.tsv
cd $hgcA_folder
#awk -F '\t' '{ print ">"$1"\n"$2 }' working_directory/hgcA_GN_blast.tsv > hgcA_GN.faa
awk -F '\t' '{ print $3 }' working_directory/hgcA_GN_blast.tsv | \
  awk -F '|ref|' '{ print $2 }' | \
  sed 's/|//g' | \
  sort | uniq > working_directory/hgcA_GN_refseq_list.txt
awk -F '\t' '{ print $1"\t"$2"\t"$3 }' working_directory/hgcA_GN_blast.tsv > hgcA_GN_blast_data.tsv

epost -db protein -input working_directory/hgcA_GN_refseq_list.txt | \
    esummary | \
    xtract -pattern DocumentSummary -element AccessionVersion,Title,TaxId > hgcA_GN_blast_refseq_data.tsv

# Files to bring to computer:
#hgcA_geneNeighborhood_raw.gff
#hgcA_GN_kofamscan_qualityHits_clean.tsv
#hgcA_GN_blast_data.tsv
#hgcA_GN_blast_refseq_data.tsv
#hgcB_clean.faa
# Download to dataEdited/ABA/hgcA/GN



############################################
############################################
# Phylogenetic analysis of hgcA
############################################
############################################

#########################
# Generate good tree with Hg-MATE seqs using RAxML
#########################

# Done locally
BLiMMP
cp /Users/benjaminpeterson/Documents/research/Hg_MATE/versions/v1.01142021/Hg-MATE-Db.v1.01142021_ISOCELMAG_Hgc.fas \
    references
rm -f dataEdited/hgcA_analysis/phylogeny/raxml/HgMate_reference_seqs_to_use.faa
cut -d"_" -f1-3 dataEdited/hgcA_analysis/phylogeny/raxml/reference_names_to_use.txt | while read reference_name
do
  grep -A 1 \
    $reference_name \
    references/Hg-MATE-Db.v1.01142021_ISOCELMAG_Hgc.fas \
    >> dataEdited/hgcA_analysis/phylogeny/raxml/HgMate_reference_seqs_to_use.faa
done


# Done on GLBRC
screen -S BLI_hgcA_tree
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=""
PERL5LIB=''
cd ~/BLiMMP/dataEdited/hgcA_analysis/phylogeny/raxml
# Upload needed references
cat 5M_bin_seqs.faa \
    HgMate_reference_seqs_to_use.faa \
    jones_hgcA_seqs.faa \
    ../hgcA_for_phylogeny.faa \
    hgcA_paralogs_for_rooting.faa \
    > hgcA_for_tree.faa

# Generate alignment
muscle -in hgcA_for_tree.faa \
        -out hgcA_for_tree.afa

# Upload masked alignment
# Then run RAxML to generate tree
raxml=/opt/bifxapps/raxml-8.2.11/raxmlHPC-PTHREADS
$raxml -f a \
        -p 283976 \
        -m PROTGAMMAAUTO \
        -N autoMRE \
        -x 2381 \
        -T 20 \
        -s hgcA_for_tree_masked.afa \
        -n hgcA


# Check out sequence similarity of study seqs to references
cd ~/BLiMMP/dataEdited/hgcA_analysis
mkdir derep_references
cdhit=~/programs/cdhit-master
$cdhit/cd-hit -g 1 \
              -i phylogeny/raxml/hgcA_for_tree.faa \
              -o derep_references/hgcA_ref.faa \
              -c 0.97 \
              -n 5 \
              -d 0
$cdhit/cd-hit -g 1 \
              -i phylogeny/raxml/hgcA_for_tree.faa \
              -o derep_references/hgcA_ref_80.faa \
              -c 0.80 \
              -n 5 \
              -d 0
