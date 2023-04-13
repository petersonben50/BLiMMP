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
#cd HomeBio
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
FastTree working_directory/hgcA_for_phylogeny_with_refs.afa \
    > hgcA_phylogeny.tree
# Download this to my local computer.



############################################
############################################
# Classify hgcA seqs with pplacer workflow
############################################
############################################

mkdir ~/BLiMMP/references
mkdir ~/BLiMMP/references/hgcA

screen -S BLI_hgcA_pplacer
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=""
PERL5LIB=""
references=~/BLiMMP/references/hgcA
workingDirectory=~/BLiMMP/dataEdited/hgcA_analysis
scripts=~/BLiMMP/code/HomeBio/fasta_manipulation/
mkdir $workingDirectory/classification

# Generate alignment of sequences of interest
cd $references/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.refpkg
python $scripts/convert_stockhold_to_fasta.py Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.stockholm
mv Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.afackholm $workingDirectory/classification/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.afa
cd $workingDirectory
$scripts/cleanFASTA.py classification/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.afa
sed -e '/^[A-Z]/s/-//g' classification/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.afa_tempCleanedFile > classification/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.faa
sed -i -e '/^-/s/-//g' classification/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.faa
rm -f classification/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.afa_tempCleanedFile

cat identification/hgcA_good.faa \
    classification/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.faa \
    > classification/hgcA_for_classification.faa
cd classification/
sed -i 's/\*//' hgcA_for_classification.faa
muscle -super5 hgcA_for_classification.faa \
        -output hgcA_for_classification.afa
# Convert to stockholm format
python $scripts/convert_fasta_to_stockholm.py hgcA_for_classification.afa
conda deactivate

# Run pplacer
conda activate hgcA_classifier
pplacer -p \
        --keep-at-most 1 \
        --max-pend 1 \
        -c $references/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.refpkg/ \
        hgcA_for_classification.sto

# Make sqlite database of reference
rppr prep_db \
      --sqlite Hg_MATE_classify \
      -c $references/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.refpkg

# Generate taxonomic assignments using guppy
guppy classify -c $references/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.refpkg \
               --pp \
               --sqlite Hg_MATE_classify \
               hgcA_for_classification.jplace

# Save out this data to csv
guppy to_csv --point-mass \
              --pp \
              -o hgcA_classification.csv \
              hgcA_for_classification.jplace

# Visualize placements on FastTree
guppy tog --pp \
          -o hgcA_classification.nwk \
          hgcA_for_classification.jplace






############################################
############################################
# Genomic context for hgcA
############################################
############################################

######################
# First pull out hgcA+ scaffolds
######################
screen -S BLI_hgcA_context
cd ~/BLiMMP/dataEdited/hgcA_analysis
# rm -r scaffolds
mkdir scaffolds

# Pull out scaffold FNA seqences and GFF entries
rm -f scaffolds/hgcA_scaffolds.fna scaffolds/hgcA_scaffolds.gff
awk -F '_' '{ print $1"_"$2"_"$3 }' identification/hgcA_good.txt | while read scaffold
do
  assemblyName=$(echo $scaffold | awk -F '_' '{ print $1"_"$2 }')
  echo "Pulling out" $scaffold "from" $assemblyName
  grep -A 1 $scaffold\$ ~/BLiMMP/dataEdited/assemblies/scaffolds/$assemblyName\_assembly.fna \
      >> scaffolds/hgcA_scaffolds.fna
  echo "Pulling out" $scaffold "GFF entries"
  awk -v scaffold="$scaffold" '{ if ($1 == scaffold) print }' ~/BLiMMP/dataEdited/assemblies/ORFs/$assemblyName.gff \
      >> scaffolds/hgcA_scaffolds.gff
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
cd ~/BLiMMP/dataEdited/hgcA_analysis
mkdir hgcB
rm -f hgcB/downstream_gene_list.txt
python $scripts/retrieve_downstream_gene_name.py \
          identification/hgcA_good.txt \
          scaffolds/hgcA_scaffolds.gff \
          hgcB/downstream_gene_list.txt
conda deactivate

# Extract downstream amino acid sequences
rm -f hgcB/downstream_genes.faa
cat hgcB/downstream_gene_list.txt | while read gene
do
  assemblyName=$(echo $gene | cut -d "_" -f1,2)
  echo "Pulling out" $gene "faa entry from" $assemblyName
  grep -A 1 \
        $gene$ \
        $ORFs/$assemblyName.faa \
      >> hgcB/downstream_genes.faa
done
grep '>' hgcB/downstream_genes.faa | \
  sed 's/>//' > hgcB/downstream_genes_present.txt

# Search adjacent genes with hgcB HMM
cd ~/BLiMMP/dataEdited/hgcA_analysis/hgcB
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
python $scripts/convert_stockhold_to_fasta.py hgcB.sto
python $scripts/cleanFASTA.py hgcB.afa
mv -f hgcB.afa_tempCleanedFile hgcB.afa
cat hgcB.afa | while read line; do echo ${line^^} >> hgcB_clean.afa; done
sed -i 's/_ASSEMBLY/_assembly/' hgcB_clean.afa
sed -i 's/_COASSEMBLY/_coassembly/' hgcB_clean.afa
# Check sequences in alignment_hgcB.R.

# Check GFF entries for cut hgcB seqs
cd ~/BLiMMP/dataEdited/hgcA_analysis/scaffolds
cat hgcA_scaffolds.gff

# They all check out, so I'll just use the hgcB.faa
# file as our final set of hgcB seqs.
cd ~/BLiMMP/dataEdited/hgcA_analysis/hgcB
grep '>' hgcB.faa | \
  sed 's/>//' > hgcB.txt

# Pull out non-hgcB downstream genes
cd ~/BLiMMP/dataEdited/hgcA_analysis/hgcB
cat downstream_genes_present.txt hgcB.txt | \
  sort | \
  uniq -u | while read gene
  do
    grep -A 1 $gene downstream_genes.faa
  done




######################
# Isolate gene neighborhoods
######################

screen -S BLI_hgcA_gene_neighborhood
cd ~/BLiMMP/dataEdited/hgcA_analysis
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate py_viz
PYTHONPATH=''
scripts=~/BLiMMP/code

cat identification/hgcA_good.txt | while read hgcA_id
do
  scaffold_id=$(echo $hgcA_id | cut -d '_' -f 1-3)
  echo "Working on" $hgcA_id", on scaffold" $scaffold_id
  awk -F '\t' -v scaffold_id="$scaffold_id" '$1 == scaffold_id { print $0 }' scaffolds/hgcA_scaffolds.gff > scaffolds/temp_scaffolds.gff
  gene_id=$(echo $hgcA_id | \
              cut -d '_' -f 3-4 | \
              sed 's/^0*//g')
  echo "Searching for" $gene_id
  python $scripts/gene_neighborhood_extraction.py scaffolds/temp_scaffolds.gff \
                                                  scaffolds/hgcA_scaffolds.fna \
                                                  $gene_id \
                                                  5000 \
                                                  scaffolds/temp_$gene_id

  rm -f scaffolds/temp_scaffolds.gff
done

cd scaffolds
cat temp_*.gff > hgcA_geneNeighborhood_raw.gff
cat temp_*.fna > hgcA_geneNeighborhood_raw.fna
rm -f *_neighborhood.*





############################################
############################################
# Metatranscriptome counts
############################################
############################################
MT_depth=~/BLiMMP/dataEdited/metatranscriptomes/alignment
cd ~/BLiMMP/dataEdited/hgcA_analysis
mkdir MT_abundance
cat identification/hgcA_good.txt | while read hgcA_ID
do
  assemblyID=`echo $hgcA_ID | cut -d'_' -f1,2`
  echo "Pulling out MT reads for" $hgcA_ID "from" $assemblyID
  grep $hgcA_ID $MT_depth/BLI*_MT_*_to_$assemblyID\_abundance.tsv >> MT_abundance/hgcA_MT_hits.tsv
done
cd MT_abundance
sed 's/\/home\/glbrc.org\/bpeterson26\/BLiMMP\/dataEdited\/metatranscriptomes\/alignment\///' hgcA_MT_hits.tsv | \
  sed -e 's/_to_BLI2[0-1]_[a-z]*[0-9]*_abundance.tsv:/\t/' > hgcA_MT_hits_clean.tsv


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
