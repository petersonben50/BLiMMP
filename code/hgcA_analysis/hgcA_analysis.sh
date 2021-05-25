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
# Identify hgcA sequences
############################################
############################################

######################
# Identify putative hgcA genes with HMM
######################

screen -S BLI_hgcA_search
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''
scripts=~/BLiMMP/code
ORFs=~/BLiMMP/dataEdited/assemblies/ORFs
mkdir ~/BLiMMP/dataEdited/hgcA_analysis
mkdir ~/BLiMMP/dataEdited/hgcA_analysis/identification
cd ~/BLiMMP/dataEdited/hgcA_analysis/identification
mkdir 2020
cd ~/BLiMMP/dataEdited

cat ~/BLiMMP/metadata/assembly_list.txt | while read line
do
  assembly=`echo $line | awk -F '\t' '{ print $1 }'`
  year=`echo $line | awk -F '\t' '{ print $2 }'`
  echo "Searching for hgcA in" $assembly "from" $year
  if [ -e $ORFs/$assembly.faa ]; then
    if [ ! -e hgcA_analysis/identification/$year/$assembly\_hgcA_report.txt ]; then
      echo "Search for hgcA in" $assembly
      hmmsearch --tblout hgcA_analysis/identification/$year/$assembly\_hgcA.out \
                --cpu 4 \
                --cut_tc \
                ~/references/hgcA/hgcA.hmm \
                $ORFs/$assembly.faa \
                > hgcA_analysis/identification/$year/$assembly\_hgcA_report.txt
      python $scripts/extract_protein_hitting_HMM.py \
              hgcA_analysis/identification/$year/$assembly\_hgcA.out \
              $ORFs/$assembly.faa \
              hgcA_analysis/identification/$year/$assembly\_hgcA.faa
    else
      echo "Search is already done in" $assembly
    fi
  else
    echo "Genes aren't predicted for" $assembly
  fi
done

######################
# Concatenate and align all hgcA seqs for curation
######################
chmod +x /home/GLBRCORG/bpeterson26/BLiMMP/code/convert_stockhold_to_fasta.py
cd ~/BLiMMP/dataEdited/hgcA_analysis/identification/
cat 202*/*_hgcA.faa > hgcA_raw.faa
# Align seqs
hmmalign -o hgcA_raw.sto \
            ~/references/hgcA/hgcA.hmm \
            hgcA_raw.faa
$scripts/convert_stockhold_to_fasta.py hgcA_raw.sto

# Curate the hgcA list
# See notes in md file.
# Then:
cd ~/BLiMMP/dataEdited/hgcA_analysis/identification/
grep '>' hgcA_good.afa | \
    sed 's/>//' \
    > hgcA_good.txt
sed 's/-//g' hgcA_good.afa > hgcA_good.faa






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
scripts=~/BLiMMP/code
mkdir $workingDirectory/classification

# Generate fasta file of reference alignment
cp -avr ~/Everglades/references/hgcA/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.refpkg $references/
cd $references/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.refpkg
$scripts/convert_stockhold_to_fasta.py Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.stockholm
mv Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.afackholm $workingDirectory/classification/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.afa

# Generate alignment of sequences of interest
cd $workingDirectory
muscle -in identification/hgcA_good.faa \
        -out classification/hgcA_muscle.afa
cd classification/

# Combine the alignments of seqs from this study and references
muscle -profile -in1 hgcA_muscle.afa \
        -in2 Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.afa \
        -out hgcA_for_classification.afa
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

# Search adjacent genes with hgcB HMM
cd ~/BLiMMP/dataEdited/hgcA_analysis
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''
scripts=~/BLiMMP/code
hmmsearch --tblout hgcB/hgcB.out \
          --cpu 4 \
          -T 30 \
          ~/references/hgcA/hgcB_5M.HMM \
          hgcB/downstream_genes.faa \
          > hgcB/hgcB_report.txt

# Check the hgcB hits
rm -f hgcB/hgcB.faa
grep -v "#" hgcB/hgcB.out | awk '{ print $1 }' | while read geneID
do
  grep -A 1 $geneID$ hgcB/downstream_genes.faa >> hgcB/hgcB.faa
done
grep '>' hgcB/hgcB.faa | wc -l
# Align sequences to HMM
hmmalign -o hgcB/hgcB.sto \
            ~/references/hgcA/hgcB_5M.HMM \
            hgcB/hgcB.faa
# Convert alignment to fasta format
$scripts/convert_stockhold_to_fasta.py hgcB/hgcB.sto


# Check sequences in Geneious.

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
cat downstream_gene_list.txt hgcB.txt | \
  sort | \
  uniq -u
grep -A 1 'fall2017cluster6_000000000428_49' downstream_genes.faa
grep -A 1 'fall2017coassembly_000000448578_3' downstream_genes.faa
grep -A 1 'fall2017coassembly_000001334838_1' downstream_genes.faa
grep -A 1 'HC18HY300_000000013755_6' downstream_genes.faa
grep -A 1 'HC18ME02_000000018262_3' downstream_genes.faa
grep -A 1 'KMBP004F_000000216801_1' downstream_genes.faa
grep -A 1 'KMBP009B_000000084934_2' downstream_genes.faa




######################
# Isolate gene neighborhoods
######################

screen -S BLI_hgcA_gene_neighborhood
cd ~/BLiMMP/dataEdited/hgcA_analysis
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate py_viz
PYTHONPATH=''
IFS=$'\n'
scripts=~/BLiMMP/code/generalUse

grep '>' identification/hgcA_raw.faa | \
    sed 's/>//' | \
    while read hgcA_id
    do
      echo "Working on" $hgcA_id
      scaffold_id=$(echo $hgcA_id | cut -d '_' -f 1-2)
      awk -F '\t' -v scaffold_id="$scaffold_id" '$1 == scaffold_id { print $0 }' scaffolds/hgcA_scaffolds.gff > scaffolds/temp_scaffolds.gff
      gene_id=$(echo $hgcA_id | \
                  cut -d '_' -f 2-3 | \
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
rm -f hgcA_geneNeighborhood.gff hgcA_geneNeighborhood.fna
cat temp_*.gff > hgcA_geneNeighborhood_raw.gff
cat temp_*.fna > hgcA_geneNeighborhood_raw.fna
rm -f *_neighborhood.*
