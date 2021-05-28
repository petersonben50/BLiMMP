#!/bin/sh

#########################
# code/metabolic_analyses/metabolic_hmms.sh
# Benjamin D. Peterson
#########################


############################################
############################################
# Identification of metabolic proteins
############################################
############################################

screen -S BLI_metabolic_HMMs

cd ~/BLiMMP/dataEdited/
mkdir metabolic_analyses
mkdir metabolic_analyses/identification

scripts=~/BLiMMP/code
metabolic_HMMs=~/BLiMMP/references/metabolic_HMMs
ORFs=~/BLiMMP/dataEdited/assemblies/ORFs
assembly_list=~/BLiMMP/metadata/assembly_list.txt
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''

tail -n +2 $metabolic_HMMs.csv | awk -F ',' '{ print $2 }' | while read HMM
do
  cd ~/BLiMMP/dataEdited/metabolic_analyses
  geneName=$(awk -F ',' -v HMM="$HMM" '$2 == HMM { print $1 }' $metabolic_HMMs.csv)

  if [ ! -e identification/$geneName.afa ]; then
    echo "Searching for" $geneName
    mkdir identification/$geneName

    awk -F '\t' '{ print $1 }' $assembly_list | while read assembly
    do
      year=$(awk -v assembly="$assembly" ' $1 == assembly { print $2 }' $assembly_list)
      if [ ! -d identification/$geneName/$year ]; then
        mkdir identification/$geneName/$year
      fi

      if [ ! -e identification/$geneName/$year/$assembly.out ]; then
        echo "Searching for" $geneName "in" $assembly
        hmmsearch --tblout identification/$geneName/$year/$assembly\_$geneName.out \
                  --cpu 4 \
                  --cut_tc \
                  $metabolic_HMMs/$HMM \
                  $ORFs/$assembly.faa \
                  > identification/$geneName/$year/$assembly\_$geneName.txt
        lineCount=`wc -l < identification/$geneName/$year/$assembly\_$geneName.out`
        if [ $lineCount -eq 13 ]; then
          echo "No" $gene "hits in" $geneName
        else
          echo "Pulling" $geneName "sequences out of" $assembly
          python $scripts/extract_protein_hitting_HMM.py \
                  identification/$geneName/$year/$assembly\_$geneName.out \
                  $ORFs/$assembly.faa \
                  identification/$geneName/$year/$assembly\_$geneName.faa
        fi
      else
        echo "Search for" $geneName "is already done in" $assembly
      fi
    done

    # Aggregate all sequences and align to HMM
    cd identification
    awk -F '\t' '{ print $2 }' ~/BLiMMP/metadata/assembly_list.txt | \
      sort | uniq | \
      while read year
      do
        shopt -s nullglob
        for i in $geneName/$year/*$geneName.faa; do FOUND=$i;break;done

        if [ ! -z $FOUND ]; then
          echo "Concatenating and aligning" $geneName
          cat $geneName/$year/*$geneName.faa \
              > $geneName\_from_$year\_all.faa
          hmmalign -o $geneName\_from_$year.sto \
                      $metabolic_HMMs/$HMM \
                      $geneName\_from_$year\_all.faa
          $scripts/convert_stockhold_to_fasta.py $geneName\_from_$year.sto
          grep '>' $geneName\_from_$year\_all.faa | \
            sed 's/>//' \
            > $geneName\_from_$year\_all_list.txt
        else
          echo "No" $geneName "sequences found at all :("
        fi
    done
    FOUND=""
  else
    echo "Already pulled out" $geneName "sequences"
  fi
done

conda deactivate



############################################
############################################
# Dereplicate sequences
############################################
############################################

cd ~/BLiMMP/dataEdited/metabolic_analyses
mkdir dereplication
cdhit=~/programs/cdhit-master

echo "geneID,geneName" > metabolic_gene_key.csv

tail -n +2 $metabolic_HMMs.csv | awk -F ',' '{ print $1 }' | while read geneName
do
  cd ~/BLiMMP/dataEdited/metabolic_analyses
    awk -F '\t' '{ print $2 }' ~/BLiMMP/metadata/assembly_list.txt | \
      sort | uniq | \
      while read year
      do
        if [ -e identification/$geneName\_from_$year\_all.faa ]; then
          echo "Dereplicating" $geneName
          $cdhit/cd-hit -g 1 \
                        -i identification/$geneName\_from_$year\_all.faa \
                        -o dereplication/$geneName\_from_$year.faa \
                        -c 0.97 \
                        -n 5 \
                        -d 0
          $cdhit/clstr2txt.pl dereplication/$geneName\_from_$year.faa.clstr \
            > dereplication/$geneName\_from_$year.tsv
        else
          echo "There were no hits for" $geneName
        fi
    done

  cd dereplication
  if [ -e $geneName\_from_*.faa ]; then

    cat $geneName\_from_*.faa > $geneName\_derep.faa
    grep '>' $geneName\_derep.faa | \
      sed 's/>//' \
      > $geneName\_derep_list.txt
    awk -v geneName="$geneName" '{ print $0","geneName }' $geneName\_derep_list.txt >> ../metabolic_gene_key.csv
  else
    echo "Still no hits for" $geneName
  fi
done


############################################
############################################
# Classify dsrA genes
############################################
############################################

screen -S BLI_dsrA_tree
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=""
PERL5LIB=''
scripts=~/BLiMMP/code
# Set up directory
cd ~/BLiMMP/dataEdited/metabolic_analyses
mkdir sulfur
mkdir sulfur/dsrA

#########################
# Generate dsrA alignment
#########################

# Copy over my sequences and align them
sed 's/*//' dereplication/dsrA_derep.faa > sulfur/dsrA/dsrA.faa
grep '>' sulfur/dsrA/dsrA.faa | \
  sed 's/>//' \
  > sulfur/dsrA/dsrA_list.txt
metabolic_HMMs=~/BLiMMP/references/metabolic_HMMs
cd sulfur/dsrA
hmmalign -o dsrA.sto \
            $metabolic_HMMs/TIGR02064.HMM \
            dsrA.faa
$scripts/convert_stockhold_to_fasta.py dsrA.sto

# Copy in reference sequences
cp ~/references/metabolicProteins/sulfur/dsrA/dsrA_karthik_clean.afa \
    dsrA_karthik_clean.afa

# Align seqs
muscle -profile \
        -in1 dsrA_karthik_clean.afa \
        -in2 dsrA.afa \
        -out dsrA_phylogeny.afa
# Trim alignment
trimal -in dsrA_phylogeny.afa \
        -out dsrA_phylogeny_trimmed.afa \
        -gt 0.5


#########################
# Generate ML tree
#########################

cd ~/BLiMMP/dataEdited/metabolic_analyses/sulfur/dsrA/
python $scripts/cleanFASTA.py dsrA_phylogeny_trimmed.afa
mv -f dsrA_phylogeny_trimmed.afa_temp.fasta dsrA_phylogeny_trimmed.afa
FastTree dsrA_phylogeny_trimmed.afa \
    > dsrA_phylogeny_trimmed.tree


############################################
############################################
# Identify potential PCCs
############################################
############################################

#########################
# First identify MHCs
#########################

screen -S BLI_MHCs
cd ~/BLiMMP/dataEdited/metabolic_analyses
mkdir MHCs
cd MHCs
scripts=~/BLiMMP/code
assembly_list=~/BLiMMP/metadata/assembly_list.txt
ORFs=~/BLiMMP/dataEdited/assemblies/ORFs
chmod +x $scripts/Find_multiheme_protein.py
# Let's look for proteins with at least 3 haem-binding sites
awk -F '\t' '{ print $1 }' $assembly_list | while read assembly
do
  year=$(awk -v assembly="$assembly" ' $1 == assembly { print $2 }' $assembly_list)
  if [ ! -d $year ]; then
    mkdir $year
  fi
  echo "Searching for MHCs in" $assembly
  $scripts/Find_multiheme_protein.py $ORFs/$assembly.faa 3
  # Move the scripts to the correct directory
  mv $ORFs/$assembly\_3_heme* $year/
done

cat */*3_heme_list* > MHC_list.txt


#########################
# Pull out adjacent genes
#########################

# Pull out names of adjacent genes
cd ~/BLiMMP/dataEdited/metabolic_analyses
mkdir BBOMP
rm -f BBOMP/adjacent_genes_all_list.txt
ORFs=~/BLiMMP/dataEdited/assemblies/ORFs
cat MHCs/MHC_list.txt | while read gene
do
  echo "Working on" $gene
  scaffold=$(echo $gene | rev | cut -d"_" -f2- | rev)
  assembly=$(echo $gene | rev | cut -d"_" -f3- | rev)
  ORFnumber=$(echo $gene | rev | cut -d"_" -f1 | rev)
  preceedingORFnumber=$(expr $ORFnumber - 1)
  followingORFnumber=$(expr $ORFnumber + 1)
  echo $scaffold"_"$preceedingORFnumber >> BBOMP/adjacent_genes_all_list.txt
  echo $scaffold"_"$followingORFnumber >> BBOMP/adjacent_genes_all_list.txt
done

# Find unique gene names
wc -l BBOMP/adjacent_genes_all_list.txt
sort BBOMP/adjacent_genes_all_list.txt | \
  uniq \
  > BBOMP/adjacent_genes_unique_list.txt

# Pull out adjacent genes
rm -f BBOMP/adjacent_genes.faa
cat BBOMP/adjacent_genes_unique_list.txt | while read geneID
do
  assembly=$(echo $geneID | rev | cut -d"_" -f3- | rev)
  echo "Looking for" $geneID "in" $assembly
  grep -A 1 -m 1 $geneID$ $ORFs/$assembly.faa >> BBOMP/adjacent_genes.faa
done


#########################
# Search adjacent genes for BBOMP
#########################
cd ~/BLiMMP/dataEdited/metabolic_analyses/BBOMP
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=""
PERL5LIB=""

# Run the BBOMP HMM
pcc_omp_HMM=~/BLiMMP/references/metabolicProteins/EET/pcc_omp.HMM
hmmsearch --tblout pcc_omp_custom.out \
          -T 50 \
          $pcc_omp_HMM \
          adjacent_genes.faa \
          > pcc_omp_custom_output.txt
# Pull out the sequences of interest
scripts=~/BLiMMP/code
python $scripts/extract_protein_hitting_HMM.py \
        pcc_omp_custom.out \
        adjacent_genes.faa \
        pcc_omp_custom.faa
# Dereplicate sequences
cd-hit -i pcc_omp_custom.faa \
        -o pcc_omp_custom_derep.faa \
        -g 1 \
        -c 0.97 \
        -n 5 \
        -d 0
clstr2txt.pl pcc_omp_custom_derep.faa.clstr > pcc_omp_custom_derep.tsv
grep '>' pcc_omp_custom_derep.faa | sed 's/>//' > pcc_omp_custom_derep_list.txt

# Align sequences to HMM
hmmalign -o pcc_omp_custom_derep.sto \
            $pcc_omp_HMM \
            pcc_omp_custom_derep.faa
# Convert alignment to fasta format
scripts=~/BLiMMP/code
$scripts/convert_stockhold_to_fasta.py pcc_omp_custom_derep.sto
python $scripts/cleanFASTA.py pcc_omp_custom_derep.afa
mv -f pcc_omp_custom_derep.afa_temp.fasta pcc_omp_custom_derep.afa


#########################
# Prepare reference sequences
#########################
cd ~/BLiMMP/dataEdited/metabolic_analyses/BBOMP
mkdir references
# Find references from RefSeq
blastp -query pcc_omp_custom_derep.faa \
        -db ~/references/ncbi_db/refseq/refseq_protein \
        -evalue 0.005 \
        -outfmt '6 sseqid sseq' \
        -max_target_seqs 5 \
        -num_threads 3 \
        -out references/refseq_pcc_omp.txt

# Pull out amino acid sequences and dereplicate them
cd references
awk -F '\t' '{ print ">"$1"\n"$2 }' refseq_pcc_omp.txt > refseq_pcc_omp.faa
sed -i 's/ref|//' refseq_pcc_omp.faa
sed -i 's/|//' refseq_pcc_omp.faa
sed -i 's/-//g' refseq_pcc_omp.faa
cd-hit -g 1 \
        -i refseq_pcc_omp.faa \
        -o refseq_pcc_omp_derep.faa \
        -c 0.97 \
        -n 5 \
        -d 0
rm -f refseq_pcc_omp_clean.faa
grep '>' refseq_pcc_omp_derep.faa | sort | uniq | while read fastaHeader
do
  grep -A 1 -m 1 $fastaHeader refseq_pcc_omp_derep.faa >> refseq_pcc_omp_clean.faa
done

# Dereplicate refseq sequences against reference dataset
cd-hit-2d -i ~/BLiMMP/references/metabolicProteins/EET/pcc_omp.faa \
          -i2 refseq_pcc_omp_clean.faa \
          -o blast_pcc_omp_uniq.faa \
          -c 0.97 \
          -n 5

# Get final list of refseq references
grep '>' blast_pcc_omp_uniq.faa | \
  sed 's/>//' \
  > blast_pcc_omp_uniq_list.txt

# Retrieve reference information on local computer
BLiMMP
cd dataEdited/metabolic_analyses/PCC
epost -db protein -input blast_pcc_omp_uniq_list.txt | \
    esummary | \
    xtract -pattern DocumentSummary -element AccessionVersion,Organism > refseq_bbomp_metadata.tsv


#########################
# Generate a phylogeny
#########################
# Concatenate sequences
cd ~/BLiMMP/dataEdited/metabolic_analyses/BBOMP
mkdir phylogeny
cat references/blast_pcc_omp_uniq.faa \
    ~/BLiMMP/references/metabolicProteins/EET/pcc_omp.faa \
    pcc_omp_custom_derep.faa \
    > phylogeny/bbomp_phylogeny.faa
# Generate alignment
cd phylogeny
muscle -in bbomp_phylogeny.faa \
        -out bbomp_phylogeny.afa
# Trim the alignment
trimal -in bbomp_phylogeny.afa \
        -out bbomp_phylogeny_trimmed.afa \
        -gt 0.5
FastTree bbomp_phylogeny_trimmed.afa > bbomp.tree


#########################
# Add sequence IDs to metabolic gene key
#########################
cd ~/BLiMMP/dataEdited/metabolic_analyses/BBOMP
cat pcc_omp_custom_derep_list.txt | while read geneID
do
  echo $geneID",BBOMP" >> ~/BLiMMP/dataEdited/metabolic_analyses/metabolic_gene_key.csv
done

############################################
############################################
# Check nitrate reductase genes
############################################
############################################

# Get set up
mkdir ~/BLiMMP/references/custom_hmms
cp ~/5M/references/MoORs/MoOR.HMM ~/BLiMMP/references/custom_hmms
cp ~/5M/references/MoORs/MoOR_final_reference.afa ~/BLiMMP/references/custom_hmms/MoOR_reference.afa
screen -S BLI_MoORs
cd ~/BLiMMP/dataEdited/metabolic_analyses
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''
scripts=~/BLiMMP/code

# Copy in nitrate reductase sequences
mkdir MoORs
cat dereplication/narG_derep.faa \
    dereplication/napA_derep.faa \
    > MoORs/nitrate_reductases.faa
sed -i 's/*//' MoORs/nitrate_reductases.faa
# Generate alignment of nitrate reductases
cd MoORs
hmmalign -o nitrate_reductases.sto \
          ~/BLiMMP/references/custom_hmms/MoOR.HMM \
          nitrate_reductases.faa
$scripts/convert_stockhold_to_fasta.py nitrate_reductases.sto

# Generate consensus alignment
muscle -profile \
        -in1 nitrate_reductases.afa \
        -in2 ~/BLiMMP/references/custom_hmms/MoOR_reference.afa \
        -out nitrate_reductases_references.afa

python $scripts/cleanFASTA.py nitrate_reductases_references.afa
mv -f nitrate_reductases_references.afa_temp.fasta nitrate_reductases_references.afa

# Mask the alignment at 50% gaps
trimal -in nitrate_reductases_references.afa \
        -out nitrate_reductases_references_masked.afa \
        -gt 0.5
FastTree nitrate_reductases_references_masked.afa > nitrate_reductases_moor.tree


############################################
############################################
# Extract depths of all scaffolds
############################################
############################################

# Generate list of all scaffolds that have potential metabolic genes on them
cd ~/BLiMMP/dataEdited/metabolic_analyses
mkdir depth
tail -n +2 metabolic_gene_key.csv | \
  cut -d"_" -f1-3 | \
  sort | \
  uniq \
  > depth/scaffold_all_list.txt

# Pull out all depths
mkdir ~/BLiMMP/reports/aggregate_depth_proteins
cd ~/BLiMMP/reports/aggregate_depth_proteins
mkdir outs errs logs
cd ~/BLiMMP/dataEdited/metabolic_analyses/depth/
awk -F '\t' '{ print $1 }' ~/BLiMMP/metadata/metagenome_list.txt > metagenomes_to_use.txt

chmod +x /home/GLBRCORG/bpeterson26/BLiMMP/code/aggregate_depth_proteins.sh
condor_submit /home/GLBRCORG/bpeterson26/BLiMMP/code/aggregate_depth_proteins.sub
