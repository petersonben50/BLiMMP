#!/bin/sh

##########################
# code/binning/binning_workflow.sh
# Benjamin D. Peterson

# This workflow will generate our set of
# hgcA+ genomes, as well as a set of
# uncurated bins.
##########################

####################################################
####################################################
# Prepare scaffolds and mapping files
####################################################
####################################################


##########################
# Filter out short scaffolds
##########################

screen -S BLI_binning
mkdir ~/BLiMMP/dataEdited/binning
cd ~/BLiMMP/dataEdited
mkdir binning/scaffolds
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate anvio6.2
PYTHONPATH=''
assembly_list=~/BLiMMP/metadata/assembly_list.txt


awk -F '\t' '{ print $1 }' $assembly_list | while read assembly
do
  if [ ! -e binning/scaffolds/$assembly\_filtered_scaffolds.fna ]; then
    echo "Processing" $assembly "scaffolds for binning"
    anvi-script-reformat-fasta assemblies/scaffolds/$assembly\_assembly.fna \
                              -o binning/scaffolds/$assembly\_filtered_scaffolds.fna \
                              -l 2000
  else
    echo $assembly "scaffolds already processed for binning"
  fi
done


##########################
# Map reads to filtered scaffolds
##########################
cd /home/GLBRCORG/bpeterson26/BLiMMP/reports/
mkdir binning_mapping
cd binning_mapping
rm -f outs/*_binningMapping.out \
      errs/*_binningMapping.err \
      logs/*_binningMapping.log
mkdir outs errs logs

cd /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/binning
mkdir mapping
mkdir mapping/indices
awk -F '\t' '{ print $1 }' ~/BLiMMP/metadata/assembly_list.txt > assembly_list_to_use.txt

cd /home/GLBRCORG/bpeterson26/BLiMMP/code/
chmod +x binning_mapping.sh
condor_submit binning_mapping.sub


####################################################
####################################################
# Run automatic binning algorithms
####################################################
####################################################
mkdir ~/BLiMMP/dataEdited/binning/autoBinning
cd ~/BLiMMP/dataEdited/binning/autoBinning
mkdir metabat2 maxbin2 dasTool finalBins
cd ~/BLiMMP/reports
mkdir autoBinning
cd autoBinning
mkdir outs errs logs

cd /home/GLBRCORG/bpeterson26/BLiMMP/code/
chmod +x automatic_binning.sh
condor_submit automatic_binning.sub



####################################################
####################################################
# Prepare anvi'o databases for manual binning
####################################################
####################################################
mkdir /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/binning/manualBinning
mkdir /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/binning/manualBinning/anvioDBs
mkdir ~/BLiMMP/reports/anvioDBprep
cd ~/BLiMMP/reports/anvioDBprep
mkdir outs errs logs


##########################
# Generate contig databases
##########################
cd /home/GLBRCORG/bpeterson26/BLiMMP/code/
chmod +x anvio_DB_prep.sh
condor_submit anvio_DB_prep.sub


##########################
# Generate read profiles
##########################
mkdir ~/BLiMMP/reports/anvioProfiling
cd ~/BLiMMP/reports/anvioProfiling
mkdir outs errs logs
cd /home/GLBRCORG/bpeterson26/BLiMMP/code/
chmod +x anvio_profiling.sh
condor_submit anvio_profiling.sub


##########################
# Run CONCOCT binning
##########################
mkdir /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/binning/manualBinning/anvioDBs/original_summaries
mkdir ~/BLiMMP/reports/concoctBinning
cd ~/BLiMMP/reports/concoctBinning
mkdir outs errs logs
cd /home/GLBRCORG/bpeterson26/BLiMMP/code/
chmod +x generate_large_bin_clusters_anvio.sh
condor_submit generate_large_bin_clusters_anvio.sub


##########################
# Search bins for hgcA
##########################
screen -S BLI_binning
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate anvio6.2
PYTHONPATH=''
assembly_list=~/BLiMMP/metadata/assembly_list.txt
cd /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/binning/manualBinning/anvioDBs
mkdir original_summaries/hgcA_search
echo -e "assembly\tbin" > original_summaries/hgcA_search/original_hgcA_bin_list.txt

awk -F '\t' '{ print $1 }' $assembly_list | while read assembly
do
  if [ ! -e original_summaries/uncurated_bins_from_$assembly ]; then
    echo "Summarizing binning from" $assembly
    anvi-summarize -c $assembly.db \
                    -p $assembly.merged/PROFILE.db \
                    -C CONCOCT \
                    -o original_summaries/uncurated_bins_from_$assembly

    ls original_summaries/uncurated_bins_from_$assembly/bin_by_bin | sed 's/\///' \
        > original_summaries/$assembly\_original_bin_list.txt
    cat original_summaries/$assembly\_original_bin_list.txt | while read bin
    do
      if [ -s original_summaries/uncurated_bins_from_$assembly/bin_by_bin/$bin/$bin-hgcaAnvio-hmm-sequences.txt ]; then
        echo $assembly$'\t'$bin >> original_summaries/hgcA_search/original_hgcA_bin_list.txt
      fi
    done
  else
    echo $assembly "already summarized"
  fi
done



####################################################
####################################################
# Add automatic binning information into anvio databases
####################################################
####################################################

screen -S BLI_anvioDBs_add_autoBins
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate anvio6.2
PYTHONPATH=""

# Set up variables
metabat2=~/BLiMMP/dataEdited/binning/autoBinning/metabat2
maxbin2=~/BLiMMP/dataEdited/binning/autoBinning/maxbin2
dasTool=~/BLiMMP/dataEdited/binning/autoBinning/finalBins
anvioDB=~/BLiMMP/dataEdited/binning/manualBinning/anvioDBs
assembly_list=~/BLiMMP/metadata/assembly_list.txt
mkdir $anvioDB/S2B_files
mkdir $anvioDB/binning_collections

# Add in automatic binning data to anvio DB.
awk -F '\t' '{ print $1 }' $assembly_list | while read assembly
do
  echo "Adding scaffold to bin files for" $assembly
  # Copy all S2B files to one folder
  cp $metabat2/$assembly\_metabat_S2B.tsv \
      $maxbin2/$assembly\_maxbin_S2B.tsv \
      $dasTool/$assembly\_dasTool_S2B.tsv \
      $anvioDB/S2B_files
  anvi-import-collection $anvioDB/S2B_files/$assembly\_metabat_S2B.tsv \
                           -c $anvioDB/$assembly.db \
                           -p $anvioDB/$assembly.merged/PROFILE.db \
                           -C metabat2 \
                           --contigs-mode
  anvi-import-collection $anvioDB/S2B_files/$assembly\_maxbin_S2B.tsv \
                          -c $anvioDB/$assembly.db \
                          -p $anvioDB/$assembly.merged/PROFILE.db \
                          -C maxbin2 \
                          --contigs-mode
  anvi-import-collection $anvioDB/S2B_files/$assembly\_dasTool_S2B.tsv \
                          -c $anvioDB/$assembly.db \
                          -p $anvioDB/$assembly.merged/PROFILE.db \
                          -C dasTool \
                          --contigs-mode
  anvi-script-merge-collections -c $anvioDB/$assembly.db \
                                -i $anvioDB/S2B_files/$assembly*S2B.tsv \
                                -o $anvioDB/binning_collections/$assembly\_collections.tsv
done


################################################
################################################
# Manually bin hgcA+ bins
################################################
################################################
screen -S BLI_anvioDBs_binning
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate anvio6.2
PYTHONPATH=""
cd ~/BLiMMP/dataEdited/binning/manualBinning/
# Copy the database folder
#cp -avr anvioDBs anvioDBs_modified
#cat /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/binning/manualBinning/anvioDBs/original_summaries/hgcA_search/original_hgcA_bin_list.txt
assembly=BLI20_coassembly
bin=Bin_6
anvi-refine -p anvioDBs_modified/$assembly.merged/PROFILE.db \
            -c anvioDBs_modified/$assembly.db \
            -C CONCOCT \
            -b $bin \
            -A anvioDBs_modified/binning_collections/$assembly\_collections.tsv \
            --taxonomic-level "t_phylum"


################################################
################################################
# Summarize/export curated bins
################################################
################################################

##########################
# Rename and summarize bins
##########################
screen -S BLI_binsPostProcessing
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate anvio6.2
PYTHONPATH=""
PERL5LIB=""
assembly_list=~/BLiMMP/metadata/assembly_list.txt
cd ~/BLiMMP/dataEdited/binning/manualBinning/anvioDBs_modified

# Summarize them
awk -F '\t' '{ print $1 }' $assembly_list | while read assembly
do
  # If old summary exists that we want to delete, uncomment
  # the following:
  # rm -r $assembly.summary.curated
  if [ ! -d $assembly.curated.summary ]; then
    echo "Summarizing bins for" $assembly
    anvi-summarize -c $assembly.db \
                    -p $assembly.merged/PROFILE.db \
                    -C CONCOCT \
                    -o $assembly.curated.summary
  else
    echo "We already summarized the curated bins for" $assembly
  fi
done
conda deactivate


##########################
# Pull out DNA files from hgcA+ bins
##########################
# First need to set up new directory
cd ~/BLiMMP/dataEdited/binning/manualBinning/
mkdir binsRaw
mkdir binsRaw/DNA
binsRaw=~/BLiMMP/dataEdited/binning/manualBinning/binsRaw

awk -F '\t' '{ print $1 }' $assembly_list | while read assembly
do
  binSummary=~/BLiMMP/dataEdited/binning/manualBinning/anvioDBs_modified/$assembly.curated.summary
  if [ -e $binSummary ]; then
    if [ ! -e $binsRaw/DNA/$assembly* ]; then
      cd $binSummary/bin_by_bin
      ls | sed 's/\///' | while read bin
      do
        isThereHgcA=`cat $bin/$bin\-hgcaAnvio-hmm-sequences.txt | wc -l`
        if [ ! $isThereHgcA -eq 0 ]; then
          echo "Copying" $bin "to binsRaw folder"
          cp $bin/$bin-contigs.fa $binsRaw/DNA/$bin.fna
        else
          echo "No hgcA in" $bin
        fi
      done
    else
      echo "Hey, there are some bins from" $assembly "already in here"
      echo "You might wanna check that out before you start overwriting stuff"
    fi
  else
    echo "Summarize anvioDB for" $assembly", dummy."
  fi
done

# Generate list of hgcA+ bins
cd $binsRaw/DNA
ls *.fna | \
  sed 's/.fna//' \
  > binsRaw_hgcA_list.txt


####################################################
####################################################
# Check quality of bins
####################################################
####################################################

##########################
# Completeness/redundancy estimates from anvio
##########################
binsRaw=~/BLiMMP/dataEdited/binning/manualBinning/binsRaw
mkdir $binsRaw/anvio_data
mkdir $binsRaw/anvio_data/completeness_redundancy

# Copy summary files into a single folder.
awk -F '\t' '{ print $1 }' $assembly_list | while read assembly
do
  binSummary=~/BLiMMP/dataEdited/binning/manualBinning/anvioDBs_modified/$assembly.curated.summary
  if [ -e $binSummary/bins_summary.txt ]; then
    cp $binSummary/bins_summary.txt $binsRaw/anvio_data/completeness_redundancy/$assembly\_bins_summary.txt
  else
    echo $assembly "has not been summarized."
  fi
done

# Concatenate summaries into a single file.
cd $binsRaw/anvio_data/completeness_redundancy
head -n 1 fall2017cluster1_bins_summary.txt > bins_summary_all.txt
ls *_bins_summary.txt | while read file
do
  tail -n +2 $file >> bins_summary_all.txt
done

# Only keep the summaries for the hgcA+ bins.
rm -f bins_summary_hgcA.txt
cat $binsRaw/DNA/binsRaw_hgcA_list.txt | while read hgcA_bin
do
  grep $hgcA_bin bins_summary_all.txt >> bins_summary_hgcA.txt
done

# Which bins have greater than 50% C and less than 10% R
head -n 1 bins_summary_hgcA.txt > bins_summary_hgcA_good.txt
awk '{ if (($7 > 50) && ($8 < 10)) print $0 }' bins_summary_hgcA.txt >> bins_summary_hgcA_good.txt
tail -n +2 bins_summary_hgcA_good.txt | \
  awk '{ print $1 }' \
  > bins_list_hgcA_good.txt




##########################
# Completeness/redundancy estimates from CheckM
##########################

screen -S BLI_checkM
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
PYTHONPATH=""
conda activate bioinformatics
binsRaw=~/BLiMMP/dataEdited/binning/manualBinning/binsRaw

cd $binsRaw
if [ -d checkM ]; then
  echo "Removing old checkM folder"
  rm -rf checkM
fi
mkdir checkM
checkm lineage_wf \
      -x .fna \
      -t 16 \
      DNA \
      checkM
checkm qa checkM/lineage.ms \
        checkM \
        -o 2 \
        -f checkM/checkm.out \
        --tab_table
awk -F '\t' \
  -v OFS=',' \
  '{ print $1,$6,$7,$8,$9,$11,$13,$15,$17,$19,$23 }' \
  checkM/checkm.out \
  > checkM/checkM_stats.csv
# Download checkM/checkM_stats.csv to local computer:
# dataEdited/binning/rawBins/bin_quality

cd ~/BLiMMP/dataEdited/binning/manualBinning
mkdir binsGood
mkdir binsGood/DNA
mkdir binsGood/checkM
awk -F ',' '{ if (($2 > 50) && ($3 < 10)) print $0 }' \
  binsRaw/checkM/checkM_stats.csv \
  > binsGood/checkM/good_bins_data.txt
awk -F ',' '{ print $1 }' binsGood/checkM/good_bins_data.txt \
  > binsGood/checkM/good_bins_list.txt


##########################
# Final bin list
##########################
cd ~/BLiMMP/dataEdited/binning/manualBinning
cat binsGood/checkM/good_bins_list.txt binsRaw/anvio_data/completeness_redundancy/bins_list_hgcA_good.txt | \
  sort | uniq \
  > goodBins_list.txt

# Copy over good bins
cat goodBins_list.txt | while read binsGood
do
  echo "Copying" $binsGood
  cp binsRaw/DNA/$binsGood.fna binsGood/DNA
done

cd binsGood/DNA
scripts=~/BLiMMP/code
ls *fna | while read fna
do
  echo "Cleaning" $fna
  python $scripts/cleanFASTA.py $fna
  mv -f $fna\_temp.fasta $fna
done


####################################################
####################################################
# Check out taxonomy of bins with GTDB
####################################################
####################################################

screen -S BLI_GTDB
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
PYTHONPATH=""
conda activate gtdbtk
cd ~/BLiMMP/dataEdited/binning/manualBinning/binsGood
rm -rf taxonomy
mkdir taxonomy

gtdbtk classify_wf \
        --cpus 16 \
        --extension fna \
        --genome_dir ./DNA \
        --out_dir taxonomy
# Summarize them
cd taxonomy
awk -F '\t' '{ print $1"\t"$2 }' gtdbtk.*.summary.tsv \
        > taxonomy_summary.txt
conda deactivate



####################################################
####################################################
# Pull out coverage information from anvio
####################################################
####################################################

# Coverage data from anvio
binsGood=~/BLiMMP/dataEdited/binning/manualBinning/binsGood
assembly_list=~/BLiMMP/metadata/assembly_list.txt
mkdir $binsGood/coverageAnvio
cd $binsGood/coverageAnvio

# Pull out coverage information
awk -F '\t' '{ print $1 }' $assembly_list | while read assembly
do
  summary=~/BLiMMP/dataEdited/binning/manualBinning/anvioDBs_modified/$assembly.curated.summary/bins_across_samples
  if [ -e $summary/mean_coverage_Q2Q3.txt ]; then
    if [ ! -e $binsGood/coverageAnvio/$assembly\_coverage.txt ]; then
      cp -f $summary/mean_coverage_Q2Q3.txt $binsGood/coverageAnvio/$assembly\_coverage.txt
    fi
  else
    echo $assembly "has not been summarized."
  fi
done

# Combine coverage information
awk -F '\t' '{ print $2 }' $assembly_list | \
  sort | uniq | \
  while read year
do
  initialAssembly=`awk -F '\t' -v year="$year" '$2 == year { print $1 }' $assembly_list | head -n 1`
  head -n 1 $initialAssembly\_coverage.txt > coverage_$year.txt

  awk -F '\t' -v year="$year" '$2 == year { print $1 }' $assembly_list | while read assembly
  do
    tail -n +2 $assembly\_coverage.txt >> coverage_$year.txt
  done

  head -n 1 coverage_$year.txt > coverage_goodBins_$year.txt
  grep -f ~/BLiMMP/dataEdited/binning/manualBinning/goodBins_list.txt coverage_$year.txt >> coverage_goodBins_$year.txt
  rm -f coverage_$year.txt
done
rm -f *_coverage.txt
# Download the "coverage_goodBins_202*.txt" files to my computer



####################################################
####################################################
# Get ORFs for bins
####################################################
####################################################

screen -S BLI_binsORFS
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=""
scripts=/home/GLBRCORG/bpeterson26/BLiMMP/code
binsGood=~/BLiMMP/dataEdited/binning/manualBinning/binsGood
binsGoodList=~/BLiMMP/dataEdited/binning/manualBinning/goodBins_list.txt
cd $binsGood
mkdir ORFs

cat $binsGoodList | while read bin
do
  if [ ! -e ORFs/$bin.gff ]; then
    echo "Predicting proteins for" $bin
    prodigal -i DNA/$bin.fna \
              -o ORFs/$bin.gff \
              -f gff \
              -a ORFs/$bin.faa \
              -d ORFs/$bin.fna \
              -p single
  else
    echo $bin "already run."
  fi
done

# Clean up the gene sequence files
cd ORFs
cat $binsGoodList | while read bin
do
  echo "Cleaning" $bin
  python $scripts/cleanFASTA.py $bin.fna
  mv -f $bin.fna_temp.fasta $bin.fna
  python $scripts/cleanFASTA.py $bin.faa
  mv -f $bin.faa_temp.fasta $bin.faa
done

# Combine gene sequence files and generate S2B/G2B files
cd $binsGood
$scripts/Fasta_to_Scaffolds2Bin.sh -e fna \
                                    -i DNA \
                                    > binsGood_S2B.tsv
cat DNA/*.fna > DNA.fna
$scripts/Fasta_to_Scaffolds2Bin.sh -e faa \
                                    -i ORFs \
                                    > binsGood_G2B.tsv
cat ORFs/*.faa > ORFs.faa



####################################################
####################################################
# Run ANI comparisons on good bins
####################################################
####################################################
# The following workflow is the code to run
# Sarah Stevens's ANI calculator on a
# folder full of bins.
# Details: https://github.com/sstevens2/ani_compare_dag


mkdir ~/BLiMMP/dataEdited/binning/manualBinning/binsGood/ANI_comparison
cd ~/BLiMMP/dataEdited/binning/manualBinning/binsGood/ANI_comparison
wget https://ani.jgi-psf.org/download_files/ANIcalculator_v1.tgz
tar -xzvf ANIcalculator_v1.tgz
git clone https://github.com/sstevens2/ani_compare_dag.git
mv ani_compare_dag BLI_bins_ANI
cd BLI_bins_ANI/
mkdir goodBins



cp ~/BLiMMP/dataEdited/binning/manualBinning/binsGood/ORFs/*fna goodBins/

echo 'goodBins' > groupslist.txt
# Change path of executable and transfer_input_files location in the group.sub file
# Stored here: /Users/benjaminpeterson/Documents/research/BLiMMP/code/binning/group.sub
# transfer_input_files lines
#executable = /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/binning/manualBinning/binsGood/ANI_comparison/BLI_bins_ANI/group.sh
#transfer_input_files = /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/binning/manualBinning/binsGood/ANI_comparison/ANIcalculator_v1/ANIcalculator,/home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/binning/manualBinning/binsGood/ANI_comparison/ANIcalculator_v1/nsimscan,$(spllist),$(totransfer)
cd ~/BLiMMP/dataEdited/binning/manualBinning/binsGood/ANI_comparison/BLI_bins_ANI
condor_submit_dag runAllANIcompare.dag
# Download output file (goodBins.all.ani.out.cleaned)
# to my computer:
#

cd ~/BLiMMP/dataEdited/binning/manualBinning
mkdir binsFinal
# Upload "final_bins.txt" list to binsFinal folder
mkdir binsFinal/DNA
mkdir binsFinal/ORFs
cat binsFinal/final_bins.txt | while read bin
do
  cp binsGood/DNA/$bin.fna binsFinal/DNA
  cp binsGood/ORFs/$bin* binsFinal/ORFs
done

# Aggregate ORFs
cd ~/BLiMMP/dataEdited/binning/manualBinning/binsFinal
$scripts/Fasta_to_Scaffolds2Bin.sh -e faa \
                                    -i ORFs \
                                    > binsFinal_G2B.tsv
cat ORFs/*faa > ORFs.faa

cd ~/BLiMMP/dataEdited/binning/manualBinning/binsFinal
$scripts/Fasta_to_Scaffolds2Bin.sh -e fna \
                                    -i DNA \
                                    > binsFinal_S2B.tsv
cat DNA/*fna > DNA.fna




####################################################
####################################################
# Metabolic analyses
####################################################
####################################################

##########################
# Custom set of metabolic HMMs
##########################
screen -S BLI_metabolic_HMMs
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate batch_HMMs
PYTHONPATH=''
PERL5LIB=''
binsFinal=~/BLiMMP/dataEdited/binning/manualBinning/binsFinal
scripts=~/BLiMMP/code
metabolic_HMMs=~/BLiMMP/references/metabolic_HMMs
cd $binsFinal
mkdir metabolism
chmod +x $scripts/batch_HMMs.py

python $scripts/batch_HMMs.py --orf_file $binsFinal/ORFs.faa \
                              --g2b $binsFinal/binsFinal_G2B.tsv \
                              --hmm_folder $metabolic_HMMs \
                              --hmm_csv $metabolic_HMMs\_bins.csv \
                              --output $binsFinal/metabolism/batch_HMMs
conda deactivate


##########################
# Search for MHCs
##########################
screen -S MHCs
mkdir $binsFinal/metabolism/MHCs
cd $binsFinal/metabolism/MHCs
scripts=~/BLiMMP/code
binsFinal=~/BLiMMP/dataEdited/binning/manualBinning/binsFinal
$scripts/Find_multiheme_protein.py $binsFinal/ORFs.faa 3
mv $binsFinal/ORFs_3_heme* .

echo -e "binID\tgeneID\themeCount" > heme_count_bins.tsv
tail -n +2 ORFs_3_heme_count.txt | awk -F '\t' '{ print $1 }' | while read geneID
do
  binID=`awk -F '\t' -v geneID="$geneID" '{ if ($1 == geneID) print $2 }' $binsFinal/binsFinal_G2B.tsv`
  hemeCount=`awk -F '\t' -v geneID="$geneID" '{ if ($1 == geneID) print $2 }' ORFs_3_heme_count.txt`
  echo -e $binID"\t"$geneID"\t"$hemeCount
  echo -e $binID"\t"$geneID"\t"$hemeCount >> heme_count_bins.tsv
done

##########################
# Search for MHCs with 10 heme-binding sites
##########################
mkdir $binsFinal/metabolism/MHCs_10
cd $binsFinal/metabolism/MHCs_10
$scripts/Find_multiheme_protein.py $binsFinal/ORFs.faa 10
mv $binsFinal/ORFs_10_heme* .

echo -e "binID\tgeneID\themeCount" > heme_count_bins_10.tsv
tail -n +2 ORFs_10_heme_count.txt | awk -F '\t' '{ print $1 }' | while read geneID
do
  binID=`awk -F '\t' -v geneID="$geneID" '{ if ($1 == geneID) print $2 }' $binsFinal/binsFinal_G2B.tsv`
  hemeCount=`awk -F '\t' -v geneID="$geneID" '{ if ($1 == geneID) print $2 }' ORFs_10_heme_count.txt`
  echo -e $binID"\t"$geneID"\t"$hemeCount
  echo -e $binID"\t"$geneID"\t"$hemeCount >> heme_count_bins.tsv
done

##########################
# Search for BBOMPs
##########################
# Pull out names of adjacent genes
binsFinal=~/BLiMMP/dataEdited/binning/manualBinning/binsFinal
cd $binsFinal/metabolism
mkdir PCC
mkdir PCC/list
rm -f $binsFinal/metabolism/PCC/adjacent_genes_all_list.txt

cat MHCs/ORFs_3_heme_list.txt | while read gene
do
  echo "Working on" $gene
  scaffold=$(echo $gene | rev | cut -d"_" -f2- | rev)
  assembly=$(echo $gene | rev | cut -d"_" -f3- | rev)
  ORFnumber=$(echo $gene | rev | cut -d"_" -f1 | rev)
  preceedingORFnumber=$(expr $ORFnumber - 1)
  followingORFnumber=$(expr $ORFnumber + 1)
  echo $scaffold"_"$preceedingORFnumber >> PCC/list/adjacent_genes_all_list.txt
  echo $scaffold"_"$followingORFnumber >> PCC/list/adjacent_genes_all_list.txt
done

# Find unique gene names
cd PCC
wc -l list/adjacent_genes_all_list.txt
sort list/adjacent_genes_all_list.txt | \
  uniq \
  > list/adjacent_genes_unique_list.txt

# Pull out adjacent genes
rm -f adjacent_genes.faa
cat list/adjacent_genes_unique_list.txt | while read geneID
do
  assembly=$(echo $geneID | rev | cut -d"_" -f3- | rev)
  echo "Looking for" $geneID "in" $assembly
  grep -A 1 -m 1 $geneID$ $binsFinal/ORFs.faa >> adjacent_genes.faa
done


#########################
# Search adjacent genes for BBOMPs
#########################
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=""
PERL5LIB=""
cd $binsFinal/metabolism/PCC

# Run the PCC HMM
pcc_omp_HMM=~/BLiMMP/references/metabolicProteins/EET/pcc_omp.HMM
hmmsearch --tblout pcc_omp_custom.out \
          -T 40 \
          $pcc_omp_HMM \
          adjacent_genes.faa \
          > pcc_omp_custom_output.txt


#########################
# Run METABOLIC
#########################
mkdir ~/BLiMMP/dataEdited/binning/manualBinning/binsFinal/metabolism/forMETABOLIC
cd ~/BLiMMP/dataEdited/binning/manualBinning/binsFinal
cat final_bins.txt | while read bin
do
  cp DNA/$bin.fna metabolism/forMETABOLIC/$bin.fasta
done
cd metabolism/forMETABOLIC
condor_submit ~/BLiMMP/code/METABOLIC_processing.sub



#########################
# Run Kofamscan on ORFs
#########################
cp -avr ~/references/kofamscan_DBs/kofam_scan-1.3.0 ~/BLiMMP/code/
# Update config.yml file
binsFinal=~/BLiMMP/dataEdited/binning/manualBinning/binsFinal
mkdir $binsFinal/metabolism/KOFAM_output

# Run kofamscan on each set of bins
screen -S kofamscan
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate kofamscan
PYTHONPATH=""
PERL5LIB=""
binsFinal=~/BLiMMP/dataEdited/binning/manualBinning/binsFinal
KOFAM_output=/home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/binning/manualBinning/binsFinal/metabolism/KOFAM_output
ORFs=/home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/binning/manualBinning/binsFinal/ORFs
cd /home/GLBRCORG/bpeterson26/BLiMMP/code/kofam_scan-1.3.0

rm -f $KOFAM_output/output_forDecoder.tsv
cat $binsFinal/final_bins.txt | while read binID
do
  if [ ! -d $KOFAM_output/$binID ]; then
    echo "Running KOFAMscan on" $binID
    ./exec_annotation -f detail-tsv \
                      -o $KOFAM_output/$binID.tsv \
                      $ORFs/$binID.faa

    echo "Cleaning up results from" $binID
    grep '*' $KOFAM_output/$binID.tsv > $KOFAM_output/$binID\_qualityHits.tsv
    binIDclean=`echo $binID | sed 's/_//g'`
    awk -F '\t' -v binIDclean="$binIDclean" '{ print binIDclean"_"$2"\t"$3 }' $KOFAM_output/$binID\_qualityHits.tsv >> $KOFAM_output/output_forDecoder.tsv
  fi
done
conda deactivate


#########################
# Run KEGG-decoder on output
#########################
screen -S kegg_decoder
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate kegg_decoder
PYTHONPATH=""
PERL5LIB=""
KOFAM_output=/home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/binning/manualBinning/binsFinal/metabolism/KOFAM_output
KEGG-decoder --input $KOFAM_output/output_forDecoder.tsv \
              --output $KOFAM_output/bin_metabolism_decoded_html.txt \
              --vizoption interactive


#########################
# Identify CAZymes in each bin
#########################
#http://bcb.unl.edu/dbCAN2/blastation.php?jobid=20210602130442


#########################
# Characterize MoORs in bins
#########################
# Search bins for MoORs
screen -S BLI_MoORs
cd ~/BLiMMP/dataEdited/binning/manualBinning/binsFinal
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''
scripts=~/BLiMMP/code/

mkdir metabolism/MoORs
hmmsearch --tblout metabolism/MoORs/MoOR.out \
          -T 50 \
          ~/BLiMMP/references/custom_hmms/MoOR.HMM \
          ORFs.faa

# Pull out the gene names, use this to find correct portion of G2B file.
cd metabolism/MoORs/
rm -f MoOR_G2B.tsv
rm -f putative_MoORs.faa
grep -v '#' MoOR.out | \
  awk '{ print $1 }' | while read gene
  do
    awk -F '\t' -v gene="$gene" '$1 == gene { print $0 }' ../../binsFinal_G2B.tsv >> MoOR_G2B.tsv
    grep -A 1 $gene\$ ../../ORFs.faa >> putative_MoORs.faa
    tail -n 1 MoOR_G2B.tsv
  done
sed -i 's/*//g' putative_MoORs.faa

# Align the putative MoORs to the references
hmmalign -o putative_MoORs.sto \
          ~/BLiMMP/references/custom_hmms/MoOR.HMM \
          putative_MoORs.faa
$scripts/convert_stockhold_to_fasta.py putative_MoORs.sto
# Check this in Geneious
# Removed sequences with fewer than 400 aligned residues
muscle -profile \
        -in1 putative_MoORs_trimmed.afa \
        -in2 ~/BLiMMP/references/custom_hmms/MoOR_reference.afa \
        -out putative_MoORs_ref_1.afa

python $scripts/cleanFASTA.py putative_MoORs_ref_1.afa
mv -f putative_MoORs_ref_1.afa_temp.fasta putative_MoORs_ref_1.afa
#sed 's/-//g' putative_MoORs_ref_1.afa > putative_MoORs_ref_1.faa

# Mask the alignment at 50% gaps
trimal -in putative_MoORs_ref_1.afa \
        -out putative_MoORs_ref_1_masked.afa \
        -gt 0.5
FastTree putative_MoORs_ref_1_masked.afa > putative_MoORs_1.tree


          ####################################################
          ####################################################
          # hgcA analysis
          ####################################################
          ####################################################

          screen -S BLI_hgcA_bins
          source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
          conda activate bioinformatics
          PYTHONPATH=''
          PERL5LIB=''
          binsGood=~/BLiMMP/dataEdited/binning/manualBinning/binsGood
          hgcaAssembly=~/BLiMMP/dataEdited/binning/manualBinning/binsGood

          scripts=~/BLiMMP/code/generalUse/
          cd $binsGood
          mkdir hgcA

          hmmsearch --tblout hgcA/hgcA.out \
                    --cpu 4 \
                    --cut_tc \
                    ~/references/hgcA/hgcA.hmm \
                    ORFs.faa \
                    > hgcA/hgcA_report.txt
          cd hgcA
          grep -v "#" hgcA.out | awk '{ print $1 }' > hgcA_list.txt

          # Get hgcA to bin file
          echo -e 'hgcA_ID\tbinID' > hgcA_to_bin.tsv
          cat hgcA_list.txt | while read hgcA
          do
            awk -F '\t' -v hgcA="$hgcA" '$1 == hgcA { print $0 }' $binsGood/binsFinal_G2B.tsv >> hgcA_to_bin.tsv
          done



####################################################
####################################################
# Pull out coverage information from bam files
####################################################
####################################################

binsFinal=~/BLiMMP/dataEdited/binning/manualBinning/binsFinal
cd $binsFinal
mkdir depth
awk -F '\t' '{ print $1 }' /home/GLBRCORG/bpeterson26/BLiMMP/metadata/metagenome_list.txt > metagenomes_to_use.txt
cd ~/BLiMMP/reports
mkdir aggregate_depth_bins
cd aggregate_depth_bins
mkdir outs logs errs
cd /home/GLBRCORG/bpeterson26/BLiMMP/code
chmod +x aggregate_depth_bins.sh
chmod +x calculate_depth_of_bins.py
condor_submit aggregate_depth_bins.sub
