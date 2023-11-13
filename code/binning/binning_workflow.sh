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
rm -rf binning_mapping
mkdir binning_mapping
cd binning_mapping
mkdir outs errs logs

cd /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/binning
mkdir mapping
mkdir mapping/indices
cp ~/BLiMMP/metadata/assembly_list.txt assembly_list_to_use.txt

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
assembly=BLI21_coassembly
bin=Bin_66
anvi-refine -p anvioDBs_modified/$assembly.merged/PROFILE.db \
            -c anvioDBs_modified/$assembly.db \
            -C CONCOCT \
            -b $bin \
            -A anvioDBs_modified/binning_collections/$assembly\_collections.tsv \
            --taxonomic-level "t_phylum"


################################################
################################################
# Aggregate all bins and identify final bins
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
# Pull out DNA files from all bins
##########################
# First need to set up new directory
cd ~/BLiMMP/dataEdited/binning
mkdir binsRaw
mkdir binsRaw/DNA
binsRaw=~/BLiMMP/dataEdited/binning/binsRaw
echo -e "assemblyID\tnewBinName\toldBinName" > $binsRaw/renaming_file.tsv

# Copy over bins from manual binning
awk -F '\t' '{ print $1 }' $assembly_list | while read assembly
do
  i=1
  binSummary=~/BLiMMP/dataEdited/binning/manualBinning/anvioDBs_modified/$assembly.curated.summary
  if [ ! -e $binsRaw/DNA/$assembly* ]; then
    if [ -e $binSummary/bin_by_bin ]; then
      cd $binSummary/bin_by_bin
      ls | sed 's/\///' | while read bin
      do
        if [ -e $bin/$bin-contigs.fa ]; then
          newFile="$(printf "$assembly\_anvio_bin_%04d.fna" "$i" | sed 's/\\//g')"
          echo "Copying" $bin "to" $newFile "in binsRaw folder"
          cp $bin/$bin-contigs.fa $binsRaw/DNA/$newFile
          echo -e $assembly"\t"$newFile"\t"$bin >> $binsRaw/renaming_file.tsv
          i=$((i+1))
        fi
      done
    fi
  else
    echo "Hey, there are some bins from" $assembly "already in here"
  fi
done

# Copy over bins from auto binning
dasTool=/home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/binning/autoBinning/dasTool
awk -F '\t' '{ print $1 }' $assembly_list | while read assembly
do
  if [ -e $dasTool/$assembly\_output/$assembly\_bins_DASTool_bins ]; then
    cd $dasTool/$assembly\_output/$assembly\_bins_DASTool_bins
    i=1
    ls *.fa | while read file
    do
      newFile="$(printf "$assembly\_dasTool_bin_%04d.fna" "$i" | sed 's/\\//g')"
      echo "Moving" "$file" "to" $newFile
      echo -e $assembly"\t"$newFile"\t"$file >> $binsRaw/renaming_file.tsv
      cp $file $binsRaw/DNA/$newFile
      i=$((i+1))
    done
  fi
done


##########################
# Completeness/redundancy estimates from CheckM
##########################

screen -S BLI_checkM
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=""
binsRaw=~/BLiMMP/dataEdited/binning/binsRaw

cd $binsRaw
if [ -d checkM ]; then
  echo "Removing old checkM folder"
  rm -rf checkM
fi
mkdir checkM
checkm lineage_wf \
      -x .fna \
      -t 24 \
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

##########################
# Keep bins above 50% complete and below 10% redundancy
##########################

cd ~/BLiMMP/dataEdited/binning
mkdir binsGood
mkdir binsGood/DNA
awk -F ',' '{ if (($2 > 50) && ($3 < 10)) print $0 }' \
  binsRaw/checkM/checkM_stats.csv \
  > binsGood/checkM_stats.txt
awk -F ',' '{ print $1 }' binsGood/checkM_stats.txt | while read binsGood
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

