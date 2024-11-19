####################################################
####################################################
# Screen bins for hgcA
####################################################
####################################################
# Get set up
cd /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/
mkdir bin_based_analyses
mkdir bin_based_analyses/bin_screen_for_hgcA
mkdir bin_based_analyses/bin_screen_for_hgcA/DNA
cp -avr binning/binsGood/DNA/*fna bin_based_analyses/bin_screen_for_hgcA/DNA
cp binning/binsGood/checkM_stats.txt bin_based_analyses/bin_screen_for_hgcA

# Get ORFs for bins
screen -S BLI_binsORFS
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=""
scripts=/home/GLBRCORG/bpeterson26/BLiMMP/code
bin_screen=/home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/bin_based_analyses/bin_screen_for_hgcA
mkdir $bin_screen/ORFs

cd $bin_screen/DNA
ls *.fna | sed 's/.fna//' | while read bin
do
  if [ ! -e ORFs/$bin.gff ]; then
    echo "Predicting proteins for" $bin
    prodigal -i $bin.fna \
              -o $bin_screen/ORFs/$bin.gff \
              -f gff \
              -a $bin_screen/ORFs/$bin.faa \
              -d $bin_screen/ORFs/$bin.fna \
              -p single
  else
    echo $bin "already run."
  fi
done

# Clean up the gene sequence files
cd $bin_screen/ORFs
ls *.fna | sed 's/.fna//' | while read bin
do
  echo "Cleaning" $bin
  python $scripts/cleanFASTA.py $bin.fna
  mv -f $bin.fna_temp.fasta $bin.fna
  python $scripts/cleanFASTA.py $bin.faa
  mv -f $bin.faa_temp.fasta $bin.faa
done

# Combine gene sequence files and generate S2B/G2B files
cd ~/BLiMMP/dataEdited/bin_based_analyses/bin_screen_for_hgcA/
$scripts/Fasta_to_Scaffolds2Bin.sh -e fna \
                                    -i DNA \
                                    > bins_DNA_S2B.tsv
cat DNA/*.fna > DNA.fna

$scripts/Fasta_to_Scaffolds2Bin.sh -e faa \
                                    -i ORFs \
                                    > bins_ORF_G2B.tsv
cat ORFs/*.faa > ORFs.faa
exit


# Check for hgcA gene
screen -S BLI_hgcA_bins
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''
bin_screen=/home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/bin_based_analyses/bin_screen_for_hgcA
scripts=~/BLiMMP/code/

cd $bin_screen
mkdir hgcA

# Search for hgcA sequences
hmmsearch --tblout hgcA/hgcA.out \
        --cpu 4 \
        --cut_tc \
        ~/references/hgcA/hgcA.hmm \
        ORFs.faa \
        > hgcA/hgcA_report.txt
cd hgcA
grep -v "#" hgcA.out | awk '{ print $1 }' | sort | uniq > hgcA_list.txt

# Pull out hgcA sequences for alignment and manual curation.
cat hgcA_list.txt | while read hgcA
do
  echo "Working on" $hgcA
  grep -A 1 -m 1 $hgcA$ ../ORFs.faa >> hgcA.faa
done

muscle -align hgcA.faa \
        -output hgcA.afa

# Manually inspected these in NCBI MSA Viewer. Remove three:
# Two missing cap helix domain:
#BLI21_coassembly_000000261085_4
#BLI21_assembly106_000000095479_4
# One is just too truncated:
#BLI21_coassembly_000000028244_10
echo BLI21_coassembly_000000261085_4 >> hgcA_to_remove.txt
echo BLI21_assembly106_000000095479_4 >> hgcA_to_remove.txt
echo BLI21_coassembly_000000028244_10 >> hgcA_to_remove.txt

grep -v -f hgcA_to_remove.txt hgcA_list.txt > hgcA_clean_list.txt





####################################################
####################################################
# Move hgcA+ bins over
####################################################
####################################################
# Set up directories
cd /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/bin_based_analyses
mkdir hgcA_bins
mkdir hgcA_bins/DNA
mkdir hgcA_bins/ORFs
mkdir hgcA_bins/final_bin_data
cp bin_screen_for_hgcA/hgcA/hgcA_clean_list.txt hgcA_bins
grep -wFf hgcA_bins/hgcA_clean_list.txt bin_screen_for_hgcA/bins_ORF_G2B.tsv > hgcA_bins/final_bin_data/hgcA_to_bin.tsv

# Move bins over
awk -F '\t' '{ print $2 }' hgcA_bins/final_bin_data/hgcA_to_bin.tsv | while read bin
do
  if [ ! -e hgcA_bins/DNA/$bin.fna ]; then
    echo "Working on" $bin
    cp bin_screen_for_hgcA/DNA/$bin.fna hgcA_bins/DNA
    cp bin_screen_for_hgcA/ORFs/$bin.* hgcA_bins/ORFs
    grep $bin bin_screen_for_hgcA/checkM_stats.txt >> hgcA_bins/final_bin_data/checkM_stats.tsv
  else
    echo $bin "already moved."
  fi
done
cat hgcA_bins/final_bin_data/checkM_stats.tsv | sort | uniq > hgcA_bins/final_bin_data/checkM_stats_clean.tsv
rm -f hgcA_bins/final_bin_data/checkM_stats.tsv

# Generate ORFs- and scaffolds-to-bin files
cd /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/bin_based_analyses/hgcA_bins
scripts=/home/GLBRCORG/bpeterson26/BLiMMP/code
$scripts/Fasta_to_Scaffolds2Bin.sh -e fna \
                                    -i DNA \
                                    > final_bin_data/hgcA_bins_S2B.tsv
$scripts/Fasta_to_Scaffolds2Bin.sh -e faa \
                                    -i ORFs \
                                    > final_bin_data/hgcA_bins_G2B.tsv
cat ORFs/*faa > ORFs.faa

####################################################
####################################################
# Check out taxonomy of bins with GTDB
####################################################
####################################################
screen -S BLI_GTDB
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
PYTHONPATH=""
conda activate gtdbtk
cd /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/bin_based_analyses/hgcA_bins
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
        > ../final_bin_data/taxonomy_summary.txt
conda deactivate


####################################################
####################################################
# Run ANI comparisons
####################################################
####################################################
# The following workflow is the code to run
# Sarah Stevens's ANI calculator on a
# folder full of bins.
# Details: https://github.com/sstevens2/ani_compare_dag

mkdir /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/bin_based_analyses/hgcA_bins/ANI_comparison
cd /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/bin_based_analyses/hgcA_bins/ANI_comparison
wget https://ani.jgi-psf.org/download_files/ANIcalculator_v1.tgz
tar -xzvf ANIcalculator_v1.tgz
git clone https://github.com/sstevens2/ani_compare_dag.git
mv ani_compare_dag BLI_bins_ANI
cd BLI_bins_ANI/
mkdir goodBins
cp /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/bin_based_analyses/hgcA_bins/ORFs/*.fna goodBins/
echo 'goodBins' > groupslist.txt
# Change path of executable and transfer_input_files location in the group.sub file
sed -i "s/executable = group.sh/executable = \/home\/GLBRCORG\/bpeterson26\/BLiMMP\/dataEdited\/bin_based_analyses\/hgcA_bins\/ANI_comparison\/BLI_bins_ANI\/group.sh/" group.sub
sed -i "s/transfer_input_files = \/home\/sstevens2\/ANIcalculator_v1\/ANIcalculator,\/home\/sstevens2\/ANIcalculator_v1\/nsimscan,\$(spllist),\$(totransfer)/transfer_input_files = \/home\/GLBRCORG\/bpeterson26\/BLiMMP\/dataEdited\/bin_based_analyses\/hgcA_bins\/ANI_comparison\/ANIcalculator_v1\/ANIcalculator,\/home\/GLBRCORG\/bpeterson26\/BLiMMP\/dataEdited\/bin_based_analyses\/hgcA_bins\/ANI_comparison\/ANIcalculator_v1\/nsimscan,\$(spllist),\$(totransfer)/" group.sub
condor_submit_dag runAllANIcompare.dag
# Save output 
mv goodBins.all.ani.out.cleaned /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/bin_based_analyses/hgcA_bins/final_bin_data/ANI_comparison.txt


####################################################
####################################################
# Calculate metagenomic coverage of bins
####################################################
####################################################

# Need a bin to assembly tsv file:
cd /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/bin_based_analyses/hgcA_bins/final_bin_data
awk -F '\t' '{ print $2 }' hgcA_to_bin.tsv | while read binID
do
  assemblyID=$(echo $binID | cut -d"_" -f1,2)
  echo -e $binID"\t"$assemblyID
done | sort | uniq > bin_to_assembly.tsv

HomeBio=/home/GLBRCORG/bpeterson26/BLiMMP/code/HomeBio
cd
python $HomeBio/homebio_workflows/calculate_coverage_of_bin_set.py --bam_folder testing_bin_cov/mapping \
                                                                    --s2b_file testing_bin_cov/bins/hgcA_bins_S2B.tsv \
                                                                    --b2a_file testing_bin_cov/bins/bin_to_assembly.tsv \
                                                                    --output_file testing_bin_cov/test_output.tsv \
                                                                    --exclude_bases 150 \
                                                                    --cores 24




####################################################
####################################################
# Identify housekeeping gene gyrB in bins
####################################################
####################################################
screen -S BLI_gyrB_bins
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''
hgcA_bins=/home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/bin_based_analyses/hgcA_bins
scripts=~/BLiMMP/code/

cd $hgcA_bins
mkdir gyrB

# Search for hgcA sequences
hmmsearch --tblout gyrB/gyrB.out \
        --cpu 4 \
        --cut_tc \
        ~/references/TIGR01059.HMM \
        ORFs.faa \
        > gyrB/gyrB_report.txt
cd gyrB
grep -v "#" gyrB.out | awk '{ print $1 }' | sort | uniq > gyrB_list.txt

# Pull out bin info:
grep -w -f gyrB_list.txt ../final_bin_data/hgcA_bins_G2B.tsv > ../final_bin_data/gyrB_to_bin.tsv




##################################################
##################################################
# Pseudoalign RNA reads to bins
##################################################
##################################################
mkdir ~/BLiMMP/dataEdited/metatranscriptomes/alignment_bins
mkdir ~/BLiMMP/dataEdited/metatranscriptomes/alignment_bins/ORFs
mkdir ~/BLiMMP/dataEdited/metatranscriptomes/alignment_bins/index
mkdir ~/BLiMMP/dataEdited/metatranscriptomes/alignment_bins/reports
mkdir /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/metatranscriptomes/alignment_bins/output

screen -S kallisto_indexing
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
PYTHONPATH=""
referenceDirectory=/home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/bin_based_analyses/hgcA_bins/ORFs/
indexFolder=/home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/metatranscriptomes/alignment_bins/index
tempRefFolder=/home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/metatranscriptomes/alignment_bins/ORFs
outputFolder=/home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/metatranscriptomes/alignment_bins/output
readsLocation=/home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/metatranscriptomes/mRNA_reads
mtList=/home/GLBRCORG/bpeterson26/BLiMMP/metadata/metatranscriptome_list.txt

conda activate kallisto

# Generated needed index files
cd $indexFolder
cat ~/BLiMMP/metadata/assembly_list.txt | while read assemblyID
do
  if test -n "$(find $referenceDirectory -maxdepth 1 -name $assemblyID\_anvio*.fna -print -quit)"; then
    cat $referenceDirectory/$assemblyID\_anvio*.fna > $tempRefFolder/$assemblyID\_anvio_ORFs.fna
    sed -i 's/\*//' $tempRefFolder/$assemblyID\_anvio_ORFs.fna
    echo "Indexing" $assemblyID "for anvio"
    kallisto index -i $assemblyID\_anvio.idx $tempRefFolder/$assemblyID\_anvio_ORFs.fna
  else
    echo "No anvio ORFs for" $assemblyID
  fi
  if test -n "$(find $referenceDirectory -maxdepth 1 -name $assemblyID\_dasTool*.fna -print -quit)"; then
    cat $referenceDirectory/$assemblyID\_dasTool*.fna > $tempRefFolder/$assemblyID\_dasTool_ORFs.fna
    sed -i 's/\*//' $tempRefFolder/$assemblyID\_anvio_ORFs.fna
    echo "Indexing" $assemblyID "for Das Tool"
    kallisto index -i $assemblyID\_dasTool.idx $tempRefFolder/$assemblyID\_dasTool_ORFs.fna
  else
    echo "No Das Tool ORFs for" $assemblyID
  fi
done
conda deactivate

# Run pseudomapping
# Submit kallisto analyses
mkdir /home/GLBRCORG/bpeterson26/BLiMMP/reports/metatranscriptomes_bins
cd /home/GLBRCORG/bpeterson26/BLiMMP/reports/metatranscriptomes_bins
mkdir outs errs logs
chmod +x /home/GLBRCORG/bpeterson26/BLiMMP/code/kallisto_pseudoalignment_bins.sh
condor_submit /home/GLBRCORG/bpeterson26/BLiMMP/code/kallisto_pseudoalignment_bins.sub

# Aggregate data
screen -S MT_aggregate
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
PYTHONPATH=""
conda activate bioinformatics
python kallisto_aggregate_data.py


##################################################
##################################################
# Extract gene neighborhood information
##################################################
##################################################









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
scripts=~/BLiMMP/code
cd /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/bin_based_analyses/hgcA_bins/
mkdir -p metabolism
cd metabolism

# Split up the anvio and DasTool ORFs to analyze separately
ORFs=~/BLiMMP/dataEdited/bin_based_analyses/hgcA_bins/ORFs
cd /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/bin_based_analyses/hgcA_bins/metabolism
mkdir ORFs_for_metabolism
mkdir ORFs_for_metabolism/anvio
mkdir ORFs_for_metabolism/DasTool
cp $ORFs/*anvio*.faa ORFs_for_metabolism/anvio
cp $ORFs/*dasTool*.faa ORFs_for_metabolism/DasTool
# Generate S2B files
$scripts/Fasta_to_Scaffolds2Bin.sh -e faa \
                                    -i ORFs_for_metabolism/anvio \
                                    > ORFs_for_metabolism/anvio_bins_G2B.tsv
$scripts/Fasta_to_Scaffolds2Bin.sh -e faa \
                                    -i ORFs_for_metabolism/DasTool \
                                    > ORFs_for_metabolism/DasTool_bins_G2B.tsv
# Concatenate the ORFs
cat ORFs_for_metabolism/anvio/*.faa > ORFs_for_metabolism/anvio.faa
cat ORFs_for_metabolism/DasTool/*.faa > ORFs_for_metabolism/DasTool.faa

scripts=~/BLiMMP/code
metabolic_HMMs=/home/GLBRCORG/bpeterson26/BLiMMP/references/metabolic_HMMs
chmod +x $scripts/batch_HMMs.py
rm -rf batch_HMMs
mkdir batch_HMMs_processing
#bin_g2b=/home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/bin_based_analyses/hgcA_bins/final_bin_data/hgcA_bins_G2B.tsv
python $scripts/batch_HMMs.py --orf_file ORFs_for_metabolism/anvio.faa \
                              --g2b ORFs_for_metabolism/anvio_bins_G2B.tsv \
                              --hmm_folder $metabolic_HMMs \
                              --hmm_csv $metabolic_HMMs\_bins.csv \
                              --output batch_HMMs_anvio

python $scripts/batch_HMMs.py --orf_file ORFs_for_metabolism/DasTool.faa \
                              --g2b ORFs_for_metabolism/DasTool_bins_G2B.tsv \
                              --hmm_folder $metabolic_HMMs \
                              --hmm_csv $metabolic_HMMs\_bins.csv \
                              --output batch_HMMs_DasTool

cat batch_HMMs_anvio/bin_counts/all_bin_hits.tsv > ../final_bin_data/batch_HMMs_bin_hits.tsv
tail -n +2 batch_HMMs_DasTool/bin_counts/all_bin_hits.tsv >> ../final_bin_data/batch_HMMs_bin_hits.tsv

conda deactivate




#########################
# Run Kofamscan on ORFs
#########################
#cp -avr ~/references/kofamscan_DBs/kofam_scan-1.3.0 ~/BLiMMP/code/
# https://github.com/takaram/kofam_scan
# Update ~/BLiMMP/code/kofam_scan-1.3.0/config.yml file to the following:

'
Path to your KO-HMM database
# A database can be a .hmm file, a .hal file or a directory in which
# .hmm files are. Omit the extension if it is .hal or .hmm file
profile: /home/GLBRCORG/bpeterson26/references/kofamscan_DBs/profiles

# Path to the KO list file
ko_list: /home/GLBRCORG/bpeterson26/references/kofamscan_DBs/ko_list

# Path to an executable file of hmmsearch
# You do not have to set this if it is in your $PATH
# hmmsearch: /usr/local/bin/hmmsearch

# Path to an executable file of GNU parallel
# You do not have to set this if it is in your $PATH
# parallel: /usr/local/bin/parallel

# Number of hmmsearch processes to be run parallelly
cpu: 8
'

# Run kofamscan on each set of bins
screen -S kofamscan
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate kofamscan
PYTHONPATH=""
PERL5LIB=""
b2a_file=/home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/bin_based_analyses/hgcA_bins/final_bin_data/bin_to_assembly.tsv
KOFAM_output=/home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/bin_based_analyses/hgcA_bins/metabolism/KOFAM_output
ORFs=/home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/bin_based_analyses/hgcA_bins/ORFs/
final_bin_data=/home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/bin_based_analyses/hgcA_bins/final_bin_data
cd /home/GLBRCORG/bpeterson26/BLiMMP/code/kofam_scan-1.3.0
mkdir $KOFAM_output

awk -F '\t' '{ print $1 }' $b2a_file | while read binID
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


cd $KOFAM_output
ls *qualityHits.tsv | while read file
do
  binID=`echo $file | sed 's/_qualityHits.tsv//'`
  echo "Working on" $binID
  awk -F '\t' -v binID="$binID" '{ print binID"\t"$0 }' $file >> $final_bin_data/kofam_data.tsv
done



##########################
# Search for MHCs
##########################
screen -S MHCs
scripts=~/BLiMMP/code
MHC_output=/home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/bin_based_analyses/hgcA_bins/metabolism/MHCs
mkdir $MHC_output
cd $MHC_output
ORFs_for_metabolism=/home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/bin_based_analyses/hgcA_bins/metabolism/ORFs_for_metabolism

# Run the MHC-finding script
$scripts/Find_multiheme_protein.py $ORFs_for_metabolism/DasTool.faa 3
$scripts/Find_multiheme_protein.py $ORFs_for_metabolism/anvio.faa 3

# Move the results to the MHCs folder
mv $ORFs_for_metabolism/*_heme_* .

# Concatenate the data into a single file
echo -e "binID\tgeneID\themeCount" > heme3_count_bins.tsv
tail -n +2 DasTool_3_heme_count.txt | awk -F '\t' '{ print $1 }' | while read geneID
do
  binID=`awk -F '\t' -v geneID="$geneID" '{ if ($1 == geneID) print $2 }' $ORFs_for_metabolism/DasTool_bins_G2B.tsv`
  hemeCount=`awk -F '\t' -v geneID="$geneID" '{ if ($1 == geneID) print $2 }' DasTool_3_heme_count.txt`
  echo -e $binID"\t"$geneID"\t"$hemeCount
  echo -e $binID"\t"$geneID"\t"$hemeCount >> heme3_count_bins.tsv
done
tail -n +2 anvio_3_heme_count.txt | awk -F '\t' '{ print $1 }' | while read geneID
do
  binID=`awk -F '\t' -v geneID="$geneID" '{ if ($1 == geneID) print $2 }' $ORFs_for_metabolism/anvio_bins_G2B.tsv`
  hemeCount=`awk -F '\t' -v geneID="$geneID" '{ if ($1 == geneID) print $2 }' anvio_3_heme_count.txt`
  echo -e $binID"\t"$geneID"\t"$hemeCount
  echo -e $binID"\t"$geneID"\t"$hemeCount >> heme3_count_bins.tsv
done
mv heme3_count_bins.tsv ../../final_bin_data/heme3_count_bins.tsv
