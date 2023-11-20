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




