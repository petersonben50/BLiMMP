screen -S BLI_dsr_bins

finalBins=~/BLiMMP/dataEdited/binning/autoBinning/finalBins
cd $finalBins
mkdir metabolism
mkdir metabolism/dsr

dsrA_list=~/BLiMMP/dataEdited/metabolic_analyses/identification/dsrA_from_2020_all_list.txt
rm -f metabolism/dsr/dsrA_bin_S2B.tsv
cat $dsrA_list | while read dsrA
do
  scaffoldID=`echo $dsrA | cut -d"_" -f1-3`
  echo $dsrA "from" $scaffoldID
  grep -h $scaffoldID binning/*_S2B.tsv >> metabolism/dsr/dsrA_bin_S2B.tsv
done


dsrD_list=~/BLiMMP/dataEdited/metabolic_analyses/identification/dsrD_from_2020_all_list.txt
rm -f metabolism/dsr/dsrD_bin_S2B.tsv
cat $dsrD_list | while read dsrD
do
  scaffoldID=`echo $dsrD | cut -d"_" -f1-3`
  echo $dsrD "from" $scaffoldID
  grep -h $scaffoldID binning/*_S2B.tsv >> metabolism/dsr/dsrD_bin_S2B.tsv
done


cd metabolism/dsr/
awk -F '\t' '{ print $2 }' dsrA_bin_S2B.tsv | while read binID
do
  grep $binID $finalBins/taxonomy/taxonomy_summary.txt
done



############################
# Mapping to 2017 reference bins
############################

scripts=/home/GLBRCORG/bpeterson26/BLiMMP/code/

cd /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/binning
mkdir mapping_to_2017_MGs
mkdir mapping_to_2017_MGs/DNA
cp ~/5M/dataEdited/binAnalysis/bins_processed/*.fna mapping_to_2017_MGs/DNA
rm -rf mapping_to_2017_MGs/DNA/KIR_0008_original.fna
cd mapping_to_2017_MGs
cat DNA/*fna > bin_scaffolds.fna
cd DNA
$scripts/Fasta_to_Scaffolds2Bin.sh -e fna > ../bin_scaffolds_S2B.tsv

cd /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/binning/mapping_to_2017_MGs
mkdir mapping/
mkdir mapping/indices

# Build indicies
if [ ! -e mapping/indices/bin_scaffolds.1.bt2 ]; then
  /opt/bifxapps/bowtie2-2.2.2/bowtie2-build bin_scaffolds.fna \
                                            mapping/indices/bin_scaffolds
else
  echo "Already made index for" $assembly
fi


awk -F '\t' '{ print $1 }' ~/BLiMMP/metadata/metagenome_list.txt > metagenome_list_to_use.txt


cd /home/GLBRCORG/bpeterson26/BLiMMP/reports/
mkdir binning2017_mapping
cd binning2017_mapping
rm -fr outs errs logs
mkdir outs errs logs

cd /home/GLBRCORG/bpeterson26/BLiMMP/code/
chmod +x binning2017_mapping.sh
condor_submit binning2017_mapping.sub




############################
# Pull out coverage information from bam files
############################

finalBins=~/BLiMMP/dataEdited/binning/autoBinning/finalBins
cd $finalBins
mkdir depth
awk -F '\t' '{ print $1 }' /home/GLBRCORG/bpeterson26/BLiMMP/metadata/metagenome_list.txt > metagenomes_to_use.txt
cd ~/BLiMMP/reports
mkdir aggregate_depth_bins2017
cd aggregate_depth_bins2017
mkdir outs logs errs
cd /home/GLBRCORG/bpeterson26/BLiMMP/code
chmod +x aggregate_depth_bins2017.sh
chmod +x calculate_depth_of_bins.py
condor_submit aggregate_depth_bins2017.sub
