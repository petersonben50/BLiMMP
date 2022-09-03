

############################################
############################################
# Metagenome assembly
############################################
############################################

######################
# Assemblies by metaSPAdes
######################
mkdir ~/BLiMMP/dataEdited/assemblies
mkdir ~/BLiMMP/dataEdited/assemblies/assembly_files

screen -S BLiMMP_metagenome_coassembly
cd ~/BLiMMP/dataEdited/assemblies
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=""
PERL5LIB=""

code=~/BLiMMP/code
assembly_grouping=~/BLiMMP/metadata/assembly_key.csv
read_storage=~/BLiMMP/dataEdited/metagenomes
output=~/BLiMMP/dataEdited/assemblies/assembly_files

#chmod +x $code/assembly_by_group.py

awk -F '\t' '{ print $1 }' ~/BLiMMP/metadata/assembly_list.txt | while read assembly
do
  if [ ! -d $output/$assembly ]; then
    mkdir $output/$assembly
  fi
  if [ ! -e $output/$assembly/scaffolds.fasta ]; then
    echo "Assembling" $assembly
    python $code/assembly_by_group.py $assembly \
                                      $assembly_grouping \
                                      $read_storage \
                                      $output/$assembly
    # To continue a paused run:
    # assembly=XXXXXX
    # metaspades.py --continue -o assembly
  else
    echo $assembly "already assembled"
  fi
done


######################
# Clean up assemblies
######################
screen -S HCC_clean_metagenome_assembly
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate anvio6.2
PYTHONPATH=''
cd ~/BLiMMP/dataEdited/assemblies
mkdir scaffolds
mkdir scaffolds/renaming_reports

awk -F '\t' '{ print $1 }' ~/BLiMMP/metadata/assembly_list.txt | while read assembly
do
  if [ -e assembly_files/$assembly/scaffolds.fasta ]; then
    echo $assembly "is done assembling. Let's clean it."
    if [ -e scaffolds/$assembly\_assembly.fna ]; then
      echo "Dude, you already cleaned the assembly for" $assembly". Relax."
    else
      echo "Cleaning the assembly for" $assembly
      anvi-script-reformat-fasta assembly_files/$assembly/scaffolds.fasta \
                                  -o scaffolds/$assembly\_assembly.fna \
                                  -l 1000 \
                                  --simplify-names \
                                  --prefix $assembly \
                                  --report-file scaffolds/renaming_reports/$assembly\_report_file.txt
    fi
  else
    echo "Yo, you gotta go back and assemble" $assembly
  fi
done
conda deactivate


############################################
############################################
# Predict ORFs
############################################
############################################

# Get set up
#cd /home/GLBRCORG/bpeterson26/BLiMMP/
#cp /home/GLBRCORG/bpeterson26/HellsCanyon/code/generalUse/cleanFASTA.py code/cleanFASTA.py
#mkdir reports
#cd reports
#mkdir outs errs logs


mkdir /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/assemblies/ORFs
chmod +x /home/GLBRCORG/bpeterson26/BLiMMP/code/ORF_prediction_assemblies.sh
condor_submit /home/GLBRCORG/bpeterson26/BLiMMP/code/ORF_prediction_assemblies.sub

# Clean up ORF files
cd ~/BLiMMP/dataEdited/assemblies/ORFs
code=~/BLiMMP/code
awk -F '\t' '{ print $1 }' ~/BLiMMP/metadata/assembly_list.txt | while read assembly
do
  echo "Working on cleaning ORFs from" $assembly
  python $code/cleanFASTA.py $assembly.faa
  mv -f $assembly.faa_temp.fasta $assembly.faa
  python $code/cleanFASTA.py $assembly.fna
  mv -f $assembly.fna_temp.fasta $assembly.fna
done


############################################
############################################
# Get assembly stats
############################################
############################################

screen -S BLI_assembly_stats
cd ~/BLiMMP/dataEdited/assemblies
mkdir reports
code=~/BLiMMP/code
cp ~/HellsCanyon/code/assembly_processing/abyss-fac.pl $code/

# Use abyss-fac to check stats
awk -F '\t' '{ print $1 }' ~/BLiMMP/metadata/assembly_list.txt | while read assembly
do
  if [ -e scaffolds/$assembly\_assembly.fna ]; then
    echo $assembly "has been cleaned. Let's get some stats on it."
    if [ ! -e reports/$assembly\_report.txt ]; then
      perl $code/abyss-fac.pl scaffolds/$assembly\_assembly.fna \
          > reports/$assembly\_report.txt
    else
      echo "Already checked" $assembly
    fi
  else
    echo $assembly "has not been cleaned. Do that first."
  fi
done

# Aggregate stats
cd ~/BLiMMP/dataEdited/assemblies/reports
echo -e "n\tn:200\tL50\tmin\tN80\tN50\tN20\tmax\tsum\tassemblyID" > all_assemblies_stats.txt
for file in *report.txt
do
  tail -n +2 $file >> all_assemblies_stats.txt
done
rm BLI*_report.txt

# Clean up the report file
awk -F '\t' '{ print $1 }' ~/BLiMMP/metadata/assembly_list.txt | while read assembly
do
  sed "s/scaffolds\/$assembly\_assembly.fna/$assembly/" all_assemblies_stats.txt \
    > all_assemblies_stats.txt_edited
  mv -f all_assemblies_stats.txt_edited all_assemblies_stats.txt
done


# Count ORFs
cd ~/BLiMMP/dataEdited/assemblies
echo -e 'assemblyID\tORF_count' > reports/ORF_counts.tsv
awk -F '\t' '{ print $1 }' ~/BLiMMP/metadata/assembly_list.txt | while read assembly
do
  if [ -e ORFs/$assembly.faa ]; then
    echo "Count ORFs in" $assembly
    ORF_count=$(grep '>' ORFs/$assembly.faa | wc -l)
    echo -e "$assembly\t$ORF_count" >> reports/ORF_counts.tsv
  fi
done
