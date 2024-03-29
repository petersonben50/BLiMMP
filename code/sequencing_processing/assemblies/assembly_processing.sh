#!/bin/sh

#########################
# code/assemblies/assembly_processing.sh
# Benjamin D. Peterson
#########################


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
mkdir ~/BLiMMP/reports/assemblies
mkdir ~/BLiMMP/dataEdited/assemblies/scaffolds
mkdir ~/BLiMMP/dataEdited/assemblies/scaffolds/renaming_reports
tail -n +2 ~/BLiMMP/metadata/assembly_key.csv | \
  awk -F ',' '{ print $1 }' | \
  sort | \
  uniq \
  > ~/BLiMMP/metadata/assembly_list.txt
chmod +x /home/GLBRCORG/bpeterson26/BLiMMP/code/assembly_by_group.py
chmod +x /home/GLBRCORG/bpeterson26/BLiMMP/code/assembly_by_group_execute.sh
condor_submit /home/GLBRCORG/bpeterson26/BLiMMP/code/assembly_by_group_submit.sub


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


############################################
############################################
# Get assembly stats
############################################
############################################
screen -S BLI_assembly_stats
cd ~/BLiMMP/dataEdited/assemblies
code=~/BLiMMP/code
$cp ~/HellsCanyon/code/assembly_processing/abyss-fac.pl $code/

# Use abyss-fac to check stats
cat ~/BLiMMP/metadata/assembly_list.txt | while read assembly
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
rm -f BLI*_report.txt
# Clean up the report file
cat ~/BLiMMP/metadata/assembly_list.txt | while read assembly
do
  sed -i "s/scaffolds\/$assembly\_assembly.fna/$assembly/" all_assemblies_stats.txt
done


############################################
############################################
# Count ORFs
############################################
############################################
cd ~/BLiMMP/dataEdited/assemblies
echo -e 'assemblyID\tORF_count' > reports/ORF_counts.tsv
cat ~/BLiMMP/metadata/assembly_list.txt | while read assembly
do
  if [ -e ORFs/$assembly.faa ]; then
    echo "Count ORFs in" $assembly
    ORF_count=$(grep '>' ORFs/$assembly.faa | wc -l)
    echo -e "$assembly\t$ORF_count" >> reports/ORF_counts.tsv
  fi
done
