#!/bin/sh


######################
# Assemblies by metaSPAdes
######################
cd /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/assemblies
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=""
PERL5LIB=""


if [ ! -d $output/$assembly ]; then
  mkdir $output/$assembly
fi
if [ ! -e $output/$assembly/scaffolds.fasta ]; then
  echo "Assembling" $assembly
  python $code/assembly_by_group.py --groupName $assembly \
                                    --groupSpreadsheetName $assembly_grouping \
                                    --readLocation $read_storage \
                                    --mergedReads yes \
                                    --output $output/$assembly \
                                    --threads 36 \
                                    --memory 1400
else
  echo $assembly "already assembled"
fi
conda deactivate



######################
# Clean up assemblies
######################
conda activate anvio6.2
PYTHONPATH=''

if [ -e scaffolds/$assembly\_assembly.fna ]; then
  echo "Assembly for" $assembly "is already cleaned."
else
  echo "Cleaning the assembly for" $assembly
  anvi-script-reformat-fasta assembly_files/$assembly/scaffolds.fasta \
                              -o scaffolds/$assembly\_assembly.fna \
                              -l 1000 \
                              --simplify-names \
                              --prefix $assembly \
                              --report-file scaffolds/renaming_reports/$assembly\_report_file.txt
fi
