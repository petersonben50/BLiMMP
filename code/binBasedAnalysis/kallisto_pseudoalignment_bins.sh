#!/bin/sh

##############################
# code/executables/kallisto_pseudoalignment.sh
# Benjamin D. Peterson

# This is a script that will pseudoalign
# reads from a transcriptome to a specified
# set of ORFs that have been indexed.
##############################

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
PYTHONPATH=""
conda activate kallisto

# For anvio bins
if [ -f $indexFolder/$assemblyID\_anvio.idx ]; then

    echo "Pseudoaligning" $mtID "to" $assemblyID "for anvio"
    kallisto quant \
            -i $indexFolder/$assemblyID\_anvio.idx \
            -o $outputFolder/$mtID\_to_$assemblyID\_anvio_output \
            -b 100 \
            -t 7 \
            $readsLocation/$mtID\_mRNA_R1.fastq \
            $readsLocation/$mtID\_mRNA_R2.fastq
else
echo "No anvio ORFs for" $assemblyID
fi
# Gather and rename abundance files
cp $outputFolder/$mtID\_to_$assemblyID\_anvio_output/abundance.tsv $outputFolder/$mtID\_to_$assemblyID\_kallisto_anvio.tsv
cp $outputFolder/$mtID\_to_$assemblyID\_anvio_output/run_info.json $reportFolder/$mtID\_to_$assemblyID\_run_info_dasTool.json
rm -rf $outputFolder/$mtID\_to_$assemblyID\_anvio_output



# For Das Tools bins
if [ -f $indexFolder/$assemblyID\_dasTool.idx ]; then

    echo "Pseudoaligning" $mtID "to" $assemblyID "for Das Tools"
    kallisto quant \
            -i $indexFolder/$assemblyID\_dasTool.idx \
            -o $outputFolder/$mtID\_to_$assemblyID\_dasTool_output \
            -b 100 \
            -t 7 \
            $readsLocation/$mtID\_mRNA_R1.fastq \
            $readsLocation/$mtID\_mRNA_R2.fastq
else
echo "No Das Tool ORFs for" $assemblyID
fi
# Gather and rename abundance files
cp $outputFolder/$mtID\_to_$assemblyID\_dasTool_output/abundance.tsv $outputFolder/$mtID\_to_$assemblyID\_kallisto_dasTool.tsv
cp $outputFolder/$mtID\_to_$assemblyID\_dasTool_output/run_info.json $reportFolder/$mtID\_to_$assemblyID\_run_info_dasTool.json
rm -rf $outputFolder/$mtID\_to_$assemblyID\_dasTool_output
