#!/bin/sh

######################
# code/assemblies/mapping_clean_metagenome.sh
# Benjamin D. Peterson
######################

if [ ! -e $output_folder/$metagenomeID\_R1.fastq.gz ]; then
  echo "Processing" $metagenomeID "reads for mapping:"
  ls $input_folder/$metagenomeID\_*R1*fastq.gz
  ls $input_folder/$metagenomeID\_*R2*fastq.gz
  $fastp --in1 $input_folder/$metagenomeID\_*R1*fastq.gz \
          --in2 $input_folder/$metagenomeID\_*R2*fastq.gz \
          --out1 $output_folder/$metagenomeID\_R1.fastq.gz \
          --out2 $output_folder/$metagenomeID\_R2.fastq.gz \
          --unpaired1 $output_folder/$metagenomeID\_single.fastq.gz \
          --unpaired2 $output_folder/$metagenomeID\_single.fastq.gz \
          --failed_out $report_folder/$metagenomeID\_failed.fastq.gz \
          --html $report_folder/$metagenomeID\_report.html \
          --detect_adapter_for_pe \
          --cut_tail \
          --cut_tail_window_size 10 \
          --cut_tail_mean_quality 20 \
          --length_required 100
else
  echo "Already processed" $metagenomeID
fi
