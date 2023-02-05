#!/bin/sh

echo "Running fastp on" $mtID
cd $transcriptomeDirectory/workingDirectory
echo "Working in" $(pwd)
echo ""
echo "Reading in" $originalReadsDirectory/$mtID\_R1.fastq.gz
echo "and" $originalReadsDirectory/$mtID\_R2.fastq.gz
echo ""
echo "Report saved to:" $transcriptomeDirectory/reports/$mtID\_report.html
echo "Failed reads saved to:" $transcriptomeDirectory/reports/$mtID\_failed.fastq.gz
echo ""
echo "Output sent to standard output and split into files with" $linesToCut "lines"


$fastp --in1 $originalReadsDirectory/$mtID\_R1.fastq.gz \
       --in2 $originalReadsDirectory/$mtID\_R2.fastq.gz \
       --stdout \
       --failed_out $transcriptomeDirectory/reports/$mtID\_failed.fastq.gz \
       --html $transcriptomeDirectory/reports/$mtID\_report.html \
       --detect_adapter_for_pe \
       --cut_tail \
       --cut_tail_window_size 10 \
       --cut_tail_mean_quality 20 \
       --length_required 100 | \
    split -l $linesToCut - $mtID\_splitFiles_
