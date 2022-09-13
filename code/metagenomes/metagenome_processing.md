### Metagenome processing for BLiMMP project

This document is meant to explain the scripts found in `code/metagenomes/metagenome_processing.sh`.
These scripts are used to process the metagenomes from the BLiMMP project.


**Trim metagenomes**

Next I needed to quality trim the metagenomes.
I used fastp for this.
I had been having issues with using the conda version of fastp, so I downloaded the binary instead.
It was downloaded 2021-04-29, version 0.21.0.
I processed the reads in the `dataRaw/metagenomes` folder and saved the output in `dataEdited/metagenomes` with the original names of the samples that I had generated.
I cut the tails of the reads at a mean quality of 20 over a 10 bp sliding window, then removed any reads that had less than 100 bp in them.
I saved all the newly unpaired reads that occurred after trimming.
I merged reads that were overlapping.
I also automatically detected adaptors and had them trimmed off.



**Check size of metagenomes**

Next I wanted to generate a couple of metrics for the different metagenomes.
First I counted the reads in both the trimmed and untrimmed metagenomes using `zgrep` to count the lines starting with "@".
These will give us some metrics for how much sequencing data we actually have and how much we lost due to trimming.


**Data aggregation and visualization**


I then downloaded these to my local computer and started to look at them.
Data here: `code/metagenomes/data_aggregation_metagenomes.R`



**Run Mash on metagenomes**

I used Mash version 2.3.
First I sketched the metagenomes using `mash sketch`.
I did each metagenome individually, so that later on I could analyze them in any combination that I liked, rather than generating one sketch file.
I ran it with 12 threads (`-p 12`).
Send output to `sketch_files` folder, with metagenomeID as prefix (`-o sketch_files/$metagenomeID`).
I used a kmer size of 21 (`-k 21`).
I used a sketch size of 100,000 (`-s 100000`).
The default is much lower (1000), but I think with metagenomes as opposed to genomes that might not be enough resolution.
Started with a seed score of 50 (`-S 50`).
Specified that our input is a read set (`-r`).
I also filtered out unique kmers by requiring that kmers be present at least twice (`-m 2`).
I only did the sketch on the forward reads for now.
After the sketches were complete, I combined them into one sketch (`BLI_MG_sketches.msh`) and calculated the distance between the metagenomes
