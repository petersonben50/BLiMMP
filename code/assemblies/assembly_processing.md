### Assembly processing for BLiMMP metagenomes

This document is meant to explain the scripts found in `code/assemblies/assembly_processing.sh`.
These scripts are used to process and assemble the metagenomes from the BLiMMP project.


**Upload list files**

First I uploaded the files that I needed down the line to the metadata folder:
- assembly_key.csv
- assembly_list.txt
- metagenome_list.txt

The `assembly_key.csv` file has two columns: groupID and metagenomeID.
GroupID should really be assemblyID, as this file identifies the metagenomes that are needed for each assembly.
This will be used by the `assembly_by_group.py` script for generating the assemblies.


**Transfer data**

Then I needed to retrieve the metagenomic data.

*Transfer 2020 data*

I started processing this after the 2020 season, then will just add the 2021 data later.
So, first I downloaded the data from QB3 using the `lftp` protocol they outlined for us (email on 2021-04-26).
This set of metagenomes included 5 from the BLiMMP project (BLI30_MG_001 through 005), but also included 5 metagenomes for Rachel, Elizabeth, and Charles.
The naming system that QB3 attached to the libraries was different than what I had named them, so I generated a key to link the two names here: `~/BLiMMP/dataEdited/metagenomes/reports/naming_key.tsv`.
These were all transferred into the `dataRaw/metagenomes` folder.


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
I then used the `kseq_fastq_base` script to calculate the coverage over the pre- and post-trimmed metagenomes.
These will give us some metrics for how much sequencing data we actually have and how much we lost due to trimming.


**Metagenome assembly**

*Assemblies by metaSPAdes*

I did all the assemblies for this project using metaSPAdes.
I used SPAdes v3.14.1.
Because I wanted to do a coassembly and since I hope to use duplicate metagenomes in 2021, I used a script (`assembly_by_group.py`) to run the assemblies using the needed metagenomes.
As usual, I used these kmers for assembly: 21,33,55,77,99,127.
I ran it using 24 cores with a memory limit of 1 Tb.

*Clean up assemblies*

I didn't want to include any scaffolds under 1000 bp in length, so I used anvi'o to strip out the sequences lower than 1000 bp in length.
I also simplified the scaffold names.


**Predict ORFs**

I then used a submission script to run the ORF predictions, using the same methods as used for 5M, EG, and HCC studies.
I used Prodigal (V2.6.3) on the `meta` mode and saved out fna, faa, and gff files.


**Get assembly stats**

Finally, I took a look at the assemblies.
I used the `abyss-fac.pl` script, written by [Shaun Jackman](sjackman@bcgsc.ca), to check the length and N50 values of the assemblies.
I aggregated those stats.
I then also counted the number of predicted ORFs for each assembly.
