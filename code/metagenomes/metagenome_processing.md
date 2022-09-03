### Metagenome processing for BLiMMP project

This document is meant to explain the scripts found in `code/metagenomes/metagenome_processing.sh`.
These scripts are used to process the metagenomes from the BLiMMP project.


**Transfer data**

First I needed to retrieve the metagenomic data.

*Transfer 2020 data*

I started processing this after the 2020 season, then will just add the 2021 data later.
So, first I downloaded the data from QB3 using the `lftp` protocol they outlined for us (email on 2021-04-26).
This set of metagenomes included 5 from the BLiMMP project (BLI30_MG_001 through 005), but also included 5 metagenomes for Rachel, Elizabeth, and Charles.
The naming system that QB3 attached to the libraries was different than what I had named them, so I generated a key to link the two names here: `~/BLiMMP/dataEdited/metagenomes/reports/naming_key.tsv`.
These were all transferred into the `dataRaw/metagenomes` folder.
Later, I changed this so that the metagenomes got imported into a `dataRaw/metagenomes/2020_data` folder.
Then, I moved the metagenomes up a folder.


*Transfer 2021 data*

First downloaded the data from QB3 using the `lftp` protocol they outlined for us (email on 2021-07-22).
Email here: `dataRaw/metagenomes/2021_data_email.pdf`.
For this, we generated 12 metagenomes, duplicate metagenomes from the six sites where we did incubations.
Transferred them into `2021_data`.
Looks like we also got a bunch of empty files, all the ones with the "L001".
These files aren't listed in the md5 files they sent along.
So, let's remove them, only work with the "L002" files.


**Generate and upload the naming file**

The names that we submit to QB3 and the ones they ultimately use in their sequencing don't match up with our naming convention.
Thus, I generated a two-column table that contains the sequencing ID used by QB3 in the first column and the desired metagenome ID that we want in the second column.
This is found here: `dataEdited/metagenomes/reports/naming_key.tsv`
Upload to GLBRC here: `~/BLiMMP/dataEdited/metagenomes/reports/naming_key.tsv`.
Note that 5 of the metagenomes included here with the 2020 sequencing were not from this project.
They will not be processed further in this workflow.


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

I then downloaded these to my local computer and started to look at them.


**Run Mash on metagenomes**
