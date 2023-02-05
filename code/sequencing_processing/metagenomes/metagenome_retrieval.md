### Metagenome retrieval for BLiMMP project

This document is meant to explain the scripts found in `code/metagenomes/metagenome_retrieval.sh`.
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
I used this to rename the metagenomes.
Note that 5 of the metagenomes included here with the 2020 sequencing were not from this project.
They will not be processed further in this workflow.


**Clean up!**

I downloaded the additional data contained in the `2020_data` and `2021_data` folders.
I then deleted those additional files and moved the metagenomes for other people's projects to a different folder (`~/collaborators/2020_other_metagenomes`).

I then uploaded the metagenomes to NCBI, as described here: `code/metagenomes/metagenome_upload_to_NCBI.docx`
