### Assembly processing for BLiMMP metagenomes

This document is meant to explain the scripts found in `code/assemblies/assembly_processing.sh`.
These scripts are used to assemble the metagenomes from the BLiMMP project.


**Prepare for assembly**

The first task I needed to do was assign each metagenome to an assembly.
I chose to assemble the metagenomes by sample location, meaning that the two duplicate metagenomes from the 2021 samples were coassembled, while the 2020 metagenomes were assembled individually.
Coassemblies were generated from both years.
This information is stored in a csv file here: `metadata/lists/assembly_key.csv`.
This file has two columns: groupID and metagenomeID.
GroupID should really be assemblyID, as this file identifies the metagenomes that are needed for each assembly.
This will be used by the `assembly_by_group.py` script for generating the assemblies.


**Metagenome assembly**

*Assemblies by metaSPAdes*

I did all the assemblies for this project using metaSPAdes (SPAdes v3.14.1).
I used a python script (`assembly_by_group.py`) to run the assemblies using the needed metagenomes by group.
This was submitted to a condor cluster and is called by the `assembly_by_group_execute.sh` execution file and the `assembly_by_group_submit.sub` submission file.
As usual, I used these kmers for assembly: 21,33,55,77,99,127.
I ran it using 24 cores with a memory limit of 1 Tb, but ended up having to do a separate submission for the 2021 coassembly with 40 CPUs and 1.5 Tb.

*Clean up assemblies*

I didn't want to include any scaffolds under 1000 bp in length, so I used anvi'o to strip out the sequences lower than 1000 bp in length.
I also simplified the scaffold names to use the assembly name as a pre-fix.
This was all included in the `assembly_by_group_execute.sh` file.


**Predict ORFs**

I then used a submission/execution script (prefix: `ORF_prediction_assemblies`) to run the ORF predictions, using the same methods as used for 5M, EG, and HCC studies.
I used Prodigal (V2.6.3) on the `meta` mode and saved out fna, faa, and gff files.
I then cleaned up the ORF files (the faa and fna ones) using the `cleanFASTA.py` script (also in the `ORF_prediction_assemblies.sh` executable).


**Get assembly stats**

Finally, I took a look at the assemblies.
I used the `abyss-fac.pl` script, written by [Shaun Jackman](sjackman@bcgsc.ca), to check the length and N50 values of the assemblies.
I aggregated those stats (`all_assemblies_stats.txt`).
I then also counted the number of predicted ORFs for each assembly (`ORF_counts.tsv`).
I downloaded these files to my local computer (`dataEdited/assemblies/all_assemblies_stats.txt`).
