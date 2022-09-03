### Assembly processing for BLiMMP metagenomes

This document is meant to explain the scripts found in `code/assemblies/assembly_processing.sh`.
These scripts are used to assemble the metagenomes from the BLiMMP project.



- assembly_key.csv
- assembly_list.txt
The `assembly_key.csv` file has two columns: groupID and metagenomeID.
GroupID should really be assemblyID, as this file identifies the metagenomes that are needed for each assembly.
This will be used by the `assembly_by_group.py` script for generating the assemblies.







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
I then cleaned up the ORF files (the faa and fna ones) using the `cleanFASTA.py` script.


**Get assembly stats**

Finally, I took a look at the assemblies.
I used the `abyss-fac.pl` script, written by [Shaun Jackman](sjackman@bcgsc.ca), to check the length and N50 values of the assemblies.
I aggregated those stats.
I then also counted the number of predicted ORFs for each assembly.
