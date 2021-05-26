### Protocol to identify and quantify SCGs in Hells Canyon metagenomes

**Get set up**

First I needed to set up the directories for the data and for the reports.

**Submit jobs**

I wrote these analyses as submission files.
Each job goes after a different gene.
First it pulls the copies of the gene out of all the different assemblies, then dereplicates those genes at a 97% identity cutoff.
Finally, it calculates the coverage of all the identified genes in all the different metagenomes.

**Concatenate files**

Concatenate all the coverage files into one file.

**Clean up**

Remove the assembly/metagenome lists that I generated in this step.
