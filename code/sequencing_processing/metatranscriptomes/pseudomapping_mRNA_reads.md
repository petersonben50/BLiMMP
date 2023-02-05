### Protocol for mapping to metatranscriptomes to ORF fna files

This is the protocol I'm using for mapping metatranscriptome reads to the assembled ORFs for both 2020 and 2021.

**Pair assembly ORFs and metatranscriptomes for analysis**

Like with the metagenome assemblies, we'll do metatranscriptome mapping only within the same year.
In the `code/metatranscriptomes/housekeeping_metatranscriptomes.R` R script, I generated a table that matched the metatranscriptomes to the assemblies for mapping purposes.
Mapping will be done by year, so all 2020 metatranscriptomes are mapped to 2020 assembly ORFs and 2021 metatranscriptomes to 2021 assembly ORFs.
TSV file saved out here: `dataEdited/metatranscriptomes/reports/pseudomapping_key.tsv`
I then uploaded this file to GLBRC.


**Pseudoalign RNA reads to references**

Next step was to pseudoalign the RNA reads to the ORFs using the kallisto software.
This was done in a unique conda environment.

First step was to index the ORF files.
I ran this as a loop on GLBRC 9.
Still running in morning, should have submitted it.

Running the pseudoalignment as a submitted file.
