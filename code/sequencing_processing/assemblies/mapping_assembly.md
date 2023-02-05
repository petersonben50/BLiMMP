### Protocol for assembly mapping to metagenomes

This is the protocol I'm using for mapping reads to my metagenome assemblies for both 2020 and 2021.

**Generate cleaned metagenomes for mapping**

In cleaning the reads for assembly, I merged overlapping paired reads.
However, for mapping purposes this doesn't work well because it eliminates one copy of the residues in the overlapping parts.
If several genomes have been shredded more than other genomes resulting in shorter DNA segments, that could influence the coverage.
While this is pretty unlikely, I'll just avoid the issue by generating cleaned metagenomes with no merging.
This was done with fastp, a version that I downloaded and compiled for this project, with the standard read processing thresholds I've used elsewhere:

```
--detect_adapter_for_pe \
--cut_tail \
--cut_tail_window_size 10 \
--cut_tail_mean_quality 20 \
--length_required 100
```


**Pick assemblies and metagenomes for analysis**

In the `code/assemblies/mapping_stats.R` R script, I generated a table that matched the metagenomes to the assemblies for mapping purposes.
Mapping will be done by year, so all 2020 metagenomes are mapped to 2020 assemblies and 2021 metagenomes to 2021 assemblies.
TSV file saved out here: `dataEdited/assemblies/reports/mapping_key.tsv`
I then uploaded this file to GLBRC.


**Calculate total coverage for each metagenome**

I wanted to calculate the total number of nucleotides from each metagenome.
For this, I used the [readfq program](https://github.com/billzt/readfq), counting the number of nucleotides in the forward, reverse and single read files, all of which will be used for mapping.
These got saved to this file: `metagenome_coverage.tsv`.


**Map reads and process output**

First I had to build an index of each assembly using bowtie.
I did this in screen, since it doesn't take too long to complete.
Then, I checked to see which of the mapping steps I had already done (since this analysis was done in stages to a certain extent).
I saved out this list of the mapping pairs that needed to be completed yet.

Then I mapped the paired-end reads, the single reads, and the merged reads to the indices using bowtie2.
Once the mapping was done, I used samtools to convert the files to BAM files, then sorted and indexed the files.






**Rerunning mapping**

When analyzing the hgcA sequences, I noticed that the mapping from BLI20_MG_003 to the coassembly and from BLI20_MG_005 to BLI20_assembly005 was not completed.
Looking at the reports, it seems that there was not enough space for it.
So, I'm going to up the memory and run it again on glbrc9.
Turns out BLI20_MG_002 wasn't mapped to BLI20_assembly003 either, but didn't notice because there are no hgcA seqs from assembly003.
I upped the memory requested to 240GB and reran it.
