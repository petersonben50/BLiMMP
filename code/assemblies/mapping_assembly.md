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
