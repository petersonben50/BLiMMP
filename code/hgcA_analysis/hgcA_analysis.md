#### hgcA analysis

This is the markdown script walking through the analysis of the *hgcA* content of the metagenomes.


**Identify hgcA sequences**

*Identify putative hgcA genes with HMM*

I used the HMM I built for the 5M project to search through the open reading frames from the assemblies for HgcA-like sequences.
These results were stored in folders by year.
At this point, it's just for 2020.
I then pulled out the amino acid sequences for each of the putative HgcA sequences.

*Concatenate and align all hgcA seqs for curation*

I then collected all the HgcA amino acid sequences and aligned them against the HgcA HMM.
This alignment I downloaded to my computer and inspected in Geneious.
To include the sequence in analysis, it needed to have the cap helix domain and predicted transmembrane domains, or at least the length at the C-terminus end that normally is occupied by the transmembrane domains.

Removal notes:
- 001: Removed two seqs, one for missing cap helix and one for missing TM domains. Two more seqs are truncated at N-terminus, but included them here anyways.
- 002: Removed two seqs for missing cap helix. Removed 4 for not having the TM domain. There's one that's pretty truncated at the C-terminus, only one TM domain. Few others that appear truncated at N-terminus.
- 003: One hit, and it didn't have any predicted TM domains, so cut it.
- 004: One seq removed for missing predicted TM domains.
- 005: Two seqs removed for missing predicted TM domains.
- coassembly: Five seqs removed for missing predicted TM domains, and three removed for missing cap helix domain. Two sequences were truncated but included (one truncated at C-terminus and one at N-terminus).

I then saved out this curated alignment as `hgcA_good.faa`.



**Classify hgcA seqs with pplacer workflow**

I used a workflow developed by Caitlin Gionfriddo to classify the *hgcA* genes, so all credit to her for this really nice workflow.
A very helpful tutorial can be found here: https://caitlingio.com/tutorial-for-hgcab-amplicon-sequencing-data/.
It relies on the Hg-MATE database (I used version 1 for this analysis).
First is the generation of an alignment of the HgcA sequences from this study using MUSCLE.
Then, consensus alignment is made between this alignment and the alignment of the Hg-MATE database.
The sequences are placed on a pre-generated tree using pplacer.
The program rppr is used to make a sqlite database of the taxonomy of the reference sequences.
Guppy was then used to classify the sequences and visualize the placements on a tree.

I then downloaded the files I needed (`Hg_MATE_classify`, `hgcA_for_classification.jplace`) to my local computer and used the R script `clean_hgcA_classification.R` (also adapted from Caitlin's well-annotated workflow) to clean up the classification and save out a csv file with that information.



**Genomic context for hgcA**

*First pull out hgcA+ scaffolds*

Next I identified the scaffolds that contained *hgcA* and pulled out the fna files and the GFF entries.

*Search downstream genes for hgcB*

Using the GFF entries, I identified the ORF downstream from each hgcA sequence (if there was one on the scaffold).
I pulled out all these ORFs, then searched through them using an HMM I built from hgcB genes for the 5M project to see if they were hgcB.
I pulled out the hits and aligned them to the HMM.
I downloaded the alignment to my local computer and inspected it in Geneious.

There were 38 hits, out of 42 identified downstream genes.
All but 3 of these hits had the expected C(M/I)ECGAC motif that confirmed it was hgcB.
The other three were cut off before that motif.
These sequences were all truncated though, according to the GFFs, so we'll include them as true hgcB seqs since the alignment at the N-terminus end was pretty good.
Two of these genes, BLI20_assembly004_000000016181_3 and BLI20_coassembly_000000077989_5 had a long trailing sequence at the C-terminus.
The trailing sequence was identical between the two and are from two different assemblies, so it might be from the same organism.
This trailing sequence includes a segment that had a hit for a thioredoxin motif using MotifFinder, but the score is very low.
Not sure what this is here, but we'll include it for now.
