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

Okay, we'll just keep all the hgcB hits.
Now, let's take a closer look at the genes that didn't have a hit downstream.
- BLI20_assembly004_000000034765_5: Looks like this might actually be an hgcB sequence, as it has the CIGCGMCA motif. Yeah, going to include this as likely hgcB. Very short.
- BLI20_assembly004_000000075580_2: No motifs found. This gene is on the opposite strand from the hgcA, and is truncated. Ends 13 bp down from the end of the hgcA. Classified as no hgcB.
- BLI20_assembly004_000000098923_1: Also very short (MSPLPTLRYIETAVTLQLDP) and has no hits in MOTIF. Overlaps with hgcA. Not out of the question that this is hgcB actually based on the similarity to the hgcB alignments. I would say this one is small enough that we could say the downstream gene is truncated and therefore can't comment on presence of hgcB.
- BLI20_coassembly_000000053122_2: Big gene here. Elongation factor Tu GTP binding domain. Starts 200 bp downstream of hgcA. Unlikely that there's an unseen hgcB in there, but worth looking.

*Isolate gene neighborhoods*

Now, let's look at the gene neighborhoods, see if we can learn anything more from that.
I isolated the scaffold containing hgcA with 5000 bp upstream of the start and 5000 bp downstream of the end using a custom script that I adopted from a script that Tyler Barnum sent me.
I downloaded the fasta file and GFF and loaded it into Geneious.
Let's take a look.

First I looked at the four genes we looked at above, that had a gene downstream that wasn't a hit for hgcB.
BLI20_assembly004_000000034765 does indeed end with the downstream gene being cut off, as does BLI20_assembly004_000000098923.

BLI20_assembly004_000000075580 does a gene that's reversed.
I did a nucleotide blast of the downstream DNA using blastn on the NCBI website.
Shares ~75-78% identity with Dehalobacter or Geobacter genomes, but no specific genes.
Tried running blastx, using the following as an entry:
- AAGCTTTCATTACCTGCGTAATGGTTCGAGCCTGGTGGTGGATGCCGGGCGCTGTACTGGCTGCTTGGCCTGTCTGGAGGTGTGCCCGCACGGCGTTCTGGCGGCGGATGGGCCGCCGGTTGCCGGTGGTGGTTCCGGTAGCCGGTTGGCGGTGCTGGTAGCCGACCGGCCGGCCTGCATGGAATGCGGCGCCTGCGCCCGTAATTGTCCCGCCGGGGCGATTAGTGTGCGCAGCGGTGTGGGCTGCGCGTACGCCATTATCCGGGGCAAGTTGCGCGGTACGGCTCCGGATTGTTCGTGCGGCTGCGGTACGCAGTCCGGTTGTTGGTGACTGGTGGCCCCAGGCGGCGCGCTTG

This actually hit against some ferrodoxins.
In Geneious, I extracted this sequence.
I then did a 6 frame translation and aligned all those translations to the hgcB alignment.
One of them matched well (BLI20_assembly004_000000075580 extraction translation frame 2).
It has the expected motif, so yes, I included this one as having a downstream hgcB.

For BLI20_coassembly_000000053122, I skipped straight to doing the translation of the downstream nucleotides in Geneious.
Surprisingly enough, this also spawned a hit to the hgcB alignment, including the necessary motif (BLI20_coassembly_000000053122 extraction translation frame 1).

So, I'm just going to count all 4 of these as likely having a downstream hgcB gene.
So that's 42 with a downstream gene.
The other 15 hgcA seqs didn't have a downstream gene, so we'll assume those scaffolds were truncated.
While I will list all 42 as having downstream hgcB genes, I'll only include the hgcB sequences from those that were properly predicted as hgcB in the final hgcB.faa file that gets posted.
To make this easier, I'll use the `downstream_genes_present.txt` as the list for hgcB genes.

**Dereplication**

Next I wanted to dereplicate these sequences across the different assemblies.
For now, I only have data from 2020, but eventually I'll be analyzing 2020 and 2021 together.
So, I'll try to write it so that I can easily transition to that when I have to redo the analysis, but for now will focus on doing the analysis for 2020.
I used CD-HIT (v4.6) to cluster the hgcA sequences at 97% identity, with the default cutoff of 50% alignment.
Using a cut-off of 95% also results in 32 clusters (out of 57 sequences), but dropping to 93% drops it to 30 clusters.
We'll stick with the 97% cutoff here.

I then downloaded the `hgcA_good_acrossYear.tsv` file to my local computer and read it into R (`code/hgcA_analysis/hgcA_dereplication.R`).
Here, I combined the hgcB data, the classification info, and the dereplication info.

**When I do this with the 2021 data for a final version of this to be published, I'll need to look more closely at which hgcA seqs we're using for phylogenetic analysis. For now, we'll roll with the auto derep data, to save time that is going to be duplicated down the line.**

I uploaded `hgcA_final_list.txt` and `hgcA_final_abundance_list.txt` to the main `hgcA_analysis` folder on GLBRC.
We'll use these for the phylogenetic tree and the abundance calculation, respectively.
