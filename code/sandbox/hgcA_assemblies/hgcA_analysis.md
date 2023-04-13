#### hgcA analysis

This is the markdown script walking through the analysis of the *hgcA* content of the metagenomes.


**Identify hgcA sequences**

*Identify putative hgcA genes with HMM*

I used the HMM I built for the 5M project to search through the open reading frames from the assemblies for HgcA-like sequences.
These results were stored in folders by year, since I only mapped metagenomes to assemblies from the same year.
Thus, depth calculations will be by year.
I then pulled out the amino acid sequences for each of the putative HgcA sequences.
No reads were identified from BLI21_assembly101.

*Concatenate and align all hgcA seqs for curation*

I then collected all the HgcA amino acid sequences and aligned them against the HgcA HMM.
This alignment I downloaded to my computer and inspected using Geneious.
Don't have the paid version, so this is a little wonky.
To include the sequence in analysis, it needed to have the cap helix domain and predicted transmembrane domains, or at least the length at the C-terminus end that normally is occupied by the transmembrane domains.
I saved my notes on this here: `dataEdited/hgcA_analysis/identification/hgcA_raw_inclusion_notes.xlsx`.
The sequence IDs of the sequences to be removed were copied here: `dataEdited/hgcA_analysis/identification/seqs_to_remove.txt`.
I then uploaded this file to GLBRC.

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



**Pull out depth of hgcA+ scaffolds**

I then pulled out the depths of the hgcA+ scaffolds.
This was set up to run with multiple years if need be, so once we get 2021 metagenomes I'll be set to run this.
I calculated the depth as the average coverage over each nucleotide in the scaffold except for the nucleotides within 150 bp of either end of the scaffold.

I downloaded the aggregated depths to my local computer.
Then I aggregated the data using a R script: `code/hgcA_analysis/clean_hgcA_abundance.R`.
This data was normalized to the SCG coverage I calculated for each metagenome.


**Genomic context for hgcA**

*First pull out hgcA+ scaffolds*

Next I identified the scaffolds that contained *hgcA* and pulled out the fna files and the GFF entries.

*Search downstream genes for hgcB*

Using the GFF entries, I identified the ORF downstream from each hgcA sequence (if there was one on the scaffold).
I pulled out all these ORFs, then searched through them using an HMM I built from hgcB genes for the 5M project to see if they were hgcB.
I pulled out the hits and aligned them to the HMM.
I downloaded the alignment to my local computer and inspected using the `alignment_hgcB.R` script.

There were 148 hits, out of 169 identified downstream genes.
All but 6 of these hits had the expected C(M/I)ECGAC motif that confirmed it was hgcB.
The other six were cut off before that motif.
These sequences were all at the end of the scaffold according to the GFFs, suggesting that they are just truncated.
Thus, we'll include them as true hgcB seqs since the alignment at the N-terminus end was pretty good.
Three of these genes, BLI21_assembly106_000000092827_1, BLI20_assembly004_000000016181_3 and BLI20_coassembly_000000077989_5 had a long trailing sequence at the C-terminus.
The trailing sequence was identical between the two that were from BLI20 and are from two different assemblies, so it might be from the same organism.
The sequence from BLI21 is not quite as long as the ones from the BLI20 assembly.
According to the auto-classification, they're associated with the PVC superphylum.
This trailing sequence includes a segment that had a hit for a thioredoxin motif using MotifFinder, but the score is very low.
Not sure what this is here, but we'll include it for now.
In all, we just kept all the hgcB hits.

We'll take a closer look at the downstream genes that didn't hit the hgcB HMM later.
I inspected them manually and ran them through [MOTIF](https://www.genome.jp/tools/motif/).
Notes on this here: ``



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


**Phylogenetic analysis of hgcA**

Next I wanted to generate a good phylogenetic tree.

*Generate rough tree with Hg-MATE seqs using FastTree*

First, I wanted to just get a sense for it using the Hg-MATE database with FastTree.
I pulled out the hgcA representatives from my data set, aligned them, then generated a consensus alignment with the Hg-MATE database.
I used FastTree with the default settings to make a ML tree.

I downloaded this tree to my local computer, then looked at in R.
Scripts are here: `tree_hgcA_FastTree.R`.
Wow, this is really dominated by PVC microbes, many of them close to what we found last time.
They're even more numerically dominant than in the 2017 data, I think.
So cool!


*Generate good tree with Hg-MATE seqs using RAxML*

They're concentrated enough into a few smaller clades that we can probably just manually select sequences to use for a RAxML tree.
I'll also include all the hgcA sequences from the initial 2017 study (so I won't pull any of those out of the FastTree tree).
Hg-MATE names of references to be used can be found here: `dataEdited/hgcA_analysis/phylogeny/reference_names_to_use.txt`.
I used this list to pull out the sequences I needed, locally.
Needed to do a couple of iterations, since there were a few sequences that were in the refpackage but not in the database.
I think Caitlin had updated the refpackage separately or something.
Either way, good to go now.
I uploaded the references to GLBRC.

I also pulled the references from my study and from Jones et al, 2019.
I had to do a bit of curation on my sequences, man was I bad at keeping records.
I got the Jones sequences from the 5M study references.
I also added two hgcA paralogs, to use in rooting the tree.
I concatanated all these sequences and aligned them using MUSCLE.
I downloaded the alignment and inspected it in Geneious to ensure it looks good.
Which it does.
I then masked the alignment at 50% gaps, which seems to work well.
I exported it (`hgcA_for_tree_masked.afa`) and uploaded it to the GLBRC.
I then ran RAxML (v8.2.11) on it to generate a ML tree.
I used rapid bootstrapping with automatic detection of limits and autodetection of the mutation model.
