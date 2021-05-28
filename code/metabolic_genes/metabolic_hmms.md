### Protocol notes for analyzing metabolic genes in assemblies

**Identification of metabolic proteins**

First I searched through the ORFs for the metabolic genes of interest.
I have accumulated these HMMs through PFAM, KEGG, and a few custom ones/other sources.
I set it up so that it groups hits within the year of the analysis, so when we get metagenomes from 2021, I can just add them in to the analysis.
I also used hmmalign here to generate alignments in case I wanted to check the quality of the hits.


**Dereplicate sequences**

I then dereplicated all the sequences.
These were done within each year since I only do the mapping within a single year of analysis.
Not sure why I feel that this is the better approach, but I do.
Again, this is set up to function this way once the 2021 metagenomes are added.


**Classify dsrA genes**

I wanted to identify which of the *dsrA* genes in the assemblies were reductive *dsrA* as opposed to reverse *dsrA* genes.
To do this, I generated a phylogenetic tree using the sequences from Karthik's dissimilatory sulfate reduction paper.


*Generate dsrA alignment*

I pulled out the *dsrA* genes I identified and aligned them together.
I then generated a consensus alignment between this alignment and an alignment of all the genes that I got from Karthik.
I trimmed this consensus alignment using `trimal`, cutting off any residue with 50% gaps.

*Generate ML tree*

I then generated a tree from this alignment using FastTree with default settings.
I downloaded this file to my local computer and looked at it in R: `code/metabolic_analyses/dsr/dsr_tree.R`.
I read in the tree and saved out a PDF of the image of the unrooted tree.
I looked for the root, which I had previously set using Candidatus_Hydrothermarchaeota_archaeon_JdFR_18_JGI24020J35080_1000005, Aigarchaeota_candidate_division_pSL4_archaeon_ASPF01000004, and Caldivirga_sp._MU80 as outgroups, to mimic the structure of the tree from Karthik's paper.
Aigarchaeota_candidate_division_pSL4_archaeon_ASPF01000004 was just outside of the cluster, so I used the other two as outgroup markers.
Node 922 was basal to these two and close to the division of Aigarchaeota_candidate_division_pSL4_archaeon_ASPF01000004, so I used that to root.
This rooted tree looked like what I was looking for.

Now I looked for the rdsr branch.
As I've done with previous projects, I looked for the "Nitrospirae_bacterium_RBG_19FT_COMBO_42_15" and "RBG_16_scaffold_151951" nodes, which marked the rdsrA branch.
Using this, I was able to find the node corresponding to the rdsrA branch in R (node 625).
I then pulled out the *dsrA* sequences from this study that were in that group.
This included 13 of the 15 *dsrA* genes that I identified, which matches up nicely with that fact that I found two *dsrD* genes.


**Identify potential PCCs**

I also wanted to search for potential porin-cytochrome c complex clusters within these metagenomes.

*First identify MHCs*

I started by looking for multiheme cytochrome c (MHC) proteins in all the assemblies.
I defined MHC here as a protein with at least 3 heme-binding sites.
I identified these using a python script that Shaomei wrote.

*Pull out adjacent genes*

The general PCC structure is that there is a MHC next to a beta-barrel outer membrane protein (BB-OMP), so I pulled out all the ORFs to either side of each MHC gene.
I first got the ID of the possible adjacent sequences, then pulled out those amino acid sequences from the ORF.faa file.

*Search adjacent genes for BBOMP*

I then used a custom HMM that I built to search these adjacent genes for BB-OMPs.
I used a cutoff of 50 for this, which is lower than expected for BB-OMPs but I'll build a phylogeny later.
I then extracted the protein sequences of the hits and dereplicated them at 97% identity.
Finally, I aligned all these hits to the custom HMM and inspected the alignment.
Alignment looks good, only BLI20_assembly004_000000000702_1 that might be a little wonky, but alignment with other seqs looks good at the N-terminus end.

*Prepare reference sequences*

I then searched for references in the NCBI RefSeq database.
I cleaned up the fasta files and dereplicated the sequences.
I then dereplicated the sequences against the BB-OMP references that I've manually curated.
Finally, I retrieved the metadata for the RefSeq sequences I ended up using.

*Generate a phylogeny*

I then aligned all the references and identified sequences together, trimmed the alignment at 50% gaps, and generated a ML tree using FastTree.
Finally, I inspected the tree in R: `code/binning/metabolism/bbomp_tree.R`.
This showed that four of the sequences form a cluster of their own, and the fifth is somewhat near a BBOMP from *Thioalkalivibrio paradoxus*.
Not super confident in these as EET mediators, but I added them to the metabolic_gene_key and thus to the abundance calculations just to see the abundance patterns.


**Check nitrate reductase genes**

I also wanted to confirm the nitrate reductase genes we identified, both the *napA* and the *narG* genes.
These are both molybdopterin oxidoreductases, so I built a tree with the MoOR reference sequences.

*Check nitrate reductase genes*

I concatenated the genes identified as *narG* and *napA* and aligned them to the MoOR HMM.
I then generated a consensus alignment with the MoOR references I gathered for the 5M project.
This alignment I masked at 50% gaps using Trimal, then generated a ML tree using FastTree on the default settings.
I downloaded the `nitrate_reductases_moor.tree` and investigated this tree using R: `code/metabolic_analyses/moor_tree.R`.
All the identified sequences fall within the expected branches on the tree, so I think we're good to move forward with what the HMMs identified.


**Clean abundance info**

Once I finished verifying the above genes, I calculated the abundance of each scaffold that included at least one metabolic gene using a submission script.
This script submitted each metagenome as a separate job.
I then continued cleaning the abundance information for the metabolic genes.
This was done in R: `code/metabolic_analyses/metabolic_proteins_depth_aggregate.R`.
This script normalized the gene coverage to the SCG coverage calculated earlier.
