#### code/binning/metabolism/MoORs.R ####
# Benjamin D. Peterson


#### Get set up ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP")
library(dplyr)
library(ggimage)
library(ggtree)
library(pdftools)
library(phangorn)
library(readxl)
library(treeio)

#### Generate bin/gene conversion vectors ####
G2B <- read.table("dataEdited/binning/metabolism/MoORs/MoOR_G2B.tsv",
                  stringsAsFactors = FALSE,
                  sep = "\t",
                  header = FALSE)
names(G2B) <- c("geneID", "binID")
# Gene to bin vector
gene.to.bin.vector <- G2B$binID
names(gene.to.bin.vector) <- G2B$geneID

# Bin to gene vector
bin.to.gene.vector <- G2B$geneID
names(bin.to.gene.vector) <- G2B$binID

rm(G2B)


#### Read in tree ####
MoOR.tree <- read.newick("dataEdited/binning/metabolism/MoORs/putative_MoORs_1.tree")



#### Generate taxonomy label for each bin ####
bin.taxonomy <- read.table("dataEdited/binning/taxonomy/gtdbtk.bac120.summary.tsv",
                           sep = "\t",
                           header = TRUE) %>%
  filter(user_genome %in% gene.to.bin.vector[MoOR.tree$tip.label]) %>%
  rename(binID = user_genome) %>%
  mutate(phylum = classification %>%
           strsplit("p__") %>% sapply("[", 2) %>%
           strsplit(";c__") %>% sapply("[", 1)) %>%
  select(binID, phylum)
bin.taxonomy.vector <- bin.taxonomy$phylum
names(bin.taxonomy.vector) <- bin.taxonomy$binID
rm(bin.taxonomy)



#### Generate renaming vector ####
renaming.vector <- MoOR.tree$tip.label
names(renaming.vector) <- MoOR.tree$tip.label

# For my bins
my.bin.index <- which(renaming.vector %in% bin.to.gene.vector)
renaming.vector[my.bin.index] <- paste(renaming.vector[my.bin.index], " - ",
                                       gene.to.bin.vector[renaming.vector[my.bin.index]],
                                       " (", bin.taxonomy.vector[gene.to.bin.vector[renaming.vector[my.bin.index]]], ")",
                                       sep = "")



#### Generate metadata vector ####
metadata <- read.table("references/MoORs/MoOR_key.tsv",
                       stringsAsFactors = FALSE,
                       sep = '\t',
                       header = FALSE)
names(metadata) <- c("geneID", "name")
ref.ID.to.name.vector <- make.unique(metadata$name)
names(ref.ID.to.name.vector) <- metadata$geneID
ref.name.to.ID.vector <- make.unique(metadata$geneID)
names(ref.name.to.ID.vector) <- metadata$name
#rm(metadata)

reference.index <- which(renaming.vector %in% names(ref.ID.to.name.vector))
renaming.vector[reference.index] <- ref.ID.to.name.vector[renaming.vector[reference.index]]

rm(reference.index)

#### Make color vector ####
color.vector <- rep("black", length(MoOR.tree$tip.label))
color.vector[my.bin.index] <- "red"


#### Write out new tree ####
MoOR.tree.for.printing <- MoOR.tree
MoOR.tree.for.printing$tip.label <- renaming.vector[MoOR.tree$tip.label]
write.tree(MoOR.tree.for.printing,
           "dataEdited/binning/metabolism/MoORs/FastTree_MoOR_edited.tree")


#### Save PDF of tree ####
pdf("dataEdited/binning/metabolism/MoORs/MoOR.tree.pdf",
    height = 30,
    width = 8)
ggtree(MoOR.tree.for.printing,
       aes(x = 0, xend = 6)) + 
  geom_tiplab(size = 2,
              col = color.vector)
dev.off()

