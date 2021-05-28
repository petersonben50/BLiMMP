#### code/metabolic_analyses/moor_tree.R ####
# Benjamin D. Peterson


#### Clean up on aisle R ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP")
library(ape)
library(ggtree)
library(phangorn)
library(readxl)
library(tidyverse)
library(treeio)



#### Read in tree ####

moor.tree.unrooted <- read.newick(file = "dataEdited/metabolic_analyses/moor/nitrate_reductases_moor.tree")
moor.tree <- midpoint(moor.tree.unrooted,
                      node.labels = "support")



#### Rename references ####

# Read in metadata and prepare vectors
metadata <- read.table("references/metabolic_genes/MoOR_key.tsv",
                       stringsAsFactors = FALSE,
                       sep = '\t',
                       header = FALSE)
names(metadata) <- c("geneID", "name")
ref.ID.to.name.vector <- make.unique(metadata$name)
names(ref.ID.to.name.vector) <- metadata$geneID
rm(metadata)

# Renaming references in tree
reference.index <- which(moor.tree$tip.label %in% names(ref.ID.to.name.vector))
moor.tree$tip.label[reference.index] <- ref.ID.to.name.vector[moor.tree$tip.label[reference.index]]
rm(reference.index)



#### Make color vector ####
color.vector <- rep("black", length(moor.tree$tip.label))
my.bin.index <- grep("BLI20", moor.tree$tip.label)
color.vector[my.bin.index] <- "red"



#### Visualize rooted tree ####
pdf("dataEdited/metabolic_analyses/moor/moor_tree_original.pdf",
    width = 10,
    height = 30)
ggtree(moor.tree) +
  geom_tiplab(size = 2,
              color = color.vector) +
  geom_text2(aes(subset=!isTip, label=node))
dev.off()
