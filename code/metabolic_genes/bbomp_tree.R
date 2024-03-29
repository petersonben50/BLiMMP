#### code/binning/metabolism/bbomp_tree.R ####
# Benjamin D. Peterson

#### Clean up on aisle R ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP")
library(ape)
library(ggtree)
library(patchwork)
library(phangorn)
library(readxl)
library(tidyverse)
library(treeio)
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")


#### Read in reference data ####
ref.key <- read_xlsx("/Users/benjaminpeterson/Box/references/proteins/EET/PCC/pcc_omp_key.xlsx")
rename.vector <- ref.key$organismName
rename.vector[which(!is.na(ref.key$geneName))] <- paste(rename.vector[which(!is.na(ref.key$geneName))],
                                                        " - ",
                                                        ref.key$geneName[which(!is.na(ref.key$geneName))],
                                                        sep = "")
names(rename.vector) <- ref.key$accessionID
# Refseq references
refseq.metadata <- read.table("dataEdited/metabolic_analyses/PCC/refseq_bbomp_metadata.tsv",
                              sep = '\t',
                              stringsAsFactors = FALSE)
refseq.metadata.vector <- refseq.metadata$V2
names(refseq.metadata.vector) <- refseq.metadata$V1
# Clean up
rm(ref.key, refseq.metadata)


# #### Read in sequence IDs ####
# list.of.seqs <- readLines("dataEdited/binning/metabolism/PCC/pcc_omp_custom_output.txt")
# 

#### Read in tree file ####
BBOMP.tree.unrooted <- read.newick("dataEdited/metabolic_analyses/PCC/bbomp.tree")
BBOMP.tree <- midpoint(BBOMP.tree.unrooted,
                       node.labels = "support")

#### Get indices ####
ref.index <- which(BBOMP.tree$tip.label %in% names(rename.vector))
refseq.index <- which(BBOMP.tree$tip.label %in% names(refseq.metadata.vector))
seq.index <- which(BBOMP.tree$tip.label == "KMBP009B_000000035879_4")


#### Rename references ####
BBOMP.tree$tip.label[ref.index] <- rename.vector[BBOMP.tree$tip.label[ref.index]]
BBOMP.tree$tip.label[refseq.index] <- refseq.metadata.vector[BBOMP.tree$tip.label[refseq.index]]


#### Color vector ####
tree.color.vector <- rep("grey", length(BBOMP.tree$tip.label))
tree.color.vector[seq.index] <- cb.translator["bluishgreen"]
tree.color.vector[ref.index] <- cb.translator["black"]
tree.color.vector[refseq.index] <- cb.translator["orange"]


#### Plot tree ####
BBOMP.tree.viz <- ggtree(BBOMP.tree,
                         aes(x = 0,
                             xend = 6)) +
  geom_tiplab(col = tree.color.vector)
pdf("dataEdited/metabolic_analyses/PCC/bbomp_original_tree.pdf",
    height = 10,
    width = 5)
BBOMP.tree.viz
dev.off()
