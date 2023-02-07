#### code/hgcA_analysis/tree_hgcA_FastTree.R ####
# Benjamin D. Peterson

# This script will look at the ML tree
# we generated using FastTree with the
# Hg-MATE database.
# It will save out a PDF of that for
# inspection.



#### Always start with a clean slate ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP/")
library(ggtree)
library(phangorn)
library(readxl)
library(tidyverse)
library(treeio)
cb.translator <- readRDS("~/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")



#### Check out FastTree of final alignment ####

# Read in tree
hgcA.tree.name <- "dataEdited/hgcA_analysis/phylogeny/rough_hgcA.tree"
hgcA.tree.unrooted <- read.newick(hgcA.tree.name)
rm(hgcA.tree.name)

# Check out unrooted tree
pdf("dataEdited/hgcA_analysis/phylogeny/hgcA_tree_FastTree_unrooted.pdf",
    height = 120,
    width = 8)
ggtree(hgcA.tree.unrooted,
       aes(x = 0,
           xend = 9)) +
  geom_tiplab(size=2.5, align = TRUE) +
  geom_text2(aes(subset=!isTip, label=node),
             size = 1.5)
dev.off()

# Branch leading to paralogs is 1320



#### Root tree ####
hgcA.tree <- root(hgcA.tree.unrooted,
                  node = 1320,
                  edgelabel = TRUE)



#### List of hgcA sequences ####
hgcA.list <- readLines("dataEdited/hgcA_analysis/dereplication/hgcA_final_list.txt")
this.study.indices <- which(hgcA.tree$tip.label %in% hgcA.list)



#### Set color vector ####
color.vector <- rep("black", length(hgcA.tree$tip.label))
color.vector[this.study.indices] <- "red"

#### Visualize tree ####
pdf("dataEdited/hgcA_analysis/phylogeny/hgcA_tree_FastTree.pdf",
    height = 120,
    width = 8)
ggtree(hgcA.tree) + 
  geom_tiplab(size=2.5,
              align = TRUE,
              col = color.vector) +
  geom_text2(aes(subset = !isTip,
                 label = node),
             size = 1.5)
dev.off()
