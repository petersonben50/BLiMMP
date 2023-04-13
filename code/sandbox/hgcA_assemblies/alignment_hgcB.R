#### code/hgcA_analysis/hgcB_alignment.R ####
# Benjamin D. Peterson

# This script will generate images of
# the phylogenetic trees for the PVC bins.


#### Always start with a clean slate ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP/")
library(Biostrings)
library(DECIPHER)
library(tidyverse)


#### Investigate alignment for hgcA hits ####
readAAStringSet("dataEdited/hgcA_analysis/hgcB/hgcB_clean.afa",
                format = "fasta",
                nrec = -1L,
                skip = 0L,
                seek.first.rec = FALSE,
                use.names = TRUE) %>%
  BrowseSeqs(htmlFile = "dataEdited/hgcA_analysis/hgcB/hgcB_raw_hits.html")
