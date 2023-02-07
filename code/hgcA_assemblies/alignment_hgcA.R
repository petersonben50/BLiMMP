#### code/binning/PVC_details/hgcB_alignment.R ####
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
readAAStringSet("dataEdited/hgcA_analysis/hgcA_raw.afa",
                format = "fasta",
                nrec = -1L,
                skip = 0L, seek.first.rec=FALSE, use.names=TRUE) %>%
  BrowseSeqs(htmlFile = "dataEdited/hgcA_analysis/hgcA_raw_hits.html")



#### Investigate alignment for final hgcA set ####
readAAStringSet("dataEdited/hgcA_analysis/identification/hgcA_muscle.afa",
                format = "fasta",
                nrec = -1L,
                skip = 0L, seek.first.rec=FALSE, use.names=TRUE) %>%
  BrowseSeqs(htmlFile = "dataEdited/hgcA_analysis/hgcA_raw_hits.html")

