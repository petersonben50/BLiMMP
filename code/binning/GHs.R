#### code/binning/metabolism/GHs.R ####
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("~/Documents/research/BLiMMP")
library(lubridate)
library(readxl)
library(tidyverse)


#### Read in GH data ####
GH.data <- read.table("dataEdited/binning/metabolism/GHs/GHs_bins.tsv",
                      sep = "\t",
                      header = TRUE) %>%
  rename(geneID = GeneID,
         GH_ID = HMMER) %>%
  select(geneID, GH_ID) %>%
  mutate(GH_class = substr(GH_ID, 1, 2))


#### Read in G2B file ####
G2B.df <- read.table("dataEdited/binning/binsFinal_G2B.tsv",
                     col.names = c("geneID", "binID"))


#### Combine data and aggregate counts ####
GH.count.data <- GH.data %>%
  left_join(G2B.df) %>%
  group_by(binID, GH_class) %>%
  summarize(GH_count = n()) %>%
  ungroup()


#### Add taxonomy data ####
tax.data <- read.table("dataEdited/binning/taxonomy/gtdbtk.bac120.summary.tsv",
                       sep = '\t',
                       header = TRUE) %>%
  rename(binID = user_genome) %>%
  select(binID, classification)
GH.count.data <- GH.count.data %>%
  left_join(tax.data) %>%
  spread(key = GH_class,
         value = GH_count)
write.csv(GH.count.data,
          file = "dataEdited/binning/metabolism/GHs/clean_cazyme_data.csv",
          row.names = FALSE)
