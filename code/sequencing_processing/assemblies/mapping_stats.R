#### code/assemblies/mapping_stats.R ####
# Benjamin D. Peterson

# This file contains the code to generate the mapping
# key file and aggregate the mapping statistics.


#### Always start with a clean slate, ya bum ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP/")
library(tidyverse)


#### Generate mapping key file ####
# Map metagenomes to assemblies from the same year.
# For now.
metagenomes <- read.table("dataEdited/metagenomes/reports/naming_key.tsv",
                          sep = '\t') %>%
  rename(metagenomeID = V2) %>%
  select(metagenomeID) %>%
  filter(grepl("BLI", metagenomeID)) %>%
  mutate(year = metagenomeID %>%
           strsplit("_") %>% sapply("[", 1) %>%
           gsub("BLI", "20", .))
assemblies <- read.table("dataEdited/assemblies/reports/all_assemblies_stats.txt",
                         header = TRUE) %>%
  select(assemblyID) %>%
  mutate(year = assemblyID %>%
           strsplit("_") %>% sapply("[", 1) %>%
           gsub("BLI", "20", .))
mapping.key <- full_join(metagenomes,
                         assemblies) %>%
  select(-year)
write.table(mapping.key,
            file = "metadata/lists/mapping_key.tsv",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
