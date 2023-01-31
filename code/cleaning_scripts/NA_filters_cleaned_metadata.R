#### code/metagenome_metadata.R ####
# Benjamin D. Peterson

# This file contains the code to generate nice tables with
# information on the metagenome sequencing and assembly.
# Also generated R objects for use in other analyses.


#### Always start with a clean slate, ya bum ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP/")
library(dplyr)
library(gridExtra)
library(readxl)


#### Generate mapping key ####
MG.list <- read.table("metadata/lists/metagenome_list.txt",
                      header = FALSE,
                      col.names = c("metagenomeID", "samplingYear"))
assembly.list <- read.table("metadata/lists/assembly_list.txt",
                            header = FALSE,
                            col.names = c("assemblyID", "samplingYear"))
mapping.key <- full_join(MG.list,
                         assembly.list) %>%
  select(metagenomeID, assemblyID, samplingYear)

write.table(mapping.key,
            "metadata/lists/mapping_key.tsv",
            col.names = FALSE,
            row.names = FALSE,
            sep = "\t",
            quote = FALSE)




#### Generate table of MG sample site information ####

# Read in filter metadata
filter.data <- read_xlsx("metadata/raw_metadata/NA_IDs.xlsx") %>%
  select(filterID, sampleID, volumeFiltered, flashFrozenOnLN2, notes) %>%
  filter(!is.na(sampleID)) %>%
  # Read in sample metadata
  left_join(read_xlsx("metadata/raw_metadata/2_sample_IDs.xlsx") %>%
              select(sampleID, tripID, depth)) %>%
  # Read in trip metadata
  left_join(read_xlsx("metadata/raw_metadata/1_trip_IDs.xlsx") %>%
              select(tripID, startDate)) %>%
  # Read in metagenome metadata
  select(filterID, sampleID, startDate, depth, volumeFiltered, flashFrozenOnLN2, notes)


#### Save out metadata for analysis ####
filter.data %>%
  write.csv("metadata/processedMetadata/filter_metadata.csv",
            row.names = FALSE)

# right_join(read_xlsx("dataEdited/dnaExtractions/DNA_extractions.xlsx") %>%
#              filter(filterID != "NA") %>%
#              select(filterID, extractionID)) %>%
