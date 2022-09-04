#### code/metagenomes/metagenome_metadata.R ####
# Benjamin D. Peterson

# This file contains the code to generate the combined
# metadata for the metagenome sequencing and assembly.


#### Always start with a clean slate, ya bum ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP/")
library(dplyr)
library(readxl)


#### Generate table of MG sample site information ####
# Read in filter metadata
all.data <- read_xlsx("metadata/raw_metadata/NA_IDs.xlsx") %>%
  select(filterID, sampleID, volumeFiltered) %>%
  # Read in extraction metadata
  right_join(read_xlsx("dataEdited/dnaExtractions/DNA_extractions.xlsx") %>%
               filter(filterID != "NA") %>%
               select(filterID, extractionID)) %>%
  # Read in sample metadata
  left_join(read_xlsx("metadata/raw_metadata/2_sample_IDs.xlsx") %>%
              select(sampleID, tripID, depth)) %>%
  # Read in trip metadata
  left_join(read_xlsx("metadata/raw_metadata/1_trip_IDs.xlsx") %>%
              select(tripID, startDate)) %>%
  # Read in metagenome metadata
  right_join(rbind(read_xlsx("dataEdited/dnaSequencing/2020/samplePrep/KMBP010_dilutions.xlsx") %>%
                     select(extractionID, metagenomeID),
                   read_xlsx("dataEdited/dnaSequencing/2021/samplePrep/BLI21_MG_dilutions.xlsx") %>%
                     rename(metagenomeID = `original metagenomeID`) %>%
                     select(extractionID, metagenomeID))) %>%
  select(metagenomeID, sampleID, startDate, depth, volumeFiltered) %>%
  as.data.frame()


#### Save out metadata for analysis ####
all.data %>%
  write.csv("metadata/metagenome_metadata.csv",
            row.names = FALSE)
