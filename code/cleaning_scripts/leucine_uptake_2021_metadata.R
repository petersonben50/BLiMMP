#### code/cleaning_scripts/clean_leucine_uptake_2020.R ####
# Benjamin D. Peterson

#### Clean up ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP")
library(readxl)
library(tidyverse)


#### Read in metadata for incubations ####

# Water column data
leucine.incubations <- read_xlsx("metadata/LEU.xlsx",
                                 sheet = "Leu_incubations_2021") %>%
  select(-notes) %>%
  filter(!is.na(sampleID))
leucine.syringes <- read_xlsx("metadata/LEU.xlsx",
                              sheet = "Leu_syringes_2021") %>%
  select(-notes) %>%
  filter(!is.na(leu_inc_ID))
sampleID.metadata <- read_xlsx("metadata/2_sample_IDs.xlsx") %>%
  select(-notes)
tripID.metadata <- read_xlsx("metadata/1_trip_IDs.xlsx") %>%
  select(-notes)

# Clean up water column data
metadata.WC <- leucine.incubations %>%
  left_join(leucine.syringes) %>%
  left_join(sampleID.metadata) %>%
  left_join(tripID.metadata) %>%
  select(uptakeID, sampleID, depth, startDate, treatment, timePoint)
write.csv(metadata.WC,
          "metadata/processedMetadata/LEU_2021.csv",
          row.names = FALSE)
