#### code/cleaning_scripts/geochem_metadata.R ####
# Benjamin Peterson

# This script will concatenate the needed metadata
# for the geochemistry samples



#### Get ready, and fly! ####

rm(list = ls())
setwd("~/Documents/research/BLiMMP/")
library(dplyr)
library(lubridate)
library(tidyr)
library(readxl)



#### Prepare metadata ####
S.metadata <- read_xlsx("metadata/chem_S.xlsx") %>%
  select(-notes)
sample_IDs <- read_xlsx("metadata/2_sample_IDs.xlsx")
incubation_IDs <- read_xlsx("metadata/4_MA_ID.xlsx") %>%
  select(-notes)
trip_IDs <- read_xlsx("metadata/1_trip_IDs.xlsx") %>%
  select(-notes)



#### Save out water column data ####
WC.metadata <- S.metadata %>%
  filter(!(sampleID == "NA")) %>%
  left_join(sample_IDs) %>%
  left_join(trip_IDs) %>%
  select(sulfurID, sampleID, tripID, depth, startDate)
write.csv(WC.metadata,
          file = "metadata/processedMetadata/sulfide_WC.csv",
          row.names = FALSE,
          quote = FALSE)



#### Save out incubation data ####
MA.metadata <- S.metadata %>%
  filter(!(incubationID == "NA")) %>%
  select(sulfurID, incubationID, incubationTimePoint) %>%
  left_join(incubation_IDs) %>%
  left_join(sample_IDs) %>%
  left_join(trip_IDs) %>%
  select(sulfurID, incubationID, sampleID, tripID, depth, startDate, dateCollected, filtered, amendment, incubationTimePoint)

# Add treatment column with all treatment information
filtered.vector <- c("filtered", "unfiltered")
names(filtered.vector) <- c("yes", "no")
MA.metadata <- MA.metadata %>%
  mutate(treatment = paste(filtered.vector[filtered],
                           amendment,
                           sep = "-")) %>%
  as.data.frame()
write.csv(MA.metadata,
          file = "metadata/processedMetadata/sulfide_MA.csv",
          row.names = FALSE,
          quote = FALSE)

rm(list = ls())




#### Prepare metadata for 2020 metals ####
filtered.metals.metadata <- read_xlsx("metadata/chem_FM.xlsx") %>%
  select(-notes) %>%
  mutate(filteredInField = "yes") %>%
  rename(metalID = MFID) %>%
  filter(grepl("BLI20",
               metalID))
unfiltered.metals.metadata <- read_xlsx("metadata/chem_UM.xlsx") %>%
  select(-notes) %>%
  mutate(filteredInField = "no") %>%
  rename(metalID = MUID) %>%
  filter(grepl("BLI20",
               metalID))
metals.metadata <- rbind(unfiltered.metals.metadata,
                         filtered.metals.metadata)
rm(unfiltered.metals.metadata,
   filtered.metals.metadata)

sample_IDs <- read_xlsx("metadata/2_sample_IDs.xlsx")
incubation_IDs <- read_xlsx("metadata/4_MA_ID.xlsx") %>%
  select(-notes)
trip_IDs <- read_xlsx("metadata/1_trip_IDs.xlsx") %>%
  select(-notes)



#### Save out water column data for 2020 ####
WC.metals.metadata <- metals.metadata %>%
  filter(!(sampleID == "NA")) %>%
  left_join(sample_IDs) %>%
  left_join(trip_IDs) %>%
  select(metalID, sampleID, tripID, startDate, depth, tare,
         mass, preservativeID, preservativeVol, filteredInField)
write.csv(WC.metals.metadata,
          file = "metadata/processedMetadata/metals_WC.csv",
          row.names = FALSE,
          quote = FALSE)



#### Save out incubation data for 2020 ####
MA.metadata <- metals.metadata %>%
  filter(!(incubationID == "NA")) %>%
  select(metalID, incubationID, incubationTimePoint, tare,
         mass, preservativeID, preservativeVol, filteredInField) %>%
  left_join(incubation_IDs %>%
              select(incubationID, sampleID, filtered, amendment, dateCollected)) %>%
  left_join(sample_IDs) %>%
  left_join(trip_IDs) %>%
  select(metalID, incubationID, sampleID, tripID, depth, startDate,
         dateCollected, filtered, amendment, incubationTimePoint, tare,
         mass, preservativeID, preservativeVol, filteredInField)

# Add treatment column with all treatment information
filtered.vector <- c("filtered", "unfiltered")
names(filtered.vector) <- c("yes", "no")
MA.metadata <- MA.metadata %>%
  mutate(treatment = paste(filtered.vector[filtered],
                           amendment,
                           sep = "-")) %>%
  as.data.frame()
write.csv(MA.metadata,
          file = "metadata/processedMetadata/metals_MA.csv",
          row.names = FALSE,
          quote = FALSE)

