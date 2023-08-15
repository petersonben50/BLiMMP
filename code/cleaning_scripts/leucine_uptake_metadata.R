#### code/cleaning_scripts/clean_leucine_uptake_2020.R ####
# Benjamin D. Peterson

#### Clean up ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP")
library(readxl)
library(tidyverse)


#### Read in metadata for incubations ####

# Incubation data
leucine_incubations_2020 <- read_xlsx("metadata/raw_metadata/LEU.xlsx",
                                 sheet = "LEU_2019_2020") %>%
  select(uptakeID, sampleID) %>%
  filter(grepl("BLI20", uptakeID)) %>%
  mutate(protocol = "direct",
         # The "direct" protocol means the sample in the syringe was collected directly from the water column, not incubated in a serum bottle first.
         treatment = NA,
         timePoint = NA) %>%
  select(sampleID, uptakeID, protocol, treatment, timePoint)
leucine_syringes_2021 <- read_xlsx("metadata/raw_metadata/LEU.xlsx",
                                   sheet = "Leu_syringes_2021") %>%
  select(-notes) %>%
  filter(!is.na(leu_inc_ID))
leucine_incubations_2021 <- read_xlsx("metadata/raw_metadata/LEU.xlsx",
                                 sheet = "Leu_incubations_2021") %>%
  select(-notes) %>%
  filter(!is.na(sampleID)) %>%
  full_join(leucine_syringes_2021,
            by = "leu_inc_ID",
            multiple = "all") %>%
  mutate(protocol = "indirect") %>% # The "indirect" protocol means the sample was collected into a serum bottle and incubated in the lab before testing leucine uptake rates.
  select(sampleID, uptakeID, protocol, treatment, timePoint)

rm(leucine_syringes_2021)
leucine_incubations <- rbind(leucine_incubations_2020,
                             leucine_incubations_2021)
rm(leucine_incubations_2020,
   leucine_incubations_2021)

# Add sample metadata
sampleID_metadata <- read_xlsx("metadata/raw_metadata/2_sample_IDs.xlsx") %>%
  select(-notes)
tripID_metadata <- read_xlsx("metadata/raw_metadata/1_trip_IDs.xlsx") %>%
  select(-notes)
leu_metadata <- leucine_incubations %>%
  left_join(sampleID_metadata) %>%
  left_join(tripID_metadata) %>%
  select(uptakeID, sampleID, depth, startDate, protocol, treatment, timePoint) %>%
  filter(startDate %in% c("2020-09-02", "2020-10-10",
                          "2021-09-10", "2021-10-14"))


# Save out metadata
write.csv(leu_metadata,
          "metadata/processedMetadata/LEU_metadata.csv",
          row.names = FALSE)
