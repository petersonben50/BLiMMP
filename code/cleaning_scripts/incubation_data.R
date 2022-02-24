#### code/cleaning_scripts/incubation_data.R ####
# Written for BLiMMP project
# Benjamin D. Peterson


#### Prep workspace ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP")
library(readxl)
library(tidyverse)



#### Read in data ####
incubation.metadata <- read.csv("metadata/processedMetadata/incubation_metadata.csv")
MeHg.data <- read_xlsx("dataRaw/incubations/MeHg/MeHg_BENDOTA_SUMMARY.xlsx",
                       col_names = c("bottleID", "MeHg_ambient_ppt", "MeHg_198_ppt", "MeHg_204_ppt", "MeHg_excess_DDL", "MeHg_DOA"),
                       skip = 1)
HgT.data <- read_xlsx("dataRaw/incubations/HgT/HgT_BENDOTA_SUMMARY.xlsx",
                      col_names = c("bottleID", "HgT_ambient_ppt", "HgT_198_ppt", "HgT_204_ppt", "HgT_excess_DDL", "HgT_DOA"),
                      skip = 1)


#### Prepare MeHg data ####
MeHg.data.clean <- right_join(incubation.metadata %>%
                                filter(constituent == "MeHg"),
                              MeHg.data) %>%
  # Identify which samples measured Me198Hg above detection
  mutate(MeHg_198_above_DDL = (MeHg_198_ppt >= MeHg_excess_DDL)) %>%
  # Round off measurements to 3 decimal places
  mutate(MeHg_ambient_ppt = round(MeHg_ambient_ppt, 3),
         MeHg_198_ppt = round(MeHg_198_ppt, 3),
         MeHg_204_ppt = round(MeHg_204_ppt, 3)) %>%
  select(-c(bottleID, constituent, MeHg_DOA))
rm(MeHg.data)

#### Prepare HgT data ####
HgT.data.clean <- right_join(incubation.metadata %>%
                               filter(constituent == "HgT"),
                             HgT.data) %>%
  select(-c(bottleID, constituent, HgT_DOA))
rm(HgT.data, incubation.metadata)


#### Combine data ####
Hg.incubation.data <- full_join(MeHg.data.clean,
                                HgT.data.clean)
rm(MeHg.data.clean,
   HgT.data.clean)


#### Read out test data for 2022 ####
write.csv(Hg.incubation.data %>% filter(year(startDate) == 2021), file = "~/Downloads/Hg_inc_data.csv")
