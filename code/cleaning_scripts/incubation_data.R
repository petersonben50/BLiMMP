#### code/cleaning_scripts/incubation_data.R ####
# Written for BLiMMP project
# Benjamin D. Peterson


#### Prep workspace ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP")
library(lubridate)
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
  mutate(bottleID_MeHg = bottleID) %>%
  select(-c(bottleID, constituent))
rm(MeHg.data)

#### Prepare HgT data ####
HgT.data.clean <- right_join(incubation.metadata %>%
                               filter(constituent == "HgT"),
                             HgT.data) %>%
  mutate(bottleID_HgT = bottleID) %>%
  select(-c(bottleID, constituent))
rm(HgT.data, incubation.metadata)


#### Combine data ####
Hg.incubation.data <- full_join(MeHg.data.clean,
                                HgT.data.clean)
rm(MeHg.data.clean,
   HgT.data.clean)


#### Calculate percent MeHg for ambient, 198, and 204
Hg.incubation.data <- Hg.incubation.data %>%
  mutate(percent_amb_MeHg = MeHg_ambient_ppt / HgT_ambient_ppt * 100,
         percent_198_MeHg = MeHg_198_ppt / HgT_198_ppt * 100,
         percent_204_MeHg = MeHg_204_ppt / HgT_204_ppt * 100)

#### Read out all data ####
write.csv(Hg.incubation.data,
          file = "dataEdited/incubations/clean_incubation_data.csv",
          row.names = FALSE)



#### Calculate Kmet from t1 ####
kmet.data <- Hg.incubation.data %>%
  select(sampleID, incubationID,
         startDate, depth,
         t, treatment, durationInDays,
         MeHg_198_ppt, HgT_198_ppt) %>%
  filter(treatment %in% c("unfiltered-unamended",
                          "unfiltered-molybdate",
                          "filtered-unamended"),
         t == "t1") %>%
  mutate(Kmet = (MeHg_198_ppt / HgT_198_ppt) / durationInDays)


#### Read out Kmet data ####
write.csv(kmet.data,
          file = "dataEdited/incubations/clean_incubation_Kmet_data.csv",
          row.names = FALSE)
