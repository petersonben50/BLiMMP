#### code/cleaning_scripts/clean_THg_2019_data.R ####
# Written for BLiMMP project
# Benjamin D. Peterson




#### Every last drop ####

rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP")
#library(readODS)
library(lubridate)
library(readxl)
library(tidyverse)




#### Generate needed metadata ####

trip.id <- read_xlsx("metadata/1_trip_IDs.xlsx") %>%
  select(tripID, startDate)
sample.ids <- read_xlsx("metadata/2_sample_IDs.xlsx") %>%
  select(sampleID, tripID, depth)
incubation.ids <- read_xlsx("metadata/4_MA_ID.xlsx") %>%
  select(incubationID, sampleID, filtered, amendment, sampleVolume, dateSpike, timeSpike)
hg.ids <- read_xlsx("metadata/5_MA_Hg_samples.xlsx") %>%
  select(bottleID, incubationID, dateKilled, timeKilled, t, notes)
all.metadata <- hg.ids %>%
  filter(bottleID != "") %>%
  left_join(incubation.ids,
            by = "incubationID") %>%
  left_join(sample.ids) %>%
  left_join(trip.id)

# Add treatment column with all treatment information
filtered.vector <- c("filtered", "unfiltered")
names(filtered.vector) <- c("yes", "no")
all.metadata <- all.metadata %>%
  mutate(treatment = paste(filtered.vector[filtered],
                           amendment,
                           sep = "-"))

rm(hg.ids, incubation.ids, sample.ids, trip.id, filtered.vector)





#### THg data, first run ####

# Jake ran these samples for me, and aggregated all
# the data into one file
THg.data <- read_xlsx("dataRaw/incubations/THg/BENDOTA SUMMARY.xlsx") %>%
  rename(bottleID = `Bottle ID`,
         amb_Hg = `Amb THg (from 202) (ng/L)`,
         THg_198 = `Excess 198THg (ng/L)`,
         THg_204 = `Excess 204THg (ng/L)`,
         DDL_excess = `Excess DDL (ng/L)`,
         DDL_amb = `Amb DDL (ng/L)`) %>%
  select(bottleID, amb_Hg, THg_198,
         THg_204, DDL_excess, DDL_amb)






#### Process excess THg data ####

THg.data.final <- THg.data %>%
  left_join(all.metadata) %>%
  select(bottleID, incubationID, sampleID, tripID, startDate, depth, treatment, t, amb_Hg, THg_198, THg_204, DDL_excess, DDL_amb)
rm(THg.data, all.metadata)






#### Write out data ####

write.csv(THg.data.final,
          "dataEdited/incubations/THg/incubations2019_THg.csv",
          row.names = FALSE)
