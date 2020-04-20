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




#### Generate processing metadata ####

processing.data <- all.metadata %>%
  # filter(tripID != "BLiMMP_trip_003") %>%
  mutate(date_time_killed = ymd_hm(paste(dateKilled,
                                         timeKilled),
                                   tz = "US/Central")) %>%
  mutate(date_time_spiked = ymd_hm(paste(dateSpike,
                                         timeSpike),
                                   tz = "US/Central")) %>%
  select(incubationID, date_time_spiked, date_time_killed, t) %>%
  spread(key = t,
         value = date_time_killed) %>%
  mutate(spike_to_kill_0_time = as.duration(date_time_spiked %--% t0),
         t0_to_t1_time = as.duration(t0 %--% t1)) %>%
  select(incubationID, spike_to_kill_0_time, t0_to_t1_time)





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
         THg_204, DDL_excess, DDL_amb) %>%
  mutate(amb_Hg = round(amb_Hg, 3),
         THg_198 = round(THg_198, 3),
         THg_204 = round(THg_204, 3))






#### Process excess THg data ####

THg.data.final <- THg.data %>%
  left_join(all.metadata) %>%
  select(bottleID, incubationID, sampleID, tripID, startDate, depth, treatment, t, amb_Hg, THg_198, THg_204, DDL_excess, DDL_amb)
rm(THg.data, all.metadata)






#### Write out data ####

write.csv(THg.data.final,
          "dataEdited/incubations/THg/incubations2019_THg.csv",
          row.names = FALSE)




#### Calculate change in 204Hg in each bag ####

# Calculate change
THg.data.final.delta.204 <- THg.data.final %>%
  select(-c(THg_198, DDL_excess, DDL_amb, amb_Hg, bottleID)) %>%
  spread(key = t,
         value = THg_204) %>%
  mutate(per_change_204 = ((t1 - t0) / t0)*100) 
write.csv(THg.data.final.delta.204,
          "dataEdited/incubations/THg/incubations2019_204Hg_delta.csv",
          row.names = FALSE)

# Calculate normalized change
normalization_time <- 60 * 60 * 24
THg.data.final.delta.204.norm <- THg.data.final.delta.204 %>%
  full_join(processing.data) %>%
  mutate(per_change_204_norm = round(((per_change_204 * as.numeric(t0_to_t1_time)) / normalization_time),
                                     3)) %>%
  select(-c(spike_to_kill_0_time,
            t0_to_t1_time,
            per_change_204))
write.csv(THg.data.final.delta.204.norm,
          "dataEdited/incubations/THg/incubations2019_204Hg_delta_norm.csv",
          row.names = FALSE)





#### Calculate change in 198Hg in each bag ####

# Calculate change
THg.data.final.delta.198 <- THg.data.final %>%
  select(-c(THg_204, DDL_excess, DDL_amb, amb_Hg, bottleID)) %>%
  spread(key = t,
         value = THg_198) %>%
  mutate(per_change_198 = ((t1 - t0) / t0)*100) 
write.csv(THg.data.final.delta.198,
          "dataEdited/incubations/THg/incubations2019_198Hg_delta.csv",
          row.names = FALSE)

# Calculate normalized change
normalization_time <- 60 * 60 * 24
THg.data.final.delta.198.norm <- THg.data.final.delta.198 %>%
  full_join(processing.data) %>%
  mutate(per_change_198_norm = round(((per_change_198 * as.numeric(t0_to_t1_time)) / normalization_time),
                                     3)) %>%
  select(-c(spike_to_kill_0_time,
            t0_to_t1_time,
            per_change_198))
write.csv(THg.data.final.delta.198.norm,
          "dataEdited/incubations/THg/incubations2019_198Hg_delta_norm.csv",
          row.names = FALSE)



