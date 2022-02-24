#### code/cleaning_scripts/incubation_metadata.R ####
# Written for BLiMMP project
# Benjamin D. Peterson



#### Prep workspace ####

rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP")
library(dplyr)
library(lubridate)
library(readxl)
library(tidyr)



#### Aggregate all metadata ####
# Read in trip IDs
trip.id <- read_xlsx("metadata/raw_metadata/1_trip_IDs.xlsx") %>%
  select(tripID, startDate)
# Read in sample IDs
sample.ids <- read_xlsx("metadata/raw_metadata/2_sample_IDs.xlsx") %>%
  select(sampleID, tripID, depth)
# Read in incubations IDs
incubation.ids <- read_xlsx("metadata/raw_metadata/4_MA_ID.xlsx") %>%
  select(incubationID, sampleID, filtered, amendment, sampleVolume, dateSpike, timeSpike)

# Read in incubation IDs for 2019 and 2020
hg.ids.2019.2020 <- read_xlsx("metadata/raw_metadata/5_MA_Hg_samples.xlsx") %>%
  filter(!is.na(bottleID)) %>%
  # In 2021, we had separate sample bottles for MeHg and HgT,
  # so we need to accomodate that here as well
  rename(bottleIdMeHg = bottleID) %>%
  mutate(bottleIdHgT = bottleIdMeHg) %>%
  select(bottleIdMeHg, bottleIdHgT, incubationID, dateKilled, timeKilled, t, notes)
hg.ids.2021 <- read_xlsx("metadata/5_MA_Hg_samples.xlsx",
                                 sheet = "2021_MeHg") %>%
  rename(bottleIdMeHg = bottleID) %>%
  select(-c(notes, volCollected)) %>%
  full_join(read_xlsx("metadata/raw_metadata/5_MA_Hg_samples.xlsx",
                      sheet = "2021_HgT")) %>%
  rename(bottleIdHgT = bottleID) %>%
  # Don't keep the bottle ID info here, we just
  # want the data for when we collected samples
  select(bottleIdMeHg, bottleIdHgT, incubationID, dateKilled, timeKilled, t, notes)
hg.ids <- rbind(hg.ids.2021,
                hg.ids.2019.2020)


#### Combine all the data ####
all.metadata <- hg.ids %>%
  left_join(incubation.ids,
            by = "incubationID") %>%
  left_join(sample.ids) %>%
  left_join(trip.id)
rm(hg.ids, hg.ids.2021, hg.ids.2019.2020,
   incubation.ids, sample.ids, trip.id)


#### Add treatment column with all treatment information ####
filtered.vector <- c("filtered", "unfiltered")
names(filtered.vector) <- c("yes", "no")
all.metadata <- all.metadata %>%
  mutate(treatment = paste(filtered.vector[filtered],
                           amendment,
                           sep = "-"))


#### Make incubation time column ####
incubationTime <- all.metadata %>%
  mutate(dateTime = ymd_hm(paste(dateKilled,
                                 timeKilled,
                                 sep = " "))) %>%
  select(incubationID, dateTime, t) %>%
  spread(key = t,
         value = dateTime) %>%
  mutate(t1 = as.duration(t0 %--% t1),
         t2 = as.duration(t0 %--% t2),
         t0 = 0) %>%
  gather(key = t,
         value = durationInDays,
         -1) %>%
  filter(!is.na(durationInDays)) %>%
  mutate(durationInDays = round(durationInDays / 60 / 60 / 24, 6))
# Add to rest of metadata
all.metadata <- all.metadata %>%
  left_join(incubationTime)
rm(incubationTime)


#### Make incubation time since collection column ####
incubationTimeSinceT <- all.metadata %>%
  mutate(dateTime = ymd_hm(paste(dateKilled,
                                 timeKilled,
                                 sep = " "))) %>%
  select(incubationID, dateTime, t) %>%
  spread(key = t,
         value = dateTime) %>%
  mutate(t2 = as.duration(t1 %--% t2),
         t1 = as.duration(t0 %--% t1),
         t0 = 0) %>%
  gather(key = t,
         value = durationSinceTimepointInDays,
         -1) %>%
  filter(!is.na(durationSinceTimepointInDays)) %>%
  mutate(durationSinceTimepointInDays = round(durationSinceTimepointInDays / 60 / 60 / 24, 6))
all.metadata <- all.metadata %>%
  left_join(incubationTimeSinceT)
rm(incubationTimeSinceT)


#### Make the data long with respect to the HgT and MeHg barcodes ####
all.metadata.long <- all.metadata %>%
  gather(key = constituent,
         value = bottleID,
         c(1, 2)) %>%
  mutate(constituent = gsub("bottleId", "", constituent))


#### Read out metadata ####
metadata.to.save.out <- all.metadata.long %>%
  select(bottleID, constituent, incubationID, sampleID, tripID, dateKilled,
         timeKilled, depth, t, durationInDays, durationSinceTimepointInDays,
         filtered, amendment, treatment) %>%
  arrange(incubationID, constituent, t)


#### Read out metadata ####
write.csv(metadata.to.save.out,
          "metadata/processedMetadata/incubation_metadata.csv",
          row.names = FALSE)
