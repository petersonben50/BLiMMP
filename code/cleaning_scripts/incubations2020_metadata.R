#### code/cleaning_scripts/clean_incubation_metadata_2020.R ####
# Written for BLiMMP project
# Benjamin D. Peterson



#### Prep workspace ####

rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP")
library(dplyr)
library(lubridate)
library(readxl)
library(tidyr)



#### Generate correct metadata ####

trip.id <- read_xlsx("metadata/1_trip_IDs.xlsx") %>%
  select(tripID, startDate)
sample.ids <- read_xlsx("metadata/2_sample_IDs.xlsx") %>%
  select(sampleID, tripID, depth)
incubation.ids <- read_xlsx("metadata/4_MA_ID.xlsx") %>%
  select(incubationID, sampleID, filtered, amendment, sampleVolume, dateSpike, timeSpike)
hg.ids <- read_xlsx("metadata/5_MA_Hg_samples.xlsx") %>%
  select(bottleID, incubationID, dateKilled, timeKilled, volCollected, t, notes)
all.metadata.2020 <- hg.ids %>%
  filter(bottleID != "") %>%
  filter(year(dateKilled) == 2020) %>%
  left_join(incubation.ids,
            by = "incubationID") %>%
  left_join(sample.ids) %>%
  left_join(trip.id)



# Add treatment column with all treatment information
filtered.vector <- c("filtered", "unfiltered")
names(filtered.vector) <- c("yes", "no")
all.metadata.2020 <- all.metadata.2020 %>%
  mutate(treatment = paste(filtered.vector[filtered],
                           amendment,
                           sep = "-"))


#### Make incubation time column ####

incubationTime <- all.metadata.2020 %>%
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

all.metadata.2020 <- all.metadata.2020 %>%
  left_join(incubationTime)


#### Make incubation time since collection column ####

incubationTimeSinceT <- all.metadata.2020 %>%
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

all.metadata.2020 <- all.metadata.2020 %>%
  left_join(incubationTimeSinceT)




#### Read out metadata ####
metadata.to.save.out <- all.metadata.2020 %>%
  select(bottleID, incubationID, sampleID, tripID, dateKilled,
         timeKilled, volCollected, depth, t, durationInDays,
         durationSinceTimepointInDays,filtered, amendment, treatment) %>%
  arrange(bottleID)

write.csv(metadata.to.save.out,
          "metadata/processedMetadata/incubation_Hg_metadata_2020.csv",
          row.names = FALSE,
          quote = FALSE)

rm(hg.ids, incubation.ids, incubationTime, sample.ids, trip.id, filtered.vector, all.metadata.2020)



#### Read out metadata for USGS team ####
USGS.metadata <- metadata.to.save.out %>%
  mutate(`Site Name` = "MENDOTA - DEEP HOLE") %>%
  mutate(`Site ID` = "") %>%
  mutate(Date = dateKilled) %>%
  mutate(Time = timeKilled) %>%
  mutate(Depth = depth) %>%
  mutate(Length = "-999") %>%
  mutate(Rep = "") %>%
  mutate(`Sample Comments` = "") %>%
  mutate(`Container ID` = bottleID) %>%
  mutate(Medium = "WSQ") %>%
  mutate(Analysis = "") %>%
  mutate(Isotope = "") %>%
  mutate(Filter = "") %>%
  mutate(`Filter Vol` = "") %>%
  mutate(Preservation = "ACID") %>%
  mutate(Acid = "HCL") %>%
  mutate(`Acid Vol` = as.numeric(volCollected) * 0.02) %>%
  mutate(Comments = paste("incubationID=", incubationID, ";",
                          "sampleID=", sampleID, ";",
                          "filtered=", filtered, ";",
                          "amendment=", amendment, ";",
                          "t=", t, ";",
                          "durationInDays=", durationInDays, ";",
                          "volCollected=", volCollected,
                          sep = "")) %>%
  select(`Site Name`, `Site ID`, Date, Time, Depth, Length, Rep,
         `Sample Comments`, `Container ID`, Medium, Analysis,
         Isotope, Filter, `Filter Vol`, Preservation, Acid,
         `Acid Vol`, Comments)

# The medium values need to be WSQ for all of them except the first one for each date, depth, site.
USGS.metadata.uniq.date.depth <- USGS.metadata %>%
  mutate(date.depth = paste(Date, Depth)) %>%
  arrange(Time) %>%
  filter(!duplicated(date.depth,)) %>%
  select(`Container ID`) %>%
  unlist(use.names = FALSE)
USGS.metadata[which(USGS.metadata$`Container ID` %in% USGS.metadata.uniq.date.depth), "Medium"] <- "WS"


write.csv(USGS.metadata,
          file = "metadata/processedMetadata/incubation_Hg_metadata_2020_USGS.csv",
          row.names = FALSE)
