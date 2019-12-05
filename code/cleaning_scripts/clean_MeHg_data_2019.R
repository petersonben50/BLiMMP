#### code/cleaning_scripts/clean_Hg_2019_data.R ####
# Written for BLiMMP project
# Benjamin D. Peterson


#### Prep workspace ####
rm(list = ls())
setwd("~/Box/BLiMMP/")
library(dplyr)
library(readODS)
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
  select(bottleID, incubationID, dateKilled, timeKilled, t, notes)
final.data <- hg.ids %>%
  filter(bottleID != "") %>%
  left_join(incubation.ids,
            by = "incubationID") %>%
  left_join(sample.ids) %>%
  left_join(trip.id)
data.for.checking <- final.data %>%
  select(bottleID, incubationID, tripID, dateKilled, timeKilled, depth, t, filtered, amendment)

rm(hg.ids, incubation.ids, sample.ids, trip.id)

write.csv(data.for.checking,
          "metadata/processedMetadata/incubation_Hg_metadata.csv",
          row.names = FALSE,
          quote = FALSE)



#### Calculate time gaps in processing (spiking or sample collection) ####

processing.times <- final.data %>%
  filter(tripID != "BLiMMP_trip_003") %>%
  mutate(date_time_killed = ymd_hm(paste(dateKilled,
                                         timeKilled),
                                   tz = "US/Central")) %>%
  mutate(date_time_spiked = ymd_hm(paste(dateSpike,
                                         timeSpike),
                                   tz = "US/Central")) %>%
  select(incubationID, date_time_spiked, date_time_killed, t) %>%
  spread(key = t,
         value = date_time_killed) %>%
  rename(t0 = `0`,
         t24 = `24`) %>%
  mutate(spike_to_kill_0_time = as.duration(date_time_spiked %--% t0),
         t0_to_t24_time = as.duration(t0 %--% t24)) %>%
  select(incubationID, spike_to_kill_0_time, t0_to_t24_time)


#### Methylmercury data, first run ####
file.name <- "dataRaw/incubations/MeHg/I080719 BENDOTA.xlsx"
file.data <- read_xlsx(file.name,
                         sheet = "Summary",
                         skip = 6)
# This whole run looked good, so we'll keep all of it. 
file.data.clean <- file.data %>%
  select(1, 3, 4, 5, 9) %>%
  as.data.frame()
colnames(file.data.clean) <- c("bottleID",
                                 "amb_MeHg_ng.L",
                                 "excess_MeHg_198_ng.L",
                                 "excess_MeHg_204_ng.L",
                                 "excess_DDL")
file.data.clean <- file.data.clean[1:20, ]

# Remove spikes
file.data.clean <- file.data.clean %>%
  filter(!duplicated(bottleID))

# Combine with metadata
file.data.clean <- file.data.clean %>%
  left_join(data.for.checking)

# Generate single treatment column
filter.vector <- c("unfiltered", "filtered")
names(filter.vector) <- c("no", "yes")

file.data.clean <- file.data.clean %>%
  mutate(filtered = filter.vector[filtered])

unamended.index <- file.data.clean$treatment == ""
file.data.clean$treatment[unamended.index] <- "unamended"

file.data.clean <- file.data.clean %>%
  mutate(treatment = paste(filtered,
                           treatment,
                           sep = "-")) %>%
  select(-filtered)

# Read out data
write.csv(file.data.clean,
          "dataEdited/incubations/MeHg/incubations2019_20190807.csv",
          row.names = FALSE,
          quote = FALSE)
rm(list = ls(pattern = "file.*"))


#### Methylmercury data, 2019-09-24 ####

file.name <- "dataRaw/incubations/MeHg/I092419 BENDOTA_unchecked.xlsx"
file.data <- read_xlsx(file.name,
                         sheet = "Summary",
                         skip = 6)
file.data.clean <- file.data %>%
  select(1, 3, 4, 5, 9) %>%
  as.data.frame()
colnames(file.data.clean) <- c("bottleID",
                                 "amb_MeHg_ng.L",
                                 "excess_MeHg_198_ng.L",
                                 "excess_MeHg_204_ng.L",
                                 "excess_DDL")
file.data.clean <- file.data.clean[1:36, ]

# Remove spikes
file.data.clean <- file.data.clean %>%
  filter(!duplicated(bottleID))

# Combine with metadata
file.data.clean <- file.data.clean %>%
  left_join(data.for.checking)

# Generate single treatment column
filter.vector <- c("unfiltered", "filtered")
names(filter.vector) <- c("no", "yes")

file.data.clean <- file.data.clean %>%
  mutate(filtered = filter.vector[filtered])

unamended.index <- file.data.clean$treatment == ""
file.data.clean$treatment[unamended.index] <- "unamended"

file.data.clean <- file.data.clean %>%
  mutate(treatment = paste(filtered,
                           treatment,
                           sep = "-")) %>%
  select(-filtered)

write.csv(file.data.clean,
          "dataEdited/incubations/MeHg/incubations2019_20190924.csv",
          row.names = FALSE,
          quote = FALSE)
rm(list = ls(pattern = "file.*"))


#### Methylmercury data, 2019-09-25 ####

file.name <- "dataRaw/incubations/MeHg/I092519 BENDOTA, HCC ISCO_unchecked.xlsx"
file.data <- read_xlsx(file.name,
                       sheet = "Summary",
                       skip = 6)
file.data.clean <- file.data %>%
  select(1, 3, 4, 5, 9) %>%
  as.data.frame()
colnames(file.data.clean) <- c("bottleID",
                               "amb_MeHg_ng.L",
                               "excess_MeHg_198_ng.L",
                               "excess_MeHg_204_ng.L",
                               "excess_DDL")
file.data.clean <- file.data.clean[1:36, ]

# Remove spikes
file.data.clean <- file.data.clean %>%
  filter(!duplicated(bottleID))

# Combine with metadata
file.data.clean <- file.data.clean %>%
  left_join(data.for.checking) %>%
  filter(!is.na(incubationID))

# Generate single treatment column
filter.vector <- c("unfiltered", "filtered")
names(filter.vector) <- c("no", "yes")

file.data.clean <- file.data.clean %>%
  mutate(filtered = filter.vector[filtered])

unamended.index <- file.data.clean$treatment == ""
file.data.clean$treatment[unamended.index] <- "unamended"

file.data.clean <- file.data.clean %>%
  mutate(treatment = paste(filtered,
                           treatment,
                           sep = "-")) %>%
  select(-filtered)

write.csv(file.data.clean,
          "dataEdited/incubations/MeHg/incubations2019_20190925.csv",
          row.names = FALSE,
          quote = FALSE)
rm(list = ls(pattern = "file.*"))









#### THg data, first run ####
THg.1 <- "dataRaw/incubations/THg/I060519 BENDOTA WATERS.xlsx"
THg.1.data <- read_xlsx(THg.1,
                        sheet = "Summary",
                        skip = 10) %>%
  select(1, 3, 4, 5, 6, 10)
colnames(THg.1.data) <- c("sampleType",
                          "barcode",
                          "amb_THg_ng.L",
                          "excess_THg_198_ng.L",
                          "excess_THg_204_ng.L",
                          "excess_DDL")
THg.1.data.clean <- THg.1.data %>%
  filter(sampleType == "SAM") %>%
  select(-sampleType) %>%
  as.data.frame()

write.csv(THg.1.data.clean,
          "dataEdited/incubations/THg/incubations2018_I060519.csv",
          row.names = FALSE,
          quote = FALSE)

rm(list = ls(pattern = "THg.1*"))

#### THg data, second run ####
THg.2 <- "dataRaw/incubations/THg/I061019 BENDOTA WATERS.xlsx"
THg.2.data <- read_xlsx(THg.2,
                        sheet = "Summary",
                        skip = 10) %>%
  select(1, 3, 4, 5, 6, 10)
colnames(THg.2.data) <- c("sampleType",
                          "barcode",
                          "amb_THg_ng.L",
                          "excess_THg_198_ng.L",
                          "excess_THg_204_ng.L",
                          "excess_DDL")
THg.2.data.clean <- THg.2.data %>%
  filter(sampleType == "SAM") %>%
  select(-sampleType) %>%
  as.data.frame()

write.csv(THg.2.data.clean,
          "dataEdited/incubations/THg/incubations2018_I061019.csv",
          row.names = FALSE,
          quote = FALSE)

rm(list = ls(pattern = "THg.2*"))

#### Read in all MeHg data ####
list.o.files <- list.files(path = "dataEdited/incubations/MeHg",
                           pattern = "incubations2018")
MeHg.data.list <- list()
for (file.number in 1:length(list.o.files)) {
  MeHg.data.list[[file.number]] <- read.csv(paste("dataEdited/incubations/MeHg/",
                                                  list.o.files[file.number],
                                                  sep = ""),
                                            stringsAsFactors = FALSE)
}
MeHg.data <- do.call(rbind,
                     MeHg.data.list)
rm(list.o.files,
   MeHg.data.list,
   file.number)

#### Combine MeHg data with metadata ####
MeHg.data.final <- left_join(MeHg.data,
                             metadata) %>%
  arrange(timePoint) %>%
  arrange(bag)


#### Read in all THg data ####
list.o.files <- list.files(path = "dataEdited/incubations/THg/",
                           pattern = "incubations2018")
THg.data.list <- list()
for (file.number in 1:length(list.o.files)) {
  THg.data.list[[file.number]] <- read.csv(paste("dataEdited/incubations/THg/",
                                                  list.o.files[file.number],
                                                  sep = ""),
                                            stringsAsFactors = FALSE)
}
THg.data <- do.call(rbind,
                     THg.data.list)
rm(list.o.files,
   THg.data.list,
   file.number)

#### Combine THg data with metadata ####
THg.data.final <- left_join(THg.data,
                             metadata) %>%
  arrange(timePoint) %>%
  arrange(bag)

rm(THg.data)
