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
all.metadata <- hg.ids %>%
  filter(bottleID != "") %>%
  left_join(incubation.ids,
            by = "incubationID") %>%
  left_join(sample.ids) %>%
  left_join(trip.id) 

metadata.to.save.out <- all.metadata %>%
  select(bottleID, incubationID, sampleID, tripID, dateKilled,
         timeKilled, depth, t, filtered, amendment)

write.csv(metadata.to.save.out,
          "metadata/processedMetadata/incubation_Hg_metadata.csv",
          row.names = FALSE,
          quote = FALSE)

rm(hg.ids, incubation.ids, sample.ids, trip.id, metadata.to.save.out)



#### Calculate time gaps in processing (spiking or sample collection) ####

processing.data <- all.metadata %>%
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
  mutate(spike_to_kill_0_time = as.duration(date_time_spiked %--% t0),
         t0_to_t1_time = as.duration(t0 %--% t1)) %>%
  select(incubationID, spike_to_kill_0_time, t0_to_t1_time)



#### Cleaning function ####

# file.name.input <- "dataRaw/incubations/MeHg/I080719 BENDOTA.xlsx"
# output.file.name <- "dataEdited/incubations/MeHg/incubations2019_20190807.csv"

clean.MeHg.data.from <- function(file.name.input,
                                 output.file.name,
                                 metadata.df = all.metadata,
                                 barcodes.to.remove = NULL) {
  
  file.name <- file.name.input
  file.data <- read_xlsx(file.name,
                         sheet = "Summary",
                         skip = 6)
  
  # Select needed columns and rename them.
  file.data.clean <- file.data %>%
    select(1, 3, 4, 5, 9) %>%
    as.data.frame()
  colnames(file.data.clean) <- c("bottleID",
                                 "amb_MeHg_ng.L",
                                 "excess_MeHg_198_ng.L",
                                 "excess_MeHg_204_ng.L",
                                 "excess_DDL")
  file.data.clean <- file.data.clean %>%
    filter(bottleID %in% metadata.df$bottleID)
  
  # Remove unwanted samples by barcode
  if (!is.null(barcodes.to.remove)) {
    
    file.data.clean <- file.data.clean %>%
      filter(!(bottleID %in% barcodes.to.remove))
    
  }
  
  # Remove spikes
  file.data.clean <- file.data.clean %>%
    filter(!duplicated(bottleID)) 
  
  # Round off numbers
  file.data.clean[, -1] <- file.data.clean[, -1] %>%
    mutate_all(as.numeric) %>% 
    round(digits = 3)
  
  # Write out data
  write.csv(file.data.clean,
            output.file.name,
            row.names = FALSE,
            quote = FALSE)
}



#### Clean files ####

clean.MeHg.data.from("dataRaw/incubations/MeHg/I080719 BENDOTA.xlsx",
                     "dataEdited/incubations/MeHg/cleaned/incubations2019_20190807.csv")

clean.MeHg.data.from("dataRaw/incubations/MeHg/I092419 BENDOTA_unchecked.xlsx",
                     "dataEdited/incubations/MeHg/cleaned/incubations2019_20190924.csv")

clean.MeHg.data.from("dataRaw/incubations/MeHg/I092519 BENDOTA, HCC ISCO_unchecked.xlsx",
                     "dataEdited/incubations/MeHg/cleaned/incubations2019_20190925.csv")

clean.MeHg.data.from("dataRaw/incubations/MeHg/I120419 BENDOTA.xlsx",
                     "dataEdited/incubations/MeHg/cleaned/incubations2019_20191204.csv",
                     barcodes.to.remove = "MSC576AR")

clean.MeHg.data.from("dataRaw/incubations/MeHg/I120519 BENDOTA.xlsx",
                     "dataEdited/incubations/MeHg/cleaned/incubations2019_20191205.csv")



#### Combine MeHg data ####

rm(clean.MeHg.data.from)

# Load up all incubation data
list.o.results <- list.files(path = "dataEdited/incubations/MeHg/cleaned",
                             pattern = "incubations2019",
                             full.names = TRUE)

for (file.name in list.o.results) {
  if (file.name == list.o.results[1]) {
    MeHg.results <- read.csv(file.name,
                             stringsAsFactors = FALSE)
  } else {
    MeHg.results <- rbind(MeHg.results,
                          read.csv(file.name,
                                   stringsAsFactors = FALSE))
  }
  
}

rm(list.o.results,
   file.name)



#### Normalize t1 values to 24 hours ####

MeHg.198.production.t1 <- MeHg.results %>%
  select(bottleID, excess_MeHg_198_ng.L) %>%
  left_join(all.metadata) %>%
  filter(t == "t1") %>%
  full_join(processing.data) %>%
  mutate(t0_to_t1_time = (as.numeric(t0_to_t1_time) / (24*60*60))) %>%
  arrange(incubationID) %>%
  mutate(t0_to_t1_time = replace_na(t0_to_t1_time, 1)) %>%
  mutate(excess_MeHg_198_ng.L_normalized = round((excess_MeHg_198_ng.L / t0_to_t1_time), 3)) %>%
  select(bottleID, incubationID, sampleID, tripID, startDate, depth, t, excess_MeHg_198_ng.L_normalized) %>%
  rename(excess_MeHg_198_ng.L = excess_MeHg_198_ng.L_normalized)


#### Check out spike-to-kill times and if it influences t0 values ####

check_t0_timing <- MeHg.results  %>%
  left_join(all.metadata) %>%
  filter(t == "t0")%>%
  full_join(processing.data) %>%
  filter(!is.na(spike_to_kill_0_time)) %>%
  select(excess_MeHg_198_ng.L, spike_to_kill_0_time)
plot(x = check_t0_timing$spike_to_kill_0_time,
     y = check_t0_timing$excess_MeHg_198_ng.L,
     pch = 18,
     xlab = "Spike time to kill time (sec)",
     ylab = "Excess Me198Hg in t0 sample")
# The high MeHg samples in the t0 samples do not seem to be 
# related to longer wait times between spiking and killing t=0.
# More likely to be due to contamination of some sort. 
# We'll keep the t0 time points as they are. 

MeHg.198.t0 <- MeHg.results  %>%
  left_join(all.metadata) %>%
  filter(t == "t0")%>%
  full_join(processing.data) %>%
  select(bottleID, incubationID, sampleID, tripID, startDate, depth, t, excess_MeHg_198_ng.L)



#### Combine Me198Hg data and save it out ####

MeHg.198.data <- rbind(MeHg.198.t0, MeHg.198.production.t1)

write.csv(MeHg.198.data,
          "dataEdited/incubations/MeHg/incubations2019_Me198Hg.csv",
          row.names = FALSE)
