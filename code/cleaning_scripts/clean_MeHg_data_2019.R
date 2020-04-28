#### code/cleaning_scripts/clean_MeHg_2019_data.R ####
# Written for BLiMMP project
# Benjamin D. Peterson



#### Prep workspace ####

rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP")
library(dplyr)
library(lubridate)
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

# Add treatment column with all treatment information
filtered.vector <- c("filtered", "unfiltered")
names(filtered.vector) <- c("yes", "no")
all.metadata <- all.metadata %>%
  mutate(treatment = paste(filtered.vector[filtered],
                           amendment,
                           sep = "-"))


metadata.to.save.out <- all.metadata %>%
  select(bottleID, incubationID, sampleID, tripID, dateKilled,
         timeKilled, depth, t, filtered, amendment, treatment) %>%
  filter()

write.csv(metadata.to.save.out,
          "metadata/processedMetadata/incubation_Hg_metadata_2019.csv",
          row.names = FALSE,
          quote = FALSE)

rm(hg.ids, incubation.ids, sample.ids, trip.id, metadata.to.save.out)



#### Calculate time gaps in processing (spiking or sample collection) ####

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
    select(1, 3, 4, 5, 17) %>%
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
# Checked this one manually. All looks good to go ahead.

clean.MeHg.data.from("dataRaw/incubations/MeHg/I092419 BENDOTA.xlsx",
                     "dataEdited/incubations/MeHg/cleaned/incubations2019_20190924.csv")
# Checked this one manually. All looks good to go ahead.

clean.MeHg.data.from("dataRaw/incubations/MeHg/I092519 BENDOTA.xlsx",
                     "dataEdited/incubations/MeHg/cleaned/incubations2019_20190925.csv")

# clean.MeHg.data.from("dataRaw/incubations/MeHg/I120419 BENDOTA.xlsx",
#                      "dataEdited/incubations/MeHg/cleaned/incubations2019_20191204.csv",
#                      barcodes.to.remove = "MSC576AR")
# This one had a super high blank. Not gonna keep it. Re-running these samples
# in February 2020.

clean.MeHg.data.from("dataRaw/incubations/MeHg/I120519 BENDOTA.xlsx",
                     "dataEdited/incubations/MeHg/cleaned/incubations2019_20191205.csv")
# Checked this one manually. All looks good to go ahead.


clean.MeHg.data.from("dataRaw/incubations/MeHg/I020520 BENDOTA.xlsx",
                     "dataEdited/incubations/MeHg/cleaned/incubations2019_20200205.csv",
                     barcodes.to.remove = "MSC576AR")
# MSC576AR was run on 2019-12-05 as well, and they had similar results.
# So, removing this on here.
# Checked this one manually. All looks good to go ahead.



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

write.csv(MeHg.results,
          "dataEdited/incubations/MeHg/incubations2019_MeHg_all.csv",
          row.names = FALSE)

rm(list.o.results,
   file.name)



#### Clean Me198Hg data ####

# Check out spike-to-kill times and if it influences t0 values

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
rm(check_t0_timing)

# Pull out needed data from t0 samples
MeHg.198.t0 <- MeHg.results  %>%
  left_join(all.metadata) %>%
  filter(t == "t0")%>%
  left_join(processing.data) %>%
  select(bottleID, incubationID, sampleID, tripID, startDate, depth, treatment, t, excess_MeHg_198_ng.L, excess_DDL)

# Pull out t1 data

MeHg.198.production.t1 <- MeHg.results %>%
  select(bottleID, excess_MeHg_198_ng.L, excess_DDL) %>%
  left_join(all.metadata) %>%
  filter(t == "t1") %>%
  select(bottleID, incubationID, sampleID, tripID, startDate, depth, treatment, t, excess_MeHg_198_ng.L, excess_DDL)


# Combine Me198Hg data

MeHg.198.data <- rbind(MeHg.198.t0, MeHg.198.production.t1) %>%
  arrange(t) %>%
  arrange(treatment) %>%
  arrange(sampleID)

# Which ones are above detection?

MeHg.198.data <- MeHg.198.data %>%
  mutate(above_DDL = (excess_MeHg_198_ng.L > excess_DDL)) %>%
  select(-excess_DDL)



# Write out data

write.csv(MeHg.198.data,
          "dataEdited/incubations/MeHg/incubations2019_Me198Hg.csv",
          row.names = FALSE)
rm(MeHg.198.data)







#### Calculate fraction Me198Hg in each sample collected ####

# Read in THg data
THg.data <- read.csv("dataEdited/incubations/THg/incubations2019_THg.csv",
                     stringsAsFactors = FALSE)
MeHg.198.data <- read.csv("dataEdited/incubations/MeHg/incubations2019_Me198Hg.csv",
                     stringsAsFactors = FALSE)
# Calculate percent MeHg
per.MeHg.data <- full_join(MeHg.198.data,
                           THg.data) %>%
  select(tripID, incubationID, startDate, depth, treatment, t, excess_MeHg_198_ng.L, THg_198) %>%
  mutate(per_MeHg_198 = (excess_MeHg_198_ng.L / THg_198) * 100) %>%
  select(-c(excess_MeHg_198_ng.L, THg_198))

# Write out percentage MeHg data
write.csv(per.MeHg.data,
          "dataEdited/incubations/MeHg/incubations2019_Me198Hg_percent.csv",
          row.names = FALSE)









#### Calculate change in percent Me198Hg ####

# Join percent data with processing time
per.MeHg.data.change <- per.MeHg.data %>% 
  spread(key = t,
         value = per_MeHg_198) %>%
  left_join(processing.data) %>%
  # Calculate fraction of 24 hours that the samples were incubated for.
  mutate(t0_to_t1_time.fraction = (as.numeric(t0_to_t1_time) / (24*60*60))) %>%
  arrange(incubationID) %>%
  mutate(change_in_per_MeHg = (t1 - t0)) %>%
  mutate(change_in_per_MeHg_norm = round((change_in_per_MeHg / t0_to_t1_time.fraction), 3)) %>%
  select(-c(change_in_per_MeHg, t0_to_t1_time.fraction, t0_to_t1_time, spike_to_kill_0_time, t0, t1))

# Write out percentage MeHg data
write.csv(per.MeHg.data.change,
          "dataEdited/incubations/MeHg/incubations2019_Me198Hg_change_in_percent.csv",
          row.names = FALSE)

rm(list = ls(pattern = "per."),
   THg.data)








#### Calculate Me198Hg production over assay ####

MeHg.delta <- MeHg.198.data %>%
  select(-c(bottleID, above_DDL)) %>%
  spread(t, excess_MeHg_198_ng.L) %>%
  left_join(processing.data) %>%
  # Calculate fraction of 24 hours that the samples were incubated for.
  mutate(t0_to_t1_time.fraction = (as.numeric(t0_to_t1_time) / (24*60*60))) %>%
  arrange(incubationID) %>%
  # mutate(t0_to_t1_time = replace_na(t0_to_t1_time, 1)) %>%
  # Calculate change in MeHg normalized to incubation time
  mutate(Me198Hg_production = round(((t1 - t0) / t0_to_t1_time.fraction), 3)) %>%
  select(incubationID, sampleID, tripID, startDate, depth, treatment, Me198Hg_production)

# Write out data
write.csv(MeHg.delta,
          "dataEdited/incubations/MeHg/incubations2019_Me198Hg_production.csv",
          row.names = FALSE)


# Clean up

rm(list = ls(pattern = "MeHg.198"),
   MeHg.delta)




#### Clean Me204Hg data ####

# Clean t1 sample values

MeHg.204.production.t1 <- MeHg.results %>%
  select(bottleID, excess_MeHg_204_ng.L, excess_DDL) %>%
  left_join(all.metadata) %>%
  filter(t == "t1") %>%
  select(bottleID, incubationID, sampleID, tripID, startDate, depth, treatment, t, excess_MeHg_204_ng.L, excess_DDL)


# Check out spike-to-kill times and if it influences t0 values

check_t0_timing <- MeHg.results  %>%
  left_join(all.metadata) %>%
  filter(t == "t0")%>%
  full_join(processing.data) %>%
  filter(!is.na(spike_to_kill_0_time)) %>%
  select(excess_MeHg_204_ng.L, spike_to_kill_0_time)
par(mfrow = c(1,1))
plot(x = check_t0_timing$spike_to_kill_0_time,
     y = check_t0_timing$excess_MeHg_204_ng.L,
     pch = 18,
     xlab = "Spike time to kill time (sec)",
     ylab = "Excess Me204Hg in t0 sample")
# It looks like there might actually be some rapid demethylation occurring.
# Hard to tell for sure, but there does seem to be a clear trend in
# decreasing Me204Hg as time to spike increased.
rm(check_t0_timing)

# Pull out needed data from t0 samples
MeHg.204.t0 <- MeHg.results  %>%
  left_join(all.metadata) %>%
  filter(t == "t0")%>%
  left_join(processing.data) %>%
  select(bottleID, incubationID, sampleID, tripID, startDate, depth, treatment, t, excess_MeHg_204_ng.L, excess_DDL)


# Combine Me204Hg data

MeHg.204.data <- rbind(MeHg.204.t0, MeHg.204.production.t1) %>%
  arrange(t) %>%
  arrange(treatment) %>%
  arrange(sampleID)

# Which ones are above detection?

MeHg.204.data <- MeHg.204.data %>%
  mutate(above_DDL = (excess_MeHg_204_ng.L > excess_DDL)) %>%
  select(-excess_DDL)



# Write out data

write.csv(MeHg.204.data,
          "dataEdited/incubations/MeHg/incubations2019_Me204Hg.csv",
          row.names = FALSE)

# Clean up
rm(list = ls(pattern = "MeHg"))

















#### Calculate percent Me204Hg in each sample collected ####

# Read in THg data
THg.data <- read.csv("dataEdited/incubations/THg/incubations2019_THg.csv",
                     stringsAsFactors = FALSE)
# Read in MeHg data
MeHg.data <- read.csv("dataEdited/incubations/MeHg/incubations2019_Me204Hg.csv",
                     stringsAsFactors = FALSE)

# Calculate percent MeHg
per.MeHg.data <- full_join(MeHg.data,
                           THg.data) %>%
  select(tripID, incubationID, startDate, depth, treatment, t, excess_MeHg_204_ng.L, THg_204) %>%
  mutate(per_MeHg_204 = (excess_MeHg_204_ng.L / THg_204) * 100) %>%
  select(-c(excess_MeHg_204_ng.L, THg_204))

# Write out percentage MeHg data
write.csv(per.MeHg.data,
          "dataEdited/incubations/MeHg/incubations2019_Me204Hg_percent.csv",
          row.names = FALSE)










#### Calculate change in percent Me204Hg ####

# Join percent data with processing time
per.MeHg.data.change <- per.MeHg.data %>% 
  spread(key = t,
         value = per_MeHg_204) %>%
  left_join(processing.data) %>%
  # Calculate fraction of 24 hours that the samples were incubated for.
  mutate(t0_to_t1_time.fraction = (as.numeric(t0_to_t1_time) / (24*60*60))) %>%
  arrange(incubationID) %>%
  mutate(change_in_per_MeHg = ((t1 - t0) / t0) * 100) %>%
  mutate(change_in_per_MeHg_norm = round((change_in_per_MeHg / t0_to_t1_time.fraction), 3)) %>%
  select(-c(change_in_per_MeHg, t0_to_t1_time.fraction, t0_to_t1_time, spike_to_kill_0_time, t0, t1))
  
# Write out percentage MeHg data
write.csv(per.MeHg.data.change,
          "dataEdited/incubations/MeHg/incubations2019_Me204Hg_percent_change.csv",
          row.names = FALSE)

rm(list = ls(pattern = "per."),
   THg.data,
   MeHg.data)
  
  
  
  
  






#### Clean ambient Hg data ####

MeHg.results <- read.csv("dataEdited/incubations/MeHg/incubations2019_MeHg_all.csv",
                         stringsAsFactors = FALSE)

ambient.MeHg <- MeHg.results %>%
  left_join(all.metadata) %>%
  select(bottleID, incubationID, sampleID, tripID, startDate, depth, treatment, t, amb_MeHg_ng.L) %>%
  arrange(t) %>%
  arrange(treatment) %>%
  arrange(sampleID)
write.csv(ambient.MeHg,
          "dataEdited/incubations/MeHg/incubations2019_MeHg_ambient.csv",
          row.names = FALSE)


#### Predicted water column MeHg ####
WC.MeHg <- ambient.MeHg %>%
  filter(t == "t0",
         treatment == "unfiltered-unamended") %>%
  group_by(tripID, sampleID, startDate, depth) %>%
  summarise(amb_MeHg_ng.L = mean(amb_MeHg_ng.L))
write.csv(WC.MeHg,
          "dataEdited/incubations/MeHg/incubations2019_MeHg_ambient_WC.csv",
          row.names = FALSE)



#### Which samples haven't been analyzed? ####

unanalyzed.samples <- all.metadata$bottleID[which(!(all.metadata$bottleID %in% ambient.MeHg$bottleID))]
