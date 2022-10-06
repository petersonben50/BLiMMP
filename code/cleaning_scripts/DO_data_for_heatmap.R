#### code/cleaning_scripts/DO_data_for_heatmap.R ####
# Benjamin D. Peterson

# This script prepares a set of DO data to generate
# a heatmap of DO over 2020 and 2021 summers.
# For 2020, this data was collected by Mark Gahler.
# In 2021, I used the Mendota Microbial Observatory
# data through August 2nd, when the sonde broke.
# I'm really just interested in the onset of anoxia,
# so that will be sufficient for my needs for now.


#### This is your mess, you need to clean it up ####
rm(list = ls())
setwd("~/Documents/research/BLiMMP")
library(dplyr)
library(lubridate)
library(readxl)



#### Read in and combine data ####
sonde.files <- list.files("dataRaw/DO_profiles",
                          pattern = "_ME_profile.csv")

for (position.of.interest in 1:length(sonde.files)) {
  file.of.interest <- sonde.files[[position.of.interest]]
  if (file.of.interest == sonde.files[1]) {
    exo.data <- read.csv(paste("dataRaw/DO_profiles/",
                               file.of.interest,
                               sep = "")) %>%
      select(sampledate, depth, do_raw)
  } else {
    exo.data <- rbind(exo.data,
                      read.csv(paste("dataRaw/DO_profiles/",
                                     file.of.interest,
                                     sep = "")) %>%
                        select(sampledate, depth, do_raw))
  }
}


#### Fix up sample date column ####
exo.data <- exo.data %>%
  rename(sampleDate = sampledate) %>%
  mutate(sampleDate = mdy(sampleDate))


#### Remove negative DO values ####
# Have a few DO values that are below 0, which is instrument
# error. Change those values to 0.
exo.data$do_raw[which(exo.data$do_raw < 0)] <- 0


#### Only keep the values down to 20 m ####
exo.data.2020 <- exo.data %>%
  filter(depth <= 20)

rm(exo.data,
   file.of.interest,
   position.of.interest,
   sonde.files)





#### Read in and combine MeMO data from 2021 ####
sonde.files <- list.files("dataRaw/DO_profiles/2021",
                          pattern = "YSI_MO_")

for (position.of.interest in 1:length(sonde.files)) {
  file.of.interest <- sonde.files[[position.of.interest]]
  if (file.of.interest == sonde.files[1]) {
    MMO.data <- read.csv(paste("dataRaw/DO_profiles/2021/",
                               file.of.interest,
                               sep = "")) %>%
      mutate(sampleDate = mdy_hms(Timestamp) %>% as.Date(),
             depth = as.numeric(Folder),
             do_raw = Dissolved.Oxygen..mg.L.) %>%
      select(sampleDate, depth, do_raw)
  } else {
    MMO.data <- rbind(MMO.data,
                      read.csv(paste("dataRaw/DO_profiles/2021/",
                                     file.of.interest,
                                     sep = "")) %>%
                        mutate(sampleDate = mdy_hms(Timestamp) %>% as.Date(),
                               depth = as.numeric(Folder),
                               do_raw = Dissolved.Oxygen..mg.L.) %>%
                        select(sampleDate, depth, do_raw))
  }
}


#### Remove all data after 2021-08-02 ####
# The sonde crapped out after this date
MMO.data <- MMO.data %>%
  filter(sampleDate <= as.Date("2021-08-02"))



#### All data ####
all.data <- rbind(MMO.data,
                  exo.data.2020)

#### Save out data ####
saveRDS(all.data,
        "dataEdited/DO_profiles/DO_profiles_data.rds")
