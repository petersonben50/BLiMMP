#### code/BGC_profiles/heatmaps_sampling_locations.R ####
# Benjamin D. Peterson

# This script will generate heatmaps of the temperature
# and DO of the lake over the course of the ice-free
# season and include markers for our sampling locations.


#### Get set up ####
rm(list = ls())
setwd("~/Documents/research/BLiMMP")
library(lubridate)
library(readxl)
library(tidyverse)


#### Download data ####
# Go for the hourly data.
# Find it here: https://lter.limnology.wisc.edu/node/56136/data_form
# Download all columns, from 2020-04-01 to 2020-12-31.
# Save it here: ~/Documents/research/BLiMMP/dataRaw/buoy/sensor_mendota_lake_watertemp_hourly.csv
# Got the 2021 hourly data from Mark Gahler on 2022-10-06

#### Set variables of interest for testing ####
buoy.data.2020 <- read.csv(file = "dataRaw/buoy/sensor_mendota_lake_watertemp_hourly.csv") %>%
  select(sampledate, hour, depth, wtemp, flag_wtemp)
buoy.data.2021 <- read.csv(file = "dataRaw/buoy/me_wtemp_2021.csv") %>%
  mutate(hour = gsub(":", "", sampletime) %>%
           substr(1, 4) %>%
           as.integer()) %>%
  select(sampledate, hour, depth, wtemp, flag_wtemp)
data.to.use <- rbind(buoy.data.2020,
                     buoy.data.2021)


#### Save out data ####
write.csv(data.to.use,
          file = "dataFinal/temperature_buoy_data.csv",
          row.names = FALSE)
