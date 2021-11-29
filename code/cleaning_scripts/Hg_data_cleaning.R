#### code/cleaning_scripts/Hg_data_cleaning.R ####

# This script used to read in the Hg data that I received from the 
# USGS Hg team for the 2020 sampling trips and combine it
# into one data set. Now I've changed it, I downloaded all the Hg
# data related to my project from Merlins into one big file. Everything
# from the Lake Mendota project after 2020-07-01. This is saved here:
# dataRaw/Hg/USGS_results_2021-11-09.xlsx

# This file includes the incubation results as well.


#### Always get a clean start ####
rm(list = ls())
setwd("~/Documents/research/BLiMMP")
library(dplyr)
library(lubridate)
library(tidyr)
library(readxl)


#### Set variables ####
dataSheet <- "dataRaw/Hg/USGS_results_2021-11-09.xlsx"
waterDepth <- 24
output <- "dataEdited/Hg/Hg_data_clean.csv"


#### Read in data ####
Hg.data.raw <- read_xlsx(dataSheet) %>%
  rename(bottleID = bottle,
         depthOriginal = depth,
         concentration = result_value,
         sampleTime = sample_time) %>%
  mutate(corewater = depthOriginal < 0,
         depth = depthOriginal + (waterDepth * corewater),
         sampleDate = ymd(sample_date)) %>%
  select(sampleDate, depth, sampleTime, depthOriginal, concentration,
         unit, constituent, detection_flag, corewater)


#### Keep only water column data ####
# No corewater or incubation data
Hg.data <- Hg.data.raw %>%
  filter(!(constituent %in% c("UMHG", "UTHG")),
         !(sampleDate %in% c(as.Date("2020-08-14"),
                             as.Date("2020-09-24"))))
rm(Hg.data.raw)

#### Prep inorganic data ####
Hg.data <- Hg.data %>%
  mutate(constituent = paste(constituent, "_", unit,
                             sep = "") %>%
           gsub("/", ".", .)) %>%
  select(-unit)

iHg.data <- Hg.data %>%
  select(-detection_flag) %>%
  filter(constituent %in% c("FTHG_NG.L", "FMHG_NG.L",
                            "PTHG_NG.L", "PMHG_NG.L")) %>%
  spread(key = constituent,
         value = concentration) %>%
  mutate(FiHg_NG.L = FTHG_NG.L - FMHG_NG.L,
         PiHg_NG.L = PTHG_NG.L - PMHG_NG.L) %>%
    gather(key = constituent,
         value = concentration,
         c(6:11)) %>%
  filter(constituent %in% c("FiHg_NG.L",
                            "PiHg_NG.L")) %>%
  mutate(detection_flag = "NA") %>%
  select(sampleDate, depth, sampleTime, depthOriginal,
         concentration, constituent, detection_flag, corewater)

Hg.data <- rbind(Hg.data,
                 iHg.data)
rm(iHg.data)


#### Fix surface samples from 2021-09-10/11 ####
# We have two surface samples listed... not sure why,
# and not sure why one is listed on 2021-09-11. Let's 
# just remove that one, we only need one here.
Hg.data <- Hg.data %>%
  filter(!(sampleDate == "2021-09-11" & depth == 0.0 & sampleTime == "1300"))

#### Read out data ####
write.csv(Hg.data,
          file = output,
          row.names = FALSE)
