#### code/cleaning_scripts/Hg_data_cleaning.R ####

# This script used to read in the Hg data that I received from the 
# USGS Hg team for the 2020 sampling trips and combine it
# into one data set. Now I've changed it, I downloaded all the Hg
# data related to my project from Merlins into one big file. Everything
# from the Lake Mendota project after 2020-07-01. This is saved here:
# dataRaw/Hg/USGS_results_2022-09-07.xlsx

# This file includes the incubation results as well.


#### Always get a clean start ####
rm(list = ls())
setwd("~/Documents/research/BLiMMP")
library(dplyr)
library(lubridate)
library(tidyr)
library(readxl)


#### Set variables ####
dataSheet <- "dataRaw/Hg/USGS_results_2022-09-07.xlsx"
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


#### Remove unneeded data ####
# No corewater or incubation data
Hg.data <- Hg.data.raw %>%
  filter(
    # No incubation data
    !(constituent %in% c("UMHG", "UTHG")),
    # No data from corewater sampling trips
    !(sampleDate %in% c(as.Date("2020-08-14"),
                        as.Date("2020-09-24"))),
    # We have two surface samples listed... not sure why,
    # and not sure why one is listed on 2021-09-11. Let's 
    # just remove that one, we only need one here.
    !(sampleDate == "2021-09-11" & depth == 0.0 & sampleTime == "1300")
    )
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



#### Calculate % MeHg for filtered and particulate ####
perMEHG.data <- Hg.data %>%
  filter(constituent %in% c("FMHG_NG.L", "FTHG_NG.L",
                            "PMHG_NG.L", "PTHG_NG.L"),
         detection_flag == "NONE") %>%
  select(-detection_flag) %>%
  spread(key = constituent,
         value = concentration) %>%
  mutate(perFMHG = FMHG_NG.L / FTHG_NG.L * 100,
         perPMHG = PMHG_NG.L / PTHG_NG.L * 100) %>%
  gather(key = constituent,
         value = concentration,
         -c(1:5)) %>%
  filter(!is.na(concentration),
         constituent %in% c("perFMHG", "perPMHG")) %>%
  mutate(detection_flag = "NONE") %>%
  select(sampleDate, depth, sampleTime, depthOriginal,
         concentration, constituent, detection_flag, corewater)
Hg.data <- rbind(Hg.data,
                 perMEHG.data)
rm(perMEHG.data)


#### Calculate solids concentrations for MeHg and HgT ####
solids.HG.data <- Hg.data %>%
  filter(constituent %in% c("SPM_MG.L", "PTHG_NG.L",
                            "PMHG_NG.L", "PiHg_NG.L")) %>%
  select(-detection_flag) %>%
  spread(key = constituent,
         value = concentration) 
# One of the duplicates is missing an SPM values.
# Let's just use the SPM value from the other dup
# to calculate the per g values.
solids.HG.data[(solids.HG.data$sampleDate == "2020-10-10" &
                  solids.HG.data$sampleTime == "1800"), "SPM_MG.L"] <- solids.HG.data[(solids.HG.data$sampleDate == "2020-10-10" &
                                                                                         solids.HG.data$sampleTime == "1820"), "SPM_MG.L"]
solids.HG.data <- solids.HG.data %>%
  mutate(PMHG_NG.G = PMHG_NG.L / SPM_MG.L * 1000,
         PTHG_NG.G = PTHG_NG.L / SPM_MG.L * 1000,
         PiHG_NG.G = PiHg_NG.L / SPM_MG.L * 1000) %>%
  gather(key = constituent,
         value = concentration,
         -c(1:5)) %>%
  filter(!is.na(concentration),
         constituent %in% c("PMHG_NG.G", "PTHG_NG.G", "PiHG_NG.G")) %>%
  mutate(detection_flag = "NONE") %>%
  select(sampleDate, depth, sampleTime, depthOriginal,
         concentration, constituent, detection_flag, corewater)
# Add solids concentrations to dataset
Hg.data <- rbind(Hg.data,
                 solids.HG.data)


#### Change dates so data from a single trip are listed on one day ####
# This is how the other geochemical data is stored.
unique(Hg.data$sampleDate)
Hg.data$sampleDate[which(Hg.data$sampleDate == "2021-09-11")] <- "2021-09-10"
Hg.data$sampleDate[which(Hg.data$sampleDate == "2020-09-03")] <- "2020-09-02"



#### Read out data ####
write.csv(Hg.data,
          file = output,
          row.names = FALSE)
