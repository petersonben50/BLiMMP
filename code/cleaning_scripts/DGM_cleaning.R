#### code/cleaning_scripts/DGM_cleaning.R ####
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("~/Documents/research/BLiMMP/")
library(lubridate)
library(readxl)
library(tidyverse)
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")
source("code/BLiMMP_functions.R")

#### Add metadata ####
sample.metadata <- read_xlsx("metadata/1_trip_IDs.xlsx") %>%
  select(tripID, startDate) %>%
  left_join(read_xlsx("metadata/2_sample_IDs.xlsx")) %>%
  select(sampleID, depth, startDate) %>%
  filter(year(startDate) == 2021)
DGM.metadata <- read_xlsx("metadata/Hg_DGM.xlsx") %>%
  mutate(trapNumber = paste("Trap ", trapNumber, sep = "")) %>%
  select(sampleID, trapNumber) %>%
  left_join(sample.metadata) %>%
  mutate(depth = as.numeric(depth))
rm(sample.metadata)


#### Read in September data ####
DGM.data.Sept <- read_xlsx("dataEdited/Hg/DGM_20210916.xlsx",
                      skip = 23)
names(DGM.data.Sept) <- c("sampleID", "peakArea", "mass",
                     "volumeAnalyzed", "concentration_ng.L",
                     "percentRecovery", "notes")
DGM.data.Sept.clean <- DGM.data.Sept %>%
  select(sampleID, concentration_ng.L) %>%
  filter(!grepl("STD", sampleID),
         !grepl("Blank", sampleID),
         !grepl("QCS", sampleID)) %>%
  mutate(trapNumber = sampleID %>%
           strsplit(" - ") %>% sapply("[", 1)) %>%
  select(trapNumber, concentration_ng.L) %>%
  left_join(DGM.metadata %>%
              filter(month(startDate) == 9)) %>%
  select(startDate, depth, concentration_ng.L)


#### Read in October data ####
DGM.data.Oct <- read_xlsx("dataEdited/Hg/DGM_20211018.xlsx",
                      skip = 23)
names(DGM.data.Oct) <- c("trapNumber", "peakArea", "mass",
                         "volumeAnalyzed", "concentration_ng.L",
                         "percentRecovery", "notes")
DGM.data.Oct.clean <- DGM.data.Oct %>%
  select(trapNumber, concentration_ng.L) %>%
  filter(!grepl("ul", trapNumber),
         !grepl("QCS", trapNumber)) %>%
  left_join(DGM.metadata %>%
              filter(month(startDate) == 10)) %>%
  select(startDate, depth, concentration_ng.L)


#### Join metadata with data ####
all.data <- rbind(DGM.data.Sept.clean,
                  DGM.data.Oct.clean)



#### Save out data ####
saveRDS(all.data,
        "dataEdited/Hg/DGM_2021_data.rds")
