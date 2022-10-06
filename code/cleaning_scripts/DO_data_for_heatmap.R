#### code/cleaning_scripts/clean_exo_data.R ####
# Benjamin D. Peterson

# This scripts reads in the sonde data collected on the Exo
# sonde. The data was collected as part of the BLiMMP project.
# These samples were collected from Lake Mendota


#### This is your mess, you need to clean it up ####
rm(list = ls())
setwd("~/Documents/research/BLiMMP")
library(dplyr)
library(lubridate)
library(readxl)



#### Read in and combine data ####
sonde.files <- list.files("dataRaw/MeMO_exo",
                          pattern = "_ME_profile.csv")

for (position.of.interest in 1:length(sonde.files)) {
  file.of.interest <- sonde.files[[position.of.interest]]
  # date.of.interest <- file.of.interest %>%
  #   strsplit("_") %>% sapply("[", 1) %>%
  #   ymd() %>% as.character()
  # 
  if (file.of.interest == sonde.files[1]) {
    exo.data <- read.csv(paste("dataRaw/MeMO_exo/",
                               file.of.interest,
                               sep = "")) %>%
      select(sampledate, depth, wtemp, do_raw, turb_fnu)
  } else {
    exo.data <- rbind(exo.data,
                      read.csv(paste("dataRaw/MeMO_exo/",
                                     file.of.interest,
                                     sep = "")) %>%
                        select(sampledate, depth, wtemp, do_raw, turb_fnu))
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
exo.data <- exo.data %>%
  filter(depth <= 20)

#### Save out data ####
saveRDS(exo.data,
        "dataEdited/MeMO_exo/MeMO_exo_data.rds")
