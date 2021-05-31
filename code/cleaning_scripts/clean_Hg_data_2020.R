# code/cleaning_scripts/clean_Hg_data_2020.R

# This scripts reads in the Hg data that I received from the 
# USGS Hg team for the 2020 sampling trips and combines it
# into one data set. The data is then saved as a csv.

#### Always get a clean start ####

rm(list = ls())
setwd("~/Documents/research/BLiMMP")
library(dplyr)
library(lubridate)
library(tidyr)
library(readxl)


#### Function to clean individual data sheets ####

dataSheet <- "dataRaw/Hg/MENDOTA_results_2020-11-03.csv"
waterDepth = 24
output = "dataEdited/Hg/2020/clean_Hg_data_20201103.csv"

clean.Hg.data.for <- function(dataSheet,
                              waterDepth = 24,
                              output) {
  
  raw.Hg.data <- read.csv(dataSheet,
                          stringsAsFactors = FALSE)
  
  clean.Hg.data <- raw.Hg.data %>%
    rename(depthOriginal = Depth,
           concentration = Result.Value,
           constituent = Constituent,
           sampleTime = Sample.Time,
           unit = Unit) %>%
    mutate(corewater = depthOriginal < 0,
           depth = depthOriginal + (waterDepth * corewater),
           sampleDate = mdy(Sample.Date)) %>%
    select(sampleDate, depth, sampleTime, depthOriginal, concentration,
           unit, constituent, DDL, corewater)
  
  write.csv(clean.Hg.data,
            output,
            row.names = FALSE)
  }




#### Clean individual sheets ####

clean.Hg.data.for(dataSheet = "dataRaw/Hg/MENDOTA_results_2020-10-09.csv",
                  waterDepth = 24,
                  output = "dataEdited/Hg/2020/clean_Hg_data_20201009.csv")
clean.Hg.data.for(dataSheet = "dataRaw/Hg/MENDOTA_results_2020-11-03.csv",
                  waterDepth = 24,
                  output = "dataEdited/Hg/2020/clean_Hg_data_20201103.csv")
clean.Hg.data.for(dataSheet = "dataRaw/Hg/MENDOTA_results_2021-01-19.csv",
                  waterDepth = 24,
                  output = "dataEdited/Hg/2020/clean_Hg_data_20210119.csv")
clean.Hg.data.for(dataSheet = "dataRaw/Hg/MENDOTA_results_2021-05-06_editedHeaders_noIncubationSamples.csv",
                  waterDepth = 24,
                  output = "dataEdited/Hg/2020/clean_Hg_data_20210506.csv")



rm(list = ls())


#### Aggregate all cleaned data files for 2020 ####
list.o.files <- list.files(path = "./dataEdited/Hg/2020",
                           pattern = "clean_Hg_data",
                           full.names = TRUE)
Hg.data <- do.call(rbind,
                   lapply(list.o.files,
                          function(x) {
                            read.csv(x,
                                     stringsAsFactors = FALSE)
                          }))


#### Prep inorganic Hg data ####
clean.Hg.data <- Hg.data %>%
  mutate(constituent = paste(constituent, "_", unit,
                             sep = "") %>%
           gsub("/", ".", .)) %>%
  group_by(sampleDate, depth, constituent, corewater) %>%
  summarize(concentration = mean(concentration)) %>%
  spread(key = constituent,
         value = concentration) %>%
  mutate(FiHg_NG.L = FTHG_NG.L - FMHG_NG.L,
         PiHg_NG.L = PTHG_NG.L - PMHG_NG.L) %>%
  gather(key = constituent,
         value = concentration,
         -c(1:3)) %>%
  filter(!is.na(concentration))
  




#### Read out data ####
write.csv(Hg.data,
          file = "dataEdited/Hg/Hg_data_2020.csv",
          row.names = FALSE)


