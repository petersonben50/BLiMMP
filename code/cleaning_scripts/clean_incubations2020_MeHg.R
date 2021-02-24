#### code/cleaning_scripts/clean_incubations2020_MeHg.R ####
# Written for BLiMMP project
# Benjamin D. Peterson


#### Prep workspace ####

rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP")
library(dplyr)
library(lubridate)
library(readxl)
library(tidyr)
source("code/cleaning_scripts/cleaning_functions.R")


#### Read in metadata ####
all.metadata.2020 <- read.csv("metadata/processedMetadata/incubation_Hg_metadata_2020.csv",
                              stringsAsFactors = FALSE)


#### Cleaning data ####

# BLiMMP Lab Notebook 2 pg 64
clean.MeHg.data.from("dataRaw/incubations/MeHg/2020_incubations/I102920 BENDOTA.xlsx",
                     "dataEdited/incubations/MeHg/cleaned/incubations2020_20201029.csv",
                     date.of.analysis = "2020-10-29")
# BLiMMP Lab Notebook 2 pg 65
clean.MeHg.data.from("dataRaw/incubations/MeHg/2020_incubations/I103020 BENDOTA.xlsx",
                     "dataEdited/incubations/MeHg/cleaned/incubations2020_20201030.csv",
                     barcodes.to.remove = c("MSC260AW",
                                            "MSC262AW",
                                            "MSC227AW"),
                     date.of.analysis = "2020-10-30")
# BLiMMP Lab Notebook 2 pg 69
clean.MeHg.data.from("dataRaw/incubations/MeHg/2020_incubations/I111220 BENDOTA.xlsx",
                     "dataEdited/incubations/MeHg/cleaned/incubations2020_20201112.csv",
                     date.of.analysis = "2020-11-12")

# BLiMMP Lab Notebook 2 pg 75
# clean.MeHg.data.from("dataRaw/incubations/MeHg/2020_incubations/I120220 BENDOTA.xlsx",
#                      "dataEdited/incubations/MeHg/cleaned/incubations2020_20201202.csv")
# 
# 
# # This is shit don't keep this
# clean.MeHg.data.from("dataRaw/incubations/MeHg/2020_incubations/I120320 BENDOTA.xlsx",
#                      "dataEdited/incubations/MeHg/cleaned/incubations2020_20201203.csv")
# clean.MeHg.data.from("dataRaw/incubations/MeHg/2020_incubations/suspect/I121120 BENDOTA.xlsx",
#                      "dataEdited/incubations/MeHg/cleaned/incubations2020_20201211.csv")

# BLiMMP Lab Notebook 2 pg 
clean.MeHg.data.from("dataRaw/incubations/MeHg/2020_incubations/I011321 BENDOTA.xlsx",
                     "dataEdited/incubations/MeHg/cleaned/incubations2020_20210113.csv",
                     date.of.analysis = "2021-01-13")

# BLiMMP Lab Notebook 2 pg 
clean.MeHg.data.from("dataRaw/incubations/MeHg/2020_incubations/I011421 BENDOTA.xlsx",
                     "dataEdited/incubations/MeHg/cleaned/incubations2020_20210114.csv",
                     barcodes.to.remove = "MSC956AR",
                     date.of.analysis = "2021-01-14")
# MSC956AR had only ~30ml of distillate. Re-ran this sample on 2021-01-15
# BLiMMP Lab Notebook 2 pg 
clean.MeHg.data.from("dataRaw/incubations/MeHg/2020_incubations/I011521 BENDOTA.xlsx",
                     "dataEdited/incubations/MeHg/cleaned/incubations2020_20210115.csv",
                     date.of.analysis = "2021-01-15")



#### Combine MeHg data ####

rm(clean.MeHg.data.from)

# Load up all incubation data
list.o.results <- list.files(path = "dataEdited/incubations/MeHg/cleaned",
                             pattern = "incubations2020",
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


#### Combine with metadata ####

MeHg.meta.data <- left_join(MeHg.results,
                            all.metadata.2020) %>%
  select(bottleID, incubationID, tripID, depth, t, durationInDays,
         dateKilled, timeKilled, filtered, amendment, treatment,
         amb_MeHg_ng.L, excess_MeHg_198_ng.L, excess_MeHg_204_ng.L,
         above_DDL_MeHg, dateAnalyzed_MeHg)



#### Write out data ####
write.csv(MeHg.meta.data,
          file = "dataEdited/incubations/MeHg/incubations2020_MeHg.csv",
          row.names = FALSE)

# Which samples need to be analyzed yet?
samples.to.be.run <- full_join(MeHg.results,
                               all.metadata.2020) %>%
  filter(is.na(amb_MeHg_ng.L)) %>%
  arrange(bottleID)
write.csv(samples.to.be.run,
          file = "dataEdited/incubations/MeHg/incubations2020_MeHg_samplesToBeRun.csv",
          row.names = FALSE)
