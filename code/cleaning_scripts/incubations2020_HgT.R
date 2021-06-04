#### code/cleaning_scripts/clean_incubations2020_THg.R ####
# Written for BLiMMP project
# Benjamin D. Peterson


#### Clean up! ####

rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP")
library(lubridate)
library(readxl)
library(tidyverse)
source("code/cleaning_scripts/cleaning_functions.R")




#### Generate needed metadata ####
all.metadata.2020 <- read.csv("metadata/processedMetadata/incubation_Hg_metadata_2020.csv")



#### Clean data ####
# See notes in Obsedian file: `BLiMMP 2020 incubation HgT analyses`
clean.HgT.data.from(file.name.input = "dataRaw/incubations/HgT/2020_incubations/I020321 BENDOTA WATER.xlsx",
                     output.file.name = "dataEdited/incubations/HgT/2020_incubations/20210203.csv",
                     metadata.df = all.metadata.2020,
                     date.of.analysis = "2021-02-03",
                     barcodes.to.remove = NULL)

clean.HgT.data.from(file.name.input = "dataRaw/incubations/HgT/2020_incubations/I020421 BENDOTA WATER.xlsx",
                    output.file.name = "dataEdited/incubations/HgT/2020_incubations/20210204.csv",
                    metadata.df = all.metadata.2020,
                    date.of.analysis = "2021-02-04",
                    barcodes.to.remove = NULL)

clean.HgT.data.from(file.name.input = "dataRaw/incubations/HgT/2020_incubations/I020521 BENDOTA WATER.xlsx",
                    output.file.name = "dataEdited/incubations/HgT/2020_incubations/20210205.csv",
                    metadata.df = all.metadata.2020,
                    date.of.analysis = "2021-02-05",
                    barcodes.to.remove = NULL)

clean.HgT.data.from(file.name.input = "dataRaw/incubations/HgT/2020_incubations/I021821 BENDOTA WATER.xlsx",
                    output.file.name = "dataEdited/incubations/HgT/2020_incubations/20210218.csv",
                    metadata.df = all.metadata.2020,
                    date.of.analysis = "2021-02-18",
                    barcodes.to.remove = NULL)

clean.HgT.data.from(file.name.input = "dataRaw/incubations/HgT/2020_incubations/I021921 BENDOTA WATER.xlsx",
                    output.file.name = "dataEdited/incubations/HgT/2020_incubations/20210219.csv",
                    metadata.df = all.metadata.2020,
                    date.of.analysis = "2021-02-19",
                    barcodes.to.remove = NULL)

clean.HgT.data.from(file.name.input = "dataRaw/incubations/HgT/2020_incubations/I022021 BENDOTA WATER.xlsx",
                    output.file.name = "dataEdited/incubations/HgT/2020_incubations/20210220.csv",
                    metadata.df = all.metadata.2020,
                    date.of.analysis = "2021-02-20",
                    barcodes.to.remove = NULL)

clean.HgT.data.from(file.name.input = "dataRaw/incubations/HgT/2020_incubations/I022221 BENDOTA WATER.xlsx",
                    output.file.name = "dataEdited/incubations/HgT/2020_incubations/20210222.csv",
                    metadata.df = all.metadata.2020,
                    date.of.analysis = "2021-02-22",
                    barcodes.to.remove = NULL)

clean.HgT.data.from(file.name.input = "dataRaw/incubations/HgT/2020_incubations/I022321 BENDOTA WATER.xlsx",
                    output.file.name = "dataEdited/incubations/HgT/2020_incubations/20210223.csv",
                    metadata.df = all.metadata.2020,
                    date.of.analysis = "2021-02-23",
                    barcodes.to.remove = NULL)


#### Combine HgT data ####

# rm(clean.MeHg.data.from)

# Load up all incubation data
list.o.results <- list.files(path = "dataEdited/incubations/HgT/2020_incubations",
                             pattern = ".csv",
                             full.names = TRUE)

for (file.name in list.o.results) {
  if (file.name == list.o.results[1]) {
    HgT.results <- read.csv(file.name,
                            stringsAsFactors = FALSE)
  } else {
    HgT.results <- rbind(HgT.results,
                         read.csv(file.name,
                                  stringsAsFactors = FALSE))
  }
  
}


#### Look for samples below DL ####
# We really shouldn't have any of these...
HgT.results.below.DDL <- HgT.results %>%
  filter(above_DDL_HgT_198 == FALSE |
           above_DDL_HgT_204 == FALSE)

#### Check duplicated samples ####
HgT.results %>%
  filter(bottleID %in% HgT.results$bottleID[duplicated(HgT.results$bottleID)]) %>%
  arrange(bottleID)
# These look pretty much the same between runs. Going to just stick with
# the first of each run.
HgT.results <- HgT.results %>%
  filter(!duplicated(bottleID))


#### Write out data ####
write.csv(HgT.results,
          file = "dataEdited/incubations/HgT/incubations2020_HgT.csv",
          row.names = FALSE)


#### Read in MeHg data ####
MeHg.data <- read.csv("dataEdited/incubations/MeHg/incubations2020_MeHg.csv")


#### Combine with MeHg and metadata ####
Hg.results.metadata <- left_join(MeHg.data,
                                 HgT.results)


#### Check fraction of HgT as MeHg ####
Hg.results.metadata <- Hg.results.metadata %>%
  mutate(fraction_MeHg_amb = amb_MeHg_ng.L / amb_HgT_ng.L,
         fraction_MeHg_198 = excess_MeHg_198_ng.L / excess_HgT_198_ng.L,
         fraction_MeHg_204 = excess_MeHg_204_ng.L / excess_HgT_204_ng.L)


#### Check spike samples ####
HgT.results %>%
  filter(!(bottleID %in% MeHg.data$bottleID))


#### Remove errant incubations from October metalimnion dataset ####
Hg.results.metadata <- Hg.results.metadata %>%
  filter(!(incubationID %in% c("BLI20_MA_050",
                               "BLI20_MA_051")))


#### Write out Hg data ####
write.csv(Hg.results.metadata,
          "dataEdited/incubations/2020incubations_Hg_data.csv",
          row.names = FALSE)

# Which samples need to be analyzed yet?
samples.to.be.run <- Hg.results.metadata %>%
  filter(is.na(excess_HgT_198_ng.L)) %>%
  arrange(bottleID)
write.csv(samples.to.be.run,
          file = "dataEdited/incubations/HgT/incubations2020_HgT_samplesToBeRun.csv",
          row.names = FALSE)
