#### code/cleaning_scripts/clean_NA_data.R ####
# Benjamin D. Peterson


#### Keep the order in the kingdom ####
rm(list = ls())
setwd("~/Box/BLiMMP/")
library(readxl)
library(tidyverse)


#### Read in needed metadata ####

NA.metadata <- read_xlsx("metadata/NA_IDs.xlsx")
sample_IDs <- read_xlsx("metadata/2_sample_IDs.xlsx")
trip_IDs <- read_xlsx("metadata/1_trip_IDs.xlsx") %>%
  select(-notes)


#### Join needed metadata ####

NA.metadata.all <- NA.metadata %>%
  left_join(sample_IDs) %>%
  left_join(trip_IDs) %>% 
  select(filterID, sampleID, tripID, startDate,
         depth, replicate, volumeFiltered)
  

#### Write out data ####

write.csv(NA.metadata.all,
          "metadata/processedMetadata/NA_metadata.csv",
          row.names = FALSE)

