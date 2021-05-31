# code/waterChemistry/clean_sulfate_data.R

# This script will clean up the sulfate data we
# generated in WSEL on the IC and combine it with
# the metadata

#### Always start with a clean slate ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP/")
library(dplyr)
library(readxl)
library(tidyr)
source("code/cleaning_scripts/cleaning_functions.R")


#### Prepare metadata ####
S.metadata <- read_xlsx("metadata/chem_S.xlsx") %>%
  select(-notes) %>%
  filter(!(sulfurID %in% c("BLI20_TS_043",
                           "BLI20_TS_044",
                           "BLI20_TS_045")))


#### Generate needed processing metadata ####
processing.metadata <- S.metadata %>%
  select(sulfurID, mass, tare, preservativeVol) %>%
  filter(grepl("BLI2", sulfurID)) %>%
  mutate_at(.vars = c("mass", "tare", "preservativeVol"),
            as.numeric)
rm(S.metadata)


#### Run cleaning function ####
clean.sulfate.data.from(data.file.name = "dataRaw/waterChemistry/sulfate/20210316_BENDOTA.xlsx",
                        processing.metadata.df = processing.metadata,
                        output.file.name = "dataEdited/waterChemistry/sulfate/20210316_BENDOTA.csv")
clean.sulfate.data.from(data.file.name = "dataRaw/waterChemistry/sulfate/20210423_BENDOTA.xlsx",
                        processing.metadata.df = processing.metadata,
                        samples.to.remove = c("BLI20_TS_087", "BLI20_TS_088", "BLI20_TS_089", "BLI20_TS_090",
                                              "BLI20_TS_091", "BLI20_TS_092", "BLI20_TS_092_duplicate",
                                              "BLI20_TS_093", "BLI20_TS_094"),
                        output.file.name = "dataEdited/waterChemistry/sulfate/20210423_BENDOTA.csv")


#### Clean up before combining all samples ####

rm(list = ls())


# Load up all incubation data
list.o.results <- list.files(path = "dataEdited/waterChemistry/sulfate/",
                             pattern = "_BENDOTA.csv",
                             full.names = TRUE)
for (file.name in list.o.results) {
  if (file.name == list.o.results[1]) {
    S.results <- read.csv(file.name,
                          stringsAsFactors = FALSE)
  } else {
    S.results <- rbind(S.results,
                       read.csv(file.name,
                                stringsAsFactors = FALSE))
  }
  
}



#### Combine all water column data ####
WC.metadata <- read.csv("metadata/processedMetadata/sulfide_WC.csv",
                        stringsAsFactors = FALSE)
waterDepth = 24

WC.results <- S.results %>%
  inner_join(WC.metadata) %>%
  arrange(sulfurID) %>%
  mutate(depthOriginal = depth) %>%
  mutate(corewater = grepl(pattern = "-",
                           x = depthOriginal))
WC.results[WC.results$corewater, ] <- WC.results[WC.results$corewater, ] %>%
  mutate(depth = paste("-",
                       strsplit(depth, "-") %>% sapply("[", 2),
                       sep = "") %>%
           as.numeric() / 100) %>%
  mutate(depth = depth + (waterDepth * corewater))

# Save out data
write.csv(WC.results,
          "dataEdited/waterChemistry/sulfate/WC_data.csv",
          row.names = FALSE,
          quote = FALSE)

