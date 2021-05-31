
##############################
##############################
##############################
##############################
##############################
##############################
##############################
##############################
##############################
##############################
##############################
##############################
##############################
##############################
##############################
##############################

# Make sure I get the tare values on the tubes

##############################
##############################
##############################
##############################
##############################
##############################
##############################
##############################
##############################
##############################
##############################
##############################
##############################
##############################
##############################
##############################
##############################
##############################

rm(list = ls())
setwd("~/Documents/research/BLiMMP/")
library(readxl)
library(tidyverse)


#### Read in metadata ####
WC.metadata <- read.csv("metadata/processedMetadata/metals_WC.csv") %>%
  select(metalID, sampleID, tripID, startDate, filteredInField,
         depth, tare, mass, replicate, preservativeVol)


#### Read in WC ICP-MS data ####
icpms.wc.data <- read_xlsx("dataEdited/waterChemistry/ICP/unprocessed/Bendota_ICPMS_results043021.xlsx",
                           sheet = "metals_WC") %>%
  filter(!grepl(" dup", metalID)) %>%
  rename(Mn_ppb = `[Mn] ppb`,
         Fe_ppb = `[Fe] ppb`,
         Co_ppb = `[Co] ppb`) %>%
  select(metalID, Mn_ppb, Fe_ppb, Co_ppb)


#### Join all data ####
all.data <- full_join(WC.metadata,
                      icpms.wc.data) 
rm(icpms.wc.data)


#### Convert cobalt below DL to 0 ####
all.data$Co_ppb[grep("<", all.data$Co_ppb)] <- 0
all.data$Co_ppb <- as.numeric(all.data$Co_ppb)


#### Adjust concentrations for dilution with preservative ####
all.data <- all.data %>%
  select(-metalID) %>%
  gather(key = constituent,
         value = concentration,
         c(Mn_ppb, Fe_ppb, Co_ppb)) %>%
  rename(concentration_diluted = concentration) %>%
  mutate(concentration = (concentration_diluted * (mass - tare)) / (mass - tare - preservativeVol)) %>%
  select(-c(tare, mass, preservativeVol))


#### Check out duplicates ####
dup.samples <- all.data %>%
  filter(replicate == "yes") %>%
  select(sampleID) %>%
  unlist(use.names = FALSE) %>%
  unique()
all.data %>%
  filter(sampleID %in% dup.samples)
# These are all pretty close. We'll just keep rep #1
all.data <- all.data %>%
  filter(is.na(replicate)) %>%
  select(-replicate)


#### Calculate particulate Mn ####
all.data.clean <- all.data %>%
  select(-concentration_diluted) %>%
  spread(key = filteredInField,
         value = concentration) %>%
  mutate(particulate = no - yes) %>%
  rename(dissolved = yes) %>%
  select(-no) %>%
  gather(key = state,
         value = concentration,
         c(6,7)) %>%
  mutate(constituent = paste(state, "_", constituent,
                             sep = "")) %>%
  select(-state)


#### Prepare depth for corewater samples ####
waterDepth <- 24
all.data.clean <- all.data.clean %>%
  mutate(depthOriginal = depth) %>%
  mutate(corewater = grepl(pattern = "-",
                           x = depthOriginal))
all.data.clean[all.data.clean$corewater, ] <- all.data.clean[all.data.clean$corewater, ] %>%
  mutate(depth = paste("-",
                       strsplit(depth, "-") %>% sapply("[", 2),
                       sep = "") %>%
           as.numeric() / 100) %>%
  mutate(depth = depth + (waterDepth * corewater))
rm(waterDepth)


#### Save out data ####
saveRDS(all.data.clean,
        "dataEdited/waterChemistry/ICP/2020_WC_data.rds")
