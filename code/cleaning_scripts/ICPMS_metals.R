#### code/cleaning_scripts/ICPMS_metals.R ####
# Benjamin D. Peterson


#### Set it up ####
rm(list = ls())
setwd("~/Documents/research/BLiMMP/")
library(readxl)
library(tidyverse)


#### Read in metadata ####
metadata <- read.csv("metadata/processedMetadata/metals_WC.csv") %>%
  select(metalID, sampleID, tripID, startDate, filteredInField,
         depth, tare, mass, replicate, preservativeVol)


#### 2021 data: Read in ICP-MS data ####
icpms.data <- rbind(read_xlsx("dataEdited/waterChemistry/ICP/2021_ICPMS_data.xlsx",
                              sheet = "2022-02-10_conc_edited"),
                    read_xlsx("dataEdited/waterChemistry/ICP/2021_ICPMS_data.xlsx",
                              sheet = "2022-02-11_conc_edited")) %>%
  rename(metalID = `Sample Name`) %>%
  filter(grepl("BLI21", metalID)) %>%
  select(metalID, Mn_ppb, Fe_ppb)


#### 2021 data: Add metadata to data ####
all.data <- inner_join(metadata,
                       icpms.data) %>%
  select(sampleID, startDate, depth, filteredInField,
         tare, mass, preservativeVol, replicate,
         Mn_ppb, Fe_ppb) %>%
  mutate(tare = 6.75)



#### 2021 data: Adjust concentrations for dilution with preservative ####
all.data <- all.data %>%
  gather(key = constituent,
         value = concentration,
         c(Mn_ppb, Fe_ppb)) %>%
  mutate(concentration = as.numeric(gsub("<", 0, concentration))) %>%
  rename(concentration_diluted = concentration) %>%
  mutate(concentration = (concentration_diluted * (mass - tare)) / (mass - tare - preservativeVol)) %>%
  select(-c(tare, mass, preservativeVol))


#### 2021 data: Check out duplicates ####
dup.samples <- all.data %>%
  filter(replicate == "yes") %>%
  select(sampleID) %>%
  unlist(use.names = FALSE) %>%
  unique()
all.data %>%
  filter(sampleID %in% dup.samples) %>%
  arrange(constituent, filteredInField, startDate, depth)
# These are all pretty close. We'll just keep rep #1
all.data <- all.data %>%
  group_by(startDate, depth, filteredInField, constituent) %>%
  summarize(concentration = mean(concentration))


#### 2021 data: Calculate particulate Mn ####
all.data.clean <- all.data %>%
  spread(key = filteredInField,
         value = concentration) %>%
  mutate(particulate = no - yes) %>%
  rename(dissolved = yes) %>%
  select(-no) %>%
  gather(key = state,
         value = concentration,
         c(dissolved, particulate)) %>%
  rename(date = startDate) %>%
  mutate(depth = as.numeric(depth))
rm(all.data, icpms.data)


#### 2020 data: Read in data ####
icpms.wc.data <- read_xlsx("dataEdited/waterChemistry/ICP/unprocessed/Bendota_ICPMS_results043021.xlsx",
                           sheet = "metals_WC") %>%
  filter(!grepl(" dup", metalID)) %>%
  rename(Mn_ppb = `[Mn] ppb`,
         Fe_ppb = `[Fe] ppb`) %>%
  select(metalID, Mn_ppb, Fe_ppb)


#### 2020 data: Add metadata ####
all.data.2020 <- inner_join(metadata,
                            icpms.wc.data) 
rm(icpms.wc.data)


#### 2020 data: Adjust concentrations for dilution with preservative ####
all.data.2020 <- all.data.2020 %>%
  select(-metalID) %>%
  gather(key = constituent,
         value = concentration,
         c(Mn_ppb, Fe_ppb)) %>%
  rename(concentration_diluted = concentration) %>%
  mutate(concentration = (concentration_diluted * (mass - tare)) / (mass - tare - preservativeVol)) %>%
  select(-c(tare, mass, preservativeVol))


#### 2020 data: Check out duplicates ####
dup.samples <- all.data.2020 %>%
  filter(replicate == "yes") %>%
  select(sampleID) %>%
  unlist(use.names = FALSE) %>%
  unique()
all.data.2020 %>%
  filter(sampleID %in% dup.samples) %>%
  arrange(constituent, filteredInField, startDate, depth)
# Let's average these values
all.data.2020 <- all.data.2020 %>%
  group_by(startDate, depth, filteredInField, constituent) %>%
  summarize(concentration = mean(concentration))



#### 2020 data: Calculate particulate Mn ####
all.data.2020.clean <- all.data.2020 %>%
  spread(key = filteredInField,
         value = concentration) %>%
  filter(!grepl("-", depth)) %>%
  mutate(particulate = no - yes) %>%
  rename(dissolved = yes) %>%
  select(-no) %>%
  gather(key = state,
         value = concentration,
         c(dissolved, particulate)) %>%
  rename(date = startDate) %>%
  mutate(depth = as.numeric(depth))


#### Combine both years of data ####
all.data.both.years <- rbind(all.data.2020.clean,
                             all.data.clean)

#### Save out data ####
saveRDS(all.data.both.years,
        "dataEdited/waterChemistry/ICP/2020_2021_Fe_Mn_data.rds")
