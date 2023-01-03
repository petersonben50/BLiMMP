#### code/cleaning_scripts/aggregate_water_chem_data.R ####
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("~/Documents/research/BLiMMP/")
library(tidyverse)


#### Set detection flags ####
dl.vector <- c("NONE", "<")
names(dl.vector) <- c(FALSE, TRUE)


#### Read in sulfide data ####
sulfide.data <- read.csv("dataEdited/waterChemistry/sulfide/WC_data.csv") %>%
  rename(containerID = sulfurID,
         date = startDate,
         concentration = S_conc_uM) %>%
  mutate(constituent = "sulfide_uM") %>%
  select(date, depth, corewater, replicate, constituent, concentration) %>%
  # Detection limit ~1 µM
  mutate(detection_flag = unname(dl.vector[as.character(concentration < 1.0)])) %>%
  filter(year(date) != 2019)


#### Read in sulfate data ####
sulfate.data <- read.csv("dataEdited/waterChemistry/sulfate/WC_data.csv") %>%
  select(-concentration) %>%
  rename(date = startDate,
         concentration = sulfate_uM) %>%
  mutate(constituent = "sulfate_uM") %>%
  select(date, depth, corewater, replicate, constituent, concentration) %>%
  # Detection limit ~1 µM
  mutate(detection_flag = unname(dl.vector[as.character(concentration < 1.0)]))


#### Read in metals data ####
metals.data <- readRDS("dataEdited/waterChemistry/ICP/2020_2021_Fe_Mn_data.rds") %>%
  mutate(corewater = FALSE,
         constituent = paste(substr(constituent, 1, 2),
                             substr(state, 1, 4),
                             substr(constituent, 4, 6),
                             sep = "_"),
         concentration = round(concentration, 3)) %>%
  select(date, depth, corewater, replicate, constituent, concentration) %>%
  # Detection limit ~0.02 ppm
  mutate(detection_flag = unname(dl.vector[as.character(concentration < 0.02)])) %>%
  as.data.frame()


#### Read in Hg data ####
Hg.data <- read.csv("dataEdited/Hg/Hg_data_clean.csv") %>%
  rename(date = sampleDate) %>%
  select(date, depth, corewater, replicate, constituent, concentration, detection_flag)
  

#### Read in DGM data ####
DGM.data <- readRDS("dataEdited/Hg/DGM_2021_data.rds") %>%
  rename(date = startDate,
         concentration = concentration_ng.L) %>%
  mutate(corewater = as.logical("FALSE"),
         constituent = "DGM_NG.L",
         detection_flag = as.logical("NA")) %>%
  select(date, depth, corewater, replicate, constituent, concentration, detection_flag) %>%
  mutate(concentration = round(concentration, 3))


#### Combine data ####
all.data <- metals.data %>%
  rbind(sulfate.data) %>%
  rbind(sulfide.data) %>%
  rbind(Hg.data) %>%
  rbind(DGM.data) %>%
  ungroup() %>%
  filter(!is.na(concentration))


#### Anything under 0, set to 0 ####
all.data$concentration[all.data$concentration < 0] <- 0


#### Make data wide ####
all.data.wide <- all.data %>%
  group_by(date, depth, corewater, constituent) %>%
  mutate(replicate = row_number()) %>%
  gather(key = type,
         value = assigned.value,
         c(concentration, detection_flag)) %>%
  mutate(type = gsub("concentration", "",
                     type),
         type = gsub("detection_flag", "_df",
                     type),
         constituent = paste(constituent, type,
                             sep = "")) %>%
  select(-type) %>%
  spread(key = constituent,
         value = assigned.value) %>%
  filter(!(date %in% c("2020-09-17",
                       "2021-08-31")))


#### Save data ####
write.csv(all.data.wide,
          file = "dataFinal/water_chem_data.csv",
          row.names = FALSE)
