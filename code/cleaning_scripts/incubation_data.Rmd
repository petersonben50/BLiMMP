#### code/cleaning_scripts/incubation_data.R ####
# Written for BLiMMP project
# Benjamin D. Peterson



#### Prep workspace ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP")
library(lubridate)
library(readxl)
library(tidyverse)



#### Read in data ####
incubation_metadata <- read.csv("metadata/processedMetadata/incubation_metadata.csv")
MeHg_data_2020 <- read_xlsx("dataRaw/incubations/MeHg/2020_incubations/2020_BENDOTA_aggregated_data.xlsx",
                            sheet = "keeper-data",
                            skip = 1,
                            col_names = c("bottleID", "MeHg_ambient_ppt", "MeHg_198_ppt", "MeHg_204_ppt", "MeHg_amb_DDL", "MeHg_excess_DDL", "MeHg_DOA"),
                            col_types = c("text", rep("numeric", 5), "text", rep("skip", 6))) %>%
  mutate(MeHg_DOA = as.Date(MeHg_DOA))
MeHg_data_2021 <- read_xlsx("dataRaw/incubations/MeHg/2021_incubations/2021_BENDOTA_aggregated_data.xlsx",
                            sheet = "keeper-data",
                            skip = 1,
                            col_names = c("bottleID", "MeHg_ambient_ppt", "MeHg_198_ppt", "MeHg_204_ppt", "MeHg_amb_DDL", "MeHg_excess_DDL", "MeHg_DOA"),
                            col_types = c("text", rep("numeric", 5), "text", rep("skip", 6))) %>%
  mutate(MeHg_DOA = as.Date(MeHg_DOA))
MeHg_data <- rbind(MeHg_data_2020,
                   MeHg_data_2021)
rm(MeHg_data_2020, MeHg_data_2021)
HgT_data <- read_xlsx("dataRaw/incubations/HgT/HgT_BENDOTA_SUMMARY.xlsx",
                      col_names = c("bottleID", "HgT_ambient_ppt", "HgT_198_ppt", "HgT_204_ppt", "HgT_excess_DDL", "HgT_DOA"),
                      skip = 1)



#### Prepare MeHg data ####
MeHg_data_clean <- right_join(incubation_metadata %>%
                                filter(constituent == "MeHg"),
                              MeHg_data) %>%
  # Identify which samples measured Me198Hg above detection
  mutate(MeHg_198_above_DDL = (MeHg_198_ppt >= MeHg_excess_DDL)) %>%
  # Round off measurements to 3 decimal places
  mutate(MeHg_ambient_ppt = round(MeHg_ambient_ppt, 3),
         MeHg_198_ppt = round(MeHg_198_ppt, 3),
         MeHg_204_ppt = round(MeHg_204_ppt, 3)) %>%
  mutate(bottleID_MeHg = bottleID) %>%
  select(-c(bottleID, constituent))
rm(MeHg_data)



#### Prepare HgT data ####
HgT_data_clean <- right_join(incubation_metadata %>%
                               filter(constituent == "HgT"),
                             HgT_data) %>%
  mutate(bottleID_HgT = bottleID) %>%
  select(-c(bottleID, constituent)) %>%
  filter(year(startDate) %in% c(2020, 2021))
rm(HgT_data, incubation_metadata)



#### Combine data ####
Hg_incubation_data <- full_join(MeHg_data_clean,
                                HgT_data_clean)
rm(MeHg_data_clean,
   HgT_data_clean)



#### Calculate percent MeHg for ambient, 198, and 204
Hg_incubation_data <- Hg_incubation_data %>%
  mutate(percent_amb_MeHg = MeHg_ambient_ppt / HgT_ambient_ppt * 100,
         percent_198_MeHg = MeHg_198_ppt / HgT_198_ppt * 100,
         percent_204_MeHg = MeHg_204_ppt / HgT_204_ppt * 100)



#### Read out all data ####
write.csv(Hg_incubation_data,
          file = "dataFinal/incubation_Hg_conc_data.csv",
          row.names = FALSE)



#### Calculate Kmet and 198HgT loss ####
kmet_data <- Hg_incubation_data %>%
  select(sampleID, incubationID,
         startDate, depth,
         t, treatment, durationInDays,
         MeHg_198_ppt, HgT_198_ppt) %>%
  filter(treatment %in% c("unfiltered-unamended",
                          "unfiltered-molybdate",
                          "filtered-unamended")) %>%
  gather(key = constituent, value = concentration, c(MeHg_198_ppt, HgT_198_ppt, durationInDays)) %>%
  mutate(const_time = paste(constituent, t,
                            sep = "_")) %>%
  select(sampleID, incubationID, startDate, depth,
         treatment,const_time, concentration) %>%
  spread(key = const_time,
         value = concentration) %>%
  mutate(Kmet_t1 = (MeHg_198_ppt_t1 - MeHg_198_ppt_t0) / HgT_198_ppt_t1 / durationInDays_t1,
         Kmet_t2 = (MeHg_198_ppt_t2 - MeHg_198_ppt_t1) / HgT_198_ppt_t2 / (durationInDays_t2 - durationInDays_t1),
         Kmet_total = (MeHg_198_ppt_t2 - MeHg_198_ppt_t0) / ((HgT_198_ppt_t2 + HgT_198_ppt_t1) / 2) / (durationInDays_t2 - durationInDays_t0),
         HgT_198_daily_percent_loss_t1 = (HgT_198_ppt_t0 - HgT_198_ppt_t1) / HgT_198_ppt_t0 / durationInDays_t1 * 100,
         HgT_198_daily_percent_loss_t2 = (HgT_198_ppt_t1 - HgT_198_ppt_t2) / HgT_198_ppt_t1 / (durationInDays_t2 - durationInDays_t1) * 100,
         HgT_198_percent_loss_total = (HgT_198_ppt_t0 - HgT_198_ppt_t2) / HgT_198_ppt_t0 * 100) %>%
  select(sampleID, incubationID, startDate, depth, treatment,
         Kmet_t1, Kmet_t2, Kmet_total,
         HgT_198_daily_percent_loss_t1,
         HgT_198_daily_percent_loss_t2,
         HgT_198_percent_loss_total)



#### Calculate Kdem and 204HgT loss ####
kdem_data <- Hg_incubation_data %>%
  select(sampleID, incubationID,
         startDate, depth,
         t, treatment, durationInDays,
         MeHg_204_ppt, HgT_204_ppt) %>%
  filter(treatment %in% c("unfiltered-unamended",
                          "unfiltered-molybdate",
                          "filtered-unamended")) %>%
  gather(key = constituent, value = concentration, c(MeHg_204_ppt, HgT_204_ppt, durationInDays)) %>%
  mutate(const_time = paste(constituent, t,
                            sep = "_")) %>%
  select(sampleID, incubationID, startDate, depth,
         treatment,const_time, concentration) %>%
  spread(key = const_time,
         value = concentration) %>%
  mutate(Kdem_t1 = -((MeHg_204_ppt_t1 / HgT_204_ppt_t1) - (MeHg_204_ppt_t0 / HgT_204_ppt_t0)) / durationInDays_t1,
         Kdem_t2 = -((MeHg_204_ppt_t2 / HgT_204_ppt_t2) - (MeHg_204_ppt_t1 / HgT_204_ppt_t1)) / (durationInDays_t2 - durationInDays_t1),
         Kdem_total = -((MeHg_204_ppt_t2 / HgT_204_ppt_t2) - (MeHg_204_ppt_t0 / HgT_204_ppt_t0)) / (durationInDays_t2 - durationInDays_t0),
         HgT_204_daily_percent_loss_t1 = (HgT_204_ppt_t0 - HgT_204_ppt_t1) / HgT_204_ppt_t0 / durationInDays_t1 * 100,
         HgT_204_daily_percent_loss_t2 = (HgT_204_ppt_t1 - HgT_204_ppt_t2) / HgT_204_ppt_t1 / (durationInDays_t2 - durationInDays_t1) * 100,
         HgT_204_percent_loss_total = (HgT_204_ppt_t0 - HgT_204_ppt_t2) / HgT_204_ppt_t0 * 100) %>%
  select(sampleID, incubationID, startDate, depth, treatment,
         Kdem_t1, Kdem_t2, Kdem_total,
         HgT_204_daily_percent_loss_t1,
         HgT_204_daily_percent_loss_t2,
         HgT_204_percent_loss_total)



#### Combine all rate data ####
Hg_rate_data <- full_join(kmet_data,
                          kdem_data)



#### Read out Kmet data ####
write.csv(Hg_rate_data,
          file = "dataFinal/incubation_Hg_rate_data.csv",
          row.names = FALSE)

