#### code/cleaning_scripts/clean_exo_data.R ####
# Benjamin D. Peterson

# This scripts reads in the sonde data collected on the Exo
# sonde. The data was collected as part of the BLiMMP project.
# These samples were collected from Lake Mendota


#### This is your mess, you need to clean it up ####
rm(list = ls())
setwd("~/Documents/research/BLiMMP")
library(dplyr)
library(readxl)


#### 2019-08-21 ####
# Read in data
exo.data <- read_xlsx("dataRaw/exo/EXO_SD_FLAME_13E101468_082119_125251.xlsx",
                      skip = 22) %>%
  as.data.frame()

# Clean up datafile
# Read in custom headers that I made (found in sonde_profile_headers.csv)
sonde.names <- c("date", "time", "time_fraction", "site_name", "fault_code", "Battery_V",
                 "Cable_Pwr_V", "Chlorophyll_RFU", "Chlorophyll_µg.L", "BGA-PC_RFU", "BGA-PC_µg.L",
                 "ODO_sat", "ODO_conc", "fDOM_RFU", "fDOM_QSU", "Turbidity_FNU", "TSS_mg.L",
                 "Temp_C", "Cond_µS.cm", "SpCond_µS.cm", "Sal_psu", "nLF_Cond_µS.cm", "TDS_mg.L",
                 "Press_psi_a", "Depth_m")
# Add names, clean the data
colnames(exo.data) <- sonde.names
# Start at 1 meters
exo.data <- exo.data[which(exo.data$Depth_m >= 1), ]
# End at maximum depth
exo.data <- exo.data[1:which.max(exo.data$Depth_m), ]
# Order by depth and time
exo.data <- exo.data[order(exo.data$Depth_m, exo.data$time), ]
# Round depths to nearest 0.1m
exo.data$Depth_m <- round(exo.data$Depth_m,
                                       digits = 1)
# Keep only columns of interest
exo.data <- exo.data %>%
  select(Temp_C, Cond_µS.cm, SpCond_µS.cm, Sal_psu,
         ODO_sat, ODO_conc, fDOM_RFU,
         Turbidity_FNU, Chlorophyll_RFU, Depth_m)
# Aggregate data by depth
exo.data <- aggregate(exo.data,
                      by = list(exo.data$Depth_m),
                      FUN = mean)

# Write out csv of data
write.csv(exo.data,
          file = "dataEdited/exo/2019-08-21_profile.csv",
          row.names = FALSE)
rm(exo.data)

#### 2019-08-30 ####
# Read in data
exo.data.2 <- read_xlsx("dataRaw/exo/EXO_SD_FLAME_13E101468_083019_191108.xlsx",
                        skip = 23) %>%
  as.data.frame()

# Clean up datafile
# Read in custom headers that I made (found in sonde_profile_headers.csv)
sonde.names <- c("date", "time", "time_fraction", "site_name", "fault_code", "Battery_V",
                 "Cable_Pwr_V", "Temp_C", "Cond_µS.cm", "SpCond_µS.cm", "Sal_psu",
                 "nLF_Cond_µS.cm", "TDS_mg.L", "ODO_sat", "ODO_conc", "fDOM_RFU", "fDOM_QSU",
                 "pH", "pH_mV", "ORP_mV", "ORP_raw_mV", "Turbidity_FNU", "TSS_mg.L",
                 "Chlorophyll_RFU", "Chlorophyll_µg.L", "BGA-PC_RFU", "BGA-PC_µg.L",
                 "Press_psi_a", "Depth_m")
# Add names, clean the data
colnames(exo.data.2) <- sonde.names
# Start at 1 meters
exo.data.2 <- exo.data.2[which(exo.data.2$Depth_m >= 1), ]
# End at maximum depth
exo.data.2 <- exo.data.2[1:which.max(exo.data.2$Depth_m), ]
# Order by depth and time
exo.data.2 <- exo.data.2[order(exo.data.2$Depth_m, exo.data.2$time), ]
# Round depths to nearest 0.1m and average these measurements
exo.data.2$Depth_m <- round(exo.data.2$Depth_m,
                            digits = 1)
exo.data.2 <- exo.data.2 %>%
  select(Temp_C, Cond_µS.cm, SpCond_µS.cm, Sal_psu,
         ODO_sat, ODO_conc, fDOM_RFU, pH, ORP_mV,
         Turbidity_FNU, Chlorophyll_RFU, Depth_m)
exo.data.2 <- aggregate(exo.data.2,
                        by = list(exo.data.2$Depth_m),
                        FUN = mean)

# Write out csv of data
write.csv(exo.data.2,
          file = "dataEdited/exo/2019-08-30_profile.csv",
          row.names = FALSE)

#### 2019-10-08 ####
# Read in data
exo.data.3 <- read_xlsx("dataRaw/exo/EXO_SD_FLAME_13E101468_100819_183725.xlsx",
                        skip = 24) %>%
  as.data.frame()

# Clean up datafile
# Read in custom headers that I made (found in sonde_profile_headers.csv)
sonde.names <- c("date", "time", "time_fraction", "site_name", "fault_code", "Battery_V",
                 "Cable_Pwr_V", "Temp_C", "Cond_µS.cm", "SpCond_µS.cm", "Sal_psu",
                 "nLF_Cond_µS.cm", "TDS_mg.L", "ODO_sat", "ODO_conc", "fDOM_RFU", "fDOM_QSU",
                 "pH", "pH_mV", "ORP_mV", "ORP_raw_mV", "Turbidity_FNU", "TSS_mg.L",
                 "Chlorophyll_RFU", "Chlorophyll_µg.L", "BGA-PC_RFU", "BGA-PC_µg.L",
                 "Press_psi_a", "Depth_m")
# Add names, clean the data
colnames(exo.data.3) <- sonde.names
# Start at 1 meters
exo.data.3 <- exo.data.3[which(exo.data.3$Depth_m >= 1), ]
# End at maximum depth
exo.data.3 <- exo.data.3[1:which.max(exo.data.3$Depth_m), ]
# Order by depth and time
exo.data.3 <- exo.data.3[order(exo.data.3$Depth_m, exo.data.3$time), ]
# Round depths to nearest 0.1m and average these measurements
exo.data.3$Depth_m <- round(exo.data.3$Depth_m,
                            digits = 1)
exo.data.3 <- exo.data.3 %>%
  select(Temp_C, Cond_µS.cm, SpCond_µS.cm, Sal_psu,
         ODO_sat, ODO_conc, fDOM_RFU, pH, ORP_mV,
         Turbidity_FNU, Chlorophyll_RFU, Depth_m)
exo.data.3 <- aggregate(exo.data.3,
                        by = list(exo.data.3$Depth_m),
                        FUN = mean) %>%
  select(-Group.1)

# Write out csv of data
write.csv(exo.data.3,
          file = "dataEdited/exo/20191008_profile.csv",
          row.names = FALSE)












#### 2020 data ####
rm(list = ls())
source("code/cleaning_scripts/cleaning_functions.R")

clean.exo.data(file.name = "exo_profile_ME20200807a.xlsx",
               output.name = "2020-08-07_profile.csv")
clean.exo.data(file.name = "exo_profile_ME20200814.xlsx",
               output.name = "2020-08-14_profile.csv")
clean.exo.data(file.name = "exo_profile_ME20200902.xlsx",
               output.name = "2020-09-02_profile.csv",
               depth.adjustment = 1)
clean.exo.data(file.name = "exo_profile_ME20200917.xlsx",
               output.name = "2020-09-17_profile.csv")
clean.exo.data(file.name = "exo_profile_ME20200924.xlsx",
               output.name = "2020-09-24_profile.csv")
clean.exo.data(file.name = "exo_profile_ME20201010_profile.xlsx",
               output.name = "2020-10-10_profile.csv",
               pH.sensor.present = FALSE,
               depth.adjustment = 0.6)




#### 2021 data ####
rm(list = ls())
source("code/cleaning_scripts/cleaning_functions.R")

custom.header.names.vector = c("date", "time", "time_fraction", "site_name", "fault_code", "Battery_V",
                               "Cable_Pwr_V", "ODO_sat", "ODO_conc", "fDOM_RFU", "fDOM_QSU", "Temp_C",
                               "Cond_µS.cm", "SpCond_µS.cm", "Sal_psu", "nLF_Cond_µS.cm", "TDS_mg.L", 
                               "Turbidity_FNU", "TSS_mg.L", "Chlorophyll_RFU", "Chlorophyll_µg.L",
                               "BGA-PC_RFU", "BGA-PC_µg.L", "Press_psi_a", "depth")
clean.exo.data(file.name = "EXO_SD_FLAME_13E101468_091021_150339.xlsx",
               output.name = "2021-09-10_profile.csv",
               custom.header.names = custom.header.names.vector,
               pH.sensor.present = FALSE)
