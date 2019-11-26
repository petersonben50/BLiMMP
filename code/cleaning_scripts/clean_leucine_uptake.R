#### code/cleaning_scripts/clean_leucine_uptake.R ####
# Benjamin D. Peterson

#### Your mom isn't here, clean this shit up ####
rm(list = ls())
setwd("~/Box/BLiMMP/")
library(readxl)
library(tidyverse)


#### Read in metadata ####

# Water column data
leucine.metadata <- read.csv("metadata/LEU.csv",
                             stringsAsFactors = FALSE)
sampleID.metadata <- read_xlsx("metadata/2_sample_IDs.xlsx")
tripID.metadata <- read_xlsx("metadata/1_trip_IDs.xlsx")

metadata.WC <- leucine.metadata %>%
  left_join(sampleID.metadata) %>%
  left_join(tripID.metadata) %>%
  select(uptakeID, sampleID, depth, startDate)

# Incubation data
incubationID.metadata <- read_xlsx("metadata/4_MA_ID.xlsx")

filtered.vector <- c("filtered", "unfiltered")
names(filtered.vector) <- c("yes", "no")

metadata.MA <- leucine.metadata %>%
  filter(!is.na(incubationID)) %>%
  select(-sampleID) %>%
  left_join(incubationID.metadata) %>%
  left_join(sampleID.metadata) %>%
  left_join(tripID.metadata,
            by = "tripID") %>%
  mutate(condition = paste(amendment,
                           " - ",
                           filtered.vector[filtered],
                           sep = "")) %>%
  select(uptakeID, incubationID, startDate, depth, condition)

rm(list = ls(pattern = ".metadata"),
   filtered.vector)


#### Calculate conversion factor ####

efficiency.data <- read_xlsx("dataRaw/leucineUptake/20191009_leucine.xlsx") %>%
  filter(Condition == "efficiency")

# 2µM leucine solution -> 2pmol/ul
# 1µl volume
leucine.pmol.eff <- 2*1

efficiency.data <- efficiency.data %>%
  mutate(CF_A = CPMA / leucine.pmol.eff) %>%
  mutate(CF_B = CPMB / leucine.pmol.eff) %>%
  select(CF_A, CF_B) %>%
  colMeans()

rm(leucine.pmol.eff)

#### Read in and clean data: 2019-10-14 ####

incubation.time.hours.temp <- 90/60
incubation.volume.L.temp <- 0.001

leucine.data.temp <- read_xlsx("dataRaw/leucineUptake/20191014_leucine.xlsx") 

blank.data.temp <- leucine.data.temp %>%
  filter(uptakeID == "negative_control") %>%
  select(CPMA, CPMB) %>%
  colMeans()
  
incubation.data.20191014 <- leucine.data.temp %>%
  mutate(Leu_uptake_pM_A = (CPMA - blank.data.temp["CPMA"]) / (incubation.time.hours.temp * incubation.volume.L.temp * efficiency.data["CF_A"])) %>%
  mutate(Leu_uptake_pM_B = (CPMB - blank.data.temp["CPMB"]) / (incubation.time.hours.temp * incubation.volume.L.temp * efficiency.data["CF_B"])) %>%
  mutate(ratio_A_B = Leu_uptake_pM_A / Leu_uptake_pM_B) %>%
  filter(Condition == "ambient") %>%
  select(c(uptakeID, Leu_uptake_pM_A))

rm(list = ls(pattern = ".temp"))

# Combine data with metadata

incubation.data.20191014 <- incubation.data.20191014 %>%
  left_join(metadata.MA)





#### Read in and clean data: 2019-10-09 ####

incubation.time.hours.temp <- 1
incubation.volume.L.temp <- 0.001

leucine.data.temp <- read_xlsx("dataRaw/leucineUptake/20191009_leucine.xlsx") %>%
  left_join(metadata.WC) %>%
  as.data.frame()


# Generate blank data
blank.data.temp <- leucine.data.temp %>%
  filter(Condition == "killed") %>%
  select(CPMA, CPMB, sampleID) %>%
  group_by(sampleID) %>%
  summarise_all(mean) %>%
  rename(blank_CPMA = CPMA,
         blank_CPMB = CPMB)



WC.data.20191009 <- leucine.data.temp %>%
  filter(Condition == "ambient") %>%
  left_join(blank.data.temp) %>%
  mutate(Leu_uptake_pM_A = (CPMA - blank_CPMA) / (incubation.time.hours.temp * incubation.volume.L.temp * efficiency.data["CF_A"])) %>%
  mutate(Leu_uptake_pM_B = (CPMB - blank_CPMB) / (incubation.time.hours.temp * incubation.volume.L.temp * efficiency.data["CF_B"])) %>%
  mutate(ratio_A_B = Leu_uptake_pM_A / Leu_uptake_pM_B) %>%
  select(c(uptakeID, sampleID, startDate, depth, Leu_uptake_pM_A))


rm(list = ls(pattern = ".temp"))

# Combine data with metadata

WC.data.20191009 <- WC.data.20191009 %>%
  left_join(metadata.WC)





#### Read in and clean data: 2019-08-31 ####

incubation.time.hours.temp <- 1
incubation.volume.L.temp <- 0.001

leucine.data.temp <- read_xlsx("dataRaw/leucineUptake/20190831_leucine.xlsx") %>%
  left_join(metadata.WC) %>%
  as.data.frame()

# Generate blank data
blank.data.temp <- leucine.data.temp %>%
  filter(Condition == "killed") %>%
  select(CPMA, CPMB, sampleID) %>%
  group_by(sampleID) %>%
  summarise_all(mean) %>%
  rename(blank_CPMA = CPMA,
         blank_CPMB = CPMB)


WC.data.20190831 <- leucine.data.temp %>%
  filter(Condition == "ambient") %>%
  left_join(blank.data.temp) %>%
  mutate(Leu_uptake_pM_A = (CPMA - blank_CPMA) / (incubation.time.hours.temp * incubation.volume.L.temp * efficiency.data["CF_A"])) %>%
  mutate(Leu_uptake_pM_B = (CPMB - blank_CPMB) / (incubation.time.hours.temp * incubation.volume.L.temp * efficiency.data["CF_B"])) %>%
  mutate(ratio_A_B = Leu_uptake_pM_A / Leu_uptake_pM_B) %>%
  select(c(uptakeID, sampleID, startDate, depth, Leu_uptake_pM_A))

rm(list = ls(pattern = ".temp"))





#### Write out water colum data ####

WC.data <- rbind(WC.data.20190831,
                 WC.data.20191009)

write.csv(WC.data,
          "dataEdited/leucineUptake/water_column_uptake.csv",
          quote = FALSE,
          row.names = FALSE)






#### Write out incubation data ####

incubation.data <- rbind(incubation.data.20191014)

write.csv(incubation.data,
          "dataEdited/leucineUptake/incubation_uptake.csv",
          quote = FALSE,
          row.names = FALSE)
