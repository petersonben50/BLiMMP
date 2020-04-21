#### code/cleaning_scripts/clean_leucine_uptake.R ####
# Benjamin D. Peterson

#### Your mom isn't here, and even if she was she shouldn't have to clean up after you ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP")
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
  select(uptakeID, incubationID, sampleID, startDate, depth, condition)

rm(list = ls(pattern = ".metadata"),
   filtered.vector)


#### Calculate conversion factor ####

efficiency.data <- read_xlsx("dataRaw/leucineUptake/20191009_leucine.xlsx") %>%
  filter(Condition == "efficiency")

# We used 1µl of a 2µM leucine solution for our efficiency test.
# Therefore, we'll divide the number of counts that we got to 
# calculate the overall efficiency of the counts, which in this
# case means the counts that we get per leucine in solution.
# This will provide a direct link between the counts we got and
# the amount of leucine taken up.

# 2µM leucine solution -> 2pmol/ul
# 1µl volume
leucine.pmol <- 2*1

efficiency.data.count.per.pmol <- efficiency.data %>%
  mutate(CF_A = CPMA / leucine.pmol) %>%
  mutate(CF_B = CPMB / leucine.pmol) %>%
  select(CF_A, CF_B) %>%
  colMeans()

rm(efficiency.data, leucine.pmol)



#### Read in and clean data: 2019-10-14 ####

# I'll walk through the first one closely with notes

# First we set the time that each of these incubated for.
incubation.time.hours.temp <- 90/60

# We also need the total volume of the incubation in L.
incubation.volume.mL.temp <- 1
incubation.volume.L.temp <- incubation.volume.mL.temp / 1000
rm(incubation.volume.mL.temp)

# I'm assuming that the ID of leucine is 2, as described by Simon and Azam, 1989.
leucine.ID <- 2

# Set the molar ratio of leucine in protein as 7.3%, according to Exercise 19 in
# Limnological Analyses by Wetzel and Likens. This ratio was determined in Simon
# and Azam, 1989
Leu.per.prot.molar.ratio <- 0.073

# Set molar mass leucine
molar.mass.leucine <- 131.2

# Set cellular carbon per protein
cell_C_to_protein <- 0.86

# Then we'll read in the raw data
leucine.data.temp <- read_xlsx("dataRaw/leucineUptake/20191014_leucine.xlsx") 

# Calculate the blank values
blank.data.temp <- leucine.data.temp %>%
  filter(uptakeID == "negative_control") %>%
  select(CPMA, CPMB) %>%
  colMeans()

# Now for the calculations.
# For both A and B counting channels, we'll first take the counts and subtract out the blank values.
# Then we'll divide the counts by the efficiency data, to calculate the pmol of leucine that have been
# taken up. Then we'll divide by the incubation volume to get pmol/L. And then we'll divide it by
# the number of hours of incubation, to normalize production to 1 hour.
# This gives us the pM of leucine incorporated per hour.

incubation.data.20191014 <- leucine.data.temp %>%
  mutate(Leu_uptake_pM_A = (CPMA - blank.data.temp["CPMA"]) / (efficiency.data.count.per.pmol["CF_A"] * incubation.volume.L.temp * incubation.time.hours.temp)) %>%
  mutate(Leu_uptake_pM_B = (CPMB - blank.data.temp["CPMB"]) / (efficiency.data.count.per.pmol["CF_B"] * incubation.volume.L.temp * incubation.time.hours.temp)) %>%
  #mutate(ratio_A_B = Leu_uptake_pM_A / Leu_uptake_pM_B) %>%
  
  # To convert this to overall pM of protein produced per house, we divide the pM of leucine by the
  # molar ratio of leucine in protein, which is approximately 7.3% and multiply by the molar mass of 
  # leucine and the leucine isotope dilution, which here we are assuming is 2.
  # I think this is a BS way of doing it, since the 7.3% value is a molar ratio, not a mass ratio.
  # But, this is how people have been calculating it forever, so I'm gonna stick with it. At some point 
  # I should use Simon and Azam's tables to calculate the mass ratio.
  mutate(µgBPP_per_L_hr = round(((Leu_uptake_pM_A / Leu.per.prot.molar.ratio) * molar.mass.leucine * leucine.ID) / 1000000, 3)) %>%
  
  # The last conversion to make is from protein biomass to overall carbon biomass.
  mutate(µgBCP_per_L_hr = µgBPP_per_L_hr * cell_C_to_protein) %>%
  filter(Condition == "ambient") %>%
  select(c(uptakeID, Leu_uptake_pM_A, µgBPP_per_L_hr, µgBCP_per_L_hr))

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
  mutate(Leu_uptake_pM_A = (CPMA - blank_CPMA) / (efficiency.data.count.per.pmol["CF_A"] * incubation.volume.L.temp * incubation.time.hours.temp)) %>%
  mutate(Leu_uptake_pM_B = (CPMB - blank_CPMB) / (efficiency.data.count.per.pmol["CF_B"] * incubation.volume.L.temp * incubation.time.hours.temp)) %>%
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
  mutate(Leu_uptake_pM_A = (CPMA - blank_CPMA) / (efficiency.data.count.per.pmol["CF_A"] * incubation.volume.L.temp * incubation.time.hours.temp)) %>%
  mutate(Leu_uptake_pM_B = (CPMB - blank_CPMB) / (efficiency.data.count.per.pmol["CF_B"] * incubation.volume.L.temp * incubation.time.hours.temp)) %>%
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
