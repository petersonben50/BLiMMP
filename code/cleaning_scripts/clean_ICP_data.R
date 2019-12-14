#### code/cleaning_scripts/clean_ICP_data.R ####
# Benjamin D. Peterson

# This script will clean up the data generated from
# the ICP in water sciences. Most of the QC and 
# filtering of samples is done within the program,
# so most of the work here is just adjusting for
# the sample dilution during collection and 
# preservation, and saving out the data. 

# I will have inspected the quality of the calibration
# curve, which is not passed if not greater than 0.999.
# I will also have manually inspected the RSD values for
# all the samples, and will have removed any bad ones.
# Continuing calibration blanks and checks will be included 
# in the run, and the run discarded if they do not match.
# We'll also manually inspect method blanks, which need
# to be below detection.

#### All the world's a stage, but it needs to be cleaned sometimes ####
rm(list = ls())
setwd("~/Box/BLiMMP/")
library(readxl)
library(tidyverse)


#### Prepare metadata ####
# 
# UM.metadata <- read_xlsx("metadata/chem_UM.xlsx") %>%
#   filter(!is.na(sampleID)) %>%
#   select(MUID, sampleID, tare, mass, preservativeVol) %>%
#   rename(metalID = MUID)
# FM.metadata <- read_xlsx("metadata/chem_FM.xlsx") %>%
#   filter(!is.na(sampleID)) %>%
#   select(MFID, sampleID, tare, mass, preservativeVol) %>%
#   rename(metalID = MFID)
# 
# sample_IDs <- read_xlsx("metadata/2_sample_IDs.xlsx")
# trip_IDs <- read_xlsx("metadata/1_trip_IDs.xlsx") %>%
#   select(-notes)
# 
# # Join all metadata
# metal.metadata <- rbind(UM.metadata,
#                         FM.metadata) %>%
#   left_join(sample_IDs) %>%
#   left_join(trip_IDs) %>%
#   select(metalID, sampleID, tripID, depth, startDate, tare, mass, preservativeVol)
# 
# # Save out metadata
# write.csv(metal.metadata,
#           file = "metadata/processedMetadata/metal_metadata.csv",
#           row.names = FALSE,
#           quote = FALSE)
# 
# rm(FM.metadata,
#    UM.metadata,
#    sample_IDs,
#    trip_IDs)



#### Read in data ####

data.file.name <- "dataEdited/waterChemistry/ICP/unprocessed/ICP_20191213.csv"


ICP.data <- read.csv(data.file.name,
                     skip = 2,
                     stringsAsFactors = FALSE) %>%
  select(Label, Element.Label, Concentration, Unit) %>%
  rename(metalID = Label,
         element = Element.Label,
         conc_ppm = Concentration)



#### Establish unit of measurement ####

unit.of.measurement <- ICP.data$Unit[1]
ICP.data <- ICP.data %>%
  select(-Unit)

#### Check that technical replicates are within 15% ####

rep.data <- ICP.data[grep("rep", ICP.data$metalID), ]
rep.data$metalID <- strsplit(rep.data$metalID, "_rep") %>%
  sapply("[", 1)
rep.data <- rep.data %>%
  group_by(metalID, element) %>%
  summarise(replicate.mean = mean(conc_ppm),
            replicate.RSD = sd(conc_ppm) / mean(conc_ppm) * 100)
# Replicates look good. Now average them
ICP.data <- ICP.data %>%
  mutate(metalID = strsplit(metalID, "_rep") %>%
           sapply("[", 1)) %>%
  group_by(metalID, element) %>%
  summarise(conc_ppm = mean(conc_ppm))


#### Join with metadata ####

# Remove standards
all.data <- ICP.data %>%
  inner_join(metal.metadata) %>%
  mutate(corr_conc_ppm = round(((mass - tare) / (mass - tare - preservativeVol)) * conc_ppm,
                               2)) %>%
  select(metalID, element,sampleID, tripID, depth, startDate, conc_ppm, corr_conc_ppm) %>%
  as.data.frame()

# Read out data
write.csv(all.data,
          "dataEdited/waterChemistry/ICP/metal_data.csv",
          row.names = FALSE)
