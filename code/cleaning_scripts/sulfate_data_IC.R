# code/waterChemistry/clean_sulfate_data.R

# This script will clean up the sulfate data we
# generated in WSEL on the IC and combine it with
# the metadata

#### Always start with a clean slate ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP/")
library(lubridate)
library(readxl)
library(tidyverse)


#### Prepare metadata ####
S.metadata <- read_xlsx("metadata/raw_metadata/chem_S.xlsx") %>%
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



#### Define sulfate data cleaning function ####
data.file.name = "dataRaw/waterChemistry/sulfate/20211203_BENDOTA.xlsx"
processing.metadata.df = processing.metadata
output.file.name = "dataEdited/waterChemistry/sulfate/20210316_BENDOTA.csv"
samples.to.remove = NULL

clean.sulfate.data.from <- function(data.file.name,
                                    processing.metadata.df = processing.metadata,
                                    samples.to.remove = NULL,
                                    output.file.name,
                                    final.say = "NO") {
  # Read in headers for Excel dataframe
  sulfate.headers <- read_xlsx(data.file.name,
                               sheet = "headers",
                               n_max = 1,
                               col_names = FALSE) %>%
    unlist(use.names = FALSE)
  
  # Read in sulfate data
  sulfate.data <- read_xlsx(data.file.name,
                            sheet = "raw_results",
                            skip = 4,
                            col_names = sulfate.headers) %>%
    filter(injectionID != "NA",
           concentration != "n.a.") %>%
    rename(sulfurID = injectionID,
           rawConcentration = concentration) %>%
    mutate(rawConcentration = as.numeric(rawConcentration)) %>%
    select(sulfurID, type, rawConcentration)
  
  # Check check standards
  check.standards <-sulfate.data %>%
    filter(grepl("ppm", sulfurID) &
             type == "Unknown") %>%
    mutate(expected_concentration = sulfurID %>%
             strsplit("ppm") %>% sapply("[", 1) %>%
             as.numeric(),
           recovery = rawConcentration / expected_concentration *100)
  print(check.standards)
  
  # Check blanks
  check.blanks <-sulfate.data %>%
    filter(grepl("BLANK", sulfurID))
  print(check.blanks)
  check.blanks$rawConcentration*6
  
  # Read in analytical prep data and adjust for dilution
  prep.data <- read_xlsx(data.file.name,
                         sheet = "sample_prep") %>%
    rename(sulfurID = injectionID) %>%
    right_join(sulfate.data) %>%
    mutate(concentration = rawConcentration * dilution) %>%
    select(-rawConcentration)
  # Check duplicates
  dup.ids <- prep.data %>%
    filter(grepl("_dup", sulfurID)) %>%
    mutate(sulfurID = sulfurID %>%
             strsplit("_dup") %>% sapply("[", 1)) %>%
    select(sulfurID) %>%
    unlist(use.names = FALSE)
  dup.data <- prep.data %>%
    mutate(dup = grepl("_dup", sulfurID)) %>%
    mutate(sulfurID = sulfurID %>%
             strsplit("_dup") %>% sapply("[", 1)) %>%
    filter(sulfurID %in% dup.ids) %>%
    spread(key = dup,
           value = concentration) %>%
    mutate(percent.difference = (`FALSE` - `TRUE`)/`FALSE` * 100)
  print(dup.data)
  
  if (final.say == "NO") {
    return("Shut it down")
  } else {
    # Finalize data 
    prep.data.final <- prep.data %>%
      mutate(sulfurID = sulfurID %>%
               strsplit("_dup") %>% sapply("[", 1)) %>%
      group_by(sulfurID) %>%
      summarise(concentration = mean(concentration)) %>%
      inner_join(processing.metadata) %>%
      mutate(preservativeDilutionFactor = ((mass - tare) / (mass - tare - preservativeVol))) %>%
      mutate(concentration = concentration * preservativeDilutionFactor) %>%
      select(sulfurID, concentration) 
    
    if (!is.null(samples.to.remove)) {
      prep.data.final <- prep.data.final %>%
        filter(!(sulfurID %in% samples.to.remove))
      }
    # Write out data
    write.csv(prep.data.final,
              output.file.name,
              row.names = FALSE,
              quote = FALSE)
    
  }
}



#### Run cleaning function ####
# clean.sulfate.data.from(data.file.name = "dataRaw/waterChemistry/sulfate/20210316_BENDOTA.xlsx",
#                         processing.metadata.df = processing.metadata,
#                         output.file.name = "dataEdited/waterChemistry/sulfate/20210316_BENDOTA.csv",
#                         final.say = "YES")
# clean.sulfate.data.from(data.file.name = "dataRaw/waterChemistry/sulfate/20210423_BENDOTA.xlsx",
#                         processing.metadata.df = processing.metadata,
#                         samples.to.remove = c("BLI20_TS_087", "BLI20_TS_088", "BLI20_TS_089", "BLI20_TS_090",
#                                               "BLI20_TS_091", "BLI20_TS_092", "BLI20_TS_092_duplicate",
#                                               "BLI20_TS_093", "BLI20_TS_094"),
#                         output.file.name = "dataEdited/waterChemistry/sulfate/20210423_BENDOTA.csv")
# clean.sulfate.data.from(data.file.name = "dataRaw/waterChemistry/sulfate/20211203_BENDOTA.xlsx",
#                         processing.metadata.df = processing.metadata,
#                         output.file.name = "dataEdited/waterChemistry/sulfate/20211203_BENDOTA.csv",
#                         samples.to.remove = c(#"BLI21_TS_015", "BLI21_TS_014",
#                                               "BLI21_TS_015_dup"),
#                         final.say = "YES")


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
  left_join(WC.metadata) %>%
  arrange(sulfurID) %>%
  mutate(depthOriginal = depth) %>%
  mutate(corewater = grepl(pattern = "-",
                           x = depthOriginal)) %>%
  mutate(sulfate_uM = round(concentration / 96.07 * 1000,
                            1)) %>%
  filter(!is.na(sampleID))
WC.results[WC.results$corewater, ] <- WC.results[WC.results$corewater, ] %>%
  mutate(depth = paste("-",
                       strsplit(depth, "-") %>% sapply("[", 2),
                       sep = "") %>%
           as.numeric() / 100) %>%
  mutate(depth = depth + (waterDepth * corewater))


#### Add replicate numbers ####
WC.results <- WC.results %>%
  arrange(sulfurID) %>%
  group_by(startDate, depth, corewater) %>%
  mutate(replicate = row_number()) %>%
  ungroup()


#### Save out data ####
write.csv(WC.results,
          "dataEdited/waterChemistry/sulfate/WC_data.csv",
          row.names = FALSE,
          quote = FALSE)


#### WC samples to analyze yet ####

WC.samples.to.analyze <- WC.metadata %>%
  left_join(S.results) %>%
  filter(is.na(concentration)) %>%
  filter(year(startDate) == 2021)
