#### code/cleaning_scripts/clean_leucine_uptake_2020.R ####
# Benjamin D. Peterson

#### Clean up ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP")
library(readxl)
library(tidyverse)


#### Read in metadata ####

leu.metadata <- read.csv("metadata/processedMetadata/LEU_2021.csv")



#### Function to calculate C analyzed ####
# testdatafile = "dataRaw/leucineUptake/20210916_leucine.xlsx"

# Default conditions for calculations. You can modify these with the inputs.

# I'm assuming that the ID of leucine is 2, as described by Simon and Azam, 1989.

# Set the molar ratio of leucine in protein as 7.3%, according to Exercise 19 in
# Limnological Analyses by Wetzel and Likens. This ratio was determined in Simon
# and Azam, 1989

convert_leuUptake_to_C <- function(datafile = testdatafile,
                                   multiple.amendments = FALSE,
                                   incubation_minutes = 60,
                                   volume_mL = 1,
                                   leucine_ID = 2,
                                   Leu_per_prot_molar_ratio = 0.073,
                                   cell_C_to_protein = 0.86,
                                   metadata.of.interest = leu.metadata,
                                   output.file) {
  
  
  
  # Set some standards
  
  # We set the time that each of these incubated for in hr
  incubation_hrs <- as.numeric(incubation_minutes) / 60
  # We also need the total volume of the incubation in L.
  volume_L <- as.numeric(volume_mL) / 1000
  # Set molar mass leucine
  molar.mass.leucine <- 131.2
  
  
  # # Calculate blank values
  # blank.values <- read_xlsx(datafile,
  #                           sheet = "leucine_standards_counts") %>%
  #   filter(grepl("neg", leucineID)) %>%
  #   select(leucineID, CPMA, CPMB) %>%
  #   group_by(leucineID) %>%
  #   summarise(CPMA_blank = mean(CPMA),
  #             CPMB_blank = mean(CPMB)) %>%
  #   select(CPMA_blank, CPMB_blank) %>%
  #   unlist()
  
  # Read in standards data
  spike.info <- read_xlsx(datafile,
                          sheet = "leucine_standards_counts") %>%
    select(`Leucine mass (pmoles)`, CPMA) %>%
    rename(leuMass = `Leucine mass (pmoles)`)

  # Normally the amount of leucine taken up is calculated by
  # first getting a measure of the counting efficiency.
  # However, it is just as easy for me to make a standard
  # curve of the concentration of the leucine against the
  # number of counts. So, we'll do that here.
  std.curve.regression <- lm(spike.info$CPMA ~ spike.info$leuMass)
  
  
  
  # Then we'll read in the raw data and calculate the pmol
  # of leucine
  leucine.data <- read_xlsx(datafile,
                                 sheet = "results") %>%
    mutate(leu_pmol = (CPMA - std.curve.regression$coefficients[1]) / std.curve.regression$coefficients[2]) %>%
    # Then we'll divide by the incubation volume and incubation hours to get the right value
    mutate(leu_pM_per_hour = leu_pmol / (incubation_hrs * volume_L)) %>%
    # To convert this to overall pM of protein produced per house, we divide the pM of leucine by the
    # molar ratio of leucine in protein, which is approximately 7.3% and multiply by the molar mass of 
    # leucine and the leucine isotope dilution, which here we are assuming is 2.
    # I think this is a BS way of doing it, since the 7.3% value is a molar ratio, not a mass ratio.
    # But, this is how people have been calculating it, so I'm gonna stick with it. At some
    # point I should use Simon and Azam's tables to calculate the mass ratio.
    mutate(µgBPP_per_L_hr = round(((leu_pM_per_hour / Leu_per_prot_molar_ratio) * molar.mass.leucine * leucine_ID) / 1000000, 3)) %>%
    
    # The last conversion to make is from protein biomass to overall carbon biomass.
    mutate(µgBCP_per_L_hr = µgBPP_per_L_hr * cell_C_to_protein) %>%
    select(c(uptakeID, leu_pM_per_hour, µgBPP_per_L_hr, µgBCP_per_L_hr))

  # Combine data with metadata
  C.incorporation <- leucine.data %>%
    left_join(metadata.of.interest)

  # Write out csv file
  write.csv(C.incorporation,
            output.file,
            row.names = FALSE)  

}




#### Convert leucine uptake to C growth ####
convert_leuUptake_to_C("dataRaw/leucineUptake/20210916_leucine.xlsx",
                       multiple.amendments = TRUE,
                       output.file = "dataEdited/leucineUptake/20210916_leucine.csv")
convert_leuUptake_to_C("dataRaw/leucineUptake/20210917_leucine.xlsx",
                       multiple.amendments = TRUE,
                       output.file = "dataEdited/leucineUptake/20210917_leucine.csv")
