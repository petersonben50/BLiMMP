#### code/cleaning_scripts/clean_leucine_uptake_2020.R ####
# Benjamin D. Peterson

#### Clean up ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP")
library(readxl)
library(tidyverse)


#### Read in metadata ####

# Water column data
leucine.metadata <- read_xlsx("metadata/LEU.xlsx") %>%
  select(-notes)
sampleID.metadata <- read_xlsx("metadata/2_sample_IDs.xlsx") %>%
  select(-notes)
tripID.metadata <- read_xlsx("metadata/1_trip_IDs.xlsx") %>%
  select(-notes)

# Clean up water column data
metadata.WC <- leucine.metadata %>%
  left_join(sampleID.metadata) %>%
  left_join(tripID.metadata) %>%
  select(uptakeID, sampleID, depth, startDate)
write.csv(metadata.WC,
          "metadata/clean_LEU_metadata.csv",
          row.names = FALSE)



#### Function to calculate C analyzed ####
# testdatafile = "dataRaw/leucineUptake/20201011_leucine.xlsx"

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
                                   metadata.of.interest = metadata.WC,
                                   output.file) {
  
  
  
  # Set some standards
  
  # We set the time that each of these incubated for in hr
  incubation_hrs <- as.numeric(incubation_minutes) / 60
  # We also need the total volume of the incubation in L.
  volume_L <- as.numeric(volume_mL) / 1000
  # Set molar mass leucine
  molar.mass.leucine <- 131.2
  
  
  # Calculate blank values
  blank.values <- read_xlsx(datafile,
                            sheet = "leucine_standards_counts") %>%
    filter(grepl("neg", leucineID)) %>%
    select(leucineID, CPMA, CPMB) %>%
    group_by(leucineID) %>%
    summarise(CPMA_blank = mean(CPMA),
              CPMB_blank = mean(CPMB)) %>%
    select(CPMA_blank, CPMB_blank) %>%
    unlist()
  
  # Read in spiking info
  spike.info <- read_xlsx(datafile,
                          sheet = "leucine_solutions") %>%
    select(leucineID, final_leu_stock_uM)
  
  # Read in counts data for standards
  efficiency.data <- read_xlsx(datafile,
                               sheet = "leucine_standards_counts") %>%
    filter(grepl("spike", leucineID)) %>%
    mutate(CPMA = CPMA - blank.values["CPMA_blank"],
           CPMB = CPMB - blank.values["CPMB_blank"])
  

  # We'll divide the number of counts that we got to calculate 
  # the overall efficiency of the counts, which in this case 
  # means the counts that we get per leucine in solution.
  # This counting efficienty (CF) will provide a direct link
  # between the counts we got and the amount of leucine taken up.
  # This CF value will be in units of counts per pmol of leucine.
  
  CF.calcs <- full_join(spike.info,
            efficiency.data) %>%
    mutate(CF_A = CPMA / (final_leu_stock_uM * volume_ul)) %>%
    mutate(CF_B = CPMB / (final_leu_stock_uM * volume_ul)) %>%
    group_by(leucineID) %>%
    summarise(CF_A = median(CF_A),
              CF_B = median(CF_B))
  

  
  
  # Then we'll read in the raw data and join it with the CF
  # calculations
  leucine.data <- read_xlsx(datafile,
                                 sheet = "results") %>%
    mutate(CPMA = CPMA - blank.values["CPMA_blank"],
           CPMB = CPMB - blank.values["CPMB_blank"]) %>%
    full_join(CF.calcs)
  
 
  # Now for the calculations.
  # For both A and B counting channels, we'll first take the counts and subtract out the blank values.
  # Then we'll divide the counts by the efficiency data, to calculate the pmol of leucine that have been
  # taken up. Then we'll divide by the incubation volume to get pmol/L. And then we'll divide it by
  # the number of hours of incubation, to normalize production to 1 hour.
  # This gives us the pM of leucine incorporated per hour.
  
  C.incorporation <- leucine.data %>%
    mutate(Leu_uptake_pM_A = CPMA / (CF_A * incubation_hrs * volume_L)) %>%
    mutate(Leu_uptake_pM_B = CPMB / (CF_B * incubation_hrs * volume_L)) %>%
    #mutate(ratio_A_B = Leu_uptake_pM_A / Leu_uptake_pM_B) %>%
    
    # To convert this to overall pM of protein produced per house, we divide the pM of leucine by the
    # molar ratio of leucine in protein, which is approximately 7.3% and multiply by the molar mass of 
    # leucine and the leucine isotope dilution, which here we are assuming is 2.
    # I think this is a BS way of doing it, since the 7.3% value is a molar ratio, not a mass ratio.
    # But, this is how people have been calculating it, so I'm gonna stick with it. At some
    # point I should use Simon and Azam's tables to calculate the mass ratio.
    mutate(µgBPP_per_L_hr = round(((Leu_uptake_pM_A / Leu_per_prot_molar_ratio) * molar.mass.leucine * leucine_ID) / 1000000, 3)) %>%
    
    # The last conversion to make is from protein biomass to overall carbon biomass.
    mutate(µgBCP_per_L_hr = µgBPP_per_L_hr * cell_C_to_protein) %>%
    select(c(uptakeID, Leu_uptake_pM_A, µgBPP_per_L_hr, µgBCP_per_L_hr))

  # Combine data with metadata
  C.incorporation <- C.incorporation %>%
    left_join(metadata.of.interest)
  
  # Read in amendment info
  amendment.info <- read_xlsx(datafile,
                              sheet = "spiking_sheet") %>%
    select(uptakeID, killed, amendmentID)
  
  if (multiple.amendments == FALSE) {
    # Add amendment info
    amendment.info = amendment.info %>%
      select(-amendmentID)
  }
  
  C.incorporation <- C.incorporation %>%
    left_join(amendment.info)

  # Write out csv file
  write.csv(C.incorporation,
            output.file,
            row.names = FALSE)  

}




#### Convert leucine uptake to C growth ####
convert_leuUptake_to_C("dataRaw/leucineUptake/20200917_leucine.xlsx",
                       multiple.amendments = TRUE,
                       output.file = "dataEdited/leucineUptake/20200917_leucine.csv")

convert_leuUptake_to_C("dataRaw/leucineUptake/20201011_leucine.xlsx",
                       output.file = "dataEdited/leucineUptake/20201011_leucine.csv")
