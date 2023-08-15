#### code/cleaning_scripts/clean_leucine_uptake_2020.R ####
# Benjamin D. Peterson

#### Clean up ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP")
library(readxl)
library(tidyverse)


#### Set conditions for calculations ####
incubation_minutes <- 60
volume_mL <- 1
molar.mass.leucine <- 131.2
cell_C_to_protein <- 0.86
# I'm assuming that the ID of leucine is 2, as described by Simon and Azam, 1989.
leucine_ID <- 2
# Set the molar ratio of leucine in protein as 7.3%, according to Exercise 19 in
# Limnological Analyses by Wetzel and Likens. This ratio was determined in Simon
# and Azam, 1989
Leu_per_prot_molar_ratio <- 0.073



#### Read in metadata ####
leu_metadata <- read.csv("metadata/processedMetadata/LEU_metadata.csv")



#### Define functions for each year to calculate C incorporated ####
convert_leuUptake_to_C_2020 <- function(datafile = testdatafile,
                                        multiple.amendments = FALSE,
                                        metadata.of.interest = leu_metadata,
                                        output.file) {

  # We set the time that each of these incubated for in hr
  incubation_hrs <- as.numeric(incubation_minutes) / 60
  # We also need the total volume of the incubation in L.
  volume_L <- as.numeric(volume_mL) / 1000  
  
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
    mutate(Leu_uptake_pM_per_hour = CPMA / (CF_A * incubation_hrs * volume_L)) %>%
    mutate(Leu_uptake_pM_B = CPMB / (CF_B * incubation_hrs * volume_L)) %>%
    #mutate(ratio_A_B = Leu_uptake_pM_A / Leu_uptake_pM_B) %>%
    
    # To convert this to overall pM of protein produced per house, we divide the pM of leucine by the
    # molar ratio of leucine in protein, which is approximately 7.3% and multiply by the molar mass of 
    # leucine and the leucine isotope dilution, which here we are assuming is 2.
    # I think this is a BS way of doing it, since the 7.3% value is a molar ratio, not a mass ratio.
    # But, this is how people have been calculating it, so I'm gonna stick with it. At some
    # point I should use Simon and Azam's tables to calculate the mass ratio.
    mutate(µgBPP_per_L_hr = ((Leu_uptake_pM_per_hour / Leu_per_prot_molar_ratio) * molar.mass.leucine * leucine_ID) / 1000000) %>%
    
    # The last conversion to make is from protein biomass to overall carbon biomass.
    mutate(µgBCP_per_L_hr = µgBPP_per_L_hr * cell_C_to_protein) %>%
    select(c(uptakeID, Leu_uptake_pM_per_hour, µgBPP_per_L_hr, µgBCP_per_L_hr))
  # 
  # # Combine data with metadata
  # C.incorporation <- C.incorporation %>%
  #   left_join(metadata.of.interest)
  
  # Read in amendment info
  amendment.info <- read_xlsx(datafile,
                              sheet = "spiking_sheet") %>%
    select(uptakeID, killed, amendmentID)
  
  if (multiple.amendments == FALSE) {
    # Add amendment info
    amendment.info = amendment.info %>%
      select(-amendmentID)
    C.incorporation <- C.incorporation %>%
      mutate(amendmentID = NA)
  }
  
  C.incorporation <- C.incorporation %>%
    left_join(amendment.info) %>%
    select(c(uptakeID, Leu_uptake_pM_per_hour, µgBPP_per_L_hr, µgBCP_per_L_hr, killed, amendmentID))
  
  # Write out csv file
  write.csv(C.incorporation,
            output.file,
            row.names = FALSE)  
  
}

convert_leuUptake_to_C_2021 <- function(datafile = testdatafile,
                                        multiple.amendments = FALSE,
                                        metadata.of.interest = leu.metadata,
                                        output.file) {
  
  # We set the time that each of these incubated for in hr
  incubation_hrs <- as.numeric(incubation_minutes) / 60
  # We also need the total volume of the incubation in L.
  volume_L <- as.numeric(volume_mL) / 1000
  # Set molar mass leucine  
 
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
    mutate(Leu_uptake_pM_per_hour = leu_pmol / (incubation_hrs * volume_L)) %>%
    # To convert this to overall pM of protein produced per house, we divide the pM of leucine by the
    # molar ratio of leucine in protein, which is approximately 7.3% and multiply by the molar mass of 
    # leucine and the leucine isotope dilution, which here we are assuming is 2.
    # I think this is a BS way of doing it, since the 7.3% value is a molar ratio, not a mass ratio.
    # But, this is how people have been calculating it, so I'm gonna stick with it. At some
    # point I should use Simon and Azam's tables to calculate the mass ratio.
    mutate(µgBPP_per_L_hr = ((Leu_uptake_pM_per_hour / Leu_per_prot_molar_ratio) * molar.mass.leucine * leucine_ID) / 1000000) %>%
    
    # The last conversion to make is from protein biomass to overall carbon biomass.
    mutate(µgBCP_per_L_hr = µgBPP_per_L_hr * cell_C_to_protein) %>%
    select(c(uptakeID, Leu_uptake_pM_per_hour, µgBPP_per_L_hr, µgBCP_per_L_hr)) %>%
    mutate(killed = NA,
           amendmentID = NA)


  # Write out csv file
  write.csv(leucine.data,
            output.file,
            row.names = FALSE)  

}



#### Convert leucine uptake to C growth ####
convert_leuUptake_to_C_2020("dataRaw/leucineUptake/20200917_leucine.xlsx",
                            multiple.amendments = TRUE,
                            output.file = "dataEdited/leucineUptake/20200917_leucine.csv")

convert_leuUptake_to_C_2020("dataRaw/leucineUptake/20201011_leucine.xlsx",
                            output.file = "dataEdited/leucineUptake/20201011_leucine.csv")
convert_leuUptake_to_C_2021("dataRaw/leucineUptake/20210916_leucine.xlsx",
                            multiple.amendments = TRUE,
                            output.file = "dataEdited/leucineUptake/20210916_leucine.csv")
convert_leuUptake_to_C_2021("dataRaw/leucineUptake/20210917_leucine.xlsx",
                            multiple.amendments = TRUE,
                            output.file = "dataEdited/leucineUptake/20210917_leucine.csv")
convert_leuUptake_to_C_2021("dataRaw/leucineUptake/20211018_leucine.xlsx",
                            multiple.amendments = TRUE,
                            output.file = "dataEdited/leucineUptake/20211018_leucine.csv")



#### Combine all leucine uptake data ####
rm(list = ls())
# Read in metadata
leu_metadata <- read.csv("metadata/processedMetadata/LEU_metadata.csv")
# Read in data
list_of_data_files <- list.files("dataEdited/leucineUptake/",
                                 pattern = "_leucine.csv",
                                 full.names = TRUE)
leucine_uptake_data <- do.call(rbind,
                               lapply(list_of_data_files,
                                      function(csv_datafile) {
                                        read.csv(csv_datafile)
                                        })
                               )
rm(list_of_data_files)

# Merge data and clean it up
leucine_data_metadata <- inner_join(leu_metadata,
                                    leucine_uptake_data)

treatment_vector <- rep(NA, dim(leucine_data_metadata)[1])
treatment_vector[which(leucine_data_metadata$treatment == "molybdate")] <- "molybdate"
treatment_vector[which(leucine_data_metadata$treatment == "control")] <- "ambient"
treatment_vector[which(leucine_data_metadata$treatment == "filtered")] <- "control"
treatment_vector[which(leucine_data_metadata$killed == "no")] <- "ambient"
treatment_vector[which(leucine_data_metadata$killed == "yes")] <- "control"
leucine_data_metadata <- leucine_data_metadata %>%
  mutate(treatment = treatment_vector) %>%
  select(uptakeID, sampleID, startDate, depth, treatment, protocol, timePoint,
         Leu_uptake_pM_per_hour, µgBPP_per_L_hr, µgBCP_per_L_hr) %>%
  mutate(Leu_uptake_pM_per_hour = round(Leu_uptake_pM_per_hour, 0),
         µgBPP_per_L_hr = round(µgBPP_per_L_hr, 2),
         µgBCP_per_L_hr = round(µgBCP_per_L_hr, 2)) %>%
  filter(treatment != "control")
