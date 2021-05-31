
# file.name.input <- "dataRaw/incubations/HgT/2020_incubations/I020321 BENDOTA WATER.xlsx"
# output.file.name <- "dataEdited/incubations/HgT/2020_incubations/20210203.csv"
# all.metadata.2020 <- read.csv("metadata/processedMetadata/incubation_Hg_metadata_2020.csv",
#                               stringsAsFactors = FALSE)
#


#### Clean HgT data ####
clean.HgT.data.from <- function(file.name.input,
                                output.file.name,
                                metadata.df = all.metadata.2020,
                                date.of.analysis,
                                barcodes.to.remove = NULL) {

  file.name <- file.name.input
  file.data <- read_xlsx(file.name,
                         sheet = "Summary",
                         skip = 10)

  # Select needed columns and rename them.
  file.data.clean <- file.data %>%
    select(1, 3, 4, 5, 6, 10) %>%
    as.data.frame()
  colnames(file.data.clean) <- c("sampleType",
                                 "bottleID",
                                 "amb_HgT_ng.L",
                                 "excess_HgT_198_ng.L",
                                 "excess_HgT_204_ng.L",
                                 "excess_HgT_DDL")
  file.data.clean <- file.data.clean %>%
    filter(sampleType == "SAM") %>%
    select(-sampleType) %>%
    mutate(dateAnalyzed_HgT = date.of.analysis)

  # Remove unwanted samples by barcode
  if (!is.null(barcodes.to.remove)) {

    file.data.clean <- file.data.clean %>%
      filter(!(bottleID %in% barcodes.to.remove))

  }

  # Remove entries lacking results
  file.data.clean <- file.data.clean %>%
    filter(!is.na(excess_HgT_198_ng.L),
           !is.na(excess_HgT_204_ng.L))

  # Round off numbers
  file.data.clean[, c("amb_HgT_ng.L", "excess_HgT_198_ng.L", "excess_HgT_204_ng.L", "excess_HgT_DDL")] <- file.data.clean[, c("amb_HgT_ng.L", "excess_HgT_198_ng.L", "excess_HgT_204_ng.L", "excess_HgT_DDL")] %>%
    mutate_all(as.numeric) %>%
    round(digits = 3)

  file.data.clean <- file.data.clean %>%
    mutate(above_DDL_HgT_198 = (excess_HgT_198_ng.L > excess_HgT_DDL)) %>%
    mutate(above_DDL_HgT_204 = (excess_HgT_204_ng.L > excess_HgT_DDL))


  # Write out data
  write.csv(file.data.clean,
            output.file.name,
            row.names = FALSE,
            quote = FALSE)
}







#### Clean MeHg data ####
clean.MeHg.data.from <- function(file.name.input,
                                 output.file.name,
                                 metadata.df = all.metadata.2020,
                                 date.of.analysis,
                                 barcodes.to.remove = NULL) {

  file.name <- file.name.input
  file.data <- read_xlsx(file.name,
                         sheet = "Summary",
                         skip = 6)

  # Select needed columns and rename them.
  file.data.clean <- file.data %>%
    select(1, 3, 4, 5, 9) %>%
    as.data.frame()
  colnames(file.data.clean) <- c("bottleID",
                                 "amb_MeHg_ng.L",
                                 "excess_MeHg_198_ng.L",
                                 "excess_MeHg_204_ng.L",
                                 "excess_MeHg_DDL")
  file.data.clean <- file.data.clean %>%
    filter(bottleID %in% metadata.df$bottleID)

  # Remove unwanted samples by barcode
  if (!is.null(barcodes.to.remove)) {

    file.data.clean <- file.data.clean %>%
      filter(!(bottleID %in% barcodes.to.remove))

  }

  # Remove spikes
  file.data.clean <- file.data.clean %>%
    filter(!duplicated(bottleID))

  # Remove entries lacking results
  file.data.clean <- file.data.clean %>%
    filter(!is.na(excess_MeHg_198_ng.L),
           !is.na(excess_MeHg_204_ng.L))

  # Round off numbers
  file.data.clean[, -1] <- file.data.clean[, -1] %>%
    mutate_all(as.numeric) %>%
    round(digits = 3)

  file.data.clean <- file.data.clean %>%
    mutate(above_DDL_MeHg_198 = (excess_MeHg_198_ng.L > excess_MeHg_DDL),
           above_DDL_MeHg_204 = (excess_MeHg_204_ng.L > excess_MeHg_DDL)) %>%
    mutate(dateAnalyzed_MeHg = date.of.analysis)


  # Write out data
  write.csv(file.data.clean,
            output.file.name,
            row.names = FALSE,
            quote = FALSE)
}










#### Clean sulfate data ####
clean.sulfate.data.from <- function(data.file.name,
                                    processing.metadata.df = processing.metadata,
                                    samples.to.remove = NULL,
                                    output.file.name) {
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
    mutate(rawConcentration = as.numeric(rawConcentration))

  # Read in analytical prep data
  prep.data <- read_xlsx(data.file.name,
                         sheet = "sample_prep") %>%
    rename(sulfurID = injectionID)

  # Combine data and adjust concentration for dilution
  sample.data <- inner_join(processing.metadata,
                            sulfate.data %>% select(sulfurID, rawConcentration)) %>%
    left_join(prep.data) %>%
    mutate(concentration = rawConcentration * dilution) %>%
    mutate(preservativeDilutionFactor = ((mass - tare) / (mass - tare - preservativeVol))) %>%
    mutate(concentration = concentration * preservativeDilutionFactor) %>%
    select(sulfurID, concentration) 
  if (!is.null(samples.to.remove)) {
    sample.data <- sample.data %>%
      filter(!(sulfurID %in% samples.to.remove))
    
  }
  


  # Write out data
  write.csv(sample.data,
            output.file.name,
            row.names = FALSE,
            quote = FALSE)

}



#### Clean exo data ####
clean.exo.data <- function(file.name = file.name.of.interest,
                           output.name = output.name.of.interest,
                           pH.sensor.present = TRUE,
                           depth.adjustment = 0) {
  # Read in data
  exo.data <- read_xlsx(paste("dataRaw/exo/",
                              file.name,
                              sep = ""),
                        skip = 24) %>%
    as.data.frame()

  # Clean up datafile
  # Read in custom headers that I made (found in sonde_profile_headers.csv)
  sonde.names <- c("date", "time", "time_fraction", "site_name", "fault_code", "Battery_V",
                   "Cable_Pwr_V", "Temp_C", "Cond_µS.cm", "SpCond_µS.cm", "Sal_psu",
                   "nLF_Cond_µS.cm", "TDS_mg.L", "ODO_sat", "ODO_conc", "fDOM_RFU", "fDOM_QSU",
                   "pH", "pH_mV", "ORP_mV", "ORP_raw_mV", "Turbidity_FNU", "TSS_mg.L",
                   "Chlorophyll_RFU", "Chlorophyll_µg.L", "BGA-PC_RFU", "BGA-PC_µg.L",
                   "Press_psi_a", "depth")
  if (pH.sensor.present == FALSE) {
    sonde.names <- c("date", "time", "time_fraction", "site_name", "fault_code", "Battery_V",
                     "Cable_Pwr_V", "Temp_C", "Cond_µS.cm", "SpCond_µS.cm", "Sal_psu",
                     "nLF_Cond_µS.cm", "TDS_mg.L", "ODO_sat", "ODO_conc", "fDOM_RFU", "fDOM_QSU",
                     "Turbidity_FNU", "TSS_mg.L", "Chlorophyll_RFU", "Chlorophyll_µg.L",
                     "BGA-PC_RFU", "BGA-PC_µg.L", "Press_psi_a", "depth")
  }
  # Add names, clean the data
  colnames(exo.data) <- sonde.names

  # Factor in depth adjustment
  exo.data$depth <- exo.data$depth - depth.adjustment

  # Start at 1 meters
  exo.data <- exo.data[which(exo.data$depth >= 1), ]
  # End at maximum depth
  exo.data <- exo.data[1:which.max(exo.data$depth), ]
  # Order by depth and time
  exo.data <- exo.data[order(exo.data$depth, exo.data$time), ]
  # Round depths to nearest 0.1m and average these measurements
  exo.data$depth <- round(exo.data$depth,
                          digits = 1)

  if (pH.sensor.present == TRUE) {

    exo.data <- exo.data %>%
      select(Temp_C, Cond_µS.cm, SpCond_µS.cm, Sal_psu,
             ODO_sat, ODO_conc, fDOM_RFU, pH, ORP_mV,
             Turbidity_FNU, Chlorophyll_RFU, depth)
  } else {
    exo.data <- exo.data %>%
      select(Temp_C, Cond_µS.cm, SpCond_µS.cm, Sal_psu,
             ODO_sat, ODO_conc, fDOM_RFU, Turbidity_FNU,
             Chlorophyll_RFU, depth)
  }

  exo.data <- aggregate(exo.data,
                        by = list(exo.data$depth),
                        FUN = mean) %>%
    select(-Group.1)

  # Write out csv of data
  write.csv(exo.data,
            file = paste("dataEdited/exo/",
                         output.name,
                         sep = ""),
            row.names = FALSE)

}
