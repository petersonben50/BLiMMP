#### code/cleaning_scripts/clean_sulfide_data.R ####
# Benjamin D. Peterson

#### Get ready, and fly! ####
rm(list = ls())
setwd("~/Box/BLiMMP/")
library(dplyr)
library(lubridate)
library(tidyr)
library(readxl)

#### Prepare metadata ####
S.metadata <- read_xlsx("metadata/chem_S.xlsx")
sampleID <- read_xlsx("metadata/2_sample_IDs.xlsx")
incubationID <- read_xlsx("metadata/4_MA_ID.xlsx")
tripID <- read_xlsx("metadata/1_trip_IDs.xlsx")

#### Save out water column data ####

WC.metadata <- S.metadata %>%
  filter(!(sampleID == "NA")) %>%
  left_join(sampleID) %>%
  left_join(tripID) %>%
  select(sulfurID, sampleID, depth, startDate)
write.csv(WC.metadata,
          file = "metadata/processedMetadata/sulfide_WC.csv",
          row.names = FALSE,
          quote = FALSE)

#### Save out incubation data ####

MA.metadata <- S.metadata %>%
  filter(!(incubationID == "NA")) %>%
  select(sulfurID, incubationID) %>%
  left_join(incubationID) %>%
  left_join(sampleID) %>%
  left_join(tripID) %>%
  select(sulfurID, incubationID, depth, dateCollected, filtered, amendment)
write.csv(MA.metadata,
          file = "metadata/processedMetadata/sulfide_MA.csv",
          row.names = FALSE,
          quote = FALSE)

#### Generate needed processing metadata ####

processing.metadata <- S.metadata %>%
  select(sulfurID, mass, tare, preservativeVol) %>%
  mutate_at(.vars = c("mass", "tare", "preservativeVol"),
            as.numeric)
  
  
# Combine incubation data

rm(S.metadata,
   sampleID,
   incubationID,
   tripID)


#### Define output locations ####
report.location <- "dataEdited/waterChemistry/sulfide/reports/"
bad.data.location <- "dataEdited/waterChemistry/sulfide/badData/"
good.data.location <- "dataEdited/waterChemistry/sulfide/dataForReview"

#### Data processing function ####
data_file <- "dataRaw/waterChemistry/sulfide/sulfide_20191023.xlsx"

# Define the function
data_processing_function <- function(data_file,
                                     override = "off") {

  keeper.status <- "pass"
  
  #### Read in needed data ####
  
  # Read in raw data
  data.spreadsheet.raw <- read_xlsx(data_file,
                                sheet = "raw_results") %>%
    as.data.frame()
  
  # Read in run info
  run.info <- read_xlsx(data_file,
                        sheet = "runInfo",
                        col_names = c("variable", "value")) %>%
    as.data.frame()
  run.info.vector <- run.info$value
  names(run.info.vector) <- run.info$variable
  rm(run.info)
  
  # Read in info on standards
  standard.info <- read_xlsx(data_file,
                             sheet = paste(run.info.vector["standard_curve"],
                                           "_curve",
                                           sep = "")) %>%
    select(sulfurID,
           sulfide_uM)
  
  # Pull out date of analysis
  date.of.analysis <- data_file %>%
    strsplit("sulfide_") %>%
    sapply("[", 2) %>%
    strsplit(".csv") %>%
    sapply("[", 1) %>%
    ymd()
                              

  #### Prep blank data ####
  blank.data <- data.spreadsheet.raw %>%
    filter(analysisType == "BLK") %>%
    select(Absorbance_667) %>%
    unlist(use.names = FALSE) 
  blank.value <- blank.data %>%
    mean() %>%
    round(3)
  blank.sd <- blank.data %>%
    sd() %>%
    round(3)
  
  
  #### Generate standard curve ####
  curve.data <-  data.spreadsheet.raw %>%
    filter(analysisType == "STD") %>%
    # filter(sulfurID != "H_6") %>%
    left_join(standard.info) %>%
    mutate(Ab667 = Absorbance_667 - blank.value) %>%
    select(sulfide_uM,
           Ab667,
           sulfurID)
  std.curve.regression <- lm(curve.data$Ab667 ~ curve.data$sulfide_uM)
  
  #### QC: Check quality of curve ####
  summary(std.curve.regression)
  adj.r.squared <- round(summary(std.curve.regression)$adj.r.squared, 5)
  if (adj.r.squared < 0.995) {
    std.curve.qc <- "fail"
  } else {
    std.curve.qc <- "pass"
  }
  
  # Read out image of curve
  pdf.of.curve <- paste(report.location,
                        date.of.analysis,
                        "_curve.pdf",
                        sep = "")
  pdf(pdf.of.curve,
      height = 5,
      width = 5)
  par(mar = c(4.5, 4.5, 3, 1))
  plot(x = curve.data$sulfide_uM,
       y = curve.data$Ab667,
       ylab = "Absorbance (667nm)",
       xlab = "Sulfide (µM)",
       main = date.of.analysis)
  abline(std.curve.regression)
  dev.off()
  
  #### QC: Check blanks ####
  
  CCB.values <- data.spreadsheet.raw %>%
    filter(analysisType == "CCB")
  CCB.values
  data.spreadsheet.raw <- data.spreadsheet.raw %>%
    filter(analysisType != "CCB")
  
  
  #### Calculate concentrations using standard curve ####
  
  curve.coeff <- signif(std.curve.regression$coefficients[2], 4)
  inter.coeff <- signif(std.curve.regression$coefficients[1], 4)
  data.spreadsheet <- data.spreadsheet.raw %>%
    mutate(S_conc_uM = ((Absorbance_667 - blank.value) - inter.coeff) / curve.coeff) %>%
    filter(analysisType != "STD") %>%
    filter(analysisType != "BLK") 
  
  #### Calculate detection limit: NOT DONE FOR SULFIDE ####
  
  # LOD.fluorescence <- inter.coeff + 3*summary(std.curve.regression)$coefficients[2,2]
  # LOD.mass <- (LOD.fluorescence - std.curve.regression$coefficients[1]) / std.curve.regression$coefficients[2]
  
  
  #### QC: Ensure CCVs are within 10% ####
  # Known concentration is 6mg/L
  CCV.values <- data.spreadsheet %>%
    filter(analysisType == "CCV") %>%
    left_join(curve.data) %>%
    select(sulfurID,
           S_conc_uM,
           sulfide_uM) %>%
    mutate(deviation = (S_conc_uM - sulfide_uM)/S_conc_uM*100)
  data.spreadsheet <- data.spreadsheet %>%
    filter(analysisType != "CCV")
  if (any(CCV.values$deviation > 10)) {
    ccv.qc <- "fail"
  } else {
    ccv.qc <- "pass"
  }
  
  #### QC: Check matrix spike ####
  
  # Pull out needed data
  spiked.sample <- data.spreadsheet[which(data.spreadsheet$analysisType == "MS"), "sulfurID"]
  MS.samples <- data.spreadsheet %>%
    filter(sulfurID == spiked.sample)
  
  # Calculate contribution from sample
  contribution.from.sample.nmol <- MS.samples %>%
    filter(analysisType == "SAM") %>%
    select(S_conc_uM) %>%
    unlist(use.names = FALSE) * (as.numeric(run.info.vector["sample_added_ul"]) - as.numeric(run.info.vector["spike_volume_ul"])) / 1000
  
  # Calculate total mass in spiked solution
  spiked.solution.nmol <- MS.samples %>%
    filter(analysisType == "MS") %>%
    select(S_conc_uM)%>%
    unlist(use.names = FALSE)
  
  # Calculate added mass
  added.mass.nmol <- as.numeric(run.info.vector["spike_concentration_uM"]) * as.numeric(run.info.vector["spike_volume_ul"]) / 1000
  
  # Calculate percent recovery
  
  MS.recovery <- signif(((spiked.solution.nmol - contribution.from.sample.nmol) / added.mass.nmol) * 100, 3)
  MS.recovery
  
  # Remove MS sample 
  data.spreadsheet <- data.spreadsheet %>%
    filter(analysisType != "MS")
  
  # Check replicates are within 15%
  if (MS.recovery > 115 | MS.recovery < 85) {
    ms.qc <- "fail"
  } else {
    ms.qc <- "pass"
  }
  
  
  #### QC: Check on triplicates ####
  
  # Find triplicates
  triplicated.samples <- data.spreadsheet %>%
    filter(analysisType == "TRIP1") %>%
    select(sulfurID) %>%
    unlist(use.names = FALSE) %>%
    unique()
  
  # Check triplicates
  rep.data <- data.spreadsheet %>%
    filter(sulfurID %in% triplicated.samples) %>%
    group_by(sulfurID) %>%
    summarise(mean.rep = mean(S_conc_uM),
              rsd.rep = sd(S_conc_uM)/mean(S_conc_uM)*100)
  
  # Check replicates are within 15%
  if (any(rep.data$rsd.rep > 15)) {
    rep.qc <- "fail"
  } else {
    rep.qc <- "pass"
  }
  
  # Take mean of replicates
  data.spreadsheet <- data.spreadsheet %>%
    group_by(sulfurID) %>%
    summarise(S_conc_uM = mean(S_conc_uM))
  
  
  #### Overall QC validation ####
  all.qc <- c(std.curve.qc,
              ccv.qc,
              ms.qc,
              rep.qc)
  if (any(all.qc == "fail") && override == "off") {
    keeper.status = "fail"
  } else if (override == "pass") {
    keeper.status = "pass"
  } else if (override == "fail") {
    keeper.status = "fail"
  }
  
  
  
  #### Adjust for dilution with ZnAc ####
  data.spreadsheet <- data.spreadsheet %>%
    left_join(processing.metadata) %>%
    mutate(S_conc_uM = S_conc_uM * ((mass - tare)/(mass - tare - preservativeVol))) %>%
    select(sulfurID, S_conc_uM) %>%
    mutate(S_conc_uM = signif(S_conc_uM, 3)) %>%
    as.data.frame()
  
  
  #### Read out data files ####
  
  good_data_location = paste(good.data.location,
                             date.of.analysis,
                             "_sulfide.csv",
                             sep = "")
  bad_data_location = paste(bad.data.location,
                            date.of.analysis,
                            "_sulfide_bad.csv",
                            sep = "")
  
  if (run.info.vector["standard_curve"] == "high" && keeper.status == "pass") {
    
    print("Standard is high, and we're all good on QC.")
    
    data.spreadsheet.passing <- data.spreadsheet %>%
      mutate(in_curve = S_conc_uM > 40) %>% 
      filter(in_curve == TRUE) %>%
      select(-in_curve)
    write.csv(file = good_data_location,
              x = data.spreadsheet.passing,
              row.names = FALSE,
              quote = FALSE)
    
    if (dim(data.spreadsheet.passing)[1] != dim(data.spreadsheet)[1]) {
      data.spreadsheet.failing <- data.spreadsheet %>%
        mutate(in_curve = S_conc_uM > 40) %>% 
        filter(in_curve == FALSE) %>%
        select(-in_curve)
      write.csv(file = bad_data_location,
                x = data.spreadsheet.failing,
                row.names = FALSE,
                quote = FALSE)     
    }
  } else if (run.info.vector["standard_curve"]  == "low" && keeper.status == "pass") {
    
    print("Standard is low, and we're all good on QC.")
    
    data.spreadsheet.passing <- data.spreadsheet %>%
      mutate(in_curve = S_conc_uM < 40) %>% 
      filter(in_curve == TRUE) %>%
      select(-in_curve)
    write.csv(file = good_data_location,
              x = data.spreadsheet.passing,
              row.names = FALSE,
              quote = FALSE)
    
    if (dim(data.spreadsheet.passing)[1] != dim(data.spreadsheet)[1]) {
      data.spreadsheet.failing <- data.spreadsheet %>%
        mutate(in_curve = S_conc_uM < 40) %>% 
        filter(in_curve == FALSE) %>%
        select(-in_curve)
      write.csv(file = bad_data_location,
                x = data.spreadsheet.failing,
                row.names = FALSE,
                quote = FALSE)     
    }
  } else if (keeper.status == "fail") {
    
    print("QC SUCKS.")
    
    write.csv(file = bad_data_location,
              x = data.spreadsheet,
              row.names = FALSE,
              quote = FALSE)  
  }
  
  
  
  #### Read out report file ####
  
  report_file <- paste(report.location,
                           date.of.analysis,
                           "_sulfide_data_report.txt",
                           sep = "")
  
  fileConn <- file(report_file)
  writeLines(c(paste("Analysis run on: ",
                     date.of.analysis),
               paste(""),
               paste("Standard Curve Info"),
               paste("Slope of curve:", curve.coeff),
               paste("Intercept of curve:", inter.coeff),
               paste("R^2 of Curve:", adj.r.squared),
               paste("QC on R^2:", std.curve.qc),
               paste(""),
               paste(""),
               # paste("Limit of Detection:", LOD.mass),
               paste("Calculated from intercept of curve and standard error of regression"),
               paste(""),
               paste(""),
               paste("Blank value =", blank.value),
               paste("Blank SD =", blank.sd),
               paste(""),
               paste(""),
               paste("Check measurements"),
               paste(""),
               paste("CCB area:"),
               paste(CCB.values$Absorbance_667),
               paste(""),
               paste("CCV % deviation:"),
               paste(round(CCV.values$deviation, 2)),
               paste("QC on CCV:", ccv.qc),
               paste(""),
               paste("Matrix Spike recovery:", MS.recovery),
               paste(""),
               paste(""),
               paste("Replicate data:"),
               paste(rep.data),
               paste("QC on replicates:", rep.qc),
               paste("")),
             fileConn)
  close(fileConn)
  write.table(rep.data,
              file = report_file,
              append = TRUE,
              row.names = FALSE)
  sink(file = report_file,
       append = TRUE)
  cat("\nClean Data\n")
  sink()
  write.table(data.spreadsheet,
              file = report_file,
              append = TRUE,
              row.names = FALSE)
  
}


#### Process it ####

data_processing_function("dataRaw/waterChemistry/sulfide/sulfide_20191023.xlsx")

data_processing_function("dataRaw/waterChemistry/sulfide/sulfide_20191029.xlsx")
# We'll need to re-run this on the low curve. 




#### Clean up before combining all samples ####

rm(list = ls())


# Load up all incubation data
list.o.results <- list.files(path = "dataEdited/waterChemistry/sulfide",
                             pattern = "_sulfide.csv",
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


#### Combine all incubation data ####

# Read in incubation metadata
MA.metadata <- read.csv("metadata/processedMetadata/sulfide_MA.csv",
                        stringsAsFactors = FALSE)
incubation.results <- S.results %>%
  inner_join(MA.metadata)

write.csv(incubation.results,
          "dataEdited/waterChemistry/sulfide/MA_data.csv",
          row.names = FALSE,
          quote = FALSE)

#### Incubation samples to analyze yet ####
unanalyzed.incubation.samples <- MA.metadata %>%
  filter(!(sulfurID %in% incubation.results$sulfurID)) %>%
  select(sulfurID) %>%
  unlist(use.names = FALSE) %>%
  write(file = "protocols/sulfideCline_protocol/incubation_samples_to_analyze.txt")



#### Combine all WC data ####
WC.metadata <- read.csv("metadata/processedMetadata/sulfide_WC.csv",
                        stringsAsFactors = FALSE)
WC.results <- S.results %>%
  inner_join(WC.metadata)

write.csv(WC.results,
          "dataEdited/waterChemistry/sulfide/WC_data.csv",
          row.names = FALSE,
          quote = FALSE)