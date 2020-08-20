#### code/cleaning_scripts/clean_sulfide_data.R ####
# Benjamin D. Peterson



#### Get ready, and fly! ####

rm(list = ls())
setwd("~/Documents/research/BLiMMP/")
library(dplyr)
library(lubridate)
library(tidyr)
library(readxl)



#### Prepare metadata ####

S.metadata <- read_xlsx("metadata/chem_S.xlsx") %>%
  select(-notes)
sample_IDs <- read_xlsx("metadata/2_sample_IDs.xlsx")
incubation_IDs <- read_xlsx("metadata/4_MA_ID.xlsx") %>%
  select(-notes)
trip_IDs <- read_xlsx("metadata/1_trip_IDs.xlsx") %>%
  select(-notes)



#### Save out water column data ####

WC.metadata <- S.metadata %>%
  filter(!(sampleID == "NA")) %>%
  left_join(sample_IDs) %>%
  left_join(trip_IDs) %>%
  select(sulfurID, sampleID, tripID, depth, startDate)
write.csv(WC.metadata,
          file = "metadata/processedMetadata/sulfide_WC.csv",
          row.names = FALSE,
          quote = FALSE)



#### Save out incubation data ####

MA.metadata <- S.metadata %>%
  filter(!(incubationID == "NA")) %>%
  select(sulfurID, incubationID) %>%
  left_join(incubation_IDs) %>%
  left_join(sample_IDs) %>%
  left_join(trip_IDs) %>%
  select(sulfurID, incubationID, sampleID, tripID, depth, startDate, dateCollected, filtered, amendment)

# Add treatment column with all treatment information
filtered.vector <- c("filtered", "unfiltered")
names(filtered.vector) <- c("yes", "no")
MA.metadata <- MA.metadata %>%
  mutate(treatment = paste(filtered.vector[filtered],
                           amendment,
                           sep = "-")) %>%
  as.data.frame()
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
   incubation_IDs,
   sample_IDs,
   trip_IDs)



#### Define output locations ####
report.location <- "dataEdited/waterChemistry/sulfide/reports/"
bad.data.location <- "dataEdited/waterChemistry/sulfide/badData/"
good.data.location <- "dataEdited/waterChemistry/sulfide/dataForReview/"



#### Data processing function ####

#data_file <- "dataRaw/waterChemistry/sulfide/sulfide_20200818.xlsx"
#override = "pass"

# Define the function
data_processing_function <- function(data_file,
                                     override = "off",
                                     accept.shitty.curve = FALSE) {

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
    strsplit(".xlsx") %>%
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
  
  if (std.curve.qc == "pass" | accept.shitty.curve == TRUE) {
    
    print("Curve looks good")
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
    
    # If sample calculated to be less than 0µM, make it 0
    data.spreadsheet$S_conc_uM[which(data.spreadsheet$S_conc_uM < 0)] <- 0
    
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
    CCV.values
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
    if (MS.recovery > 120 | MS.recovery < 80) {
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
    if (any(abs(rep.data$rsd.rep) < 15 | is.na(rep.data$rsd.rep))) {
      rep.qc <- "pass"
    } else {
      rep.qc <- "fail"
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
        mutate(in_curve = S_conc_uM >= 20) %>% 
        filter(in_curve == TRUE) %>%
        select(-in_curve)
      write.csv(file = good_data_location,
                x = data.spreadsheet.passing,
                row.names = FALSE,
                quote = FALSE)
      
      if (dim(data.spreadsheet.passing)[1] != dim(data.spreadsheet)[1]) {
        data.spreadsheet.failing <- data.spreadsheet %>%
          mutate(in_curve = S_conc_uM >= 20) %>% 
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
        mutate(in_curve = S_conc_uM <= 50) %>% 
        filter(in_curve == TRUE) %>%
        select(-in_curve)
      write.csv(file = good_data_location,
                x = data.spreadsheet.passing,
                row.names = FALSE,
                quote = FALSE)
      
      if (dim(data.spreadsheet.passing)[1] != dim(data.spreadsheet)[1]) {
        data.spreadsheet.failing <- data.spreadsheet %>%
          mutate(in_curve = S_conc_uM <= 50) %>% 
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
  } else {
    print("Curve is shit. You ain't going anywhere")
  }
 
}


#### Process data ####

data_processing_function("dataRaw/waterChemistry/sulfide/sulfide_20191023.xlsx")
# This one looks good.
# Copy it into the good_data folder

data_processing_function("dataRaw/waterChemistry/sulfide/sulfide_20191029.xlsx",
                         override = "fail")
# We'll need to re-run this on the low curve. 

data_processing_function("dataRaw/waterChemistry/sulfide/sulfide_20191105.xlsx")
# Bad curve. Will need re-running.

data_processing_function("dataRaw/waterChemistry/sulfide/sulfide_20191119.xlsx")
# This one looks good.
# Copy it into the good_data folder

data_processing_function("dataRaw/waterChemistry/sulfide/sulfide_20191126.xlsx")
# This one looks good. 
# Copy it into the good_data folder

data_processing_function("dataRaw/waterChemistry/sulfide/sulfide_20191203.xlsx",
                         override = "pass")
# The MS QC failed on this one, it was too high. The original measured curve was 
# pretty low (Ab667 of 0.069). Let's let this one pass, since everything else
# looks so good. 
# Copy it into the good_data folder
# Anna zeroed the machine with 1% ZnAc instead of water, which shifted the raw 
# values slightly. However, it does not seem to have negatively impacted the 
# data processing, as we have a good standard curve. We kept this anyways. 

data_processing_function("dataRaw/waterChemistry/sulfide/sulfide_20191209.xlsx")
# This one looks good. 
# All the data is within the curve. 
# BLiMMP_S_0080 is a field blank that I set up a while back. 
# It was an unwashed 30ml HDPE bottle that we filled with 
# nanopore water and ZnAc. No detectable sulfide in it. Good!

data_processing_function("dataRaw/waterChemistry/sulfide/sulfide_20191209B.xlsx")
# Standard curve looks like shit. This needs to be re-run. 

data_processing_function("dataRaw/waterChemistry/sulfide/sulfide_20191210.xlsx")
# BLiMMP_S_0067 is below the curve (is 0) so was removed automatically.
# However, we'll just include it, and will manually add it to the file 
# that will be transferred to the good_data folder. 


data_processing_function("dataRaw/waterChemistry/sulfide/sulfide_20200818.xlsx")
# Curve looks good
# Two samples were below 20µM and will need to be run on the low curve
# when I get a chance. Other than that, looks great.


#### Clean up before combining all samples ####

rm(list = ls())


# Load up all incubation data
list.o.results <- list.files(path = "dataEdited/waterChemistry/sulfide/goodData/",
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
MA.results <- S.results %>%
  inner_join(MA.metadata) %>%
  arrange(sulfurID)

write.csv(MA.results,
          "dataEdited/waterChemistry/sulfide/MA_data.csv",
          row.names = FALSE,
          quote = FALSE)



#### Combine all WC data ####

WC.metadata <- read.csv("metadata/processedMetadata/sulfide_WC.csv",
                        stringsAsFactors = FALSE)
WC.results <- S.results %>%
  inner_join(WC.metadata) %>%
  arrange(sulfurID)

write.csv(WC.results,
          "dataEdited/waterChemistry/sulfide/WC_data.csv",
          row.names = FALSE,
          quote = FALSE)



#### Samples to analyze yet ####

# Make a vector with unanalyzed samples that 
# will not be analyzed
samples.to.skip <- c("BLiMMP_INC_S_007",
                     "BLiMMP_INC_S_008",
                     "BLiMMP_INC_S_009")

unanalyzed.MA.samples <- MA.metadata %>%
  filter(!(sulfurID %in% MA.results$sulfurID)) %>%
  filter(!(sulfurID %in% samples.to.skip))
write.csv(unanalyzed.MA.samples,
          file = "protocols/sulfideCline_protocol/MA_incubation_samples_to_analyze.csv",
          row.names = FALSE)

unanalyzed.WC.samples <- WC.metadata %>%
  filter(!(sulfurID %in% WC.results$sulfurID)) %>%
  filter(!(sulfurID %in% samples.to.skip))
write.csv(unanalyzed.WC.samples,
          file = "protocols/sulfideCline_protocol/WC_incubation_samples_to_analyze.csv",
          row.names = FALSE)

