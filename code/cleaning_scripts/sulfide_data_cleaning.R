#### code/cleaning_scripts/clean_sulfide_data.R ####
# Benjamin D. Peterson


#### Get ready, and fly! ####
rm(list = ls())
setwd("~/Documents/research/BLiMMP/")
library(dplyr)
library(lubridate)
library(tidyverse)
library(readxl)


#### Prepare metadata ####
S.metadata <- read_xlsx("metadata/raw_metadata/chem_S.xlsx",
                        sheet = "Sheet1") %>%
  select(-notes) %>%
  filter(!is.na(sampleID) & !is.na(incubationID))


#### Generate needed processing metadata ####
processing.metadata <- S.metadata %>%
  select(sulfurID, mass, tare, preservativeVol) %>%
  mutate_at(.vars = c("mass", "tare", "preservativeVol"),
            as.numeric)
# rm(S.metadata)


#### Define output locations ####
report.location <- "dataEdited/waterChemistry/sulfide/reports/"
bad.data.location <- "dataEdited/waterChemistry/sulfide/badData/"
good.data.location <- "dataEdited/waterChemistry/sulfide/dataForReview/"


#### Data processing function ####
# data_file <- "dataRaw/waterChemistry/sulfide/sulfide_20211123a.xlsx"
# override = "off"
# accept.shitty.curve = FALSE
# remove.these.samples = c("H_1")
# output.file.name = NULL
# Define the function
data_processing_function <- function(data_file,
                                     override = "off",
                                     accept.shitty.curve = FALSE,
                                     remove.these.samples = NULL,
                                     output.file.name = NULL,
                                     high_curve_cutoff = 20,
                                     low_curve_cutoff = 30) {

  keeper.status <- "pass"
  
  #### Read in needed data ####
  
  # Read in raw data
  data.spreadsheet.raw <- read_xlsx(data_file,
                                    sheet = "raw_results") %>%
    filter(!is.na(sulfurID)) %>%
    mutate(Absorbance_667 = Absorbance_667 / dilutionFactor) %>%
    as.data.frame()
  
  if (length(remove.these.samples) > 0) {
    data.spreadsheet.raw <- data.spreadsheet.raw %>%
      filter(!(sulfurID %in% remove.these.samples))
    print(paste("Removing samples", remove.these.samples))
  }
  
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
  
  if (accept.shitty.curve == TRUE) {
    std.curve.qc <- "crappy_but_we_passing"
  }
  
  # Read out image of curve
  if (is.null(output.file.name)) {
    pdf.of.curve <- paste(report.location,
                          date.of.analysis,
                          "_curve.pdf",
                          sep = "")
  } else {
    pdf.of.curve <- paste(report.location,
                          output.file.name,
                          "_curve.pdf",
                          sep = "")
  }
  
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
  
  if (std.curve.qc == "pass" | accept.shitty.curve == TRUE) {
    
    print("Curve looks good")
    
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
    
    #### Calculate detection limit ####
    
    LOD.fluorescence <- inter.coeff + 3*summary(std.curve.regression)$coefficients[2,2]
    LOD.mass <- (LOD.fluorescence - std.curve.regression$coefficients[1]) / std.curve.regression$coefficients[2]
    
    
    #### QC: Ensure CCVs are within 10% ####
    
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
    
    # Check MS is within 20%
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
    data.spreadsheet <- processing.metadata %>%
      right_join(data.spreadsheet) %>%
      mutate(S_conc_uM = S_conc_uM * ((mass - tare)/(mass - tare - preservativeVol))) %>%
      select(sulfurID, S_conc_uM) %>%
      mutate(S_conc_uM = signif(S_conc_uM, 3)) %>%
      as.data.frame()
    
    
    #### Read out data files ####
    
    if (is.null(output.file.name)) {
      good_data_location = paste(good.data.location,
                                 date.of.analysis,
                                 "_sulfide.csv",
                                 sep = "")
      bad_data_location = paste(bad.data.location,
                                date.of.analysis,
                                "_sulfide_bad.csv",
                                sep = "")
      
      report_file = paste(report.location,
                          date.of.analysis,
                          "_sulfide_data_report.txt",
                          sep = "")
    } else {
      print(paste("Supplied the output file name prefix as",
                  output.file.name))
      good_data_location = paste(good.data.location,
                                 output.file.name,
                                 "_sulfide.csv",
                                 sep = "")
      bad_data_location = paste(bad.data.location,
                                output.file.name,
                                "_sulfide_bad.csv",
                                sep = "")
      report_file = paste(report.location,
                          output.file.name,
                          "_sulfide_data_report.txt",
                          sep = "")
    }
    
    #### Check if samples are within the cut-offs we set ####
    if (run.info.vector["standard_curve"] == "high" && keeper.status == "pass") {
      
      print("Standard is high, and we're all good on QC.")
      
      data.spreadsheet.passing <- data.spreadsheet %>%
        mutate(in_curve = S_conc_uM >= high_curve_cutoff) %>% 
        filter(in_curve == TRUE) %>%
        select(-in_curve)
      write.csv(file = good_data_location,
                x = data.spreadsheet.passing,
                row.names = FALSE,
                quote = FALSE)
      
      if (dim(data.spreadsheet.passing)[1] != dim(data.spreadsheet)[1]) {
        data.spreadsheet.failing <- data.spreadsheet %>%
          mutate(in_curve = S_conc_uM >= high_curve_cutoff) %>% 
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
        mutate(in_curve = S_conc_uM <= low_curve_cutoff) %>% 
        filter(in_curve == TRUE) %>%
        select(-in_curve)
      write.csv(file = good_data_location,
                x = data.spreadsheet.passing,
                row.names = FALSE,
                quote = FALSE)
      
      if (dim(data.spreadsheet.passing)[1] != dim(data.spreadsheet)[1]) {
        data.spreadsheet.failing <- data.spreadsheet %>%
          mutate(in_curve = S_conc_uM <= low_curve_cutoff) %>% 
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
                 paste("Limit of Detection:", LOD.mass),
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
                 paste("QC on replicates:", rep.qc)),
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

# data_processing_function("dataRaw/waterChemistry/sulfide/sulfide_20191023.xlsx")
# # This one looks good.
# # Copy it into the good_data folder
# 
# data_processing_function("dataRaw/waterChemistry/sulfide/sulfide_20191029.xlsx",
#                          override = "fail")
# # We'll need to re-run this on the low curve. 
# 
# data_processing_function("dataRaw/waterChemistry/sulfide/sulfide_20191105.xlsx")
# # Bad curve. Will need re-running.
# 
# data_processing_function("dataRaw/waterChemistry/sulfide/sulfide_20191119.xlsx")
# # This one looks good.
# # Copy it into the good_data folder
# 
# data_processing_function("dataRaw/waterChemistry/sulfide/sulfide_20191126.xlsx")
# # This one looks good. 
# # Copy it into the good_data folder
# 
# data_processing_function("dataRaw/waterChemistry/sulfide/sulfide_20191203.xlsx",
#                          override = "pass")
# # The MS QC failed on this one, it was too high. The original measured curve was 
# # pretty low (Ab667 of 0.069). Let's let this one pass, since everything else
# # looks so good. 
# # Copy it into the good_data folder
# # Anna zeroed the machine with 1% ZnAc instead of water, which shifted the raw 
# # values slightly. However, it does not seem to have negatively impacted the 
# # data processing, as we have a good standard curve. We kept this anyways. 
# 
# data_processing_function("dataRaw/waterChemistry/sulfide/sulfide_20191209.xlsx")
# # This one looks good. 
# # All the data is within the curve. 
# # BLiMMP_S_0080 is a field blank that I set up a while back. 
# # It was an unwashed 30ml HDPE bottle that we filled with 
# # nanopore water and ZnAc. No detectable sulfide in it. Good!
# 
# data_processing_function("dataRaw/waterChemistry/sulfide/sulfide_20191209B.xlsx")
# Standard curve looks like shit. This needs to be re-run.
# 
# data_processing_function("dataRaw/waterChemistry/sulfide/sulfide_20191210.xlsx")
# # BLiMMP_S_0067 is below the curve (is 0) so was removed automatically.
# # However, we'll just include it, and will manually add it to the file 
# # that will be transferred to the good_data folder. 
# 
# 
# data_processing_function("dataRaw/waterChemistry/sulfide/sulfide_20200818.xlsx")
# # Curve looks good
# # Two samples were below 20µM and will need to be run on the low curve
# # when I get a chance. Other than that, looks great.
# 
# data_processing_function("dataRaw/waterChemistry/sulfide/sulfide_20200909.xlsx",
#                          remove.these.samples = "H_2")
# 
# data_processing_function(data_file = "dataRaw/waterChemistry/sulfide/sulfide_20200914a.xlsx",
#                          output.file.name = "2020-09-14a",
#                          remove.these.samples = "H_1")
# 
# data_processing_function(data_file = "dataRaw/waterChemistry/sulfide/sulfide_20200914b.xlsx",
#                          output.file.name = "2020-09-14b")
# # Looking good here!
#  
# data_processing_function(data_file = "dataRaw/waterChemistry/sulfide/sulfide_20200917a.xlsx",
#                          output.file.name = "2020-09-17a",
#                          override = "pass")
# # The replicates are failing here because the samples are so low.
# # We're gonna go ahead and pass it.
# 
# data_processing_function(data_file = "dataRaw/waterChemistry/sulfide/sulfide_20200917b.xlsx",
#                          output.file.name = "2020-09-17b")
# # This looks good
# 
# data_processing_function(data_file = "dataRaw/waterChemistry/sulfide/sulfide_20200917c.xlsx",
#                          output.file.name = "2020-09-17c")
# # This looks good
# 
# data_processing_function(data_file = "dataRaw/waterChemistry/sulfide/sulfide_20201019.xlsx",
#                          output.file.name = "2020-10-19",
#                          remove.these.samples = "H_1")
# # Starting to look like H_1 is not linear with the rest of the bunch in the high curve.
# # Wonder if it's saturating at that high of an absorbance.
# # Let's refrain from including this one in the future. 
# 
# data_processing_function(data_file = "dataRaw/waterChemistry/sulfide/sulfide_20201021a.xlsx",
#                          output.file.name = "2020-10-21a",
#                          remove.these.samples = "H_1")
# # This one looks good. Kept H_1 out as decided before.
# # Two samples below curve, will need to rerun those on low curve.
# 
# data_processing_function(data_file = "dataRaw/waterChemistry/sulfide/sulfide_20201021b.xlsx",
#                          output.file.name = "2020-10-21b",
#                          remove.these.samples = "H_1")
# # Looks good! All samples within curve.
# 
# data_processing_function(data_file = "dataRaw/waterChemistry/sulfide/sulfide_20201021c.xlsx",
#                          output.file.name = "2020-10-21c",
#                          remove.these.samples = "H_1",
#                          override = "pass")
# # CCV is -10.17%. I'm going to count that as a pass, since the other checks look good.
# # These are all incubation samples anyways.
# 
# 
# 
# data_processing_function(data_file = "dataRaw/waterChemistry/sulfide/sulfide_20201027a.xlsx",
#                          output.file.name = "2020-10-27a")
# # Passes with flying colors
# 
# data_processing_function(data_file = "dataRaw/waterChemistry/sulfide/sulfide_20201027b.xlsx",
#                          output.file.name = "2020-10-27b",
#                          remove.these.samples = "L_1")
# # Used L_3 for CCV, had to manually change the sulfurID
# # in the raw data file.
# # This run also looks great.
# 
# data_processing_function(data_file = "dataRaw/waterChemistry/sulfide/sulfide_20201027c.xlsx",
#                          output.file.name = "2020-10-27c",
#                          remove.these.samples = c("H_1"),
#                          override = "pass")
# # Forgot to add matrix spike to this sample...
# # Other QC looks good though, so we'll pass it.
# 
# data_processing_function(data_file = "dataRaw/waterChemistry/sulfide/sulfide_20201110c.xlsx",
#                          output.file.name = "2020-11-10c",
#                          remove.these.samples = c("L_1"))
# # Looks good here
# 
# data_processing_function(data_file = "dataRaw/waterChemistry/sulfide/sulfide_20201110a.xlsx",
#                          output.file.name = "2020-11-10a",
#                          remove.these.samples = c("H_1"))
# # Looks good here



# data_processing_function(data_file = "dataRaw/waterChemistry/sulfide/sulfide_20211110a.xlsx",
#                          output.file.name = "2021-11-10a",
#                          remove.these.samples = c("H_1", "H_2"))
# Not great, will redo with new standards

# data_processing_function(data_file = "dataRaw/waterChemistry/sulfide/sulfide_20211110b.xlsx",
#                          output.file.name = "2021-11-10b",
#                          remove.these.samples = c("L_1"))
# Not great, will redo with new standards

# data_processing_function(data_file = "dataRaw/waterChemistry/sulfide/sulfide_20211115a.xlsx",
#                          output.file.name = "2021-11-15a")
# Looks good! I fucked up the QCS on this one, added 100µM standard instead of 1mM, so that's not great.
# Still passed with it though, so we'll go with it. Not sure the matrix spike is too important for this analysis anyways.
# Copy to goodData.
# data_processing_function(data_file = "dataRaw/waterChemistry/sulfide/sulfide_20211115b.xlsx",
#                          output.file.name = "2021-11-15b",
#                          override = "pass")
# Looks good! Only blemish is the matrix spike. For this, I accidently added the 1 mM standard,
# rather than the 100 µM. I think this is too high for this Cline's reagent, so just going to ignore it.
# Copy to goodData.

# data_processing_function(data_file = "dataRaw/waterChemistry/sulfide/sulfide_20211117a.xlsx",
#                          output.file.name = "2021-11-17a")
# Looks good enough!
# Copy to goodData

# data_processing_function(data_file = "dataRaw/waterChemistry/sulfide/sulfide_20211117b.xlsx",
#                          output.file.name = "2021-11-17b")
# # Good, copy to goodData
# 
# data_processing_function(data_file = "dataRaw/waterChemistry/sulfide/sulfide_20211117c.xlsx",
#                          output.file.name = "2021-11-17c")
# # Good, copy to goodData
# 
# data_processing_function(data_file = "dataRaw/waterChemistry/sulfide/sulfide_20211118a.xlsx",
#                          output.file.name = "2021-11-18a")
# # Good, copy to goodData
# 
# data_processing_function(data_file = "dataRaw/waterChemistry/sulfide/sulfide_20211118b.xlsx",
#                          output.file.name = "2021-11-18b")
# # Good, copy to goodData

# data_processing_function(data_file = "dataRaw/waterChemistry/sulfide/sulfide_20211118c.xlsx",
#                          output.file.name = "2021-11-18c",
#                          remove.these.samples = c("H_1", "H_6"))
# H_1 looked terrible. But, the samples weren't that high on the curve, and the other standards
# looked good and matched up with previous results. So, we'll remove that standard and continue
# on. Remake H_1 before continuing more analyses.
# The 11.9 m samples, from the t0 point, were all below 20.
# Also going to remove the 20 µM standard, after seeing the ones below, I'm not sure that we should
# include that low of a standard for the high curve... Seems like we're losing linearity.
# Good, copy to goodData

# data_processing_function(data_file = "dataRaw/waterChemistry/sulfide/sulfide_20211119a.xlsx",
#                          output.file.name = "2021-11-19a",
#                          remove.these.samples = c("H_6"))
# CCV has 10.5% RSD. But, I think this is due to pipetting, not drift, since I reran the original
# H_4 sample and it was close to the value measured first. I think we'll pass this, since everything
# else looks good.
# Removing the low standard, since I think that will be best after seeing other data. Seems to start to
# lose linearity less than 0.2 absorbance

# data_processing_function(data_file = "dataRaw/waterChemistry/sulfide/sulfide_20211119b.xlsx",
#                          output.file.name = "2021-11-19b",
#                          remove.these.samples = c("H_1", "H_6"))
# Pulled out the top and bottom samples. Ah, these standards don't look bad. Nick prepped this one, need to make 
# pipetting is going carefully here.


# data_processing_function(data_file = "dataRaw/waterChemistry/sulfide/sulfide_20211123a.xlsx",
#                          output.file.name = "2021-11-23a",
#                          override = "pass")
# These look great except for the matrix spike... We'll stick with it.

# data_processing_function(data_file = "dataRaw/waterChemistry/sulfide/sulfide_20211123b.xlsx",
#                          output.file.name = "2021-11-23b")
# All looks great! Keep it.
# Used curve from sulfide_20211123a.xlsx to stay consistent with all high curve samples.

# data_processing_function(data_file = "dataRaw/waterChemistry/sulfide/sulfide_20211123c.xlsx",
#                          output.file.name = "2021-11-23c")
# Used curve from sulfide_20211123a.xlsx to stay consistent with all high curve samples.

# data_processing_function(data_file = "dataRaw/waterChemistry/sulfide/sulfide_20211123d.xlsx",
#                          output.file.name = "2021-11-23d")
# Looks great! Used curve from sulfide_20211123a.xlsx to stay consistent with all high curve samples.
# Only made difference of 4 µM.

# data_processing_function(data_file = "dataRaw/waterChemistry/sulfide/sulfide_20211123e.xlsx",
#                          output.file.name = "2021-11-23e")
# PERFECT

# data_processing_function(data_file = "dataRaw/waterChemistry/sulfide/sulfide_20211124a.xlsx",
#                          output.file.name = "2021-11-24a",
#                          override = "pass")
# Spike recovery way off, but also used the wrong spike. Everything else looks spot on, so
# we'll accept it.
# 
# data_processing_function(data_file = "dataRaw/waterChemistry/sulfide/sulfide_20211124b.xlsx",
#                          output.file.name = "2021-11-24b",
#                          override = "pass")
# CCV is a little off, but we also used the second lowest standard, which wasn't a good idea...
# too much variation. Also used the wrong spike here. The data looks good overall, so we'll
# stick with it.

# data_processing_function(data_file = "dataRaw/waterChemistry/sulfide/sulfide_20211124c.xlsx",
#                          output.file.name = "2021-11-24c",
#                          override = "pass")
# Wrong spike. Looks great other than that.
# The spikes were super consistent across the three runs from 2021-11-24.


 

#### Clean up before combining all samples ####

rm(list = ls())


# Load up all incubation data
list.o.results <- list.files(path = "dataEdited/waterChemistry/sulfide/goodData",
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
# Read in timing information (from Hg incubation data)
timing.data <- read.csv("metadata/processedMetadata/incubation_metadata.csv") %>%
  spread(key = constituent, value = bottleID) %>%
  select(incubationID, t, durationInDays) %>%
  rename(incubationTimePoint = t)
# Combine data
MA.results <- S.results %>%
  inner_join(MA.metadata) %>%
  left_join(timing.data) %>%
  arrange(sulfurID)

# Save out data
write.csv(MA.results,
          "dataEdited/incubation_sulfide_data.csv",
          row.names = FALSE,
          quote = FALSE)



#### Combine all WC data ####

WC.metadata <- read.csv("metadata/processedMetadata/sulfide_WC.csv",
                        stringsAsFactors = FALSE)


#### Calculate depth of corewater samples ####
waterDepth = 24
WC.results <- S.results %>%
  inner_join(WC.metadata) %>%
  arrange(sulfurID) %>%
  mutate(depthOriginal = depth) %>%
  mutate(corewater = grepl(pattern = "-",
                           x = depthOriginal))
WC.results[WC.results$corewater, ] <- WC.results[WC.results$corewater, ] %>%
  mutate(depth = paste("-",
                       strsplit(depth, "-") %>% sapply("[", 2),
                       sep = "") %>%
           as.numeric() / 100) %>%
  mutate(depth = depth + (waterDepth * corewater))


#### Add replicate info ####
WC.results <- WC.results %>%
  group_by(startDate, depth) %>%
  mutate(replicate = row_number()) %>%
  ungroup() %>%
  as.data.frame()


#### Save out data ####
write.csv(WC.results,
          "dataEdited/waterChemistry/sulfide/WC_data.csv",
          row.names = FALSE,
          quote = FALSE)



#### Samples to analyze yet ####

# Make a vector with unanalyzed samples that 
# will not be analyzed
samples.to.skip <- c("BLiMMP_INC_S_007",
                     "BLiMMP_INC_S_008",
                     "BLiMMP_INC_S_009",
                     "BLI20_TS_001",
                     "BLI20_TS_002",
                     "BLI20_TS_003",
                     "BLI20_TS_043",
                     "BLI20_TS_044",
                     "BLI20_TS_045",
                     "BLI21_TS_005")

unanalyzed.MA.samples <- MA.metadata %>%
  filter(!(sulfurID %in% MA.results$sulfurID)) %>%
  filter(!(sulfurID %in% samples.to.skip)) %>%
  arrange(sulfurID)
write.csv(unanalyzed.MA.samples,
          file = "dataEdited/waterChemistry/sulfide/MA_incubation_samples_to_analyze.csv",
          row.names = FALSE)

unanalyzed.WC.samples <- WC.metadata %>%
  filter(!(sulfurID %in% WC.results$sulfurID)) %>%
  filter(!(sulfurID %in% samples.to.skip)) %>%
  arrange(sulfurID)
write.csv(unanalyzed.WC.samples,
          file = "dataEdited/waterChemistry/sulfide/WC_incubation_samples_to_analyze.csv",
          row.names = FALSE)

