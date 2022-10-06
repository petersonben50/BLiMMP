#### code/BGC_profiles/all_profiles.R ####
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("~/Documents/research/BLiMMP/")
library(lubridate)
library(readxl)
library(tidyverse)
source("code/BLiMMP_functions.R")
detection.vector <- c(4, 16)


#### Read in data ####
sulfide.data <- read.csv("dataEdited/waterChemistry/sulfide/WC_data.csv")
sulfate.data <- read.csv("dataEdited/waterChemistry/sulfate/WC_data.csv")
metals.data <- readRDS("dataEdited/waterChemistry/ICP/2020_2021_Fe_Mn_data.rds")
Hg.data <- read.csv("dataEdited/Hg/Hg_data_clean.csv")



#### Define function: Plot exo data ####
plot.exo.data <- function(date.of.sampling,
                          points.of.sampling = NULL,
                          main.title = NULL,
                          legend.location = "topright") {
  # Read in exo data
  exo.data.file.name <- paste("dataEdited/exo/",
                              date.of.sampling,
                              "_profile.csv",
                              sep = "")
  exo.data <- read.csv(exo.data.file.name,
                       stringsAsFactors = FALSE)

  # Plot the temp points and set up graph
  plot(x = (exo.data$Temp_C-5)/2.5,
       y = -exo.data$depth,
       type = "l",
       lwd = 3,
       xaxt = "n",
       xlab = "",
       xlim = c(0, 10),
       yaxt = "n",
       ylab = "",
       ylim = c(-25, 2),
       col = cb.translator["blue"])
  # Plot DO
  DO.fudge.factor <- 15
  lines(x = exo.data$ODO_sat/DO.fudge.factor,
        y = -exo.data$depth,
        col = cb.translator["black"],
        lwd = 3,
        lty = 4)
  # Plot the turbidity values
  turb.fudge.factor <- 0.5
  lines(x = exo.data$Turbidity_FNU/turb.fudge.factor,
        y = -exo.data$depth,
        col = cb.translator["reddishpurple"],
        lwd = 4,
        lty = 3)
  # Add y-axis with depth measurements
  axis(2,
       at = seq(0, -25, by = -5),
       labels = seq(0, 25, by = 5))
  # Add axis for temperature
  axis(1,
       at = seq(0, 10, by = 2),
       labels = seq(0, 10, by = 2)*3)
  # Add axis for DO values
  axis(1,
       at = seq(0, 10, by = 2),
       labels = seq(0, 10, by = 2)*DO.fudge.factor,
       line = 3,
       par(mgp = c(3, 0.5, 0)))
  # Add axis for turbidity values
  axis(1,
       at = seq(0, 10, by = 2),
       labels = seq(0, 10, by = 2)*turb.fudge.factor,
       line = 6,
       par(mgp = c(3, 0.5, 0)))
  
  
  # Add titles for depth
  title(xlab = "Temp (C)",
        ylab = "Depth (m)",
        line = 1.5)
  # Add title for dissolved oxygen
  title(xlab = "ODO (% sat)",
        line = 4.5)
  # Add title for turbidity
  title(xlab = "Turbidity (FNU)",
        line = 7.5)
  
  # if (is.null(main.title)) {
  #   title(main = date.of.sampling)
  # } else {
  #   title(main = main.title)
  # }
  # 
  
  # OPTIONAL: Add in sampling points
  if (!is.null(points.of.sampling)) {
    points(x = rep(9, length(points.of.sampling)),
           y = -points.of.sampling,
           pch = 16,
           cex = 2,
           col = cb.translator["bluishgreen"])
    legend(legend.location,
           legend = c("Temperature", "DO", "Turbidity", "Incubations"),
           col = c(cb.translator["blue"], cb.translator["black"], cb.translator["reddishpurple"], cb.translator["bluishgreen"]),
           lwd = c(3, 3, 3, NA),
           lty = c(1, 4, 3, NA),
           pt.cex = c(NA, NA, NA, 1.75),
           pch = c(NA, NA, NA, 16),
           bg = "white")
  } else {
    legend(legend.location,
           legend = c("Temperature", "DO", "Turbidity"),
           col = c(cb.translator["blue"], cb.translator["black"], cb.translator["reddishpurple"]),
           lwd = 3,
           lty = c(1, 4, 3, NA),
           pch = NA,
           bg = "white")
  }
}


#### Define function: Plot metals data ####
plot.metals.data <- function(icp.data.to.use,
                             constituent.of.interest,
                             color.vector.to.use,
                             point.vector.to.use,
                             date.of.interest,
                             xlim.to.use = c(0, 0.3),
                             xlab.to.use = "Metals (ppm)",
                             legend.location = "topright",
                             say.anything.for.two.legends = NULL) {
  
  icp.data.to.use.filtered <- icp.data.to.use %>%
    filter(date == date.of.interest)
  
  # Metals data
  plot(x = icp.data.to.use.filtered$concentration,
       y = as.numeric(icp.data.to.use.filtered$depth),
       ylim = c(25, -2),
       xlim = xlim.to.use,
       col = color.vector.to.use[icp.data.to.use.filtered$constituent],
       pch = point.vector.to.use[icp.data.to.use.filtered$state],
       xlab = "",
       ylab = "")
  
  for (metal in names(color.vector.to.use)) {
    line.data.p <- icp.data.to.use.filtered %>%
      filter(constituent == metal,
             state == "particulate") %>%
      arrange(depth)
    lines(x = line.data.p$concentration,
          y = line.data.p$depth,
          col = color.vector.to.use[line.data.p$constituent],
          lty = 2)
    line.data.f <- icp.data.to.use.filtered %>%
      filter(constituent == metal,
             state == "dissolved") %>%
      arrange(depth)
    lines(x = line.data.f$concentration,
          y = line.data.f$depth,
          col = color.vector.to.use[line.data.p$constituent])
  }
  
  if (is.null(say.anything.for.two.legends)) {
    legend(legend.location,
           legend = c(names(color.vector.to.use), names(point.vector.to.use)),
           pch = c(rep(15, length(color.vector.to.use)), point.vector.to.use),
           lty = c(rep(NA, length(color.vector.to.use)), 1, 2),
           col= color.vector.to.use, rep("darkgrey", length(point.vector.to.use)))
    
  } else {
    
    legend(x = 0.125,
           y = 3.5,
           legend = c("Mn (ppm)",
                      "Fe (ppm)"),
           pch = rep(15, length(color.vector.to.use)),
           col= color.vector.to.use,
           bg = "white")
    legend("topright",
           legend = c("Particulate",
                      "Dissolved"),
           pch = point.vector.to.use,
           lty = c(1, 2),
           col= rep("darkgrey", length(point.vector.to.use)),
           bg = "white")
   
  }

  title(xlab = xlab.to.use,
        ylab = "Depth (m)",
        line = 1.5)
}


#### Define function: Plot sulfide and sulfate data ####
plot.sulfide.sulfate.data <- function(sulfide.data.to.use,
                                      sulfate.data.to.use,
                                      xlim.to.use = c(0, 250),
                                      legend.location = "topright") {
  # Sulfide data
  plot(x = sulfide.data.to.use$S_conc_uM,
       y = sulfide.data.to.use$depth,
       pch = 16,
       ylim = c(25, -2),
       xlim = xlim.to.use,
       col = cb.translator["blue"],
       xlab = "",
       ylab = "")
  sulfide.data.to.use.lines <- sulfide.data.to.use %>%
    group_by(depth) %>%
    summarise(S_conc_uM = mean(S_conc_uM)) %>%
    arrange(depth)
  lines(x = sulfide.data.to.use.lines$S_conc_uM,
        y = sulfide.data.to.use.lines$depth,
        col = cb.translator["blue"])
  
  # Sulfate data
  points(x = sulfate.data.to.use$sulfate_uM,
         y = sulfate.data.to.use$depth,
         pch = 1,
         col = cb.translator["blue"])
  sulfate.data.to.use.lines <- sulfate.data.to.use %>%
    group_by(depth) %>%
    summarise(sulfate_uM = mean(sulfate_uM)) %>%
    arrange(depth)
  lines(x = sulfate.data.to.use.lines$sulfate_uM,
        y = sulfate.data.to.use.lines$depth,
        col = cb.translator["blue"],
        lty = 2)
  
  legend(legend.location,
         legend = c("Sulfate (µM)", "Sulfide (µM)"),
         pch = c(1, 16),
         lty = c(1, 2),
         col = cb.translator["blue"],
         bg = "white")
  title(xlab = "Sulfide/sulfate (µM)",
        ylab = "Depth (m)",
        line = 1.5)
  
  }


#### Function to plot Hg data ####
plot.Hg.profile <- function(date.to.use = "2021-09-10",
                            Hg.data.to.use = Hg.data,
                            phase = "filtered") {
  
  # Prepare Hg data
  Hg.data.date <- Hg.data.to.use %>%
    filter(sampleDate == date.to.use)
  
  if (phase == "filtered") {
    # Split out values
    FMHG.values <- Hg.data.date %>%
      filter(constituent == "FMHG_NG.L") %>%
      arrange(depth) %>%
      select(depth, concentration)
    FMHG.values.no.dup <- FMHG.values %>%
      group_by(depth) %>%
      summarise(concentration = mean(concentration))
    FiHG.values <- Hg.data.date %>%
      filter(constituent == "FiHg_NG.L") %>%
      arrange(depth) %>%
      select(depth, concentration)
    FiHG.values.no.dup <- FiHG.values %>%
      group_by(depth) %>%
      summarise(concentration = mean(concentration))
    PerFMHG.values <- Hg.data.date %>%
      filter(constituent == c("perFMHG"))
    PerFMHG.values.no.dup <- PerFMHG.values %>%
      group_by(depth) %>%
      summarise(concentration = mean(concentration))
  }
  
  if (phase == "filtered") {
    Hg.adjustment.factor <- 1
    Hg.label <- "Dissolved Hg (ng/L)"
    # Split out values
    MHG.values <- Hg.data.date %>%
      filter(constituent == "FMHG_NG.L") %>%
      arrange(depth) %>%
      select(depth, concentration)
    MHG.values.no.dup <- MHG.values %>%
      group_by(depth) %>%
      summarise(concentration = mean(concentration))
    iHG.values <- Hg.data.date %>%
      filter(constituent == "FiHg_NG.L") %>%
      arrange(depth) %>%
      select(depth, concentration)
    iHG.values.no.dup <- iHG.values %>%
      group_by(depth) %>%
      summarise(concentration = mean(concentration))
    PerMHG.values <- Hg.data.date %>%
      filter(constituent == c("perFMHG"))
    PerMHG.values.no.dup <- PerMHG.values %>%
      group_by(depth) %>%
      summarise(concentration = mean(concentration))
  } else if (phase == "particulate") {
    # Split out values
    Hg.adjustment.factor <- 1
    Hg.label <- "Particulate Hg (ng/L)"
    MHG.values <- Hg.data.date %>%
      filter(constituent == "PMHG_NG.L") %>%
      arrange(depth) %>%
      select(depth, concentration)
    MHG.values.no.dup <- MHG.values %>%
      group_by(depth) %>%
      summarise(concentration = mean(concentration))
    iHG.values <- Hg.data.date %>%
      filter(constituent == "PiHg_NG.L") %>%
      arrange(depth) %>%
      select(depth, concentration)
    iHG.values.no.dup <- iHG.values %>%
      group_by(depth) %>%
      summarise(concentration = mean(concentration))
    PerMHG.values <- Hg.data.date %>%
      filter(constituent == c("perPMHG"))
    PerMHG.values.no.dup <- PerMHG.values %>%
      group_by(depth) %>%
      summarise(concentration = mean(concentration))
    } else if (phase == "particulate-SPM") {
      Hg.adjustment.factor <- 0.005
      Hg.label <- "Particulate Hg (ng/g)"
      MHG.values <- Hg.data.date %>%
        filter(constituent == "PMHG_NG.G") %>%
        arrange(depth) %>%
        select(depth, concentration) %>%
        mutate(concentration = concentration)
      MHG.values.no.dup <- MHG.values %>%
        group_by(depth) %>%
        summarise(concentration = mean(concentration))
      iHG.values <- Hg.data.date %>%
        filter(constituent == "PiHG_NG.G") %>%
        arrange(depth) %>%
        select(depth, concentration) %>%
        mutate(concentration = concentration)
      iHG.values.no.dup <- iHG.values %>%
        group_by(depth) %>%
        summarise(concentration = mean(concentration))
      PerMHG.values <- Hg.data.date %>%
        filter(constituent == c("perPMHG"))
      PerMHG.values.no.dup <- PerMHG.values %>%
        group_by(depth) %>%
        summarise(concentration = mean(concentration))
    
  }
  
  # Plot MeHg values
  plot(x = MHG.values$concentration*Hg.adjustment.factor,
       y = MHG.values$depth,
       ylim = c(25, -2),
       xlim = c(0, 1.5),
       pch = 17,
       cex = 1.5,
       col = cb.translator["vermillion"],
       xaxt = "n",
       xlab = "",
       yaxt = "n",
       ylab = "")
  lines(MHG.values.no.dup$concentration*Hg.adjustment.factor,
        MHG.values.no.dup$depth,
        col = cb.translator["vermillion"])
  
  # Plot iHg values
  points(x = iHG.values$concentration*Hg.adjustment.factor,
         y = iHG.values$depth,
         pch = 18,
         cex = 1.5,
         col = cb.translator["bluishgreen"])
  lines(x = iHG.values.no.dup$concentration*Hg.adjustment.factor,
        y = iHG.values.no.dup$depth,
        col = cb.translator["bluishgreen"])
  
  # Plot percent MeHg values
  points(x = PerMHG.values$concentration*0.015,
         y = PerMHG.values$depth,
         pch = 16,
         cex = 1.5,
         col = cb.translator["black"])
  lines(x = PerMHG.values.no.dup$concentration*0.015,
        y = PerMHG.values.no.dup$depth,
        col = cb.translator["black"])
  
  # Add y-axis with depth measurements
  axis(2,
       at = seq(0, 25, by = 5),
       labels = seq(0, 25, by = 5))
  # Add axis for Hg
  axis(1,
       at = seq(0, 1.5, by = 0.3),
       labels = seq(0, 1.5, by = 0.3)/Hg.adjustment.factor)
  axis(1,
       line = 3,
       at = seq(0, 1.5, by = 0.375),
       labels = seq(0, 1.5, by = 0.375)/1.5)
  
  # Add labels
  # Add label for percent MeHg
  title(xlab = "Fraction MeHg",
        line = 4.5)
  # Add label for depth
  title(ylab = "Depth (m)",
        line = 1.5)
  # Add label for sulfide
  title(xlab = Hg.label,
        line = 1.5)
  
  
  # Add legend
  legend("topright",
         legend = c("iHg",
                    "MeHg",
                    "% MeHg"),
         col = c(cb.translator["bluishgreen"],
                 cb.translator["vermillion"],
                 cb.translator["black"]),
         pch = c(18, 17, 16),
         pt.cex = c(1.75, 1.5, 1.75),
         bg = "white")
}




#### Define function: To combine plots for a specific date ####
generate_geochem_plots <- function(date.of.interest = date.to.use,
                                   sampling.depths.of.interest = NULL) {
  

  #### Plot sonde casts ####
  plot.exo.data(date.of.interest,
                main.title = "",
                legend.location = "bottomright")
  text(x = 4,
       y = 1.5,
       cex = 1.2,
       labels = paste(month(date.of.interest, label = TRUE), " ", day(date.of.interest), ", ", year(date.of.interest),
                      sep = ""))
  if (!is.null(sampling.depths.of.interest)) {
    for (depth.of.interest in sampling.depths.of.interest) {
      abline(h = -depth.of.interest)
    }
  }
  
  
  #### Plot metals data ####
  color.vector <- c(cb.translator["orange"],
                    cb.translator["yellow"])
  names(color.vector) <- c("Mn_ppm", "Fe_ppm")
  point.vector <- c(1, 16)
  names(point.vector) <- c("particulate", "dissolved")
  
  plot.metals.data(icp.data.to.use = metals.data,
                   constituent.of.interest = point.vector,
                   color.vector.to.use = color.vector,
                   point.vector.to.use = point.vector,
                   xlim.to.use = c(0, 0.4),
                   date.of.interest = date.of.interest,
                   say.anything.for.two.legends = "okay")
  
  
  #### Plot sulfide data
  plot.sulfide.sulfate.data(sulfide.data.to.use = sulfide.data %>%
                              filter(startDate == date.of.interest),
                            sulfate.data.to.use = sulfate.data %>%
                              filter(startDate == date.of.interest),
                            legend.location = "topright")
  
  
  plot.Hg.profile(date.to.use = date.of.interest)
  # plot.Hg.profile(date.to.use = date.of.interest,
  #                 phase = "particulate")
  plot.Hg.profile(date.to.use = date.of.interest,
                  phase = "particulate-SPM")
  
  
  
  
  
}


pdf("results/BGC_profiles/BGC_profiles.pdf",
    width = 8,
    height = 15)
par(mfrow = c(4, 5),
    mar = c(9, 3, 2.5, 1),
    mgp = c(1.5, 0.4, 0),
    tck = -0.008)
generate_geochem_plots("2020-09-02")
generate_geochem_plots("2020-10-10")
generate_geochem_plots("2021-09-10")
generate_geochem_plots("2021-10-14")
dev.off()
