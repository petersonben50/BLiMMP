#### code/BGC_profiles/TEAP_profiles_2021.R ####
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("~/Documents/research/BLiMMP/")
library(readxl)
library(tidyverse)
source("code/BLiMMP_functions.R")
detection.vector <- c(4, 16)


#### Read in data ####
sulfide.data <- read.csv("dataEdited/waterChemistry/sulfide/WC_data.csv")
sulfate.data <- read.csv("dataEdited/waterChemistry/sulfate/WC_data.csv")

metals.data.2021 <- readRDS("dataEdited/waterChemistry/ICP/2020_2021_Fe_Mn_data.rds")


#### Function to combine plots ####
generate_TEAP_plots <- function(date.of.interest = date.to.use,
                                sampling.depths.of.interest = NULL) {
  
  pdf(paste("results/BGC_profiles/TEAP_profile_", date.of.interest, ".pdf",
            sep = ""),
      width = 7.5,
      height = 5)
  par(mfrow = c(1, 3))
  
  # Sonde casts
  plot.exo.data(date.of.interest,
                main.title = "",
                legend.location = "bottomright")
  text(x = 2.5,
       y = -0.1,
       cex = 1.5,
       labels = paste(month(date.of.interest, label = TRUE), " ", day(date.of.interest), ", ", year(date.of.interest),
                      sep = ""))
  if (!is.null(sampling.depths.of.interest)) {
    for (depth.of.interest in sampling.depths.of.interest) {
      abline(h = -depth.of.interest)
    }
  }
  
  
  # Metals data
  color.vector <- c(cb.translator["orange"],
                    cb.translator["yellow"])
  names(color.vector) <- c("Mn_ppb", "Fe_ppb")
  point.vector <- c(1, 16)
  names(point.vector) <- c("particulate", "dissolved")
  plot.metals.data(icp.data.to.use = metals.data.2021,
                   constituent.of.interest = point.vector,
                   color.vector.to.use = color.vector,
                   point.vector.to.use = point.vector,
                   xlim.to.use = c(0, 400),
                   date.of.interest = date.of.interest,
                   legend.location = "bottomright")
  if (!is.null(sampling.depths.of.interest)) {
    for (depth.of.interest in sampling.depths.of.interest) {
      abline(h = depth.of.interest)
    }
  }
  
  # Sulfide data
  plot.sulfide.sulfate.data(sulfide.data.to.use = sulfide.data %>%
                              filter(startDate == date.of.interest),
                            sulfate.data.to.use = sulfate.data %>%
                              filter(startDate == date.of.interest),
                            legend.location = "bottomright")
  if (!is.null(sampling.depths.of.interest)) {
    for (depth.of.interest in sampling.depths.of.interest) {
      abline(h = depth.of.interest)
    }
  }
  
  
  dev.off()
  
}



#### Plots: run function ####
generate_TEAP_plots("2020-09-02",
                    sampling.depths.of.interest = c(11, 15.5, 20.7))
generate_TEAP_plots("2020-10-10",
                    sampling.depths.of.interest = c(15.7, 20.9))

generate_TEAP_plots("2021-09-10",
                    sampling.depths.of.interest = c(10.8, 11.9, 19.9))
generate_TEAP_plots("2021-10-14",
                    sampling.depths.of.interest = c(14.2, 15.2, 19.9))