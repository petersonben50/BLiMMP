#### code/BGC_profiles/TEAP_profiles_2021.R ####
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("~/Documents/research/BLiMMP/")
library(readxl)
library(tidyverse)
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")
source("code/BLiMMP_functions.R")
detection.vector <- c(4, 16)


#### Plotting functions ####
plot.sulfide.sulfate.data <- function(sulfide.data.to.use,
                                      sulfate.data.to.use,
                                      xlim.to.use = c(0, 250)) {
  # Sulfide data
  plot(x = sulfide.data.to.use$S_conc_uM,
       y = sulfide.data.to.use$depth,
       pch = 16,
       ylim = c(25, 0),
       xlim = xlim.to.use,
       col = cb.translator["blue"],
       xlab = "Sulfide/Sulfate (uM)",
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
  
  legend(x = sum(xlim.to.use)/5,
         y = 0,
         legend = c("Sulfate (µM)", "Sulfide (µM)"),
         pch = c(1, 16),
         lty = c(1, 2),
         col= cb.translator["blue"])
}



#### Read in data ####
sulfide.data <- read.csv("dataEdited/waterChemistry/sulfide/WC_data.csv")
sulfate.data <- read.csv("dataEdited/waterChemistry/sulfate/WC_data.csv")



#### Plots for September 2020 ####
date.to.use <- "2020-09-02"
pdf(paste("results/BGC_profiles/TEAP_profile_", date.to.use, ".pdf",
          sep = ""),
    width = 7,
    height = 7)
par(mfrow = c(1, 2))

# Sonde casts
plot.exo.data(date.to.use,
              main.title = "",
              legend.location = "topleft")
text(x = 6,
     y = -24,
     cex = 1,
     labels = "September 2nd, 2020")

# Sulfide and sulfate data
plot.sulfide.sulfate.data(sulfide.data.to.use = sulfide.data %>%
                            filter(startDate == date.to.use),
                          sulfate.data.to.use = sulfate.data %>%
                            filter(startDate == date.to.use))
dev.off()



#### Plots for October 2020 ####
date.to.use <- "2020-10-10"
pdf(paste("results/BGC_profiles/TEAP_profile_", date.to.use, ".pdf",
          sep = ""),
    width = 7,
    height = 7)
par(mfrow = c(1, 2))

# Sonde casts
plot.exo.data(date.to.use,
              main.title = "",
              legend.location = "topleft")
text(x = 6,
     y = -24,
     cex = 1,
     labels = "October 10th, 2020")

# Sulfide and sulfate data
plot.sulfide.sulfate.data(sulfide.data.to.use = sulfide.data %>%
                            filter(startDate == date.to.use),
                          sulfate.data.to.use = sulfate.data %>%
                            filter(startDate == date.to.use))
dev.off()




#### Plots for September 2021 ####
date.to.use <- "2021-09-10"
pdf(paste("results/BGC_profiles/TEAP_profile_", date.to.use, ".pdf",
          sep = ""),
    width = 7,
    height = 7)
par(mfrow = c(1, 2))

# Sonde casts
plot.exo.data(date.to.use,
              main.title = "",
              legend.location = "topleft")
text(x = 6,
     y = -24,
     cex = 1,
     labels = "September 10th, 2021")

# Sulfide data
plot.sulfide.sulfate.data(sulfide.data.to.use = sulfide.data %>%
                    filter(startDate == date.to.use),
                  sulfate.data.to.use = sulfate.data %>%
                    filter(startDate == date.to.use))
dev.off()



#### Plots for October 2021 ####
date.to.use <- "2021-10-14"
pdf(paste("results/BGC_profiles/TEAP_profile_", date.to.use, ".pdf",
          sep = ""),
    width = 7,
    height = 7)
par(mfrow = c(1, 2))

# Sonde casts
plot.exo.data(date.to.use,
              main.title = "",
              legend.location = "topleft")
text(x = 6,
     y = -24,
     cex = 1,
     labels = "October 14th, 2021")

# Sulfide data
plot.sulfide.sulfate.data(sulfide.data %>%
                            filter(startDate == date.to.use),
                          sulfate.data %>%
                            filter(startDate == date.to.use))
dev.off()
