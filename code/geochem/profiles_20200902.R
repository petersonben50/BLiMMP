#### code/incubations/profiles_20200902.R ####
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("~/Documents/research/BLiMMP/")
library(tidyverse)
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")
source("code/BLiMMP_functions.R")



#### Generate geochem profiles ####
pdf("results/geochem/profile_20200902.pdf",
    height = 4.5,
    width = 6)
par(mfrow = c(1,3))
plot.exo.data("2020-09-02",
              points.of.sampling = c(11, 15.5, 20.7))
plot.redox.profile(trip = "BLiMMP_trip_010",
                   sulfide.data.location = "dataEdited/waterChemistry/sulfide/WC_data.csv",
                   sulfate.data.location = "dataEdited/waterChemistry/sulfate/WC_data.csv",
                   ICP.data.name = "dataEdited/waterChemistry/ICP/2020_WC_data.rds")
plot.Hg.profile(tripDates = c("2020-09-02", "2020-09-03"),
                Hg.data.location = "dataEdited/Hg/Hg_data_2020.csv")
dev.off()
