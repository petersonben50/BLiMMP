

#### Get set up ####
rm(list = ls())
setwd("~/Documents/research/BLiMMP")
library(tidyverse)
library(gtools)
library(LakeMetabolizer)
library(rLakeAnalyzer)
library(stringr)
library(RcppRoll)
library(lubridate)
library(viridis)
cb.translator <- readRDS("~/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")


#### Download data ####
# Go for the hourly data.
# Find it here: https://lter.limnology.wisc.edu/node/54845/data_form
# Download all columns, from 2020-04-01 to 2020-12-31.
# Save it here: ~/Documents/research/BLiMMP/dataRaw/buoy/sensor_mendota_lake_watertemp_hourly.csv

#### Set variables of interest for testing ####

data.file.of.interest <- "dataRaw/buoy/sensor_mendota_lake_watertemp_hourly.csv"
year.of.interest <- 2020
starting.date <- "2020-04-30"
ending.date <- "2020-11-05"

data.to.use <- read.csv(file = data.file.of.interest,
                        stringsAsFactors = FALSE)
#### Plot out heat maps ####

pdf("results/manuscript_figs/2020_heatmap_thermocline.pdf",
    width = 7.5,
    height = 5)
temp.profile.thermocline(data.of.interest = data.to.use,
                         year.of.interest = 2020,
                         starting.date = "2020-04-30",
                         ending.date = "2020-11-05")

mtext("˚C",
      side = 4,
      line = -1,
      cex = 1.2)
dev.off()
