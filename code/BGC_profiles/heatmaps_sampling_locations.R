#### code/BGC_profiles/heatmaps_sampling_locations.R ####
# Benjamin D. Peterson

# This script will generate heatmaps of the temperature
# and DO of the lake over the course of the ice-free
# season and include markers for our sampling locations.


#### Get set up ####
rm(list = ls())
setwd("~/Documents/research/BLiMMP")
library(tidyverse)
# library(gtools)
# library(LakeMetabolizer)
library(rLakeAnalyzer)
library(stringr)
# library(RcppRoll)
library(readxl)
library(lubridate)
library(viridis)
source("code/BLiMMP_functions.R")


#### Download data ####
# Go for the hourly data.
# Find it here: https://lter.limnology.wisc.edu/node/56136/data_form
# Download all columns, from 2020-04-01 to 2021-12-31.
# Save it here: ~/Documents/research/BLiMMP/dataRaw/buoy/sensor_mendota_lake_watertemp_hourly.csv


#### Set variables of interest for testing ####
data.to.use <- read.csv(file = "dataRaw/buoy/sensor_mendota_lake_watertemp_hourly.csv",
                        stringsAsFactors = FALSE)
sampling.location.info <- read_xlsx("dataEdited/incubations/incubation_sites_notes.xlsx")




#### Define function: Plot temp profile over year ####
temp.profile.thermocline <- function(data.of.interest,
                                     year.of.interest,
                                     starting.date,
                                     ending.date,
                                     temp.scale = c(4, 30),
                                     sampling.date.location.df.to.use = NULL) {
  
  if (!exists("data.of.interest")) {
    print("Gotta supply a data file.")
    break()
  }
  
  if (!exists("year.of.interest")) {
    print("Gotta set year of interest.")
    break()
  }
  
  # Read in the data and clean the timestamps
  LTERdf <- data.of.interest %>%
    filter(year4 == year.of.interest) %>%
    mutate(hour = str_pad(hour, 4, pad = "0")) %>%
    mutate(date_time = ymd_hm(paste(sampledate, hour))) %>%
    select(date_time, depth, wtemp) %>%
    spread(key = depth,
           value = wtemp) %>%
    mutate(`24` = `20`) %>%
    gather(key = depth,
           value = wtemp,
           -1) %>%
    mutate(depth = as.numeric(depth)) %>%
    arrange(date_time)
  
  if (exists("starting.date")) {
    print(paste("Cutting off dates before ",
                starting.date,
                ".",
                sep = ""
    ))
    LTERdf <- LTERdf %>%
      filter(as.Date(date_time) >= as.Date(starting.date))
  }
  if (exists("ending.date")) {
    print(paste("Cutting off dates after ",
                ending.date,
                ".",
                sep = ""
    ))
    LTERdf <- LTERdf %>%
      filter(as.Date(date_time) <= ending.date)
  }
  
  
  
  # Set up depth and date vectors
  depths <- unique(LTERdf$depth)
  dates <- unique(LTERdf$date_time)
  depths.to.plot <- seq(25, 0, by = -5)
  dates.to.plot <- seq(ceiling_date(min(dates), "months"),
                       floor_date(max(dates), "months"),
                       by='months')
  
  
  LTERdf <- LTERdf %>%
    spread(key = depth,
           value = wtemp,
           fill = NA)
  
  
  # Find dates with NAs and exclude them
  # LTERdf.missing.vector <- LTERdf %>%
  #   select(-date_time) %>%
  #   rowSums() %>%
  #   is.na()
  # if (length(which(LTERdf.missing.vector)) > 0) {
  #   print(paste("NAs present on ",
  #               LTERdf[which(LTERdf.missing.vector), "date_time"],
  #               sep = ""))
  #   LTERdf <- LTERdf[!LTERdf.missing.vector, ]
  # }
  # 
  
  # Set up matrix
  LTERmt <- LTERdf %>%
    select(-date_time) %>%
    as.matrix()
  # LTERmt_top<-LTERmt[,1:15]
  
  
  
  #### Calculate thermocline depth ####
  names(LTERdf)[-1] <- paste("var_", names(LTERdf)[-1], sep = "")
  thermocline.depth.ts <- LTERdf %>%
    rename(datetime = date_time) %>%
    ts.thermo.depth() %>%
    group_by(as.Date(datetime)) %>%
    summarize(thermo.depth = median(thermo.depth,
                                    na.rm = TRUE)) %>%
    rename(datetime = `as.Date(datetime)`)
  
  
  if (is.null(sampling.date.location.df.to.use)) {
    sampling.date.location.df.to.use = as.data.frame(cbind(NA, NA))
    names(sampling.date.location.df.to.use) <- c("date", "depth")
  } else {
    sampling.date.location.df.to.use <- sampling.date.location.df.to.use %>%
      filter(as.Date(date) >= as.Date(starting.date)) %>%
      filter(as.Date(date) <= as.Date(ending.date))
    
  }
  
  # Generate plot
  par(mar = c(3, 3.5, 2, 1),
      mgp=c(1.5,0.4,0),
      tck=-0.01)
  filled.contour(x = dates,
                 y = depths,
                 z = LTERmt,
                 ylim = c(max(depths), 0),
                 zlim = temp.scale,
                 nlevels = 40,
                 color.palette = colorRampPalette(c("violet", "blue", "cyan", "green3", "yellow", "orange", "red"),
                                                  bias = 1, space = "rgb"),
                 # plot.title = { title(main = paste("Temperature profile and thermocline depth for ",
                 #                                   year.of.interest,
                 #                                   sep = "")) },
                 plot.axes = { axis(1,
                                    at = dates.to.plot,
                                    labels = month(dates.to.plot, label = TRUE));
                   axis(2,
                        at = depths.to.plot);
                   lines(as.POSIXct(thermocline.depth.ts$datetime),
                         thermocline.depth.ts$thermo.depth,
                         lwd = 1,
                         col = 'black')
                   points(x = as.POSIXct(sampling.date.location.df.to.use$date),
                          y = sampling.date.location.df.to.use$depth,
                          cex = 1.4,
                          pch = 18,
                          col = cb.translator['vermillion']) }
  )
  
}


#### Plot out heat maps ####
temp.profile.thermocline(data.of.interest = data.to.use,
                         year.of.interest = 2020,
                         starting.date = "2020-04-30",
                         ending.date = "2020-11-05",
                         sampling.date.location.df.to.use = sampling.location.info)


temp.profile.thermocline(data.of.interest = data.to.use,
                         year.of.interest = 2021,
                         starting.date = "2021-04-30",
                         ending.date = "2021-11-05",
                         sampling.date.location.df.to.use = sampling.location.info)
