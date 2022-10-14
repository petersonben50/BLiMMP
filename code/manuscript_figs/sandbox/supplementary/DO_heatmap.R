#### code/waterChem/DO_heatmap.R ####
# Benjamin D. Peterson


#### Set it up ####
rm(list = ls())
setwd("~/Documents/research/BLiMMP")
library(akima)
# library(colorspace)
# library(ggplot2)
library(fields)
library(ggimage)
library(lubridate)
library(tidyverse)
library(viridis)
cb.translator <- readRDS("~/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")


#### Read data ####
sonde.data <- readRDS("dataEdited/MeMO_exo/MeMO_exo_data.rds")


#### Dates ####
date.vector <- c("May", "Jun", "Jul", "Aug",
                 "Sep", "Oct", "Nov")
date.naming.vector <- yday(mdy(paste(date.vector,
                                     "-01-2020",
                                     sep = "")))
names(date.naming.vector) <- date.vector


#### Make common limits for both years ####
# Start May 1st
start.date <- yday("2020-04-30")
end.date <- yday("2020-11-05")



#### Define function ####
# year.of.interest <- 2020
DO.heatmap <- function(year.of.interest) {
  sonde.data.year <- sonde.data %>%
    filter(year(sampleDate) == year.of.interest) %>%
    mutate(number.day.of.year = yday(sampleDate)) %>%
    filter(!is.na(do_raw))
  
  # Extend down to 24 m, use values from 20 m. 
  heatmap.data <- sonde.data.year %>%
    select(number.day.of.year, depth, do_raw) %>%
    spread(key = depth,
           value = do_raw) %>%
    mutate(`24` = `20`) %>%
    gather(key = depth,
           value = do_raw,
           -1) %>%
    arrange(number.day.of.year)
    
  
  
  heatmap.data.interp <- interp(x = heatmap.data$number.day.of.year,
                                y = -as.numeric(heatmap.data$depth),
                                z = heatmap.data$do_raw)
  
  par(mar = c(3, 3.5, 2, 1),
      mgp=c(1.5,0.4,0),
      tck=-0.01)
  image.plot(heatmap.data.interp,
             axes = F,
             col = viridis(20),
             xlim = c(start.date,
                      end.date),
             zlim = c(0,15))

  axis(side = 2,
       at = seq(-25, 0,
                by = 5),
       labels = seq(25, 0,
                    by = -5))
  
  axis(side = 1,
       at = date.naming.vector,
       labels = names(date.naming.vector),
       las = 1)
  
  title(main = paste("DO profile for ",
                     year.of.interest,
                     sep = ""),
        ylab = "Depth (m)")
  
  mtext("mg/L",
        side=4,
        line=3.3,
        cex=1.2)
  
  points(x = c(rep(yday("2020-09-02,2200"), 3),
               rep(yday("2020-10-10 2200"), 2)),
         y = -c(20.7, 15.5, 11.0,
               20.9, 15.7),
         cex = 1.4,
         pch = 18,
         col = cb.translator['vermillion'])
  
  
}


#### Save out image of heatmap ####
pdf("results/manuscript_figs/2020_DO_heatmap.pdf",
    width = 7.5,
    height = 5)
par(mfrow = c(1, 1))
DO.heatmap(2020)
dev.off()
