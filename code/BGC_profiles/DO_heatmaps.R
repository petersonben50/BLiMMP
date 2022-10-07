#### code/BGC_profiles/DO_heatmaps.R ####
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("~/Documents/research/BLiMMP")
library(fields)
# library(gplots)
library(lubridate)
library(RColorBrewer)
library(readxl)
library(tidyverse)
library(viridis)
cb.translator <- readRDS("references/colorblind_friendly_colors.rds")


#### Read in data ####
DO.profile <- readRDS("dataEdited/DO_profiles/DO_profiles_data.rds")


#### Make date vector ####

date.vector <- c("May 1st", "Jun 1st", "July 1st", "August 1st")
date.naming.vector <- yday(mdy(paste(date.vector,
                                     ", 2017",
                                     sep = "")))
names(date.naming.vector) <- date.vector




#### Function to plot out heatmap ####
DO.heatmap <- function(year.of.interest,
                       date.limits = c("2020-05-01", "2020-08-01")) {
  # Generate a matrix of the DO values
  clean.DO.profile <- DO.profile %>%
    filter(year(sampleDate) == year.of.interest) 
  
  clean.DO.profile.matrix <- clean.DO.profile%>%
    select(-sampleDate) %>%
    as.matrix() %>%
    image.smooth()
  
  # Add daynum and depth as x and y, respectively
  clean.DO.profile.matrix$x <- yday(clean.DO.profile$sampleDate)
  clean.DO.profile.matrix$y <- names(clean.DO.profile)[-1] %>% as.numeric()
  
  depths.to.plot.on.axes <- seq(20, 0, by = -5)
  dates.to.plot.on.axes <- seq(from = as.Date(date.limits[1]),
                               to = as.Date(date.limits[2]),
                               by = 'months') 
  
  
  #### Plotting function ####
  par(mar = c(3, 3.5, 2, 1),
      mgp=c(1.5,0.4,0),
      tck=-0.01)
  filled.contour(x = clean.DO.profile.matrix$x,
                 y = clean.DO.profile.matrix$y,
                 z = clean.DO.profile.matrix$z,
                 xlim = date.limits %>% yday(),
                 ylim = c(20, 0),
                 zlim = c(0, 14),
                 nlevels = 40,
                 color.palette = colorRampPalette(c(cb.translator["black"], cb.translator["blue"], cb.translator["bluishgreen"], cb.translator["orange"]),
                                                  bias = 1, space = "rgb"),
                 plot.axes = { 
                   axis(1,
                        at = dates.to.plot.on.axes %>% as.Date() %>% yday(),
                        labels = paste(month(dates.to.plot.on.axes, label = TRUE, abbr = FALSE),
                                       " 1st"));
                   axis(2,
                        at = depths.to.plot.on.axes);
                   
                 },
                 key.title = {par(cex.main=1);title(main = "DO (mg/L)")},
                 key.axes = axis(4, seq(0, 15, by = 3),asp=1)
  )
}



#### Generate the heatmap ####
pdf("results/BGC_profiles/DO_heatmaps.pdf",
    width = 7.5,
    height = 4)
DO.heatmap(year.of.interest = "2020",
           date.limits = c("2020-05-01", "2020-08-01"))
DO.heatmap(year.of.interest = "2021",
           date.limits = c("2021-05-01", "2021-08-01"))
dev.off()
