#### code/incubations/incubation_analysis_2019_204.R ####
# Written for BLiMMP project
# Benjamin D. Peterson



#### Prep workspace ####

rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP")
library(dplyr)
library(ggplot2)
library(lubridate)
library(patchwork)
library(readxl)
library(tidyr)





#### Read in MeHg data ####

MeHg.data <- read.csv("dataEdited/incubations/MeHg/incubations2019_MeHg_ambient.csv",
                      stringsAsFactors = FALSE)





#### Read in colors ####

cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")




#### Function to generate empty plot ####

empty.plot <- function(x) {
  plot(x = 0,
       cex = 0,
       xaxt = "n",
       yaxt = "n",
       xlab = "",
       bty = "n")
}




#### Function to plot change in ambient MeHg for each depth on a given trip ####

plot.MeHg.time.course.depth <- function(selected.depth,
                                           trip_data,
                                           color.vector.input = NULL) {
  
  par(mar=c(3,3,2.5,1), mgp=c(1.5,0.4,0), tck=-0.008)
  
  #y_max <- ceiling(max(trip_data$amb_MeHg_ng.L)*100)/100 + 0.01
  y_max <- 1
  y_min <- floor(min(trip_data$amb_MeHg_ng.L)*100)/100
  
  # Isolate by depth
  trip_data_depth <- trip_data %>%
    filter(depth == selected.depth)
  
  #### Make a color vector ####
  if (is.null(color.vector.input)) {
    color.vector <- c(cb.translator["black"],
                      cb.translator["bluishgreen"],
                      cb.translator["vermillion"])
    names(color.vector) <- c("unfiltered-molybdate",
                             "unfiltered-unamended",
                             "filtered-unamended")
  } else {
    color.vector <- color.vector.input
  }
  
  ####  Make a time vector #### 
  time.vector <- c(0, 1)
  names(time.vector) <- c("t0", "t1")
  

  plot(x = time.vector[trip_data_depth$t],
       y = trip_data_depth$amb_MeHg_ng.L,
       ylim = c(y_min, y_max),
       xlim = c(-0.5, 1.5),
       pch = 18,
       col = color.vector[trip_data_depth$treatment],
       cex = 1.25,
       xaxt = "n",
       xlab = "",
       ylab = "",
       main = paste(unique(trip_data_depth$startDate),
                    ": ",
                    selected.depth,
                    "m",
                    sep = ""))
  axis(1,
       at = c(0, 1),
       labels = c("t0", "t1"))
  
  #### Add segments to link samples within an incubation ####
  
  trip_data_depth_segments <- trip_data_depth %>%
    select(-c(depth)) %>%
    spread(key = t,
           value = amb_MeHg_ng.L)
  
  for (row.num in 1:nrow(trip_data_depth_segments)) {
    segments(0, trip_data_depth_segments$t0[row.num],
             1, trip_data_depth_segments$t1[row.num],
             col = color.vector[trip_data_depth_segments$treatment[row.num]])
  }
  
  title(ylab = "Ambient MeHg (ng/L)",
        line = 1.75)
  
  legend("topright",
         legend =  names(color.vector)[length(color.vector):1],
         text.col = color.vector[length(color.vector):1],
         bty = "n")
}






#### Function to run plotting function over each depth ####

plot.MeHg.changes <- function(trip.of.interest,
                              color.input = NULL) {
  trip_data_df <- MeHg.data %>%
    filter(tripID == trip.of.interest) %>%
    select(incubationID, depth, startDate, t, treatment, amb_MeHg_ng.L)
  sapply(unique(trip_data_df$depth) %>% sort(),
         function(depth) {
           plot.MeHg.time.course.depth(selected.depth = depth,
                                          trip_data = trip_data_df,
                                          color.vector.input = color.input)
         }
  )
}







#### Plot out change in MeHg depths ####

# For 006
color.vector.006 <- c(cb.translator["black"],
                      cb.translator["bluishgreen"],
                      cb.translator["vermillion"],
                      cb.translator["orange"],
                      cb.translator["blue"],
                      cb.translator["reddishpurple"])
names(color.vector.006) <- c("unfiltered-molybdate",
                             "unfiltered-unamended",
                             "filtered-unamended",
                             "unfiltered-starch",
                             "unfiltered-starch-molybdate",
                             "unfiltered-algal")

png("results/incubations/ambient_Hg/ambient_MeHg_incubations.png",
    units = "in",
    res = 240,
    width = 10,
    height = 12)

par(mfrow = c(3, 3))
plot.MeHg.changes(trip.of.interest = "BLiMMP_trip_003")
empty.plot()
plot.MeHg.changes(trip.of.interest = "BLiMMP_trip_005")
plot.MeHg.changes(trip.of.interest = "BLiMMP_trip_006",
                  color.input = color.vector.006)
empty.plot()

dev.off()


rm(amb_MeHg_ng.L,
   MeHg.data,
   plot.MeHg.changes)
