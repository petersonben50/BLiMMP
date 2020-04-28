#### code/incubations/incubation_analysis_2019_THg.R ####
# Benjamin D. Peterson
# Written for BLiMMP project

#### Prep workspace ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP/")
library(dplyr)
library(lubridate)
library(readxl)
library(tidyr)



#### Read in MeHg data ####

THg.data <- read.csv("dataEdited/incubations/THg/incubations2019_THg.csv",
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

#### Function to plot change in THg for each depth on a given trip ####

plot.THg.time.course.depth <- function(selected.depth,
                                       trip_data,
                                       color.vector.input = NULL,
                                       isotope = 198) {
  
  if (isotope == 198) {
    trip_data <- trip_data %>%
      rename(THg_ng_L = THg_198)
  } else if (isotope == 204) {
    trip_data <- trip_data %>%
      rename(THg_ng_L = THg_204)
  } else {
    break
  }
  
  trip_data <- trip_data %>%
    select(incubationID, startDate, depth, treatment, t, THg_ng_L)
  
  par(mar=c(3,3,2.5,1), mgp=c(1.5,0.4,0), tck=-0.008)
  
  y_max <- ceiling(max(trip_data$THg_ng_L)*100)/100 + 0.01
  y_min <- floor(min(trip_data$THg_ng_L)*100)/100
  
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
  
  #### Make plot ####
  plot(x = time.vector[trip_data_depth$t],
       y = trip_data_depth$THg_ng_L,
       ylim = c(y_min, y_max),
       xlim = c(-0.5, 1.5),
       pch = 18,
       col = color.vector[trip_data_depth$treatment],
       cex = 1.25,
       xaxt = "n",
       xlab = "",
       ylab = "",
       main = paste(selected.depth,
                    "m",
                    sep = ""))
  axis(1,
       at = c(0, 1),
       labels = c("t0", "t1"))
  
  #### Add segments to link samples within an incubation ####
  
  trip_data_depth_segments <- trip_data_depth %>%
    select(-c(depth)) %>%
    spread(key = t,
           value = THg_ng_L)
  
  for (row.num in 1:nrow(trip_data_depth_segments)) {
    segments(0, trip_data_depth_segments$t0[row.num],
             1, trip_data_depth_segments$t1[row.num],
             col = color.vector[trip_data_depth_segments$treatment[row.num]])
  }
  
  title(ylab = paste("Excess T",
                     isotope,
                     "Hg (ng/L)",
                     sep = ""),
        line = 1.75)
  
  legend("topleft",
         legend =  names(color.vector)[length(color.vector):1],
         text.col = color.vector[length(color.vector):1],
         bty = "n")
}







#### Plot 198 data for each depth within a trip ####

png("results/incubations/THg_198_BLiMMP_incubations.png",
    units = "in",
    res = 240,
    height = 10,
    width = 8)
par(mfrow = c(3, 3))



color.vector.003 <- c(cb.translator["bluishgreen"],
                      cb.translator["vermillion"])
names(color.vector.003) <- c("unfiltered-unamended",
                             "filtered-unamended")
trip.of.interest <- "BLiMMP_trip_003"
trip_data_df <- THg.data %>%
  filter(tripID == trip.of.interest)
sapply(X = sort(unique(trip_data_df$depth)),
       function(x) {
         plot.THg.time.course.depth(x,
                                    trip_data = trip_data_df,
                                    color.vector.input = color.vector.003)
       }
)
empty.plot()



# Trip 005

color.vector.005 <- c(cb.translator["vermillion"],
                      cb.translator["bluishgreen"],
                      cb.translator["black"])
names(color.vector.005) <- c("filtered-unamended",
                             "unfiltered-unamended",
                             "unfiltered-molybdate")

trip.of.interest <- "BLiMMP_trip_005"
trip_data_df <- THg.data %>%
  filter(tripID == trip.of.interest)
sapply(X = sort(unique(trip_data_df$depth)),
       function(x) {
         plot.THg.time.course.depth(x,
                                    trip_data = trip_data_df,
                                    color.vector.input = color.vector.005)
       }
)


color.vector.006 <- c(cb.translator["vermillion"],
                      cb.translator["bluishgreen"],
                      cb.translator["black"],
                      cb.translator["orange"],
                      cb.translator["blue"],
                      cb.translator["reddishpurple"])
names(color.vector.006) <- c("filtered-unamended",
                             "unfiltered-unamended",
                             "unfiltered-molybdate",
                             "unfiltered-starch",
                             "unfiltered-starch-molybdate",
                             "unfiltered-algal")

trip.of.interest <- "BLiMMP_trip_006"
trip_data_df <- THg.data %>%
  filter(tripID == trip.of.interest)

sapply(X = sort(unique(trip_data_df$depth)),
       function(x) {
         plot.THg.time.course.depth(x,
                                    trip_data = trip_data_df,
                                    color.vector.input = color.vector.006)
       }
)

dev.off()





















#### Plot 204 data for each depth within a trip ####

png("results/incubations/THg_204_BLiMMP_incubations.png",
    units = "in",
    res = 240,
    height = 10,
    width = 8)
par(mfrow = c(3, 3))



color.vector.003 <- c(cb.translator["bluishgreen"],
                      cb.translator["vermillion"])
names(color.vector.003) <- c("unfiltered-unamended",
                             "filtered-unamended")
trip.of.interest <- "BLiMMP_trip_003"
trip_data_df <- THg.data %>%
  filter(tripID == trip.of.interest)
sapply(X = sort(unique(trip_data_df$depth)),
       function(x) {
         plot.THg.time.course.depth(x,
                                    trip_data = trip_data_df,
                                    color.vector.input = color.vector.003,
                                    isotope = 204)
       }
)
empty.plot()



# Trip 005

color.vector.005 <- c(cb.translator["vermillion"],
                      cb.translator["bluishgreen"],
                      cb.translator["black"])
names(color.vector.005) <- c("filtered-unamended",
                             "unfiltered-unamended",
                             "unfiltered-molybdate")

trip.of.interest <- "BLiMMP_trip_005"
trip_data_df <- THg.data %>%
  filter(tripID == trip.of.interest)
sapply(X = sort(unique(trip_data_df$depth)),
       function(x) {
         plot.THg.time.course.depth(x,
                                    trip_data = trip_data_df,
                                    color.vector.input = color.vector.005,
                                    isotope = 204)
       }
)


color.vector.006 <- c(cb.translator["vermillion"],
                      cb.translator["bluishgreen"],
                      cb.translator["black"],
                      cb.translator["orange"],
                      cb.translator["blue"],
                      cb.translator["reddishpurple"])
names(color.vector.006) <- c("filtered-unamended",
                             "unfiltered-unamended",
                             "unfiltered-molybdate",
                             "unfiltered-starch",
                             "unfiltered-starch-molybdate",
                             "unfiltered-algal")

trip.of.interest <- "BLiMMP_trip_006"
trip_data_df <- THg.data %>%
  filter(tripID == trip.of.interest)

sapply(X = sort(unique(trip_data_df$depth)),
       function(x) {
         plot.THg.time.course.depth(x,
                                    trip_data = trip_data_df,
                                    color.vector.input = color.vector.006,
                                    isotope = 204)
       }
)

dev.off()


