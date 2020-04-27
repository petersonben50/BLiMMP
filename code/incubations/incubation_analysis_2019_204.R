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

MeHg.data <- read.csv("dataEdited/incubations/MeHg/incubations2019_Me204Hg.csv",
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




#### Function to plot change in Me204Hg for each depth on a given trip ####

plot.Me204Hg.time.course.depth <- function(selected.depth,
                                           trip_data,
                                           color.vector.input = NULL) {
  
  par(mar=c(3,3,2.5,1), mgp=c(1.5,0.4,0), tck=-0.008)
  
  y_max <- ceiling(max(trip_data$excess_MeHg_204_ng.L)*100)/100 + 0.01
  y_min <- floor(min(trip_data$excess_MeHg_204_ng.L)*100)/100
  
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
  
  #### Add column for above DDL ####
  trip_data_depth[, "DDL_shape"] <- rep(4, length(trip_data_depth$above_DDL))
  trip_data_depth$DDL_shape[trip_data_depth$above_DDL] <- 18
  
  
  plot(x = time.vector[trip_data_depth$t],
       y = trip_data_depth$excess_MeHg_204_ng.L,
       ylim = c(y_min, y_max),
       xlim = c(-0.5, 1.5),
       pch = trip_data_depth$DDL_shape,
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
    select(-c(above_DDL, DDL_shape, depth)) %>%
    spread(key = t,
           value = excess_MeHg_204_ng.L)
  
  for (row.num in 1:nrow(trip_data_depth_segments)) {
    segments(0, trip_data_depth_segments$t0[row.num],
             1, trip_data_depth_segments$t1[row.num],
             col = color.vector[trip_data_depth_segments$treatment[row.num]])
  }
  
  title(ylab = "Excess Me204Hg (ng/L)",
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
    select(incubationID, depth, t, treatment, excess_MeHg_204_ng.L, above_DDL)
  sapply(unique(trip_data_df$depth) %>% sort(),
         function(depth) {
           plot.Me204Hg.time.course.depth(selected.depth = depth,
                                          trip_data = trip_data_df,
                                          color.vector.input = color.input)
         }
  )
}







#### Plot out change in Me204Hg depths ####

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

png("results/incubations/demethylation/demeth_raw_values.png",
    units = "in",
    res = 240,
    width = 10,
    height = 12)

par(mfrow = c(3, 3))
plot.MeHg.changes("BLiMMP_trip_003")
empty.plot()
plot.MeHg.changes("BLiMMP_trip_005")
plot.MeHg.changes("BLiMMP_trip_006",
                  color.input = color.vector.006)
empty.plot()

dev.off()


rm(plot.Me204Hg.time.course.depth,
   MeHg.data,
   plot.MeHg.changes)
















#### Function to plot percent Me204Hg at each time point for each depth on a given trip ####

per.MeHg.data <- read.csv("dataEdited/incubations/MeHg/incubations2019_Me204Hg_percent.csv",
                          stringsAsFactors = FALSE)

plot.per.Me204Hg.time.course.depth <- function(selected.depth,
                                           trip_data,
                                           color.vector.input = NULL) {
  
  par(mar=c(3,3,2.5,1), mgp=c(1.5,0.4,0), tck=-0.008)
  
  y_max <- ceiling(max(trip_data$per_MeHg_204)*100)/100 + 0.01
  y_min <- floor(min(trip_data$per_MeHg_204)*100)/100
  
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
       y = trip_data_depth$per_MeHg_204,
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
    select(-depth) %>%
    spread(key = t,
           value = per_MeHg_204)
  
  for (row.num in 1:nrow(trip_data_depth_segments)) {
    segments(0, trip_data_depth_segments$t0[row.num],
             1, trip_data_depth_segments$t1[row.num],
             col = color.vector[trip_data_depth_segments$treatment[row.num]])
  }
  
  title(ylab = "Fraction MeHg for 204Hg",
        line = 1.75)
  
  legend("topright",
         legend =  names(color.vector)[length(color.vector):1],
         text.col = color.vector[length(color.vector):1],
         bty = "n")
}






#### Function to run plotting function over each depth ####

trip.of.interest <- "BLiMMP_trip_003"
plot.per.MeHg.changes <- function(trip.of.interest,
                              color.input = NULL) {
  trip_data_df <- per.MeHg.data %>%
    filter(tripID == trip.of.interest) %>%
    select(incubationID, depth, startDate, t, treatment, per_MeHg_204)
  sapply(unique(trip_data_df$depth) %>% sort(),
         function(depth) {
           plot.per.Me204Hg.time.course.depth(selected.depth = depth,
                                              trip_data = trip_data_df,
                                              color.vector.input = color.input)
         }
  )
}




png("results/incubations/demethylation/demeth_percent.png",
    units = "in",
    res = 240,
    width = 10,
    height = 12)

par(mfrow = c(3, 3))
plot.per.MeHg.changes("BLiMMP_trip_003")
empty.plot()
plot.per.MeHg.changes("BLiMMP_trip_005")
plot.per.MeHg.changes("BLiMMP_trip_006",
                      color.input = color.vector.006)
empty.plot()

dev.off()

rm(list = ls(pattern = "per"),
   color.vector.006)




#### Change in MeHg over incubation ####

Me204Hg.per.change <- read.csv("dataEdited/incubations/MeHg/incubations2019_Me204Hg_percent_change.csv",
                               stringsAsFactors = FALSE)


# Generate needed color vector
color.vector.003 <- c(cb.translator["vermillion"],
                      cb.translator["bluishgreen"])
names(color.vector.003) <- c("filtered\nunamended",
                         "unfiltered\nunamended")

color.vector.005 <- c(cb.translator["vermillion"],
                      cb.translator["bluishgreen"],
                      cb.translator["black"])
names(color.vector.005) <- c("filtered\nunamended",
                             "unfiltered\nunamended",
                             "unfiltered\nmolybdate")

color.vector.006 <- c(cb.translator["black"],
                      cb.translator["bluishgreen"],
                      cb.translator["vermillion"],
                      cb.translator["orange"],
                      cb.translator["blue"],
                      cb.translator["reddishpurple"])
names(color.vector.006) <- c("unfiltered\nmolybdate",
                             "unfiltered\nunamended",
                             "filtered\nunamended",
                             "unfiltered\nstarch",
                             "unfiltered\nstarch\nmolybdate",
                             "unfiltered\nalgal")


plot.date.data <- function(date.of.interest,
                           color.vector) {
  
  per.change.MeHg.data.date <- Me204Hg.per.change %>%
    filter(startDate == date.of.interest)
  
  for (depth.number in 1:length(unique(per.change.MeHg.data.date$depth))) {
    
    depth.vector <- unique(per.change.MeHg.data.date$depth) %>% sort()
    depth.of.interest <- depth.vector[depth.number]
    
    depth.of.interest.m <- paste(depth.of.interest, "m", sep = "")
    
    plot.list[[depth.number]] <- per.change.MeHg.data.date %>%
      filter(depth == depth.of.interest) %>%
      mutate(treatment = gsub("-", "\n", treatment)) %>%
      mutate(treatment = factor(treatment, levels = names(color.vector))) %>%
      ggplot(aes(x = treatment,
                 y = change_in_per_MeHg_norm,
                 color = treatment)) +
      geom_jitter(position = position_jitter(0.1)) +
      scale_color_manual(values = color.vector) +
      theme_classic() +
      theme(legend.position = "none") +
      ylim(min(per.change.MeHg.data.date$change_in_per_MeHg_norm)-5,
           max(per.change.MeHg.data.date$change_in_per_MeHg_norm)+5) +
      xlab("") +
      ggtitle(paste(date.of.interest,
                    ": ",
                    depth.of.interest.m,
                    sep = ""))
    
  }
  
  if (length(plot.list) == 2) {
    plot.list[[3]] <- ggplot() + theme_void()
  }
  
  plot.list[[1]] + plot.list[[2]] + plot.list[[3]]
  
}

plot.list <- list()
plot.date.data("2019-07-30",
               color.vector = color.vector.003)

plot.list <- list()
plot.date.data("2019-08-30",
               color.vector = color.vector.005)

plot.list <- list()
plot.date.data("2019-10-08",
               color.vector = color.vector.006)

