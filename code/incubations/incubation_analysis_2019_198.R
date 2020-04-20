#### code/incubations/incubation_analysis_2019.R ####
# Written for BLiMMP project
# Benjamin D. Peterson

#### Prep workspace ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP/")
library(dplyr)
library(lubridate)
library(readxl)
library(tidyr)



#### Read in MeHg data ####

MeHg.data <- read.csv("dataEdited/incubations/MeHg/incubations2019_Me198Hg.csv",
                      stringsAsFactors = FALSE)


#### Read in colors ####

cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")


#### Function to plot exo data ####

plot.exo.data <- function(date.of.sampling,
                          points.of.sampling = NULL) {
  
  # Read in exo data
  exo.data.file.name <- paste("dataEdited/exo/",
                              date.of.sampling,
                              "_profile.csv",
                              sep = "")
  exo.data <- read.csv(exo.data.file.name,
                       stringsAsFactors = FALSE)
  
  # Generate plot
  par(mar=c(9,3,2.5,1), mgp=c(1.5,0.4,0), tck=-0.008)
  
  # Plot the temp points and set up graph
  plot(x = (exo.data$Temp_C-5)/2.5,
       y = -exo.data$Depth_m,
       type = "l",
       lwd = 3,
       xaxt = "n",
       xlab = "",
       xlim = c(0, 10),
       yaxt = "n",
       ylab = "",
       ylim = c(-25, 0),
       col = cb.translator["blue"])
  # Plot DO
  DO.fudge.factor <- 15
  lines(x = exo.data$ODO_sat/DO.fudge.factor,
        y = -exo.data$Depth_m,
        col = cb.translator["black"],
        lwd = 3)
  # Plot the turbidity values
  turb.fudge.factor <- 1.5
  lines(x = exo.data$Turbidity_FNU/turb.fudge.factor,
        y = -exo.data$Depth_m,
        col = cb.translator["reddishpurple"],
        lwd = 3)
  # Add y-axis with depth measurements
  axis(2,
       at = seq(0, -25, by = -5),
       labels = seq(0, 25, by = 5))
  # Add axis for temperature
  axis(1,
       at = seq(0, 10, by = 2),
       labels = seq(0, 10, by = 2)*3)
  # Add axis for DO values
  axis(1,
       at = seq(0, 10, by = 2),
       labels = seq(0, 10, by = 2)*DO.fudge.factor,
       line = 3,
       par(mgp = c(3, 0.5, 0)))
  # Add axis for turbidity values
  axis(1,
       at = seq(0, 10, by = 2),
       labels = seq(0, 10, by = 2)*turb.fudge.factor,
       line = 6,
       par(mgp = c(3, 0.5, 0)))
  
  
  # Add titles for depth
  title(xlab = "Temp (C)",
        ylab = "Depth (m)",
        line = 1.5)
  # Add title for dissolved oxygen
  title(xlab = "ODO (% sat)",
        line = 4.5)
  # Add title for turbidity
  title(xlab = "Turbidity (FNU)",
        line = 7.5)
  title(main = date.of.sampling)
  
  # OPTIONAL: Add in sampling points
  if (!is.null(points.of.sampling)) {
    points(x = rep(5, length(points.of.sampling)),
           y = points.of.sampling,
           pch = 16,
           col = cb.translator["bluishgreen"])
  }
  
  legend("bottomright",
         legend = c("Temperature", "DO", "Turbidity", "Sampling Points"),
         col = c(cb.translator["blue"], cb.translator["black"], cb.translator["reddishpurple"], cb.translator["bluishgreen"]),
         lwd = c(3, 3, 3, NA),
         pch = c(NA, NA, NA, 16),
         bty = "n")
  
}

#### Function to plot redox profiles ####

plot.redox.profile <- function(trip) {
  
  # Read in sulfide data
  sulfide.data.name <- "dataEdited/waterChemistry/sulfide/WC_data.csv"
  sulfide.data.date <- read.csv(sulfide.data.name,
                           stringsAsFactors = FALSE) %>%
    filter(tripID == trip)
  
  # Read in ICP data
  ICP.data.name <- "dataEdited/waterChemistry/ICP/metal_data.csv"
  Mn.data.date <- read.csv(ICP.data.name,
                                stringsAsFactors = FALSE) %>%
    filter(tripID == trip,
           element == "Mn")
  Fe.data.date <- read.csv(ICP.data.name,
                           stringsAsFactors = FALSE) %>%
    filter(tripID == trip,
           element == "Fe")
  
  
  
  # Generate plot, add sulfide data points
  sulfide.adjustment.for.plotting <- 0.1
  plot(x = sulfide.data.date$S_conc*sulfide.adjustment.for.plotting,
       y = -sulfide.data.date$depth,
       ylim = c(-25, 0),
       xlim = c(-1, 10),
       pch = 16,
       col = cb.translator["blue"],
       xaxt = "n",
       xlab = "",
       yaxt = "n",
       ylab = "")
  
  # Add Mn data points
  Mn.adjustment.for.plotting <- 10
  points(x = Mn.data.date$diss_conc_ppm*Mn.adjustment.for.plotting,
         y = -Mn.data.date$depth,
         pch = 18,
         col = cb.translator["orange"])
  
  # Add y-axis with depth measurements
  axis(2,
       at = seq(0, -25, by = -5),
       labels = seq(0, 25, by = 5))

  # Add axis for sulfide
  axis(1,
       at = seq(0, 10, by = 2),
       labels = seq(0, 10, by = 2)/sulfide.adjustment.for.plotting)

  # Add axis for Mn
  axis(1,
       line = 3,
       at = seq(0, 10, by = 2),
       labels = seq(0, 10, by = 2)/sulfide.adjustment.for.plotting)
  
  # Add label for depth
  title(ylab = "Depth (m)",
        line = 1.5)
    
  # Add label for sulfide
  title(xlab = "Sulfide (uM)",
        line = 1.5)
  
  # Add label for Mn
  title(xlab = "Mn (ppm)",
        line = 4.5)
  
  # Add legend
  legend("topright",
         legend = c("Sulfide",
                    "Manganese"),
         col = c(cb.translator["blue"],
                 cb.translator["orange"]),
         pch = 18)
}

#### Function to plot bacterial data ####

plot.bacterial.profile <- function(date.of.sampling) {

  leu.data.name <- "dataEdited/leucineUptake/water_column_uptake.csv"
  leu.data.date <- read.csv(leu.data.name,
                            stringsAsFactors = FALSE) %>%
    filter(startDate == date.of.sampling)
  
  leucine.adj.for.plotting <- 100
  
  plot(x = NA,
       y = NA,
       ylim = c(25, 0),
       xlim = c(-1, 10),
       xaxt = "n",
       xlab = "",
       yaxt = "n",
       ylab = "")
  
  points(x = leu.data.date$Leu_uptake_pM_A/leucine.adj.for.plotting,
         y = leu.data.date$depth,
         pch = 18,
         col = cb.translator["bluishgreen"])
  
  # Add y-axis with depth measurements
  axis(2,
       at = seq(0, 25, by = 5),
       labels = seq(0, 25, by = 5))
  
  # Add axis for sulfide
  axis(1,
       at = seq(0, 10, by = 2),
       labels = seq(0, 10, by = 2)*leucine.adj.for.plotting)
  
  # Add label for depth
  title(ylab = "Depth (m)",
        line = 1.5)
  
  # Add label for leucine
  title(xlab = "Leucine uptake (pM)",
        line = 1.5)
  }
  

#### Function to generate empty plot ####
empty.plot <- function(x) {
  plot(x = 0,
       cex = 0,
       xaxt = "n",
       yaxt = "n",
       xlab = "",
       bty = "n")
}

#### Function to plot change in Me198Hg for each depth on a given trip ####

plot.MeHg.time.course.depth <- function(selected.depth,
                                        trip_data,
                                        color.vector.input = NULL) {

  par(mar=c(3,3,2.5,1), mgp=c(1.5,0.4,0), tck=-0.008)
  
  y_max <- ceiling(max(trip_data$excess_MeHg_198_ng.L)*100)/100 + 0.01
  y_min <- floor(min(trip_data$excess_MeHg_198_ng.L)*100)/100
  
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
       y = trip_data_depth$excess_MeHg_198_ng.L,
       ylim = c(y_min, y_max),
       xlim = c(-0.5, 1.5),
       pch = trip_data_depth$DDL_shape,
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
    select(-c(above_DDL, DDL_shape, depth)) %>%
    spread(key = t,
           value = excess_MeHg_198_ng.L)
  
  for (row.num in 1:nrow(trip_data_depth_segments)) {
    segments(0, trip_data_depth_segments$t0[row.num],
             1, trip_data_depth_segments$t1[row.num],
             col = color.vector[trip_data_depth_segments$treatment[row.num]])
  }
  
  title(ylab = "MeHg production (ng/L per day)",
        line = 1.75)
  
  legend("topleft",
         legend =  names(color.vector)[length(color.vector):1],
         text.col = color.vector[length(color.vector):1],
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
       main = paste(selected.depth,
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
  
  title(ylab = "MeHg production (ng/L per day)",
        line = 1.75)
  
  legend("topleft",
         legend =  names(color.vector)[length(color.vector):1],
         text.col = color.vector[length(color.vector):1],
         bty = "n")
}









#### Plot all data using functions ####
trip.of.interest <- "BLiMMP_trip_006"
plot.incubation.data <- function(trip.of.interest,
                                 plot.exo = TRUE,
                                 should.i.plot.geochem.data = "yes",
                                 should.i.plot.bacterial.data = "yes",
                                 empty.plots = NULL,
                                 plot.Hg.data = TRUE,
                                 color.vector.for.Hg = NULL) {
  
  date.of.interest <- read_xlsx("metadata/1_trip_IDs.xlsx") %>%
    filter(tripID == trip.of.interest) %>%
    select(startDate) %>%
    unlist(use.names = FALSE)
  
  trip_data_df <- MeHg.data %>%
    filter(tripID == trip.of.interest) %>%
    select(incubationID, depth, t, treatment, excess_MeHg_198_ng.L, above_DDL)
  
  # Sampling locations:
  sampling.depths <- unique(trip_data_df$depth)
  
  # Plot exo data. Option to not plot it if we don't have the data. 
  if (plot.exo == TRUE) {
    plot.exo.data(date.of.interest,
                  points.of.sampling = -sampling.depths)
  }
  
  # Plot redox data
  
  if (should.i.plot.geochem.data == "yes") {
    plot.redox.profile(trip.of.interest)
  }
  
  if (should.i.plot.bacterial.data == "yes") {
    plot.bacterial.profile(date.of.interest)
    }
  
  if (!is.null(empty.plots)) {
    for (number in 1:empty.plots) {
      empty.plot()
    }
  }
  
  # Plot all the incubation information
  if (plot.Hg.data == TRUE) {
    sapply(X = sort(unique(trip_data_df$depth)),
           function(x) {
             plot.MeHg.time.course.depth(x,
                                         trip_data = trip_data_df,
                                         color.vector = color.vector.for.Hg)
           }
    )
  }
  
}



#### Generate plots ####
color.vector.003 <- c(cb.translator["bluishgreen"],
                      cb.translator["vermillion"])
names(color.vector.003) <- c("unfiltered-unamended",
                             "filtered-unamended")

# png("results/incubations/BLiMMP_trip_003_incubations.png",
#     units = "in",
#     res = 240,
#     height = 10,
#     width = 6)
pdf("results/incubations/BLiMMP_trip_003_incubations.pdf",
    height = 6,
    width = 6,
    useDingbats = FALSE)
par(mfrow = c(1, 2))
plot.incubation.data(trip.of.interest = "BLiMMP_trip_003",
                     plot.exo = FALSE,
                     should.i.plot.geochem.data = "no",
                     should.i.plot.bacterial.data = "no",
                     color.vector.for.Hg = color.vector.003)
dev.off()


# png("results/incubations/BLiMMP_trip_005_incubations.png",
#     units = "in",
#     res = 240,
#     height = 10,
#     width = 6)
pdf("results/incubations/BLiMMP_trip_005_incubations.pdf",
    height = 10,
    width = 6,
    useDingbats = FALSE)
par(mfrow = c(2, 3))
plot.incubation.data(trip.of.interest = "BLiMMP_trip_005",
                     plot.exo = TRUE,
                     empty.plots = NULL)
dev.off()


# Set color vector
# unique(plotting.data$treatment)
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

# png("results/incubations/BLiMMP_trip_006_incubations.png",
#     units = "in",
#     res = 240,
#     height = 10,
#     width = 6)
pdf("results/incubations/BLiMMP_trip_006_incubations.pdf",
    height = 10,
    width = 6,
    useDingbats = FALSE)
par(mfrow = c(2, 3))
par(mfrow = c(2, 3))
plot.incubation.data(trip.of.interest = "BLiMMP_trip_006",
                     plot.exo = TRUE,
                     empty.plots = NULL,
                     color.vector.for.Hg = color.vector.006)
dev.off()
