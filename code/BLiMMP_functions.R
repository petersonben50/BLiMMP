
library(fields)
library(lubridate)
library(readxl)
library(rLakeAnalyzer)
library(stringr)
library(tidyverse)
library(viridis)

cb.translator <- readRDS("references/colorblind_friendly_colors.rds")
detection.vector <- c(4, 16)



#### Plot temperature from buoy over year as heatmap ####
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
    filter(year(sampledate) == year.of.interest) %>%
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
  par(mar = c(3, 3.5, 3.5, 1),
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
  title(main = year.of.interest)
  
}


#### Plot dissolved oxygen (from MeMO data and Mark Gahler's exo profiles) over year as a heatmap ####
DO.heatmap <- function(DO.profile.data = DO.data,
                       year.of.interest,
                       date.limits = c("2020-05-01", "2020-08-01")) {
  # Generate a matrix of the DO values
  clean.DO.profile <- DO.profile.data %>%
    filter(year(sampleDate) == year.of.interest) 
  
  clean.DO.profile.matrix <- clean.DO.profile%>%
    select(-sampleDate) %>%
    as.matrix() %>%
    image.smooth()
  
  # Add daynum and depth as x and y, respectively
  clean.DO.profile.matrix$x <- yday(clean.DO.profile$sampleDate)
  clean.DO.profile.matrix$y <- names(clean.DO.profile)[-1] %>% gsub("X", "", .) %>% as.numeric()
  
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
                 zlim = c(0, 15),
                 nlevels = 40,
                 color.palette = colorRampPalette(c(cb.translator["black"], cb.translator["blue"], cb.translator["bluishgreen"], cb.translator["orange"]),
                                                  bias = 1, space = "rgb"),
                 plot.axes = { 
                   axis(1,
                        at = dates.to.plot.on.axes %>% as.Date() %>% yday(),
                        labels = paste(month(dates.to.plot.on.axes, label = TRUE, abbr = FALSE),
                                       "1st"));
                   axis(2,
                        at = depths.to.plot.on.axes);
                   
                 },
                 key.title = {par(cex.main=1);title(main = "DO (mg/L)")},
                 key.axes = axis(4, seq(0, 15, by = 3),asp=1)
  )
  title(main = year(date.limits[1]))
}


#### Define function: Plot exo data from this study ####
plot.exo.data <- function(exo.data.to.use,
                          date.of.sampling,
                          points.of.sampling = NULL,
                          # main.title = NULL,
                          legend.location = "topright") {
  # Isolate needed exo data
  exo.data.for.profile <- exo.data.to.use %>%
    filter(date == date.of.sampling)
  
  # Plot the temp points and set up graph
  plot(x = (exo.data.for.profile$Temp_C-5)/2.5,
       y = -exo.data.for.profile$depth,
       type = "l",
       lwd = 3,
       xaxt = "n",
       xlab = "",
       xlim = c(0, 10),
       yaxt = "n",
       ylab = "",
       ylim = c(-25, 2),
       col = cb.translator["blue"])
  # Plot DO
  DO.fudge.factor <- 15
  lines(x = exo.data.for.profile$ODO_sat/DO.fudge.factor,
        y = -exo.data.for.profile$depth,
        col = cb.translator["black"],
        lwd = 3,
        lty = 4)
  # Plot the turbidity values
  turb.fudge.factor <- 0.5
  lines(x = exo.data.for.profile$Turbidity_FNU/turb.fudge.factor,
        y = -exo.data.for.profile$depth,
        col = cb.translator["reddishpurple"],
        lwd = 4,
        lty = 3)
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
  
  # if (is.null(main.title)) {
  #   title(main = date.of.sampling)
  # } else {
  #   title(main = main.title)
  # }
  # 
  
  # OPTIONAL: Add in sampling points
  if (!is.null(points.of.sampling)) {
    points(x = rep(9, length(points.of.sampling)),
           y = -points.of.sampling,
           pch = 16,
           cex = 2,
           col = cb.translator["bluishgreen"])
    legend(legend.location,
           legend = c("Temperature", "DO", "Turbidity", "Incubations"),
           col = c(cb.translator["blue"], cb.translator["black"], cb.translator["reddishpurple"], cb.translator["bluishgreen"]),
           lwd = c(3, 3, 3, NA),
           lty = c(1, 4, 3, NA),
           pt.cex = c(NA, NA, NA, 1.75),
           pch = c(NA, NA, NA, 16),
           bg = "white")
  } else {
    legend(legend.location,
           legend = c("Temperature", "DO", "Turbidity"),
           col = c(cb.translator["blue"], cb.translator["black"], cb.translator["reddishpurple"]),
           lwd = 3,
           lty = c(1, 4, 3, NA),
           pch = NA,
           bg = "white")
  }
  text(x = 4,
       y = 1.5,
       cex = 1.2,
       labels = paste(month(date.of.sampling, label = TRUE), " ",
                      day(date.of.sampling), ", ",
                      year(date.of.sampling),
                      sep = ""))
  
}


#### Define function: Plot metals data ####
plot.chem.data <- function(icp.data.to.use,
                             color.vector.to.use,
                             point.vector.to.use,
                             line.vector.to.use,
                             naming.vector.to.use,
                             date.of.sampling,
                             xlim.to.use = c(0, 0.3),
                             xlab.to.use = "Metals (ppm)",
                             legend.location = "topright") {
  
  # Set up data
  icp.data.to.use.filtered <- icp.data.to.use %>%
    filter(date == date.of.sampling)
  icp.data.to.use.filtered <- icp.data.to.use.filtered[, c("depth", names(color.vector.to.use))]
  icp.data.to.use.filtered <- icp.data.to.use.filtered %>%
    gather(key = constituent,
           value = concentration,
           -1) %>%
    filter(!is.na(concentration))
  
  # Metals data
  plot(x = icp.data.to.use.filtered$concentration,
       y = as.numeric(icp.data.to.use.filtered$depth),
       ylim = c(25, -2),
       xlim = xlim.to.use,
       col = color.vector.to.use[icp.data.to.use.filtered$constituent],
       pch = point.vector.to.use[icp.data.to.use.filtered$constituent],
       xlab = "",
       ylab = "")
  
  for (metal in names(color.vector.to.use)) {
    line.data <- icp.data.to.use.filtered %>%
      filter(constituent == metal) %>%
      group_by(depth) %>%
      summarise(concentration = mean(concentration)) %>%
      arrange(depth)
    lines(x = line.data$concentration,
          y = line.data$depth,
          col = color.vector.to.use[metal],
          lty = line.vector.to.use[metal])
  }
  
  legend(legend.location,
         legend = naming.vector.to.use,
         pch = point.vector.to.use,
         lty = line.vector.to.use,
         col= color.vector.to.use)
  

  title(xlab = xlab.to.use,
        ylab = "Depth (m)",
        line = 1.5)
}


#### Define function: Plot sulfide and sulfate data ####
plot.sulfide.sulfate.data <- function(sulfide.data.to.use,
                                      sulfate.data.to.use,
                                      xlim.to.use = c(0, 250),
                                      legend.location = "topright") {
  
  waterChem.data
  # Sulfide data
  plot(x = sulfide.data.to.use$S_conc_uM,
       y = sulfide.data.to.use$depth,
       pch = 16,
       ylim = c(25, -2),
       xlim = xlim.to.use,
       col = cb.translator["blue"],
       xlab = "",
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
  
  legend(legend.location,
         legend = c("Sulfate (µM)", "Sulfide (µM)"),
         pch = c(1, 16),
         lty = c(1, 2),
         col = cb.translator["blue"],
         bg = "white")
  title(xlab = "Sulfide/sulfate (µM)",
        ylab = "Depth (m)",
        line = 1.5)
  
}


#### Function to plot Hg data ####
plot.Hg.profile <- function(date.to.use = "2021-09-10",
                            Hg.data.to.use = Hg.data,
                            phase = "filtered") {
  
  # Prepare Hg data
  Hg.data.date <- Hg.data.to.use %>%
    filter(sampleDate == date.to.use)
  
  if (phase == "filtered") {
    # Split out values
    FMHG.values <- Hg.data.date %>%
      filter(constituent == "FMHG_NG.L") %>%
      arrange(depth) %>%
      select(depth, concentration)
    FMHG.values.no.dup <- FMHG.values %>%
      group_by(depth) %>%
      summarise(concentration = mean(concentration))
    FiHG.values <- Hg.data.date %>%
      filter(constituent == "FiHg_NG.L") %>%
      arrange(depth) %>%
      select(depth, concentration)
    FiHG.values.no.dup <- FiHG.values %>%
      group_by(depth) %>%
      summarise(concentration = mean(concentration))
    PerFMHG.values <- Hg.data.date %>%
      filter(constituent == c("perFMHG"))
    PerFMHG.values.no.dup <- PerFMHG.values %>%
      group_by(depth) %>%
      summarise(concentration = mean(concentration))
  }
  
  if (phase == "filtered") {
    Hg.adjustment.factor <- 1
    Hg.label <- "Dissolved Hg (ng/L)"
    # Split out values
    MHG.values <- Hg.data.date %>%
      filter(constituent == "FMHG_NG.L") %>%
      arrange(depth) %>%
      select(depth, concentration)
    MHG.values.no.dup <- MHG.values %>%
      group_by(depth) %>%
      summarise(concentration = mean(concentration))
    iHG.values <- Hg.data.date %>%
      filter(constituent == "FiHg_NG.L") %>%
      arrange(depth) %>%
      select(depth, concentration)
    iHG.values.no.dup <- iHG.values %>%
      group_by(depth) %>%
      summarise(concentration = mean(concentration))
    PerMHG.values <- Hg.data.date %>%
      filter(constituent == c("perFMHG"))
    PerMHG.values.no.dup <- PerMHG.values %>%
      group_by(depth) %>%
      summarise(concentration = mean(concentration))
  } else if (phase == "particulate") {
    # Split out values
    Hg.adjustment.factor <- 1
    Hg.label <- "Particulate Hg (ng/L)"
    MHG.values <- Hg.data.date %>%
      filter(constituent == "PMHG_NG.L") %>%
      arrange(depth) %>%
      select(depth, concentration)
    MHG.values.no.dup <- MHG.values %>%
      group_by(depth) %>%
      summarise(concentration = mean(concentration))
    iHG.values <- Hg.data.date %>%
      filter(constituent == "PiHg_NG.L") %>%
      arrange(depth) %>%
      select(depth, concentration)
    iHG.values.no.dup <- iHG.values %>%
      group_by(depth) %>%
      summarise(concentration = mean(concentration))
    PerMHG.values <- Hg.data.date %>%
      filter(constituent == c("perPMHG"))
    PerMHG.values.no.dup <- PerMHG.values %>%
      group_by(depth) %>%
      summarise(concentration = mean(concentration))
  } else if (phase == "particulate-SPM") {
    Hg.adjustment.factor <- 0.005
    Hg.label <- "Particulate Hg (ng/g)"
    MHG.values <- Hg.data.date %>%
      filter(constituent == "PMHG_NG.G") %>%
      arrange(depth) %>%
      select(depth, concentration) %>%
      mutate(concentration = concentration)
    MHG.values.no.dup <- MHG.values %>%
      group_by(depth) %>%
      summarise(concentration = mean(concentration))
    iHG.values <- Hg.data.date %>%
      filter(constituent == "PiHG_NG.G") %>%
      arrange(depth) %>%
      select(depth, concentration) %>%
      mutate(concentration = concentration)
    iHG.values.no.dup <- iHG.values %>%
      group_by(depth) %>%
      summarise(concentration = mean(concentration))
    PerMHG.values <- Hg.data.date %>%
      filter(constituent == c("perPMHG"))
    PerMHG.values.no.dup <- PerMHG.values %>%
      group_by(depth) %>%
      summarise(concentration = mean(concentration))
    
  }
  
  # Plot MeHg values
  plot(x = MHG.values$concentration*Hg.adjustment.factor,
       y = MHG.values$depth,
       ylim = c(25, -2),
       xlim = c(0, 1.5),
       pch = 17,
       cex = 1.5,
       col = cb.translator["vermillion"],
       xaxt = "n",
       xlab = "",
       yaxt = "n",
       ylab = "")
  lines(MHG.values.no.dup$concentration*Hg.adjustment.factor,
        MHG.values.no.dup$depth,
        col = cb.translator["vermillion"])
  
  # Plot iHg values
  points(x = iHG.values$concentration*Hg.adjustment.factor,
         y = iHG.values$depth,
         pch = 18,
         cex = 1.5,
         col = cb.translator["bluishgreen"])
  lines(x = iHG.values.no.dup$concentration*Hg.adjustment.factor,
        y = iHG.values.no.dup$depth,
        col = cb.translator["bluishgreen"])
  
  # Plot percent MeHg values
  points(x = PerMHG.values$concentration*0.015,
         y = PerMHG.values$depth,
         pch = 16,
         cex = 1.5,
         col = cb.translator["black"])
  lines(x = PerMHG.values.no.dup$concentration*0.015,
        y = PerMHG.values.no.dup$depth,
        col = cb.translator["black"])
  
  # Add y-axis with depth measurements
  axis(2,
       at = seq(0, 25, by = 5),
       labels = seq(0, 25, by = 5))
  # Add axis for Hg
  axis(1,
       at = seq(0, 1.5, by = 0.3),
       labels = seq(0, 1.5, by = 0.3)/Hg.adjustment.factor)
  axis(1,
       line = 3,
       at = seq(0, 1.5, by = 0.375),
       labels = seq(0, 1.5, by = 0.375)/1.5)
  
  # Add labels
  # Add label for percent MeHg
  title(xlab = "Fraction MeHg",
        line = 4.5)
  # Add label for depth
  title(ylab = "Depth (m)",
        line = 1.5)
  # Add label for sulfide
  title(xlab = Hg.label,
        line = 1.5)
  
  
  # Add legend
  legend("topright",
         legend = c("iHg",
                    "MeHg",
                    "% MeHg"),
         col = c(cb.translator["bluishgreen"],
                 cb.translator["vermillion"],
                 cb.translator["black"]),
         pch = c(18, 17, 16),
         pt.cex = c(1.75, 1.5, 1.75),
         bg = "white")
}











#### Function to plot change in Me198Hg for each depth on a given trip ####
plot.Hg.time.course.depth <- function(selected.depth,
                                      trip_data,
                                      parameter.of.interest,
                                      color.vector.input = NULL,
                                      treatment.names.vector = NULL,
                                      DDL.column = NULL,
                                      y.label.to.use = NULL,
                                      ylim.to.use = NULL,
                                      legend.position = "topleft") {
  
  par(mar=c(3,3,2.5,1), mgp=c(1.5,0.4,0), tck=-0.008)
  
  if (is.null(y.label.to.use)) {
    y.label.to.use <- parameter.of.interest
  }
  
  trip_data[, "parameter.to.use"] <- trip_data[, parameter.of.interest]
  trip_data <- trip_data %>%
    filter(!is.na(parameter.to.use))
  
  if (is.null(ylim.to.use)) {
    y_max <- ceiling(max(trip_data$parameter.to.use)*100)/100 + 0.01
    y_min <- floor(min(trip_data$parameter.to.use)*100)/100
    ylim.to.use <- c(y_min, y_max)
  }

  
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
  
  
  #### Generate shape vector ####
  if (!is.null(DDL.column)) {
    trip_data_depth[, "DDL_shape"] <- rep(4, length(trip_data_depth$parameter.to.use))
    trip_data_depth$DDL_shape[trip_data_depth[, DDL.column]] <- 18
  } else {
    trip_data_depth[, "DDL_shape"] <- rep(18, length(trip_data_depth$parameter.to.use))
  }
  
  
  plot(x = trip_data_depth$durationInDays,
       y = trip_data_depth$parameter.to.use,
       ylim = ylim.to.use,
       xlim = c(-0.5, 4),
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
       at = c(0, 1, 3.5),
       labels = c("t0", "24hr", "84hr"))
  
  #### Add segments to link samples within an incubation ####
  
  trip_data_depth_segments <- trip_data_depth %>%
    select(c(incubationID, treatment, t, parameter.to.use)) %>%
    spread(key = t,
           value = parameter.to.use) %>%
    left_join(trip_data_depth %>%
                select(c(incubationID, treatment, t, durationInDays)) %>%
                mutate(t = paste(t, "_duration", sep = "")) %>%
                spread(key = t,
                       value = durationInDays))
  
  for (row.num in 1:nrow(trip_data_depth_segments)) {
    segments(trip_data_depth_segments$t0_duration[row.num],
             trip_data_depth_segments$t0[row.num],
             trip_data_depth_segments$t1_duration[row.num],
             trip_data_depth_segments$t1[row.num],
             col = color.vector[trip_data_depth_segments$treatment[row.num]])
    if (!is.null(trip_data_depth_segments$t2_duration)) {
      segments(trip_data_depth_segments$t1_duration[row.num],
               trip_data_depth_segments$t1[row.num],
               trip_data_depth_segments$t2_duration[row.num],
               trip_data_depth_segments$t2[row.num],
               col = color.vector[trip_data_depth_segments$treatment[row.num]])
    }
    
  }
  
  title(ylab = y.label.to.use,
        line = 1.75)
  
  if (is.null(treatment.names.vector)) {
    legend(legend.position,
           legend =  names(color.vector)[length(color.vector):1],
           text.col = color.vector[length(color.vector):1],
           bty = "n")
  } else {
    legend(legend.position,
           legend =  treatment.names.vector[names(color.vector)[length(color.vector):1]],
           text.col = color.vector[length(color.vector):1],
           bty = "n")
  }
  
}



#### Function to plot change in parameter over time ####
plot.production.of.parameter <- function(dataFile.of.interest,
                                         tripID.of.interest,
                                         t.of.interest,
                                         named.color.vector,
                                         parameter.of.interest,
                                         treatment.names.vector = NULL,
                                         plot.title = "",
                                         ylims.to.use = NULL,
                                         y.label.to.use = NULL) {
  
  dataFile.of.interest[, "parameter.to.plot"] <- dataFile.of.interest[, parameter.of.interest]
  
  if (is.null(y.label.to.use)) {
    y.label.to.use = parameter.of.interest
  }
  
  plot.to.plot <- dataFile.of.interest %>%
    filter(tripID == tripID.of.interest,
           t == t.of.interest) %>%
    mutate(treatment = fct_relevel(treatment, names(named.color.vector)),
           depth = paste(depth, "meters")) %>%
    ggplot(aes(x = treatment,
               y = parameter.to.plot,
               color = treatment)) +
    geom_point(stat = "identity") +
    scale_color_manual(values = color.vector,
                       labels = treatment.names.vector) +
    facet_grid(~depth) +
    scale_x_discrete(labels = gsub(pattern = " ",
                                   "\n",
                                   treatment.names.vector[names(color.vector)])) +
    ylab(y.label.to.use) +
    xlab("") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust=1)) +
    theme(legend.position = "none") +
    ggtitle(plot.title)
  
  if (!is.null(ylims.to.use)) {
    plot.to.plot <- plot.to.plot +
      ylim(ylims.to.use)
  }
  plot.to.plot
}




#### Define function to plot multiple proteins ####
plot.profile.for.multiple.genes <- function(marker.depth.df,
                                            genesOfInterest,
                                            #yearOfInterest,
                                            monthOfInterest,
                                            gene.name.column = "geneName",
                                            xlab.to.use = "Gene coverage normalized\nto SCG coverage (100X)",
                                            depth_limits = c(25, 0),
                                            coverage_limits = NULL,
                                            show.mean.coverage = TRUE,
                                            color.vector.to.use = NULL,
                                            point.vector.to.use = NULL,
                                            legend.position.to.use = "default",
                                            legend.title.to.use = element_blank(),
                                            titleToUse = NULL) {
  
  marker.depth.df[, "gene.names.to.use.for.filtering"] <- marker.depth.df[, gene.name.column]
  
  clean.coverage <- marker.depth.df %>%
    filter(gene.names.to.use.for.filtering %in% genesOfInterest) %>%
    # filter(year(ymd(startDate)) == yearOfInterest) %>%
    filter(month(ymd(startDate)) == monthOfInterest) %>%
    group_by(gene.names.to.use.for.filtering, depth) %>%
    summarise(coverage = sum(coverage)) %>%
    ungroup() %>%
    arrange(depth)
  
  graph.of.interest <- clean.coverage %>%
    ggplot(aes(x = depth,
               y = coverage)) +
    theme_classic() +
    scale_x_reverse(limits = depth_limits)  +
    theme(axis.text.x = element_text(colour="black"),
          axis.text.y = element_text(colour="black"))
  
  #### Add labels and title ####
  if (is.null(titleToUse)) {
    titleToUse = paste("Coverage of",
                       paste(genesOfInterest, collapse = ","))
    
  }
  graph.of.interest <- graph.of.interest +
    labs(title = titleToUse,
         y = xlab.to.use,
         x = "Depth (m)")
  
  #### Add colors if defined in the call ####
  if (!is.null(color.vector.to.use)) {
    graph.of.interest <- graph.of.interest +
      geom_line(aes(group = gene.names.to.use.for.filtering,
                    color = gene.names.to.use.for.filtering,
                    linetype = gene.names.to.use.for.filtering)) +
      geom_point(aes(color = gene.names.to.use.for.filtering,
                     shape = gene.names.to.use.for.filtering)) +
      scale_colour_manual(values=color.vector.to.use)
  } else {
    graph.of.interest <- graph.of.interest +
      geom_point(shape = gene.names.to.use.for.filtering) +
      geom_line(aes(group = gene.names.to.use.for.filtering,
                    linetype = gene.names.to.use.for.filtering))
  }
  
  #### Add point shapes if defined in the call ####
  if (!is.null(point.vector.to.use)) {
    graph.of.interest <- graph.of.interest +
      scale_shape_manual(values = point.vector.to.use)
  } 
  #### Constrain x-axis if defined in the call ####
  if (!is.null(coverage_limits)) {
    graph.of.interest <- graph.of.interest +
      scale_y_continuous(limits = coverage_limits)
  } else {
    max.coverage <- clean.coverage %>%
      select(coverage) %>%
      unlist() %>%
      max()
    graph.of.interest <- graph.of.interest +
      scale_y_continuous(limits = c(0, max.coverage))
    
  }
  
  if (show.mean.coverage == TRUE) {
    graph.of.interest <- graph.of.interest +
      stat_summary(geom = "point", fun = "mean",
                   col = "black", fill = "red",
                   size = 3, shape = 24)
  }
  
  #### Set up legend position ####
  if (legend.position.to.use[1] != "default") {
    graph.of.interest <- graph.of.interest +
      theme(legend.position = legend.position.to.use,
            legend.title = legend.title.to.use,
            legend.box.background = element_rect(colour = "black"))
  }
  
  graph.of.interest + coord_flip()
  
}


