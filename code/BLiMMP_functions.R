






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
       y = -exo.data$depth,
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
        y = -exo.data$depth,
        col = cb.translator["black"],
        lwd = 3,
        lty = 4)
  # Plot the turbidity values
  turb.fudge.factor <- 0.5
  lines(x = exo.data$Turbidity_FNU/turb.fudge.factor,
        y = -exo.data$depth,
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
  title(main = date.of.sampling)
  
  # OPTIONAL: Add in sampling points
  if (!is.null(points.of.sampling)) {
    points(x = rep(9, length(points.of.sampling)),
           y = -points.of.sampling,
           pch = 16,
           cex = 2,
           col = cb.translator["bluishgreen"])
    legend("topright",
           legend = c("Temperature", "DO", "Turbidity", "Incubations"),
           col = c(cb.translator["blue"], cb.translator["black"], cb.translator["reddishpurple"], cb.translator["bluishgreen"]),
           lwd = c(3, 3, 3, NA),
           lty = c(1, 4, 3, NA),
           pt.cex = c(NA, NA, NA, 1.75),
           pch = c(NA, NA, NA, 16),
           bg = "white")
  } else {
    legend("topright",
           legend = c("Temperature", "DO", "Turbidity"),
           col = c(cb.translator["blue"], cb.translator["black"], cb.translator["reddishpurple"]),
           lwd = 3,
           lty = c(1, 4, 3, NA),
           pch = NA,
           bg = "white")
  }
}



#### Function to plot geochem data ####
plot.redox.profile <- function(trip,
                               sulfide.data.location,
                               sulfate.data.location,
                               ICP.data.name,
                               plot.fe = "yes") {
  
  # Read in sulfide data
  sulfide.data.date <- read.csv(sulfide.data.location,
                                stringsAsFactors = FALSE) %>%
    filter(tripID == trip) %>%
    mutate(depth = as.numeric(depth)) %>%
    arrange(depth)
  
  # Read in sulfate data
  sulfate.data.date <- read.csv(sulfate.data.location,
                                stringsAsFactors = FALSE) %>%
    filter(tripID == trip) %>%
    mutate(depth = as.numeric(depth)) %>%
    arrange(depth)
  
  # Read in ICP data
  Mn.data.date <- readRDS(ICP.data.name) %>%
    filter(tripID == trip,
           constituent == "dissolved_Mn_ppb") %>%
    mutate(diss_conc_uM = concentration/ 54.938) %>%
    mutate(depth = as.numeric(depth)) %>%
    arrange(depth)
  Mn.part.data.date <- readRDS(ICP.data.name) %>%
    filter(tripID == trip,
           constituent == "particulate_Mn_ppb") %>%
    mutate(part_conc_uM = concentration/ 54.938) %>%
    mutate(depth = as.numeric(depth)) %>%
    arrange(depth)

  Fe.data.date <- readRDS(ICP.data.name) %>%
    filter(tripID == trip,
           constituent == "dissolved_Fe_ppb") %>%
    mutate(diss_conc_uM = concentration/ 55.85) %>%
    mutate(depth = as.numeric(depth)) %>%
    arrange(depth)
  Fe.part.data.date <- readRDS(ICP.data.name) %>%
    filter(tripID == trip,
           constituent == "particulate_Fe_ppb") %>%
    mutate(diss_conc_uM = concentration/ 55.85) %>%
    mutate(depth = as.numeric(depth)) %>%
    arrange(depth)
  
  # Generate plot, add sulfide data points
  sulfide.adjustment.for.plotting <- 0.1
  plot(x = sulfide.data.date$S_conc*sulfide.adjustment.for.plotting,
       y = -sulfide.data.date$depth,
       ylim = c(-25, 0),
       xlim = c(-1, 23),
       pch = 16,
       cex = 1.8,
       col = cb.translator["blue"],
       xaxt = "n",
       xlab = "",
       yaxt = "n",
       ylab = "")
  lines(x = sulfide.data.date$S_conc*sulfide.adjustment.for.plotting,
        y = -sulfide.data.date$depth,
        col = cb.translator["blue"])

  # Add sulfate data points
  sulfate.adjustment.for.plotting <- 0.1
  points(x = sulfate.data.date$concentration/96.06*1000*sulfate.adjustment.for.plotting,
         y = -as.numeric(sulfate.data.date$depth),
         pch = 1,
         cex = 1.5,
         col = cb.translator["blue"])
  lines(x = sulfate.data.date$concentration/96.06*1000*sulfate.adjustment.for.plotting,
        y = -as.numeric(sulfate.data.date$depth),
        col = cb.translator["blue"])
  
  # Add dissolved Mn data points
  Mn.adjustment.for.plotting <- 3
  points(x = Mn.data.date$diss_conc_uM*Mn.adjustment.for.plotting,
         y = -as.numeric(Mn.data.date$depth),
         pch = 18,
         cex = 2,
         col = cb.translator["orange"])
  lines(x = Mn.data.date$diss_conc_uM*Mn.adjustment.for.plotting,
        y = -as.numeric(Mn.data.date$depth),
        col = cb.translator["orange"])
  
  # Add particulate Mn data points
  points(x = Mn.part.data.date$part_conc_uM*Mn.adjustment.for.plotting,
         y = -as.numeric(Mn.part.data.date$depth),
         pch = 5,
         cex = 1.5,
         col = cb.translator["orange"])
  lines(x = Mn.part.data.date$part_conc_uM*Mn.adjustment.for.plotting,
        y = -as.numeric(Mn.part.data.date$depth),
        col = cb.translator["orange"])

  # # Add dissolved Fe data points
  # Fe.adjustment.for.plotting <- 1
  # points(x = Fe.part.data.date$diss_conc_uM*Fe.adjustment.for.plotting,
  #        y = -as.numeric(Fe.part.data.date$depth),
  #        pch = 15,
  #        cex = 1.5,
  #        col = cb.translator["reddishpurple"])
  # lines(x = Fe.part.data.date$diss_conc_uM*Fe.adjustment.for.plotting,
  #       y = -as.numeric(Fe.part.data.date$depth),
  #       col = cb.translator["reddishpurple"])
  # 
  # # Add particulate Fe data points
  # Fe.adjustment.for.plotting <- 1
  # points(x = Fe.data.date$diss_conc_uM*Fe.adjustment.for.plotting,
  #        y = -as.numeric(Fe.data.date$depth),
  #        pch = 0,
  #        cex = 1.5,
  #        col = cb.translator["reddishpurple"])
  # lines(x = Fe.data.date$diss_conc_uM*Fe.adjustment.for.plotting,
  #       y = -as.numeric(Fe.data.date$depth),
  #       col = cb.translator["reddishpurple"])
  
  
  # Add y-axis with depth measurements
  axis(2,
       at = seq(0, -25, by = -5),
       labels = seq(0, 25, by = 5))
  
  # Add axis for sulfide
  axis(1,
       at = seq(0, 20, by = 5),
       labels = seq(0, 20, by = 5)/sulfide.adjustment.for.plotting)
  
  # Add axis for Mn
  axis(1,
       line = 3,
       at = seq(0, 21, by = 3),
       labels = seq(0, 21, by = 3)/Mn.adjustment.for.plotting)

  # Add label for depth
  title(ylab = "Depth (m)",
        line = 1.5)
  
  # Add label for sulfide
  title(xlab = "Sulfide/sulfate (uM)",
        line = 1.5)
  
  # Add label for Mn
  title(xlab = "Mn (µM)",
        line = 4.5)
  
  # Add legend
  legend("topright",
         legend = c("Sulfide",
                    "Sulfate",
                    "Dissolved Mn",
                    "Particulate Mn"),
         col = c(cb.translator["blue"],
                 cb.translator["blue"],
                 cb.translator["orange"],
                 cb.translator["orange"]),
         pch = c(16, 1, 18, 5),
         pt.cex = c(1.75, 1.75, 1.75, 1.25),
         bty = "n")
}


#### Function to plot Hg data ####
plot.Hg.profile <- function(tripDates,
                            Hg.data.location) {
  
  # Read in Hg data
  Hg.data.date <- read.csv(Hg.data.location) %>%
    filter(sampleDate %in% tripDates)
  
  # Split out values
  FMHG.values <- Hg.data.date %>%
    filter(constituent == "FMHG_NG.L") %>%
    arrange(depth) %>%
    select(depth, concentration)
  FiHG.values <- Hg.data.date %>%
    filter(constituent == "FiHg_NG.L") %>%
    arrange(depth) %>%
    select(depth, concentration)
  

  percent.MeHg <- Hg.data.date %>%
    spread(key = constituent,
           value = concentration) %>%
    mutate(perMeHg = FMHG_NG.L / FTHG_NG.L) %>%
    arrange(depth)

  # Plot MeHg values
  plot(x = FMHG.values$concentration,
       y = FMHG.values$depth,
       ylim = c(25, 0),
       xlim = c(0, 1.5),
       pch = 17,
       cex = 1.5,
       col = cb.translator["vermillion"],
       xaxt = "n",
       xlab = "",
       yaxt = "n",
       ylab = "")
  lines(FMHG.values$concentration,
        FMHG.values$depth,
        col = cb.translator["vermillion"])
  
  # Plot iHg values
  points(x = FiHG.values$concentration,
         y = FiHG.values$depth,
         pch = 18,
         cex = 1.5,
         col = cb.translator["bluishgreen"])
  lines(x = FiHG.values$concentration,
        y = FiHG.values$depth,
        col = cb.translator["bluishgreen"])
  
  #### Plot percent MeHg values ####
  points(x = percent.MeHg$perMeHg*1.5,
         y = percent.MeHg$depth,
         pch = 16,
         cex = 1.5,
         col = cb.translator["black"])
  lines(x = percent.MeHg$perMeHg*1.5,
        y = percent.MeHg$depth,
        col = cb.translator["black"])
  
  # Add y-axis with depth measurements
  axis(2,
       at = seq(0, 25, by = 5),
       labels = seq(0, 25, by = 5))
    # Add axis for Hg
  axis(1,
       at = seq(0, 1.5, by = 0.3),
       labels = seq(0, 1.5, by = 0.3))
  axis(1,
       line = 3,
       at = seq(0, 1.5, by = 0.375),
       labels = seq(0, 1.5, by = 0.375)/1.5)
  
  #### Add labels ####
  # Add label for percent MeHg
  title(xlab = "Fraction MeHg",
        line = 4.5)
  # Add label for depth
  title(ylab = "Depth (m)",
        line = 1.5)
  # Add label for sulfide
  title(xlab = "Hg (ng/L)",
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
         bty = "n")
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


