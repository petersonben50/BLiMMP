
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
