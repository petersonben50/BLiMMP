#### code/figures/BGC_profiles_main_plotting_functions.R ####
# Benjamin D. Peterson



library(fields)
library(lubridate)
library(readxl)
library(rLakeAnalyzer)
library(stringr)
library(tidyverse)
library(viridis)

cb.translator <- readRDS("references/colorblind_friendly_colors.rds")

#### Global variables ####
detection.vector <- c(4, 16)
depths_to_use <- c(25, 0)
xlim_setting <- c(-0.5, 12)
swi_line <- 24
axis_tick_label_size <- 7/12
axis_label_size <- 8/12


#### Define function: Plot exo data from this study ####
plot.exo.data <- function(exo.data.to.use = exo.data,
                          incubation.depths.to.use = incubation.depths,
                          date.of.sampling = "2021-09-10",
                          plot_legend = NULL) {
  par(mar = c(0, 0, 0, 0),
      tck = -0.008)
  #### Exo variables ####
  DO_fudge_factor <- 12
  turb_fudge_factor <- 0.5
  
  #### Isolate needed exo data ####
  exo.data.for.profile <- exo.data.to.use %>%
    filter(date == date.of.sampling)
  if (!is.null(incubation.depths.to.use)) {
    sampling.depths.of.interest <- incubation.depths.to.use %>%
      filter(startDate == date.of.sampling) %>%
      select(depth) %>%
      unlist(use.names = FALSE)
  }
  
  #### Plot the temp points and set up graph ####
  plot(x = (exo.data.for.profile$Temp_C-5)/2.5,
       y = exo.data.for.profile$depth,
       type = "l",
       lwd = 3,
       xlim = xlim_setting,
       xaxs = "i",
       yaxs = "i",
       xaxt = "n",
       yaxt = "n",
       xlab = "",
       ylab = "",
       ylim = depths_to_use,
       col = cb.translator["blue"])
  # Plot DO
  lines(x = exo.data.for.profile$ODO_sat/DO_fudge_factor,
        y = exo.data.for.profile$depth,
        col = cb.translator["black"],
        lwd = 3,
        lty = 2)
  # Plot the turbidity values
  lines(x = exo.data.for.profile$Turbidity_FNU/turb_fudge_factor,
        y = exo.data.for.profile$depth,
        col = cb.translator["reddishpurple"],
        lwd = 4,
        lty = 3)
  if (!is.null(incubation.depths.to.use)) {
    # Plot the incubation depths
    points(x = rep(9, length(sampling.depths.of.interest)),
           y = sampling.depths.of.interest,
           pch = 21,
           cex = 1,
           bg = cb.translator["yellow"],
           col = cb.translator["black"])
  }

  
  #### Add scales ####
  par(mgp = c(1.5, 0.2, 0))
  axis(2,
       at = seq(0, 25, by = 5),
       labels = seq(0, 25, by = 5),
       cex.axis = axis_tick_label_size,
       las = 1)
  par(mgp = c(1.5, -0.1, 0))
  axis(1,
       at = seq(0, 10, by = 2),
       labels = seq(0, 10, by = 2)*2.5+5,
       cex.axis = axis_tick_label_size,
       gap.axis = -1,
       padj = 0)
  axis(1,
       line = 2,
       at = seq(0, 10, by = 2.5),
       labels = seq(0, 10, by = 2.5)*DO_fudge_factor,
       cex.axis = axis_tick_label_size,
       gap.axis = -1,
       padj = 0)
  axis(1,
       line = 4,
       at = seq(0, 10, by = 2),
       labels = seq(0, 10, by = 2)*turb_fudge_factor,
       cex.axis = axis_tick_label_size,
       gap.axis = -1,
       padj = 0)
  
  #### Add axis labels ####
  mtext("Depth (m)",
        side = 2,
        line = 0.75,
        cex = axis_label_size*1.5)
  mtext(paste("Temp (", "\U02DA","C)",
              sep = ""),
        side = 1,
        line = 0.75,
        cex = axis_label_size)
  mtext(paste("DO (%)",
              sep = ""),
        side = 1,
        line = 2.75,
        cex = axis_label_size)
  mtext(paste("Turbidity (FNU)",
              sep = ""),
        side = 1,
        line = 4.75,
        cex = axis_label_size)

  #### Line for SWI ####
  abline(h = swi_line,
         lty = 3)
  
  #### Legend ####
  if (!is.null(plot_legend)) {
    legend("topleft",
           legend = c("Temp.", "DO", "Turb.", "Inc."),
           col = c(cb.translator["blue"], cb.translator["black"], cb.translator["reddishpurple"], cb.translator["black"]),
           pt.bg = c(NA, NA, NA, cb.translator["yellow"]),
           cex = 7/12,
           lwd = c(1.5, 1.5, 1.5, NA),
           # pt.lwd = c(6, 6, 6, NA),
           lty = c(1, 2, 3, NA),
           pt.cex = c(NA, NA, NA, 1),
           pch = c(NA, NA, NA, 21),
           bty = "n")
  }
}



#### Define function: Plot redox data ####
plot_redox_data <- function(date.of.sampling,
                            legend_location = NULL) {
  par(mar = c(0, 0, 0, 0),
      tck = -0.008)
  #### Set up vectors for redox data ####
  color.vector <- c(cb.translator["orange"], cb.translator["orange"],
                    cb.translator["blue"], cb.translator["blue"])
  names(color.vector) <- c("Mn_part_ppm", "Mn_diss_ppm",
                           "sulfate_ppm", "sulfide_ppm")
  point.vector <- c(1, 16, 1, 16)
  names(point.vector) <- names(color.vector)
  line.vector <- c(2, 1, 2, 1)
  names(line.vector) <- names(color.vector)
  naming.vector <- c("p.Mn (mg/L)", "f.Mn (mg/L)",
                     "Sulfate (mg/L)", "Sulfide (mg/L)")
  names(line.vector) <- names(color.vector)
  
  fudge_vector <- c(rep(20, 2), 0.5, 2)
  names(fudge_vector) <- names(color.vector)
  
  
  #### Set up data ####
  waterChem.data.filtered <- waterChem.data %>%
    filter(date == date.of.sampling)
  waterChem.data.filtered <- waterChem.data.filtered[, c("depth", names(color.vector))]
  waterChem.data.filtered <- waterChem.data.filtered %>%
    gather(key = constituent,
           value = concentration,
           -1) %>%
    filter(!is.na(concentration))
  
  
  #### Generate empty plot ####
  plot(NULL,
       ylim = depths_to_use,
       xlim = xlim_setting,
       xaxs = "i",
       yaxs = "i",
       xaxt = "n",
       yaxt = "n",
       xlab = "",
       ylab = "")
  
  #### Add points and lines ####
  for (species in names(color.vector)) {
    point.data <- waterChem.data.filtered %>%
      filter(constituent == species)
    points(x = point.data$concentration * fudge_vector[species],
           y = point.data$depth,
           col = color.vector[species],
           pch = point.vector[species])
    line.data <- point.data %>%
      group_by(depth) %>%
      summarise(concentration = mean(concentration)) %>%
      arrange(depth)
    lines(x = line.data$concentration * fudge_vector[species],
          y = line.data$depth,
          col = color.vector[species],
          lty = line.vector[species])
  }
    
  #### Add scales ####
  par(mgp = c(1.5, -0.1, 0))
  axis(1,
       at = seq(0, 9, by = 3),
       labels = seq(0, 9, by = 3)/fudge_vector["Mn_part_ppm"],
       cex.axis = axis_tick_label_size,
       gap.axis = -1,
       padj = 0)
  axis(1,
       line = 2,
       at = seq(0, 10, by = 2.4),
       labels = seq(0, 10, by = 2.4)/fudge_vector["sulfide_ppm"],
       cex.axis = axis_tick_label_size,
       gap.axis = -1,
       padj = 0)
  axis(1,
       line = 4,
       at = seq(0, 12, by = 3),
       labels = seq(0, 12, by = 3)/fudge_vector["sulfate_ppm"],
       cex.axis = axis_tick_label_size,
       gap.axis = -1,
       padj = 0)
  
  #### Add axis labels ####
  mtext(paste("Manganese (mg/L)",
              sep = ""),
        side = 1,
        line = 0.75,
        cex = axis_label_size)
  mtext(paste("Sulfide (mg/L)",
              sep = ""),
        side = 1,
        line = 2.75,
        cex = axis_label_size)
  mtext(paste("Sulfate (mg/L)",
              sep = ""),
        side = 1,
        line = 4.75,
        cex = axis_label_size)
  
  #### Add SWI line ####
  abline(h = swi_line,
         lty = 3)
  
  #### Legend ####
  if (!is.null(legend_location)) {
    legend(legend_location,
           legend = naming.vector,
           col = color.vector,
           cex = 7/12,
           lwd = 1.5,
           lty = line.vector,
           pt.cex = 1,
           pch = point.vector,
           bty = "n")
  }
  }


#### Function to plot Hg data ####
plot.Hg.profile <- function(date.of.sampling = "2021-09-10",
                            legend_location = NULL) {
  
  par(mar = c(0, 0, 0, 0),
      tck = -0.008)
  #### Set up vectors for redox data ####
  color.vector <- c(cb.translator["bluishgreen"], cb.translator["vermillion"], 
                    cb.translator["black"])
  names(color.vector) <- c("FiHg_NG.L", "FMHG_NG.L", "perFMHG")
  point.vector <- c(1, 16, 19)
  names(point.vector) <- names(color.vector)
  line.vector <- c(2, 1, 4)
  names(line.vector) <- names(color.vector)
  naming.vector <- c("f.Hg(II)", "f.MeHg", "% MeHg")
  names(line.vector) <- names(color.vector)
  
  fudge_vector <- c(rep(10, 2), 0.1)
  names(fudge_vector) <- names(color.vector)
  
  
  #### Set up data ####
  waterChem.data.filtered <- waterChem.data %>%
    filter(date == date.of.sampling)
  waterChem.data.filtered <- waterChem.data.filtered[, c("depth", names(color.vector))]
  waterChem.data.filtered <- waterChem.data.filtered %>%
    gather(key = constituent,
           value = concentration,
           -1) %>%
    filter(!is.na(concentration))
  
  
  #### Generate empty plot ####
  plot(NULL,
       ylim = depths_to_use,
       xlim = xlim_setting,
       xaxs = "i",
       yaxs = "i",
       xaxt = "n",
       yaxt = "n",
       xlab = "",
       ylab = "")
  
  #### Add points and lines ####
  for (species in names(color.vector)) {
    point.data <- waterChem.data.filtered %>%
      filter(constituent == species)
    points(x = point.data$concentration * fudge_vector[species],
           y = point.data$depth,
           col = color.vector[species],
           pch = point.vector[species])
    line.data <- point.data %>%
      group_by(depth) %>%
      summarise(concentration = mean(concentration)) %>%
      arrange(depth)
    lines(x = line.data$concentration * fudge_vector[species],
          y = line.data$depth,
          col = color.vector[species],
          lty = line.vector[species])
  }
  
  #### Add scales ####
  par(mgp = c(1.5, -0.1, 0))
  axis(1,
       at = seq(0, 10, by = 2.5),
       labels = seq(0, 10, by = 2.5)/fudge_vector["FiHg_NG.L"],
       cex.axis = axis_tick_label_size,
       gap.axis = -1,
       padj = 0)
  axis(1,
       line = 2,
       at = seq(0, 10, by = 2.4),
       labels = seq(0, 10, by = 2.4)/fudge_vector["perFMHG"],
       cex.axis = axis_tick_label_size,
       gap.axis = -1,
       padj = 0)
  
  #### Add axis labels ####
  mtext(paste("Hg(II), MeHg (ng/L)",
              sep = ""),
        side = 1,
        line = 0.75,
        cex = axis_label_size)
  mtext(paste("Percent MeHg",
              sep = ""),
        side = 1,
        line = 2.75,
        cex = axis_label_size)
  
  #### Add SWI line ####
  abline(h = swi_line,
         lty = 3)
  
  
  #### Legend ####
  if (!is.null(legend_location)) {
    legend(legend_location,
           legend = naming.vector,
           col = color.vector,
           cex = 7/12,
           lwd = 1.5,
           lty = line.vector,
           pt.cex = 1,
           pch = point.vector,
           bty = "n")
  }
}
