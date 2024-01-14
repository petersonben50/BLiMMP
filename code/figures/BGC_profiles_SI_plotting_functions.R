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
axis_tick_label_size <- 6/12
axis_label_size <- 7/12

axis_line_1 <- 0
axis_line_2 <- 1.25
axis_label_line_1 <- 0.25
axis_label_line_2 <- 1.5


#### Define function: Plot exo data from this study ####
plot.exo.data <- function(exo.data.to.use = exo.data,
                          date.of.sampling = "2021-09-10",
                          plot_legend_only = "NO") {
  par(mar = c(0, 0, 0, 0),
      tck = -0.008)
  if (plot_legend_only == "NO") {
    #### Exo variables ####
    DO_fudge_factor <- 15
    turb_fudge_factor <- 0.5
    
    #### Isolate needed exo data ####
    exo.data.for.profile <- exo.data.to.use %>%
      filter(date == date.of.sampling)
    
    
    #### Plot the temp points and set up graph ####
    plot(x = NULL,
         y = NULL,
         xlim = xlim_setting,
         xaxs = "i",
         yaxs = "i",
         xaxt = "n",
         yaxt = "n",
         xlab = "",
         ylab = "",
         ylim = depths_to_use)
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
    
    #### Add scales ####
    par(mgp = c(1.5, 0.15, 0))
    axis(2,
         at = seq(0, 25, by = 5),
         labels = seq(0, 25, by = 5),
         cex.axis = axis_tick_label_size,
         las = 1)
    par(mgp = c(1.5, -0.3, 0))
    axis(1,
         line = axis_line_1,
         at = seq(0, 10, by = 5),
         labels = seq(0, 10, by = 5)*DO_fudge_factor,
         cex.axis = axis_tick_label_size,
         gap.axis = -1,
         padj = 0)
    axis(1,
         line = axis_line_2,
         at = seq(0, 12, by = 4),
         labels = seq(0, 12, by = 4)*turb_fudge_factor,
         cex.axis = axis_tick_label_size,
         gap.axis = -1,
         padj = 0)
    
    #### Add axis labels ####
    mtext("Depth (m)",
          side = 2,
          line = 0.75,
          cex = axis_label_size)
    mtext(paste("DO (%)",
                sep = ""),
          side = 1,
          line = axis_label_line_1,
          cex = axis_label_size)
    mtext(paste("Turbidity (FNU)",
                sep = ""),
          side = 1,
          line = axis_label_line_2,
          cex = axis_label_size)
    
    #### Line for SWI ####
    abline(h = swi_line,
           lty = 3)
  } else if (plot_legend_only == "YES") {
    #### Legend ####
    plot.new()
    legend("bottom",
           legend = c("DO", "Turb."),
           col = c(cb.translator["black"], cb.translator["reddishpurple"]),
           cex = 7/12,
           lwd = c(1.5, 1.5),
           lty = c(2, 3),
           pt.cex = c(NA, NA),
           pch = c(NA, NA),
           bty = "n")
  }
}



#### Define function: Plot metals data ####
plot_metal_data <- function(date.of.sampling,
                            legend_location = NULL,
                            plot_legend_only = "NO") {
  par(mar = c(0, 0, 0, 0),
      tck = -0.010)
  
  #### Set up vectors for redox data ####
  color.vector <- c(cb.translator["orange"], cb.translator["orange"],
                    cb.translator["yellow"], cb.translator["yellow"])
  names(color.vector) <- c("Mn_part_ppm", "Mn_diss_ppm",
                           "Fe_part_ppm", "Fe_diss_ppm")
  point.vector <- c(1, 16, 1, 16)
  names(point.vector) <- names(color.vector)
  line.vector <- c(2, 1, 2, 1)
  names(line.vector) <- names(color.vector)
  naming.vector <- c("Part. Mn (ppm)", "Diss. Mn (ppm)",
                     "Part. Fe (ppm)", "Diss. Fe (ppm)")
  names(line.vector) <- names(color.vector)
  
  fudge_vector <- c(rep(25, 4))
  names(fudge_vector) <- names(color.vector)
  
  if (plot_legend_only == "NO") {
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
    par(mgp = c(1.5, -0.3, 0))
    axis(1,
         at = seq(0, 12, by = 5),
         labels = seq(0, 12, by = 5)/fudge_vector["Mn_part_ppm"],
         cex.axis = axis_tick_label_size,
         gap.axis = -1,
         padj = 0)
  
    #### Add axis labels ####
    mtext(paste("Metals (ppm)",
                sep = ""),
          side = 1,
          line = axis_label_line_1,
          cex = axis_label_size)
  
    #### Add SWI line ####
    abline(h = swi_line,
           lty = 3)
  } else if (plot_legend_only == "YES") {
    #### Legend ####
    plot.new()
    legend("bottom",
           legend = naming.vector,
           col = color.vector,
           cex = 6/12,
           lwd = 1.5,
           lty = line.vector,
           pt.cex = 0.8,
           pch = point.vector,
           bty = "n")
    }
  }


#### Define function: Plot sulfur data ####
plot_sulfur_data <- function(date.of.sampling,
                             plot_legend_only = "NO") {
  par(mar = c(0, 0, 0, 0),
      tck = -0.008)
  #### Set up vectors for redox data ####
  color.vector <- c(cb.translator["blue"], cb.translator["blue"])
  names(color.vector) <- c("sulfate_uM", "sulfide_uM")
  point.vector <- c(1, 16)
  names(point.vector) <- names(color.vector)
  line.vector <- c(2, 1)
  names(line.vector) <- names(color.vector)
  naming.vector <- c("Sulfate (µM)", "Sulfide (µM)")
  names(line.vector) <- names(color.vector)
  
  fudge_vector <- c(rep(0.04, 2))
  names(fudge_vector) <- names(color.vector)
  
  if (plot_legend_only == "NO") {
      
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
    par(mgp = c(1.5, -0.3, 0))
    axis(1,
         line = axis_line_1,
         at = seq(0, 9, by = 3),
         labels = seq(0, 9, by = 3)/fudge_vector["sulfate_uM"],
         cex.axis = axis_tick_label_size,
         gap.axis = -1,
         padj = 0)
    
    #### Add axis labels ####
    mtext(paste("S species (µM)",
                sep = ""),
          side = 1,
          line = axis_label_line_1,
          cex = axis_label_size)
    
    #### Add SWI line ####
    abline(h = swi_line,
           lty = 3)
  } else if (plot_legend_only == "YES") {
    plot.new()
    #### Legend ####
    legend("bottom",
           legend = naming.vector,
           col = color.vector,
           cex = 6/12,
           lwd = 1.5,
           lty = line.vector,
           pt.cex = 0.8,
           pch = point.vector,
           bty = "n")
  }
}


#### Function to plot filter-passing Hg data ####
plot_filt_Hg_profile <- function(date.of.sampling = "2021-09-10",
                                 plot_legend_only = "NO") {
  
  par(mar = c(0, 0, 0, 0),
      tck = -0.008)
  #### Set up vectors for redox data ####
  color.vector <- c(cb.translator["bluishgreen"], cb.translator["vermillion"], 
                    cb.translator["black"])
  names(color.vector) <- c("FTHG_NG.L", "FMHG_NG.L", "perFMHG")
  point.vector <- c(1, 16, 19)
  names(point.vector) <- names(color.vector)
  line.vector <- c(2, 1, 4)
  names(line.vector) <- names(color.vector)
  naming.vector <- c("FP iHg", "FP MeHg", "% MeHg")
  names(line.vector) <- names(color.vector)
  
  fudge_vector <- c(rep(5, 2), rep(0.1, 2))
  names(fudge_vector) <- names(color.vector)
  
  if (plot_legend_only == "NO") {
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
    par(mgp = c(1.5, -0.3, 0))
    axis(1,
         line = axis_line_1,
         at = seq(0, 9, by = 3),
         labels = seq(0, 9, by = 3)/fudge_vector["FTHG_NG.L"],
         cex.axis = axis_tick_label_size,
         gap.axis = -1,
         padj = 0)
    axis(1,
         line = axis_line_2,
         at = seq(0, 10, by = 2.4),
         labels = seq(0, 10, by = 2.4)/fudge_vector["perFMHG"],
         cex.axis = axis_tick_label_size,
         gap.axis = -1,
         padj = 0)
    
    #### Add axis labels ####
    mtext(paste("FP Hg (ng/L)",
                sep = ""),
          side = 1,
          line = axis_label_line_1,
          cex = axis_label_size)
    mtext(paste("FP MeHg (%)",
                sep = ""),
          side = 1,
          line = axis_label_line_2,
          cex = axis_label_size)
    
    #### Add SWI line ####
    abline(h = swi_line,
           lty = 3)
  } else if (plot_legend_only == "YES") {
    #### Legend ####
    plot.new()
    legend("bottom",
           legend = naming.vector,
           col = color.vector,
           cex = 6/12,
           lwd = 1.5,
           lty = line.vector,
           pt.cex = 0.8,
           pch = point.vector,
           bty = "n")
  }
}


#### Function to plot particulate Hg data ####
plot_part_Hg_profile <- function(date.of.sampling = "2021-09-10",
                                 plot_legend_only = "NO") {
  
  par(mar = c(0, 0, 0, 0),
      tck = -0.008)
  #### Set up vectors for redox data ####
  color.vector <- c(cb.translator["bluishgreen"], cb.translator["vermillion"], 
                    cb.translator["black"])
  names(color.vector) <- c("PTHG_NG.L", "PMHG_NG.L", "perPMHG")
  point.vector <- c(1, 16, 19)
  names(point.vector) <- names(color.vector)
  line.vector <- c(2, 1, 4)
  names(line.vector) <- names(color.vector)
  naming.vector <- c("Part. iHg", "Part. MeHg", "% Part. MeHg")
  names(line.vector) <- names(color.vector)
  
  fudge_vector <- c(rep(5, 2), rep(0.1, 2))
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
  
  if (plot_legend_only == "NO") {
      
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
    par(mgp = c(1.5, -0.3, 0))
    axis(1,
         line = axis_line_1,
         at = seq(0, 9, by = 3),
         labels = seq(0, 9, by = 3)/fudge_vector["PTHG_NG.L"],
         cex.axis = axis_tick_label_size,
         gap.axis = -1,
         padj = 0)
    axis(1,
         line = axis_line_2,
         at = seq(0, 10, by = 2.4),
         labels = seq(0, 10, by = 2.4)/fudge_vector["perPMHG"],
         cex.axis = axis_tick_label_size,
         gap.axis = -1,
         padj = 0)
    
    #### Add axis labels ####
    mtext(paste("Part. Hg (ng/L)",
                sep = ""),
          side = 1,
          line = axis_label_line_1,
          cex = axis_label_size)
    mtext(paste("Part. MeHg (%)",
                sep = ""),
          side = 1,
          line = axis_label_line_2,
          cex = axis_label_size)
    
    #### Add SWI line ####
    abline(h = swi_line,
           lty = 3)
  } else if (plot_legend_only == "YES") {
    plot.new()
    #### Legend ####
    legend("bottom",
           legend = naming.vector,
           col = color.vector,
           cex = 6/12,
           lwd = 1.5,
           lty = line.vector,
           pt.cex = 0.8,
           pch = point.vector,
           bty = "n")
  }
}



#### Function to plot particulate partitioning Hg data ####
plot_part_part_Hg_profile <- function(date.of.sampling = "2021-09-10",
                                      plot_legend_only = "NO") {
  
  par(mar = c(0, 0, 0, 0),
      tck = -0.008)
  #### Set up vectors for redox data ####
  color.vector <- c(cb.translator["bluishgreen"], cb.translator["vermillion"])
  names(color.vector) <- c("PTHG_NG.G", "PMHG_NG.G")
  point.vector <- c(1, 16)
  names(point.vector) <- names(color.vector)
  line.vector <- c(2, 1)
  names(line.vector) <- names(color.vector)
  naming.vector <- c("Part. iHg", "Part. MeHg")
  names(line.vector) <- names(color.vector)
  fudge_vector <- 0.04
  
  #### Set up data ####
  waterChem.data.filtered <- waterChem.data %>%
    filter(date == date.of.sampling)
  waterChem.data.filtered <- waterChem.data.filtered[, c("depth", names(color.vector))]
  waterChem.data.filtered <- waterChem.data.filtered %>%
    gather(key = constituent,
           value = concentration,
           -1) %>%
    filter(!is.na(concentration))
  
  if (plot_legend_only == "NO") {
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
      points(x = point.data$concentration * fudge_vector,
             y = point.data$depth,
             col = color.vector[species],
             pch = point.vector[species])
      line.data <- point.data %>%
        group_by(depth) %>%
        summarise(concentration = mean(concentration)) %>%
        arrange(depth)
      lines(x = line.data$concentration * fudge_vector,
            y = line.data$depth,
            col = color.vector[species],
            lty = line.vector[species])
    }
    
    #### Add scales ####
    par(mgp = c(1.5, -0.3, 0))
    axis(1,
         line = axis_line_1,
         at = seq(0, 9, by = 3),
         labels = seq(0, 9, by = 3)/fudge_vector,
         cex.axis = axis_tick_label_size,
         gap.axis = -1,
         padj = 0)
    
    #### Add axis labels ####
    mtext(paste("Part. Hg (ng/g)",
                sep = ""),
          side = 1,
          line = axis_label_line_1,
          cex = axis_label_size)
    
    #### Add SWI line ####
    abline(h = swi_line,
           lty = 3)
  } else if (plot_legend_only == "YES") {
    plot.new()
    legend("bottom",
           legend = naming.vector,
           col = color.vector,
           cex = 6/12,
           lwd = 1.5,
           lty = line.vector,
           pt.cex = 0.8,
           pch = point.vector,
           bty = "n")
  }
}



#### Function to plot DOC and SPM data ####
plot_C_profile <- function(date.of.sampling = "2021-09-10",
                           plot_legend_only = "NO") {
  
  par(mar = c(0, 0, 0, 0),
      tck = -0.008)
  #### Set up vectors for redox data ####
  color.vector <- c(cb.translator["skyblue"], cb.translator["reddishpurple"])
  names(color.vector) <- c("DOC_MG.L", "SPM_MG.L")
  point.vector <- c(1, 16)
  names(point.vector) <- names(color.vector)
  line.vector <- c(2, 1)
  names(line.vector) <- names(color.vector)
  naming.vector <- c("DOC (mg/L)", "SPM (mg/L)")
  names(line.vector) <- names(color.vector)
  fudge_vector <- 2
  
  if (plot_legend_only == "NO") {
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
      points(x = point.data$concentration * fudge_vector,
             y = point.data$depth,
             col = color.vector[species],
             pch = point.vector[species])
      line.data <- point.data %>%
        group_by(depth) %>%
        summarise(concentration = mean(concentration)) %>%
        arrange(depth)
      lines(x = line.data$concentration * fudge_vector,
            y = line.data$depth,
            col = color.vector[species],
            lty = line.vector[species])
    }
    
    #### Add scales ####
    par(mgp = c(1.5, -0.3, 0))
    axis(1,
         line = axis_line_1,
         at = seq(0, 9, by = 3),
         labels = seq(0, 9, by = 3)/fudge_vector,
         cex.axis = axis_tick_label_size,
         gap.axis = -1,
         padj = 0)
    
    #### Add axis labels ####
    mtext(paste("DOC/SPM (mg/L)",
                sep = ""),
          side = 1,
          line = axis_label_line_1,
          cex = axis_label_size)
    
    #### Add SWI line ####
    abline(h = swi_line,
           lty = 3)
  } else if (plot_legend_only == "YES") {
    # plot.new()
    legend("bottom",
           legend = naming.vector,
           col = color.vector,
           cex = 6/12,
           lwd = 1.5,
           lty = line.vector,
           pt.cex = 0.8,
           pch = point.vector,
           bty = "n")
  }
}


#### Function to plot DGM data ####
plot_DGM_profile <- function(date.of.sampling = "2021-09-10",
                             plot_legend_only = "NO") {
  
  par(mar = c(0, 0, 0, 0),
      tck = -0.008)
  #### Set up vectors for redox data ####
  color.vector <- c(cb.translator["orange"])
  names(color.vector) <- c("DGM_NG.L")
  point.vector <- c(16)
  names(point.vector) <- names(color.vector)
  line.vector <- c(1)
  names(line.vector) <- names(color.vector)
  naming.vector <- c("DGM (ng/L)")
  names(line.vector) <- names(color.vector)
  fudge_vector <- 300
  
  if (plot_legend_only == "NO") {
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
      points(x = point.data$concentration * fudge_vector,
             y = point.data$depth,
             col = color.vector[species],
             pch = point.vector[species])
      line.data <- point.data %>%
        group_by(depth) %>%
        summarise(concentration = mean(concentration)) %>%
        arrange(depth)
      lines(x = line.data$concentration * fudge_vector,
            y = line.data$depth,
            col = color.vector[species],
            lty = line.vector[species])
    }
    
    #### Add scales ####
    par(mgp = c(1.5, -0.3, 0))
    axis(1,
         line = axis_line_1,
         at = seq(0, 12, by = 3),
         labels = seq(0, 12, by = 3)/fudge_vector*1000,
         cex.axis = axis_tick_label_size,
         gap.axis = -1,
         padj = 0)
    
    #### Add axis labels ####
    mtext(paste("DGM (pg/L)",
                sep = ""),
          side = 1,
          line = axis_label_line_1,
          cex = axis_label_size)
    
    #### Add SWI line ####
    abline(h = swi_line,
           lty = 3)
  } else if (plot_legend_only == "YES") {
    plot.new()
    legend("bottom",
           legend = naming.vector,
           col = color.vector,
           cex = 6/12,
           lwd = 1.5,
           lty = line.vector,
           pt.cex = 0.8,
           pch = point.vector,
           bty = "n")
  }
}
