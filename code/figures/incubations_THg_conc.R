#### code/figures/incubations_MeHg_conc.R ####
# Written for BLiMMP project
# Benjamin D. Peterson



#### Prep workspace ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP")
library(ggpubr)
library(lubridate)
library(readxl)
library(tidyverse)
cb.translator <- readRDS("references/colorblind_friendly_colors.rds")


#### Treatment renaming vector ####
color.vector <- c(cb.translator["vermillion"],
                  cb.translator["bluishgreen"],
                  cb.translator["black"])
names(color.vector) <- c("filtered-unamended",
                         "unfiltered-unamended",
                         "unfiltered-molybdate")
point.vector.to.use <- c(16, 4)
names(point.vector.to.use) <- c(TRUE, FALSE)


#### Read in data ####
Hg_inc_data <- read.csv("dataFinal/incubation_Hg_conc_data.csv") %>%
  filter(treatment %in% names(color.vector))
# Hg_inc_data[which(Hg_inc_data$MeHg_ambient_ppt < 0), "MeHg_ambient_ppt"] <- min(abs(Hg_inc_data$MeHg_ambient_ppt))

#### Function: Plot of Me198Hg given date and depth ####
plot.Hg.species.over.incubation <- function(date_to_use,
                                            depth_to_use,
                                            constituent = "HgT_ambient_ppt",
                                            conc_range,
                                            day.range,
                                            xlab.to.use = NULL,
                                            ylab.to.use = NULL) {
  plotting.data <- Hg_inc_data %>%
    filter(startDate == date_to_use,
           depth == depth_to_use)
  plotting.data[, "plotting_data_column"] <- plotting.data[, constituent]
  
  plot.to.show <- plotting.data %>%
    ggplot(aes(x = durationInDays,
               y = plotting_data_column,
               color = treatment)) +
    geom_point() +
    geom_line(aes(group = incubationID)) +
    scale_color_manual(values = color.vector) +
    scale_y_continuous(limits = conc_range,
                       transform = "log10") +
    xlim(day.range) +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black"),
          axis.title.x = element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.text.x =  element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          axis.ticks.x = element_line(colour = "black"),
          axis.ticks.y = element_line(colour = "black"),
          axis.line = element_line(colour = "black")) +
    labs(subtitle = paste(date_to_use, "-", depth_to_use, "m"))
  
  if (is.null(xlab.to.use)) {
    plot.to.show <- plot.to.show +
      theme(axis.title.x = element_blank())
  } else {
    plot.to.show <- plot.to.show +
      xlab(xlab.to.use)
  }
  
  if (is.null(ylab.to.use)) {
    plot.to.show <- plot.to.show +
      theme(axis.title.y = element_blank())
  } else {
    plot.to.show <- plot.to.show +
      ylab(ylab.to.use)
  }
  
  plot.to.show
}



#### Generate plots for each date/depth ####
concentration.to.use <- c(0.35, 10)
constituent.ylab <- "THg (ng/L)"
day.range.to.use <- c(0, 4)

sept.2020.11m <- plot.Hg.species.over.incubation("2020-09-02", 11, conc_range = concentration.to.use, ylab.to.use = constituent.ylab, day.range = day.range.to.use)
sept.2020.15m <- plot.Hg.species.over.incubation("2020-09-02", 15.5, conc_range = concentration.to.use, day.range = day.range.to.use)
sept.2020.20m <- plot.Hg.species.over.incubation("2020-09-02", 20.7, conc_range = concentration.to.use, day.range = day.range.to.use)

oct.2020.15m <- plot.Hg.species.over.incubation("2020-10-10", 15.7, conc_range = concentration.to.use, ylab.to.use = constituent.ylab, day.range = day.range.to.use)
empty.plot <- ggplot() + theme_void()
oct.2020.20m <- plot.Hg.species.over.incubation("2020-10-10", 20.9, conc_range = concentration.to.use, day.range = day.range.to.use)

sept.2021.10m <- plot.Hg.species.over.incubation("2021-09-10", 10.8, conc_range = concentration.to.use, ylab.to.use = constituent.ylab, day.range = day.range.to.use)
sept.2021.11m <- plot.Hg.species.over.incubation("2021-09-10", 11.9, conc_range = concentration.to.use, day.range = day.range.to.use)
sept.2021.19m <- plot.Hg.species.over.incubation("2021-09-10", 19.9, conc_range = concentration.to.use, day.range = day.range.to.use)

oct.2021.14m <- plot.Hg.species.over.incubation("2021-10-14", 14.2, conc_range = concentration.to.use, xlab = "Days", ylab.to.use = constituent.ylab, day.range = day.range.to.use)
oct.2021.15m <- plot.Hg.species.over.incubation("2021-10-14", 15.2, conc_range = concentration.to.use, xlab = "Days", day.range = day.range.to.use)
oct.2021.19m <- plot.Hg.species.over.incubation("2021-10-14", 19.9, conc_range = concentration.to.use, xlab = "Days", day.range = day.range.to.use)


#### Write out plots to PDF ####
pdf("results/figures/incubations_THg_conc.pdf",
    width = 7,
    height = 7)
ggarrange(sept.2020.11m, sept.2020.15m, sept.2020.20m,
          oct.2020.15m, empty.plot, oct.2020.20m,
          sept.2021.10m, sept.2021.11m, sept.2021.19m,
          oct.2021.14m, oct.2021.15m, oct.2021.19m,
          nrow = 4, ncol = 3)
dev.off()
