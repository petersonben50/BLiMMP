#### code/incubations/HgT_in_incubations.R ####
# Written for BLiMMP project
# Benjamin D. Peterson


#### Clean up ####
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
Hg.incubation.data <- read.csv("dataEdited/incubation_Hg_data.csv") %>%
  filter(treatment %in% names(color.vector),
         year(startDate) %in% c(2020, 2021)) %>%
  mutate(treatment = fct_relevel(treatment,
                                 names(color.vector)))



#### Function: Plot of Hg species given date and depth ####
plot.Hg.species.over.incubation <- function(date.to.use,
                                            depth.to.use,
                                            constituent,
                                            conc.range,
                                            day.range,
                                            xlab.to.use = NULL,
                                            ylab.to.use = NULL) {
  plotting.data <- Hg.incubation.data %>%
    filter(startDate == date.to.use,
           depth == depth.to.use)
  plotting.data[, "plotting_data_column"] <- plotting.data[, constituent]
  
  plot.to.show <- plotting.data %>%
    ggplot(aes(x = durationInDays,
               y = plotting_data_column,
               color = treatment)) +
    geom_point() +
    geom_line(aes(group = incubationID)) +
    scale_color_manual(values = color.vector) +
    ylim(conc.range) +
    xlim(day.range) +
    theme_classic() +
    theme(legend.position = "none") +
    labs(subtitle = paste(date.to.use, "-", depth.to.use, "m"))
  
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



#### Function: All plots for a given constituent ####
make.all.plots.for.constituent <- function(constituent.to.use,
                                           concentration.to.use,
                                           constituent.ylab = NULL,
                                           day.range.to.use = c(0, 3.6)) {

  sept.2020.11m <- plot.Hg.species.over.incubation("2020-09-02", 11, constituent = constituent.to.use, conc.range = concentration.to.use, ylab.to.use = constituent.ylab, day.range = day.range.to.use)
  sept.2020.15m <- plot.Hg.species.over.incubation("2020-09-02", 15.5, constituent = constituent.to.use, conc.range = concentration.to.use, day.range = day.range.to.use)
  sept.2020.20m <- plot.Hg.species.over.incubation("2020-09-02", 20.7, constituent = constituent.to.use, conc.range = concentration.to.use, day.range = day.range.to.use)
  
  oct.2020.15m <- plot.Hg.species.over.incubation("2020-10-10", 15.7, constituent = constituent.to.use, conc.range = concentration.to.use, ylab.to.use = constituent.ylab, day.range = day.range.to.use)
  empty.plot <- ggplot() + theme_void()
  oct.2020.20m <- plot.Hg.species.over.incubation("2020-10-10", 20.9, constituent = constituent.to.use, conc.range = concentration.to.use, day.range = day.range.to.use)

  sept.2021.10m <- plot.Hg.species.over.incubation("2021-09-10", 10.8, constituent = constituent.to.use, conc.range = concentration.to.use, ylab.to.use = constituent.ylab, day.range = day.range.to.use)
  sept.2021.11m <- plot.Hg.species.over.incubation("2021-09-10", 11.9, constituent = constituent.to.use, conc.range = concentration.to.use, day.range = day.range.to.use)
  sept.2021.19m <- plot.Hg.species.over.incubation("2021-09-10", 19.9, constituent = constituent.to.use, conc.range = concentration.to.use, day.range = day.range.to.use)

  oct.2021.14m <- plot.Hg.species.over.incubation("2021-10-14", 14.2, constituent = constituent.to.use, conc.range = concentration.to.use, xlab = "Days", ylab.to.use = constituent.ylab, day.range = day.range.to.use)
  oct.2021.15m <- plot.Hg.species.over.incubation("2021-10-14", 15.2, constituent = constituent.to.use, conc.range = concentration.to.use, xlab = "Days", day.range = day.range.to.use)
  oct.2021.19m <- plot.Hg.species.over.incubation("2021-10-14", 19.9, constituent = constituent.to.use, conc.range = concentration.to.use, xlab = "Days", day.range = day.range.to.use)
  
  ggarrange(sept.2020.11m, sept.2020.15m, sept.2020.20m,
            oct.2020.15m, empty.plot, oct.2020.20m,
            sept.2021.10m, sept.2021.11m, sept.2021.19m,
            oct.2021.14m, oct.2021.15m, oct.2021.19m,
            nrow = 4, ncol = 3)  
}




#### Ambient Hg species plots ####
make.all.plots.for.constituent("HgT_ambient_ppt", concentration.to.use = c(0, 7), constituent.ylab = "HgT (ppt)")
make.all.plots.for.constituent("MeHg_ambient_ppt", concentration.to.use = c(0, 4.5), constituent.ylab = "MeHg (ppt)")



#### Save out PDF with HgT, MeHg, and %MeHg ####
pdf("results/incubations/timecourse_data.pdf",
    width = 7.2,
    height = 10)

make.all.plots.for.constituent("HgT_198_ppt", concentration.to.use = c(0, 1.2), constituent.ylab = "198HgT (ppt)")
make.all.plots.for.constituent("HgT_204_ppt", concentration.to.use = c(0, 1.4), constituent.ylab = "204HgT (ppt)")

make.all.plots.for.constituent("MeHg_198_ppt", concentration.to.use = c(-0.1, 0.6), constituent.ylab = "Me198Hg (ppt)")
make.all.plots.for.constituent("MeHg_204_ppt", concentration.to.use = c(0, 1.6), constituent.ylab = "Me204Hg (ppt)")

make.all.plots.for.constituent("percent_198_MeHg", concentration.to.use = c(-4, 65), constituent.ylab = "% Me198Hg")
make.all.plots.for.constituent("percent_204_MeHg", concentration.to.use = c(0, 160), constituent.ylab = "% Me204Hg")

dev.off()




#### Plot: Percent MeHg with just t1 ####
make.all.plots.for.constituent("percent_204_MeHg", concentration.to.use = c(0, 150), constituent.ylab = "% Me204Hg", day.range.to.use = c(0, 1.2))
make.all.plots.for.constituent("percent_198_MeHg", concentration.to.use = c(-4, 25), constituent.ylab = "% Me198Hg", day.range.to.use = c(0, 1.2))
