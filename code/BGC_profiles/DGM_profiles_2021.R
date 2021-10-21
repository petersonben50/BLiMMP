#### code/incubations/profiles_20201010.R ####
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("~/Documents/research/BLiMMP/")
library(readxl)
library(tidyverse)
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")
source("code/BLiMMP_functions.R")


#### Read in data ####
DGM.data <- readRDS("dataEdited/Hg/DGM_2021_data.rds")
DGM.data.sept <- DGM.data %>%
  filter(month(startDate) == 9)
DGM.data.oct <- DGM.data %>%
  filter(month(startDate) == 10)




pdf("~/Downloads/DGM_temp_DO_turb.pdf",
    height = 8,
    width = 8)
par(mfrow = c(1, 2))
plot.exo.data("2021-09-10",
              main.title = "",
              legend.location = "topleft")
points(x = DGM.data.sept$concentration_ng.L*1000 / 5,
       y = -DGM.data.sept$depth,
       pch = 16,
       col = "orange", cex = 1.75)

DGM.line <- DGM.data.sept %>%
  group_by(depth) %>%
  summarise(concentration_ng.L = mean(concentration_ng.L))

lines(x = DGM.line$concentration_ng.L*1000 / 5,
      y = -DGM.line$depth,
      col = "orange", lwd = 3.5)
axis(3,
     at = seq(0, 10, by = 2),
     labels = seq(0, 10, by = 2)*5,
     par(mgp=c(1.4, 0.15, 0), tck=-0.008),)
title(main = "DGM (pg/L)",
      col.main = "orange",
      line = 1.5)
text(x = 7,
     y = -24,
     cex = 1.3,
     labels = "September 10th, 2021")
# dev.off()


plot.exo.data("2021-10-14",
              main.title = "",
              legend.location = "topleft")
points(x = DGM.data.oct$concentration_ng.L*1000 / 5,
       y = -DGM.data.oct$depth,
       pch = 16,
       col = "orange", cex = 1.75)

DGM.line <- DGM.data.oct %>%
  group_by(depth) %>%
  summarise(concentration_ng.L = mean(concentration_ng.L))

lines(x = DGM.line$concentration_ng.L*1000 / 5,
      y = -DGM.line$depth,
      col = "orange", lwd = 3.5)
axis(3,
     at = seq(0, 10, by = 2),
     labels = seq(0, 10, by = 2)*5,
     par(mgp=c(1.4, 0.15, 0), tck=-0.008),)
title(main = "DGM (pg/L)",
      col.main = "orange",
      line = 1.5)
text(x = 7,
     y = -24,
     cex = 1.3,
     labels = "October 14th, 2021")

dev.off()