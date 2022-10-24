#### code/BGC_profiles/DGM_profiles.R ####
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("~/Documents/research/BLiMMP/")
library(readxl)
library(tidyverse)
cb.translator <- readRDS("references/colorblind_friendly_colors.rds")


#### Read in data ####
DGM.data <- readRDS("dataEdited/Hg/DGM_2021_data.rds")


#### Split up by date ####
DGM.data.Sept <- DGM.data %>%
  filter(startDate == "2021-09-10")
DGM.data.Oct <- DGM.data %>%
  filter(startDate == "2021-10-14")


#### Set up data for lines ####
DGM.line.data.Sept <- DGM.data.Sept %>%
  group_by(depth) %>%
  summarise(concentration_ng.L = mean(concentration_ng.L))
DGM.line.data.Oct <- DGM.data.Oct %>%
  group_by(depth) %>%
  summarise(concentration_ng.L = mean(concentration_ng.L))



#### Generate plot ####
pdf("results/BGC_profiles/DGM_profiles.pdf",
    height = 4.5,
    width = 6)
par(mfrow = c(1, 2),
    mar = c(3, 3, 2, 1),
    mgp = c(1.5, 0.4, 0),
    tck = -0.008)
plot(x = DGM.data.Sept$concentration_ng.L*1000,
     y = DGM.data.Sept$depth,
     xlim = c(0, 50),
     ylim = c(24, 0),
     xlab = "DGM (pg/L)",
     ylab = "Depth (m)",
     pch = 16,
     col = "orange",
     cex = 1.75)
lines(x = DGM.line.data.Sept$concentration_ng.L*1000,
      y = DGM.line.data.Sept$depth,
      col = "orange", lwd = 3.5)
title(main = "2021-09-10",
      line = 0.5)

plot(x = DGM.data.Oct$concentration_ng.L*1000,
     y = DGM.data.Oct$depth,
     xlim = c(0, 50),
     ylim = c(24, 0),
     xlab = "DGM (pg/L)",
     ylab = "Depth (m)",
     pch = 16,
     col = "orange",
     cex = 1.75)
lines(x = DGM.line.data.Oct$concentration_ng.L*1000,
      y = DGM.line.data.Oct$depth,
      col = "orange", lwd = 3.5)
title(main = "2021-10-14",
      line = 0.5)
dev.off()
