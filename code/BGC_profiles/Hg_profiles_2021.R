#### code/BGC_profiles/Hg_profiles_2021.R ####
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("~/Documents/research/BLiMMP/")
library(readxl)
library(tidyverse)
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")
source("code/BLiMMP_functions.R")
detection.vector <- c(4, 16)
names(detection.vector) <- c("<", "NONE") 

#### Read in data ####
DGM.data <- readRDS("dataEdited/Hg/DGM_2021_data.rds")
DGM.data.sept <- DGM.data %>%
  filter(month(startDate) == 9)
DGM.data.oct <- DGM.data %>%
  filter(month(startDate) == 10)

Hg.data <- read.csv("dataEdited/Hg/Hg_data_clean.csv")


#### Function to plot Hg data from USGS ####
plot.Hg.data <- function(Hg.data.to.use.1,
                         Hg.data.to.use.2 = NULL,
                         color.vector.to.use,
                         range.to.plot = c(0, 1)) {
  # iHg plots
  plot(x = Hg.data.to.use.1$concentration,
       y = Hg.data.to.use.1$depth,
       xlab = "iHg (ng/L)",
       pch = detection.vector[Hg.data.to.use.1$detection_flag],
       ylim = c(25, 0),
       xlim = range.to.plot,
       col = color.vector.to.use[1])
  # iHg lines
  Hg.data.to.use.1.line <- Hg.data.to.use.1 %>%
    group_by(sampleDate, depth) %>% summarise(concentration = mean(concentration)) %>%
    arrange(depth)
  lines(x = Hg.data.to.use.1.line$concentration,
        y = Hg.data.to.use.1.line$depth,
        col = color.vector.to.use[1],
        lwd = 2)
  
  # MeHg plots
  points(x = Hg.data.to.use.2$concentration,
         y = Hg.data.to.use.2$depth,
         pch = detection.vector[Hg.data.to.use.2$detection_flag],
         xlab = "iHg (ng/L)",
         ylim = c(25, 0),
         col = color.vector.to.use[2])
  # MeHg lines
  Hg.data.to.use.2.line <- Hg.data.to.use.2 %>%
    group_by(sampleDate, depth) %>% summarise(concentration = mean(concentration)) %>%
    arrange(depth)
  lines(x = Hg.data.to.use.2.line$concentration,
        y = Hg.data.to.use.2.line$depth,
        col = color.vector.to.use[2],
        lwd = 2)
  legend("topright",
         legend = names(color.vector.to.use),
         fill = color.vector.to.use)
}




#### Plot Sept 2020 data ####

# Dissolved iHg
Sept.2020.FiHg.data <- Hg.data %>%
  filter(sampleDate %in% c("2020-09-02",
                           "2020-09-03"),
         constituent == "FiHg_NG.L") %>%
  select(sampleDate, depth, sampleTime, constituent, concentration) %>%
  # Use the HgT detection flag for iHg.
  left_join(Hg.data %>%
              filter(sampleDate %in% c("2020-09-02",
                                       "2020-09-03"),
                     constituent == "FTHG_NG.L") %>%
              select(sampleDate, depth, sampleTime, detection_flag))

# Dissolved MeHg
Sept.2020.FMHG.data <- Hg.data %>%
  filter(sampleDate %in% c("2020-09-02",
                           "2020-09-03"),
         constituent == "FMHG_NG.L")

# Particulate iHg
Sept.2020.PiHg.data <- Hg.data %>%
  filter(sampleDate %in% c("2020-09-02",
                           "2020-09-03"),
         constituent == "PiHg_NG.L") %>%
  select(sampleDate, depth, sampleTime, constituent, concentration) %>%
  # Use the HgT detection flag for iHg.
  left_join(Hg.data %>%
              filter(sampleDate %in% c("2020-09-02",
                                       "2020-09-03"),
                     constituent == "PTHG_NG.L") %>%
              select(sampleDate, depth, sampleTime, detection_flag))

# Dissolved MeHg
Sept.2020.PMHG.data <- Hg.data %>%
  filter(sampleDate %in% c("2020-09-02",
                           "2020-09-03"),
         constituent == "PMHG_NG.L")

# Color vector
color.vector.Hg.data <- c(cb.translator["vermillion"], cb.translator["bluishgreen"])
names(color.vector.Hg.data) <- c("FiHg", "FMHg")

pdf("results/BGC_profiles/Hg_profile_2020-09-02.pdf",
    height = 5,
    width = 6)
par(mfrow = c(1, 3))
plot.exo.data("2020-09-02",
              main.title = "",
              legend.location = "topleft")
text(x = 6,
     y = -24,
     cex = 1.3,
     labels = "September 10th, 2021")

plot.Hg.data(Hg.data.to.use.1 = Sept.2020.FiHg.data,
             Hg.data.to.use.2 = Sept.2020.FMHG.data,
             color.vector.to.use = color.vector.Hg.data)
plot.Hg.data(Hg.data.to.use.1 = Sept.2020.PiHg.data,
             Hg.data.to.use.2 = Sept.2020.PMHG.data,
             color.vector.to.use = color.vector.Hg.data,
             )


dev.off()







#### Plot October 2020 data ####

# Dissolved iHg
Oct.2020.FiHg.data <- Hg.data %>%
  filter(sampleDate %in% c("2020-10-10"),
         constituent == "FiHg_NG.L") %>%
  select(sampleDate, depth, sampleTime, constituent, concentration) %>%
  # Use the HgT detection flag for iHg.
  left_join(Hg.data %>%
              filter(sampleDate %in% c("2020-10-10"),
                     constituent == "FTHG_NG.L") %>%
              select(sampleDate, depth, sampleTime, detection_flag))

# Dissolved MeHg
Oct.2020.FMHG.data <- Hg.data %>%
  filter(sampleDate %in% c("2020-10-10"),
         constituent == "FMHG_NG.L")

# Particulate iHg
Oct.2020.PiHg.data <- Hg.data %>%
  filter(sampleDate %in% c("2020-10-10"),
         constituent == "PiHg_NG.L") %>%
  select(sampleDate, depth, sampleTime, constituent, concentration) %>%
  # Use the HgT detection flag for iHg.
  left_join(Hg.data %>%
              filter(sampleDate %in% c("2020-10-10"),
                     constituent == "PTHG_NG.L") %>%
              select(sampleDate, depth, sampleTime, detection_flag))

# Particulate MeHg
Oct.2020.PMHG.data <- Hg.data %>%
  filter(sampleDate %in% c("2020-10-10"),
         constituent == "PMHG_NG.L")

# Color vector
color.vector.Hg.data <- c(cb.translator["vermillion"], cb.translator["bluishgreen"])
names(color.vector.Hg.data) <- c("FiHg", "FMHg")

pdf("results/BGC_profiles/Hg_profile_2020-10-10.pdf",
    height = 5,
    width = 6)
par(mfrow = c(1, 3))
plot.exo.data("2020-10-10",
              main.title = "",
              legend.location = "topleft")
text(x = 6,
     y = -24,
     cex = 1.3,
     labels = "October 10th, 2020")

plot.Hg.data(Hg.data.to.use.1 = Oct.2020.FiHg.data,
             Hg.data.to.use.2 = Oct.2020.FMHG.data,
             color.vector.to.use = color.vector.Hg.data)
plot.Hg.data(Hg.data.to.use.1 = Oct.2020.PiHg.data,
             Hg.data.to.use.2 = Oct.2020.PMHG.data,
             color.vector.to.use = color.vector.Hg.data)


dev.off()








#### Plot Sept 2021 data ####

# Dissolved iHg
Sept.FiHg.data <- Hg.data %>%
  filter(sampleDate %in% c("2021-09-10",
                           "2021-09-11"),
         constituent == "FiHg_NG.L") %>%
  select(sampleDate, depth, sampleTime, constituent, concentration) %>%
  # Use the HgT detection flag for iHg.
  left_join(Hg.data %>%
              filter(sampleDate %in% c("2021-09-10",
                                       "2021-09-11"),
                     constituent == "FTHG_NG.L") %>%
              select(sampleDate, depth, sampleTime, detection_flag))

# Dissolved MeHg
Sept.FMHG.data <- Hg.data %>%
  filter(sampleDate %in% c("2021-09-10",
                           "2021-09-11"),
         constituent == "FMHG_NG.L")

# Color vector
color.vector.Hg.data <- c(cb.translator["vermillion"], cb.translator["bluishgreen"])
names(color.vector.Hg.data) <- c("FiHg", "FMHg")

pdf("results/BGC_profiles/Hg_profile_2021-09-10.pdf",
    height = 5,
    width = 6)
par(mfrow = c(1, 3))
plot.exo.data("2021-09-10",
              main.title = "",
              legend.location = "topleft")
text(x = 6,
     y = -24,
     cex = 1.3,
     labels = "September 10th, 2021")

plot.Hg.data(Hg.data.to.use.1 = Sept.FiHg.data,
             Hg.data.to.use.2 = Sept.FMHG.data,
             color.vector.to.use = color.vector.Hg.data)


# September DGM
plot(x = DGM.data.sept$concentration_ng.L,
     y = DGM.data.sept$depth,
     pch = 16,
     xlab = "DGM (ng/L)",
     ylim = c(25, 0),
     xlim = c(0, 0.035),
     col = "orange",
     cex = 1.75)
DGM.line <- DGM.data.sept %>%
  group_by(depth) %>%
  summarise(concentration_ng.L = mean(concentration_ng.L))
lines(x = DGM.line$concentration_ng.L,
      y = DGM.line$depth,
      col = "orange", lwd = 3.5)

dev.off()





#### October 2021 Hg profiles ####
pdf("results/BGC_profiles/Hg_profile_2021-10-10.pdf",
    height = 5,
    width = 6)
par(mfrow = c(1, 3))
plot.exo.data("2021-10-14",
              main.title = "",
              legend.location = "topleft")
# October DGM
plot(x = DGM.data.oct$concentration_ng.L,
     y = DGM.data.oct$depth,
     pch = 16,
     xlab = "DGM (ng/L)",
     xlim = c(0, 0.035),
     ylim = c(25, 0),
     col = "orange", cex = 1.75)
DGM.line.oct <- DGM.data.oct %>%
  group_by(depth) %>%
  summarise(concentration_ng.L = mean(concentration_ng.L))
lines(x = DGM.line.oct$concentration_ng.L,
      y = DGM.line.oct$depth,
      col = "orange", lwd = 3.5)


dev.off()
