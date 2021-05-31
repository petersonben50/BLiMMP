#### code/incubations/2020/HgT_change_across_incubation.R ####
# Written for BLiMMP project
# Benjamin D. Peterson

# Notes in the Obsidian:
# 2020 incubations - Hg data processing

#### Prep workspace ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP/")
library(dplyr)
library(lubridate)
library(patchwork)
library(tidyverse)
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")


#### Treatment renaming vector ####
color.vector <- c(cb.translator["vermillion"],
                  cb.translator["bluishgreen"],
                  cb.translator["black"])
names(color.vector) <- c("filtered-unamended",
                         "unfiltered-unamended",
                         "unfiltered-molybdate")


#### Read in data ####
Hg.data <- read.csv("dataEdited/incubations/2020incubations_Hg_data.csv") %>%
  mutate(depth = paste(depth, "m", sep = "")) %>%
  filter(treatment %in% names(color.vector))


#### Prepare data for plotting ####
Hg.data$monthOfIncubations <- "October"
Hg.data$monthOfIncubations[which(Hg.data$tripID == "BLiMMP_trip_010")] = "September"
Hg.data.plotting <- Hg.data %>%
  mutate(sampleInfo = paste(monthOfIncubations,
                            depth, sep = "-"),
         sampleInfo = fct_relevel(sampleInfo,
                                  c("September-11m", "September-15.5m", "September-20.7m",
                                    "October-15.7m", "October-20.9m")),
         treatment = fct_relevel(treatment,
                                 names(color.vector))) %>%
  select(sampleInfo, incubationID, treatment, t, durationSinceTimepointInDays, amb_MeHg_ng.L)


#### Plot ambient MeHg over incubation ####
pdf("results/incubations/2020incubations/ambient_Hg/MeHg_levels_over_incubation.pdf",
    width = 8,
    height = 6)
Hg.data.plotting %>%
  ggplot(aes(x = t,
             y = amb_MeHg_ng.L,
             color = treatment)) +
  geom_point() +
  geom_line(aes(group = incubationID)) +
  scale_color_manual(values = color.vector) +
  facet_wrap(~sampleInfo, nrow = 2) +
  theme_classic() +
  ylim(c(0, 2.2))
dev.off()


#### Calculate loss of ambient MeHg over each time point ####
amb.MeHg.delta <- Hg.data.plotting %>%
  select(incubationID, t, amb_MeHg_ng.L) %>%
  spread(key = t,
         value = amb_MeHg_ng.L) %>%
  mutate(t2 = t2 - t1,
         t1 = t1 - t0) %>%
  select(incubationID, t1, t2) %>%
  gather(key = t,
         value = net_amb_MeHg_ng.L,
         -1) %>%
  left_join(Hg.data.plotting) %>%
  mutate(net_amb_MeHg_ng.L.day = round(net_amb_MeHg_ng.L / durationSinceTimepointInDays, 3)) %>%
  select(incubationID, sampleInfo, treatment, t, net_amb_MeHg_ng.L, net_amb_MeHg_ng.L.day) %>%
  filter(!is.na(net_amb_MeHg_ng.L.day))



#### Plot rates of loss of ambient MeHg ####
pdf(file = "results/incubations/2020incubations/ambient_Hg/MeHg_loss_over_incubation.pdf",
    width = 8,
    height = 6)
amb.MeHg.delta %>%
  ggplot(aes(x = t,
             y = net_amb_MeHg_ng.L.day,
             color = treatment)) +
  geom_point() +
  geom_line(aes(group = incubationID)) +
  scale_color_manual(values = color.vector) +
  facet_wrap(~sampleInfo, nrow = 2) +
  theme_classic() +
  ylim(c(-1, 1))
dev.off()



#### Summarize the changes in ambient MeHg ####
amb.MeHg.delta %>%
  group_by(t) %>%
  summarise(net_amb_MeHg_ng.L.day.mean = mean(net_amb_MeHg_ng.L.day),
            net_amb_MeHg_ng.L.day.sd = sd(net_amb_MeHg_ng.L.day))
mean.data <- amb.MeHg.delta %>%
  group_by(t, sampleInfo, treatment) %>%
  summarise(net_amb_MeHg_ng.L.day.mean = mean(net_amb_MeHg_ng.L.day),
            net_amb_MeHg_ng.L.day.sd = sd(net_amb_MeHg_ng.L.day))
mean.data

