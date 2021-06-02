#### code/incubations/2020/initial_injection_concentration.R ####
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
# source("code/BLiMMP_functions.R")


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
Hg.data$monthOfIncubations <- "October"
Hg.data$monthOfIncubations[which(Hg.data$tripID == "BLiMMP_trip_010")] = "September"
Hg.data$monthOfIncubations <- relevel(as.factor(Hg.data$monthOfIncubations),
                                      c("September"))


#### Prepare data for plotting ####
Hg.data.plotting <- Hg.data %>%
  filter(t == "t0") %>%
  mutate(sampleInfo = paste(monthOfIncubations,
                            depth, sep = "-"),
         sampleInfo = fct_relevel(sampleInfo,
                                  c("September-11m", "September-15.5m", "September-20.7m",
                                    "October-15.7m", "October-20.9m")),
         treatment = fct_relevel(treatment,
                                 names(color.vector))) %>%
  select(sampleInfo, treatment, excess_HgT_198_ng.L, excess_HgT_204_ng.L)


#### Plot 198 amendment ####
Hg.data.plotting %>%
  ggplot(aes(x = treatment,
             y = excess_HgT_198_ng.L,
             color = treatment)) +
  geom_point() +
  facet_wrap(~sampleInfo) +
  ylim(c(0, 1.3)) +
  ylab("Excess 198HgT (ng/L)") +
  xlab(element_blank()) +
  scale_color_manual(values = color.vector) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())


#### Plot 204 amendment ####
Hg.data.plotting %>%
  ggplot(aes(x = treatment,
             y = excess_HgT_204_ng.L,
             color = treatment)) +
  geom_point() +
  facet_wrap(~sampleInfo) +
  ylim(c(0, 1.3)) +
  ylab("Excess 204HgT (ng/L)") +
  xlab(element_blank()) +
  scale_color_manual(values = color.vector) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())


#### Check mean, SD ####
max(Hg.data.plotting$excess_HgT_198_ng.L)
min(Hg.data.plotting$excess_HgT_198_ng.L)
mean(Hg.data.plotting$excess_HgT_198_ng.L)
sd(Hg.data.plotting$excess_HgT_198_ng.L)

sd(Hg.data.plotting$excess_HgT_198_ng.L) / mean(Hg.data.plotting$excess_HgT_198_ng.L)


#### Check mean, SD for each depth ####
Hg.data.plotting %>%
  group_by(sampleInfo) %>%
  summarise(mean = mean(excess_HgT_198_ng.L),
            sd = sd(excess_HgT_198_ng.L),
            min = min(excess_HgT_198_ng.L),
            max = max(excess_HgT_198_ng.L))
