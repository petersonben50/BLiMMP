#### code/microbial/leucine_uptake_profiles_2021.R ####
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("~/Documents/research/BLiMMP/")
library(ggpubr)
library(patchwork)
library(tidyverse)
cb.translator <- readRDS("~/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")


#### Color vector ####
color.vector <- cb.translator[c("bluishgreen", "black")]
names(color.vector) <- c("control", "molybdate")


#### Read in leucine uptake data from September ####
leu.uptake.sept.8hr <- read.csv("dataEdited/leucineUptake/20210916_leucine.csv",
                                stringsAsFactors = FALSE)
leu.uptake.sept.24hr <- read.csv("dataEdited/leucineUptake/20210917_leucine.csv",
                                stringsAsFactors = FALSE)




#### Plot out leucine uptake after 8hr ####
leu.uptake.8hr.plot <- leu.uptake.sept.8hr %>%
  filter(treatment != "filtered") %>%
  mutate(depth = as.factor(depth)) %>%
  ggplot(aes(y = µgBCP_per_L_hr,
             x = treatment,
             col = treatment)) +
  geom_boxplot() +
  geom_point() +
  facet_wrap(~depth) +
  scale_color_manual(values = color.vector) +
  ylim(c(0, 25)) +
  ylab("Bacterial carbon production (µgC/L/hr)") +
  theme_bw() +
  theme(legend.position = c(0.8, 0.8))




#### Plot out leucine uptake after 24hr ####
leu.uptake.24hr.plot <- leu.uptake.sept.24hr %>%
  filter(treatment != "filtered") %>%
  mutate(depth = as.factor(paste(depth,
                                 "m"))) %>%
  ggplot(aes(y = µgBCP_per_L_hr,
             x = treatment,
             col = treatment)) +
  geom_boxplot() +
  geom_point() +
  facet_wrap(~depth) +
  scale_color_manual(values = color.vector) +
  ylim(c(0, 20)) +
  ylab("Bacterial carbon production (µgC/L/hr)") +
  theme_bw() +
  theme(legend.position = c(0.8, 0.8))


#### View both plots ####
leu.uptake.8hr.plot / leu.uptake.24hr.plot
# These are very similar, will probably just
# go with the 8 hr time point for publications.
leu.uptake.8hr.plot


#### Plot out leucine uptake in October ####
# Read in leucine uptake data
leu.uptake.oct <- read.csv("dataEdited/leucineUptake/20211018_leucine.csv",
                                 stringsAsFactors = FALSE)

leu.uptake.oct.plot <- leu.uptake.oct %>%
  mutate(depth = as.factor(depth)) %>%
  ggplot(aes(y = µgBCP_per_L_hr,
             x = treatment,
             col = treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(shape = as.character(replicate)),
              width = 0.1) +
  facet_wrap(~depth) +
  scale_color_manual(values = color.vector) +
  ylim(c(0, 25)) +
  ylab("Bacterial carbon production (µgC/L/hr)") +
  theme_bw() +
  theme(legend.position = c(0.8, 0.8))
leu.uptake.oct.plot
# No real difference between the replicates.
# We'll combine them, no need to use the different
# shapes.
leu.uptake.oct.plot <- leu.uptake.oct %>%
  mutate(depth = as.factor(paste(depth,
                                 "m"))) %>%
  ggplot(aes(y = µgBCP_per_L_hr,
             x = treatment,
             col = treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1) +
  facet_wrap(~depth) +
  scale_color_manual(values = color.vector) +
  ylim(c(0, 25)) +
  ylab("Bacterial carbon production (µgC/L/hr)") +
  theme_bw() +
  theme(legend.position = c(0.8, 0.8))
leu.uptake.oct.plot


#### Visualize both sets of plots ####
pdf("results/microbial/leucine_uptake.pdf",
    width = 8,
    height = 6)
ggarrange(leu.uptake.8hr.plot,
          leu.uptake.oct.plot,
          ncol = 1,
          labels = c("A. September 10th, 2021",
                     "B. October 14th, 2021"))
dev.off()
