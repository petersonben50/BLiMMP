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
  select(sampleInfo, incubationID, treatment, t, durationSinceTimepointInDays, excess_HgT_198_ng.L, excess_HgT_204_ng.L)



#### Plot HgT for 198 at each time point ####
# Plot HgT levels for 198
pdf(file = "results/incubations/2020incubations/HgT_data/198HgT_levels.pdf",
    width = 8,
    height = 6)
Hg.data.plotting %>%
  ggplot(aes(x = t,
             y = excess_HgT_198_ng.L,
             color = treatment)) +
  geom_point() +
  geom_line(aes(group = incubationID)) +
  scale_color_manual(values = color.vector) +
  facet_wrap(~sampleInfo, nrow = 2) +
  theme_classic() +
  ylim(c(0, 1.3))
dev.off()


#### Plot HgT for 204 at each time point ####
# Plot HgT levels for 204
pdf(file = "results/incubations/2020incubations/HgT_data/204HgT_levels.pdf",
    width = 8,
    height = 6)
Hg.data.plotting %>%
  ggplot(aes(x = t,
             y = excess_HgT_204_ng.L,
             color = treatment)) +
  geom_point() +
  geom_line(aes(group = incubationID)) +
  scale_color_manual(values = color.vector) +
  facet_wrap(~sampleInfo, nrow = 2) +
  theme_classic() +
  ylim(c(0, 1.3))
dev.off()






#### Calculate loss of HgT in 198 and 204 over each time point ####
# Calculate loss of 198HgT
HgT.delta.198 <- Hg.data.plotting %>%
  select(incubationID, t, excess_HgT_198_ng.L) %>%
  spread(key = t,
         value = excess_HgT_198_ng.L) %>%
  mutate(t2 = t2 - t1,
         t1 = t1 - t0) %>%
  select(incubationID, t1, t2) %>%
  gather(key = t,
         value = net_198HgT_ng.L,
         -1) %>%
  left_join(Hg.data.plotting) %>%
  mutate(net_198HgT_ng.L.day = round(net_198HgT_ng.L / durationSinceTimepointInDays, 3)) %>%
  select(incubationID, sampleInfo, treatment, t, net_198HgT_ng.L, net_198HgT_ng.L.day) %>%
  filter(!is.na(net_198HgT_ng.L.day))
# Calculate loss of 204HgT
HgT.delta.204 <- Hg.data.plotting %>%
  select(incubationID, t, excess_HgT_204_ng.L) %>%
  spread(key = t,
         value = excess_HgT_204_ng.L) %>%
  mutate(t2 = t2 - t1,
         t1 = t1 - t0) %>%
  select(incubationID, t1, t2) %>%
  gather(key = t,
         value = net_204HgT_ng.L,
         -1) %>%
  left_join(Hg.data.plotting) %>%
  mutate(net_204HgT_ng.L.day = round(net_204HgT_ng.L / durationSinceTimepointInDays, 3)) %>%
  select(incubationID, sampleInfo, treatment, t, net_204HgT_ng.L, net_204HgT_ng.L.day) %>%
  filter(!is.na(net_204HgT_ng.L.day))
# Combine rates of loss for each isotope
all.HgT.delta.data <- full_join(HgT.delta.198,
                                HgT.delta.204)
rm(HgT.delta.198,
   HgT.delta.204)




#### Plot rates of loss of HgT ####

# Losses in 198HgT
pdf(file = "results/incubations/2020incubations/HgT_data/198HgT_loss.pdf",
    width = 8,
    height = 6)
all.HgT.delta.data %>%
  ggplot(aes(x = t,
             y = net_198HgT_ng.L.day,
             color = treatment)) +
  geom_point() +
  geom_line(aes(group = incubationID)) +
  scale_color_manual(values = color.vector) +
  facet_wrap(~sampleInfo, nrow = 2) +
  theme_classic() +
  ylim(c(-0.5, 0.2))
dev.off()

# Losses in 204HgT
pdf(file = "results/incubations/2020incubations/HgT_data/204HgT_loss.pdf",
    width = 8,
    height = 6)
all.HgT.delta.data %>%
  ggplot(aes(x = t,
             y = net_204HgT_ng.L.day,
             color = treatment)) +
  geom_point() +
  geom_line(aes(group = incubationID)) +
  scale_color_manual(values = color.vector) +
  facet_wrap(~sampleInfo, nrow = 2) +
  theme_classic() +
  ylim(c(-0.5, 0.2))
dev.off()




#### Summarize the changes in HgT ####
all.HgT.delta.data %>%
  group_by(t) %>%
  summarise(net_198HgT_ng.L.day.mean = mean(net_198HgT_ng.L.day),
            net_198HgT_ng.L.day.sd = sd(net_198HgT_ng.L.day),
            net_204HgT_ng.L.day.mean = mean(net_204HgT_ng.L.day),
            net_204HgT_ng.L.day.sd = sd(net_204HgT_ng.L.day))
mean.data <- all.HgT.delta.data %>%
  group_by(t, sampleInfo, treatment) %>%
  summarise(net_198HgT_ng.L.day.mean = mean(net_198HgT_ng.L.day),
            net_198HgT_ng.L.day.sd = sd(net_198HgT_ng.L.day),
            net_204HgT_ng.L.day.mean = mean(net_204HgT_ng.L.day),
            net_204HgT_ng.L.day.sd = sd(net_204HgT_ng.L.day))
rm(mean.data)





#### Check for differences between different variables ####
# 198HgT
anova.test.of.Hg <- aov(net_198HgT_ng.L.day ~ sampleInfo * treatment * t,
                        data = all.HgT.delta.data)
summary(anova.test.of.Hg)
rm(anova.test.of.Hg)
# 204HgT
anova.test.of.Hg <- aov(net_204HgT_ng.L.day ~ sampleInfo * treatment * t,
                        data = all.HgT.delta.data)
summary(anova.test.of.Hg)
rm(anova.test.of.Hg)

