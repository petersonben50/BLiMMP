#### code/incubations/2020/Hg_data_processing.R ####
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
# color.vector<- c(color.vector,
#                  cb.translator["yellow"], cb.translator["reddishpurple"])
# names(color.vector) <- c(names(color.vector),
#                          "unfiltered-glucose", "unfiltered-glucose+molybdate")


#### Read in data ####
Hg.data <- read.csv("dataEdited/incubations/2020incubations_Hg_data.csv") %>%
  mutate(depth = paste(depth, "meters")) %>%
  filter(treatment %in% names(color.vector))
Hg.data$monthOfIncubations <- "October"
Hg.data$monthOfIncubations[which(Hg.data$tripID == "BLiMMP_trip_010")] = "September"
Hg.data$monthOfIncubations <- relevel(as.factor(Hg.data$monthOfIncubations),
                                      c("September"))


#### Plot HgT for 198 and 204 at each time point ####
# Plot HgT levels for 198
pdf(file = "results/incubations/2020incubations/198HgT_2020incubations.pdf",
    width = 8,
    height = 6)
Hg.data %>%
  ggplot(aes(x = t,
             y = excess_HgT_198_ng.L,
             color = treatment)) +
  geom_point() +
  geom_line(aes(group = incubationID)) +
  scale_color_manual(values = color.vector) +
  facet_wrap(~monthOfIncubations + depth, nrow = 2) +
  theme_classic() +
  ylim(c(0, 1.1))
dev.off()

# Plot HgT levels for 204
pdf(file = "results/incubations/2020incubations/204HgT_levels/204HgT_2020incubations.pdf",
    width = 8,
    height = 6)
Hg.data %>%
  ggplot(aes(x = t,
             y = excess_HgT_204_ng.L,
             color = treatment)) +
  geom_point() +
  geom_line(aes(group = incubationID)) +
  scale_color_manual(values = color.vector) +
  facet_wrap(~monthOfIncubations + depth, nrow = 2) +
  theme_classic() +
  ylim(c(0, 1.4)) +
  ylab("Excess 204HgT (ng/L)") +
  xlab("Time point")
dev.off()
# Look at MeHg fraction of 204HgT 
pdf(file = "results/incubations/2020incubations/204HgT_levels/Me204Hg_fraction_2020incubations.pdf",
    width = 8,
    height = 6)
Hg.data %>%
  ggplot(aes(x = t,
             y = fraction_MeHg_204,
             color = treatment)) +
  geom_point() +
  geom_line(aes(group = incubationID)) +
  scale_color_manual(values = color.vector) +
  facet_wrap(~monthOfIncubations + depth, nrow = 2) +
  theme_classic() +
  xlab("Time point") +
  ylab("Fraction Me204Hg") +
  ylim(c(0, 1.4))
dev.off()
# Look at excess Me204Hg levels
pdf(file = "results/incubations/2020incubations/204HgT_levels/Me204Hg_2020incubations.pdf",
    width = 8,
    height = 6)
Hg.data %>%
  ggplot(aes(x = t,
             y = excess_MeHg_204_ng.L,
             color = treatment)) +
  geom_point() +
  geom_line(aes(group = incubationID)) +
  scale_color_manual(values = color.vector) +
  facet_wrap(~monthOfIncubations + depth, nrow = 2) +
  theme_bw() +
  xlab("Time point") +
  ylab("Excess Me204Hg (ng/L)") +
  ylim(c(0, 1.4))
dev.off()


#### Calculate loss of HgT in 198 and 204 over each time point ####
# Calculate loss of 198HgT
HgT.delta.198 <- Hg.data %>%
  select(incubationID, t, excess_HgT_198_ng.L) %>%
  spread(key = t,
         value = excess_HgT_198_ng.L) %>%
  mutate(t2 = t2 - t1,
         t1 = t1 - t0) %>%
  select(incubationID, t1, t2) %>%
  gather(key = t,
         value = net_198HgT_ng.L,
         -1) %>%
  left_join(Hg.data) %>%
  mutate(net_198HgT_ng.L.day = round(net_198HgT_ng.L / durationSinceTimepointInDays, 3)) %>%
  select(incubationID, monthOfIncubations, depth, treatment, t, net_198HgT_ng.L, net_198HgT_ng.L.day) %>%
  filter(!is.na(net_198HgT_ng.L.day))
# Calculate loss of 204HgT
HgT.delta.204 <- Hg.data %>%
  select(incubationID, t, excess_HgT_204_ng.L) %>%
  spread(key = t,
         value = excess_HgT_204_ng.L) %>%
  mutate(t2 = t2 - t1,
         t1 = t1 - t0) %>%
  select(incubationID, t1, t2) %>%
  gather(key = t,
         value = net_204HgT_ng.L,
         -1) %>%
  left_join(Hg.data) %>%
  mutate(net_204HgT_ng.L.day = round(net_204HgT_ng.L / durationSinceTimepointInDays, 3)) %>%
  select(incubationID, monthOfIncubations, depth, treatment, t, net_204HgT_ng.L, net_204HgT_ng.L.day) %>%
  filter(!is.na(net_204HgT_ng.L.day))
# Combine rates of loss for each isotope
all.HgT.delta.data <- full_join(HgT.delta.198,
                                HgT.delta.204)
rm(HgT.delta.198,
   HgT.delta.204)


#### Plot rates of loss of HgT ####
# Losses in 198HgT
all.HgT.delta.data %>%
  ggplot(aes(x = t,
             y = net_198HgT_ng.L.day,
             color = treatment)) +
  geom_point() +
  geom_line(aes(group = incubationID)) +
  scale_color_manual(values = color.vector) +
  facet_wrap(~monthOfIncubations + depth, nrow = 2) +
  theme_classic() +
  ylim(c(-0.5, 0.2))
# Losses in 204HgT
all.HgT.delta.data %>%
  ggplot(aes(x = t,
             y = net_204HgT_ng.L.day,
             color = treatment)) +
  geom_point() +
  geom_line(aes(group = incubationID)) +
  scale_color_manual(values = color.vector) +
  facet_wrap(~monthOfIncubations + depth, nrow = 2) +
  theme_classic() +
  ylim(c(-0.5, 0.2))


#### Summarize the changes in HgT ####
all.HgT.delta.data %>%
  group_by(t) %>%
  summarise(net_198HgT_ng.L.day.mean = mean(net_198HgT_ng.L.day),
            net_198HgT_ng.L.day.sd = sd(net_198HgT_ng.L.day),
            net_204HgT_ng.L.day.mean = mean(net_204HgT_ng.L.day),
            net_204HgT_ng.L.day.sd = sd(net_204HgT_ng.L.day))
mean.data <- all.HgT.delta.data %>%
  group_by(t, depth, treatment) %>%
  summarise(net_198HgT_ng.L.day.mean = mean(net_198HgT_ng.L.day),
            net_198HgT_ng.L.day.sd = sd(net_198HgT_ng.L.day),
            net_204HgT_ng.L.day.mean = mean(net_204HgT_ng.L.day),
            net_204HgT_ng.L.day.sd = sd(net_204HgT_ng.L.day))
rm(mean.data)



#### Check for differences between different variables ####
# 198HgT
anova.test.of.Hg <- aov(net_198HgT_ng.L.day ~ depth * treatment * t,
                        data = all.HgT.delta.data)
summary(anova.test.of.Hg)
rm(anova.test.of.Hg)
# 204HgT
anova.test.of.Hg <- aov(net_204HgT_ng.L.day ~ depth * treatment * t,
                        data = all.HgT.delta.data)
summary(anova.test.of.Hg)
rm(anova.test.of.Hg)



##########################
#### MeHg analyses ####


#### Plot MeHg at each time point ####
# Plot 198 as concentration
Hg.data %>%
  ggplot(aes(x = durationInDays,
             y = excess_MeHg_198_ng.L,
             color = treatment)) +
  geom_point() +
  geom_line(aes(group = incubationID)) +
  scale_color_manual(values = color.vector) +
  facet_wrap(~monthOfIncubations + depth, nrow = 2) +
  theme_classic() +
  ylim(c(-0.01, 0.5)) +
  xlab("Days") +
  ylab("Excess Me198Hg (ng/L)")
# Plot 198 as fraction MeHg
Hg.data %>%
  ggplot(aes(x = durationInDays,
             y = fraction_MeHg_198,
             color = treatment)) +
  geom_point() +
  geom_line(aes(group = incubationID)) +
  scale_color_manual(values = color.vector) +
  facet_wrap(~monthOfIncubations + depth, nrow = 2) +
  theme_classic() +
  ylim(c(-0.01, 0.7)) +
  xlab("Days") +
  ylab("Fraction Me198Hg")
# Plot 204 as fraction MeHg
Hg.data %>%
  ggplot(aes(x = durationInDays,
             y = fraction_MeHg_204,
             color = treatment)) +
  geom_point() +
  geom_line(aes(group = incubationID)) +
  scale_color_manual(values = color.vector) +
  facet_wrap(~monthOfIncubations + depth, nrow = 2) +
  theme_classic() +
  ylim(c(-0.01, 1.4)) +
  xlab("Days") +
  ylab("Fraction Me204Hg")



#### Calculate rates of change in 198 and 204, both concentrations and fraction ####
# Change in 198 MeHg concentration
MeHg.delta.198 <- Hg.data %>%
  select(incubationID, t, excess_MeHg_198_ng.L) %>%
  spread(key = t,
         value = excess_MeHg_198_ng.L) %>%
  mutate(t2 = t2 - t1,
         t1 = t1 - t0) %>%
  select(incubationID, t1, t2) %>%
  gather(key = t,
         value = net_Me198Hg_ng.L,
         -1) %>%
  left_join(Hg.data) %>%
  mutate(net_Me198Hg_ng.L.day = round(net_Me198Hg_ng.L / durationSinceTimepointInDays, 3)) %>%
  select(incubationID, monthOfIncubations, depth, treatment, t, net_Me198Hg_ng.L, net_Me198Hg_ng.L.day) %>%
  filter(!is.na(net_Me198Hg_ng.L.day))
# Change in 198 MeHg fraction
fraction.MeHg.delta.198 <- Hg.data %>%
  select(incubationID, t, fraction_MeHg_198) %>%
  spread(key = t,
         value = fraction_MeHg_198) %>%
  mutate(t2 = t2 - t1,
         t1 = t1 - t0) %>%
  select(incubationID, t1, t2) %>%
  gather(key = t,
         value = net_Me198Hg_fraction,
         -1) %>%
  left_join(Hg.data) %>%
  mutate(net_Me198Hg_fraction.day = round(net_Me198Hg_fraction / durationSinceTimepointInDays, 3)) %>%
  select(incubationID, monthOfIncubations, depth, treatment, t, net_Me198Hg_fraction, net_Me198Hg_fraction.day) %>%
  filter(!is.na(net_Me198Hg_fraction.day))
# Change in 204 MeHg concentration
MeHg.delta.204 <- Hg.data %>%
  select(incubationID, t, excess_MeHg_204_ng.L) %>%
  spread(key = t,
         value = excess_MeHg_204_ng.L) %>%
  mutate(t2 = t2 - t1,
         t1 = t1 - t0) %>%
  select(incubationID, t1, t2) %>%
  gather(key = t,
         value = net_Me204Hg_ng.L,
         -1) %>%
  left_join(Hg.data) %>%
  mutate(net_Me204Hg_ng.L.day = round(net_Me204Hg_ng.L / durationSinceTimepointInDays, 3)) %>%
  select(incubationID, monthOfIncubations, depth, treatment, t, net_Me204Hg_ng.L, net_Me204Hg_ng.L.day) %>%
  filter(!is.na(net_Me204Hg_ng.L.day))
# Change in 204 MeHg fraction
fraction.MeHg.delta.204 <- Hg.data %>%
  select(incubationID, t, fraction_MeHg_204) %>%
  spread(key = t,
         value = fraction_MeHg_204) %>%
  mutate(t2 = t2 - t1,
         t1 = t1 - t0) %>%
  select(incubationID, t1, t2) %>%
  gather(key = t,
         value = net_Me204Hg_fraction,
         -1) %>%
  left_join(Hg.data) %>%
  mutate(net_Me204Hg_fraction.day = round(net_Me204Hg_fraction / durationSinceTimepointInDays, 3)) %>%
  select(incubationID, monthOfIncubations, depth, treatment, t, net_Me204Hg_fraction, net_Me204Hg_fraction.day) %>%
  filter(!is.na(net_Me204Hg_fraction.day))
# Combine Hg species fraction and concentration
all.MeHg.delta.data <- full_join(MeHg.delta.198,
                           fraction.MeHg.delta.198) %>%
  full_join(MeHg.delta.204) %>%
  full_join(fraction.MeHg.delta.204)
rm(MeHg.delta.198, fraction.MeHg.delta.198,
   MeHg.delta.204, fraction.MeHg.delta.204)
# Save out data
write.csv(all.MeHg.delta.data,
          "dataEdited/incubations/2020incubations_MeHg_delta.csv",
          row.names = FALSE)



#### Plot rate of change in concentration and fraction within each incubation ####
# Plot as a point, connecting the two time points within a bag.
# MeHg concentration
all.MeHg.delta.data %>%
  ggplot(aes(x = t,
             y = net_Me198Hg_ng.L.day,
             color = treatment)) +
  geom_point() +
  geom_line(aes(group = incubationID)) +
  scale_color_manual(values = color.vector) +
  facet_wrap(~monthOfIncubations + depth, nrow = 2) +
  theme_classic() +
  ylim(c(-0.01, 0.2))
# MeHg fraction 
all.MeHg.delta.data %>%
  ggplot(aes(x = t,
             y = net_Me198Hg_fraction.day,
             color = treatment)) +
  geom_point() +
  geom_line(aes(group = incubationID)) +
  scale_color_manual(values = color.vector) +
  facet_wrap(~monthOfIncubations + depth, nrow = 2) +
  theme_classic() +
  ylim(c(-0.01, 0.25))


#### Boxplots of methylation rates by month ####
# Boxplots for September
all.MeHg.delta.data %>%
  filter(monthOfIncubations == "September") %>%
  mutate(treatment = fct_relevel(treatment,
                                 names(color.vector)[1:3])) %>%
  ggplot(aes(x = depth,
             y = net_Me198Hg_fraction.day,
             color = treatment)) +
  geom_boxplot() +
  scale_color_manual(values = color.vector) +
  facet_wrap(~t, nrow = 2) +
  theme_classic() +
  ylim(c(-0.01, 0.22))
# Boxplots for October
all.MeHg.delta.data %>%
  filter(monthOfIncubations == "October") %>%
  filter(treatment %in% names(color.vector)[1:3]) %>%
  mutate(treatment = fct_relevel(treatment,
                                 names(color.vector)[1:3])) %>%
  ggplot(aes(x = depth,
             y = net_Me198Hg_fraction.day,
             color = treatment)) +
  geom_boxplot() +
  scale_color_manual(values = color.vector) +
  facet_wrap(~t, nrow = 2) +
  theme_classic() +
  ylim(c(-0.01, 0.25))


# #### Statistical analyses ####
fraction.MeHg.delta.198.sept <- all.MeHg.delta.data %>%
  filter(monthOfIncubations == "September")
anova.sept <- aov(net_Me198Hg_fraction.day ~ treatment * depth * t,
                  data = fraction.MeHg.delta.198.sept)
summary(anova.sept)
#TukeyHSD(anova.sept)

fraction.MeHg.delta.198.oct <- all.MeHg.delta.data %>%
  filter(monthOfIncubations == "October")
anova.oct <- aov(net_Me198Hg_fraction.day ~ treatment * depth * t,
                 data = fraction.MeHg.delta.198.oct)
summary(anova.oct)
rm(fraction.MeHg.delta.198.sept, anova.sept,
   fraction.MeHg.delta.198.oct, anova.oct)


#### Calculate fraction reduction in fraction MeHg produced due to molybdate ####
MoO4.inhibition <- all.MeHg.delta.data %>%
  filter(treatment %in% names(color.vector)[2:3]) %>%
  group_by(monthOfIncubations, depth, treatment, t) %>%
  summarise(mean_methylation = mean(net_Me198Hg_ng.L.day)) %>%
  spread(key = treatment,
         value = mean_methylation) %>%
  mutate(fraction_reduction = (`unfiltered-unamended` - `unfiltered-molybdate`) / `unfiltered-unamended`)
