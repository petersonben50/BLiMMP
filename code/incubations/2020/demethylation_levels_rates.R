#### code/incubations/2020/demethylation_levels_rates.R ####
# Written for BLiMMP project
# Benjamin D. Peterson


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
  select(sampleInfo, incubationID, treatment, t, durationSinceTimepointInDays,fraction_MeHg_204)


#### Plot ambient MeHg over incubation ####
pdf("results/incubations/2020incubations/demethylation/Me204Hg_levels_over_incubation.pdf",
    width = 8,
    height = 6)
Hg.data.plotting %>%
  ggplot(aes(x = durationSinceTimepointInDays,
             y = fraction_MeHg_204,
             color = treatment)) +
  geom_point() +
  geom_line(aes(group = incubationID)) +
  scale_color_manual(values = color.vector) +
  facet_wrap(~sampleInfo, nrow = 2) +
  theme_classic() +
  ylim(c(0, 1.5))
dev.off()


#### Read in rate data ####
all.MeHg.delta.data <- read.csv("dataEdited/incubations/2020incubations_MeHg_delta.csv") %>%
  mutate(sampleInfo = paste(monthOfIncubations,
                            depth,
                            sep = "-"),
         sampleInfo = fct_relevel(sampleInfo,
                                  c("September-11m", "September-15.5m", "September-20.7m",
                                    "October-15.7m", "October-20.9m")))


#### Plot change in Me204Hg over each time point ####
all.MeHg.delta.data %>%
  ggplot(aes(x = t,
             y = net_Me204Hg_fraction.day,
             color = treatment)) +
  geom_point() +
  geom_line(aes(group = incubationID)) +
  scale_color_manual(values = color.vector) +
  facet_wrap(~sampleInfo, nrow = 2) +
  theme_classic() +
  ylim(c(-0.6, 0.4))



#### Boxplots of methylation rates by month ####
# Boxplots for September
all.MeHg.delta.data %>%
  filter(monthOfIncubations == "September") %>%
  mutate(treatment = fct_relevel(treatment,
                                 names(color.vector)[1:3])) %>%
  ggplot(aes(x = depth,
             y = net_Me204Hg_fraction.day,
             color = treatment)) +
  geom_boxplot() +
  scale_color_manual(values = color.vector) +
  facet_wrap(~t, nrow = 2) +
  theme_classic() +
  ylim(c(-0.6, 0.2))
# Boxplots for October
all.MeHg.delta.data %>%
  filter(monthOfIncubations == "October") %>%
  filter(treatment %in% names(color.vector)[1:3]) %>%
  mutate(treatment = fct_relevel(treatment,
                                 names(color.vector)[1:3])) %>%
  ggplot(aes(x = depth,
             y = net_Me204Hg_fraction.day,
             color = treatment)) +
  geom_boxplot() +
  scale_color_manual(values = color.vector) +
  facet_wrap(~t, nrow = 2) +
  theme_classic() +
  ylim(c(-0.25, 0.4))


#### Calculate methylation rate constants ####
starting.HgT <- Hg.data %>%
  filter(t == "t0") %>%
  rename(HgT_start = excess_HgT_204_ng.L) %>%
  select(incubationID, HgT_start)

meth.rate.constant <- all.MeHg.delta.data %>%
  filter(t == "t1") %>%
  select(incubationID, sampleInfo, treatment,
         net_Me198Hg_ng.L.day) %>%
  left_join(starting.HgT) %>%
  mutate(Kmet = net_Me198Hg_ng.L.day / HgT_start)


#### Plot methylation rate constants ####
pdf("results/incubations/2020incubations/methylation/Kdem.pdf",
    width = 8,
    height = 4)
meth.rate.constant %>%
  mutate(sampleInfo = fct_relevel(sampleInfo,
                                  c("September-11m", "September-15.5m", "September-20.7m",
                                    "October-15.7m", "October-20.9m")),
         treatment = fct_relevel(treatment,
                                 names(color.vector))) %>%
  ggplot(aes(x = treatment,
             y = Kmet,
             color = treatment)) +
  geom_boxplot() +
  scale_color_manual(values = color.vector) +
  ylim(c(-0.01, 0.25)) +
  facet_wrap(~sampleInfo, nrow = 2) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  xlab(element_blank()) +
  ylab("Kmet (day^-1)")
dev.off()


#### Calculate different between unamended and molybdate-inhibited methylation ####
meth.rate.constant %>%
  group_by(sampleInfo, treatment) %>%
  summarize(Kmet = mean(Kmet)) %>%
  spread(key = treatment,
         value = Kmet) %>%
  mutate(molybdateInhibitionPercent = (`unfiltered-unamended` - `unfiltered-molybdate`) / `unfiltered-unamended` * 100)


#### Statistical testing ####
testing <- aov(Kmet ~ treatment * sampleInfo,
               data = meth.rate.constant)
summary(testing)
