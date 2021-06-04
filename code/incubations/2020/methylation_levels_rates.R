#### code/incubations/2020/MeHg_levels_rates.R ####
# Written for BLiMMP project
# Benjamin D. Peterson


#### Prep workspace ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP/")
library(dplyr)
library(emmeans)
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
  select(sampleInfo, incubationID, treatment, t, durationInDays, durationSinceTimepointInDays, excess_MeHg_198_ng.L, fraction_MeHg_198, above_DDL_MeHg_198)


#### Generate vector for DDL points ####
point.vector.to.use <- c(16, 4)
names(point.vector.to.use) <- c(TRUE, FALSE)


#### Plot Me198Hg levels over incubation ####
pdf("results/incubations/2020incubations/methylation/MeHg_levels_over_incubation.pdf",
    width = 8,
    height = 6)
Hg.data.plotting %>%
  ggplot(aes(x = durationInDays,
             y = excess_MeHg_198_ng.L,
             color = treatment,
             shape = above_DDL_MeHg_198)) +
  geom_point() +
  geom_line(aes(group = incubationID)) +
  scale_color_manual(values = color.vector) +
  scale_shape_manual(values = point.vector.to.use) +
  facet_wrap(~sampleInfo, nrow = 2) +
  theme_classic() +
  ylim(c(-0.005, 0.5))
dev.off()


#### Plot Me198Hg levels over incubation, zoomed in on metalimion ####
pdf("results/incubations/2020incubations/methylation/MeHg_levels_over_incubation_metalimnion.pdf",
    width = 3,
    height = 6)
Hg.data.plotting %>%
  filter(sampleInfo %in% c("October-15.7m",
                           "September-11m")) %>%
  ggplot(aes(x = durationInDays,
             y = excess_MeHg_198_ng.L,
             color = treatment,
             shape = above_DDL_MeHg_198)) +
  geom_point() +
  geom_line(aes(group = incubationID)) +
  scale_color_manual(values = color.vector) +
  scale_shape_manual(values = point.vector.to.use) +
  facet_wrap(~sampleInfo, nrow = 2) +
  theme_classic() +
  ylim(c(-0.01, 0.025)) +
  theme(legend.position = "none")
dev.off()


#### Calculate means and SDs ####
for.inspection <- Hg.data.plotting %>%
  group_by(sampleInfo, treatment, t) %>%
  summarise(mean_198 = mean(excess_MeHg_198_ng.L),
            sd_198 = sd(excess_MeHg_198_ng.L))



#### Read in rate data ####
all.MeHg.delta.data <- read.csv("dataEdited/incubations/2020incubations_MeHg_delta.csv") %>%
  mutate(sampleInfo = paste(monthOfIncubations,
                            depth,
                            sep = "-"),
         sampleInfo = fct_relevel(sampleInfo,
                                  c("September-11m", "September-15.5m", "September-20.7m",
                                    "October-15.7m", "October-20.9m")))


#### Plot change in Me198Hg over each time point ####
all.MeHg.delta.data %>%
  ggplot(aes(x = t,
             y = net_Me198Hg_fraction.day,
             color = treatment)) +
  geom_point() +
  geom_line(aes(group = incubationID)) +
  scale_color_manual(values = color.vector) +
  facet_wrap(~sampleInfo, nrow = 2) +
  theme_classic() +
  ylim(c(-0.02, 0.3))


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


#### Calculate methylation rate constants ####
starting.HgT <- Hg.data %>%
  filter(t == "t0") %>%
  rename(HgT_start = excess_HgT_198_ng.L) %>%
  select(incubationID, HgT_start)

meth.rate.constant <- all.MeHg.delta.data %>%
  filter(t == "t1") %>%
  select(incubationID, sampleInfo, treatment,
         net_Me198Hg_ng.L.day) %>%
  left_join(starting.HgT) %>%
  mutate(Kmet = net_Me198Hg_ng.L.day / HgT_start)



#### Plot methylation rate constants ####
pdf("results/incubations/2020incubations/methylation/Kmet.pdf",
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


#### ANOVA test ####
anova.incubations <- aov(Kmet ~ treatment * sampleInfo,
                         data = meth.rate.constant)
summary(anova.incubations)


#### Planned comparisons using emmeans ####
Simple.Effects.By.Type <- emmeans(anova.incubations, ~treatment|sampleInfo)
pairs(Simple.Effects.By.Type,
      adjust = 'none')


#### Individual ANOVAs ####
date.depths <- as.character(unique(meth.rate.constant$sampleInfo))
for (date.depth in date.depths) {
  print(date.depth)
  temp.data <- meth.rate.constant %>%
    filter(sampleInfo == date.depth)
  print(summary(aov(Kmet ~ treatment,
              data = temp.data)))
}


#### Calculate different between unamended and molybdate-inhibited methylation ####
molybdate.inhibition <- meth.rate.constant %>%
  group_by(sampleInfo, treatment) %>%
  summarize(Kmet = mean(Kmet)) %>%
  spread(key = treatment,
         value = Kmet) %>%
  mutate(molybdateInhibitionPercent = (`unfiltered-unamended` - `unfiltered-molybdate`) / `unfiltered-unamended` * 100)
