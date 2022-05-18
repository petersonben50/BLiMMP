#### code/incubations/Kmet_overview.R ####
# Written for BLiMMP project
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP")
library(ggpubr)
library(lubridate)
library(readxl)
library(tidyverse)
cb.translator <- readRDS("references/colorblind_friendly_colors.rds")



#### Read in data ####
Hg.Kmet.data <- read.csv("dataEdited/incubation_Hg_rate_data.csv") %>%
  mutate(sampleInfo = paste(year(startDate), "-",
                            month(startDate,label = TRUE,abbr = FALSE), "-",
                            depth, "m", sep = "")) %>%
  filter(year(startDate) %in% c(2020, 2021)) %>%
  mutate(treatment = fct_relevel(treatment,
                                  c("unfiltered-unamended",
                                    "unfiltered-molybdate")),
         sampleInfo = fct_relevel(sampleInfo, c("2020-September-11m", "2020-September-15.5m", "2020-September-20.7m",
                                                "2020-October-15.7m", "2020-October-20.9m",
                                                "2021-September-10.8m", "2021-September-11.9m", "2021-September-19.9m",
                                                "2021-October-14.2m", "2021-October-15.2m", "2021-October-19.9m")))



#### Look for differences between Kmet calculations for ambient samples ####
ambient.conditions <- Hg.Kmet.data %>%
  filter(treatment == "unfiltered-unamended") %>%
  select(startDate, depth, sampleInfo, Kmet_t1, Kmet_t2, Kmet_total) %>%
  gather(key = calculation,
         value = Kmet,
         4:6) %>%
  ggplot(aes(x = sampleInfo,
             y = Kmet,
             color = calculation)) +
  geom_boxplot() +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        axis.title.x = element_blank()) +
  ylim(c(-0.03, 0.25)) +
  ylab("Kmet (day^-1)") +
  labs(title = "Kmet of ambient conditions")



#### Look for differences between Kmet calculations for molybdate-amended samples ####
molybdate.amended <- Hg.Kmet.data %>%
  filter(treatment == "unfiltered-molybdate") %>%
  select(startDate, depth, sampleInfo, Kmet_t1, Kmet_t2, Kmet_total) %>%
  gather(key = calculation,
         value = Kmet,
         4:6) %>%
  ggplot(aes(x = sampleInfo,
             y = Kmet,
             color = calculation)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        axis.title.x = element_blank()) +
  ylim(c(-0.03, 0.25)) +
  ylab("Kmet (day^-1)") +
  labs(title = "Kmet of molybdate-inhibited samples")


#### Save out figure ####
pdf("results/incubations/Kmet_overview.pdf",
    width = 10,
    heigh = 5)
ggarrange(ambient.conditions, molybdate.amended)
dev.off()



#### Look for differences between Kmet calculations for filtered samples ####
filtered.control <- Hg.Kmet.data %>%
  filter(treatment == "filtered-unamended") %>%
  select(startDate, depth, sampleInfo, Kmet_t1, Kmet_t2, Kmet_total) %>%
  gather(key = calculation,
         value = Kmet,
         4:6) %>%
  ggplot(aes(x = sampleInfo,
             y = Kmet,
             color = calculation)) +
  geom_boxplot() +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        axis.title.x = element_blank()) +
  ylim(c(-0.03, 0.03)) +
  labs(title = "Kmet of ambient conditions")
filtered.control
