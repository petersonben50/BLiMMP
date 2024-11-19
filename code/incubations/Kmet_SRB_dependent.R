#### code/incubations/Kmet_SRB_dependent.R ####
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



#### Read in water column sulfide data ####
sulfide.data <- read.csv("dataFinal/water_chem_data.csv") %>%
  group_by(date, depth) %>%
  summarise(sulfide_uM = mean(sulfide_uM, na.rm = TRUE))



#### Guild color vector ####
color.vector <- cb.translator[c("blue", "reddishpurple")]
names(color.vector) <- c("SRB", "non-SRB")



#### Year shape vector ####
shape.vector <- c(21, 24)
# shape.vector <- c(1, 2)
names(shape.vector) <- c(2020, 2021)



#### Read in data, calculate SRB-dependent Kmet ####
Hg.Kmet.data <- read.csv("dataEdited/incubation_Hg_rate_data.csv") %>%
  rename(date = startDate) %>%
  filter(treatment != "filtered-unamended",
         year(date) %in% c(2020, 2021)) %>%
  mutate(treatment = treatment %>%
           gsub("unfiltered-", "", .) %>%
           gsub("unamended", "ambient", .)) %>%
  select(sampleID, date, depth, treatment, Kmet_t1) %>%
  group_by(sampleID, date, depth, treatment) %>%
  summarise(Kmet_mean = mean(Kmet_t1),
            Kmet_sd = sd(Kmet_t1),
            Kmet_count = n(),
            Kmet_se = Kmet_sd / Kmet_count) %>%
  ungroup() %>%
  gather(key = constituent,
         value = value,
         -c(1:4)) %>%
  mutate(constituent = paste(treatment, "_", constituent,
                             sep = "")) %>%
  select(sampleID, date, depth, constituent, value) %>%
  spread(key = constituent,
         value = value) %>%
  mutate(SRB_Kmet_mean = ambient_Kmet_mean - molybdate_Kmet_mean,
         SRB_Kmet_se = ambient_Kmet_se + molybdate_Kmet_se) %>%
  left_join(sulfide.data) %>%
  select(date, depth, sulfide_uM, SRB_Kmet_mean, SRB_Kmet_se, molybdate_Kmet_mean, molybdate_Kmet_se) %>%
  gather(key = measurement,
         value = Kmet_value,
         -c(1:3)) %>%
  mutate(guild = measurement %>%
           strsplit("_") %>% sapply("[", 1),
         measurement = paste(measurement %>%
                               strsplit("_") %>% sapply("[", 2),
                             measurement %>%
                               strsplit("_") %>% sapply("[", 3),
                             sep = "_")) %>%
  mutate(guild = gsub("molybdate", "non-SRB", guild)) %>%
  spread(key = measurement,
         value = Kmet_value)


#### Plot: Kmet vs. sulfide ####
Kmet.plot <- Hg.Kmet.data %>%
  ggplot(aes(x = sulfide_uM,
             y = Kmet_mean,
             bg = guild,
             shape = as.character(year(date)))) +
  geom_errorbar(aes(ymin = Kmet_mean - Kmet_se,
                    ymax = Kmet_mean + Kmet_se),
                width = 1.5) +
  geom_point(size = 2.5,
             # stroke = 1.2
             ) +
  scale_fill_manual(values = color.vector, name = "Functional guild") +
  scale_shape_manual(values = shape.vector, name = "Year") +
  theme_bw() +
  ylab(bquote(K[met] ())) +
  xlab("Sulfide (µM)") +
  ylim(c(-0.005, 0.15)) +
  theme(legend.position = c(0.175, 0.75),
        legend.spacing.y = unit(1, "mm"),
        legend.key.height = unit(5, "mm"),
        # legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black",
                                             fill = "white")) +
  guides(fill = guide_legend(override.aes = list(fill = color.vector)))



#### Save out plot ####
pdf("results/incubations/Kmet_SRB_dependent.pdf",
    width = 5,
    height = 2.5)
Kmet.plot
dev.off()
