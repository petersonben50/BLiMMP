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
Hg.Kmet.data <- read.csv("dataFinal/incubation_Hg_rate_data.csv") %>%
  mutate(sampleInfo = paste(month(startDate,label = TRUE,abbr = TRUE),
                            " '", substr(year(startDate), 3, 4), "\n",
                            depth, "m", sep = "")) %>%
  filter(year(startDate) %in% c(2020, 2021)) %>%
  mutate(treatment = fct_relevel(treatment,
                                  c("unfiltered-unamended",
                                    "unfiltered-molybdate")),
         sampleInfo = fct_relevel(sampleInfo, c("Sep '20\n11m", "Sep '20\n15.5m", "Sep '20\n20.7m",
                                                "Oct '20\n15.7m", "Oct '20\n20.9m",
                                                "Sep '21\n10.8m", "Sep '21\n11.9m", "Sep '21\n19.9m",
                                                "Oct '21\n14.2m", "Oct '21\n15.2m", "Oct '21\n19.9m")))

#### Condition vector ####
condition_naming <- c("t0 to t1",
                      "t1 to t2",
                      "Slope with all\ntime points")
names(condition_naming) <- c("Kmet_t1", "Kmet_t2", "Kmet_total")
color_vector <- cb.translator[c('bluishgreen', 'skyblue', 'vermillion')]
names(color_vector) <- names(condition_naming)


#### Look for differences between Kmet calculations for ambient samples ####
ambient.conditions <- Hg.Kmet.data %>%
  filter(treatment == "unfiltered-unamended") %>%
  select(startDate, depth, sampleInfo, Kmet_t1, Kmet_t2, Kmet_total) %>%
  gather(key = calculation,
         value = Kmet,
         4:6) %>%
  ggplot(aes(x = sampleInfo,
             y = Kmet,
             color = calculation,
             shape = sampleInfo,
             group = calculation)) +
  geom_point(position = position_dodge(width = 0.5)) +
  scale_color_manual(values = color_vector,
                     labels = condition_naming,
                     name = expression('K'['met']*' calculation method')) +
  scale_shape_manual(values = c(15, 1, 16, 2, 17, 3, 18, 4, 19, 5, 20)) +
  theme_bw() +
  theme(legend.position = 'none',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",
                                    linewidth = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(colour = "black"),
        axis.text.x =  element_text(colour = "black",
                                    angle = 45,
                                    hjust = 1),
        axis.text.y = element_text(colour = "black"),
        axis.ticks.x = element_line(colour = "black"),
        axis.ticks.y = element_line(colour = "black"),
        axis.line = element_blank()) +
  guides(shape = 'none') +
  ylim(c(-0.03, 0.25)) +
  ylab(expression('K'['met']*' '('day'^-1))) +
  labs(title = expression('K'['met']*' of ambient samples'))



#### Look for differences between Kmet calculations for molybdate-amended samples ####
molybdate.amended <- Hg.Kmet.data %>%
  filter(treatment == "unfiltered-molybdate") %>%
  select(startDate, depth, sampleInfo, Kmet_t1, Kmet_t2, Kmet_total) %>%
  gather(key = calculation,
         value = Kmet,
         4:6) %>%
  ggplot(aes(x = sampleInfo,
             y = Kmet,
             color = calculation,
             shape = sampleInfo,
             group = calculation)) +
  geom_point(position = position_dodge(width = 0.5)) +
  scale_color_manual(values = color_vector,
                     labels = condition_naming,
                     name = expression('K'['met']*' calculation method')) +
  scale_shape_manual(values = c(15, 1, 16, 2, 17, 3, 18, 4, 19, 5, 20)) +
  theme_bw() +
  theme(legend.position = c(0.8, 0.95),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",
                                    linewidth = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(colour = "black"),
        axis.text.x =  element_text(colour = "black",
                                    angle = 45,
                                    hjust = 1),
        axis.text.y = element_text(colour = "black"),
        axis.ticks.x = element_line(colour = "black"),
        axis.ticks.y = element_line(colour = "black"),
        axis.line = element_blank(),
        legend.background = element_rect(color = "black")) +
  guides(shape = 'none') +
  ylim(c(-0.03, 0.1)) +
  ylab(expression('K'['met']*' '('day'^-1))) +
  labs(title = expression('K'['met']*' of molybdate-inhibited samples'))


#### Save out figure ####
ggarrange(ambient.conditions, molybdate.amended,
          nrow = 2, ncol = 1)

pdf("results/figures/incubations_Kmet_comparison.pdf",
    width = 6,
    height = 7)
ggarrange(ambient.conditions, molybdate.amended,
          nrow = 2, ncol = 1)
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

