#### code/figures/incubations_kdem_SI.R ####
# Written for BLiMMP project
# SI figure for manuscript
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP")
library(extrafont)
# library(ggpubr)
# library(lubridate)
# library(readxl)
library(tidyverse)
cb_translator <- readRDS("references/colorblind_friendly_colors.rds")


#### Read in Kmet data ####
Kdem_data <- read.csv("dataFinal/incubation_Hg_rate_data.csv") %>%
  mutate(date_depth = paste(startDate, ":", depth, "m",
                            sep = "")) %>%
  filter(date_depth != "2020-10-10:15.7m") %>%
  select(date_depth, treatment, Kdem_slope)


#### Summarize Kmet data ####
Kdem_data_summarized <- Kdem_data %>%
  group_by(date_depth, treatment) %>%
  summarize(Kdem_slope_mean = mean(Kdem_slope, na.rm = TRUE),
            Kdem_slope_sd = sd(Kdem_slope, na.rm = TRUE),
            Kdem_count = n(),
            Kdem_slope_sem = Kdem_slope_sd / sqrt(Kdem_count))

#### Read in sulfide data ####
sulfide_data <- read.csv("dataFinal/water_chem_data.csv") %>%
  group_by(date, depth) %>%
  summarize(sulfide_uM = mean(sulfide_uM, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(date_depth = paste(date, ":", depth, "m",
                            sep = "")) %>%
  select(date_depth, sulfide_uM)

#### Set up aesthetic vectors ####
naming_vector_kdem <- c("Ambient", "Molybdate-treated", "Filtered control")
names(naming_vector_kdem) <- c("unfiltered-unamended", "unfiltered-molybdate", "filtered-unamended")
color_vector_kdem <- cb_translator[c("bluishgreen", "black", "orange")]
names(color_vector_kdem) <- names(naming_vector_kdem)

date_depth_order_by_sulfide <- sulfide_data %>%
  filter(date_depth %in% Kdem_data_summarized$date_depth) %>%
  arrange(sulfide_uM) %>%
  select(date_depth) %>%
  unlist()
names(date_depth_order_by_sulfide) <- NULL


#### Set up graph ####
figure_plot <- Kdem_data_summarized %>%
  mutate(treatment = fct_relevel(treatment, names(naming_vector_kdem)),
         date_depth = as.factor(date_depth),
         date_depth = fct_relevel(date_depth, date_depth_order_by_sulfide)) %>%
  ggplot(aes(x = date_depth,
             y = Kdem_slope_mean,
             fill = treatment)) + 
  geom_bar(stat = "identity",
           color = "black", 
           position = position_dodge()) +
  scale_fill_manual(values = color_vector_kdem,
                    labels = naming_vector_kdem,
                    name = expression(italic(Treatment))) +
  geom_errorbar(aes(ymin = Kdem_slope_mean - Kdem_slope_sem,
                    ymax = Kdem_slope_mean + Kdem_slope_sem),
                width = 0.2,
                position = position_dodge(0.9)) +
  geom_hline(yintercept = 0) +
  scale_y_continuous(limits = c(-0.15, 0.5),
                     breaks = seq(-0.1, 0.5, by = 0.1)) +
  labs(x = element_blank(),
       y = expression('K'['dem'])) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =  90),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
figure_plot

#### Save out figure ####
cairo_pdf("results/figures/incubations_Kdem.pdf",
          family = "Arial",
          height = 4.5,
          width = 7.2)
figure_plot
dev.off()
