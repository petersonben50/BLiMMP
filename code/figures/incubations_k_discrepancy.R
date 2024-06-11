#### code/figures/incubations_kmet_maintext.R ####
# Written for BLiMMP project
# Supplementary figure looking at the discrepancy between
# K ratios and MeHg/Hg(II)
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP")
library(tidyverse)
cb.translator <- readRDS("references/colorblind_friendly_colors.rds")


#### Read in and prepare data ####
incubation_data <- read.csv('dataFinal/incubation_Hg_rate_data.csv') %>%
  rename(date = startDate,
         Kmet = Kmet_int_t2,
         Kdem = Kdem_slope) %>%
  select(date, depth, treatment, Kmet, Kdem) %>%
  filter(!(date == "2020-10-10" & depth == 15.7)) %>%
  group_by(date, depth, treatment) %>%
  summarise(counts = n(),
            Kmet_mean = mean(Kmet),
            Kmet_sem = sd(Kmet)/sqrt(counts),
            Kdem_mean = mean(Kdem),
            Kdem_sem = sd(Kdem)/sqrt(counts)) %>%
  ungroup() %>%
  select(-counts)
geochem_data <- read.csv("dataFinal/water_chem_data.csv") %>%
  filter(!is.na(sulfate_uM),
         !corewater) %>%
  group_by(date, depth) %>%
  summarise(FMHG_NG.L = mean(FMHG_NG.L),
            FiHg_NG.L = mean(FiHg_NG.L)) %>%
  ungroup()
all_data <- left_join(incubation_data,
                      geochem_data)
rm(geochem_data, incubation_data)


#### Calculate ratios ####
ratio_data <- all_data %>%
  filter(treatment == "unfiltered-unamended") 
ratio_data$Kdem_mean[ratio_data$Kdem_mean < 0.01] <- 0.01
ratio_data <- ratio_data %>%
  mutate(conc_ratio = FMHG_NG.L / FiHg_NG.L,
         K_ratio = Kmet_mean / Kdem_mean,
         K_discrepancy = K_ratio / conc_ratio) %>%
  select(date, depth, K_discrepancy, conc_ratio, K_ratio, FMHG_NG.L, FiHg_NG.L, Kmet_mean, Kdem_mean)


#### Generate plot ####
cairo_pdf("results/figures/incubations_K_discrepancy.pdf",
          family = "Arial",
          height = 4,
          width = 4)
par(mfrow = c(1, 1),
    mar = c(3.5, 3.5, 1, 1),
    tck = -0.008,
    mgp = c(2, 0.4, 0),
    cex.axis = 1.1,
    cex.lab = 1.2)
plot(x = ratio_data$Kmet_mean,
     y = log(ratio_data$K_discrepancy, 10),
     xlim = c(0, 0.2),
     ylim = c(-1.2, 1.2),
     xlab = expression('K'['met']*' (day'^-1*')'),
     ylab = expression('(K'['met']*' /  K'['dem']*') / (MeHg / Hg(II))'),
     pch = 16,
     cex = 1.2)
dev.off()
