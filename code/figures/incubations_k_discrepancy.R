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


#### Read in and prepare incubation data ####
incubation_data <- read.csv('dataFinal/incubation_Hg_rate_data.csv') %>%
  rename(date = startDate,
         Kmet = Kmet_int_t2,
         Kdem = Kdem_slope) %>%
  select(date, depth, treatment, Kmet, Kdem) %>%
  filter(!(date == "2020-10-10" & depth == 15.7))

incubation_data$Kdem[incubation_data$Kdem <= 0.005] <- 0.005
incubation_data$Kmet[incubation_data$Kmet <= 0.002] <- 0.002

ratio_data_rates <- incubation_data %>%
  filter(treatment == "unfiltered-unamended") %>%
  mutate(K_ratio = Kmet / Kdem) %>%
  group_by(date, depth) %>%
  summarize(K_ratio_mean = mean(K_ratio),
            K_ratio_sd = sd(K_ratio),
            K_ratio_n = n(),
            K_ratio_sem = K_ratio_sd / sqrt(K_ratio_n)) %>%
  select(date, depth, K_ratio_mean, K_ratio_sem) %>%
  ungroup()


#### Read in and prepare concentration data ####
ratio_data_conc <- read.csv("dataFinal/water_chem_data.csv") %>%
  filter(!is.na(FMHG_NG.L)) %>%
  group_by(date, depth) %>%
  summarise(FMHG_NG.L = mean(FMHG_NG.L),
            FiHg_NG.L = mean(FiHg_NG.L),
            sulfide_ppm = mean(sulfide_ppm)) %>%
  ungroup() %>%
  mutate(conc_ratio = FMHG_NG.L / FiHg_NG.L) %>%
  select(date, depth, sulfide_ppm, conc_ratio)


#### Combine ratio data ####
all_data <- left_join(ratio_data_rates,
                      ratio_data_conc)
rm(incubation_data,
   ratio_data_rates,
   ratio_data_conc)


#### Set up vectors ####
year_vector <- c(21, 24)
names(year_vector) <- c("2020", "2021")


#### Define color palette ####
# Define colour pallete
color_ramp = colorRampPalette(c("white", cb.translator["blue"]))


#### Function to plot data ####
plot_ratio_data <- function() {
  plot(y = log(all_data$conc_ratio, 10),
       x = log(all_data$K_ratio_mean, 10),
       xlim = c(-1.5, 1.5),
       ylim = c(-1.0, 1.0),
       pch = year_vector[substr(all_data$date, 1, 4)],
       col = "grey25",
       #bg = color_ramp(max(all_data$sulfide_ppm, na.rm = TRUE)*10)[all_data$sulfide_ppm*10],
       bg = "grey85",
       cex = 1.8,
       #cex = sqrt(all_data$sulfide_ppm),
       ylab = expression("Log of MeHg / Hg(II) ratio"),
       xlab = expression('Log of K'[italic('met')]*' /  K'[italic('dem')]*' ratio'))
  text(x = log(all_data$K_ratio_mean, 10),
       y = log(all_data$conc_ratio, 10),
       labels = paste(all_data$date, all_data$depth, sep = ", "),
       cex = 0.8)
  arrows(x0 = log((all_data$K_ratio_mean - all_data$K_ratio_sem), 10),
         x1 = log((all_data$K_ratio_mean + all_data$K_ratio_sem), 10),
         y0 = log(all_data$conc_ratio, 10), y1 = log(all_data$conc_ratio, 10),
         angle = 90, code = 3, length = 0.05)
  abline(a = 0, b = 1)
  legend(x = -1.35, y = 0.9,
         legend = names(year_vector),
         pch = year_vector,
         col = "grey25",
         pt.bg = "grey85")
}



#### Generate plot ####
# cairo_pdf("results/figures/incubations_K_discrepancy.pdf",
#           family = "Arial",
#           height = 4,
#           width = 4)
par(mfrow = c(1, 1),
    mar = c(3.5, 3.5, 1, 1),
    tck = -0.008,
    mgp = c(2, 0.4, 0),
    cex.axis = 1.1,
    cex.lab = 1.2)
plot_ratio_data()
# dev.off()
