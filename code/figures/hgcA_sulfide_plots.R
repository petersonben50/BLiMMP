#### code/figures/incubations_hgcA_sulfide.R ####
# Written for BLiMMP project
# Supplemental figure of manuscript
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP")
library(extrafont)
library(ggpubr)
library(tidyverse)
cb.translator <- readRDS("references/colorblind_friendly_colors.rds")


#### Set constants ####
label_size <- 1.25

year_vector <- c(21, 24)
names(year_vector) <- c("2020", "2021")


#### Prepare omic data ####
MG_metadata <- read.csv("metadata/metagenome_metadata.csv") %>%
  dplyr::rename(omicID = metagenomeID) %>%
  mutate(seqType = "MG") %>%
  select(omicID, sampleID, startDate, depth, seqType)
MT_metadata <- read.csv("metadata/metatranscriptome_metadata.csv") %>%
  dplyr::rename(omicID = metatranscriptomeID) %>%
  mutate(seqType = "MT") %>%
  select(omicID, sampleID, startDate, depth, seqType)
omic_metadata <- rbind(MG_metadata,
                       MT_metadata) %>%
  mutate(date_depth = paste(startDate, ":", depth, "m",
                            sep = ""))
rm(MG_metadata, MT_metadata)


#### Prepare hgcA data ####
hgcA_data <- read.csv("dataFinal/hgcA_data.csv")
# There were no MT hits to hgcA in BLI21_MT_003, which excluded it from the table.
hgcA_data$BLI21_MT_003 <- 0
hgcA_data_total <- hgcA_data %>%
  filter(verified_hgcA,
         clstr_rep == 1) %>%
  select(metabolic_assignment, all_of(grep("BLI2", names(hgcA_data)))) %>%
  gather(key = omicID, value = omic_coverage, -1) %>%
  left_join(omic_metadata) %>%
  group_by(seqType, omicID, date_depth) %>%
  summarise(omic_coverage = sum(omic_coverage)) %>%
  ungroup() %>%
  group_by(date_depth, seqType) %>%
  summarise(coverage_mean = mean(omic_coverage),
            coverage_sd = sd(omic_coverage),
            coverage_count = n(),
            coverage_se = coverage_sd / sqrt(coverage_count))



#### Read in sulfide data ####
sulfide_data <- read.csv("dataFinal/water_chem_data.csv") %>%
  group_by(date, depth) %>%
  summarize(sulfide_ppm = mean(sulfide_ppm, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(date_depth = paste(date, ":", depth, "m",
                            sep = ""))


#### All data ####
all_data <- left_join(hgcA_data_total, sulfide_data)



#### Function to plot hgcA vs. sulfide ####
hgcA_vs_sulfide <- function(seq_type,
                            xscale = c(0, 5),
                            xlabel = "Sulfide (mg/L)",
                            ylabel,
                            yscale) {
  plot_data <- all_data %>%
    filter(seqType == seq_type)
  
  plot(x = NULL,
       y = NULL,
       xlim = xscale,
       ylim = yscale,
       ylab = '',
       xlab = '')
  arrows(plot_data$sulfide_ppm, plot_data$coverage_mean - plot_data$coverage_se,
         plot_data$sulfide_ppm, plot_data$coverage_mean + plot_data$coverage_se,
         length = 0.05, angle = 90, code = 3)
  points(x = plot_data$sulfide_ppm,
         y = plot_data$coverage_mean,
         pch = year_vector[as.character(year(plot_data$date))],
         col = "gray25",
         bg = "gray85",
         lwd = 1.5,
         cex = 2)
  
  title(xlab = xlabel,
        cex.lab = label_size,
        line = 1.5)
  title(ylab = ylabel,
        cex.lab = label_size,
        line = 1.5)
  
}


#### Generate plots for manuscript main text figure ####
cairo_pdf("results/figures/hgcA_sulfide_plots.pdf",
          family = "Arial",
          height = 4,
          width = 7.2)
par(mfrow = c(1, 2),
    mar = c(3, 3.5, 1.5, 1),
    tck = -0.008,
    mgp = c(1.5, 0.2, 0),
    cex.axis = 1.1)
hgcA_vs_sulfide(seq_type = "MG",
                ylabel = expression(italic(hgcA)*' abundance (%)'),
                yscale = c(0, 17))
mtext("A.", at = c(-0.8))
hgcA_vs_sulfide(seq_type = "MT",
                ylabel = expression(italic(hgcA)*' transcripts (10'^6*' per L)'),
                yscale = c(0, 11))
mtext("B.", at = c(-0.8))
dev.off()


#### Linear regression for Kmet and hgcA abundance ####
regression_data <- all_data %>%
  filter(seqType == "MG")
regression_data_lm <- lm(coverage_mean ~ sulfide_ppm, data = regression_data)
summary(regression_data_lm)

#### Linear regression for Kmet and hgcA transcripts ####
regression_data <- all_data %>%
  filter(seqType == "MT")
regression_data_lm <- lm(coverage_mean ~ sulfide_ppm, data = regression_data)
summary(regression_data_lm)

