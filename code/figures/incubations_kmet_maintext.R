#### code/figures/incubations_kmet_maintext.R ####
# Written for BLiMMP project
# Figure 2 of manuscript
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP")
library(extrafont)
library(ggpubr)
# library(lubridate)
# library(readxl)
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
hgcA_abundance <- read.csv("dataFinal/hgcA_data.csv")
# There were no MT hits to hgcA in BLI21_MT_003, which excluded it from the table.
hgcA_abundance$BLI21_MT_003 <- 0
hgcA_abundance_total <- hgcA_abundance %>%
  filter(verified_hgcA,
         clstr_rep == 1) %>%
  select(metabolic_assignment, all_of(grep("BLI2", names(hgcA_abundance)))) %>%
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


#### Read in incubation rate data ####
Hg_Kmet_data <- read.csv("dataFinal/incubation_Hg_rate_data.csv") %>%
  dplyr::rename(date = startDate) %>%
  filter(treatment == "unfiltered-unamended",
         year(date) %in% c(2020, 2021)) %>%
  mutate(treatment = treatment %>%
           gsub("unfiltered-", "", .) %>%
           gsub("unamended", "ambient", .)) %>%
  select(sampleID, date, depth, treatment, Kmet_int_t2)  %>%
  group_by(sampleID, date, depth, treatment) %>%
  summarise(Kmet_mean = median(Kmet_int_t2),
            Kmet_sd = sd(Kmet_int_t2),
            Kmet_count = n(),
            Kmet_se = Kmet_sd / Kmet_count) %>%
  ungroup() %>%
  gather(key = constituent,
         value = value,
         -c(1:4)) %>%
  filter(constituent %in% c("Kmet_mean", "Kmet_se")) %>%
  mutate(date_depth = paste(date, ":", depth, "m",
                            sep = "")) %>%
  select(sampleID, date, date_depth, constituent, value) %>%
  spread(key = constituent,
         value = value) %>%
  filter(!is.na(Kmet_mean))


#### Read in sulfide data ####
sulfide_data <- read.csv("dataFinal/water_chem_data.csv") %>%
  group_by(date, depth) %>%
  summarize(sulfide_uM = mean(sulfide_uM, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(date_depth = paste(date, ":", depth, "m",
                            sep = ""))


#### Function to plot Kmet vs. sulfide ####
Kmet_vs_sulfide <- function() {
  plot_data <- Hg_Kmet_data %>%
    left_join(sulfide_data)
  plot(x = NULL,
       y = NULL,
       xlim = c(0, 150),
       ylim = c(0, 0.2),
       xlab = "",
       ylab = "")
  title(xlab = expression("Sulfide (µM)"),
        cex.lab = label_size,
        line = 2)
  title(ylab = expression('K'['met']*' (day'^-1*')'),
        cex.lab = label_size,
        line = 1.5)
  arrows(plot_data$sulfide_uM, plot_data$Kmet_mean - plot_data$Kmet_se,
         plot_data$sulfide_uM, plot_data$Kmet_mean + plot_data$Kmet_se,
         length = 0.05, angle = 90, code = 3)
  points(x = plot_data$sulfide_uM,
         y = plot_data$Kmet_mean,
         pch = year_vector[as.character(year(plot_data$date))],
         col = "gray25",
         bg = "gray85",
         lwd = 1.5,
         cex = 2)
}
Kmet_vs_hgcA <- function(seq_type, xlabel,
                         xscale) {
  plot_data <- left_join(hgcA_abundance_total %>%
                           filter(seqType == seq_type),
                         Hg_Kmet_data) %>%
    filter(!is.na(Kmet_mean))
  plot(x = NULL,
       y = NULL,
       xlim = xscale,
       ylim = c(0, 0.2),
       ylab = '',
       xlab = '')
  arrows(plot_data$coverage_mean, plot_data$Kmet_mean - plot_data$Kmet_se,
         plot_data$coverage_mean, plot_data$Kmet_mean + plot_data$Kmet_se,
         length = 0.05, angle = 90, code = 3)
  arrows(plot_data$coverage_mean - plot_data$coverage_se, plot_data$Kmet_mean,
         plot_data$coverage_mean + plot_data$coverage_se, plot_data$Kmet_mean,
         length = 0.05, angle = 90, code = 3)
  points(x = plot_data$coverage_mean,
         y = plot_data$Kmet_mean,
         pch = year_vector[as.character(year(plot_data$date))],
         col = "gray25",
         bg = "gray85",
         lwd = 1.5,
         cex = 2)
  
  title(xlab = xlabel,
        cex.lab = label_size,
        line = 2)
  title(ylab = expression('K'['met']*' (day'^-1*')'),
        cex.lab = label_size,
        line = 1.5)
  
}


#### Generate plots for manuscript main text figure ####
cairo_pdf("results/figures/incubations_kmet_maintext.pdf",
    family = "Arial",
    height = 3,
    width = 7.2)
par(mfrow = c(1, 3),
    mar = c(3, 3.5, 1.5, 1),
    tck = -0.008,
    mgp = c(1.5, 0.2, 0),
    cex.axis = 1.1)
Kmet_vs_sulfide()
mtext("A.", at = c(-30))
legend(x = 3,
       y = 0.20,
       legend = names(year_vector),
       pch = year_vector,
       col = "gray25",
       pt.bg = "gray85",
       cex = 1)
Kmet_vs_hgcA(seq_type = "MG",
             xlabel = expression(italic(hgcA)*' abundance (%)'),
             xscale = c(0, 17))
mtext("B.", at = c(-3.4))
Kmet_vs_hgcA(seq_type = "MT",
             xlabel = expression(italic(hgcA)*' transcripts (10'^6*' per L)'),
             xscale = c(0, 11))
mtext("C.", at = c(-2.2))
dev.off()
# embed_fonts("results/figures/incubations_kmet_maintext.pdf")


#### Generate plots for manuscript supplemental methods figure, part 1 ####
cairo_pdf("results/figures/methods_figure_hgcA_MG_MT_1.pdf",
          family = "Arial",
          height = 2.5,
          width = 7.5)
par(mfrow = c(1, 2),
    mar = c(1, 3.5, 1.5, 1),
    tck = -0.008,
    mgp = c(1.5, 0.2, 0),
    cex.axis = 1.1)
Kmet_vs_hgcA(seq_type = "MG",
             xlabel = "",
             xscale = c(0, 17))
Kmet_vs_hgcA(seq_type = "MT",
             xlabel = "",
             xscale = c(0, 11))
dev.off()


#### Generate plots for manuscript supplemental methods figure, part 2 ####
cairo_pdf("results/figures/methods_figure_hgcA_MG_MT_2.pdf",
          family = "Arial",
          height = 2.5,
          width = 4)
par(mfrow = c(1, 1),
    mar = c(1, 3.5, 1.5, 1),
    tck = -0.008,
    mgp = c(1.5, 0.2, 0),
    cex.axis = 1.1)
Kmet_vs_sulfide()
dev.off()
