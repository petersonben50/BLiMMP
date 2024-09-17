#### code/figures/incubation_guilds_by_sulfide_maintext.R ####
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
hgcA_abundance_guilds <- hgcA_abundance %>%
  filter(verified_hgcA,
         clstr_rep == 1) %>%
  select(metabolic_assignment, all_of(grep("BLI2", names(hgcA_abundance)))) %>%
  gather(key = omicID, value = omic_coverage, -1) %>%
  left_join(omic_metadata) %>%
  group_by(seqType, omicID, date_depth, metabolic_assignment) %>%
  summarise(omic_coverage = sum(omic_coverage)) %>%
  ungroup() %>%
  group_by(date_depth, seqType, metabolic_assignment) %>%
  summarise(coverage_mean = mean(omic_coverage),
            coverage_sd = sd(omic_coverage),
            coverage_count = n(),
            coverage_se = coverage_sd / sqrt(coverage_count)) %>%
  select(date_depth, seqType, metabolic_assignment, coverage_mean, coverage_se) %>%
  ungroup()


#### Read in incubation rate data ####
Hg_Kmet_data <- read.csv("dataFinal/incubation_Hg_rate_data.csv") %>%
  dplyr::rename(date = startDate) %>%
  filter(treatment != "filtered-unamended",
         year(date) %in% c(2020, 2021)) %>%
  mutate(treatment = treatment %>%
           gsub("unfiltered-", "", .) %>%
           gsub("unamended", "ambient", .)) %>%
  select(sampleID, date, depth, treatment, Kmet_int_t2)  %>%
  group_by(sampleID, date, depth, treatment) %>%
  summarise(Kmet_mean = mean(Kmet_int_t2),
            Kmet_sd = sd(Kmet_int_t2),
            Kmet_count = n(),
            Kmet_se = Kmet_sd / Kmet_count) %>%
  ungroup() %>%
  gather(key = constituent,
         value = value,
         -c(1:4)) %>%
  filter(constituent %in% c("Kmet_mean", "Kmet_se")) %>%
  mutate(date_depth = paste(date, ":", depth, "m",
                            sep = ""),
         constituent = paste(treatment, constituent,
                             sep = "_")) %>%
  select(sampleID, date, date_depth, constituent, value) %>%
  spread(key = constituent,
         value = value) %>%
  filter(!is.na(ambient_Kmet_mean)) %>%
  mutate(SRB_Kmet_mean = ambient_Kmet_mean - molybdate_Kmet_mean,
         SRB_Kmet_se = ambient_Kmet_se + molybdate_Kmet_se) %>%
  dplyr::rename("nonSRB_Kmet_mean" = molybdate_Kmet_mean,
                "nonSRB_Kmet_se" = molybdate_Kmet_se,
                total_Kmet_mean = ambient_Kmet_mean,
                total_Kmet_se = ambient_Kmet_se) 


#### Read in sulfide data ####
sulfide_data <- read.csv("dataFinal/water_chem_data.csv") %>%
  group_by(date, depth) %>%
  summarize(sulfide_ppm = mean(sulfide_ppm, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(date_depth = paste(date, ":", depth, "m",
                            sep = "")) %>%
  select(date, date_depth, sulfide_ppm)


#### Set up vectors for aesthetics ####
color_vector <- c(cb.translator[c("black", "blue", "yellow", "bluishgreen")], "gray")
names(color_vector) <- c("KIR", "SRB", "RESP", "FERM", "UNK")
color_vector <- color_vector[length(color_vector):1]
name_vector <- c('SRB-dependent',
                 'SRB-independent')
year_vector_hgcA <- c(16, 17)
names(year_vector_hgcA) <- c("2020", "2021")
error_bar_width <- 0.025
# Sizes
point_size <- 1.3
legend_point_size <- 1.5


#### Refactor guild order ####
hgcA_abundance_guilds <- hgcA_abundance_guilds %>%
  mutate(metabolic_assignment = fct_relevel(metabolic_assignment,
                                            names(color_vector)))


#### Function to plot Kmet against sulfide ####
plot_Kmet_vs_sulfide_guilds <- function() {
  plot_data <- left_join(Hg_Kmet_data,
                         sulfide_data)
  plot(x = NULL,
       y = NULL,
       xlim = c(0, 5),
       ylim = c(0, 0.16),
       xlab = "",
       ylab = "")
  
  title(xlab = expression("Sulfide (mg/L)"),
        cex.lab = 1.2,
        line = 2)
  title(ylab = expression('K'['met']*' (day'^-1*')'),
        cex.lab = 1.2,
        line = 1.3)
  
  arrows(plot_data$sulfide_ppm, plot_data$nonSRB_Kmet_mean - plot_data$nonSRB_Kmet_se,
         plot_data$sulfide_ppm, plot_data$nonSRB_Kmet_mean + plot_data$nonSRB_Kmet_se,
         col = "gray",
         length = error_bar_width, angle = 90, code = 3)
  arrows(plot_data$sulfide_ppm, plot_data$SRB_Kmet_mean - plot_data$SRB_Kmet_se,
         plot_data$sulfide_ppm, plot_data$SRB_Kmet_mean + plot_data$SRB_Kmet_se,
         col = color_vector["SRB"],
         length = error_bar_width, angle = 90, code = 3)
  
  points(x = plot_data$sulfide_ppm,
         y = plot_data$nonSRB_Kmet_mean,
         pch = year_vector_hgcA[as.character(year(plot_data$date))],
         cex = point_size,
         col = "gray",
         lwd = 1)
  points(x = plot_data$sulfide_ppm,
         y = plot_data$SRB_Kmet_mean,
         pch = year_vector_hgcA[as.character(year(plot_data$date))],
         cex = point_size,
         col = color_vector["SRB"],
         lwd = 1)
  
  legend(x = -3,
         y = 0.165,
         legend = name_vector,
         pch = 16,
         pt.cex = legend_point_size,
         col = c(color_vector["SRB"], "gray"),
         bty = "n")
}

#### Set up function to generate colored points ####
omics_points_function <- function(plot_data, omicType = "MT", guild = "KIR") {
  column_name <- paste(omicType, guild, "coverage_mean", sep = "_")
  plot_data_temp <- plot_data[, c("date", "sulfide_ppm", column_name)] %>%
    filter((!is.na(.data[[column_name]]))) %>%
    as.data.frame()
  points(x = plot_data_temp[, "sulfide_ppm"],
         y = plot_data_temp[, column_name],
         pch = year_vector_hgcA[as.character(year(plot_data_temp$date))],
         col = color_vector[guild],
         cex = point_size)
}

omics_points_se_function <- function(plot_data, omicType = "MT", guild = "KIR") {
  mean_column_name <- paste(omicType, guild, "coverage_mean", sep = "_")
  se_column_name <- paste(omicType, guild, "coverage_se", sep = "_")
  plot_data_temp <- plot_data[, c("date", "sulfide_ppm", mean_column_name, se_column_name)] %>%
    filter(!is.na(.data[[mean_column_name]]),
           !is.na(.data[[se_column_name]])) %>%
    as.data.frame()
  
  arrows(plot_data_temp[, "sulfide_ppm"], plot_data_temp[, mean_column_name] - plot_data_temp[, se_column_name],
         plot_data_temp[, "sulfide_ppm"], plot_data_temp[, mean_column_name] + plot_data_temp[, se_column_name],
         length = error_bar_width, angle = 90, code = 3,
         col = color_vector[guild])
}


#### hgcA vs. sulfide ####
hgcA_vs_sulfide_plot <- function(ylabel, yscale, omicType) {
  
  plot_data <- left_join(hgcA_abundance_guilds,
                         sulfide_data) %>%
    mutate(identifier = paste(seqType, metabolic_assignment,
                              sep = "_")) %>%
    select(date_depth, date, sulfide_ppm, identifier, coverage_mean, coverage_se) %>%
    gather(key = metric,
           value = value,
           -c(1:4)) %>%
    mutate(ID_metric = paste(identifier, metric, sep = "_")) %>%
    select(date_depth, date, sulfide_ppm, ID_metric, value) %>%
    spread(key = ID_metric,
           value = value)
    
  plot(x = NULL,
       y = NULL,
       xlim = c(0, 5),
       ylim = yscale,
       xlab = "",
       ylab = "")
  
  sapply(names(color_vector),
         function(guild_to_plot) {
           omics_points_se_function(plot_data, omicType, guild_to_plot)
         })
  
  sapply(names(color_vector),
         function(guild_to_plot) {
           omics_points_function(plot_data, omicType, guild_to_plot)
         })
  
  title(xlab = expression("Sulfide (mg/L)"),
        cex.lab = 1.2,
        line = 2)
  title(ylab = ylabel,
        cex.lab = 1.2,
        line = 1.3)
  
}



#### Generate plots ####
cairo_pdf("results/figures/incubations_guilds_by_sulfide_maintext.pdf",
          family = "Arial",
          height = 3,
          width = 7.2)
par(mfrow = c(1, 3),
    mar = c(3, 3, 1.5, 1),
    tck = -0.008,
    mgp = c(1.5, 0.2, 0))
plot_Kmet_vs_sulfide_guilds()
mtext("A.", at = c(-30))
hgcA_vs_sulfide_plot(ylabel = expression(italic(hgcA)*' abundance (%)'),
                     yscale = c(0, 15),
                     omicType = "MG")
mtext("B.", at = c(-30))
legend(x = 5,
       y = 15,
       legend = names(color_vector)[length(color_vector):1],
       pch = 16,
       pt.cex = legend_point_size,
       col = color_vector[length(color_vector):1],
       bty = "n")

hgcA_vs_sulfide_plot(ylabel = expression(italic(hgcA)*' transcripts (10'^6*' per L)'),
                     yscale = c(0, 6),
                     omicType = "MT")
mtext("C.", at = c(-30))
legend(x = 5,
       y = 6,
       legend = names(color_vector)[length(color_vector):1],
       pch = 16,
       pt.cex = legend_point_size,
       col = color_vector[length(color_vector):1],
       bty = "n")
dev.off()


#### Values for Kmet from different guilds ####
left_join(Hg_Kmet_data,
          sulfide_data) %>%
  mutate(percent_inhibition_nonSRB = nonSRB_Kmet_mean / total_Kmet_mean * 100,
         percent_inhibition_SRB = SRB_Kmet_mean / total_Kmet_mean * 100) %>%
  summarise(min_kmet_nonSRB = min(nonSRB_Kmet_mean),
            max_kmet_nonSRB = max(nonSRB_Kmet_mean),
            mean_kmet_nonSRB = mean(nonSRB_Kmet_mean),
            sd_kmet_nonSRB = sd(nonSRB_Kmet_mean),
            count_kmet_nonSRB = n(),
            sem_kmet_nonSRB = sd_kmet_nonSRB / sqrt(count_kmet_nonSRB),
            min_percent_nonSRB = min(percent_inhibition_nonSRB),
            max_percent_nonSRB = max(percent_inhibition_nonSRB),
            mean_percent_nonSRB = mean(percent_inhibition_nonSRB),
            sd_percent_nonSRB = sd(percent_inhibition_nonSRB),
            count_percent_nonSRB = n(),
            sem_percent_nonSRB = sd_percent_nonSRB / sqrt(count_percent_nonSRB),
            min_kmet_SRB = min(SRB_Kmet_mean),
            max_kmet_SRB = max(SRB_Kmet_mean),
            mean_kmet_SRB = mean(SRB_Kmet_mean),
            sd_kmet_SRB = sd(SRB_Kmet_mean),
            count_kmet_SRB = n(),
            sem_kmet_SRB = sd_kmet_SRB / sqrt(count_kmet_SRB),
            min_percent_SRB = min(percent_inhibition_SRB),
            max_percent_SRB = max(percent_inhibition_SRB),
            mean_percent_SRB = mean(percent_inhibition_SRB),
            sd_percent_SRB = sd(percent_inhibition_SRB),
            count_percent_SRB = n(),
            sem_percent_SRB = sd_percent_SRB / sqrt(count_percent_SRB)) %>%
  select(min_kmet_nonSRB, max_kmet_nonSRB, mean_kmet_nonSRB, sem_kmet_nonSRB,
         min_percent_nonSRB, max_percent_nonSRB, mean_percent_nonSRB, sem_percent_nonSRB,
         min_kmet_SRB, max_kmet_SRB, mean_kmet_SRB, sem_kmet_SRB,
         min_percent_SRB, max_percent_SRB, mean_percent_SRB, sem_percent_SRB) %>%
  unlist()


#### Values for Kmet from different guilds when total Kmet is above 0.002 ####
left_join(Hg_Kmet_data,
          sulfide_data) %>%
  filter(total_Kmet_mean > 0.002) %>%
  mutate(percent_inhibition_nonSRB = nonSRB_Kmet_mean / total_Kmet_mean * 100,
         percent_inhibition_SRB = SRB_Kmet_mean / total_Kmet_mean * 100) %>%
  summarise(min_kmet_nonSRB = min(nonSRB_Kmet_mean),
            max_kmet_nonSRB = max(nonSRB_Kmet_mean),
            mean_kmet_nonSRB = mean(nonSRB_Kmet_mean),
            sd_kmet_nonSRB = sd(nonSRB_Kmet_mean),
            count_kmet_nonSRB = n(),
            sem_kmet_nonSRB = sd_kmet_nonSRB / sqrt(count_kmet_nonSRB),
            min_percent_nonSRB = min(percent_inhibition_nonSRB),
            max_percent_nonSRB = max(percent_inhibition_nonSRB),
            mean_percent_nonSRB = mean(percent_inhibition_nonSRB),
            sd_percent_nonSRB = sd(percent_inhibition_nonSRB),
            count_percent_nonSRB = n(),
            sem_percent_nonSRB = sd_percent_nonSRB / sqrt(count_percent_nonSRB),
            min_kmet_SRB = min(SRB_Kmet_mean),
            max_kmet_SRB = max(SRB_Kmet_mean),
            mean_kmet_SRB = mean(SRB_Kmet_mean),
            sd_kmet_SRB = sd(SRB_Kmet_mean),
            count_kmet_SRB = n(),
            sem_kmet_SRB = sd_kmet_SRB / sqrt(count_kmet_SRB),
            min_percent_SRB = min(percent_inhibition_SRB),
            max_percent_SRB = max(percent_inhibition_SRB),
            mean_percent_SRB = mean(percent_inhibition_SRB),
            sd_percent_SRB = sd(percent_inhibition_SRB),
            count_percent_SRB = n(),
            sem_percent_SRB = sd_percent_SRB / sqrt(count_percent_SRB)) %>%
  select(min_kmet_nonSRB, max_kmet_nonSRB, mean_kmet_nonSRB, sem_kmet_nonSRB,
         min_percent_nonSRB, max_percent_nonSRB, mean_percent_nonSRB, sem_percent_nonSRB,
         min_kmet_SRB, max_kmet_SRB, mean_kmet_SRB, sem_kmet_SRB,
         min_percent_SRB, max_percent_SRB, mean_percent_SRB, sem_percent_SRB) %>%
  unlist()

