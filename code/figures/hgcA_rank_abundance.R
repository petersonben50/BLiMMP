#### code/figures/hgcA_rank_abundance.R ####
# Written for BLiMMP project
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP")
library(ggpubr)
library(readxl)
library(tidyverse)
cb.translator <- readRDS("references/colorblind_friendly_colors.rds")


#### Read in sulfide data ####
sulfide_data <- read.csv("dataFinal/water_chem_data.csv") %>%
  group_by(date, depth) %>%
  summarize(sulfide_uM = mean(sulfide_uM, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(date_depth = paste(date, ":", depth, "m",
                            sep = "")) %>%
  select(date_depth, sulfide_uM)


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
                            sep = "")) %>%
  left_join(sulfide_data) %>%
  filter(date_depth != "2020-10-10:15.7m") %>%
  arrange(sulfide_uM) %>%
  filter(sulfide_uM > 30) # This effectively cuts out the date-depths where we didn't do RNA sequencing.
rm(MG_metadata, MT_metadata, sulfide_data)


#### Set up colors ####
color_vector <- c(cb.translator[c("black", "blue", "yellow", "bluishgreen")], "gray")
names(color_vector) <- c("KIR", "SRB", "RESP", "FERM", "UNK")


#### Prepare hgcA data ####
hgcA_abundance <- read.csv("dataFinal/hgcA_data.csv")
# There were no MT hits to hgcA in BLI21_MT_003, which excluded it from the table.
hgcA_abundance$BLI21_MT_003 <- 0
hgcA_abundance <- hgcA_abundance %>%
  filter(verified_hgcA,
         clstr_rep == 1) %>%
  select(cluster_ID, taxonomic_assignment, metabolic_assignment, all_of(grep("BLI2", names(hgcA_abundance)))) %>%
  gather(key = omicID, value = omic_coverage, -c(1:3)) %>%
  left_join(omic_metadata)


#### Prepare plotting function ####
plot_rank_abund <- function(date_depth_to_use = "2020-09-02:20.7m",
                            seqType_to_use = "MG",
                            number_seqs_to_include = 20,
                            ymax_to_use = NULL) {
  plot_data <- hgcA_abundance %>%
    filter(date_depth == date_depth_to_use,
           seqType == seqType_to_use) 
  if (dim(plot_data)[1] > 0) {
    plot_data <- plot_data %>%
      group_by(cluster_ID, metabolic_assignment, date_depth) %>%
      summarise(cov_mean = mean(omic_coverage),
                cov_sd = sd(omic_coverage),
                cov_count = n(),
                cov_sem = cov_sd / sqrt(cov_count)) %>%
      arrange(desc(cov_mean))
    plot_data <- plot_data[c(1:number_seqs_to_include), ] %>%
      arrange(desc(cov_mean))
    
    # Set y max if not already:
    if (is.null(ymax_to_use)) {
      ymax_to_use <- max(plot_data$cov_mean)*1.5
    }
    
    # Set up vector for hgcA IDs
    IDs <- gsub(pattern = "BLI_hgcA_clstr_0",
                replacement = "",
                plot_data$cluster_ID)

    if (seqType_to_use == "MG" & substr(date_depth_to_use, 1, 4) == "2020") {
      barplot(plot_data$cov_mean,
              col = color_vector[plot_data$metabolic_assignment],
              ylim = c(0, ymax_to_use))
      
      mid_points <- barplot(plot_data$cov_mean[1:number_seqs_to_include],
                            plot = FALSE)
      text(x = mid_points,
           y = plot_data$cov_mean + 0.1,
           labels = IDs,
           cex = 0.75,
           adj = c(0.5, 0))
    } else {
      barplot(plot_data$cov_mean,
              col = color_vector[plot_data$metabolic_assignment],
              ylim = c(0, ymax_to_use))
      mid_points <- barplot(plot_data$cov_mean[1:number_seqs_to_include],
                            plot = FALSE)
      arrows(mid_points, plot_data$cov_mean + plot_data$cov_sem,
             mid_points, plot_data$cov_mean,
             length = 0.05, angle = 90, code = 3)
      text(x = mid_points,
           y = plot_data$cov_mean + plot_data$cov_sem + 0.1,
           labels = IDs,
           cex = 0.75,
           adj = c(0.5, 0))
    }
  } else {
    plot.new()
  }
}


### 
par(mfrow = c(length(unique(omic_metadata$date_depth)), 2),
    mar = c(1, 3, 0, 0))

for (date_depth_of_interest in unique(omic_metadata$date_depth)) {
  plot_rank_abund(date_depth_to_use = date_depth_of_interest,
                  seqType_to_use = "MG",
                  ymax_to_use = 4.5)  
  plot_rank_abund(date_depth_to_use = date_depth_of_interest,
                  seqType_to_use = "MT",
                  ymax_to_use = 1.5)
}
