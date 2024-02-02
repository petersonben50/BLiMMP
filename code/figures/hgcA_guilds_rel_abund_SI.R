#### code/figures/hgcA_guilds_rel_abund_SI.R ####
# Written for BLiMMP project
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP")
library(ggpubr)
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



#### hgcA gene statistics ####
hgcA_abundance_stats_data <- hgcA_abundance %>%
  filter(verified_hgcA,
         clstr_rep == 1)
# Counts
hgcA_abundance_stats_data %>%
  group_by(taxonomic_assignment) %>%
  summarise(count = n())

# Relative abundance metabolic groups
hgcA_abundance_stats_data_total <- hgcA_abundance_stats_data %>%
  select(all_of(grep("_M[G:T]_", colnames(hgcA_abundance_stats_data), value = TRUE))) %>%
  gather(key = omic_ID,
         value = coverage_total) %>%
  group_by(omic_ID) %>%
  summarise(coverage_total = sum(coverage_total))

hgcA_abundance_stats_data_RA <- hgcA_abundance_stats_data %>%
  select(metabolic_assignment,
         all_of(grep("_M[G:T]_", colnames(hgcA_abundance_stats_data), value = TRUE))) %>%
  gather(key = omic_ID,
         value = coverage,
         all_of(grep("_M[G:T]_", colnames(hgcA_abundance_stats_data), value = TRUE))) %>% 
  group_by(metabolic_assignment, omic_ID) %>%
  summarise(coverage = sum(coverage)) %>%
  left_join(hgcA_abundance_stats_data_total) %>%
  mutate(rel_cov = coverage / coverage_total * 100,
         omic_type = omic_ID %>%
           strsplit("_") %>% sapply("[", 2)) %>%
  filter(!is.na(rel_cov)) %>%
  group_by(metabolic_assignment, omic_type) %>%
  summarise(min_rel_cov = min(rel_cov),
            max_rel_cov = max(rel_cov),
            mean_rel_cov = mean(rel_cov),
            sd_rel_cov = sd(rel_cov),
            count_rel_cov = n(),
            sem_rel_cov = sd_rel_cov / sqrt(count_rel_cov)) %>%
  arrange(omic_type)
