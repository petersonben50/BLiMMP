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


#### Read in sulfide data ####
sulfide_data <- read.csv("dataFinal/water_chem_data.csv") %>%
  group_by(date, depth) %>%
  summarize(sulfide_uM = mean(sulfide_uM, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(date_depth = paste(date, ":", depth, "m",
                            sep = "")) %>%
  select(date, date_depth, sulfide_uM)


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
  gather(key = omicID,
         value = coverage_total) %>%
  group_by(omicID) %>%
  summarise(coverage_total = sum(coverage_total))

hgcA_abundance_data_RA <- hgcA_abundance_stats_data %>%
  select(metabolic_assignment,
         all_of(grep("_M[G:T]_", colnames(hgcA_abundance_stats_data), value = TRUE))) %>%
  gather(key = omicID,
         value = coverage,
         all_of(grep("_M[G:T]_", colnames(hgcA_abundance_stats_data), value = TRUE))) %>% 
  group_by(metabolic_assignment, omicID) %>%
  summarise(coverage = sum(coverage)) %>%
  left_join(hgcA_abundance_stats_data_total) %>%
  mutate(rel_cov = coverage / coverage_total * 100,
         omic_type = omicID %>%
           strsplit("_") %>% sapply("[", 2)) 

hgcA_abundance_data_RA %>%
  filter(!is.na(rel_cov)) %>%
  group_by(metabolic_assignment, omic_type) %>%
  summarise(min_rel_cov = min(rel_cov),
            max_rel_cov = max(rel_cov),
            mean_rel_cov = mean(rel_cov),
            sd_rel_cov = sd(rel_cov),
            count_rel_cov = n(),
            sem_rel_cov = sd_rel_cov / sqrt(count_rel_cov)) %>%
  arrange(omic_type)


#### Read in sulfide data ####
sulfide_data <- read.csv("dataFinal/water_chem_data.csv") %>%
  group_by(date, depth) %>%
  summarize(sulfide_uM = mean(sulfide_uM, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(date_depth = paste(date, ":", depth, "m",
                            sep = "")) %>%
  select(date_depth, sulfide_uM)


#### Read in metadata ####
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
  mutate(groupID = data.table::rowid(date_depth, seqType,
                                     prefix = "rep")) %>%
  arrange(date_depth)
rm(MG_metadata, MT_metadata)


#### Set up colors ####
color_vector <- c(cb.translator[c("black", "blue", "yellow", "bluishgreen")], "gray")
names(color_vector) <- c("KIR", "SRB", "RESP", "FERM", "UNK")


#### Prep plotting data ####
plot_data <- hgcA_abundance_data_RA %>%
  left_join(omic_metadata) %>%
  left_join(sulfide_data) %>%
  arrange(sulfide_uM) %>%
  select(metabolic_assignment, omic_type, date_depth, groupID, coverage, rel_cov, sulfide_uM) %>%
  mutate(date_depth = gsub(":", "\n", date_depth))
plot_data <- plot_data %>%
  mutate(date_depth = factor(date_depth,
                             levels = unique(plot_data$date_depth)),
         metabolic_assignment = factor(metabolic_assignment,
                                       levels = names(color_vector)))


#### Plotting function: Guilds ####
OMIC_TYPE = "MT"
DATE_DEPTH = "2021-10-14\n14.2m"
plot_omic_guilds <- function(OMIC_TYPE = "MG", DATE_DEPTH = "2021-10-14\n14.2m") {
  temp_data <- plot_data %>%
    filter(omic_type == OMIC_TYPE,
           date_depth == DATE_DEPTH) %>%
    select(metabolic_assignment, groupID, rel_cov) %>%
    spread(key = metabolic_assignment,
           value = rel_cov) %>%
    column_to_rownames(var = "groupID") %>%
    as.matrix()
  if (dim(temp_data)[1] == 0) {
    empty_plot_function()
  } else {
    color_vector_to_use <- rep(color_vector,
                               each = dim(temp_data)[1])
    spacing_vector = c(0, 0.1)
    barplot(temp_data,
            beside = TRUE,
            axes = FALSE,
            axisnames = FALSE,
            ylim = c(0, 100),
            col = color_vector_to_use,
            space = spacing_vector)
    
  }
}


#### Plotting function: Totals ####
OMIC_TYPE = "MT"
DATE_DEPTH = "2021-10-14\n19.9m"
plot_omics_totals <- function(OMIC_TYPE = "MG", DATE_DEPTH = "2021-10-14\n14.2m") {
  temp_data <- plot_data %>%
    filter(omic_type == OMIC_TYPE,
           date_depth == DATE_DEPTH) %>%
    group_by(groupID) %>%
    summarise(coverage = sum(coverage)) %>%
    column_to_rownames(var = "groupID") %>%
    as.matrix()
  if (OMIC_TYPE == "MT") {
    max_y = MT_totals_max
  } else if (OMIC_TYPE == "MG") {
    max_y = MG_totals_max
  }
  if (dim(temp_data)[1] == 0) {
    empty_plot_function()
  } else {
    color_vector_to_use <- rep(color_vector,
                               each = dim(temp_data)[1])
    spacing_vector = c(0, 0.1)
    barplot(temp_data,
            beside = TRUE,
            axes = FALSE,
            axisnames = FALSE,
            ylim = c(0, max_y),
            col = "gray40",
            space = spacing_vector)
    
  }
}


#### Plotting function: Sulfide ####
DATE_DEPTH = "2021-10-14\n19.9m"
plot_sulfide <- function(DATE_DEPTH = "2021-10-14\n14.2m") {
  temp_data <- sulfide_data %>%
    mutate(date_depth = gsub(":", "\n", date_depth)) %>%
    filter(date_depth == DATE_DEPTH) %>%
    select(sulfide_uM) %>%
    as.matrix()
  barplot(temp_data,
            beside = TRUE,
            axes = FALSE,
            axisnames = FALSE,
            ylim = c(0, sulfide_max),
            col = cb.translator["blue"])
    }


#### Plotting function: Empty plot ####
empty_plot_function <- function() {
  barplot(NA,
          axes = FALSE,
          axisnames = FALSE,
          ylim = c(0, 1))
}


#### Set up layout for figure ####
vert_plot_space_guilds <- 0.27
vert_plot_space_summaries <- 0.10
vert_spacer <- 0.02
title_space <- 0.05

header_vertArea <- c(1 - title_space, 1)
row_1_sulfide_vertArea <- c(header_vertArea[1] - vert_plot_space_summaries - vert_spacer,
                            header_vertArea[1] - vert_spacer)
row_2_MG_vertArea <- c(row_1_sulfide_vertArea[1] - vert_plot_space_summaries - vert_spacer,
                       row_1_sulfide_vertArea[1] - vert_spacer)
row_3_MG_guilds_vertArea <- c(row_2_MG_vertArea[1] - vert_plot_space_guilds - vert_spacer,
                              row_2_MG_vertArea[1] - vert_spacer)
row_4_MT_vertArea <- c(row_3_MG_guilds_vertArea[1] - vert_plot_space_summaries - vert_spacer,
                       row_3_MG_guilds_vertArea[1] - vert_spacer)
row_5_MT_guilds_vertArea <- c(row_4_MT_vertArea[1] - vert_plot_space_guilds - vert_spacer,
                              row_4_MT_vertArea[1] - vert_spacer)
vertical_areas <- list(header_vertArea, row_1_sulfide_vertArea, row_2_MG_vertArea, row_3_MG_guilds_vertArea, row_4_MT_vertArea, row_5_MT_guilds_vertArea)

# Horizontal areas
increment <- (0.081)
horiArea1_lab <- c(0.00,0.099)
horiArea2 <- c(horiArea1_lab[2], horiArea1_lab[2] + increment)
horiArea3 <- c(horiArea2[2], horiArea2[2] + increment)
horiArea4 <- c(horiArea3[2], horiArea3[2] + increment)
horiArea5 <- c(horiArea4[2], horiArea4[2] + increment)
horiArea6 <- c(horiArea5[2], horiArea5[2] + increment)
horiArea7 <- c(horiArea6[2], horiArea6[2] + increment)
horiArea8 <- c(horiArea7[2], horiArea7[2] + increment)
horiArea9 <- c(horiArea8[2], horiArea8[2] + increment)
horiArea10 <- c(horiArea9[2], horiArea9[2] + increment)
horiArea11 <- c(horiArea10[2], horiArea10[2] + increment)
horiArea12 <- c(horiArea11[2], horiArea11[2] + increment)
horiAreas <- list(horiArea1_lab, horiArea2, horiArea3, horiArea4, horiArea5, horiArea6, horiArea7, horiArea8, horiArea9, horiArea10, horiArea11, horiArea12)
# Set up rbind
screen_layout = NULL
for (vertical_area in vertical_areas) {
  for (horiArea in horiAreas) {
    if (is.null(screen_layout)) {
      screen_layout = c(horiArea, vertical_area)
    } else {
      screen_layout = rbind(screen_layout, c(horiArea, vertical_area))
    }
    
  }
}
screen_layout <- rbind(screen_layout, c(horiArea2[1], horiArea3[2],
                                        row_5_MT_guilds_vertArea))

#### Initializing plot ####
grDevices::cairo_pdf("results/figures/hgcA_guilds_RA_SI.pdf",
                     family = "Arial",
                     width = 7.20,
                     height = 7)


split.screen(screen_layout)
margins_to_use <- c(0, 0.5, 0, 0)

#### Plot: headers ####
i <- 2
for (date_depth_label in unique(plot_data$date_depth)) {
  screen(i)
  par(mar = margins_to_use)
  empty_plot_function()
  text(x = 0.7, y = 0.35, adj = 0.5,
       labels = date_depth_label,
       cex = 0.5)
  i <- i + 1
}


#### Plot: sulfide ####
# Sulfide y-axis
screen(i)
sulfide_max <- 150
par(mar = margins_to_use)
empty_plot_function()
axis(side = 4,
     at = seq(0, 0.9, by = 0.3),
     labels = seq(0, 0.9, by = 0.3)*sulfide_max,
     tck = 0.05,
     mgp = c(3, -0.25, 0),
     las = 2,
     cex.axis = 0.6,
     hadj = 1)
text(x = 0.5,
     y = 0.5,
     labels = "Sulfide (µM)",
     cex = 0.7,
     srt = 90)
i <- i + 1
# Sulfide barplots
for (date_depth_label in unique(plot_data$date_depth)) {
  screen(i)
  par(mar = margins_to_use)
  plot_sulfide(DATE_DEPTH = date_depth_label)
  i <- i + 1
}


#### Plot: metagenomes totals ####
# Total metagenomes y-axis
MG_totals_max <- 15
screen(i)
par(mar = margins_to_use)
empty_plot_function()
axis(side = 4,
     at = seq(0, 0.9, by = 0.3),
     labels = seq(0, 0.9, by = 0.3)*MG_totals_max,
     tck = 0.05,
     mgp = c(3, -0.25, 0),
     las = 2,
     cex.axis = 0.6,
     hadj = 1)
text(x = 0.6,
     y = 0.5,
     labels = expression('Abundance\n',
                         'of '*italic('hgcA')*' (%)'),
     cex = 0.7,
     srt = 90)
i <- i + 1
# Metagenomes totals barplots
for (date_depth_label in unique(plot_data$date_depth)) {
  screen(i)
  par(mar = margins_to_use)
  plot_omics_totals(OMIC_TYPE = "MG",
                    DATE_DEPTH = date_depth_label)
  i <- i + 1
}


#### Plot: metagenomes guilds ####
# Metagenomes guilds y-axis
screen(i)
par(mar = margins_to_use)
empty_plot_function()
axis(side = 4,
     at = seq(0, 1, by = 0.25),
     labels = seq(0, 1, by = 0.25)*100,
     tck = 0.05,
     mgp = c(3, -0.25, 0),
     las = 2,
     cex.axis = 0.6,
     hadj = 1)
text(x = 0.6,
     y = 0.5,
     labels = expression('Abundance of metabolic\n',
                         'guilds with '*italic('hgcA')*' (%)'),
     cex = 0.7,
     srt = 90)
i <- i + 1
# Metagenomes guilds barplots
for (date_depth_label in unique(plot_data$date_depth)) {
  screen(i)
  par(mar = margins_to_use)
  plot_omic_guilds(OMIC_TYPE = "MG",
                   DATE_DEPTH = date_depth_label)
  i <- i + 1
}


#### Plot: Metatranscriptomes totals ####
# Total metatranscriptomes y-axis
MT_totals_max <- 8
screen(i)
par(mar = margins_to_use)
empty_plot_function()
axis(side = 4,
     at = seq(0, 1, by = 0.25),
     labels = seq(0, 1, by = 0.25)*MT_totals_max,
     tck = 0.05,
     mgp = c(3, -0.25, 0),
     las = 2,
     cex.axis = 0.6,
     hadj = 1)
text(x = c(0.27, 0.52, 0.77),
     y = 0.5,
     labels = c('Expression',
                expression('of '*italic('hgcA')),
                expression('(10'^6*' copies/L)')),
     cex = 0.65,
     srt = 90)
i <- i + 1
# Total metatranscriptomes barplots
for (date_depth_label in unique(plot_data$date_depth)) {
  screen(i)
  par(mar = margins_to_use)
  plot_omics_totals(OMIC_TYPE = "MT",
                    DATE_DEPTH = date_depth_label)
  i <- i + 1
}


#### Plot: Metatranscriptomes guilds ####
# Metatranscriptomes guilds y-axis
screen(i)
par(mar = margins_to_use)
empty_plot_function()
axis(side = 4,
     at = seq(0, 1, by = 0.25),
     labels = seq(0, 1, by = 0.25)*100,
     tck = 0.05,
     mgp = c(3, -0.25, 0),
     las = 2,
     cex.axis = 0.6,
     hadj = 1)
text(x = c(0.35, 0.6),
     y = 0.5,
     labels = c('Relative transcription',
                expression('of '*italic('hgcA')*' by guild (%)')),
     cex = 0.7,
     srt = 90)
i <- i + 1
# Metatranscriptomes guilds barplots
for (date_depth_label in unique(plot_data$date_depth)) {
  screen(i)
  par(mar = margins_to_use)
  plot_omic_guilds(OMIC_TYPE = "MT",
                          DATE_DEPTH = date_depth_label)
  i <- i + 1
}

#### Add legend ####
screen(i)
par(mar = margins_to_use)
empty_plot_function()
legend("top",
       legend = names(color_vector),
       fill = color_vector,
       border = "black")

dev.off()
