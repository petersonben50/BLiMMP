#### analyses/hgcA_gyrB_ratio_model_selection.R ####
# Written for BLiMMP project
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP")
library(tidyverse)
library(AICcmodavg)


#### Read in data ####
bin_metadata <- readRDS(file = 'dataEdited/bin_based_analyses/bin_data_aggregate_3.rds')
bin_metadata$meta_code[bin_metadata$class == "Kiritimatiellae"] <- "KIR"
bin_metadata$meta_code <- gsub(pattern = "OB_",
                               replacement = "",
                               bin_metadata$meta_code)


#### Read in omic and sulfide metadata ####
MT_metadata <- read.csv("metadata/metatranscriptome_metadata.csv") %>%
  dplyr::rename(omicID = metatranscriptomeID) %>%
  mutate(seqType = "MT",
         date_depth = paste(startDate, ":", depth, "m",
                            sep = "")) %>%
  select(omicID, date_depth) %>%
  left_join(read.csv("dataFinal/water_chem_data.csv") %>%
              group_by(date, depth) %>%
              summarize(sulfide_uM = mean(sulfide_uM, na.rm = TRUE)) %>%
              mutate(date_depth = paste(date, ":", depth, "m",
                                        sep = "")) %>%
              ungroup()) %>%
  select(omicID, date_depth, sulfide_uM)
# sulfide_order <- MT_metadata %>%
#   arrange(sulfide_uM) %>%
#   select(date_depth) %>% unlist(use.names = FALSE)


#### Identify most expressed hgcA seqs from bins ####
total_hgcA_transc <- bin_metadata %>%
  select(hms_id, class, meta_code, all_of(grep("hgcA-B", 
                                               names(bin_metadata),
                                               value = TRUE))) %>%
  gather(key = mtID,
         value = abundance,
         -c(1:3)) %>%
  filter(!is.na(abundance)) %>%
  group_by(hms_id, class, meta_code, mtID) %>%
  summarise(abundance = mean(abundance)) %>%
  ungroup() %>%
  group_by(hms_id, class, meta_code) %>%
  summarise(total_abund = sum(abundance)) %>%
  ungroup() %>%
  arrange(desc(total_abund)) %>%
  select(hms_id, total_abund) %>%
  as.data.frame()
# Pick out bins to use:
bins_to_use <- c("BLI21_coassembly_anvio_bin_0107", "BLI21_coassembly_anvio_bin_0085", "BLI21_coassembly_anvio_bin_0103",
                 "BLI21_assembly102_anvio_bin_0033", "BLI21_coassembly_anvio_bin_0147", "BLI21_coassembly_anvio_bin_0003",
                 "BLI21_coassembly_anvio_bin_0048", "BLI21_assembly106_anvio_bin_0020", "BLI21_coassembly_anvio_bin_0005",
                 "BLI21_coassembly_anvio_bin_0019", "BLI20_coassembly_dasTool_bin_0110", "BLI21_coassembly_anvio_bin_0021",
                 "BLI21_coassembly_anvio_bin_0184", "BLI21_coassembly_anvio_bin_0128", "BLI21_coassembly_anvio_bin_0092")
bins_to_use <- bins_to_use[length(bins_to_use):1]
# Select HMSs:
hms_to_use <- bin_metadata %>%
  filter(bin_id %in% bins_to_use) %>%
  select(hms_id) %>%
  unlist(use.names = FALSE)



#### Prepare hgcA:gyrB ratio data ####
transcription_reg <- readRDS("working_directory/transcriptional_regulation_vector.rds")
hgcA_gyrB_ratio_data <- bin_metadata %>%
  filter(hms_id %in% hms_to_use) %>%
  mutate(trans_reg = transcription_reg[cluster_ID]) %>%
  select(hms_id, bin_id, trans_reg, meta_code, class, all_of(grep("-BLI",
                                                                  names(bin_metadata)))) %>%
  gather(key = code,
         value = value,
         -c(hms_id, bin_id, trans_reg, meta_code, class)) %>%
  group_by(hms_id, code, trans_reg, meta_code, class) %>%
  summarise(value = mean(value, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(gene = strsplit(code, "-") %>% sapply("[", 1),
         omicID = strsplit(code, "-") %>% sapply("[", 2)) %>%
  filter(gene %in% c("hgcA", "gyrB")) %>%
  select(-code) %>%
  spread(key = gene,
         value = value) %>%
  filter(!is.na(gyrB)) 
# hgcA_gyrB_ratio_data$hgcA[which(hgcA_gyrB_ratio_data$hgcA == 0)] <- min(hgcA_gyrB_ratio_data$hgcA[which(hgcA_gyrB_ratio_data$hgcA > 0)], na.rm = TRUE) / 2
hgcA_gyrB_ratio_data <- hgcA_gyrB_ratio_data %>%
  mutate(hgcA_gyrB_ratio = hgcA / gyrB) %>%
  left_join(MT_metadata) %>%
  filter(gyrB > 0) %>%
  filter(hgcA > 0)

rm(total_hgcA_transc, MT_metadata, bin_metadata,
   bins_to_use, hms_to_use, transcription_reg)


#### Statistical testing ####
Cand.mod <- list()
Cand.mod[["reg_only"]] <- aov(log(hgcA_gyrB_ratio, 10) ~ trans_reg,
                           data = hgcA_gyrB_ratio_data)
Cand.mod[["sulfide_only"]] <- aov(log(hgcA_gyrB_ratio, 10) ~ sulfide_uM,
                                 data = hgcA_gyrB_ratio_data)
Cand.mod[["interaction"]] <- aov(log(hgcA_gyrB_ratio, 10) ~ sulfide_uM * trans_reg,
                                 data = hgcA_gyrB_ratio_data)
Cand.mod[["no_interaction"]] <- aov(log(hgcA_gyrB_ratio, 10) ~ sulfide_uM + trans_reg,
                                   data = hgcA_gyrB_ratio_data)
aictab(cand.set = Cand.mod,modnames = names(Cand.mod))

summary(Cand.mod[["interaction"]])
summary(Cand.mod[["no_interaction"]])
summary(Cand.mod[["reg_only"]])
summary(Cand.mod[["sulfide_only"]])
