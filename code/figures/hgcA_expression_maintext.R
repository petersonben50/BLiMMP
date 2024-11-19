#### code/figures/hgcA_expression_maintext.R ####
# Written for BLiMMP project
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP")
library(gggenes)
library(ggpubr)
library(tidyverse)
cb.translator <- readRDS("references/colorblind_friendly_colors.rds")

#### Read in data ####
bin_metadata <- readRDS(file = 'dataEdited/bin_based_analyses/bin_data_aggregate_3.rds')


#### Read in omic and sulfide metadata ####
MT_metadata <- read.csv("metadata/metatranscriptome_metadata.csv") %>%
  dplyr::rename(omicID = metatranscriptomeID) %>%
  mutate(seqType = "MT",
         date_depth = paste(startDate, ":", depth, "m",
                            sep = "")) %>%
  select(omicID, date_depth) %>%
  left_join(read.csv("dataFinal/water_chem_data.csv") %>%
              group_by(date, depth) %>%
              summarize(sulfide_ppm = mean(sulfide_ppm, na.rm = TRUE)) %>%
              mutate(date_depth = paste(date, ":", depth, "m",
                                        sep = "")) %>%
              ungroup()) %>%
  select(omicID, date_depth, sulfide_ppm)
sulfide_order <- MT_metadata %>%
  arrange(sulfide_ppm) %>%
  select(date_depth) %>% unlist(use.names = FALSE) %>%
  gsub(":", "\n", .)


#### Function to plot RNA expression ####
plot_bin_MT_genes <- function(constituent_1,
                              constituent_2,
                              label_1,
                              label_2,
                              ylim_to_use = NULL,
                              legend_position = "none") {
  # Add KIR notation to meta_code
  plot_data <- bin_metadata
  plot_data$meta_code[bin_metadata$class == "Kiritimatiellae"] <- "KIR"
  plot_data$meta_code <- gsub(pattern = "OB_",
                              replacement = "",
                              plot_data$meta_code)

  plot_data_1 <- plot_data %>%
    select(meta_code, hms_id, all_of(grep("_MT_",
                                          names(plot_data),
                                          value = TRUE))) %>%
    gather(key = gene_mt_id,
           value = gene_transcription,
           all_of(grep("_MT_",
                       names(plot_data),
                       value = TRUE))) %>%
    mutate(gene = gene_mt_id %>% strsplit("-") %>% sapply("[", 1),
           mt_id = gene_mt_id %>% strsplit("-") %>% sapply("[", 2)) %>%
    select(-gene_mt_id) %>%
    filter(!is.na(gene_transcription)) %>%
    group_by(meta_code, gene, hms_id, mt_id) %>%
    summarise(gene_transcription = mean(gene_transcription)) %>%
    spread(key = gene,
           value = gene_transcription)

  if ("gyrB" %in% c(constituent_1, constituent_2)) {
    plot_data_1 <- plot_data_1 %>%
      filter(gyrB > 0)
  }
  if ("hgcA" %in% c(constituent_1, constituent_2)) {
    plot_data_1$hgcA[plot_data_1$hgcA == 0] <- 0.004
  }
  plot_data_1[, "constituent_1"] <- plot_data_1[, constituent_1]
  plot_data_1[, "constituent_2"] <- plot_data_1[, constituent_2]

  plot_data_1 %>%
    ggplot(aes(x = constituent_1,
               y = constituent_2,
               color = meta_code)) +
    geom_point() +
    scale_color_manual(values = color_vector) +
    xlab(label_1) +
    ylab(label_2) +
    theme_bw() +
    scale_y_continuous(limits = ylim_to_use,
                       transform = "log10") +
    theme(legend.position = legend_position,
          legend.title = element_blank(),
          legend.box.background = element_rect(fill = "white",
                                               colour = "black",
                                               linewidth = 1),
          axis.text = element_text(colour = "black"),
          axis.title = element_text(size = 10)
    )
}



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
  arrange(desc(total_abund)) %>%
  select(hms_id, total_abund) %>%
  as.data.frame()


#### Look at gene neighborhood of 15 most expressed hgcA seqs from bins ####
# Add KIR notation to meta_code
bin_metadata$meta_code[bin_metadata$class == "Kiritimatiellae"] <- "KIR"
bin_metadata$meta_code <- gsub(pattern = "OB_",
                               replacement = "",
                               bin_metadata$meta_code)
# Add gene neighborhood info to bin data
hgcA_GN <- bin_metadata %>%
  mutate(scaffoldID = paste(str_split(hgca_id, "_") %>% sapply("[", 1),
                            str_split(hgca_id, "_") %>% sapply("[", 2),
                            str_split(hgca_id, "_") %>% sapply("[", 3),
                            sep = "_")) %>%
  select(hms_id, bin_id, scaffoldID, meta_code, class) %>%
  left_join(read.csv("dataFinal/hgcA_GNs.csv") %>%
              select(scaffoldID, seqID, gene_start, gene_end, gene_direction, gene_info) %>%
              filter(gene_start > 2000,
                     gene_end < 9000))
hgcA_GN$gene_info[which(is.na(hgcA_GN$gene_info))] <- "NA"
# Pick out bins to use:
#view(hgcA_GN)
bin_to_use <- c("BLI21_coassembly_anvio_bin_0107", "BLI21_coassembly_anvio_bin_0085", "BLI21_coassembly_anvio_bin_0103",
                "BLI21_assembly102_anvio_bin_0033", "BLI21_coassembly_anvio_bin_0147", "BLI21_coassembly_anvio_bin_0003",
                "BLI21_coassembly_anvio_bin_0048", "BLI21_assembly106_anvio_bin_0020", "BLI21_coassembly_anvio_bin_0005",
                "BLI21_coassembly_anvio_bin_0019", "BLI20_coassembly_dasTool_bin_0110", "BLI21_coassembly_anvio_bin_0021",
                "BLI21_coassembly_anvio_bin_0184", "BLI21_coassembly_anvio_bin_0128", "BLI21_coassembly_anvio_bin_0092")
bin_to_use <- bin_to_use[length(bin_to_use):1]
hgcA_GN <- hgcA_GN %>%
  filter(bin_id %in% bin_to_use)

# Set up bin info and set order
# hms_to_bin <- paste(gsub("BLI_hgcA_", "", bin_metadata$hms_id), ":",
#                     bin_metadata$meta_code, ",",
#                     # bin_metadata$class,
#                     sep = "")
# names(hms_to_bin) <- bin_metadata$bin_id
# hms_order <- hms_to_bin[bin_to_use]
# names(hms_order) <- NULL
hgcA_GN <- hgcA_GN %>%
  # mutate(bin_info = as.factor(paste(gsub("BLI_hgcA_", "", hms_id), ":",
  #                                   meta_code, ",",
  #                                   # class,
  #                                   sep = "")),
  #        bin_info = fct_relevel(bin_info, hms_order[length(hms_order):1])) %>%
  mutate(bin_id = fct_relevel(bin_id, bin_to_use))

# Remove the unused gene annotations
hgcA_GN$gene_info[which(hgcA_GN$gene_info == "metal_efflux_RND_TolC")] <- "NA"
hgcA_GN$gene_info[which(hgcA_GN$gene_info == "transc_reg_sig54")] <- "NA"
hgcA_GN$gene_info[which(is.na(hgcA_GN$gene_info))] <- "Other"

# Set up naming and color vectors for genes
naming_vector_genes <- c(expression(italic(hgcA)), expression(italic(hgcB)),
                         expression(italic(acr3)), expression(italic(arsC)), expression(italic(arsR)*'-like'),
                         expression(italic(terC)), "other")
names(naming_vector_genes) <- c("hgcA", "hgcB",
                                "ACR3_arsenite_efflux", "ArsC_arsenate_reductase", "transc_reg_arsR",
                                "terC", "Other")

color_vector_genes <- c(cb.translator[c("black", "skyblue",
                                        "orange", "yellow", "vermillion", "reddishpurple")], "white")
names(color_vector_genes) <- names(naming_vector_genes)

# Set up color vector for metabolic groups
color_vector <- c(cb.translator[c("black", "blue", "yellow", "bluishgreen")], "gray")
names(color_vector) <- c("KIR", "SRB", "RESP", "FERM", "UNK")

# Set up color vectors for bins
bin_meta_key <- bin_metadata %>%
  mutate(color_code = color_vector[meta_code])
bin_to_hms_key <- bin_meta_key$hms_id
names(bin_to_hms_key) <- bin_meta_key$bin_id

color_vector_key <- bin_meta_key$color_code
names(color_vector_key) <- bin_meta_key$hms_id
color_vector_bins <- color_vector_key[bin_to_hms_key[bin_to_use]]
names(color_vector_bins) <- gsub("BLI_hgcA_", "", names(color_vector_bins))
rm(bin_meta_key, color_vector_key)

# Refactor the gene names
hgcA_GN$gene_info <- fct_relevel(hgcA_GN$gene_info,
                                 names(naming_vector_genes))

#### Prepare hgcA:gyrB ratio data ####
transcription_reg <- readRDS("working_directory/transcriptional_regulation_vector.rds")
# Select HMS:
hms_to_use <- bin_metadata %>%
  filter(bin_id %in% bin_to_use) %>%
  select(hms_id) %>%
  unlist(use.names = FALSE)
hgcA_gyrB_ratio_data <- bin_metadata %>%
  filter(hms_id %in% hms_to_use) %>%
  # filter(bin_id %in% bin_to_use) %>%
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
hgcA_gyrB_ratio_data$hgcA[which(hgcA_gyrB_ratio_data$hgcA == 0)] <- min(hgcA_gyrB_ratio_data$hgcA[which(hgcA_gyrB_ratio_data$hgcA > 0)], na.rm = TRUE) / 2
hgcA_gyrB_ratio_data <- hgcA_gyrB_ratio_data %>%
  mutate(hgcA_gyrB_ratio = hgcA / gyrB)





#### Generate hgcA vs. other metrics plots ####
hgcA_by_gyrB <- plot_bin_MT_genes(constituent_1 = "gyrB",
                                  constituent_2 = "hgcA",
                                  label_1 = expression(italic(gyrB)*" (10"^"6"*" transcripts/L)"),
                                  label_2 = expression(italic(hgcA)*" mRNA (10"^"6"*" copies/L)"),
                                  ylim_to_use = c(0.0039, 1.5),
                                  legend_position = c(0.2, 0.8)) +
  geom_abline() +
  scale_x_continuous(limits = c(log(0.0039, 10), 10),
                     transform = "log10",
                     breaks = c(0.01, 0.1, 1, 10),
                     labels = c("0.01", "0.10", "1.0", "10")) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank())

hgcA_by_total_RNA <- plot_bin_MT_genes(constituent_1 = "total",
                                       constituent_2 = "hgcA",
                                       label_1 = expression("mRNA (10"^"6"*" transcripts/L)"),
                                       label_2 = expression(italic(hgcA)*" mRNA (10"^"6"*" copies/L)"),
                                       ylim_to_use = c(0.0039, 1.5)) +
  scale_x_continuous(limits = c(8, 130000),
                     transform = "log10",
                     breaks = c(10, 100, 1000, 10000, 100000),
                     labels = c("10", "", "1,000", "", "100,000")) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank())

scatterplots <- ggarrange(hgcA_by_gyrB, hgcA_by_total_RNA,
                                  ncol = 2,
                                  nrow = 1)

#### Generate hgcA:gyrB ratio plot ####
color_vector_transc_reg <- cb.translator[c("vermillion",
                                           "gray25")]
names(color_vector_transc_reg) <- c("arsR", "no_reg")
transcriptional_regulators <- hgcA_gyrB_ratio_data %>%
  left_join(MT_metadata) %>%
  filter(hgcA > 0) %>%
  group_by(sulfide_ppm, date_depth, hms_id, trans_reg, class) %>%
  summarise(hgcA_gyrB_ratio = mean(hgcA_gyrB_ratio)) %>%
  mutate(date_depth = gsub(":", "\n", date_depth),
         date_depth = fct_relevel(date_depth, sulfide_order)) %>%
  ggplot(aes(x = class,
             y = log(hgcA_gyrB_ratio, 10),
             col = trans_reg)) +
  geom_boxplot() +
  scale_color_manual(values = color_vector_transc_reg,
                     labels = c(expression(italic(arsR)*"-like"),
                                "No regulator")) +
  theme_bw() +
  ylim(c(-1.75, 0.75)) +
  ylab(expression(italic(hgcA)*':'*italic(gyrB)*' ratio (log)')) +
  xlab("") +
  theme(legend.position = c(0.15, 0.85),
        legend.title = element_blank(),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 90),
        legend.box.background = element_rect(fill = "white",
                                             colour = "black",
                                             linewidth = 1))

graphical_plots <- ggarrange(scatterplots,
                             transcriptional_regulators,
                             nrow = 2, ncol = 1)



#### Gene neighborhood plot
GN_plot <- ggplot(hgcA_GN,
                  aes(xmin = gene_start,
                      xmax = gene_end,
                      y = gsub(pattern = "BLI_hgcA_", "", hms_id),
                      fill = gene_info)) +
  geom_gene_arrow() +
  scale_fill_manual(values = color_vector_genes,
                    labels = naming_vector_genes) +
  ylab("") +
  theme_bw() +
  theme(legend.position = "top",
        legend.title = element_blank(),
        axis.text.y = element_text(colour = color_vector_bins))

#### Plot up gene neighborhood ####
cairo_pdf("results/figures/hgcA_transcription_maintext.pdf",
          family = "Arial",
          width = 7.2,
          height = 6)
ggarrange(graphical_plots,
          GN_plot,
          nrow = 1,
          ncol = 2,
          widths = c(1,1))
dev.off()
