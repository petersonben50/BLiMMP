#### Clean up and prep ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP")
library(tidyverse)
cb.translator <- readRDS("references/colorblind_friendly_colors.rds")



#### Read in and normalize MT data ####
NF_vector <- readRDS("dataEdited/metatranscriptomes/normalization_vector.rds")
bin_MT_data <- read.table("dataEdited/bin_based_analyses/final_bin_data/bin_MT_data.tsv",
                          sep = '\t', header = TRUE) %>%
  gather(key = mt_id,
         value = rna_abund_per_eff_length,
         -c(1:2)) %>%
  mutate(million_counts_per_liter = rna_abund_per_eff_length*NF_vector[mt_id]/1000000) %>%
  select(-rna_abund_per_eff_length) %>%
  ungroup() %>%
  spread(key = mt_id,
         value = million_counts_per_liter)


#### Read in bin data ####
hgcA_bins <- read.csv('dataFinal/hgcA_bins.csv')
# Add KIR notation to meta_code
hgcA_bins$meta_code[grep("Kiritimatiellae", hgcA_bins$gtdb_tax)] <- "KIR"
hgcA_bins$meta_code <- gsub(pattern = "OB_",
                            replacement = "",
                            hgcA_bins$meta_code)

hgcA_bins <- hgcA_bins %>%
  mutate(scaffoldId = paste(strsplit(hgca_id, "_")[[1]][1],
                             strsplit(hgca_id, "_")[[1]][2],
                             strsplit(hgca_id, "_")[[1]][3],
                             sep = "_")) %>%
  select(bin_id, scaffoldId, gtdb_tax, meta_code, all_of(grep("hgcA", colnames(hgcA_bins), value = TRUE)))



#### Find arsR-like TR genes ####
arsR_list <- read.csv("dataEdited/bin_based_analyses/hgcA_geneNeighborhood_annotated.csv") %>%
  filter(gene_info == "transc_reg_arsR") %>%
  select(seqID) %>%
  unlist(use.names = FALSE)
arsR_transcription <- bin_MT_data %>%
  filter(gene_id %in% arsR_list) %>%
  select(-gene_id)
names(arsR_transcription)[-1] <- paste("arsR.", names(arsR_transcription)[-1], sep = "")


#### Merge arsR-like TR genes with hgcA bins ####
arsR_hgcA_data <- hgcA_bins %>%
  inner_join(arsR_transcription)



# Set up color vector for metabolic groups
color_vector <- c(cb.translator[c("black", "blue", "yellow", "bluishgreen")], "gray")
names(color_vector) <- c("KIR", "SRB", "RESP", "FERM", "UNK")



#### Plot hgcA vs. arsR-like TR genes ####
arsR_hgcA_plot <- arsR_hgcA_data %>%
  gather(key = gene_mt_id,
         value = million_counts_per_liter,
         -c(1:4)) %>%
  mutate(gene_name = gene_mt_id %>% strsplit("\\.") %>% sapply("[", 1),
         mt_id = gene_mt_id %>% strsplit("\\.") %>% sapply("[", 2)) %>%
  select(-gene_mt_id) %>% 
  group_by(scaffoldId, gtdb_tax, meta_code, gene_name, mt_id) %>%
  summarise(million_counts_per_liter = mean(million_counts_per_liter)) %>%
  spread(key = gene_name,
         value = million_counts_per_liter) %>%
  ggplot(aes(x = hgcA,
             y = arsR,
             color = meta_code)) +
  geom_point() +
  scale_color_manual(values = color_vector) +
  geom_abline(slope = 1,
              linetype = "dashed",
              color="red") +
  scale_x_continuous(limits = c(0.0008, 5),
                     trans = 'log10') +
  scale_y_continuous(limits = c(0.0008, 5),
                     trans = 'log10') +
  xlab(expression(italic(hgcA)*" mRNA (10"^"6"*" copies/L)")) +
  ylab(expression(italic(arsR)*"-like mRNA (10"^"6"*" copies/L)")) +
  theme_bw() +
  theme(legend.position = c(0.8, 0.3),
        legend.title = element_blank(),
        legend.box.background = element_rect(fill = "white",
                                             colour = "black",
                                             linewidth = 1),
        axis.text = element_text(colour = "black"),
        axis.title = element_text(size = 10))
 


#### Save out plot ####
cairo_pdf("results/figures/arsR_hgcA_plot.pdf",
          width = 4,
          height = 3.5)
arsR_hgcA_plot
dev.off()
