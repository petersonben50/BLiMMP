rm(list =ls())
library(tidyverse)
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP")
cb.translator <- readRDS("references/colorblind_friendly_colors.rds")


#### Read in omic data ####
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


#### Read in sulfide data ####
sulfide_data <- read.csv("dataFinal/water_chem_data.csv") %>%
  group_by(date, depth) %>%
  summarize(sulfide_ppm = mean(sulfide_ppm, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(date_depth = paste(date, ":", depth, "m",
                            sep = "")) %>%
  select(date_depth, sulfide_ppm)
omic_metadata <- omic_metadata %>%
  left_join(sulfide_data)


#### Read in the processed TEAP gene data ####
TEAP_data <- read.csv(file = "dataFinal/TEAP_data.csv")


#### Aggregate the data ####
TEAP_data_metadata <- TEAP_data %>%
  gather(key = omicID,
         value = coverage,
         all_of(grep("BLI2", names(TEAP_data),
                     value = TRUE))) %>%
  group_by(clstr_ID, geneName, omicID) %>%
  summarise(coverage = mean(coverage)) %>%
  ungroup() %>%
  left_join(omic_metadata) %>%
  filter(!(startDate == "2020-10-10" & depth == 15.7))
rm(TEAP_data)


#### Function to plot data ####
plot_teap_genes_by_sulfide <- function(genes_to_plot = "red_dsrA",
                                       omic_type = "MT",
                                       ylim_to_use = c(0, 2500),
                                       ylab_to_use = "Reductive dsrA transcription\n(million copies/L)",
                                       color_vector_to_use,
                                       log_scale = FALSE,
                                       legend_position = NULL) {
  
  plotting_data <- TEAP_data_metadata %>%
    filter(geneName %in% genes_to_plot,
           seqType == omic_type) %>%
    group_by(omicID, date_depth, sulfide_ppm, geneName) %>%
    summarise(coverage = sum(coverage)) %>%
    ungroup() %>%
    spread(key = geneName,
           value = coverage) %>%
    arrange(sulfide_ppm) %>%
    as.data.frame()
  
  if (log_scale == TRUE) {
    plot(x = NA,
         y = NA,
         xlim = c(0, 4.6),
         ylim = ylim_to_use,
         xlab = "Sulfide (mg/L)",
         ylab = ylab_to_use,
         log = "y")
  } else {
    plot(x = NA,
         y = NA,
         xlim = c(0, 4.6),
         ylim = ylim_to_use,
         xlab = "Sulfide (mg/L)",
         ylab = ylab_to_use)
    
  }
  
  for (gene_to_plot in genes_to_plot) {
    points(plotting_data$sulfide_ppm,
           plotting_data[, gene_to_plot],
           pch = 21,
           bg = color_vector_to_use[gene_to_plot],
           col = "black",
           cex = 1.5)
  }
  if (!is.null(legend_position)) {
    legend(legend_position,
           legend = naming_vector_TEAPs,
           pt.bg = color_vector_to_use,
           col = "black",
           pch = 21,
           pt.cex = 1.75)
  }

}


#### Generate style vectors ####
color_vector_TEAPs <- cb.translator[c("blue", "reddishpurple")]
names(color_vector_TEAPs) <- c("red_dsrA", "narG")

naming_vector_TEAPs <- c(expression("Reductive "*italic(dsrA)),
                         expression(italic(narG)))
names(naming_vector_TEAPs) <- names(color_vector_TEAPs)

#### Generate plots ####
cairo_pdf("results/figures/TEAP_genes.pdf",
          family = "Arial",
          width = 7.2,
          height = 4)
par(mfrow = c(1, 2),
    mar = c(3, 3, 1, 1),
    tck = -0.008,
    mgp = c(1.5, 0.2, 0))

plot_teap_genes_by_sulfide(genes_to_plot = names(color_vector_TEAPs),
                           omic_type = "MG",
                           ylim_to_use = c(0.1, 10),
                           ylab_to_use = "Gene abundance (%)",
                           color_vector_to_use = color_vector_TEAPs,
                           legend_position = "topright")

plot_teap_genes_by_sulfide(genes_to_plot = names(color_vector_TEAPs),
                           omic_type = "MT",
                           ylim_to_use = c(0.1, 2500),
                           ylab_to_use = "Gene transcription (million copies/L)",
                           color_vector_to_use = color_vector_TEAPs,
                           log_scale = TRUE)
dev.off()
