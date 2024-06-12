#### code/figures/hgcA_expression_maintext.R ####
# Written for BLiMMP project.
# This script is intended to visualize the phylogenetic
# distribution of the arsR-like transcriptional regulator.
# Benjamin D. Peterson

#### Clean up ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP")
library(gggenes)
library(ggnewscale)
library(ggpubr)
library(ggtree)
library(tidyverse)
library(treeio)
cb.translator <- readRDS("references/colorblind_friendly_colors.rds")



#### Read in data ####
hgcA_tree_unrooted <- read.newick("dataEdited/ABA/hgcA/RAxML_bipartitions.hgcA_for_phylogeny")
hgcA_data <- read.csv("dataFinal/hgcA_data.csv")
hgcA_reg_vector <- readRDS("working_directory/transcriptional_regulation_vector.rds")


#### Prepare key to link sequences to regulatory elements ####
hgcA_data <- hgcA_data %>%
  mutate(transc_reg = hgcA_reg_vector[cluster_ID])
hgcA_data$transc_reg[which(is.na(hgcA_data$transc_reg))] <- "unknown"
# hgcA_reg_key <- hgcA_data$transc_reg
# names(hgcA_reg_key) <- hgcA_data$seqID
# hgcA_reg_key[which(is.na(hgcA_reg_key))] <- "unknown"
# 

#### Set up color and naming vector for tree to reflect regulatory elements ####
color_vector_reg <- c(cb.translator["vermillion"], "gray25", "gray85", cb.translator["reddishpurple"])
names(color_vector_reg) <- c("arsR", "no_reg", "unknown", "marR")
naming_vector_reg <- c(expression(italic(arsR)*'-like regulator'),
                       "No known regulator identified",
                       "Unknown",
                       expression(italic(marR)*'-like regulator'))
names(naming_vector_reg) <- names(color_vector_reg)
hgcA_data$transc_reg <- fct_relevel(hgcA_data$transc_reg, names(color_vector_reg))


#### Set up color vector for tree to reflect metabolic assignments ####
color_vector_metab <- c(cb.translator[c("black", "blue", "yellow", "bluishgreen")], "gray")
names(color_vector_metab) <- c("KIR", "SRB", "RESP", "FERM", "UNK")

#### Set up color vector for tree to reflect taxonomic assignments ####
hgcA_data %>%
  filter(clstr_rep == 1) %>%
  group_by(taxonomic_assignment) %>%
  summarise(counts = n())
# Combine Geobacter, Chloroflexi, and UNK
hgcA_data$taxonomic_assignment[which(hgcA_data$taxonomic_assignment %in% c("Geobacteraceae",
                                                                           "Chloroflexi",
                                                                           "UNK"))] <- "Unknown/other"
# Combine non-Kiritimatiellaeota into Other PVC
hgcA_data$taxonomic_assignment[which(hgcA_data$taxonomic_assignment %in% c("Verrucomicrobiae",
                                                                           "Planctomycetes",
                                                                           "PVC"))] <- "Other PVC"
sort(unique(hgcA_data$taxonomic_assignment))
# Color source: https://davidmathlogic.com/colorblind/#%23648FFF-%23785EF0-%23DC267F-%23FE6100-%23FFB000-%23DCDCDC-%23000000-%23616060
color_vector_tax <- c('#DC267F', '#FE6100', '#648FFF', '#785EF0',
                      '#FFB000', '#000000', '#616060', "#DCDCDC")
names(color_vector_tax) <- c("Actinobacteria", "Bacteroidetes", "Desulfobacterota", "Desulfomonilaceae",
                             "Firmicutes", "Kiritimatiellae", "Other PVC", "Unknown/other")
  
#### Prepare heatmap data ####
hgcA_data_plotting <- hgcA_data %>%
  select(taxonomic_assignment, metabolic_assignment)
rownames(hgcA_data_plotting) <- hgcA_data$seqID


#### Generate figure ####
# Initial tree plot
hgcA_tree <- ggtree(hgcA_tree_unrooted,
                    layout = "fan") %<+% (hgcA_data %>%
                                            mutate(cluster_ID = gsub("BLI_hgcA_", "", cluster_ID))) + 
  geom_tippoint(aes(color = transc_reg),
                size = 4) +
  scale_color_manual(values = color_vector_reg,
                     labels = naming_vector_reg[names(color_vector_reg)],
                     name = "Regulators present\n(tip labels)") + 
  geom_tiplab2(aes(label = cluster_ID), align=T, linetype=NA, 
               size = 3.5, offset = 0.05, hjust = 0)

# Add the metabolic assignments
hgcA_tree_plus <- gheatmap(hgcA_tree,
                           hgcA_data_plotting %>% select(metabolic_assignment),
                           offset = 0.6, width = 0.1,
                           font.size = 3, hjust = 0.65,
                           custom_column_labels = "Met.") +
  scale_fill_manual(values = c(color_vector_metab),
                    name = "Metabolic assignment\n(inner circle)")


# Add taxonomic assignments
hgcA_tree_plus2 <- hgcA_tree_plus + new_scale_fill()
final_figure <- gheatmap(hgcA_tree_plus2,
                         hgcA_data_plotting %>% select(taxonomic_assignment),
                         offset = 0.9, width = 0.1,
                         font.size = 3, hjust = 0.35,
                         custom_column_labels = "Tax.") +
  scale_fill_manual(values = color_vector_tax,
                    name = "Taxonomic assignment\n(outer circle)") +
  theme(legend.frame = element_rect())

#### Save out figure ####
cairo_pdf("results/figures/hgcA_phylogeny.pdf",
          family = "Ariel",
          height = 10,
          width = 20)
final_figure
dev.off()

