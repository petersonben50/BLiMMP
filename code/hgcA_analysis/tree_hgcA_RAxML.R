#### code/hgcA_analysis/tree_hgcA_RAxML.R ####
# Benjamin D. Peterson

# This script will look at the ML tree
# we generated using FastTree with the
# Hg-MATE database.
# It will save out a PDF of that for
# inspection.



#### Always start with a clean slate ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP/")
library(ggtree)
library(phangorn)
library(readxl)
library(tidyverse)
library(treeio)
cb.translator <- readRDS("~/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")



#### Read in needed files ####
hgcA.list <- readLines("dataEdited/hgcA_analysis/dereplication/hgcA_final_list.txt")


#### Read in tree ####
tree.name <- "dataEdited/hgcA_analysis/phylogeny/raxml/RAxML_bipartitions.hgcA"
hgcA.tree.unrooted <- read.newick(tree.name)
rm(tree.name)


#### Root tree ####
# hgcA.tree <- root(phy = hgcA.tree.unrooted,
#                   outgroup = c("paralog_Thermosulfurimonas_dismutans",
#                                "paralog_Candidatus_Omnitrophica"),
#                   edgelabel = TRUE)
hgcA.tree <- midpoint(hgcA.tree.unrooted,
                      edgelabel = TRUE)

#### Read in hgcA to bin info ####
# hgcA.to.bin <- read.table("dataEdited/binning/hgcA/hgcA_to_bin.tsv",
#                           sep = '\t',
#                           header = TRUE)
# hgcA.to.bin.vector <- hgcA.to.bin$binID
# names(hgcA.to.bin.vector) <- hgcA.to.bin$hgcA_ID
# binned.hgcA.index <- which(names(hgcA.renaming.vector) %in% names(hgcA.to.bin.vector))
# hgcA.renaming.vector[binned.hgcA.index] <- paste(hgcA.renaming.vector[binned.hgcA.index], " - ",
#                                                  hgcA.to.bin.vector[names(hgcA.renaming.vector)[binned.hgcA.index]])


#### Read in Hg-MATE seq info ####
hgcA.tree$tip.label[grep("MATE_", hgcA.tree$tip.label)] <- 
  paste(hgcA.tree$tip.label[grep("MATE_", hgcA.tree$tip.label)] %>% strsplit("_") %>% sapply("[", 1),
        hgcA.tree$tip.label[grep("MATE_", hgcA.tree$tip.label)] %>% strsplit("_") %>% sapply("[", 2),
        hgcA.tree$tip.label[grep("MATE_", hgcA.tree$tip.label)] %>% strsplit("_") %>% sapply("[", 3),
        sep = "_")


#### Read in Hg-MATE seq info ####
hg.mate.metadata <- read_xlsx("/Users/benjaminpeterson/Documents/research/Hg_MATE/versions/v1.01142021/Hg-MATE-Db.v1.01142021_catalogue.xlsx") %>%
  rename(group = `microbial group`,
         name = `[Organism Name]_Phylum-Class`) %>%
  mutate(name = gsub("sp._", "sp.", name)) %>%
  mutate(name = paste(name %>% strsplit("_") %>% sapply("[", 1),
                      name %>% strsplit("_") %>% sapply("[", 2),
                      sep = "_"),
         treeID = paste(group,
                        " (", name, "-", MATE_id, ")",
                        sep = ""))
MATE.renaming.vector <- hg.mate.metadata$treeID
names(MATE.renaming.vector) <- hg.mate.metadata$MATE_id


#### Read in info on seqs from Jones et al, 2019 ####
jones.metadata <- read.table("references/jones_bin_names.tsv",
                             sep = '\t',
                             col.names = c("binID", "binName")) %>%
  full_join(read.table("references/jones_genes2bin.tsv",
                       sep = '\t',
                       col.names = c("geneID", "binID"))) %>%
  mutate(binName = paste(gsub("Unclassified", "", binName),
                         " (Jones et al, 2019)",
                         sep = ""))
jones.renaming.vector <- jones.metadata$binName
names(jones.renaming.vector) <- jones.metadata$geneID


#### Get indices ####
reference.indices.mate <- which(hgcA.tree$tip.label %in% names(MATE.renaming.vector))
reference.indices.jones <- which(hgcA.tree$tip.label %in% names(jones.renaming.vector))
this.study.indices <- which(hgcA.tree$tip.label %in% hgcA.list)


#### Change names ####
hgcA.tree$tip.label[reference.indices.mate] <- MATE.renaming.vector[hgcA.tree$tip.label[reference.indices.mate]]
hgcA.tree$tip.label[reference.indices.jones] <- jones.renaming.vector[hgcA.tree$tip.label[reference.indices.jones]]



#### Set color vector ####
color.vector <- rep(cb.translator["bluishgreen"], length(hgcA.tree$tip.label))
color.vector[reference.indices.mate] <- cb.translator["black"]
color.vector[this.study.indices] <- cb.translator["vermillion"]
color.vector[reference.indices.jones] <- cb.translator["skyblue"]

# Visualize tree
pdf("results/hgcA_analysis/hgcA_tree_RAxML_rooted.pdf",
    height = 16,
    width = 5)
ggtree(hgcA.tree,
       aes(x = 0,
           xend = 6)) + 
  geom_tiplab(size=2.5,
              align = TRUE,
              col = color.vector) + 
  geom_nodelab(aes(x = branch),
               vjust = -.3,
               size = 2)
dev.off()




#### Save out tree ####
saveRDS(hgcA.tree,
        "dataEdited/hgcA_analysis/phylogeny/hgcA_clean_tree.rds")


#### Save out color vector ####
saveRDS(color.vector,
        "dataEdited/hgcA_analysis/phylogeny/hgcA_clean_tree_color_vector.rds")
