#### code/hgcA_analysis/aggregate_hgcA_info.R ####
# Benjamin D. Peterson


#### Start with a clean slate ####
rm(list = ls())
setwd("~/Documents/research/BLiMMP/")
library(readxl)
library(tidyverse)
mt.NF.vector <- readRDS("dataEdited/metatranscriptomes/normalization_vector.rds")


#### Generate data frame of seqs ####
hgcA.list <- readLines("dataEdited/hgcA_analysis/identification/hgcA_good.txt")
hgcA.df <- data.frame(hgcA_ID = hgcA.list,
                      scaffoldID = paste(hgcA.list %>% strsplit("_") %>% sapply("[",1),
                                         hgcA.list %>% strsplit("_") %>% sapply("[",2),
                                         hgcA.list %>% strsplit("_") %>% sapply("[",3),
                                         sep = "_"))
rm(hgcA.list)



#### Include hgcB info ####
####################### NEED TO LOOK MANUALLY FOR HGCB ##################
# hgcB.list <- readLines("dataEdited/hgcA_analysis/hgcB/downstream_genes_present.txt")
# hgcB.df <- data.frame(hgcB_ID = hgcB.list,
#                       scaffoldID = paste(hgcB.list %>% strsplit("_") %>% sapply("[",1),
#                                          hgcB.list %>% strsplit("_") %>% sapply("[",2),
#                                          hgcB.list %>% strsplit("_") %>% sapply("[",3),
#                                          sep = "_"))
# hgcA.df <- left_join(hgcA.df,
#                      hgcB.df)
# rm(hgcB.df,
#    hgcB.list)


#### Add auto classification data ####
hgcA.classification <- read.csv("dataEdited/hgcA_analysis/classification/hgcA_taxonomy_table.csv") %>%
  mutate(classification = paste(phylum, class, order,
                                family, genus, species,
                                sep = ";"),
         hgcA_ID = seqID) %>%
  select(hgcA_ID, classification)
hgcA.df <- hgcA.df %>%
  full_join(hgcA.classification)
rm(hgcA.classification)


#### Add cluster IDs ####
clustering.info <- read.table("dataEdited/hgcA_analysis/dereplication/hgcA_good_acrossYear.tsv",
                              header = TRUE) %>%
  mutate(hgcA_ID = id) %>%
  mutate(representative = (clstr_rep == 1)) %>%
  select(hgcA_ID, clstr, representative)
hgcA.df <- hgcA.df %>%
  full_join(clustering.info)
rm(clustering.info)




#### Add depth information ####

depth.df <- read.csv("dataEdited/hgcA_analysis/depth/hgcA_coverage.csv") %>%
  select(hgcA_ID, metagenomeID, coverage) %>%
  spread(key = metagenomeID,
         value = coverage,
         fill = 0)
hgcA.df <- hgcA.df %>%
  full_join(depth.df)



#### Add hgcA MT abundance info ####
MT.data <- read.table("dataEdited/hgcA_analysis/hgcA_MT_hits_clean.tsv",
                      sep = '\t',
                      col.names = c("mtID", "hgcA_ID", "length_gene", "eff_length", "counts", "tpm")) %>%
  mutate(reads_per_base = counts / length_gene,
         copies_per_liter = round(reads_per_base * mt.NF.vector[mtID])) %>%
  select(hgcA_ID, mtID, copies_per_liter) %>%
  spread(key = mtID,
         value = copies_per_liter,
         fill = 0)
hgcA.df <- hgcA.df %>%
  full_join(MT.data)
write.csv(hgcA.df,
          "~/Downloads/hgcA_data.csv",
          row.names = FALSE)
