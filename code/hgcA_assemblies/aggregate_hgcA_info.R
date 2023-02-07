#### code/hgcA_analysis/hgcA_dereplication.R ####
# Benjamin D. Peterson


#### Start with a clean slate ####
rm(list = ls())
setwd("~/Documents/research/BLiMMP/")
library(readxl)
library(tidyverse)


#### Generate data frame of seqs ####
hgcA.list <- readLines("dataEdited/hgcA_analysis/identification/hgcA_good.txt")
hgcA.df <- data.frame(hgcA_ID = hgcA.list,
                      scaffoldID = paste(hgcA.list %>% strsplit("_") %>% sapply("[",1),
                                         hgcA.list %>% strsplit("_") %>% sapply("[",2),
                                         hgcA.list %>% strsplit("_") %>% sapply("[",3),
                                         sep = "_"))
rm(hgcA.list)



#### Include hgcB info ####
hgcB.list <- readLines("dataEdited/hgcA_analysis/hgcB/downstream_genes_present.txt")
hgcB.df <- data.frame(hgcB_ID = hgcB.list,
                      scaffoldID = paste(hgcB.list %>% strsplit("_") %>% sapply("[",1),
                                         hgcB.list %>% strsplit("_") %>% sapply("[",2),
                                         hgcB.list %>% strsplit("_") %>% sapply("[",3),
                                         sep = "_"))
hgcA.df <- left_join(hgcA.df,
                     hgcB.df)
rm(hgcB.df,
   hgcB.list)


#### Add classification data ####
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



# USE AUTO-GENERATED REPRESENTATIVES FOR NOW
hgcA.final.list <- hgcA.df %>%
  filter(representative == TRUE) %>%
  select(hgcA_ID) %>%
  unlist(use.names = FALSE)
hgcA.final.list %>%
  writeLines("dataEdited/hgcA_analysis/dereplication/hgcA_final_list.txt")

hgcA.final.list %>%
  writeLines("dataEdited/hgcA_analysis/dereplication/hgcA_final_abundance_list.txt")



# #### Read out data for dereplication ####
# write.csv(hgcA.df,
#           "dataEdited/hgcA_analysis/dereplication/hgcA_dereplication_data.csv",
#           row.names = FALSE)
# rm(hgcA.df)
# 
# 
# #### Read data back in ####
# hgcA.df <- read_xlsx("dataEdited/hgcA_analysis/hgcA_dereplication.xlsx")
# 
# hgcA.df %>%
#   filter(representative == TRUE) %>%
#   select(seqID) %>%
#   unlist(use.names = FALSE) %>%
#   writeLines("dataEdited/hgcA_analysis/hgcA_rep_list.txt")
# 
# hgcA.df %>%
#   filter(usedForAbundance == TRUE) %>%
#   select(seqID) %>%
#   unlist(use.names = FALSE) %>%
#   writeLines("dataEdited/hgcA_analysis/hgcA_repAbundance_list.txt")
# 
# 
# #### Add in other info about hgcA ####
# hgcA.manual.taxonomy <- read_xlsx("dataEdited/hgcA_analysis/phylogeny/manual_taxonomy.xlsx") %>%
#   left_join(hgcA.df %>% select(seqID, clstr)) %>%
#   select(clstr, manual_classification)
# hgcA.df <- hgcA.df %>%
#   left_join(hgcA.manual.taxonomy)
# 
# 
# #### Save out data ####
# saveRDS(hgcA.df,
#         "dataEdited/hgcA_analysis/hgcA_information.rds")
# write.csv(hgcA.df,
#           "dataEdited/hgcA_analysis/hgcA_information.csv")
