#### code/metatranscriptomes/housekeeping_metatranscriptomes.R ####
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP/")
library(tidyverse)



#### Calculate fraction of reads that were rRNA ####
rRNA.data <- read_table("dataEdited/metatranscriptomes/reports/mt_read_counts_rRNA.tsv") %>%
  mutate(total_reads = rRNA_reads + nonrRNA_reads,
         percent_rRNA = round((rRNA_reads / (total_reads) * 100),
                              1))


#### Calculate fraction of non-rRNA reads that were IS ####
IS.data <- read_table("dataEdited/metatranscriptomes/reports/mt_read_counts_IS.tsv") %>%
  mutate(total_non_rRNA = IS_reads + MT_reads,
         percent_non_rRNA_as_IS = round(IS_reads / total_non_rRNA * 100, 2))


#### Combine data ####
count.data <- full_join(rRNA.data,
                        IS.data)
# Make sure non-rRNA read counts are right
count.data %>%
  mutate(rRNA_count_discrepancy = (nonrRNA_reads != total_non_rRNA)) %>%
  select(rRNA_count_discrepancy) %>%
  unlist(use.names = FALSE) %>%
  any()

count.data <- count.data %>%
  select(mtID, MT_reads, IS_reads, percent_non_rRNA_as_IS, rRNA_reads, percent_rRNA, total_reads)


#### Calculate normalization metric using internal standard ####
volume.filtered.data <- read.csv("metadata/metatranscriptome_metadata.csv") %>%
  select(metatranscriptomeID, startDate, depth, volumeFiltered)



#### Set up pseudomapping key ####
# Map metatranscriptomes to ORFs from assemblies from the same year.
assemblies <- read.table("dataEdited/assemblies/reports/all_assemblies_stats.txt",
                         header = TRUE) %>%
  select(assemblyID) %>%
  mutate(year = assemblyID %>%
           strsplit("_") %>% sapply("[", 1) %>%
           gsub("BLI", "20", .))
metatranscriptomes <- count.data %>%
  mutate(year = mtID %>%
           strsplit("_") %>% sapply("[", 1) %>%
           gsub("BLI", "20", .)) %>%
  select(mtID, year)

# Join by year
pseudomapping.key <- full_join(metatranscriptomes,
                               assemblies) %>%
  select(-year)
write.table(mapping.key,
            file = "dataEdited/metatranscriptomes/reports/pseudomapping_key.tsv",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
