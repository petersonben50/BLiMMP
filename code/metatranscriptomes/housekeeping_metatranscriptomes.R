#### code/metatranscriptomes/housekeeping_metatranscriptomes.R ####
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP/")
library(tidyverse)



#### Calculate fraction of reads that were rRNA ####
rRNA.data <- read_table("dataEdited/metatranscriptomes/reports/mt_read_counts_rRNA.tsv") %>%
  mutate(percent_rRNA = rRNA_reads / (rRNA_reads + nonrRNA_reads) * 100)
