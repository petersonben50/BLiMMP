#### code/metagenomes/metagenome_overview.R ####
# Benjamin D. Peterson

# This file contains the code to combine the read count,
# read coverage, fastp, and Nonpareil data for each of
# the metagenomes from the BLiMMP project.


#### Clean-up ####
rm(list = ls())
setwd("~/Documents/research/BLiMMP")
library(dplyr)
library(lubridate)
library(tidyr)
library(readxl)



#### Read in coverage and read data ####
coverage.pre <- read.table("dataEdited/metagenomes/reports/metagenome_coverage_pre_trimming.tsv",
                           header = TRUE) %>%
  rename(R1_pre = R1,
         R2_pre = R2) %>%
  select(-sequencingID)
coverage.post <- read.table("dataEdited/metagenomes/reports/metagenome_coverage.tsv",
                            header = TRUE)
reads.pre <- read.table("dataEdited/metagenomes/reports/metagenome_read_count_pre_trimming.tsv",
                        header = TRUE) %>%
  rename(forwardReads_pre = forwardReads,
         reverseReads_pre = reverseReads)
reads.post <- read.table("dataEdited/metagenomes/reports/metagenome_read_count.tsv",
                         header = TRUE)


#### Check out fraction of read coverage lost due to trimming ####
full_join(coverage.pre, coverage.post) %>%
  mutate(fraction_coverage_retained = (R1 + R2 + single + merged) / (R1_pre + R2_pre)) %>%
  select(metagenomeID, fraction_coverage_retained)


#### Percentage of reads that were merged ####
full_join(reads.pre %>%
            select(metagenomeID, forwardReads_pre),
          reads.post %>%
            select(metagenomeID, mergedReads)) %>%
  mutate(fraction_merged = mergedReads / forwardReads_pre)


#### Save out pre-trimming read coverage ####
write.csv(reads.pre %>%
            rename(paired_reads_preTrim = forwardReads_pre) %>%
            select(metagenomeID, paired_reads_preTrim),
          "results/seq_data_tables/metagenome_read_counts.csv",
          row.names = FALSE)
