#### code/assemblies/assembly_stats.R ####
# Benjamin D. Peterson

# This file contains the code to generate the mapping
# key file and aggregate the mapping statistics.


#### Always start with a clean slate ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP/")
library(tidyverse)


#### Read in data ####
assembly.data <- read.table("dataEdited/assemblies/all_assemblies_stats.txt",
                            header = TRUE)
ORF.data <- read.table("dataEdited/assemblies/ORF_counts.tsv",
                       header = TRUE)
left_join(ORF.data,
          assembly.data) %>%
  rename(length_of_assembly_bp = n,
         longest_contig_length_bp = max,
         total_length_of_contigs_bp = sum) %>%
  mutate(total_length_of_contigs_Mbp = total_length_of_contigs_bp / 10^6) %>%
  select(assemblyID, length_of_assembly_bp, N50, longest_contig_length_bp,
         total_length_of_contigs_Mbp, ORF_count) %>%
  write.csv("results/seq_data_tables/assembly_data.csv",
            row.names = FALSE)
