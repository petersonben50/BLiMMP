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

#### Clean up spreadsheet and save out for reporting in manuscript ####
count.data <- count.data %>%
  select(mtID, MT_reads, IS_reads, percent_non_rRNA_as_IS, rRNA_reads, percent_rRNA, total_reads)
write.csv(count.data,
          "results/seq_data_tables/metatranscriptome_numbers.csv",
          row.names = FALSE)

#### Calculate normalization metric using internal standard ####
length.of.IS <- 1371 # nucleotides
molar.mass.IS <- 441530 # g/mol
mass.IS.added <- 12 # ng
number.of.IS.copies.added <- (mass.IS.added * 10^-9) / molar.mass.IS * (6.02*10^23)

normalization.data <- count.data %>%
  select(mtID, IS_reads) %>%
  mutate(IS_copies = number.of.IS.copies.added,
         IS_length = length.of.IS) %>%
  mutate(NF_copies_per_readsPerBp = number.of.IS.copies.added / (IS_reads / length.of.IS)) %>%
  select(mtID, NF_copies_per_readsPerBp) %>%
  left_join(read.csv("metadata/metatranscriptome_metadata.csv") %>%
              select(metatranscriptomeID, volumeFiltered) %>%
              rename(mtID = metatranscriptomeID)) %>%
  mutate(NF_counts_per_readsPerBp_per_L = NF_copies_per_readsPerBp / (volumeFiltered / 1000)) %>%
  as.data.frame()



#### Convert to vector and save out normalization data ####
normalization.vector <- normalization.data$NF_counts_per_readsPerBp_per_L
names(normalization.vector) <- normalization.data$mtID
saveRDS(normalization.vector,
        "dataEdited/metatranscriptomes/normalization_vector.rds")



#### Set up pseudomapping key ####
# Map metatranscriptomes to ORFs from assemblies from the same year.
assemblies <- read.table("dataEdited/assemblies/all_assemblies_stats.txt",
                         header = TRUE) %>%
  select(assemblyID) %>%
  unlist(use.names = FALSE)
metatranscriptomes <- count.data %>%
  select(mtID) %>%
  unlist(use.names = FALSE)

# Join by year
pseudomapping.key <- cbind(rep(metatranscriptomes, length(assemblies)),
                           rep(assemblies, length(metatranscriptomes)) %>%
                             sort())
write.table(pseudomapping.key,
            file = "dataEdited/metatranscriptomes/reports/pseudomapping_key.tsv",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
