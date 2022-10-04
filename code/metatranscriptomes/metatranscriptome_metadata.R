#### code/metatranscriptomes/metatranscriptome_metadata.R ####
# Benjamin D. Peterson

# This file contains the code to generate the combined metadata
# needed for metatranscriptome processing, storage, and analysis.


#### Always start with a clean slate, ya bum ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP/")
library(lubridate)
library(readxl)
library(tidyverse)


#### Generate table of MT sample site information ####

sample.data <-
# Read in sample metadata
  read_xlsx("metadata/raw_metadata/2_sample_IDs.xlsx") %>%
            select(sampleID, tripID, depth) %>%
  # Read in trip metadata
  left_join(read_xlsx("metadata/raw_metadata/1_trip_IDs.xlsx") %>%
              select(tripID, startDate))

all.data <-
# Read in filter metadata
read_xlsx("metadata/raw_metadata/NA_IDs.xlsx") %>%
  select(filterID, sampleID, volumeFiltered, replicate) %>%
  # Read in extraction metadata
  right_join(read_xlsx("dataEdited/rnaExtractions/RNA_extractions.xlsx") %>%
               select(filterID, extractionID)) %>%
  left_join(sample.data) %>%
  # Read in metagenome metadata
  right_join(read_xlsx("dataEdited/rnaSequencing/rna_seq_prep.xlsx") %>%
                     select(extractionID, sequencingID, metatranscriptomeID)) %>%
  select(metatranscriptomeID, sequencingID, sampleID,
         startDate, depth, replicate, volumeFiltered) %>%
  arrange(metatranscriptomeID) %>%
  as.data.frame()


#### Save out naming file to rename metatranscriptomes ####
all.data %>%
  select(sequencingID, metatranscriptomeID) %>%
  write.table("dataRaw/metatranscriptomes/reports/MT_naming_key.tsv",
              sep = "\t",
              col.names = FALSE,
              row.names = FALSE,
              quote = FALSE)


#### Save out metadata for analysis ####
all.data %>%
  select(-sequencingID) %>%
  write.csv("metadata/metatranscriptome_metadata.csv",
            row.names = FALSE)



#### Check md5 readouts ####
pre.transfer.md5 <- read.table("dataRaw/metatranscriptomes/reports/md5sum.txt") 
names(pre.transfer.md5) <- c("md5_pre", "fileName")
pre.transfer.md5 <- pre.transfer.md5 %>%
  mutate(fileName = gsub("Peterson/", "", fileName)) %>%
  filter(grepl("KMBP020", fileName))

post.transfer.md5 <- read.table("dataRaw/metatranscriptomes/reports/md5sum_post_download.txt") 
names(post.transfer.md5) <- c("md5_post", "fileName")

full_join(pre.transfer.md5,
          post.transfer.md5) %>%
  mutate(match = (md5_post == md5_pre)) %>%
  select(match) %>%
  unlist(use.names = FALSE) %>%
  all()
# All the md5 values match
rm(pre.transfer.md5, post.transfer.md5)



#### SRA entry prep ####
biosample.ids <- read.table("dataEdited/metagenomes/NCBI_upload/NCBI_info_MIMS_with_BioSample.txt",
                            sep = '\t',
                            skip = 1) %>%
  select(V1, V3)
names(biosample.ids) <- c("biosample_accession", "sample_name")
ncbi.data <- all.data %>%
  mutate(
    # Use descriptive name for sample name
    sample_name = paste("MendotaWaterColumn", "_",
                        year(startDate), "_",
                        month(startDate, label = TRUE, abbr = TRUE), "_",
                        as.character(depth), "m",
                        sep = "")) %>%
  left_join(biosample.ids) %>%
  rename(library_ID = metatranscriptomeID) %>%
  mutate(title = paste("Metatranscriptomic sequencing of microbial community from the water column of Lake Mendota on ",
                       startDate, " at ", depth, "m - biological replicate ", replicate,
                       sep = ""),
         library_strategy = "RNA-Seq",
         library_source = "METATRANSCRIPTOMIC",
         library_selection = "Inverse rRNA",
         library_layout = "paired",
         platform = "ILLUMINA",
         instrument_model = "Illumina NovaSeq 6000",
         design_description = "Library prep and sequencing done at QB3 Genomics at the University of California - Berkeley. Library prep used Kapa Biosystem Library Prep kit with a target insert length of ∼600 bp. Sequenced using the S4 chemistry. 150 bp paired-end sequencing",
         filetype = "fastq",
         filename = paste(library_ID, "_R1.fastq.gz",
                          sep = ""),
         filename2 = paste(library_ID, "_R2.fastq.gz",
                           sep = ""),
         filename3 = "",
         filename4 = "",
         assembly = "",
         fasta_file = "") %>%
  select(biosample_accession, library_ID, title, library_strategy, library_source,
         library_selection, library_layout, platform, instrument_model, design_description,
         filetype, filename, filename2, filename3, filename4, assembly, fasta_file) %>%
  arrange(biosample_accession)

write.csv(ncbi.data,
          file = "dataRaw/metatranscriptomes/temp_NCBI_info_SRA.csv",
          quote = FALSE,
          row.names = FALSE)
