#### code/metagenomes/metagenome_metadata.R ####
# Benjamin D. Peterson

# This file contains the code to generate the combined
# metadata for the metagenome sequencing and assembly.


#### Always start with a clean slate, ya bum ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP/")
library(lubridate)
library(readxl)
library(tidyverse)


#### Generate table of MG sample site information ####
sample.data <- read.csv("metadata/processedMetadata/filter_metadata.csv") %>%
  select(filterID, sampleID, startDate, depth, volumeFiltered)
all.data <- read_xlsx("dataEdited/dnaExtractions/DNA_extractions.xlsx") %>%
                        filter(!(filterID %in% c("NA", "blank"))) %>%
  left_join(sample.data) %>%
  # Read in metagenome metadata
  right_join(rbind(read_xlsx("dataEdited/dnaSequencing/2020/samplePrep/KMBP010_dilutions.xlsx") %>%
                     select(extractionID, metagenomeID),
                   read_xlsx("dataEdited/dnaSequencing/2021/samplePrep/BLI21_MG_dilutions.xlsx") %>%
                     rename(metagenomeID = `original metagenomeID`) %>%
                     select(extractionID, metagenomeID))) %>%
  select(metagenomeID, sampleID, startDate, depth, volumeFiltered) %>%
  as.data.frame()



#### Save out metadata for analysis ####
all.data %>%
  write.csv("metadata/metagenome_metadata.csv",
            row.names = FALSE)



#### BioSample entry prep ####
sampleIDs.of.interest <- unique(all.data$sampleID)
sample.data %>%
  filter(sampleID %in% sampleIDs.of.interest) %>%
  mutate(
    # Use descriptive name for sample name
    sample_name = paste("MendotaWaterColumn", "_",
                        year(startDate), "_",
                        month(startDate, label = TRUE, abbr = TRUE), "_",
                        as.character(depth), "m",
                        sep = ""),
    # Set sample title as sample ID from project
    sample_title = sampleID,
    bioproject_accession = "PRJNA876614",
    # "Organism" name selected from this list: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Tree&id=410657&lvl=3&lin=f&keep=1&srchmode=1&unlock
    # Instructions to do so here: https://www.ncbi.nlm.nih.gov/biosample/docs/organism/#metagenomes
    organism = "freshwater metagenome",
    collection_date = startDate,
    # env_broad_scale selected from this list: https://ontobee.org/ontology/ENVO?iri=http://purl.obolibrary.org/obo/ENVO_00000428
    # Instructions to do so here: https://www.gensc.org/pages/standards/all-terms.html
    env_broad_scale = "lake [ENVO:00000020]",
    # Same as above
    env_local_scale = "eutrophic lake [ENVO:01000548]|lake with an anoxic hypolimnion [ENVO:01001073]|freshwater lake [ENVO:00000021]",
    # Same as above
    env_medium = "fresh water [ENVO_00002011]",
    geo_loc_name = "USA: Lake Mendota - Deep Hole",
    lat_lon = "43.0989 N 89.4055 W") %>%
  select(sample_name, sample_title, bioproject_accession, organism, collection_date,
         depth, env_broad_scale, env_local_scale, env_medium, geo_loc_name, lat_lon) %>%
  write.csv(file = "dataEdited/metagenomes/NCBI_upload/temp_NCBI_info_MIMS.csv",
              quote = FALSE,
              row.names = FALSE)


#### SRA entry prep ####
biosample.ids <- read.table("dataEdited/metagenomes/NCBI_upload/NCBI_info_MIMS_with_BioSample.txt",
                            sep = '\t',
                            skip = 1) %>%
  select(V1, V3)
names(biosample.ids) <- c("biosample_accession", "sample_name")
all.data %>%
  mutate(
    # Use descriptive name for sample name
    sample_name = paste("MendotaWaterColumn", "_",
                        year(startDate), "_",
                        month(startDate, label = TRUE, abbr = TRUE), "_",
                        as.character(depth), "m",
                        sep = "")) %>%
  left_join(biosample.ids) %>%
  rename(library_ID = metagenomeID) %>%
  mutate(replicate = duplicated(sample_name) + 1,
         title = paste("Metagenomic sequencing of microbial community from the water column of Lake Mendota on ",
                       startDate, " at ", depth, "m - biological replicate ", replicate,
                       sep = ""),
         library_strategy = "WGS",
         library_source = "METAGENOMIC",
         library_selection = "RANDOM",
         library_layout = "paired",
         platform = "ILLUMINA",
         instrument_model = "Illumina NovaSeq 6000",
         design_description = "Library prep and sequencing done at QB3 Genomics Center at Berkeley. Library prep used Kapa Biosystem Library Prep kit with a target insert length of ∼600 bp. Sequenced using the S4 chemistry. 150 bp paired-end sequencing",
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
  arrange(biosample_accession) %>%
  write.csv(file = "dataEdited/metagenomes/NCBI_upload/temp_NCBI_info_SRA.csv",
            quote = FALSE,
            row.names = FALSE)
