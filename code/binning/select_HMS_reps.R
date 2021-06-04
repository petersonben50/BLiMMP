#### code/binning/aggregate_bin_data.R ####
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("~/Documents/research/BLiMMP")
library(igraph)
library(tidyverse)
library(vegan)



#### Calculate HMSs ####

ANI.values <- read.table(file = "dataEdited/binning/ANI/goodBins.all.ani.out.cleaned",
                         sep = "\t",
                         header = FALSE,
                         stringsAsFactors = FALSE,
                         col.names = c("GENOME1", "GENOME2",
                                       "ANI_1", "ANI_2",
                                       "AF_1", "AF_2")) %>%
  mutate(GENOME1 = gsub(".fna", "", GENOME1),
         GENOME2 = gsub(".fna", "", GENOME2))

# Check plot of values
plot(y = ANI.values$ANI_1,
     x = ANI.values$AF_1,
     xlim = c(0.2, 1.0),
     ylim = c(70, 100),
     pch = 18)
points(y = ANI.values$ANI_2,
       x = ANI.values$AF_2,
       pch = 18)
abline(h = 97)
abline(v = 0.5)


# Set cut-offs
ANI.cutoff <- 97
cov.cutoff <- 0.5

# Generate edgelist
edgelist <- ANI.values %>%
  filter((ANI_1 > ANI.cutoff & AF_1 > cov.cutoff) |
           (ANI_2 > ANI.cutoff & AF_2 > cov.cutoff)) %>%
  select(GENOME1, GENOME2)
# Generate the graph from the edgelist
adjacency.graph <- graph_from_data_frame(edgelist)
# Store the bin and HMS name in a df.
HMS.ID <- data.frame(paste("HMS.",
                           clusters(adjacency.graph)$membership,
                           sep = ""),
                     names(clusters(adjacency.graph)$membership),
                     stringsAsFactors = FALSE)
names(HMS.ID) <- c("HMS", "binID")


# Check out loner bins
bin.list <- unique(c(ANI.values$GENOME1, ANI.values$GENOME2)) %>%
  strsplit(".fna") %>%
  sapply("[", 1)
lone.bins <- data.frame(bin.list[!(bin.list %in% HMS.ID$bin)],
                        bin.list[!(bin.list %in% HMS.ID$bin)],
                        stringsAsFactors = FALSE)
names(lone.bins) <- c("HMS", "binID")
# Add loner bins to list
HMS.bin.info <- do.call("rbind",
                        list(HMS.ID, lone.bins))
rm(ANI.cutoff, cov.cutoff, HMS.ID, lone.bins,
   edgelist, ANI.values, adjacency.graph,
   bin.list)


#### Read in taxonomy data ####

taxonomy.data <- read.table("dataEdited/binning/taxonomy/gtdbtk.bac120.summary.tsv",
                            sep = "\t",
                            header = TRUE,
                            stringsAsFactors = FALSE) %>%
  rename(binID = user_genome,
         gtdb_tax = classification) %>%
  select(binID, gtdb_tax)
allData <- HMS.bin.info %>%
  left_join(taxonomy.data)
rm(taxonomy.data)


#### CheckM data ####

checkm.data <- read.csv("dataEdited/binning/quality/good_bins_data.txt",
                        stringsAsFactors = FALSE,
                        header = FALSE)
names(checkm.data) <- c("binID",
                        paste("checkM_",
                              c("completeness", "contamination", "strain_het",
                                "genome_length", "scaf_count", "N50", "mean_scaf_length",
                                "longest_scaf", "GC", "predicted_genes"),
                              sep = ""))
allData <- allData %>%
  left_join(checkm.data)
rm(checkm.data)


#### anvi'o completeness/redundancy ####

anvio.quality <- read.table("dataEdited/binning/quality/bins_summary_all.txt",
                            sep = "\t",
                            stringsAsFactors = FALSE,
                            header = FALSE) %>%
  select("V1", "V7", "V8")
names(anvio.quality) <- c("binID", "anvio_completeness", "anvio_redundancy")
allData <- allData %>%
  left_join(anvio.quality)
rm(anvio.quality)


#### Save out data ####
write.csv(allData,
          "dataEdited/binning/bin_data.csv",
          row.names = FALSE,
          quote = FALSE)


#### Save out list of good bins ####
allData %>%
  filter(!is.na(checkM_completeness)) %>%
  filter(binID != "BLI20_hgcA_031") %>%
  select(binID) %>%
  unlist(use.names = FALSE) %>%
  writeLines("dataEdited/binning/final_bins.txt")
