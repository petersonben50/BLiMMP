#### code/mash/NMDS_mash_distance.R ####
# Benjamin D. Peterson

#### Always start with a clean slate ####
rm(list = ls())
setwd("~/Documents/research/BLiMMP/")
library(lubridate)
library(tidyverse)
library(vegan)


#### Read in distance data ####

dist.data <- read.table("dataEdited/metagenomes/BLI_MG_sketches.dist",
                        stringsAsFactors = FALSE,
                        sep = "\t")
names(dist.data) <- c("refID", "queryID", "mash_dist",
                      "pvalue", "matching_hashes")
clean.dist.data <- dist.data %>%
  mutate(refID = refID %>%
           strsplit("/metagenomes/") %>% sapply("[", 2) %>%
           gsub("_R1.fastq.gz", "", .)) %>%
  mutate(queryID = queryID %>%
           strsplit("/metagenomes/") %>% sapply("[", 2) %>%
           gsub(".fastq.gz", "", .)) %>%
  select(refID, queryID, mash_dist) %>%
  spread(key = queryID,
         value = mash_dist)


#### Read in metadata ####
MG.data <- read.csv("metadata/metagenome_metadata.csv",
                    stringsAsFactors = FALSE)
MG.data.vector <- paste(
  MG.data$metagenomeID, "\n",
  MG.data$startDate, "\n",
  MG.data$depth, "m",
  sep = "")
names(MG.data.vector) <- MG.data$metagenomeID



#### Generate matrix ####

clean.dist.data.matrix <- clean.dist.data %>%
  select(-refID)
row.names(clean.dist.data.matrix) <- clean.dist.data$refID


#### Test dimensions ####
# https://ourcodingclub.github.io/tutorials/ordination/
par(mfrow = c(1,1))
NMDS.screen <- function(x) { #where x is the name of the data frame variable
  plot(rep(1, 10),
       replicate(10, metaMDS(x, autotransform = F, k = 1)$stress),
       xlim = c(1, 10),
       ylim = c(0, 0.30),
       xlab = "# of Dimensions",
       ylab = "Stress",
       main = "NMDS stress plot")
  for (i in 1:10) {
    points(rep(i + 1,10),
           replicate(10, metaMDS(x, autotransform = F, k = i + 1)$stress))
  }
}
x <- clean.dist.data.matrix
NMDS.screen(clean.dist.data.matrix)





#### Make plot ####

plot.nmds <- function(distance.matrix = clean.dist.data.matrix,
                      dimensions = 4) {
  
  # Run NMDS ordination
  nmds.sites <- metaMDS(distance.matrix,
                        k = dimensions)
  
  # Set up dataframe for plotting
  data.scores = cbind(NMDS1 = nmds.sites$points[, "MDS1"],
                      NMDS2 = nmds.sites$points[, "MDS2"],
                      site.info = MG.data.vector[rownames(nmds.sites$points)]) %>%
    as.data.frame() 
  data.scores$NMDS1 = as.numeric(data.scores$NMDS1)
  data.scores$NMDS2 = as.numeric(data.scores$NMDS2)
  # Plot data
  ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
    geom_point(size = 0) + 
    theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
          axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
          axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2)) + 
    labs(x = "NMDS1", y = "NMDS2") +
    geom_text(aes(label = site.info))
  
}


plot.nmds(clean.dist.data.matrix,
          dimensions = 3)


