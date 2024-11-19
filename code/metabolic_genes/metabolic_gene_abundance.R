#### code/metabolic_analyses/metabolic_gene_abundance_manuscript.R ####
# Benjamin D. Peterson


#### Clean up on aisle R ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP")
library(patchwork)
library(tidyverse)
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")
source("code/BLiMMP_functions.R")


#### Read in data ####
gene.depth <- readRDS("dataEdited/metabolic_analyses/depth/metabolicProtein_depth_clean.rds")


#### Set up for plots ####
par(mar=c(9,3,2.5,1),
    mgp=c(1.5,0.4,0),
    tck=-0.008)

#### September ####
color.vector <- c(cb.translator["vermillion"],
                  cb.translator["orange"])
names(color.vector) <- c("narG", "napA")
N.genes.september <- plot.profile.for.multiple.genes(marker.depth.df = gene.depth,
                                                     monthOfInterest = 9,
                                                     genesOfInterest = names(color.vector),
                                                     show.mean.coverage = FALSE,
                                                     depth_limits = c(25, 0),
                                                     coverage_limits = c(0, 55),
                                                     color.vector.to.use = color.vector,
                                                     legend.position.to.use = c(0.8, 0.8),
                                                     xlab.to.use = "Gene abundance (per 100X SCG)",
                                                     titleToUse = "A.")

color.vector <- c(cb.translator["blue"],
                  cb.translator["bluishgreen"])
names(color.vector) <- c("dsrA", "dsrD")
S.genes.september <- plot.profile.for.multiple.genes(marker.depth.df = gene.depth,
                                                     monthOfInterest = 9,
                                                     genesOfInterest = names(color.vector),
                                                     show.mean.coverage = FALSE,
                                                     depth_limits = c(25, 0),
                                                     coverage_limits = c(0, 2),
                                                     color.vector.to.use = color.vector,
                                                     legend.position.to.use = c(0.8, 0.8),
                                                     xlab.to.use = "Gene abundance (per 100X SCG)",
                                                     titleToUse = "C.")
N.genes.september + S.genes.september



#### October ####
color.vector <- c(cb.translator["vermillion"],
                  cb.translator["orange"])
names(color.vector) <- c("narG", "napA")
N.genes.october <- plot.profile.for.multiple.genes(marker.depth.df = gene.depth,
                                                   monthOfInterest = 10,
                                                   genesOfInterest = names(color.vector),
                                                   show.mean.coverage = FALSE,
                                                   depth_limits = c(25, 0),
                                                   coverage_limits = c(0, 55),
                                                   color.vector.to.use = color.vector,
                                                   legend.position.to.use = c(0.8, 0.8),
                                                   xlab.to.use = "Gene abundance (per 100X SCG)",
                                                   titleToUse = "B.")

color.vector <- c(cb.translator["blue"],
                  cb.translator["bluishgreen"])
names(color.vector) <- c("dsrA", "dsrD")
S.genes.october <- plot.profile.for.multiple.genes(marker.depth.df = gene.depth,
                                                   monthOfInterest = 10,
                                                   genesOfInterest = names(color.vector),
                                                   show.mean.coverage = FALSE,
                                                   depth_limits = c(25, 0),
                                                   coverage_limits = c(0, 2),
                                                   color.vector.to.use = color.vector,
                                                   legend.position.to.use = c(0.8, 0.8),
                                                   xlab.to.use = "Gene abundance (per 100X SCG)",
                                                   titleToUse = "D.")
pdf("results/metabolic_analyses/N_S_reduction.pdf",
    width = 7,
    height = 6)
(N.genes.september + S.genes.september) / (N.genes.october + S.genes.october)
dev.off()


#### EET genes ####
color.vector <- c(cb.translator["orange"],
                  cb.translator["reddishpurple"],
                  cb.translator["yellow"])
names(color.vector) <- c("pcc_porin", "OMP_MtrB_PioB", "BBOMP")
E.genes.september <- plot.profile.for.multiple.genes(marker.depth.df = gene.depth,
                                                   monthOfInterest = 9,
                                                   genesOfInterest = names(color.vector),
                                                   show.mean.coverage = FALSE,
                                                   depth_limits = c(25, 0),
                                                   coverage_limits = c(0, 80),
                                                   color.vector.to.use = color.vector,
                                                   legend.position.to.use = c(0.8, 0.8))
E.genes.october <- plot.profile.for.multiple.genes(marker.depth.df = gene.depth,
                                                   monthOfInterest = 10,
                                                   genesOfInterest = names(color.vector),
                                                   show.mean.coverage = FALSE,
                                                   depth_limits = c(25, 0),
                                                   coverage_limits = c(0, 50),
                                                   color.vector.to.use = color.vector,
                                                   legend.position.to.use = c(0.8, 0.8))
E.genes.september + E.genes.october
