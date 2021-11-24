#### code/incubations/incubation_ancillary_data.R ####
# Written for BLiMMP project
# Benjamin D. Peterson


#### Prep workspace ####
rm(list = ls())
setwd("~/Documents/research/BLiMMP/")
library(ggpubr)
library(tidyverse)
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")



#### Read in sulfide data from incubations, combine with metadata ####
sulfide.data.MA <- read.csv("dataEdited/waterChemistry/sulfide/MA_data.csv",
                            stringsAsFactors = FALSE)



#### Read in water column sulfide data ####
sulfide.data.WC <- read.csv("dataEdited/waterChemistry/sulfide/WC_data.csv",
                            stringsAsFactors = FALSE) %>%
  select(sampleID, S_conc_uM) %>%
  rename(original_S_conc = S_conc_uM) %>%
  group_by(sampleID) %>%
  summarise(original_S_conc = mean(original_S_conc))



#### Join sulfide data ####

sulfide.data.all <- left_join(sulfide.data.MA,
                              sulfide.data.WC) %>%
  mutate(delta_S_uM = (S_conc_uM - original_S_conc)) %>%
  select(tripID, incubationID, startDate, depth, filtered, treatment, incubationTimePoint, durationInDays, S_conc_uM, delta_S_uM, original_S_conc) %>%
  filter(incubationTimePoint != "t3") %>%
  group_by(tripID, incubationID, startDate, depth, filtered, treatment, incubationTimePoint, durationInDays) %>%
  summarize(S_conc_uM = mean(S_conc_uM),
            delta_S_uM = mean(delta_S_uM),
            original_S_conc = mean(original_S_conc)) %>%
  ungroup()

# rm(sulfide.data.WC,
#    sulfide.data.MA)

#### Treatment color vector ####
color.vector <- c(cb.translator["vermillion"],
                  cb.translator["bluishgreen"],
                  cb.translator["black"])
names(color.vector) <- c("filtered-unamended",
                         "unfiltered-unamended",
                         "unfiltered-molybdate")


#### Function to plot sulfide levels across incubation ####
start.date.of.interest <- "2021-09-10"
depth.of.interest <- 19.9
max.concentration <- 150
plot.sulfide.levels <- function(start.date.of.interest,
                                depth.of.interest,
                                max.concentration,
                                min.concentration = 0,
                                legend.position.to.use = c(0.8, 0.2),
                                color.vector.to.use = color.vector,
                                number.of.days = 2) {
  
  sulfide.data.all %>%
    filter(startDate == start.date.of.interest,
           depth == depth.of.interest) %>%
    arrange(durationInDays) %>%
    ggplot(aes(x = durationInDays,
               y = S_conc_uM,
               color = treatment,
               group = incubationID)) +
    geom_point() +
    geom_line() +
    scale_fill_manual(values = color.vector.to.use) +
    scale_color_manual(values = color.vector.to.use) +
    xlab("Timepoint (days)") +
    ylab("Sulfide (µM)") +
    xlim(c(0, number.of.days)) +
    ylim(c(min.concentration, max.concentration)) +
    theme_bw() +
    theme(legend.position = legend.position.to.use) +
    geom_hline(yintercept = sulfide.data.all %>%
                 filter(startDate == start.date.of.interest,
                        depth == depth.of.interest) %>%
                 select(original_S_conc) %>%
                 unlist(use.names = FALSE) %>% mean(),
               colour = cb.translator["reddishpurple"])
  
}



#### Plot sulfide in September incubations ####

sept.upper <- plot.sulfide.levels(start.date.of.interest = "2020-09-02",
                                  depth.of.interest = 11,
                                  max.concentration = 20,
                                  legend.position.to.use = c(0.6, 0.6),
                                  number.of.days = 4)

sept.mid <- plot.sulfide.levels(start.date.of.interest = "2020-09-02",
                                depth.of.interest = 15.5,
                                max.concentration = 50,
                                min.concentration = -5,
                                legend.position.to.use = "none",
                                number.of.days = 4)

sept.bottom <- plot.sulfide.levels(start.date.of.interest = "2020-09-02",
                                   depth.of.interest = 20.7,
                                   max.concentration = 150,
                                   legend.position.to.use = "none",
                                   number.of.days = 4)
pdf("results/incubations/sulfide_levels/2020-09_sulfide_in_incubations.pdf",
    height = 10,
    width = 6)
ggarrange(sept.upper, sept.mid, sept.bottom,
          ncol = 1)
dev.off()



#### Plot sulfide in October incubations ####

oct.upper <- plot.sulfide.levels(start.date.of.interest = "2020-10-10",
                                  depth.of.interest = 15.7,
                                  max.concentration = 20,
                                  legend.position.to.use = c(0.6, 0.6),
                                  number.of.days = 4)

oct.bottom <- plot.sulfide.levels(start.date.of.interest = "2020-10-10",
                                depth.of.interest = 20.9,
                                max.concentration = 120,
                                min.concentration = -5,
                                legend.position.to.use = "none",
                                number.of.days = 4)

pdf("results/incubations/sulfide_levels/2020-10_sulfide_in_incubations.pdf",
    height = 7,
    width = 6)
ggarrange(oct.upper, oct.bottom,
          ncol = 1)
dev.off()





#### Plot sulfide in September 2021 incubations ####

sept.upper <- plot.sulfide.levels(start.date.of.interest = "2021-09-10",
                                  depth.of.interest = 10.8,
                                  max.concentration = 10,
                                  legend.position.to.use = c(0.6, 0.6))

sept.mid <- plot.sulfide.levels(start.date.of.interest = "2021-09-10",
                                depth.of.interest = 11.9,
                                max.concentration = 75,
                                legend.position.to.use = "none")

sept.bottom <- plot.sulfide.levels(start.date.of.interest = "2021-09-10",
                                   depth.of.interest = 19.9,
                                   max.concentration = 150,
                                   legend.position.to.use = "none")
pdf("results/incubations/sulfide_levels/2021-09_sulfide_in_incubations.pdf",
    height = 10,
    width = 6)
ggarrange(sept.upper, sept.mid, sept.bottom,
          ncol = 1)
dev.off()



#### Plot sulfide in October incubations ####
oct.2021.upper <- plot.sulfide.levels(start.date.of.interest = "2021-10-14",
                                       depth.of.interest = 14.2,
                                       max.concentration = 10,
                                       legend.position.to.use = c(0.8, 0.6))

oct.2021.mid <- plot.sulfide.levels(start.date.of.interest = "2021-10-14",
                                     depth.of.interest = 15.2,
                                     max.concentration = 100,
                                     legend.position.to.use = "none")

oct.2021.bottom <- plot.sulfide.levels(start.date.of.interest = "2021-10-14",
                                        depth.of.interest = 19.9,
                                        max.concentration = 200,
                                        legend.position.to.use = "none")

pdf("results/incubations/sulfide_levels/2021-10_sulfide_in_incubations.pdf",
    height = 10,
    width = 6)
ggarrange(oct.2021.upper, oct.2021.mid, oct.2021.bottom,
          ncol = 1)

dev.off()

# 
# #### Plot difference in sulfide for specified trip ####
# 
# trip.of.interest <- "BLiMMP_trip_006"
# depth.of.interest <- 16.18
# 
# plot.S.difference <- function(trip.of.interest,
#                               depth.of.interest) {
#   
#   sulfide.data.trip <- sulfide.data.all %>%
#     filter(tripID == trip.of.interest) %>%
#     filter(depth == depth.of.interest) %>%
#     filter(filtered == "no")
#   
#   stripchart(delta_S_uM ~ amendment,
#             data = sulfide.data.trip,
#             ylab = "Change in sulfide (µM)",
#             ylim = c(-60, 0),
#             xlab = "Amendment",
#             vertical = TRUE,
#             method = "jitter",
#             pch = 18)
# 
# }
# 
# par(mfrow = c(2, 1),
#     mar = c(4.5, 4.5, 0.5, 0.5))
# plot.S.difference("BLiMMP_trip_006",
#                   16.18)
# plot.S.difference("BLiMMP_trip_006",
#                   19.3)