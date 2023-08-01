#### code/figures/incubations_sulfide.R ####
# Written for BLiMMP project
# Benjamin D. Peterson


#### Prep workspace ####
rm(list = ls())
setwd("~/Documents/research/BLiMMP/")
library(ggpubr)
library(tidyverse)
cb.translator <- readRDS("references/colorblind_friendly_colors.rds")



#### Read in sulfide data from incubations, combine with metadata ####
sulfide_data_MA <- read.csv("dataFinal/incubation_sulfide_data.csv",
                            stringsAsFactors = FALSE)



#### Read in water column sulfide data ####
sulfide_data_WC <- read.csv("dataFinal/water_chem_data.csv",
                            stringsAsFactors = FALSE) %>%
  rename(startDate = date) %>%
  group_by(startDate, depth) %>%
  summarize(original_sulfide_uM = mean(sulfide_uM, na.rm = TRUE))


#### Join sulfide data ####
sulfide_data <- left_join(sulfide_data_MA,
                          sulfide_data_WC) %>%
  mutate(delta_sulfide_uM = (S_conc_uM - original_sulfide_uM)) %>%
  select(tripID, incubationID, startDate, depth, filtered, treatment, incubationTimePoint, durationInDays, S_conc_uM, delta_sulfide_uM, original_sulfide_uM) %>%
  filter(incubationTimePoint != "t3") %>%
  group_by(tripID, incubationID, startDate, depth, filtered, treatment, incubationTimePoint, durationInDays) %>%
  summarize(S_conc_uM = mean(S_conc_uM),
            delta_sulfide_uM = mean(delta_sulfide_uM),
            original_sulfide_uM = mean(original_sulfide_uM)) %>%
  ungroup()

rm(sulfide_data_MA,
   sulfide_data_WC)

#### Treatment color vector ####
color.vector <- c(cb.translator["vermillion"],
                  cb.translator["bluishgreen"],
                  cb.translator["black"])
names(color.vector) <- c("filtered-unamended",
                         "unfiltered-unamended",
                         "unfiltered-molybdate")
naming.vector <- c("Filtered control", "Ambient conditions", "Molybdate-inhibited")


#### Keep only incubations of interest ####
sulfide_data <- sulfide_data %>%
  filter(treatment %in% names(color.vector))


#### Function to plot sulfide levels across incubation ####
# start.date.of.interest <- "2021-09-10"
# depth.of.interest <- 19.9
# max.concentration <- 150
plot.sulfide.levels <- function(start.date.of.interest,
                                depth.of.interest,
                                max.concentration,
                                min.concentration = 0,
                                legend.position.to.use = c(0.8, 0.2),
                                color.vector.to.use = color.vector,
                                number.of.days = 2,
                                conc.range = c(0, 100),
                                xlab.to.use = NULL,
                                ylab.to.use = NULL,
                                start.date.to.use,
                                title.to.use = "default") {
  
  plot.to.show <- sulfide_data %>%
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
    xlim(c(0, number.of.days)) +
    ylim(c(min.concentration, max.concentration)) +
    geom_hline(yintercept = sulfide_data %>%
                 filter(startDate == start.date.of.interest,
                        depth == depth.of.interest) %>%
                 select(original_sulfide_uM) %>%
                 unlist(use.names = FALSE) %>% mean(),
               colour = cb.translator["reddishpurple"]) +
    theme_bw() +
    theme(legend.position = legend.position.to.use,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black"),
          axis.title.x = element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.text.x =  element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          axis.ticks.x = element_line(colour = "black"),
          axis.ticks.y = element_line(colour = "black"),
          axis.line = element_line(colour = "black"))

  if (is.null(xlab.to.use)) {
    plot.to.show <- plot.to.show +
      theme(axis.title.x = element_blank())
  } else {
    plot.to.show <- plot.to.show +
      xlab(xlab.to.use)
  }
  
  if (is.null(ylab.to.use)) {
    plot.to.show <- plot.to.show +
      theme(axis.title.y = element_blank())
  } else {
    plot.to.show <- plot.to.show +
      ylab(ylab.to.use)
  }
  
  if (title.to.use == "default") {
    plot.to.show <- plot.to.show +
      labs(subtitle = paste(start.date.of.interest, ": ", depth.of.interest, "m",
                            sep = ""))
  } else if (is.null(title.to.use)) {
    plot.to.show <- plot.to.show
  } else {
    plot.to.show <- plot.to.show +
      labs(title = title.to.use)
  }
  
}



#### Plots: Generate sulfide plots for all locations ####
sept.2020.upper <- plot.sulfide.levels(start.date.of.interest = "2020-09-02",
                                       depth.of.interest = 11,
                                       max.concentration = 20,
                                       legend.position.to.use = "none",
                                       number.of.days = 4,
                                       ylab.to.use = "Sulfide (µM)",
                                       xlab.to.use = "Days")

sept.2020.mid <- plot.sulfide.levels(start.date.of.interest = "2020-09-02",
                                     depth.of.interest = 15.5,
                                     max.concentration = 50,
                                     min.concentration = -5,
                                     legend.position.to.use = "none",
                                     number.of.days = 4,
                                     xlab.to.use = "Days")

sept.2020.bottom <- plot.sulfide.levels(start.date.of.interest = "2020-09-02",
                                        depth.of.interest = 20.7,
                                        max.concentration = 150,
                                        legend.position.to.use = "none",
                                        number.of.days = 4,
                                        xlab.to.use = "Days")


oct.2020.upper <- plot.sulfide.levels(start.date.of.interest = "2020-10-10",
                                      depth.of.interest = 15.7,
                                      max.concentration = 20,
                                      legend.position.to.use = "none",
                                      number.of.days = 4,
                                      ylab.to.use = "Sulfide (µM)",
                                      xlab.to.use = "Days")
empty.plot <- sept.2020.bottom + theme_void() +
  theme(legend.position = c(0.5, 0.5)) +
  scale_color_manual(values = color.vector,
                     labels = naming.vector) +
  ylim(c(1000, 1001)) +
  labs(subtitle = "")
oct.2020.bottom <- plot.sulfide.levels(start.date.of.interest = "2020-10-10",
                                       depth.of.interest = 20.9,
                                       max.concentration = 120,
                                       min.concentration = -5,
                                       legend.position.to.use = "none",
                                       number.of.days = 4,
                                       xlab.to.use = "Days")


sept.2021.upper <- plot.sulfide.levels(start.date.of.interest = "2021-09-10",
                                       depth.of.interest = 10.8,
                                       max.concentration = 10,
                                       legend.position.to.use = "none",
                                       ylab.to.use = "Sulfide (µM)",
                                       xlab.to.use = "Days")

sept.2021.mid <- plot.sulfide.levels(start.date.of.interest = "2021-09-10",
                                     depth.of.interest = 11.9,
                                     max.concentration = 50,
                                     legend.position.to.use = "none",
                                     xlab.to.use = "Days")

sept.2021.bottom <- plot.sulfide.levels(start.date.of.interest = "2021-09-10",
                                        depth.of.interest = 19.9,
                                        max.concentration = 150,
                                        legend.position.to.use = "none",
                                        xlab.to.use = "Days")

oct.2021.upper <- plot.sulfide.levels(start.date.of.interest = "2021-10-14",
                                      depth.of.interest = 14.2,
                                      max.concentration = 10,
                                      legend.position.to.use = "none",
                                      ylab.to.use = "Sulfide (µM)",
                                      xlab.to.use = "Days")

oct.2021.mid <- plot.sulfide.levels(start.date.of.interest = "2021-10-14",
                                    depth.of.interest = 15.2,
                                    max.concentration = 100,
                                    legend.position.to.use = "none",
                                    xlab.to.use = "Days")

oct.2021.bottom <- plot.sulfide.levels(start.date.of.interest = "2021-10-14",
                                       depth.of.interest = 19.9,
                                       max.concentration = 200,
                                       legend.position.to.use = "none",
                                       xlab.to.use = "Days")



#### Save out plots as PDF ####
pdf("results/figures/incubations_sulfide.pdf",
    width = 7.2,
    height = 8)
ggarrange(sept.2020.upper, sept.2020.mid, sept.2020.bottom,
          oct.2020.upper, empty.plot, oct.2020.bottom,
          sept.2021.upper, sept.2021.mid, sept.2021.bottom,
          oct.2021.upper, oct.2021.mid, oct.2021.bottom,
          ncol = 3, nrow = 4)
dev.off()

