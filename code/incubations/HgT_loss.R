#### code/incubations/HgT_loss.R ####
# Written for BLiMMP project
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP")
library(ggpubr)
library(lubridate)
library(readxl)
library(tidyverse)
cb.translator <- readRDS("references/colorblind_friendly_colors.rds")



#### Read in data ####
HgT.data <- read.csv("dataEdited/incubation_Hg_rate_data.csv") %>%
  select(-c(Kmet_t1, Kmet_t2, Kmet_total,
            Kdem_t1, Kdem_t2, Kdem_total)) %>%
  filter(year(startDate) %in% c(2020, 2021))



#### Treatment renaming vector ####
color.vector <- c(cb.translator["vermillion"],
                  cb.translator["bluishgreen"],
                  cb.translator["black"])
names(color.vector) <- c("filtered-unamended",
                         "unfiltered-unamended",
                         "unfiltered-molybdate")
naming.vector <- c("Filtered control", "Ambient conditions", "Molybdate-inhibited")
point.vector.to.use <- c(16, 4)
names(point.vector.to.use) <- c(TRUE, FALSE)



#### Function: Plot of HgT loss ####
HgT.loss <- function(HgT.data.to.use = HgT.data,
                     date.to.use = "2020-09-02",
                     depth.to.use = 11,
                     constituents = c("HgT_198_daily_percent_loss_t1", "HgT_198_daily_percent_loss_t2"),
                     conc.range = c(0, 100),
                     xlab.to.use = NULL,
                     x.ticks,
                     ylab.to.use = NULL) {
  plotting.data <- HgT.data.to.use %>%
    filter(startDate == date.to.use,
           depth == depth.to.use) %>%
    mutate(treatment = fct_relevel(treatment, names(color.vector)))
  plotting.data <- plotting.data[, c("startDate", "depth", "treatment", constituents)]
  plotting.data <- plotting.data %>%
    gather(key = constituent,
           value = concentration,
           all_of(constituents))
  
  plot.to.show <- plotting.data %>%
    ggplot(aes(x = constituent,
               y = concentration,
               color = treatment)) +
    geom_boxplot() +
    scale_color_manual(values = color.vector) +
    ylim(conc.range) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(subtitle = paste(date.to.use, "-", depth.to.use, "m"))
  
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
  
  if (!is.null(x.ticks)) {
    plot.to.show <- plot.to.show +
      scale_x_discrete(labels = x.ticks)
    
  }
  
  plot.to.show
}




#### Function: All plots for a given constituent ####
make.all.plots.for.constituent <- function(constituent.to.use = c("HgT_198_daily_percent_loss_t1",
                                                               "HgT_198_daily_percent_loss_t2"),
                                           concentration.to.use = c(-40, 62),
                                           constituent.ylab = NULL,
                                           x.ticks.to.use = NULL) {
  
  sept.2020.11m <- HgT.loss(HgT.data.to.use = HgT.data, date.to.use = "2020-09-02", depth.to.use = 11,
                            constituent = constituent.to.use, conc.range = concentration.to.use,
                            ylab.to.use = constituent.ylab, x.ticks = x.ticks.to.use)
  sept.2020.15m <- HgT.loss(HgT.data.to.use = HgT.data, date.to.use = "2020-09-02", depth.to.use = 15.5,
                            constituent = constituent.to.use, conc.range = concentration.to.use,
                            ylab.to.use = constituent.ylab, x.ticks = x.ticks.to.use)
  sept.2020.20m <- HgT.loss(HgT.data.to.use = HgT.data, date.to.use = "2020-09-02", depth.to.use = 20.7,
                            constituent = constituent.to.use, conc.range = concentration.to.use,
                            ylab.to.use = constituent.ylab, x.ticks = x.ticks.to.use)
  
  oct.2020.15m <- HgT.loss(HgT.data.to.use = HgT.data, date.to.use = "2020-10-10", depth.to.use = 15.7,
                           constituent = constituent.to.use, conc.range = concentration.to.use,
                           ylab.to.use = constituent.ylab, x.ticks = x.ticks.to.use)
  
  empty.plot <- oct.2020.15m + theme_void() +
    theme(legend.position = c(0.5, 0.5)) +
    scale_color_manual(values = color.vector,
                       labels = naming.vector) +
    ylim(c(1000, 1001)) +
    labs(subtitle = "")
  
  oct.2020.20m <- HgT.loss(HgT.data.to.use = HgT.data, date.to.use = "2020-10-10", depth.to.use = 20.9,
                           constituent = constituent.to.use, conc.range = concentration.to.use,
                           ylab.to.use = constituent.ylab, x.ticks = x.ticks.to.use)
  
  sept.2021.10m <- HgT.loss(HgT.data.to.use = HgT.data, date.to.use = "2021-09-10", depth.to.use = 10.8,
                            constituent = constituent.to.use, conc.range = concentration.to.use,
                            ylab.to.use = constituent.ylab, x.ticks = x.ticks.to.use)
  sept.2021.11m <- HgT.loss(HgT.data.to.use = HgT.data, date.to.use = "2021-09-10", depth.to.use = 11.9,
                            constituent = constituent.to.use, conc.range = concentration.to.use,
                            ylab.to.use = constituent.ylab, x.ticks = x.ticks.to.use)
  sept.2021.19m <- HgT.loss(HgT.data.to.use = HgT.data, date.to.use = "2021-09-10", depth.to.use = 19.9,
                            constituent = constituent.to.use, conc.range = concentration.to.use,
                            ylab.to.use = constituent.ylab, x.ticks = x.ticks.to.use)
  
  oct.2021.14m <- HgT.loss(HgT.data.to.use = HgT.data, date.to.use = "2021-10-14", depth.to.use = 14.2,
                           constituent = constituent.to.use, conc.range = concentration.to.use,
                           ylab.to.use = constituent.ylab, x.ticks = x.ticks.to.use)
  oct.2021.15m <- HgT.loss(HgT.data.to.use = HgT.data, date.to.use = "2021-10-14", depth.to.use = 15.2,
                           constituent = constituent.to.use, conc.range = concentration.to.use,
                           ylab.to.use = constituent.ylab, x.ticks = x.ticks.to.use)
  oct.2021.19m <- HgT.loss(HgT.data.to.use = HgT.data, date.to.use = "2021-10-14", depth.to.use = 19.9,
                           constituent = constituent.to.use, conc.range = concentration.to.use,
                           ylab.to.use = constituent.ylab, x.ticks = x.ticks.to.use)
  
  ggarrange(sept.2020.11m, sept.2020.15m, sept.2020.20m,
            oct.2020.15m, empty.plot, oct.2020.20m,
            sept.2021.10m, sept.2021.11m, sept.2021.19m,
            oct.2021.14m, oct.2021.15m, oct.2021.19m,
            nrow = 4, ncol = 3)  
}




#### Save out PDF with HgT loss for enriched isotopes ####
pdf("results/incubations/HgT_loss_daily.pdf",
    width = 7.2,
    height = 10)
make.all.plots.for.constituent(constituent.to.use = c("HgT_198_daily_percent_loss_t1",
                                                      "HgT_198_daily_percent_loss_t2"),
                               concentration.to.use = c(-40, 62),
                               constituent.ylab = "198HgT loss (%)",
                               x.ticks.to.use = c("t1", "t2"))
make.all.plots.for.constituent(constituent.to.use = c("HgT_204_daily_percent_loss_t1",
                                                      "HgT_204_daily_percent_loss_t2"),
                               concentration.to.use = c(-40, 50),
                               constituent.ylab = "204HgT loss (%)",
                               x.ticks.to.use = c("t1", "t2"))
dev.off()


pdf("results/incubations/HgT_loss_total.pdf",
    width = 7.2,
    height = 10)
make.all.plots.for.constituent(constituent.to.use = c("HgT_198_percent_loss_total",
                                                      "HgT_204_percent_loss_total"),
                               concentration.to.use = c(-16, 62),
                               constituent.ylab = "HgT loss (%)",
                               x.ticks.to.use = c("198HgT", "204HgT"))
dev.off()
