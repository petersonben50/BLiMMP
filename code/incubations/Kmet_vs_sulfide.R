#### code/incubations/Kmet_vs_sulfide.R ####
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




#### Read in sulfide data ####
sulfide.data <- read.csv("dataEdited/waterColumn_sulfide_data.csv") %>%
  group_by(sampleID) %>%
  summarize(S_conc_uM = mean(S_conc_uM))


#### Read in data ####
Hg.Kmet.data <- read.csv("dataEdited/incubation_Hg_rate_data.csv") %>%
  filter(treatment == "unfiltered-unamended",
         year(startDate) %in% c(2020, 2021))


#### Combine data ####
all.data <- left_join(Hg.Kmet.data,
                      sulfide.data)


#### Set vectors ####
color.vector.year <- cb.translator[c("blue", "vermillion")]
names(color.vector.year) <- c("2020", "2021")
shape.vector <- c(16, 17)
names(shape.vector) <- c("September", "October")



#### Plot the Kmet vs. sulfide for t1 ####
Kmet_vs_sulfide_plot_t1 <- all.data %>%
  ggplot(aes(x = S_conc_uM,
             y = Kmet_t1,
             color = as.character(year(startDate)),
             shape = as.character(month(startDate,
                                        label = TRUE, abbr = FALSE)))) +
  geom_point(size = 3) +
  scale_color_manual(values = color.vector.year,
                     name = "Year") +
  scale_shape_manual(values = shape.vector,
                     name = "Month") +
  ylab(bquote(K[met] (day^-1))) +
  xlab("Sulfide (µM)") +
  theme_bw() +
  theme(legend.position = c(0.25, 0.72),
        legend.spacing.y = unit(0, "mm"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black",
                                             fill = "white"))


#### Plot the Kmet vs. sulfide for t2 ####
Kmet_vs_sulfide_plot_t2 <- all.data %>%
  ggplot(aes(x = S_conc_uM,
             y = Kmet_t2,
             color = as.character(year(startDate)),
             shape = as.character(month(startDate,
                                        label = TRUE, abbr = FALSE)))) +
  geom_point() +
  scale_color_manual(values = color.vector.year,
                     name = "Year") +
  scale_shape_manual(values = shape.vector,
                     name = "Month") +
  ylab("Kmet (day^-1)") +
  xlab("Sulfide (µM)") +
  theme_bw() +
  theme(legend.position = c(0.25, 0.72))


#### Plot the Kmet vs. sulfide for total incubations ####
Kmet_vs_sulfide_plot_total <- all.data %>%
  ggplot(aes(x = S_conc_uM,
             y = Kmet_total,
             color = as.character(year(startDate)),
             shape = as.character(month(startDate,
                                        label = TRUE, abbr = FALSE)))) +
  geom_point() +
  scale_color_manual(values = color.vector.year,
                     name = "Year") +
  scale_shape_manual(values = shape.vector,
                     name = "Month") +
  ylab("Kmet (day^-1)") +
  xlab("Sulfide (µM)") +
  theme_bw() +
  theme(legend.position = c(0.25, 0.72))


#### View all plots ####
ggarrange(Kmet_vs_sulfide_plot_t1 + labs(title = "Kmet t1"), Kmet_vs_sulfide_plot_t2 + labs(title = "Kmet t2"),
          Kmet_vs_sulfide_plot_total + labs(title = "Kmet total"),
          ncol = 1)

#### Save out plot ####
pdf("results/incubations/Kmet_vs_sulfide.pdf",
    width = 5,
    height = 4)
Kmet_vs_sulfide_plot_t1
dev.off()
