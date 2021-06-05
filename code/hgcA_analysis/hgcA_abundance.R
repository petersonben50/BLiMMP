#### code/hgcA_analysis/hgcA_abundance.R ####
# Written by Benjamin D. Peterson


#### Get ready. Get set ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP/")
library(gridExtra)
library(lubridate)
library(readxl)
library(tidyverse)
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")
gray.vector <- c("gray50", "gray30", "gray80")
names(gray.vector) <- c("gray50", "gray30", "gray80")
cb.translator <- c(cb.translator, gray.vector)


#### Read in hgcA classification ####
tax.data <- read_xlsx("dataEdited/hgcA_analysis/hgcA_taxonomy_assignment.xlsx") %>%
  select(seqID, manual_classification_for_colors) %>%
  rename(hgcA_ID = seqID)


#### Read in color information ####
hgcA.manual.taxonomy <- read_xlsx("dataEdited/hgcA_analysis/hgcA_taxonomy_assignment.xlsx",
                                  sheet = "colors_to_use")
color.vector.to.use <- cb.translator[hgcA.manual.taxonomy$colorsToUse]
names(color.vector.to.use) <- hgcA.manual.taxonomy$classification


#### Read in abundance data ####
coverage.data <- read.csv("dataEdited/hgcA_analysis/depth/hgcA_coverage.csv") %>%
  select(-volumeFiltered)


#### Combine data ####
all.data <- full_join(tax.data,
                      coverage.data) %>%
  mutate(depth.date = paste(month(startDate, label = TRUE), "-", depth, "m",
                            sep = ""),
         depth.date = fct_relevel(depth.date,
                                  c("Sep-11m", "Sep-15.5m", "Sep-20.7m",
                                    "Oct-15.7m", "Oct-20.9m")),
         manual_classification_for_colors = fct_relevel(manual_classification_for_colors,
                                                        names(color.vector.to.use)))


#### hgcA abundance as points ####
all.data %>%
  group_by(depth.date,
           manual_classification_for_colors) %>%
  summarize(coverage = sum(coverage)) %>%
  ggplot(aes(x = depth.date,
             y = coverage,
             color = manual_classification_for_colors)) +
  geom_point(position = position_dodge(width = .4),
             alpha = 0.8,
             size = 2.5) +
  scale_color_manual(values=color.vector.to.use) +
  theme_classic()


#### hgcA abundance as stacked bars ####
all.data %>%
  mutate(depth.date = paste(startDate, "\n", depth, " m",
                            sep = ""),
         manual_classification_for_colors = fct_relevel(manual_classification_for_colors,
                                                        names(color.vector.to.use))) %>%
  ggplot(aes(x = depth.date,
             y = coverage,
             fill = manual_classification_for_colors)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=color.vector.to.use) +
  theme_classic()



#### hgcA total abundance: September ####
hgcA.September <- all.data %>%
  filter(month(startDate) == 9) %>%
  mutate(depth.date = fct_relevel(depth.date,
                                  rev( c("Sep-11m", "Sep-15.5m", "Sep-20.7m",
                                         "Oct-15.7m", "Oct-20.9m")))) %>%
  group_by(depth.date) %>%
  summarise(coverage = sum(coverage)) %>%
  ggplot(aes(x = depth.date,
             y = coverage)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  coord_flip() +
  xlab(element_blank()) +
  ylab("Total hgcA coverage") +
  ylim(c(0, 20)) +
  theme(legend.position = "none")


#### hgcA total abundance: October ####
hgcA.October <- all.data %>%
  filter(month(startDate) == 10) %>%
  mutate(depth.date = fct_relevel(depth.date,
                                  rev( c("Sep-11m", "Sep-15.5m", "Sep-20.7m",
                                         "Oct-15.7m", "Oct-20.9m")))) %>%
  group_by(depth.date) %>%
  summarise(coverage = sum(coverage)) %>%
  ggplot(aes(x = depth.date,
             y = coverage)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  coord_flip() +
  xlab(element_blank()) +
  ylab("Total hgcA coverage") +
  ylim(c(0, 20)) +
  theme(legend.position = "none")



#### hgcA taxonomy abundance as side-by-side bars: September ####
hgcA.tax.September <- all.data %>%
  filter(month(startDate) == 9) %>%
  mutate(depth.date = fct_relevel(depth.date,
                                  rev( c("Sep-11m", "Sep-15.5m", "Sep-20.7m",
                                         "Oct-15.7m", "Oct-20.9m")))) %>%
  ggplot(aes(x = depth.date,
             y = coverage,
             fill = manual_classification_for_colors)) +
  geom_bar(stat = "identity",position = "dodge") +
  scale_fill_manual(values=color.vector.to.use) +
  theme_classic() +
  coord_flip() +
  xlab(element_blank()) +
  ylab("hgcA coverage") +
  ylim(c(0, 7)) +
  theme(legend.position = "none")


#### hgcA taxonomy abundance as side-by-side bars: October ####
hgcA.tax.October <- all.data %>%
  filter(month(startDate) == 10) %>%
  mutate(depth.date = fct_relevel(depth.date,
                                  rev( c("Sep-11m", "Sep-15.5m", "Sep-20.7m",
                                         "Oct-15.7m", "Oct-20.9m")))) %>%
  ggplot(aes(x = depth.date,
             y = coverage,
             fill = manual_classification_for_colors)) +
  geom_bar(stat = "identity",
           position = "dodge") +
  scale_fill_manual(values=color.vector.to.use,
                    guide = guide_legend(reverse = TRUE)) +
  xlab(element_blank()) +
  ylab("hgcA coverage") +
  ylim(c(0, 7)) +
  theme_classic() +
  coord_flip() +
  theme(legend.position = c(0.8, 0.5),
        legend.title = element_blank())


grid.arrange(hgcA.September, hgcA.tax.September,
             hgcA.October, hgcA.tax.October,
             heights = c(3,2),
             widths = c(1, 2))
