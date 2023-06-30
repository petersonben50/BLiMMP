#### code/figures/BGC_profiles_main.R ####
# Benjamin D. Peterson

#### Start with a clean slate ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP")
# 
source("code/figures/BGC_profiles_main_plotting_functions.R")


#### Read in data ####
incubation.depths <- read.csv("metadata/processedMetadata/incubation_metadata.csv") %>%
  select(startDate, depth) %>%
  unique()
exo.data <- read.csv("dataFinal/exo_data.csv")
waterChem.data <- read.csv("dataFinal/water_chem_data.csv")


#### Plotting functions tailored to manuscript figure ####
# testing
exo.data.to.use = exo.data
date.of.sampling = "2021-09-10"
sampling.depths.of.interest <- incubation.depths %>%
  filter(startDate == date.of.sampling) %>%
  select(depth) %>%
  unlist(use.names = FALSE)
points.of.sampling = sampling.depths.of.interest
legend.location = "topright"


#### Set up layout for figure ####
header_top_vertArea <- c(0.97, 1.0)
header_bottom_vertArea <- c(0.47, 0.50)
header_left_horiArea <- c(0.05, 0.50)
header_right_horiArea <- c(0.56, 1.0)

top_row_vertArea <- c(0.66, 0.965)
bottom_row_vertArea <- c(0.16, 0.465)
phys_plot_horiArea1 <- c(0.05,0.19)
redox_plot_horiArea1 <- c(0.19, 0.33)
hg_plot_horiArea1 <- c(0.33, 0.47)
column_2_adjustment <- 0.98 - hg_plot_horiArea1[2]

#### Save out pdf of plots ####
grDevices::cairo_pdf("results/figures/BGC_profiles_main.pdf",
                     width = 7.2,
                     height = 7.5)
# Add the headers as separate fields on the plot
split.screen(rbind(c(header_left_horiArea, header_top_vertArea), c(header_right_horiArea, header_top_vertArea), c(header_left_horiArea, header_bottom_vertArea), c(header_right_horiArea, header_bottom_vertArea),
                   c(phys_plot_horiArea1, top_row_vertArea), c(redox_plot_horiArea1, top_row_vertArea), c(hg_plot_horiArea1, top_row_vertArea),
                   c(phys_plot_horiArea1 + column_2_adjustment, top_row_vertArea), c(redox_plot_horiArea1 + column_2_adjustment, top_row_vertArea), c(hg_plot_horiArea1 + column_2_adjustment, top_row_vertArea),
                   c(phys_plot_horiArea1, bottom_row_vertArea), c(redox_plot_horiArea1, bottom_row_vertArea), c(hg_plot_horiArea1, bottom_row_vertArea),
                   c(phys_plot_horiArea1 + column_2_adjustment, bottom_row_vertArea), c(redox_plot_horiArea1 + column_2_adjustment, bottom_row_vertArea), c(hg_plot_horiArea1 + column_2_adjustment, bottom_row_vertArea)))

screen(1)
text(x = -0.035, y = 0.4, adj = c(0, 0.5),
     expression("A. September 2"^"nd"*", 2020"))
screen(2)
text(x = -0.035, y = 0.4, adj = c(0, 0.5),
     expression("B. October 10"^"th"*", 2020"))
screen(3)
text(x = -0.035, y = 0.4, adj = c(0, 0.5),
     expression("C. September 10"^"th"*", 2021"))
screen(4)
text(x = -0.035, y = 0.4, adj = c(0, 0.5),
     expression("D. October 14"^"th"*", 2021"))

screen(5)
plot.exo.data(date.of.sampling = "2020-09-02",
              plot_legend = "yes")
screen(6)
plot_redox_data(date.of.sampling = "2020-09-02",
                legend_location = "top")
screen(7)
plot.Hg.profile(date.of.sampling = "2020-09-02",
                legend_location = "top")

screen(8)
plot.exo.data(date.of.sampling = "2020-10-10")
screen(9)
plot_redox_data(date.of.sampling = "2020-10-10")
screen(10)
plot.Hg.profile(date.of.sampling = "2020-10-10")

screen(11)
plot.exo.data(date.of.sampling = "2021-09-10")
screen(12)
plot_redox_data(date.of.sampling = "2021-09-10")
screen(13)
plot.Hg.profile(date.of.sampling = "2021-09-10")

screen(14)
plot.exo.data(date.of.sampling = "2021-10-14")
screen(15)
plot_redox_data(date.of.sampling = "2021-10-14")
screen(16)
plot.Hg.profile(date.of.sampling = "2021-10-14")

dev.off()
