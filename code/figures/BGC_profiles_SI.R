#### code/figures/BGC_profiles_main.R ####
# Benjamin D. Peterson

#### Start with a clean slate ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP")
source("code/figures/BGC_profiles_SI_plotting_functions.R")


#### Read in data ####
exo.data <- read.csv("dataFinal/exo_data.csv")
waterChem.data <- read.csv("dataFinal/water_chem_data.csv")


#### Plotting functions tailored to manuscript figure ####
# testing
exo.data.to.use = exo.data
dates_of_sampling = c("2020-09-02",
                      "2020-10-10",
                      "2021-09-10",
                      "2021-10-14")
legend.location = "topright"


#### Set up layout for figure ####
vert_space <- 0.15
xlabel_space <- 0.08
header_vertArea <- c(0.920, 0.990)
row_1_vertArea <- c(header_vertArea[1] - vert_space,
                    header_vertArea[1])
row_2_vertArea <- c(row_1_vertArea[1] - xlabel_space - vert_space,
                    row_1_vertArea[1] - xlabel_space)
row_3_vertArea <- c(row_2_vertArea[1] - xlabel_space - vert_space,
                    row_2_vertArea[1] - xlabel_space)
row_4_vertArea <- c(row_3_vertArea[1] - xlabel_space - vert_space,
                    row_3_vertArea[1] - xlabel_space)

increment <- (0.1125)
horiArea1_lab <- c(0.00,0.08)
horiArea2_phy <- c(horiArea1_lab[2], horiArea1_lab[2] + increment)
horiArea3_met <- c(horiArea2_phy[2], horiArea2_phy[2] + increment)
horiArea4_sul <- c(horiArea3_met[2], horiArea3_met[2] + increment)
horiArea5_DOC <- c(horiArea4_sul[2], horiArea4_sul[2] + increment)
horiArea6_FHG <- c(horiArea5_DOC[2], horiArea5_DOC[2] + increment)
horiArea7_PHG <- c(horiArea6_FHG[2], horiArea6_FHG[2] + increment)
horiArea8_KHG <- c(horiArea7_PHG[2], horiArea7_PHG[2] + increment)
horiArea9_DGM <- c(horiArea8_KHG[2], horiArea8_KHG[2] + increment)


#### Save out pdf of plots ####
grDevices::cairo_pdf("results/figures/BGC_profiles_SI.pdf",
                     family = "Arial",
                     width = 7.20,
                     height = 8.64)
# Add the headers as separate fields on the plot
split.screen(rbind(c(horiArea1_lab, row_1_vertArea), c(horiArea1_lab, row_2_vertArea), c(horiArea1_lab, row_3_vertArea), c(horiArea1_lab, row_4_vertArea),
                   c(horiArea2_phy, header_vertArea), c(horiArea2_phy, row_1_vertArea), c(horiArea2_phy, row_2_vertArea), c(horiArea2_phy, row_3_vertArea), c(horiArea2_phy, row_4_vertArea),
                   c(horiArea3_met, header_vertArea), c(horiArea3_met, row_1_vertArea), c(horiArea3_met, row_2_vertArea), c(horiArea3_met, row_3_vertArea), c(horiArea3_met, row_4_vertArea),
                   c(horiArea4_sul, header_vertArea), c(horiArea4_sul, row_1_vertArea), c(horiArea4_sul, row_2_vertArea), c(horiArea4_sul, row_3_vertArea), c(horiArea4_sul, row_4_vertArea),
                   c(horiArea5_DOC, header_vertArea), c(horiArea5_DOC, row_1_vertArea), c(horiArea5_DOC, row_2_vertArea), c(horiArea5_DOC, row_3_vertArea), c(horiArea5_DOC, row_4_vertArea),
                   c(horiArea6_FHG, header_vertArea), c(horiArea6_FHG, row_1_vertArea), c(horiArea6_FHG, row_2_vertArea), c(horiArea6_FHG, row_3_vertArea), c(horiArea6_FHG, row_4_vertArea),
                   c(horiArea7_PHG, header_vertArea), c(horiArea7_PHG, row_1_vertArea), c(horiArea7_PHG, row_2_vertArea), c(horiArea7_PHG, row_3_vertArea), c(horiArea7_PHG, row_4_vertArea),
                   c(horiArea8_KHG, header_vertArea), c(horiArea8_KHG, row_1_vertArea), c(horiArea8_KHG, row_2_vertArea), c(horiArea8_KHG, row_3_vertArea), c(horiArea8_KHG, row_4_vertArea),
                   c(horiArea9_DGM, header_vertArea), c(horiArea9_DGM, row_1_vertArea), c(horiArea9_DGM, row_2_vertArea), c(horiArea9_DGM, row_3_vertArea), c(horiArea9_DGM, row_4_vertArea)))

i <- 1
date_labels <- c(expression("Sep. 2"^"nd"*" 2020"),
                 expression("Oct. 10"^"th"*" 2020"),
                 expression("Sep. 10"^"th"*" 2021"),
                 expression("Oct. 14"^"th"*" 2021"))
names(date_labels) <- dates_of_sampling
for (date_to_use in dates_of_sampling) {
  screen(i)
  text(x = 0.2, y = 0.25, adj = 0,
       labels = date_labels[date_to_use],
       srt = 90,
       cex = 0.8)
  i <<- i + 1
}

# Physical parameters
for (date_to_use in dates_of_sampling) {
  if (date_to_use == dates_of_sampling[1]) {
    screen(i)
    plot.exo.data(date.of.sampling = date_to_use,plot_legend_only = "YES")
    i <- i + 1
    screen(i)
    plot.exo.data(date.of.sampling = date_to_use)
    i <- i + 1
  } else {
    screen(i)
    plot.exo.data(date.of.sampling = date_to_use)
    i <- i + 1
  }
}

# Metal species
for (date_to_use in dates_of_sampling) {
  if (date_to_use == dates_of_sampling[1]) {
    screen(i)
    plot_metal_data(date.of.sampling = date_to_use,
                    plot_legend_only = "YES")
    i <- i + 1
    screen(i)
    plot_metal_data(date.of.sampling = date_to_use)
    i <- i + 1
  } else {
    screen(i)
    plot_metal_data(date.of.sampling = date_to_use)
    i <- i + 1
  }
}

# Sulfur species
for (date_to_use in dates_of_sampling) {
  if (date_to_use == dates_of_sampling[1]) {
    screen(i)
    plot_sulfur_data(plot_legend_only = "YES")
    i <- i + 1
    }
  screen(i)
  plot_sulfur_data(date.of.sampling = date_to_use)
  i <- i + 1
}

# Carbon/particulates profiles
for (date_to_use in dates_of_sampling) {
  if (date_to_use == dates_of_sampling[1]) {
    screen(i)
    plot_C_profile(plot_legend_only = "YES")
    i <- i + 1
    }
  screen(i)
  plot_C_profile(date.of.sampling = date_to_use)
  i <- i + 1
  }

# Filtered mercury
for (date_to_use in dates_of_sampling) {
  if (date_to_use == dates_of_sampling[1]) {
    screen(i)
    plot_filt_Hg_profile(plot_legend_only = "YES")
    i <- i + 1
  }
  screen(i)
  plot_filt_Hg_profile(date.of.sampling = date_to_use)
  i <- i + 1
}

# Particulate mercury
for (date_to_use in dates_of_sampling) {
  if (date_to_use == dates_of_sampling[1]) {
    screen(i)
    plot_part_Hg_profile(plot_legend_only = "YES")
    i <- i + 1
  }
  screen(i)
  plot_part_Hg_profile(date.of.sampling = date_to_use)
  i <- i + 1
}

# Particulate partitioning mercury
for (date_to_use in dates_of_sampling) {
  if (date_to_use == dates_of_sampling[1]) {
    screen(i)
    plot_part_part_Hg_profile(plot_legend_only = "YES")
    i <- i + 1
    }
  screen(i)
  plot_part_part_Hg_profile(date.of.sampling = date_to_use)
  i <- i + 1
  }


# DGM profiles
for (date_to_use in dates_of_sampling) {
  if (date_to_use %in% dates_of_sampling[1]) {
    screen(i)
    plot_DGM_profile(plot_legend_only = "YES")
    i <- i + 1
  }
  if (date_to_use %in% dates_of_sampling[3:4]) {
    screen(i)
    plot_DGM_profile(date.of.sampling = date_to_use)
  } else {
    screen(i)
    par(mar = c(0, 0, 0, 0))
    box(which = "plot")
  }
  i <- i + 1
}

dev.off()
