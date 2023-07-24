#### code/figures/temp_DO_heatmaps.R ####
# Benjamin D. Peterson


#### Set up shop ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP")
library(akima)
library(lattice)
library(readxl)
library(tidyverse)
library(viridis)



#### Read in data ####
temperature_data <- read.csv("dataFinal/temperature_buoy_data.csv")
DO_data <- read.csv("dataFinal/DO_profiles_data.csv") %>%
  mutate(sampleDate = as.Date(sampleDate))
pigment_data <- read.csv("dataFinal/pigment_profiles_data.csv")
sampling_location_info <- read_xlsx("dataEdited/incubations/incubation_sites_notes.xlsx")
cb_translator <- readRDS("references/colorblind_friendly_colors.rds")


#### Set variables ####
temp_range_to_plot <- c(5, 31)
temp_color_object <- viridis(20)
DO_range_to_plot <- c(0, 18)
DO_color_object <- colorRampPalette(c(cb_translator["black"] , cb_translator["blue"], cb_translator["bluishgreen"], cb_translator["orange"], cb_translator["orange"]),
                                    bias = 1, space = "rgb")
axis_tick_label_size <- 0.7
axis_label_size <- 0.8


#### Function to plot heatmap of temperature ####
heatmap_of_temp <- function(year_to_plot = 2020,
                            start_date = "2020-05-01",
                            end_date = "2020-11-15") {
  
  par(mar = c(0, 0, 0, 0),
      tck = -0.008)
  
  #### Prepare data ####
  data_to_use = temperature_data %>%
    filter(year(sampledate) == year_to_plot) %>%
    mutate(hour = str_pad(hour, 4, pad = "0")) %>%
    mutate(date_time = ymd_hm(paste(sampledate, hour))) %>%
    filter(as.Date(date_time) >= as.Date(start_date)) %>%
    filter(as.Date(date_time) <= as.Date(end_date)) %>%
    mutate(yday_fraction = yday(date_time) + hour(date_time)/24) %>%
    select(yday_fraction, depth, wtemp) %>%
    spread(key = depth,
           value = wtemp) %>%
    mutate(`24` = `20`) %>%
    arrange(yday_fraction) %>%
    gather(value = wtemp,
           key = depth,
           -1)
  #### Generated needed vectors ####
  dates_desired <- seq(yday(start_date),
                       yday(end_date),
                       by = 0.5)
  depths_desired <- seq(0, 24, by = 0.25)
  dates_we_have <- unique(data_to_use$yday_fraction)
  data_to_use <- data_to_use %>%
    spread(value = wtemp,
           key = yday_fraction)
  #### Interpolate data twice ####
  transformed_data <- sapply(data_to_use[, -1],
                             function(row_day) {
                               approx(data_to_use$depth,
                                      y = row_day,
                                      xout = depths_desired)$y
                             })
  retransformed_data <- sapply(t(transformed_data) %>% as.data.frame(),
                               function(depth_to_use) {
                                 approx(dates_we_have,
                                        y = depth_to_use,
                                        xout = dates_desired)$y
                               })
  #### Generate image ####
  image(z = retransformed_data,
        x = dates_desired,
        y = depths_desired,
        ylim = c(24, 0),
        zlim = temp_range_to_plot,
        col = temp_color_object,
        xaxt = 'n', yaxt = 'n',
        xlab = "", ylab = "")
  #### Calculate thermocline depth ####
  data_for_thermo_df <- temperature_data %>%
    filter(year(sampledate) == year_to_plot) %>%
    mutate(hour = str_pad(hour, 4, pad = "0")) %>%
    mutate(date_time = ymd_hm(paste(sampledate, hour))) %>%
    filter(as.Date(date_time) >= as.Date(start_date)) %>%
    filter(as.Date(date_time) <= as.Date(end_date)) %>%
    select(date_time, wtemp, depth) %>%
    spread(key = depth,
           value = wtemp) %>%
    rename(datetime = date_time)
    
  names(data_for_thermo_df)[-1] <- paste("var_", names(data_for_thermo_df)[-1], sep = "")
  thermocline_df <- data_for_thermo_df %>%
    rLakeAnalyzer::ts.thermo.depth() %>%
    group_by(as.Date(datetime)) %>%
    summarize(thermo.depth = median(thermo.depth,
                                    na.rm = TRUE)) %>%
    rename(datetime = `as.Date(datetime)`)
  lines(x = yday(thermocline_df$datetime),
        y = thermocline_df$thermo.depth)
  
  #### Add axis scales ####
  par(mgp = c(1.5, 0.2, 0))
  axis(2,
       at = seq(0, 25, by = 5),
       labels = seq(0, 25, by = 5),
       cex.axis = axis_tick_label_size,
       las = 1)
  par(mgp = c(1.5, -0.1, 0))
  dates_to_plot_on_axes <- seq(from = as.Date(start_date),
                               to = as.Date(end_date),
                               by = 'months') 
  axis(1,
       at = yday(dates_to_plot_on_axes),
       labels = month(dates_to_plot_on_axes,label = TRUE,abbr = TRUE),
       cex.axis = axis_tick_label_size,
       gap.axis = -1,
       padj = 0)
  
  #### Add axis labels ####
  mtext("Depth (m)",
        side = 2,
        line = 0.85,
        cex = axis_label_size)
  
}


#### Function to plot chlorophyll a data ####
pigment_plot <- function(start_date = "2021-05-01",
                         end_date = "2021-11-15") {
  
  par(mar = c(0, 0, 0, 0),
      tck = -0.008)
  
  pigment_data_to_plot <- pigment_data %>%
    filter(sampleDate > start_date,
           sampleDate < end_date)
  
  #### Plot data ####
  plot(x = yday(pigment_data_to_plot$sampleDate),
       y = pigment_data_to_plot$mean_phyco_rfu,
       col = "blue",
       pch = 16,
       ylim = c(-0.5, 3),
       xlim = c(yday(start_date), yday(end_date)),
       cex = 0.5,
       xaxt = 'n', yaxt = 'n',
       xaxs = 'i', yaxs = 'i',
       xlab = "", ylab = "")
  lines(x = yday(pigment_data_to_plot$sampleDate),
        y = pigment_data_to_plot$mean_phyco_rfu,
        col = "blue",
        pch = 16)
  points(x = yday(pigment_data_to_plot$sampleDate),
         y = pigment_data_to_plot$mean_chlor_rfu,
         col = "green",
         cex = 0.5,
         pch = 16)
  lines(x = yday(pigment_data_to_plot$sampleDate),
        y = pigment_data_to_plot$mean_chlor_rfu,
        col = "green",
        pch = 16)
  
  #### Add axis scales ####
  par(mgp = c(1.5, 0.2, 0))
  axis(2,
       at = seq(0, 3, by = 1),
       labels = seq(0, 3, by = 1),
       cex.axis = axis_tick_label_size,
       las = 1)
  par(mgp = c(1.5, -0.1, 0))
  dates_to_plot_on_axes <- seq(from = as.Date(start_date),
                               to = as.Date(end_date),
                               by = 'months') 
  axis(1,
       at = yday(dates_to_plot_on_axes),
       labels = month(dates_to_plot_on_axes,label = TRUE,abbr = TRUE),
       cex.axis = axis_tick_label_size,
       gap.axis = -1,
       padj = 0)
  
  #### Add axis labels ####
  mtext("Relative\nfluoresence units",
        side = 2,
        line = 0.65,
        cex = axis_label_size*0.8)
}


#### Function to plot heatmap of DO ####
# DO.profile.data = DO.data,
# year.of.interest,
# date.limits = c("2020-05-01", "2020-08-01")
DO_heatmap <- function(year_to_plot = 2021,
                       start_date = "2021-05-01",
                       end_date = "2021-11-15") {
  
  par(mar = c(0, 0, 0, 0),
      tck = -0.008)
  
  #### Generate a matrix of the DO values ####
  DO_data_year <- DO_data %>%
    filter(year(sampleDate) == year_to_plot) %>%
    arrange(sampleDate)
  dates_we_have <- unique(DO_data_year$sampleDate)
  
  DO_data_year <- DO_data_year %>%
    spread(key = sampleDate,
           value = DO)

  
  #### Interpolate data twice ####
  dates_desired <- seq(yday(start_date), yday(end_date), by = 1)
  depths_desired <- seq(0, 20, by = 0.25)
  transformed_data <- sapply(DO_data_year[, -1],
                             function(col_date) {
                               approx(DO_data_year$depth,
                                      y = col_date,
                                      xout = depths_desired)$y
                             })
  retransformed_data <- sapply(t(transformed_data) %>% as.data.frame(),
                               function(col_depth) {
                                 approx(yday(dates_we_have),
                                        y = col_depth,
                                        xout = dates_desired)$y
                               })
  
  #### Plotting function ####
  image(z = retransformed_data,
        x = dates_desired,
        y = depths_desired,
        ylim = c(20, 0),
        xlim = c(yday(start_date),
                 yday(end_date)),
        zlim = DO_range_to_plot,
        col = DO_color_object(20),
        xaxt = 'n', yaxt = 'n',
        xlab = "", ylab = "")
  
  #### Add axis scales ####
  par(mgp = c(1.5, 0.2, 0))
  axis(2,
       at = seq(0, 20, by = 5),
       labels = seq(0, 20, by = 5),
       cex.axis = axis_tick_label_size,
       las = 1)
  par(mgp = c(1.5, -0.1, 0))
  dates_to_plot_on_axes <- seq(from = as.Date(start_date),
                               to = as.Date(end_date),
                               by = 'months') 
  axis(1,
       at = yday(dates_to_plot_on_axes),
       labels = month(dates_to_plot_on_axes,label = TRUE,abbr = TRUE),
       cex.axis = axis_tick_label_size,
       gap.axis = -1,
       padj = 0)
  
  #### Add axis labels ####
  mtext("Depth (m)",
        side = 2,
        line = 0.85,
        cex = axis_label_size)
}


#### Set up plot areas ####

header_top_vertArea <- c(0.97, 1.0)
header_middle_vertArea <- c(0.59, 0.62)
header_bottom_vertArea <- c(0.33, 0.36)

header_left_horiArea <- c(0.005, 0.21)
header_right_horiArea <- c(0.565, 0.77)

top_row_vertArea <- c(0.67, 0.96)
middle_row_vertArea <- c(0.41, 0.58)
bottom_row_vertArea <- c(0.03, 0.32)

left_horiArea <- c(0.06, 0.43)
center_horiArea <- c(0.45, 0.56)
right_horiArea <- c(0.62, 0.99)


#### Generate PDF ####
cairo_pdf("results/figures/temp_DO_heatmaps.pdf",
          width = 7,
          height = 7)
split.screen(rbind(c(left_horiArea, top_row_vertArea), c(center_horiArea, top_row_vertArea), c(right_horiArea, top_row_vertArea),
                   c(left_horiArea, middle_row_vertArea), c(center_horiArea, middle_row_vertArea), c(right_horiArea, middle_row_vertArea),
                   c(left_horiArea, bottom_row_vertArea), c(center_horiArea, bottom_row_vertArea), c(right_horiArea, bottom_row_vertArea),
                   c(header_left_horiArea, header_top_vertArea),c(header_right_horiArea, header_top_vertArea),
                   c(header_left_horiArea, header_middle_vertArea), c(header_right_horiArea, header_middle_vertArea),
                   c(header_left_horiArea, header_bottom_vertArea), c(header_right_horiArea, header_bottom_vertArea)))

#### Generate temperature plots ####
screen(1)
heatmap_of_temp(2020,
                start_date = "2020-05-01",
                end_date = "2020-11-10")
screen(2)
legend_image <- as.raster(matrix(viridis(20)[20:1], ncol=1))
par(mar = c(0, 0, 0, 0))
plot(c(0,3),
     c(0,1),
     type = 'n',
     axes = F,
     xlab = '',
     ylab = '',
     main = '') 
text(x = 1.5,
     y = seq(0, 0.8, l = 6),
     labels = seq(temp_range_to_plot[1], temp_range_to_plot[2]-1, l = 6),
     cex = axis_tick_label_size)
text(x = -0.4, y = 0.9,
     pos = 4,
     labels = "Temp (˚C)",
     cex = axis_tick_label_size + 0.1)
rasterImage(legend_image, 0, 0, 1,0.8)


screen(3)
heatmap_of_temp(2021,
                start_date = "2021-05-01",
                end_date = "2021-11-10")

#### Generate pigment plots ####
screen(4)
pigment_plot(start_date = "2020-05-01",
             end_date = "2020-11-15")
screen(5)
par(mar = c(0, 0, 0, 0))
plot(c(0,3), c(0,1),
     type = 'n', axes = F,
     xlab = '', ylab = '', main = '') 
legend("bottom",
       legend = c("Phycocyanin",
                  "Chlorophyll"),
       pch = 16,
       col = c("blue", "green"),
       cex = axis_tick_label_size - 0.1,
       bty = 'n',
       pt.cex = 1.5,
       y.intersp = 2)
# text(x = -0.4, y = 0.9,
#      pos = 4,
#      labels = "Pigment\ntype",
#      cex = axis_tick_label_size + 0.1)
screen(6)
pigment_plot(start_date = "2021-05-01",
             end_date = "2021-11-15")




#### Generate DO plots ####
screen(7)
DO_heatmap(2020,
           start_date = "2020-05-01",
           end_date = "2020-11-10")

screen(8)
legend_image <- as.raster(matrix(DO_color_object(20)[20:1], ncol=1))
par(mar = c(0, 0, 0, 0))
plot(c(0,3),
     c(0,1),
     type = 'n',
     axes = F,
     xlab = '',
     ylab = '',
     main = '') 
text(x = 1.5,
     y = seq(0, 0.8, l = 7),
     labels = seq(DO_range_to_plot[1], DO_range_to_plot[2], l = 7),
     cex = axis_tick_label_size)
text(x = -0.4, y = 0.9,
     pos = 4,
     labels = "DO (mg/L)",
     cex = axis_tick_label_size + 0.1)
rasterImage(legend_image, 0, 0, 1,0.8)


screen(9)
DO_heatmap(2021,
           start_date = "2021-05-01",
           end_date = "2021-11-10")


#### Generate headers ####
header_x <- 0
header_y <- 0.3
header_adj <- c(0, 0)
screen(10)
text(x = header_x, y = header_y, adj = header_adj, "A.")
screen(11)
text(x = header_x, y = header_y, adj = header_adj, "B.")
screen(12)
text(x = header_x, y = header_y, adj = header_adj, "C.")
screen(13)
text(x = header_x, y = header_y, adj = header_adj, "D.")
screen(14)
text(x = header_x, y = header_y, adj = header_adj, "E.")
screen(15)
text(x = header_x, y = header_y, adj = header_adj, "F.")

dev.off()

