#### code/figures/bacterial_production.R ####
# Written for BLiMMP project
# Benjamin D. Peterson


#### Prep workspace ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/BLiMMP")
library(tidyverse)
source("code/BLiMMP_functions.R")


#### Read in data ####
leucine_data_metadata <- read.csv('dataFinal/bacterial_production.csv')
waterChem_data <- read.csv("dataFinal/water_chem_data.csv")
exo_data <- read.csv("dataFinal/exo_data.csv")


#### Transparent colors ####
trans_colors <- function(color_to_change,
                         percent_trans) {
  rgb_color_code <- col2rgb(color_to_change)
  ## Make new color using input color as base and alpha set by transparency
  rgb(rgb_color_code[1], rgb_color_code[2], rgb_color_code[3],
      max = 255,
      alpha = (100 - percent_trans) * 255 / 100)
  }
trans_blue <- trans_colors(cb.translator["blue"], 50)
trans_black <- trans_colors(cb.translator["black"], 50)
trans_pink <- trans_colors(cb.translator["reddishpurple"], 50)


#### Function: Plot exo data ####
plot.exo.data <- function(exo.data.to.use,
                          date.of.sampling,
                          points.of.sampling = NULL,
                          # main.title = NULL,
                          legend.location = "topright") {
  
  
  # Isolate needed exo data
  exo.data.for.profile <- exo.data.to.use %>%
    filter(date == date.of.sampling)
  
  # Plot the temp points and set up graph
  plot(x = (exo.data.for.profile$Temp_C-5)/2.5,
       y = -exo.data.for.profile$depth,
       type = "l",
       lwd = 3,
       xaxt = "n",
       xlab = "",
       xlim = c(0, 10),
       yaxt = "n",
       ylab = "",
       ylim = c(-25, 2),
       col = trans_blue)
  # Plot DO
  DO.fudge.factor <- 15
  lines(x = exo.data.for.profile$ODO_sat/DO.fudge.factor,
        y = -exo.data.for.profile$depth,
        col = trans_black,
        lwd = 3,
        lty = 4)
  # Plot the turbidity values
  turb.fudge.factor <- 0.5
  lines(x = exo.data.for.profile$Turbidity_FNU/turb.fudge.factor,
        y = -exo.data.for.profile$depth,
        col = trans_pink,
        lwd = 4,
        lty = 3)
  # Add y-axis with depth measurements
  par(mgp = c(1.5, 0.4, 0))
  axis(2,
       at = seq(0, -25, by = -5),
       labels = seq(0, 25, by = 5))
  # # Add axis for temperature
  # axis(1,
  #      at = seq(0, 10, by = 2),
  #      labels = seq(0, 10, by = 2)*3)
  # # Add axis for DO values
  # axis(1,
  #      at = seq(0, 10, by = 2),
  #      labels = seq(0, 10, by = 2)*DO.fudge.factor,
  #      line = 3,
  #      par(mgp = c(3, 0.5, 0)))
  # # Add axis for turbidity values
  # axis(1,
  #      at = seq(0, 10, by = 2),
  #      labels = seq(0, 10, by = 2)*turb.fudge.factor,
  #      line = 6,
  #      par(mgp = c(3, 0.5, 0)))
  # 
  # 
  # # Add titles for depth
  # title(xlab = "Temp (C)",
  #       ylab = "Depth (m)",
  #       line = 1.5)
  # # Add title for dissolved oxygen
  # title(xlab = "ODO (% sat)",
  #       line = 4.5)
  # # Add title for turbidity
  # title(xlab = "Turbidity (FNU)",
  #       line = 7.5)
}



#### Function: Bacterial production profile by date ####
plot_BP_profiles_by_date <- function(date_to_use) {
  
  plot.exo.data(exo.data.to.use = exo_data,
                date.of.sampling = date_to_use)
  
  data_to_plot <- leucine_data_metadata %>%
    filter(startDate == date_to_use,
           treatment == "ambient",
           (timePoint == "t1" | is.na(timePoint)))
  points(x = data_to_plot$µgBCP_per_L_hr / 2.5,
         y = -data_to_plot$depth,
         pch = 21,
         bg = cb.translator["yellow"],
         cex = 1.2)
  # data_to_plot_Mo <- leucine_data_metadata %>%
  #   filter(startDate == date_to_use,
  #          treatment == "molybdate",
  #          (timePoint == "t1" | is.na(timePoint)))
  # points(x = data_to_plot_Mo$µgBCP_per_L_hr / 2.5,
  #        y = -data_to_plot_Mo$depth,
  #        pch = 4,
  #        col = cb.translator["bluishgreen"],
  #        lwd = 1.8,
  #        cex = 1.2)
  par(mgp = c(1.5, 0.2, 0))
  axis(1, at = seq(0, 10, by = 2),
       labels = seq(0, 10, by = 2)*2.5,
       par(mgp = c(3, 0.5, 0)),
       gap.axis = -1)
  title(xlab = "BP (µgC/L/hr)",
        ylab = "Depth (m)",
        line = 1.5)
  
}



#### Function: Plot bacterial production by concentration of molybdate amendment ####
plot_BP_by_molybdate_conc <- function() {

  data_to_plot <- leucine_data_metadata %>%
    filter(startDate == "2020-09-17",
           !grepl("glucose", treatment)) %>%
    group_by(depth, treatment) %>%
    summarize(mean_µgBCP_per_L_hr = mean(µgBCP_per_L_hr),
              n_group = n(),
              se_µgBCP_per_L_hr = sd(µgBCP_per_L_hr) / sqrt(n_group)) %>%
    ungroup()
  
  # Set up means data
  data_table <- data_to_plot %>%
    select(depth, treatment, mean_µgBCP_per_L_hr) %>%
    spread(value = mean_µgBCP_per_L_hr,
           key = depth) %>%
    as.data.frame()
  rownames(data_table) <- data_table$treatment
  data_matrix <- as.matrix(data_table[, c("14.1", "20.1")])
  
  # Set up error bars
  mid_x <- barplot(data_matrix,
                   beside = TRUE,
                   plot = FALSE)
  data_table_se <- data_to_plot %>%
    select(depth, treatment, se_µgBCP_per_L_hr) %>%
    spread(value = se_µgBCP_per_L_hr,
           key = depth) %>%
    as.data.frame()
  rownames(data_table_se) <- data_table_se$treatment
  data_matrix_se <- as.matrix(data_table_se[, c("14.1", "20.1")])
  

  color_vector <- cb.translator[c("yellow", "skyblue", "bluishgreen", "blue")]
  names(color_vector) <- c("ambient", "molybdate-0.017", "molybdate-0.17", "molybdate-1.7")
  name_vector <- c("Ambient", "Molybdate (17 µM)", "Molybdate (170 µM)", "Molybdate (1.7 mM)")
  names(name_vector) <- c("ambient", "molybdate-0.017", "molybdate-0.17", "molybdate-1.7")
  
  par(mgp = c(1.5, 0.2, 0))
  barplot(data_matrix,
          xlab = "Depth (m)",
          ylab = "BP (µgC/L/hr)",
          col = color_vector[rownames(data_table)],
          ylim = c(0, 12),
          beside = TRUE) # Grouped bars
  
  arrows(x0 = mid_x, x1 = mid_x, 
         y0 = data_matrix - data_matrix_se,
         y1 = data_matrix + data_matrix_se,
         code = 3,
         angle = 90,
         length = 0.1)
  
  legend(x = 0.8,
         y = 11.75,
         legend = name_vector,
         fill = color_vector,
         cex = 0.7)
  box()
  }


#### Function: Plot bacterial production by location, with and without molybdate ####
plot_BP_molybdate_comparison <- function() {
  
  data_to_plot_2021 <- leucine_data_metadata %>%
    filter(year(startDate) == "2021") %>%
    mutate(location = paste(startDate, "\n",
                            depth, "m",
                            sep = "")) %>%
    group_by(location, treatment) %>%
    summarize(mean_µgBCP_per_L_hr = mean(µgBCP_per_L_hr),
              n_group = n(),
              se_µgBCP_per_L_hr = sd(µgBCP_per_L_hr) / sqrt(n_group)) %>%
    ungroup()
  
  # Set up means data
  data_table_mean_2021 <- data_to_plot_2021 %>%
    select(location, treatment, mean_µgBCP_per_L_hr) %>%
    spread(value = mean_µgBCP_per_L_hr,
           key = location) %>%
    as.data.frame()
  rownames(data_table_mean_2021) <- data_table_mean_2021$treatment
  data_matrix_mean_2021 <- as.matrix(data_table_mean_2021[, -1])
  
  # Set up error bars
  mid_x_positions_2021 <- barplot(data_matrix_mean_2021,
                                  beside = TRUE,
                                  plot = FALSE)
  data_table_se <- data_to_plot_2021 %>%
    select(location, treatment, se_µgBCP_per_L_hr) %>%
    spread(value = se_µgBCP_per_L_hr,
           key = location) %>%
    as.data.frame()
  rownames(data_table_se) <- data_table_se$treatment
  data_matrix_se <- as.matrix(data_table_se[, -1])
  
  color_vector <- cb.translator[c("yellow", "bluishgreen")]
  names(color_vector) <- c("ambient", "molybdate")
  name_vector <- c("Ambient", "Molybdate (170 µM)")
  names(name_vector) <- names(color_vector)
  
  par(mgp = c(1.5, 1.5, 0))
  barplot(data_matrix_mean_2021,
          ylab = "BP (µgC/L/hr)",
          col = color_vector[rownames(data_matrix_mean_2021)],
          ylim = c(0, 22),
          beside = TRUE,
          yaxt = "n")
  par(mgp = c(1.5, 0.2, 0))
  axis(2)
  
  arrows(x0 = mid_x_positions_2021, x1 = mid_x_positions_2021, 
         y0 = data_matrix_mean_2021 - data_matrix_se,
         y1 = data_matrix_mean_2021 + data_matrix_se,
         code = 3,
         angle = 90,
         length = 0.1)
  
  legend(x = 13,
         y = 19,
         legend = name_vector,
         fill = color_vector,
         cex = 0.8)
  box()
}


#### Set up plot areas ####
top_row_vertArea <- c(0.67, 1.00)
middle_row_vertArea <- c(0.33, 0.67)
bottom_row_vertArea <- c(0.00, 0.33)

left_horiArea <- c(0.00, 0.33)
center_horiArea <- c(0.33, 0.67)
right_horiArea <- c(0.67, 1.00)
right_and_center_horiArea <- c(0.33, 1.00)
all_horiArea <- c(0.00, 1.00)



#### Generate PDF ####
cairo_pdf("results/figures/BP_plots.pdf",
          width = 7,
          height = 9)
par(mar = c(3, 3, 1, 0.5),
    tck = -0.008)
split.screen(rbind(c(left_horiArea, top_row_vertArea), c(center_horiArea, top_row_vertArea), c(right_horiArea, top_row_vertArea),
                   c(left_horiArea, middle_row_vertArea), c(right_and_center_horiArea, middle_row_vertArea),
                   c(all_horiArea, bottom_row_vertArea)))

#### Plot data ####
screen(1)
# par(mgp = c(1.5, 0.2, 0))
plot_BP_profiles_by_date("2020-10-10")
mtext(text = "A. Oct. 10, 2020",
      side = 3,
      adj = 0)
legend(x = 2.2,
       y = 2,
       legend = c("Bacterial Prod.", "Temperature", "Diss. Oxygen", "Turbidity"),
       pch = c(21, NA, NA, NA),
       lty = c(NA, 1, 4, 3),
       col = c("black", trans_colors(cb.translator["blue"], 50),
               trans_colors(cb.translator["black"], 50),
               trans_colors(cb.translator["reddishpurple"], 50)),
       pt.bg = c(cb.translator["yellow"], NA, NA, NA),
       cex = 0.7,
       pt.cex = 1.2)
screen(2)
# par(mgp = c(1.5, 0.2, 0))
plot_BP_profiles_by_date("2021-09-10")
mtext(text = "B. Sept. 10, 2021",
      side = 3,
      adj = 0)
screen(3)
# par(mgp = c(1.5, 0.2, 0))
plot_BP_profiles_by_date("2021-10-14")
mtext(text = "C. Oct. 14, 2021",
      side = 3,
      adj = 0)
screen(4)
# par(mgp = c(1.5, 0.2, 0))
plot_BP_profiles_by_date("2020-09-17")
mtext(text = "D. Sept. 17, 2020",
      side = 3,
      adj = 0)

screen(5)
# par(mgp = c(1.5, 0.2, 0))
plot_BP_by_molybdate_conc()
mtext(text = "E. Molybdate treatments",
      side = 3,
      adj = 0)

screen(6)
plot_BP_molybdate_comparison()
mtext(text = "F. Molybdate effects on BP",
      side = 3,
      adj = 0)


dev.off()

