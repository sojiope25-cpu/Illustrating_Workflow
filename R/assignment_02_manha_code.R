#Contributors and Contributions:
#
# Opeyemi – Original author of the script.
#
# Manha Khan – Secondary contributor.Revised the outlier detection procedure by creating a reusable function, updated the Wilcoxon tests and boxplots to consistently use outlier-removed datasets, and added a sensitivity analysis comparing results with and without outliers.
# 01 Load libraries + datasets ----

library(tidyverse)
library(viridis)
library(ggtext)
library(ggplot2)

t_spiralis <- read_csv("data/GCA_001447595.1.csv")
t_britovi <- read_csv("data/GCA_001447585.1.csv")

# testing something now---- 
#excluding non-TE matching classes, and adding row numbers to create unique identifier for each row, enabling exclusion of outliers

ls_non_te <- c("Simple_repeat", "Low_complexity", "Satellite", "snRNA", "rRNA", "ARTEFACT", "Unknown", "simple_repeat")

t_spiralis <- t_spiralis |> mutate(row = 1:n()) |> filter(matching_class %notin% ls_non_te)
dim(t_spiralis) #1037

t_britovi <- t_britovi |> mutate(row = 1:n()) |> filter(matching_class %notin% ls_non_te)
dim(t_britovi) #1177

# 02 REVISED Identify, Plot, and Remove Outliers REVISED (using interquartile range method) ----

# Function to identify outliers and remove them

identify_outliers <- function(data) {
  Q1 <- quantile(data$perc_div, 0.25, na.rm = TRUE)
  Q3 <- quantile(data$perc_div, 0.75, na.rm = TRUE)
  IQR_value <- IQR(data$perc_div, na.rm = TRUE)
  
  lower_limit <- Q1 - 1.5 * IQR_value
  upper_limit <- Q3 + 1.5 * IQR_value
  
  # Identify outliers
  outliers <- subset(data,perc_div < lower_limit | perc_div > upper_limit)
  # Remove outliers from dataset
  cleaned_data <- data |> 
    filter(!row %in% outliers$row)
  # Return results
  list(
    outliers = outliers,
    cleaned_data = cleaned_data,
    lower_limit = lower_limit,
    upper_limit = upper_limit
  )
}
# Apply function to both species
spiralis_results <- identify_outliers(t_spiralis)
britovi_results <- identify_outliers(t_britovi)

# Extract outliers and cleaned datasets
outliers_spiralis <- spiralis_results$outliers
outliers_britovi <- britovi_results$outliers

t_spiralis_out_removed <- spiralis_results$cleaned_data
t_britovi_out_removed <- britovi_results$cleaned_data

# Check number of outliers
dim(outliers_spiralis) #7
dim(outliers_britovi) #1

# Check cleaned dataset dimensions
dim(t_spiralis_out_removed) #1030

dim(t_britovi_out_removed) #1176
# Plot outlier detection for T. spiralis
t_spiralis |>
  ggplot(aes(matching_class, y = perc_div)) +geom_boxplot() + geom_hline(aes(yintercept = spiralis_results$lower_limit,colour = "Lower Limit")) + geom_hline(aes(yintercept = spiralis_results$upper_limit,colour = "Upper Limit")) + scale_color_manual( name = "Outlier Thresholds", values = c("Lower Limit" = "red", "Upper Limit" = "blue")) +guides(x = guide_axis(angle = 45)) +labs( title = "Outlier Detection by TE Matching Class - *T. spiralis*", x = "TE Matching Class",y = "Percent Divergence") +theme(plot.title = element_markdown())

# Plot outlier detection for T. britovi
t_britovi |>
  ggplot(aes(matching_class, y = perc_div)) + geom_boxplot() + geom_hline(aes(
    yintercept = britovi_results$lower_limit, colour = "Lower Limit")) + geom_hline(aes( yintercept = britovi_results$upper_limit, colour = "Upper Limit" )) + scale_color_manual(name = "Outlier Thresholds", values = c("Lower Limit" = "red", "Upper Limit" = "blue") ) + guides(x = guide_axis(angle = 45)) + labs(title = "Outlier Detection by TE Matching Class - *T. britovi*", x = "TE Matching Class", y = "Percent Divergence" ) + theme(plot.title = element_markdown())