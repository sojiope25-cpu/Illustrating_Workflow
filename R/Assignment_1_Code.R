
# 01 Load libraries + datasets ----

library(tidyverse)
library(viridis)
library(ggtext) 
library(ggplot2)

t_spiralis <- read_csv("../data/GCA_001447595.1.csv")
t_britovi <- read_csv("../data/GCA_001447585.1.csv")

# testing something now---- 
#excluding non-TE matching classes, and adding row numbers to create unique identifier for each row, enabling exclusion of outliers

ls_non_te <- c("Simple_repeat", "Low_complexity", "Satellite", "snRNA", "rRNA", "ARTEFACT", "Unknown", "simple_repeat")

t_spiralis <- t_spiralis |> mutate(row = 1:n()) |> filter(matching_class %notin% ls_non_te)
dim(t_spiralis) #1037

t_britovi <- t_britovi |> mutate(row = 1:n()) |> filter(matching_class %notin% ls_non_te)
dim(t_britovi) #1177


# 02 Identify, Plot, and Remove Outliers (using interquartile range method) 
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

# 03 Perform Wilcoxon tests: percent divergence between species, within each TA mathcing_class ----

te_combined <- bind_rows(
  t_spiralis |> mutate(species = "spiralis"),
  t_britovi  |> mutate(species = "britovi")
) |>
  filter(
    !matching_class %in% ls_non_te,
    !is.na(matching_class),
    is.finite(perc_div)
  )

wilcox_results <- te_combined |>
  group_by(matching_class) |>
  group_modify(~ {
    x <- .x$perc_div[.x$species == "spiralis"]
    y <- .x$perc_div[.x$species == "britovi"]
    
    # Skip classes with fewer than two observations in either species
    # or where every value is identical
    can_test <- length(x) >= 2 &&
      length(y) >= 2 &&
      n_distinct(c(x, y)) > 1
    
    tibble(
      n_spiralis = length(x),
      n_britovi = length(y),
      median_spiralis = if (length(x) > 0) median(x) else NA_real_,
      median_britovi = if (length(y) > 0) median(y) else NA_real_,
      p_value = if (can_test) {
        wilcox.test(
          x, y,
          alternative = "two.sided",
          paired = FALSE,
          exact = FALSE
        )$p.value
      } else NA_real_
    )
  }) |>
  ungroup() |>
  mutate(
    p_adjusted = p.adjust(p_value, method = "BH"),
    significant = p_adjusted < 0.05
  ) |>
  arrange(p_adjusted)

print(wilcox_results, n = Inf)

# 04 Plot average percent divergence within matching_class between t_spiralis & t_britovi using boxplots ----


# find distribution of %div of classes in britovi
shared_classes <- intersect(t_spiralis$matching_class, t_britovi$matching_class)

te_shared <- bind_rows(
  t_britovi  |> 
    select(matching_class, perc_div) |> 
    mutate(species = "T. britovi"),
  t_spiralis |> 
    select(matching_class, perc_div) |> 
    mutate(species = "T. spiralis")) |>
  filter(matching_class %in% shared_classes)

te_shared |> filter(species == "T. britovi") |>
  ggplot(aes(x = perc_div, y = matching_class)) + geom_boxplot()

te_shared |> filter(species == "T. spiralis") |>
  ggplot(aes(x = perc_div, y = matching_class)) + geom_boxplot()

# combine datasets to make combined graph
ggplot(data = te_shared, mapping = aes(x = perc_div, y = matching_class, fill = species)) +
  geom_boxplot() +
  labs(
    title = "% Divergence in Shared TE Classes of T. britovi and T. spiralis",
    x = "% Divergence",
    y = "TE Class")
