library(purrr)
library(dplyr)
library(tidyr)
library(future)
library(future.apply)

raw <- fread("data/fig4.csv")
tbl <- table(raw$Taxa)
subset(tbl, tbl > 30)
# Filter each group
groups <- list(
  Mammals = raw %>% filter(Taxa == "Large mammal herbivores"),
  Seagrass = raw %>% filter(Taxa == "Seagrass"),
  Fish = raw %>% filter(Taxa == "Fish"),
  Broadleaf = raw %>% filter(Taxa == "Broadleaf"),
  Coniferous = raw %>% filter(Taxa == "Needleleaf"),
  Grassland = raw %>% filter(Taxa == "Pinales: Cupressaceae" | Taxa == "Pinales: Pinaceae"),
  MarineBacteria = raw %>% filter(Taxa == "Marine sediment bacteria")
)

# Dataset names
dataset_names <- c("Mammals", "Seagrass", "Fish", "Broadleaf", "Coniferous", "Grassland", "MarineBacteria")

# Combine and clean x_gm values
selected_xgm <- groups %>%
  map(~ as.numeric(.x$x_gm)) %>%
  unlist() %>%
  .[is.finite(.)]  # Removes NA, NaN, Inf, -Inf

# Check if any usable values remain
if (length(selected_xgm) == 0) {
  stop("No valid x_gm values found in selected datasets.")
}

# Plot histogram
hist(
  selected_xgm,
  breaks = 50,
  col = "gray",
  main = "Histogram of x_gm in Selected Datasets",
  xlab = "x_gm",
  xlim = c(0, quantile(selected_xgm, 0.99, na.rm = TRUE))
)

# Summary stats
summary_stats <- quantile(selected_xgm, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE)
print(summary_stats)

# Mean and SD
cat("Mean:", mean(selected_xgm, na.rm = TRUE), "\n")
cat("SD:", sd(selected_xgm, na.rm = TRUE), "\n")

