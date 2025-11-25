.libPaths(c(.libPaths(), "/drives/4tb/modules/R"))
library(data.table)
library(brms)
library(dplyr)
library(tidyr)
library(future)
library(future.apply)

# Set up parallel execution: 7 workers for 7 datasets
plan(multisession, workers = 7)

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

# Nonlinear model formula
nl_formula <- bf(
  growth_rate ~ r * (2 * step(a) - 1)  * (1 - (x_gm / K)^a), 
  r + K + a ~ 1, 
  nl = TRUE
)

# Adjusted wide priors (allowing negative space)
priors <- c(
  prior(normal(0, 2), nlpar = "r", lb = 0),       # Half-normal for r
  prior(normal(500, 300), nlpar = "K", lb = 1e-3),   # Log-normal to fit the distribution of x_gm
  prior(normal(0, 2), nlpar = "a")                # Normal for a
)

# Run the models in parallel
results <- future_lapply(dataset_names, function(dataset) {
  focal_data <- groups[[dataset]] %>%
    mutate(
      x_gm = as.numeric(x_gm),
      y_gm = as.numeric(y_gm),
      growth_rate = y_gm / x_gm
    ) %>%
    drop_na(x_gm, y_gm, growth_rate) %>%
    filter(is.finite(growth_rate))
  
  if (nrow(focal_data) == 0) return(NULL)  # Skip empty datasets
  
  brm(
    formula = nl_formula,
    data = focal_data,
    prior = priors,
    control = list(adapt_delta = 0.999, max_treedepth = 25),
    iter = 10000, warmup = 2000, chains = 4, cores = 4,
    threads = threading(4),
    backend = "cmdstanr",
    refresh = 500
  )
}, future.seed = TRUE)

# Save the results to disk
names(results) <- dataset_names  # Assign dataset names for easy reference
saveRDS(results, file = "rds/Hatton_v1.2.rds")
