.libPaths(c(.libPaths(), "/drives/4tb/modules/R"))
library(brms)
library(dplyr)
library(tidyr)
library(future)
library(future.apply)
library(data.table)  # For fread()

# Load merged data
raw <- fread("data/processed_microbial.csv")

# Get unique trials based on RecordID
unique_trials <- unique(raw$RecordID)

# Dynamically set parallel workers to the number of unique trials (capped if needed)
plan(multisession, workers = length(unique_trials))

# Updated nonlinear model formula
nl_formula <- bf(
  ## (2 * step(a) - 1) operates the same as sign (a); and the 'true' r= r*a
  Per_Capita_Growth_Rate ~ r * (2 * step(a) - 1) * (1 - Rescaled_Density^a), 
  r + a ~ 1, 
  nl = TRUE
)

# Updated priors
priors <- c(
  prior(normal(0, 2), nlpar = "r", lb = 0),       # Half-normal for r
  prior(normal(0, 2), nlpar = "a")                # Normal for a
)

# Run models in parallel for each unique trial
results <- future_lapply(unique_trials, function(trial_id) {
  focal_data <- raw %>%
    filter(RecordID == trial_id) %>%
    mutate(
      Rescaled_Density = as.numeric(Rescaled_Density),
      Per_Capita_Growth_Rate = as.numeric(Per_Capita_Growth)
    ) %>%
    drop_na(Rescaled_Density, Per_Capita_Growth_Rate) %>%
    filter(is.finite(Per_Capita_Growth_Rate))
  
  if (nrow(focal_data) == 0) return(NULL)  # Skip empty groups
  
  fit <- brm(
    formula = nl_formula,
    data = focal_data,
    prior = priors,
    control = list(adapt_delta = 0.999, max_treedepth = 25),
    iter = 10000, warmup = 2000, chains = 4, cores = 4,
    threads = threading(4),
    backend = "cmdstanr",
    refresh = 500
  )
  
  attr(fit, "dataset_name") <- trial_id  # Store the trial ID as metadata
  return(fit)
}, future.seed = TRUE)

# Assign model names directly from metadata
names(results) <- unique_trials  # Assign dataset names for easy reference
saveRDS(results, file = "rds/microbial_v1.2.rds")

