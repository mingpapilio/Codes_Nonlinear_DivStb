library(posterior)
library(dplyr)
library(purrr)
library(tools)

# Load your model results
rds_path <- "rds/Hatton_v1.2.rds"
results <- readRDS(rds_path)
version_tag <- file_path_sans_ext(basename(rds_path))

# Compute summary table with CIs
param_summary <- map_dfr(names(results), function(dataset) {
  model <- results[[dataset]]
  
  if (is.null(model)) {
    return(data.frame(
      Dataset = dataset,
      mean_a = NA, lower_a = NA, upper_a = NA,
      mean_r = NA, lower_r = NA, upper_r = NA,
      mean_K = NA, lower_K = NA, upper_K = NA
    ))
  }
  
  draws <- as_draws_df(model, variables = c("b_a_Intercept", "b_r_Intercept", "b_K_Intercept"))
  
  a_vals <- draws[["b_a_Intercept"]]
  r_vals <- draws[["b_r_Intercept"]]
  K_vals <- draws[["b_K_Intercept"]]
  
  valid <- is.finite(a_vals) & is.finite(r_vals) & is.finite(K_vals)
  a_vals <- a_vals[valid]
  r_vals <- r_vals[valid]
  K_vals <- K_vals[valid]
  
  r_true_vals <- r_vals / abs(a_vals)
  
  data.frame(
    Dataset = dataset,
    mean_a = mean(a_vals), lower_a = quantile(a_vals, 0.025), upper_a = quantile(a_vals, 0.975),
    mean_r = mean(r_true_vals), lower_r = quantile(r_true_vals, 0.025), upper_r = quantile(r_true_vals, 0.975),
    mean_K = mean(K_vals), lower_K = quantile(K_vals, 0.025), upper_K = quantile(K_vals, 0.975)
  )
})

# Save summary table
write.csv(param_summary, paste0("data/mean_", version_tag, ".csv"), row.names = FALSE)