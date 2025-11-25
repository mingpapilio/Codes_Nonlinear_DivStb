library(posterior)
library(dplyr)
library(purrr)

# Load your model results
rds_path <- "rds/microbial_v1.2.rds"
results <- readRDS(rds_path)

param_summary_ra <- map_dfr(names(results), function(dataset) {
  model <- results[[dataset]]
  
  if (is.null(model)) {
    return(data.frame(
      Dataset = dataset,
      mean_a = NA, lower_a = NA, upper_a = NA,
      mean_r = NA, lower_r = NA, upper_r = NA
    ))
  }
  
  draws <- as_draws_df(model)
  
  r_vals <- draws[["b_r_Intercept"]]
  a_vals <- draws[["b_a_Intercept"]]
  
  # Back-transform r to true intrinsic growth rate
  r_true_vals <- r_vals / abs(a_vals)
  
  data.frame(
    Dataset = dataset,
    mean_a = mean(a_vals), lower_a = quantile(a_vals, 0.025), upper_a = quantile(a_vals, 0.975),
    mean_r = mean(r_true_vals), lower_r = quantile(r_true_vals, 0.025), upper_r = quantile(r_true_vals, 0.975)
  )
})

version_tag <- tools::file_path_sans_ext(basename(rds_path))
write.csv(param_summary_ra, paste0("data/mean_", version_tag, ".csv"), row.names = FALSE)
