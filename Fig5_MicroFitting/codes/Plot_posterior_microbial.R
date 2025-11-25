library(ggplot2)
library(ggridges)
library(posterior)
library(dplyr)
library(purrr)

# Load results if not yet loaded
rds_path <- "rds/microbial_v1.2.rds"
results <- readRDS(rds_path)

# Extract posterior draws of 'a' parameter
posterior_a <- map_dfr(names(results), function(dataset) {
  model <- results[[dataset]]
  if (is.null(model)) return(data.frame())
  
  draws <- as_draws_df(model)
  data.frame(
    a = draws[["b_a_Intercept"]],
    Dataset = dataset
  )
})

posterior_means <- posterior_a %>%
  group_by(Dataset) %>%
  summarise(mean_a = mean(a, na.rm = TRUE))

ggplot(posterior_a, aes(x = a)) +
  geom_density(fill = "skyblue", alpha = 0.6) +
  geom_vline(data = posterior_means, aes(xintercept = mean_a), 
             color = "black", linetype = "dashed", size = 0.8) +
  geom_text(data = posterior_means, 
            aes(x = 25, y = Inf,  # Adjust x to the right edge of coord_cartesian
                label = paste0("Mean = ", round(mean_a, 2))), 
            hjust = 1, vjust = 1.5, size = 3.5, color = "black") +
  facet_wrap(~ Dataset, scales = "free_y", ncol = 3) +
  coord_cartesian(xlim = c(-5, 25)) +
  theme_minimal() +
  labs(x = "Posterior of \u03B1", y = "Density", title = "Fitting microbial monocultures")

version_tag <- tools::file_path_sans_ext(basename(rds_path))
ggsave(paste0("plots/alpha_", version_tag, ".jpg"), width = 9, height = 7)

