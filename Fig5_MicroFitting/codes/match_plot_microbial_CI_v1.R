library(brms)
library(posterior)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(data.table)
library(zoo)

# Load microbial models
rds_path <- "rds/microbial_v1.2.rds"
results <- readRDS(rds_path)
version_tag <- tools::file_path_sans_ext(basename(rds_path))

# Load empirical data
data <- fread("data/processed_microbial.csv")
groups <- split(data, data$RecordID)

# Display group map
dataset_group_map <- c(
  "FE32" = "Lab",
  "Deng_2021" = "Biomass",
  "Monod_1941" = "Biomass",
  "Varma_1994" = "Biomass",
  "Oscar_1999_batchBo54" = "Oscar (1999)",
  "Oscar_1999_batchBo62" = "Oscar (1999)",
  "Oscar_1999_batchBo63" = "Oscar (1999)",
  "Oscar_1999_batchBo92" = "Oscar (1999)",
  "Oscar_1999_batchBo93" = "Oscar (1999)",
  "Oscar_1999_batchBo99" = "Oscar (1999)"
)

color_map <- c(
  "Lab" = "red",
  "Biomass" = "#B0D0E3",
  "Oscar (1999)" = "#B8D8C0"
)

# Function to compute θ-logistic ribbon per model
get_curve_draws_rescaled <- function(model, dataset_name) {
  draws <- as_draws_df(model, variables = c("b_r_Intercept", "b_a_Intercept"))
  if (nrow(draws) == 0) return(NULL)
  
  draws <- draws %>%
    mutate(draw_id = row_number()) %>%
    filter(is.finite(b_r_Intercept), is.finite(b_a_Intercept))
  
  x_seq <- seq(0, 1, length.out = 200)
  
  expand_grid(draw_id = draws$draw_id, x = x_seq) %>%
    left_join(draws, by = "draw_id") %>%
    mutate(
      sign_a = ifelse(b_a_Intercept >= 0, 1, -1),
      growth = b_r_Intercept * sign_a * (1 - x^b_a_Intercept),
      Group = dataset_group_map[[dataset_name]],
      Dataset = dataset_name
    ) %>%
    filter(is.finite(growth), growth > -5, growth < 10) %>%
    group_by(x) %>%
    summarise(
      mean_raw = mean(growth),
      lower_raw = quantile(growth, 0.025),
      upper_raw = quantile(growth, 0.975),
      .groups = "drop"
    ) %>%
    mutate(
      Growth_Rate = rollapply(mean_raw, width = 7, FUN = mean, fill = NA, align = "center"),
      Growth_Lower = rollapply(lower_raw, width = 7, FUN = mean, fill = NA, align = "center"),
      Growth_Upper = rollapply(upper_raw, width = 7, FUN = mean, fill = NA, align = "center"),
      Rescaled_Density = x,
      Group = dataset_group_map[[dataset_name]],
      Dataset = dataset_name
    )
}

# Build CI-fitted curves
curve_data <- map_dfr(names(results), function(dataset) {
  model <- results[[dataset]]
  if (is.null(model)) return(NULL)
  if (!dataset %in% names(dataset_group_map)) return(NULL)
  get_curve_draws_rescaled(model, dataset)
})
# Add Group_Label to both datasets
curve_data <- curve_data %>%
  mutate(Group_Label = label_map[Dataset])

# Extract empirical points
point_data <- data %>%
  filter(RecordID %in% names(dataset_group_map)) %>%
  mutate(
    Dataset = RecordID,
    Group = dataset_group_map[RecordID]
  ) %>%
  filter(is.finite(Rescaled_Density), is.finite(Per_Capita_Growth)) %>%
  select(Dataset, Rescaled_Density, Per_Capita_Growth, Group)

# Plot with ribbon
ggplot() +
  geom_ribbon(data = curve_data, aes(x = Rescaled_Density, ymin = Growth_Lower, ymax = Growth_Upper, fill = Group, group = Dataset), alpha = 0) +
  geom_line(data = filter(curve_data, Group != "Lab"),
            aes(x = Rescaled_Density, y = Growth_Rate, color = Group_Label, group = Dataset),
            linewidth = 1) +
  geom_line(data = filter(curve_data, Group == "Lab"),
            aes(x = Rescaled_Density, y = Growth_Rate, color = Group),
            linewidth = 1) +
  geom_point(data = point_data, aes(x = Rescaled_Density, y = Per_Capita_Growth, color = Group, group = Dataset), 
             alpha = 0.5, size = 1.5, stroke = 0) +
  scale_color_manual(values = color_map) +
  scale_fill_manual(values = color_map) +
  coord_cartesian(ylim = c(0, 5)) +
  labs(
    x = "Rescaled Density",
    y = "Per Capita Growth Rate (1/hr)",
    title = "Microbial data fitting",
    color = "Data type",
    fill = "Data type"
  ) +
  theme_minimal() +
  theme(
    legend.position = c(1, 1),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = alpha("white", 0.8), color = NA)
  )

ggsave(paste0("plots/fit_", version_tag, ".jpg"), width = 5, height = 5)
ggsave(paste0("plots/fit_", version_tag, ".pdf"), width = 5, height = 5)

